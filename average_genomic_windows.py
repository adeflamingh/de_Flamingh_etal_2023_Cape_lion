#!/usr/bin/env python3
import sys, os, argparse, datetime, gzip, statistics

#
# Constants
#
PROG = sys.argv[0].split('/')[-1]

#
# Command line options
#
def parse_args():
    p = argparse.ArgumentParser(prog=PROG)
    p.add_argument('-f', '--fai', required=True,
                   help='Path to the genome index (FAI).')
    p.add_argument('-v', '--var-sites', required=False,
                   help='Path to MAF file with variant site positions.')
    p.add_argument('-d', '--depth', required=False,
                   help='Path to SAMtools depth file.')
    p.add_argument('-o', '--outdir', required=False, default='.',
                   help='Path to output directory.')
    p.add_argument('-m', '--min-length', required=False, default=1_000_000,
                   type=float, help='Minimum chromosome size in bp [default 1000000 (1Mbp)]')
    p.add_argument('-s', '--window-size', required=False, default=250_000,
                   type=float, help='Size of windows in bp [default 250000 (250Kbp)].')
    p.add_argument('-t', '--window-step', required=False, default=100_000,
                   type=float, help='Step of windows in bp [default 100000 (100Kbp)].')
    # Check inputs
    args = p.parse_args()
    args.outdir = args.outdir.rstrip('/')
    if not os.path.exists(args.outdir):
        sys.exit(f"Error: output directory does not exist.\n\t{args.outdir}")
    if not os.path.exists(args.fai):
        sys.exit(f"Error: genome index does not exist.\n\t{args.fai}")
    if args.var_sites is not None:
        if not os.path.exists(args.var_sites):
            sys.exit(f"Error: variant sites file does not exist.\n\t{args.var_sites}")
    if args.depth is not None:
        if not os.path.exists(args.depth):
            sys.exit(f"Error: depth file does not exist.\n\t{args.depth}")
    if not args.window_size > 0:
        sys.exit(f"Error: size of windows ({args.window_size}) must be > 0.")
    if not args.window_step > 0:
        sys.exit(f"Error: step of windows ({args.window_step}) must be > 0.")
    if not args.min_length > 0:
        sys.exit(f"Error: Min chromosome length ({args.min_length}) must be > 0.")
    if not args.window_size >= args.window_step:
        sys.exit(f"Error: Window size ({args.window_size}) must be >=  than window step ({args.window_step}).")
    # Adjust the numerica values
    args.min_length  = int(args.min_length)
    args.window_size = int(args.window_size)
    args.window_step = int(args.window_step)
    return args

#
# Print the current date
#
def date_now():
    return datetime.datetime.now().strftime("%Y-%m-%d")

#
# Print current time
#
def time_now():
    return datetime.datetime.now().strftime("%H:%M:%S")

#
# Calculate the window intervals for a given Chromosome length
#
def calculate_chr_window_intervals(chr_len, window_size=250_000, window_step=100_000):
    assert window_size >= window_step
    assert window_size > 0
    windows = list()
    window_start = 0
    window_end = window_size
    while window_end < (chr_len+window_step):
        windows.append((window_start, window_end))
        window_start += window_step
        window_end += window_step
    return windows

#
# Add a given variant site to the windows
#
def populate_variant_site_windows(variant_site, chrom_id, chr_window_boundaries, populated_window_intervals, start=0):
    assert len(chr_window_boundaries) > 0
    assert 0 <= start < len(chr_window_boundaries)
    assert type(populated_window_intervals[chrom_id]) is dict
    # Index to determine the starting index in the list of chromosome
    # window boundaries. The index increases across sequential iterations
    # to prevent looping the windows from the beginning each time.
    idx = start
    # Loop over the windows and check the sites
    for window in chr_window_boundaries[start:]:
        assert len(window) == 2
        window_sta = window[0]
        window_end = window[1]
        # Check the site:
        # Skip windows that lie before the site
        if window_end < variant_site:
            # Increase the index
            idx += 1
            continue
        # Site lands in the window
        elif window_sta <= variant_site < window_end:
            # The middle of the window serves as it main reference position
            window_mid = int(window_sta + ((window_end - window_sta)/2))
            # Initialize the output dictionary
            populated_window_intervals[chrom_id].setdefault(window_mid, 0)
            # Increase the window counter
            populated_window_intervals[chrom_id][window_mid] += 1
        # Break the loop when moving past the site
        elif window_sta > variant_site:
            break
    # Return the window with the current variant added as well as the
    # starting index for the next iteration.
    return populated_window_intervals, idx

#
# Determine the windows from the genome FAI
#
def set_windows_from_fai(fai, window_size=250_000, window_step=100_000, min_chr_size=1_000_000):
    assert window_step > 0
    genome_window_intervals = dict()
    # Set the fai file handle
    fh = None
    if fai.endswith('gz'):
        fh = gzip.open(fai, 'rt')
    else:
        fh = open(fai, 'r')
    # Parse file
    for line in fh:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')
        if not len(fields) >= 2:
            sys.exit('Error: FAI must have at least two columns, seq_id<tab>seq_len')
        seq_id = fields[0]
        if not fields[1].isnumeric():
            sys.exit('Error: column two of the FAI must be a number with the length of the sequence.')
        seq_len = int(fields[1])
        if seq_len < min_chr_size:
            continue
        windows = calculate_chr_window_intervals(seq_len, window_size, window_step)
        genome_window_intervals[seq_id] = windows
    print(f'\nGenerated window intervals for {len(genome_window_intervals):,} chromosomes/scaffolds.', flush=True)
    return genome_window_intervals

#
# Initialize dictionary of windows
#
def init_windows_dictionary(genome_window_intervals, win_type='variant'):
    assert win_type in {'variant', 'depth'}
    windows_dict = dict()
    for chrom in genome_window_intervals:
        chr_windows = genome_window_intervals[chrom]
        windows_dict.setdefault(chrom, dict())
        for window in chr_windows:
            assert len(window) == 2
            window_sta = window[0]
            window_end = window[1]
            window_mid = int(window_sta + ((window_end - window_sta)/2))
            deflt = None
            if win_type == 'variant':
                deflt = 0
            else:
                deflt = [0, 0]
            windows_dict[chrom].setdefault(window_mid, deflt)
    return windows_dict

#
# Parse the variant sites file and populate windows
#
def read_variant_sites(variant_sites_f, genome_window_intervals):
    print('\nTallying variant sites per chromosome...', flush=True)
    populated_window_intervals = init_windows_dictionary(genome_window_intervals, win_type='variant')
    chr_position_tally = dict()
    fh = None
    if variant_sites_f.endswith('gz'):
        fh = gzip.open(variant_sites_f, 'rt')
    else:
        fh = open(variant_sites_f, 'r')
    prev_chr = None
    prev_pos = None
    start = None
    genome_total = 0
    for line in fh:
        if line.startswith('#'):
            continue
        if line.startswith('chromo'):
            continue
        line = line.strip('\n')
        fields = line.split('\t')
        if not len(fields) >= 2:
            sys.exit('Error: variant sites file must have >2 columns, scaffold_id<tab>position_bp<tab>major<tab>...')

        #
        # Process the chromosome information
        #
        chr_id = fields[0]
        # Skip positions from unwanted chromosomes
        if chr_id not in genome_window_intervals:
            continue
        # Check for next chromosome and clear variables
        if chr_id != prev_chr:
            if prev_chr is None:
                prev_chr = chr_id
                prev_pos = 0
                start = 0
            else:
                # Print the previous chromosome to log
                n_sites = chr_position_tally[prev_chr]
                genome_total += n_sites
                print(f'    {prev_chr} : {n_sites:,} variant sites found.', flush=True)
                # Reset to the next
                start = 0
                prev_pos = 0
                prev_chr = chr_id
        # Tally the position from the current chromosome
        chr_position_tally.setdefault(chr_id, 0)
        chr_position_tally[chr_id] += 1

        #
        # Check the BP position
        # 
        if not fields[1].isnumeric():
            sys.exit('Error: column two of the variant sites file must be an integer position.')
        position = int(fields[1])
        # Check for sorted sequence
        if not prev_pos < position:
            sys.exit('Error: depth file must be sorted by Chr position.')

        #
        # Populate the windows
        #
        chr_window_boundaries = genome_window_intervals.get(chr_id, None)
        populated_window_intervals, start = populate_variant_site_windows(position, chr_id, chr_window_boundaries, populated_window_intervals, start)
    # Report final chromosome and total
    n_sites = chr_position_tally[prev_chr]
    genome_total += n_sites
    print(f'    {prev_chr} : {n_sites:,} variant sites found.', flush=True)
    print(f'\n    Genome-wide : {genome_total:,} total variant sites found.', flush=True)
    return populated_window_intervals

#
# Print the variant site tally file
#
def print_variant_site_tally(populated_window_intervals, window_size=250_000, outdir='.'):
    # Preparte output
    outfh = open(f'{outdir}/variant_site_tally.tsv', 'w')
    # Make header
    outfh.write('Chrom\tPos\tNumSitesWindow\tFreqSitesWindow\n')
    # Loop over dictionary
    for chr in populated_window_intervals:
        chr_windows = populated_window_intervals[chr]
        for pos in sorted(chr_windows):
            tally = chr_windows[pos]
            sites_per_window = tally/window_size
            row = f'{chr}\t{pos}\t{tally}\t{sites_per_window:.08g}\n'
            outfh.write(row)
    outfh.close()

#
# Add a given coverage to the windows
#
def populate_depth_windows(chr_id, position, depth, chr_window_boundaries, populated_depth_windows, start=0):
    assert len(chr_window_boundaries) > 0
    assert 0 <= start < len(chr_window_boundaries)
    assert type(populated_depth_windows[chr_id]) is dict
    assert depth >= 0
    # Index to determine the starting index in the list of chromosome
    # window boundaries. The index increases across sequential iterations
    # to prevent looping the windows from the beginning each time.
    idx = start
    # Loop over the windows and check the position
    for window in chr_window_boundaries[start:]:
        assert len(window) == 2
        start = window[0]
        end = window[1]
        # Skip windows that end before the position of interest
        if end < position:
            # Increase the index
            idx += 1
            continue
        # Position is in the window
        elif start <= position < end:
            # Window's middle position, serves as reference
            mid = int(start + ((end - start)/2))
            # Initialize the dictionary if empty (should not happen)
            populated_depth_windows[chr_id].setdefault(mid, [0,0])
            # Increase the tally and the sum
            populated_depth_windows[chr_id][mid][0] += 1     # How many non-zero sites
            populated_depth_windows[chr_id][mid][1] += depth # The sum of non-zero sites
        # Break when moving past the position
        elif start > position:
            break
    # Return the window with the current variant added as well as the
    # starting index for the next iteration.
    return populated_depth_windows, idx

#
# Print the per-window average coverage file
#
def print_average_coverage_file(populated_depth_windows, window_size=250_000, outdir='.'):
    # Prepare output
    outfh = open(f'{outdir}/window_average_coverage.tsv', 'w')
    outfh.write('Chrom\tPos\tNonZeroCovSites\tMeanCov\n')
    # Loop over the data
    for chr in populated_depth_windows:
        chr_windows = populated_depth_windows[chr]
        for pos in sorted(chr_windows):
            win_values = chr_windows[pos]
            assert len(win_values) == 2
            non_zero_n = win_values[0]
            non_zero_sum = win_values[1]
            # Determine a mean coverage
            mean_cov = 0
            if non_zero_n > 0:
                 mean_cov = non_zero_sum/window_size
            row = f'{chr}\t{pos}\t{non_zero_n}\t{mean_cov:.08g}\n'
            outfh.write(row)
    outfh.close()

#
# Parse the SAMtools depth file
#
def parse_depth_file(depth_f, genome_window_intervals):
    print('\nTallying depth per chromosome...', flush=True)
    populated_depth_windows = init_windows_dictionary(genome_window_intervals, win_type='depth')
    chr_depth_tally = dict()
    fh = None
    if depth_f.endswith('gz'):
        fh = gzip.open(depth_f, 'rt')
    else:
        fh = open(depth_f, 'r')
    prev_chr = None
    prev_pos = None
    start = 0
    genome_total = 0
    genome_non_zero = 0
    for line in fh:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')
        if len(fields) != 3:
            sys.exit('Error: depth file does not confirm to the SAMtools depth format: scaffold_id<tab>position_bp<tab>depth')
        if not fields[1].isnumeric() or not fields[2].isnumeric():
            sys.exit('Error: columns 2 and 3 of the depth file must be an integer.')
        chr_id = fields[0]
        position = int(fields[1])
        depth = int(fields[2])

        # Skip positions from unwanted chromosomes
        if chr_id not in genome_window_intervals:
            continue

        # When processing a brand new chromosome
        if chr_id != prev_chr:
            if prev_chr is None:
                # First line in file, prepare for incoming data
                prev_chr = chr_id
                prev_pos = 0
                start = 0
            else:
                # Print the tally to the log
                total_chr_sites = chr_depth_tally[prev_chr][0]
                nonzero_chr_sites = chr_depth_tally[prev_chr][1]
                nonzero_chr_freq = nonzero_chr_sites/total_chr_sites
                print(f'    {prev_chr} :  {total_chr_sites:,} sites seen, {nonzero_chr_sites:,} ({nonzero_chr_freq:.02%}) with non-zero depth.', flush=True)
                # Reset for the next sequence
                prev_chr = chr_id
                prev_pos = 0
                start = 0

        # Check for sorted sequence
        if not prev_pos < position:
            sys.exit('Error: depth file must be sorted by Chr position.')

        # Tally the current chromosome
        chr_depth_tally.setdefault(chr_id, [0, 0])
        chr_depth_tally[chr_id][0] += 1
        genome_total += 1
        # Initialize the window (if needed)
        populated_depth_windows.setdefault(chr_id, dict())
        chr_window_boundaries = genome_window_intervals.get(chr_id, None)
        # For efficiency, skip sites with no coverage (and tally)
        if depth < 1:
            continue
        chr_depth_tally[chr_id][1] += 1
        genome_non_zero += 1
        # Populate the windows with > 0 values only
        populated_depth_windows, start = populate_depth_windows(chr_id, position, depth, chr_window_boundaries, populated_depth_windows, start)

    # Print the tally of the last chr in file to the log
    total_chr_sites = chr_depth_tally[prev_chr][0]
    nonzero_chr_sites = chr_depth_tally[prev_chr][1]
    nonzero_chr_freq = nonzero_chr_sites/total_chr_sites
    print(f'    {prev_chr} :  {total_chr_sites:,} sites seen, {nonzero_chr_sites:,} ({nonzero_chr_freq:.02%}) with non-zero depth.', flush=True)
    # Print tally of the genome-wide stats
    genome_non_zero_f = genome_non_zero/genome_total
    print(f'\n    Genome-wide : {genome_total:,} sites seen, {genome_non_zero:,} ({genome_non_zero_f:.02%}) with non-zero depth.', flush=True)
    return populated_depth_windows

def main():
    args = parse_args()
    # Initialize script
    print(f'Started {PROG} on {date_now()} {time_now()}.')
    print(f'    Min Chrom Size (bp): {args.min_length:,}')
    print(f'    Window Size (bp): {args.window_size:,}')
    print(f'    Window Step (bp): {args.window_step:,}', flush=True)
    # Get windows from the fai
    genome_window_intervals = set_windows_from_fai(args.fai, args.window_size, args.window_step, args.min_length)
    # Parse the variant sites file and print output
    if args.var_sites is not None:
        populated_window_sites = read_variant_sites(args.var_sites, genome_window_intervals)
        print_variant_site_tally(populated_window_sites, args.window_size, args.outdir)
        del populated_window_sites
    # Parse the depth file and print output
    if args.depth is not None:
        populated_depth_windows = parse_depth_file(args.depth, genome_window_intervals)
        print_average_coverage_file(populated_depth_windows, args.window_size, args.outdir)
        del populated_depth_windows

    print(f'\nFinished {PROG} on {date_now()} {time_now()}.', flush=True)

# Run Code
if __name__ == '__main__':
    main()
