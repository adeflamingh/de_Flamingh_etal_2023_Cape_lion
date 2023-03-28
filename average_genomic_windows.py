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
    p.add_argument('-f', '--fai', required=True, help='Path to the genome index (FAI).')
    p.add_argument('-v', '--var-sites', required=False, help='Path to MAF file with variant site positions.')
    p.add_argument('-d', '--depth', required=False, help='Path to SAMtools depth file.')
    p.add_argument('-o', '--outdir', required=False, default='.', help='Path to output directory.')
    p.add_argument('-m', '--min-length', required=False, default=1_000_000, type=float, help='Minimum chromosome size in bp [default 1000000 (1Mbp)]')
    p.add_argument('-s', '--window-size', required=False, default=250_000, type=float, help='Size of windows in bp [default 250000 (250Kbp)].')
    p.add_argument('-t', '--window-step', required=False, default=100_000, type=float, help='Step of windows in bp [default 100000 (100Kbp)].')
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
def calculate_chr_window_intervals(chr_len, window_size=250_000, window_step=50_000):
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
def populate_variant_site_windows(variant_site, chrom_id, chr_window_boundaries, populated_window_intervals):
    assert len(chr_window_boundaries) > 0
    assert type(populated_window_intervals[chrom_id]) is dict
    # Loop over the windows and check the sites
    for window in chr_window_boundaries:
        assert len(window) == 2
        window_sta = window[0]
        window_end = window[1]
        # Check the site:
        # Skip windows that lie before the site
        if window_end < variant_site:
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
    return populated_window_intervals

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
    print(f'\nGenerated window intervals for {len(genome_window_intervals):,} chromosomes/scaffolds.')
    return genome_window_intervals

#
# Parse the variant sites file and populate windows
#
def read_variant_sites(variant_sites_f, genome_window_intervals):
    print('\nTallying variant sites per chromosome...')
    populated_window_intervals = dict()
    chr_position_tally = dict()
    fh = None
    if variant_sites_f.endswith('gz'):
        fh = gzip.open(variant_sites_f, 'rt')
    else:
        fh = open(variant_sites_f, 'r')
    for line in fh:
        if line.startswith('#'):
            continue
        if line.startswith('chromo'):
            continue
        line = line.strip('\n')
        fields = line.split('\t')
        if not len(fields) >= 2:
            sys.exit('Error: variant sites file must have >2 columns, scaffold_id<tab>position_bp<tab>major<tab>...')
        chr_id = fields[0]
        if not fields[1].isnumeric():
            sys.exit('Error: column two of the variant sites file must be an integer position.')
        position = int(fields[1])
        # Skip positions from unwanted chromosomes
        if chr_id not in genome_window_intervals:
            continue
        # Tally the position from the current chromosome
        chr_position_tally.setdefault(chr_id, 0)
        chr_position_tally[chr_id] += 1
        # Populate the windows
        chr_window_boundaries = genome_window_intervals.get(chr_id, None)
        populated_window_intervals.setdefault(chr_id, dict())
        populated_window_intervals = populate_variant_site_windows(position, chr_id, chr_window_boundaries, populated_window_intervals)
    genome_total = 0
    for chrom in chr_position_tally:
        n_sites = chr_position_tally[chrom]
        genome_total += n_sites
        print(f'    {chrom} : {n_sites:,} variant sites found.')
    print(f'\n    Genome-wide : {genome_total:,} total variant sites found.')
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
def populate_depth_windows(chr_id, position, depth, chr_window_boundaries, populated_depth_windows):
    assert len(chr_window_boundaries) > 0
    assert type(populated_depth_windows[chr_id]) is dict
    assert depth >= 0
    # Loop over the windows and check the position
    for window in chr_window_boundaries:
        assert len(window) == 2
        start = window[0]
        end = window[1]
        # Skip windows that end before the position of interest
        if end < position:
            continue
        # Position is in the window
        elif start <= position < end:
            # Window's middle position, serves as reference
            mid = int(start + ((end - start)/2))
            # Initialize the dictionary if empty
            populated_depth_windows[chr_id].setdefault(mid, list())
            # Append the depth to the window list
            populated_depth_windows[chr_id][mid].append(depth)
        # Break when moving past the position
        elif start > position:
            break
    return populated_depth_windows

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
            depth_list = chr_windows[pos]
            mean_cov = statistics.mean(depth_list)
            non_zero_n = len([ cov for cov in depth_list if cov > 0 ])
            row = f'{chr}\t{pos}\t{non_zero_n}\t{mean_cov:.08g}\n'
            outfh.write(row)
    outfh.close()

#
# Parse the SAMtools depth file
#
def parse_depth_file(depth_f, genome_window_intervals):
    print('\nTallying depth per chromosome...')
    populated_depth_windows = dict()
    chr_depth_tally = dict()
    fh = None
    if depth_f.endswith('gz'):
        fh = gzip.open(depth_f, 'rt')
    else:
        fh = open(depth_f, 'r')
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
        # Tally the current chromosome
        chr_depth_tally.setdefault(chr_id, [0, 0])
        chr_depth_tally[chr_id][0] += 1
        if depth > 0:
            chr_depth_tally[chr_id][1] += 1
        # Populate the windows
        chr_window_boundaries = genome_window_intervals.get(chr_id, None)
        populated_depth_windows.setdefault(chr_id, dict())
        populated_depth_windows = populate_depth_windows(chr_id, position, depth, chr_window_boundaries, populated_depth_windows)
    genome_total = 0
    genome_non_zero = 0
    for chrom in chr_depth_tally:
        # Total sites
        total_sites = chr_depth_tally[chrom][0]
        genome_total += total_sites
        # Non-zero sites
        non_zero_sites = chr_depth_tally[chrom][1]
        non_zero_freq = non_zero_sites/total_sites
        genome_non_zero += non_zero_sites
        print(f'    {chrom} :  {total_sites:,} sites seen, {non_zero_sites:,} ({non_zero_freq:.02%}) with non-zero depth.')
    genome_non_zero_f = genome_non_zero/genome_total
    print(f'\n    Genome-wide : {genome_total:,} sites seen, {genome_non_zero:,} ({genome_non_zero_f:.02%}) with non-zero depth.')
    return populated_depth_windows

def main():
    args = parse_args()
    # Initialize script
    print(f'Started {PROG} on {date_now()} {time_now()}.')
    print(f'    Min Chrom Size (bp): {args.min_length:,}')
    print(f'    Window Size (bp): {args.window_size:,}')
    print(f'    Window Step (bp): {args.window_step:,}')
    # Get windows from the fai
    genome_window_intervals = set_windows_from_fai(args.fai, args.window_size, args.window_step, args.min_length)
    # Parse the variant sites file and print output
    if args.var_sites is not None:
        populated_window_sites = read_variant_sites(args.var_sites, genome_window_intervals)
        print_variant_site_tally(populated_window_sites, args.window_size, args.outdir)
    # Parse the depth file and print output
    if args.depth is not None:
        populated_depth_windows = parse_depth_file(args.depth, genome_window_intervals)
        print_average_coverage_file(populated_depth_windows, args.window_size, args.outdir)

    print(f'\nFinished {PROG} on {date_now()} {time_now()}.')

# Run Code
if __name__ == '__main__':
    main()
