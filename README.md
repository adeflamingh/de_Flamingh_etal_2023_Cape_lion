# Cape_Lion
Scripts and code associated with Cape lion project

## Average values across genomic windows

Python3 (version 3.6 or higher) script to generate an average, per-window coverage from a *SAMtools* `depth` output file and/or generate a tally/frequency of kept variant sites per genomic windows from a `*.maf` file. Requires an index of the reference genome in `fai` format (see *SAMtools* `faidx`). The size (total span in bp) and step (how much bps the window moves forward in overlapping windows) are specified by `--window-size` and `--window-step`, respectively. The target chromosomes can be filtered using `--min-length`, or by prior filtering of the `fai` index file.

```sh
usage: average_genomic_windows.py [-h] -f FAI [-v VAR_SITES] [-d DEPTH]
                                  [-o OUTDIR] [-m MIN_LENGTH] [-s WINDOW_SIZE]
                                  [-t WINDOW_STEP]

optional arguments:
  -h, --help            show this help message and exit
  -f FAI, --fai FAI     Path to the genome index (FAI).
  -v VAR_SITES, --var-sites VAR_SITES
                        Path to file with variant site positions.
  -d DEPTH, --depth DEPTH
                        Path to SAMtools depth file.
  -o OUTDIR, --outdir OUTDIR
                        Path to output directory.
  -m MIN_LENGTH, --min-length MIN_LENGTH
                        Minimum chromosome size in bp [default 1000000 (1Mbp)]
  -s WINDOW_SIZE, --window-size WINDOW_SIZE
                        Size of windows in bp [default 250000 (250Kbp)].
  -t WINDOW_STEP, --window-step WINDOW_STEP
                        Step of windows in bp [default 100000 (100Kbp)].
```

