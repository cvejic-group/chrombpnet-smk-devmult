import gzip
import h5py
import numpy as np

BED_PKS = snakemake.input['bed_pk']
H5_SCORES = snakemake.input['h5']
BEDPE_OUT = snakemake.output['bedpe']
WINDOW_SIZE = snakemake.params['width']

# Load peaks
peak_rows = None
with open(BED_PKS, 'r') as peaks_file:
    peak_rows = peaks_file.read().splitlines()

# open BEDPE output
with gzip.open(BEDPE_OUT, 'wt') as f_out:

    # load H5 modisco
    with h5py.File(H5_SCORES, "r") as grp:

        # iterate over positive then negative groups
        for contribution_dir in ['pos', 'neg']:
            patterns_category = f'{contribution_dir}_patterns'
            if patterns_category not in grp:
              continue

            # iterate over patterns
            for (pattern_name, datasets) in grp[patterns_category].items():

                # iterate over seqlets
                for idx in range(datasets['seqlets']['start'].shape[0]):
                    seqlet_name = f'{contribution_dir}_patterns.{pattern_name}.{idx}'

                    # identify target peak
                    target_idx = datasets['seqlets']['example_idx'][idx]
                    peak_row = peak_rows[target_idx].split('\t')
                    
                    # chrom(s)
                    chrom = peak_row[0]
    
                    # positions
                    target_start = int(peak_row[1])
                    target_end = int(peak_row[2])
                    
                    target_center = (target_end + target_start) // 2
                    window_center_offset = WINDOW_SIZE // 2
                    seqlet_start_offset = datasets['seqlets']['start'][idx] + 1
                    seqlet_end_offset = datasets['seqlets']['end'][idx]
                    seqlet_start = target_center - window_center_offset + seqlet_start_offset
                    seqlet_end = target_center - window_center_offset + seqlet_end_offset
    
                    # scores
                    target_score = peak_row[4]
                    contrib_score = datasets['seqlets']['contrib_scores'][idx]
                    seqlet_score = np.sum(contrib_score**2)
    
                    # strand
                    target_strand = '*' if peak_row[5] == '.' else peak_row[5]
                    seqlet_strand = '-' if bool(datasets['seqlets']['is_revcomp'][idx]) is True else '+'
    
                    row_vars = [chrom, seqlet_start, seqlet_end, 
                                chrom, target_start, target_end, 
                                seqlet_name, seqlet_score, 
                                seqlet_strand, target_strand, 
                                target_score]
                    f_out.write("\t".join(map(str, row_vars)) + "\n")
