#!/usr/bin/env python3

import argparse
import os
import re

# Renames TrAEL index-demultiplexed fastq/QC/alignment files to include the sample names given in
# the sample sheet's Names column, matched up against the Index column (e.g. Index "5,6,7" with
# Names "Alg70bp1, 70 in ESAG8, TbCas9" means index5 -> Alg70bp1, index6 -> 70 in ESAG8, etc).
#
# This follows the same file tag convention as the existing site tool
# /bi/apps/TrAELseq/latest/TrAEL-seq/rename_trael_samples.py (same TrAEL barcode sequences, same
# renamed-tag format, same sample-name cleaning), but that tool renames every matching file in the
# current directory regardless of which original sample it came from -- fine when a folder only
# ever holds one TrAEL sample, but wrong when several samples share a folder (they can reuse the
# same TrAEL index numbers with different meanings each). This version scopes every rename to
# files that also contain the given sample's own Illumina barcode, so samples sharing a folder
# can't collide.
#
# The same renaming is used for WGS-with-indexing groups (nf_WGS_4bpIndex / nf_WGS_4bpIndexPE),
# whose files are tagged differently: WGS_index_preprocess(_PE).py inserts just the 4bp index
# sequence straight after the lane, e.g. '..._L001_AGTC_R1.fastq.gz', with no 'index<N>' part.
# --wgs_index switches to matching that tag instead (anchored on the lane so the 4bp sequence
# can't match by chance elsewhere in the filename), e.g. '_L001_AGTC_' -> '_L001_AGTC_rnh_'.

INDEX_BARCODES = {
    "1": "AGTC",
    "2": "GACT",
    "3": "CTTG",
    "4": "TCGA",
    "5": "AAGG",
    "6": "TTCC",
    "7": "GTGC",
    "8": "GCCA",
    "9": "GATG",
}


def clean_sample_name(name):
    # Same cleaning rule as rename_trael_samples.py, so names end up formatted consistently
    # regardless of which of the two tools was used.
    return re.sub(r'[^a-zA-Z0-9.\-_]+_?', '_', name.strip())


def rename_for_sample(directory, barcode, index_digits, sample_names, wgs_index=False):
    indexes = list(index_digits.strip())
    names = [n.strip() for n in sample_names.strip().split(',')]

    if len(indexes) != len(names):
        print(f"!! number of indexes ({len(indexes)}) and sample names ({len(names)}) do not "
              f"match for barcode {barcode} -- skipping renames for this sample")
        return

    barcode_upper = barcode.strip().upper()

    for index, name in zip(indexes, names):
        if index not in INDEX_BARCODES:
            print(f"!! unknown TrAEL index '{index}' -- skipping")
            continue

        if wgs_index:
            tag = rf"_L\d{{3}}_{INDEX_BARCODES[index]}"
            old_tag = f"_L<lane>_{INDEX_BARCODES[index]}_"
        else:
            tag = re.escape(f"_{INDEX_BARCODES[index]}_index{index}")
            old_tag = f"_{INDEX_BARCODES[index]}_index{index}_"
        clean_name = clean_sample_name(name)
        tag_re = re.compile(f"({tag})_")
        already_renamed_re = re.compile(f"{tag}_{re.escape(clean_name)}_")

        matched_any = False
        for fname in sorted(os.listdir(directory)):
            if barcode_upper not in fname.upper():
                continue
            if not tag_re.search(fname):
                continue
            matched_any = True
            # Already-renamed files still contain the old tag, so skip them rather than adding
            # the sample name a second time on a re-run.
            if already_renamed_re.search(fname):
                continue
            new_fname = tag_re.sub(lambda m: f"{m.group(1)}_{clean_name}_", fname, count=1)
            old_path = os.path.join(directory, fname)
            new_path = os.path.join(directory, new_fname)
            if os.path.exists(new_path):
                print(f"Already renamed: {new_fname}")
                continue
            os.rename(old_path, new_path)
            print(f"Renamed {fname} -> {new_fname}")

        if not matched_any:
            print(f"!! no files found for barcode {barcode} / index{index} "
                  f"(looked for '{barcode_upper}' and '{old_tag}')")


def main():
    parser = argparse.ArgumentParser(
        description="Rename TrAEL index-demultiplexed files with sample names from the sample "
                    "sheet, scoped to one sample's own Illumina barcode so multiple samples "
                    "sharing a folder don't collide.")
    parser.add_argument('--directory', default='.',
                         help="Directory containing the files to rename [default: current directory]")
    parser.add_argument('--barcode', required=True,
                         help="This sample's Illumina barcode (Barcode 1 from the sample sheet)")
    parser.add_argument('--index', required=True,
                         help='TrAEL index digits for this sample, e.g. "567"')
    parser.add_argument('--sample_names', required=True,
                         help='Comma-separated names in the same order as --index, '
                              'e.g. "Alg70bp1, 70 in ESAG8, TbCas9"')
    parser.add_argument('--wgs_index', action='store_true',
                         help="Files come from nf_WGS_4bpIndex(PE) ('_L001_AGTC_' tags) rather than "
                              "TrAEL ('_AGTC_index1_' tags)")
    args = parser.parse_args()

    rename_for_sample(args.directory, args.barcode, args.index, args.sample_names, args.wgs_index)


if __name__ == "__main__":
    main()
