#!/apps/python/3.12.3/bin/python

import gzip
import zlib
import gzip
import zlib
import re

VALID_SEQ = re.compile(r'^[ACGTNacgtnRYSWKMBDHVryswkmbdhv]+$')

def is_valid_record(header, seq, plus, qual):
    seq = seq.strip()
    qual = qual.strip()
    plus = plus.strip()
    if not header.startswith('@'):
        return False
    if not plus.startswith('+'):
        return False
    if len(seq) != len(qual):
        return False
    if not VALID_SEQ.match(seq):
        return False
    if not all(33 <= ord(c) <= 126 for c in qual):
        return False
    return True

def filter_fastq(infile, outfile):
    opener = gzip.open if infile.endswith('.gz') else open
    out_opener = gzip.open if outfile.endswith('.gz') else open

    skipped = 0
    kept = 0

    with opener(infile, 'rt') as fin, out_opener(outfile, 'wt') as fout:
        while True:
            try:
                header = fin.readline()
                if not header:
                    break
                seq = fin.readline()
                plus = fin.readline()
                qual = fin.readline()
            except (zlib.error, EOFError, OSError) as e:
                print(f"Decompression error after {kept} reads: {e} — stopping")
                break

            if not all([header, seq, plus, qual]):
                skipped += 1
                continue

            if is_valid_record(header, seq, plus, qual):
                fout.write(header)
                fout.write(seq)
                fout.write(plus)
                fout.write(qual)
                kept += 1
            else:
                skipped += 1

    print(f"{infile}: kept {kept}, skipped {skipped}")

# R1
filter_fastq('QC/B1201_M1_1_96_S24.R1.fq.gz', 'repaired_QC/B1201_M1_1_96_S24.R1.fq.gz')
filter_fastq('QC/B2490_F1_74_9_S2.R1.fq.gz', 'repaired_QC/B2490_F1_74_9_S2.R1.fq.gz')
filter_fastq('QC/B3015_F1_0_49_S28.R1.fq.gz', 'repaired_QC/B3015_F1_0_49_S28.R1.fq.gz')
filter_fastq('QC/B3648_M1_14_8_S19.R1.fq.gz', 'repaired_QC/B3648_M1_14_8_S19.R1.fq.gz')
filter_fastq('QC/B4092_F1_28_3_S15.R1.fq.gz', 'repaired_QC/B4092_F1_28_3_S15.R1.fq.gz')
filter_fastq('QC/B4276_F1_28_7_S20.R1.fq.gz', 'repaired_QC/B4276_F1_28_7_S20.R1.fq.gz')
filter_fastq('QC/B5053a_N1_1_14_S10.R1.fq.gz', 'repaired_QC/B5053a_N1_1_14_S10.R1.fq.gz')
filter_fastq('QC/B6565_F1_1_03_S13.R1.fq.gz', 'repaired_QC/B6565_F1_1_03_S13.R1.fq.gz')
filter_fastq('QC/B6658_L1_15_S25.R1.fq.gz', 'repaired_QC/B6658_L1_15_S25.R1.fq.gz')
filter_fastq('QC/B6852_F1_38_S1.R1.fq.gz', 'repaired_QC/B6852_F1_38_S1.R1.fq.gz')
filter_fastq('QC/B6894_F1_0_913_S23.R1.fq.gz', 'repaired_QC/B6894_F1_0_913_S23.R1.fq.gz')
filter_fastq('QC/B7002_1_F1_5_22_S9.R1.fq.gz', 'repaired_QC/B7002_1_F1_5_22_S9.R1.fq.gz')
filter_fastq('QC/B7017_M2_S22.R1.fq.gz', 'repaired_QC/B7017_M2_S22.R1.fq.gz')
filter_fastq('QC/B7024_M1_21_2_S11.R1.fq.gz', 'repaired_QC/B7024_M1_21_2_S1.R1.fq.gz')
filter_fastq('QC/B7106_F1_S21.R1.fq.gz', 'repaired_QC/B7106_F1_S21.R1.fq.gz')
filter_fastq('QC/B7115_N3_S40.R1.fq.gz', 'repaired_QC/B7115_N3_S40.R1.fq.gz')
filter_fastq('QC/B7115_N5_S41.R1.fq.gz', 'repaired_QC/B7115_N5_S41.R1.fq.gz')
filter_fastq('QC/B7116_N3_8_78_S8.R1.fq.gz', 'repaired_QC/B7116_N3_8_78_S8.R1.fq.gz')
filter_fastq('QC/B7119_F1_5_34_S5.R1.fq.gz', 'repaired_QC/B7119_F1_5_34_S5.R1.fq.gz')
filter_fastq('QC/B7203_N1_B3644_S26.R1.fq.gz', 'repaired_QC/B7203_N1_B3644_S26.R1.fq.gz')
filter_fastq('QC/B7217_F1_8_44_S22.R1.fq.gz', 'repaired_QC/B7217_F1_8_44_S22.R1.fq.gz')
filter_fastq('QC/B7229_F1_11_1_S17.R1.fq.gz', 'repaired_QC/B7229_F1_11_1_S17.R1.fq.gz')
filter_fastq('QC/B7380_F1_34_S21.R1.fq.gz', 'repaired_QC/B7380_F1_34_S21.R1.fq.gz')

#R2
filter_fastq('QC/B1201_M1_1_96_S24.R2.fq.gz', 'repaired_QC/B1201_M1_1_96_S24.R2.fq.gz')
filter_fastq('QC/B2490_F1_74_9_S2.R2.fq.gz', 'repaired_QC/B2490_F1_74_9_S2.R2.fq.gz')
filter_fastq('QC/B3015_F1_0_49_S28.R2.fq.gz', 'repaired_QC/B3015_F1_0_49_S28.R2.fq.gz')
filter_fastq('QC/B3648_M1_14_8_S19.R2.fq.gz', 'repaired_QC/B3648_M1_14_8_S19.R2.fq.gz')
filter_fastq('QC/B4092_F1_28_3_S15.R2.fq.gz', 'repaired_QC/B4092_F1_28_3_S15.R2.fq.gz')
filter_fastq('QC/B4276_F1_28_7_S20.R2.fq.gz', 'repaired_QC/B4276_F1_28_7_S20.R2.fq.gz')
filter_fastq('QC/B5053a_N1_1_14_S10.R2.fq.gz', 'repaired_QC/B5053a_N1_1_14_S10.R2.fq.gz')
filter_fastq('QC/B6565_F1_1_03_S13.R2.fq.gz', 'repaired_QC/B6565_F1_1_03_S13.R2.fq.gz')
filter_fastq('QC/B6658_L1_15_S25.R2.fq.gz', 'repaired_QC/B6658_L1_15_S25.R2.fq.gz')
filter_fastq('QC/B6852_F1_38_S1.R2.fq.gz', 'repaired_QC/B6852_F1_38_S1.R2.fq.gz')
filter_fastq('QC/B6894_F1_0_913_S23.R2.fq.gz', 'repaired_QC/B6894_F1_0_913_S23.R2.fq.gz')
filter_fastq('QC/B7002_1_F1_5_22_S9.R2.fq.gz', 'repaired_QC/B7002_1_F1_5_22_S9.R2.fq.gz')
filter_fastq('QC/B7017_M2_S22.R2.fq.gz', 'repaired_QC/B7017_M2_S22.R2.fq.gz')
filter_fastq('QC/B7024_M1_21_2_S11.R2.fq.gz', 'repaired_QC/B7024_M1_21_2_S1.R2.fq.gz')
filter_fastq('QC/B7106_F1_S21.R2.fq.gz', 'repaired_QC/B7106_F1_S21.R2.fq.gz')
filter_fastq('QC/B7115_N3_S40.R2.fq.gz', 'repaired_QC/B7115_N3_S40.R2.fq.gz')
filter_fastq('QC/B7115_N5_S41.R2.fq.gz', 'repaired_QC/B7115_N5_S41.R2.fq.gz')
filter_fastq('QC/B7116_N3_8_78_S8.R2.fq.gz', 'repaired_QC/B7116_N3_8_78_S8.R2.fq.gz')
filter_fastq('QC/B7119_F1_5_34_S5.R2.fq.gz', 'repaired_QC/B7119_F1_5_34_S5.R2.fq.gz')
filter_fastq('QC/B7203_N1_B3644_S26.R2.fq.gz', 'repaired_QC/B7203_N1_B3644_S26.R2.fq.gz')
filter_fastq('QC/B7217_F1_8_44_S22.R2.fq.gz', 'repaired_QC/B7217_F1_8_44_S22.R2.fq.gz')
filter_fastq('QC/B7229_F1_11_1_S17.R2.fq.gz', 'repaired_QC/B7229_F1_11_1_S17.R2.fq.gz')
filter_fastq('QC/B7380_F1_34_S21.R2.fq.gz', 'repaired_QC/B7380_F1_34_S21.R2.fq.gz')
