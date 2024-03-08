
from pathlib import Path
from itertools import product
from Bio import SeqIO
import gzip
import csv
import logging
logging.basicConfig(level=logging.INFO)


class Identifier:
    def __init__(self):
        self.fastq_dic = dict()
        self.primer_dic = dict()

    def add_primer_info(self, primer_type, primer_for, primer_rev):
        possible_sequences = self.generate_degenerate_sequences(primer_for)
        for sequence in possible_sequences:
            self.primer_dic.setdefault(primer_type, {}).setdefault("r1", {}).setdefault(sequence, 0)
        possible_sequences = self.generate_degenerate_sequences(primer_rev)
        for sequence in possible_sequences:
            self.primer_dic.setdefault(primer_type, {}).setdefault("r2", {}).setdefault(sequence, 0)
        return 1

    def generate_degenerate_sequences(self, degenerate_primer):
        # R = G or A
        # Y = C or T
        # M = A or C
        # K = G or T
        # S = G or C
        # W = A or T
        # B = not A (C or G or T)
        # D = not C (A or G or T)
        # H = not G (A or C or T)
        # V = not T (A or C or G)
        # N = A or C or G or T
        sequences = []
        for base in degenerate_primer:
            if base == 'N':
                sequences.append(['A', 'T', 'C', 'G'])
            elif base == 'W':
                sequences.append(['A', 'T'])
            elif base == 'Y':
                sequences.append(['C', 'T'])
            elif base == 'M':
                sequences.append(['A', 'C'])
            elif base == 'H':
                sequences.append(['A', 'C', 'T'])
            elif base == 'V':
                sequences.append(['A', 'C', 'G'])
            else:
                sequences.append([base])

        possible_combinations = product(*sequences)
        degenerate_sequences = [''.join(combination) for combination in possible_combinations]
        return degenerate_sequences

    def load_data_path(self, data_path):
        for fastq_file in data_path.glob("*.fastq.gz"):
            sample_name = fastq_file.name.split("_")[0]
            if str(fastq_file).endswith("_1.fastq.gz"):
                self.fastq_dic.setdefault(sample_name, {}).setdefault("r1", fastq_file.resolve())
            elif str(fastq_file).endswith("_2.fastq.gz"):
                self.fastq_dic.setdefault(sample_name, {}).setdefault("r2", fastq_file.resolve())
            else:
                self.fastq_dic.setdefault(sample_name, {}).setdefault("unknown", fastq_file)
        return 1

    def counting_primers(self, path_dic, read_type, read_limit):
        read_count = 0
        with gzip.open(str(path_dic[read_type]), "rt") as filehandle:
            for record in SeqIO.parse(filehandle, "fastq"):
                # counting primers
                for primer_type, info_dic in self.primer_dic.items():
                    self.primer_counting_dic.setdefault(primer_type, 0)
                    for primer_seq in info_dic[read_type]:
                        #if record.seq.startswith(primer_seq):
                        if primer_seq in record.seq:
                            self.primer_counting_dic[primer_type] += 1
                # check the limit
                read_count += 1
                if read_count in [read_limit]:
                    break
            # end biopython
        return 1

    def who_am_i(self):
        max_value = max(self.primer_counting_dic.values())
        if max_value in [0]:
            return "unknown"
        keys_with_max_value = [key for key, value in self.primer_counting_dic.items() if value == max_value]
        return ';'.join(keys_with_max_value)

    def parse_fastp_table(self, fn):
        if not Path(fn).exists():
            logging.warning(f"{fn} file is not exists.")
            return 0

        for line in open(fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ["Samples"]:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            sample_name = items[idx_dic["Samples"]]
            read_count = items[idx_dic["ReadsCount"]]
            good_rate = items[idx_dic["GoodReadRate"]]
            r1_mean_len = items[idx_dic["R1MeanLen"]]
            r2_mean_len = items[idx_dic["R2MeanLen"]]
            self.fastq_info_dic.setdefault(sample_name, {}).setdefault("read_count", read_count)
            self.fastq_info_dic.setdefault(sample_name, {}).setdefault("good_rate", good_rate)
            self.fastq_info_dic.setdefault(sample_name, {}).setdefault("r1_mean_len", r1_mean_len)
            self.fastq_info_dic.setdefault(sample_name, {}).setdefault("r2_mean_len", r2_mean_len)
        return 1

    def parse_meta_table(self, fn):
        for line in open(fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ["sampleID"]:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            sample_name = items[idx_dic["sampleID"]]
            project_id = items[idx_dic["projectID"]]
            race = items[idx_dic["Race"]]
            if race in ["uncalculated", ""]:
                race = "NA"
            continent = items[idx_dic["Continent"]]
            if continent in ["uncalculated", ""]:
                continent = "NA"
            self.meta_info_dic.setdefault(sample_name, {}).setdefault("project_id", project_id)
            self.meta_info_dic.setdefault(sample_name, {}).setdefault("race", race)
            self.meta_info_dic.setdefault(sample_name, {}).setdefault("continent", continent)
        return 1


def main():
    instance = Identifier()
    instance.load_data_path(Path("../890_Koreans"))
    instance.load_data_path(Path("../1_4464_foreigners"))
    instance.add_primer_info("V3-V4", "CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC")
    instance.add_primer_info("V4", "GTGCCAGCMGCCGCGGTAA", "GGACTACHVGGGTWTCTAAT")
    instance.add_primer_info("V1-V3", "AGAGTTTGATCMTGGCTCAG", "CTGCTGCCTYCCGTA")
    instance.add_primer_info("V4-V5", "GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTHTCTAAT")
    instance.add_primer_info("V5-V6", "GYACWCACCWGAGCTG", "CCGTCACCTTGTTACGACTT")
    instance.add_primer_info("V6-V8", "AACGCGAAGAACCTTAC", "CGGTGTGTACAAGACCC")
    instance.add_primer_info("ITS1", "CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC")
    instance.add_primer_info("ITS2", "GCATCGATGAAGAACGCAGC", "TCCTCCGCTTATTGATATGC")
    primer_set_s = ["V3-V4","V4","V1-V3","V4-V5","V5-V6","V6-V8","ITS1","ITS2"]

    instance.fastq_info_dic = dict()
    instance.parse_fastp_table("../890_Koreans/fastp.summary.tsv") # optional
    instance.parse_fastp_table("../1_4464_foreigners/fastp.summary.tsv") # optional

    instance.meta_info_dic = dict()
    instance.parse_meta_table("metadata.korea890.tsv") # optional
    instance.parse_meta_table("metadata_global4464.tsv") # optional

    outfn = "result_identifying_16S_primer.xls"
    if Path(outfn).exists():
        logging.info(f"{outfn} is exists. PASS the task.")
    else:
        logging.info(f"write the result of 16S primer identification into {outfn}")
        outfh = open(outfn, "w")
        items = ["sample_name"]
        items.extend(primer_set_s)
        items.append("most_primer")
        items.append("read_count") # optional
        items.append("good_rate") # optional
        items.append("mean_len") # optional
        items.append("project_id") # optional
        items.append("race") # optional
        items.append("continent") # optional
        outfh.write("{0}\n".format("\t".join(items)))
        for sample_name, path_dic in instance.fastq_dic.items():
            instance.primer_counting_dic = dict()
            instance.counting_primers(path_dic, 'r1', 1000)
            instance.counting_primers(path_dic, 'r2', 1000)
            items = [sample_name]
            for primer_set in primer_set_s:
                items.append(instance.primer_counting_dic[primer_set])
            items.append(instance.who_am_i())
            if sample_name in instance.fastq_info_dic:
                items.append(instance.fastq_info_dic[sample_name]['read_count'])
                items.append(instance.fastq_info_dic[sample_name]['good_rate'])
                items.append(instance.fastq_info_dic[sample_name]['r1_mean_len'])
            else:
                items.append("NA")
                items.append("NA")
                items.append("NA")
                logging.warning(f"{sample_name} has not fastp info.")
            if sample_name in instance.meta_info_dic:
                items.append(instance.meta_info_dic[sample_name]['project_id'])
                items.append(instance.meta_info_dic[sample_name]['race'])
                items.append(instance.meta_info_dic[sample_name]['continent'])
            else:
                items.append("NA")
                items.append("NA")
                items.append("NA")
                logging.warning(f"{sample_name} has not metadata info.")
            outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()


    # optional
    unusable_files = ['unusable_samples.not_v3v4.v240223']
    unusable_files.append("unusable_samples.low_read_count.v240227")
    unusable_files.append("unusable_samples.low_read_count.v240308")
    unusable_files.append("unusable_samples.global4464_race_NA.v240308")
    unusable_files.append("unusable_samples.global4464_primer_filter.v240308")
    unusable_samples = list()
    for fn in unusable_files:
        for line in open(fn):
            if line.startswith("#"):
                continue
            items = line.rstrip('\n').split('\t')
            if items[0] not in unusable_samples:
                unusable_samples.append(items[0])
        #
    logging.info(f"un-usable samples (n) : {len(unusable_samples)}")

    outfn = "samplesheet.csv"
    logging.info(f"write a samplesheet file to run nf-core/ampliseq into {outfn}")
    outfh = open(outfn, "w")
    headers = ["sampleID","forwardReads","reverseReads","run"]
    writer = csv.DictWriter(outfh, fieldnames=headers)
    writer.writeheader()
    for sample_name in sorted(instance.fastq_dic):
        if sample_name in unusable_samples:
            continue
        row_dic = dict()
        row_dic.setdefault("sampleID", sample_name)
        row_dic.setdefault("forwardReads", instance.fastq_dic[sample_name]["r1"])
        row_dic.setdefault("reverseReads", instance.fastq_dic[sample_name]["r2"])
        row_dic.setdefault("run", "1")
        writer.writerow(row_dic)
    outfh.close()

    outfn = "metadata.tsv"
    logging.info(f"write a metadata file to run nf-core/ampliseq into {outfn}")
    outfh = open(outfn, "w")
    headers = ["sampleID","projectID","Race","Continent"]
    outfh.write("{0}\n".format("\t".join(headers)))
    for sample_name in sorted(instance.fastq_dic):
        if sample_name in unusable_samples:
            continue
        items = [sample_name]
        items.append(instance.meta_info_dic[sample_name]['project_id'])
        items.append(instance.meta_info_dic[sample_name]['race'])
        items.append(instance.meta_info_dic[sample_name]['continent'])
        outfh.write("{0}\n".format("\t".join(items)))
    outfh.close()








if __name__=='__main__':
    main()



