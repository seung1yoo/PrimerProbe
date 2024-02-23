
from pathlib import Path
from itertools import product
from Bio import SeqIO
import gzip


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
                self.fastq_dic.setdefault(sample_name, {}).setdefault("r1", fastq_file)
            elif str(fastq_file).endswith("_2.fastq.gz"):
                self.fastq_dic.setdefault(sample_name, {}).setdefault("r2", fastq_file)
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



def main():
    data_path = Path("../20240206")
    #data_path = Path("../20240126_test")
    instance = Identifier()
    instance.load_data_path(data_path)
    #print(instance.fastq_dic)
    instance.add_primer_info("V3-V4", "CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC")
    instance.add_primer_info("V4", "GTGYCAGCMGCCGCGGTAA", "GACTACHVGGGTATCTAATCC")
    instance.add_primer_info("ITS1", "CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC")
    instance.add_primer_info("ITS2", "GCATCGATGAAGAACGCAGC", "TCCTCCGCTTATTGATATGC")
    #print(instance.primer_dic)
    outfn = "result_identifying_16S_primer.xls"
    outfh = open(outfn, "w")
    items = ["sample_name"]
    items.extend(["V3-V4","V4","ITS1","ITS2"])
    items.append("most_primer")
    outfh.write("{0}\n".format("\t".join(items)))
    for sample_name, path_dic in instance.fastq_dic.items():
        instance.primer_counting_dic = dict()
        instance.counting_primers(path_dic, 'r1', 1000)
        instance.counting_primers(path_dic, 'r2', 1000)
        items = [sample_name]
        items.append(instance.primer_counting_dic["V3-V4"])
        items.append(instance.primer_counting_dic["V4"])
        items.append(instance.primer_counting_dic["ITS1"])
        items.append(instance.primer_counting_dic["ITS2"])
        items.append(instance.who_am_i())
        outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
    outfh.close()






if __name__=='__main__':
    main()



