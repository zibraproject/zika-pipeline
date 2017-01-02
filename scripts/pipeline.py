#!/usr/bin/env python
import argparse, csv, subprocess

def collect_strain_mapping(samples_dir):
    '''
    return dict
    each key is string "strain"
    each value is a list of tuples ("library", "barcode")
    '''
    runs_file = samples_dir + "runs.tsv"
    mapping = {}
    with open(runs_file) as tsv:
        for row in csv.DictReader(tsv, delimiter="\t"):
            sample = row["sample_id"]
            rb_pair = (row["run_name"], row["barcode_id"])
            if sample not in mapping:
                mapping[sample] = []
            mapping[sample].append(rb_pair)
    return mapping

def construct_sample_fastas(mapping, data_dir, build_dir):
    '''
    Use nanopolish to construct a single fasta for all reads from a sample
    '''
    for sample in mapping:
        print("* Extracting " + sample)
        # nanopolish extract each run/barcode pair
        for (run, barcode) in mapping[sample]:
            input_dir = data_dir + run + "/basecalled_reads/pass_demultiplex/" + barcode
            output_file = build_dir + sample + "_" + run + "_" + barcode + ".fasta"
            f = open(output_file, "w")
            call = map(str, ['nanopolish', 'extract', '--type', '2d', input_dir])
            print(" ".join(call) + " > " + output_file)
            subprocess.call(call, stdout=f)
        # concatenate to single sample fasta
        input_file_list = [build_dir + sample + "_" + run + "_" + barcode + ".fasta"
            for (run, barcode) in mapping[sample]]
        output_file = build_dir + sample + ".fasta"
        f = open(output_file, "w")
        call = map(str, ['cat'] + input_file_list)
        print(" ".join(call) + " > " + output_file)
        subprocess.call(call, stdout=f)
        print("")

def process_sample_fastas(mapping, data_dir, build_dir):
    '''
    Run fasta_to_consensus script to construct consensus files
    '''
    for sample in mapping:
        print("* Processing " + sample)
        sample_stem = build_dir + sample
        call = map(str, ['./zibra/zika-pipeline/scripts/fasta_to_consensus.sh', '/zibra/zika-pipeline/refs/KJ776791.2.fasta', sample_stem, '/zibra/zika-pipeline/metadata/v2_500.amplicons.ver2.bed'])
        print(" ".join(call))

if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "process data")
    parser.add_argument('--data_dir', type = str, default = "/data/")
    parser.add_argument('--samples_dir', type = str, default = "/samples/")
    parser.add_argument('--build_dir', type = str, default = "/build/")
    params = parser.parse_args()

    mapping = collect_strain_mapping(params.samples_dir)
    construct_sample_fastas(mapping, params.data_dir, params.build_dir)
    process_sample_fastas(mapping, params.data_dir, params.build_dir)
