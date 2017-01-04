#!/usr/bin/env python
import argparse, csv, subprocess

def sample_to_run_data_mapping(samples_dir):
    '''
    return dict
    each key is string "sample_id"
    each value is a list of tuples ("library", "barcode")
    '''
    runs_file = samples_dir + "runs.tsv"
    sr_mapping = {}
    with open(runs_file) as tsv:
        for row in csv.DictReader(tsv, delimiter="\t"):
            sample = row["sample_id"]
            rb_pair = (row["run_name"], row["barcode_id"])
            if sample not in sr_mapping:
                sr_mapping[sample] = []
            sr_mapping[sample].append(rb_pair)
    return sr_mapping

def sample_to_metadata_mapping(samples_dir):
    '''
    return dict
    each key is string "sample_id"
    each value is a list of metadata ordered as
    ["strain", "sample_id", "collect_date", "country", "division", "location"]
    '''
    metadata_file = samples_dir + "samples.tsv"
    sm_mapping = {}
    with open(metadata_file) as tsv:
        for row in csv.DictReader(tsv, delimiter="\t"):
            sample = row["sample_id"]
            metadata = [row["strain"], row["sample_id"], row["collection_date"],
                row["country"], row["division"], row["location"]]
            sm_mapping[sample] = metadata
    return sm_mapping

def construct_sample_fastas(sr_mapping, data_dir, build_dir):
    '''
    Use nanopolish to construct a single fasta for all reads from a sample
    '''
    for sample in sr_mapping:
        print("* Extracting " + sample)
        # nanopolish extract each run/barcode pair
        for (run, barcode) in sr_mapping[sample]:
            input_dir = data_dir + run + "/basecalled_reads/pass_demultiplex/" + barcode
            output_file = build_dir + sample + "_" + run + "_" + barcode + ".fasta"
            f = open(output_file, "w")
            call = ['nanopolish', 'extract', '--type', '2d', input_dir]
            print(" ".join(call) + " > " + output_file)
            subprocess.call(call, stdout=f)
        # concatenate to single sample fasta
        input_file_list = [build_dir + sample + "_" + run + "_" + barcode + ".fasta"
            for (run, barcode) in sr_mapping[sample]]
        output_file = build_dir + sample + ".fasta"
        f = open(output_file, "w")
        call = ['cat'] + input_file_list
        print(" ".join(call) + " > " + output_file)
        subprocess.call(call, stdout=f)
        print("")

def process_sample_fastas(sm_mapping, build_dir):
    '''
    Run fasta_to_consensus script to construct consensus files
    '''
    for sample in sm_mapping:
        print("* Processing " + sample)
        # build consensus
        sample_stem = build_dir + sample
        call = ['/zibra/zika-pipeline/scripts/fasta_to_consensus.sh', '/zibra/zika-pipeline/refs/KJ776791.2.fasta', sample_stem, '/zibra/zika-pipeline/metadata/v2_500.amplicons.ver2.bed']
        print(" ".join(call))
        subprocess.call(call)
        # annotate consensus
        # >ZBRD116|ZBRD116|2015-08-28|brazil|alagoas|arapiraca
        fasta_header = ">" + "|".join(sm_mapping[sample])
        replacement = r"\~^>~s~.*~" + fasta_header + "~" # ~ rather than / to avoid conflict with strain names
        input_file = build_dir + sample + ".consensus.fasta"
        output_file = "temp.fasta"
        f = open(output_file, "w")
        call = ['sed', replacement, input_file]
        print(" ".join(call) + " > " + output_file)
        subprocess.call(call, stdout=f)
        call = ['mv', output_file, input_file]
        print(" ".join(call))
        subprocess.call(call)
        print("")

def gather_consensus_fastas(sm_mapping, build_dir, prefix):
    '''
    Gather consensus files into genomes with 'partial' (50-80% coverage)
    and good (>80% coverage) coverage
    '''
    # identify partial and good samples
    print("* Concatenating consensus fastas")
    partial_samples = []
    good_samples = []
    for sample in sm_mapping:
        consensus_file = build_dir + sample + ".consensus.fasta"
        with open(consensus_file) as f:
            lines = f.readlines()
        seq = lines[1]
        coverage = 1 - seq.count("N") / float(len(seq))
        if coverage >= 0.5 and coverage < 0.8:
            partial_samples.append(sample)
        if coverage >= 0.8:
            good_samples.append(sample)
    # sort samples
    partial_samples.sort()
    good_samples.sort()
    # concatenate partial samples
    print("Partial samples: " + " ".join(partial_samples))
    input_file_list = [build_dir + sample + ".consensus.fasta" for sample in partial_samples]
    output_file = build_dir + prefix + "_partial.fasta"
    f = open(output_file, "w")
    call = ['cat'] + input_file_list
    print(" ".join(call) + " > " + output_file)
    subprocess.call(call, stdout=f)
    # concatenate good samples
    print("Good samples: " + " ".join(good_samples))
    input_file_list = [build_dir + sample + ".consensus.fasta" for sample in good_samples]
    output_file = build_dir + prefix + "_good.fasta"
    f = open(output_file, "w")
    call = ['cat'] + input_file_list
    print(" ".join(call) + " > " + output_file)
    subprocess.call(call, stdout=f)
    print("")

if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "process data")
    parser.add_argument('--data_dir', type = str, default = "/data/")
    parser.add_argument('--samples_dir', type = str, default = "/samples/")
    parser.add_argument('--build_dir', type = str, default = "/build/")
    parser.add_argument('--prefix', type = str, default = "ZIKA")
    params = parser.parse_args()

    sr_mapping = sample_to_run_data_mapping(params.samples_dir)
    sm_mapping = sample_to_metadata_mapping(params.samples_dir)
    #construct_sample_fastas(sr_mapping, params.data_dir, params.build_dir)
    #process_sample_fastas(sm_mapping, params.build_dir)
    gather_consensus_fastas(sm_mapping, params.build_dir, params.prefix)
