SHELL=/bin/bash -o pipefail

all: $(PREFIX)-$(TPR).fasta.tagged

%.fasta.tagged: %.fasta
	tagfastas.py $(METADATA) < $(PREFIX)-tpr$(TPR).fasta > $(PREFIX)-tpr$(TPR).tagged.fasta

%.fasta: $(METADATA)
	make_stats_file.py $(METADATA) $(SET)  | intersection_vcf_stats.py /dev/stdin | awk '($$10>=$(TPR)&&$$10<=$(TPRMAX))' | grep np-new-filter_qual200-50 | cut -f1 | sort | uniq | xargs -L 1 -I '{}' margin_cons.py ../refs/EM_079517.fasta {}_hq_EM_079517_np_primer.tagged.vcf EM_079517_{}_hq_marginalign.sorted.bam all >$(PREFIX)-tpr$(TPR).fasta 2>$(PREFIX)-tpr$(TPR).stderr

