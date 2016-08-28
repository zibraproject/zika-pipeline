
all: $(PREFIX)-bigtree.pdf $(PREFIX)-bigtree.png

$(PREFIX)_aligned_trim_short_raxml.nwk $(PREFIX)_microreact_public.csv:
	go_micro.sh $(FASTA) $(PREFIX) $(THREADS)

$(PREFIX)_dist_$(DIST).clusters.txt: $(PREFIX)_aligned_trim_short_raxml.nwk $(PREFIX)_microreact_public.csv
	go_tree.sh $(PREFIX) $(PREFIX)_dist_$(DIST) $(DIST)

$(PREFIX)_dist_0.001.clusters.txt: $(PREFIX)_aligned_trim_short_raxml.nwk $(PREFIX)_microreact_public.csv
	go_tree.sh $(PREFIX) $(PREFIX)_dist_0.001 0.001

$(PREFIX)-bigtree.pdf: $(PREFIX)_dist_$(DIST).clusters.txt $(PREFIX)_dist_0.001.clusters.txt
	xvfb-run pdf_tree.py --legend --clusters $(PREFIX)_dist_0.001_clusters.txt --subclusters $(PREFIX)_dist_$(DIST)_clusters.txt $(PREFIX)_aligned_trim_short_raxml.nwk $(PREFIX)_microreact_private.csv $@

$(PREFIX)-bigtree.png: $(PREFIX)_dist_$(DIST).clusters.txt $(PREFIX)_dist_0.001.clusters.txt
	sleep 5;
	xvfb-run pdf_tree.py --legend --clusters $(PREFIX)_dist_0.001_clusters.txt --subclusters $(PREFIX)_dist_$(DIST)_clusters.txt $(PREFIX)_aligned_trim_short_raxml.nwk $(PREFIX)_microreact_private.csv $@

#find . -name "oct*.pdf" | xargs -L 1 -I '{}' convert -verbose -resize 3000 -density 1200 -trim {} -quality 100 -sharpen 0x1.0 {}.png
#pandoc -f markdown -t html report-20151004.md > report-20151004.html
#sed -i tmp 's/src/width=800 src/' report-20151005.html
sed --in-place 's/src/width=800 src//' report-20151004.html

