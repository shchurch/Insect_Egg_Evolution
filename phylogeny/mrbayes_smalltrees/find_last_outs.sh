# run this script within folders misof/oldouts/ or rainford/oldouts/
# it creates symlinks to the last output file for each taxon, so we can check
# final st dev of split frequencies

ls -1 *.gz | cut -d . -f 2 | sort | uniq | while read taxonid;
do
	lastid=$(ls -1 *.$taxonid.* | cut -d . -f 3 | sort | tail -n 1)
	mkdir -p last_outs
        cd last_outs
	ln -s ../*.$taxonid.$lastid.out.gz
	cd ..
done


