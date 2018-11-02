#get tumor and non tumor samples

touch tumorSamps
for f in $(ls *tcga-t.txt); do head -1 $f | tr "\t" "\n" | tail -n +3 >> tumorSamps; done

touch notumorSamps
for f in $(ls *tcga.txt); do head -1 $f | tr "\t" "\n" | tail -n +3 >> notumorSamps; done

