
# convert all pdf files to png files, 600 dpi
for F in *.pdf
do 
	echo $F
	convert -density 600 $F -trim -bordercolor White -border 3x3 "${F/%.pdf/.png}"
	convert "${F/%.pdf/.png}" -strip -type palette "${F/%.pdf/.tiff}"
done

