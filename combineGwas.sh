for gwas in {1..5}; do
	rm assoc.${gwas}.total
	for dir in {1..32}; do
		for chr in {1..22}; do
			cat dir${dir}/assoc.${chr}.${dir}.P${gwas}.assoc | grep -v "NA\|CHR" >> assoc.${gwas}.total
		done
	done
done
