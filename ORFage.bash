#!/bin/bash 

#### Conservation of Unannotated Gene Sequences 

# Zachary Ardern, 2017
# Contact: zachary.ardern@tum.de 
# Script Updates at https://github.com/ZacharyArdern/ORFage 


# INSTRUCTIONS: 

########### EXAMPLE: bash ORFage.sh example.bed 20 
#
########### Run on one input file at a time. 
#
########### Set up for E coli genomes, SSUNR99 SILVA 128 tree - otherwise change line 65 (origin= )
########### See note at line 216 for intergenic sequences   

#
####### INPUT ARGUMENT 1: BED file with gene positions for one chromosome 
# # e.g. first row: 	NC_002695.1	385213	385980	plus  
#	Note, BED files must have 0-based co-ordinates 
#
####### INPUT ARGUMENT 2: Number of CPU (cores) to Use 
#  
#
### Inputs 3 & 4 are obtained from SILVA database, e.g.
### https://www.arb-silva.de/fileadmin/arb_web_db/release_128/ARB_files/SSURef_NR99_128_SILVA_07_09_16_opt.arb.gz  
### and exported as trees using 'arb' program
#
####### INPUT ARGUMENT 3: Tree file in Newick format, with full names with species
#
#
####### INPUT ARGUMENT 4: Tree file in Newick format, with short names  
#
#
# OUTPUT: Summary Table, and png image with gnuplot 
#
####### PROGRAMS REQUIRED (added to PATH): 
#faidx
#pyfasta
#NCBI Blast+
#newick utils
#needleall  


#######
### Set Up 
#######

echo "Conservation of Unannotated Gene Sequences - script running" ;

input="$1" ;
bed=${input%.*} ;
cores="$2" ;
tree1="$3" ;
tree2="$4" ;  

echo -e "bed file:" "$input" "\n" "cores:" "$cores" ;

# set chromosome
chr=$(cat "$input" | awk '{if (NR==1) print $1}' ) ;

# set SILVA entry for original genome 
origin=BA000007


test ! -e LTPs128_SSU_tree.newick && wget https://www.arb-silva.de/fileadmin/silva_databases/living_tree/LTP_release_128/LTPs128_SSU/LTPs128_SSU_tree.newick ;


mkdir tmp ;
mkdir fastas ;

###############################################################################################################################

#######
####Create Input Fastas
#######


#esearch -db nucleotide -query "$chr" < /dev/null | 
#/home/zachary/edirect/efetch -format fasta > $chr.fna ; 



cat $input | awk '{if ($4=="+") print $1 "\t" $2 "\t" $3 "\t" $4 }' | faidx $chr.fna -b - | transeq -filter > "$bed"_plus.fasta ; 
cat $input | awk '{if ($4=="-") print $1 "\t" $2 "\t" $3 "\t" $4 }' | faidx $chr.fna -rcb - | transeq -filter > "$bed"_minus.fasta ; 


n=$(echo "($cores-3 ) / 4" | bc) ;
pyfasta split -n $n "$bed"_plus.fasta ;
pyfasta split -n $n "$bed"_minus.fasta ;
rm *flat ; rm *gdx ;

cat "$bed"_plus.fasta "$bed"_minus.fasta > "$bed"_all.fasta ;


# QUERY LIST

cat "$bed"_plus.fasta | awk '{ if ($1~">") print $1 }' |  sed -e 's|>||g' > query-list-"$bed"_plus.txt ;
cat "$bed"_minus.fasta | awk '{ if ($1~">") print $1 }' |  sed -e 's|>||g' > query-list-"$bed"_minus.txt ;



#######
### NCBI TBLASTN SEARCH - "refseq_genomic" Database
#######

for f in "$bed"_plus.*.fasta;  
do 
	tblastn  -query $f -db ~/EDL933/blastn/refseq_genomic/refseq_genomic -num_threads 4 -out ${f%.*}tBLASTn -evalue 0.001 -max_target_seqs 500000 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos salltitles' & 
done ;
wait ;

for f in "$bed"_minus.*.fasta;  
do 
	tblastn  -query $f -db ~/EDL933/blastn/refseq_genomic/refseq_genomic -num_threads 4 -out ${f%.*}tBLASTn -evalue 0.001 -max_target_seqs 500000 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos salltitles' & 
done ;
wait ;

cat "$bed"_plus.*tBLASTn  | awk '{print $0 "\t" NR}' > "$bed"_plus.tBLASTn.txt  ;
cat "$bed"_minus.*tBLASTn | awk '{print $0 "\t" NR}' > "$bed"_minus.tBLASTn.txt  ;


#######
###Extracting Sequences 
#######

split -dn l/$cores "$bed"_plus.tBLASTn.txt "$bed"_plus.tBLASTn.txt  ;

split -dn l/$cores "$bed"_minus.tBLASTn.txt "$bed"_minus.tBLASTn.txt  ;


# find fasta sequences for all tBLASTn hits 
# remove sequences where full range is not included in database 


for f in "$bed"_plus.tBLASTn.txt[0-9][0-9] ;
do 
	awk -F "\t" '{ split($2,a,"|"); split($1,name,"-");  													 
	if ($9<$10) { print a[4] "\t" ($9-($7*3)+3) "-" $10+((name[2])-($8*3)) "\t" "plus" "\t"  $NF } \
	else if ($10<$9 ) \
	{print a[4] "\t" $10-((name[2])-($8*3)) "-" $9+($7*3)-3 "\t" "minus" "\t" $NF}}' $f |
	awk -F "\t" '{ if (gsub(/-/, "-", $2) == 1) {print $1 "\t" $2 "\t" $3 >"tmp/'$f'.batch"} 
	else if (gsub(/-/, "-", $2) == 2) {print $NF >"tmp/'$f'.excluded"}}'  ; 

	awk 'NR==FNR {a[$1];next} !(($NF) in a)' tmp/"$f".excluded $f > tmp/$f-tBLASTn.included ;
	
	blastdbcmd -db ~/EDL933/blastn/refseq_genomic/refseq_genomic -entry_batch tmp/$f.batch -outfmt %a%f | 
	transeq -filter | awk '!/^>/ { printf "%s", $0; n = "\n" }  /^>/ { print n $0; n = "" } END { printf "%s", n } ' |
	awk '{if ($1 !~ ">") print $0}' | paste -d '\t' tmp/$f-tBLASTn.included - > tmp/bdb-batch-$f & 
done ;

wait ;


for f in "$bed"_minus.tBLASTn.txt[0-9][0-9] ;
do 	
	awk -F "\t" '{ split($2,a,"|"); split($1,name,"-");  													 
	if ($9<$10) { print a[4] "\t" ($9-($7*3)+3) "-" $10+((name[1])-($8*3)) "\t" "plus" "\t" $NF } \
	else if ($10<$9 ) \
	{print a[4] "\t" $10-((name[1])-($8*3)) "-" $9+($7*3)-3 "\t" "minus" "\t" $NF}}' $f |
	awk -F "\t" '{ if (gsub(/-/, "-", $2) == 1) {print $1 "\t" $2 "\t" $3 >"tmp/'$f'.batch"} 
	else if (gsub(/-/, "-", $2) == 2) {print $NF >"tmp/'$f'.excluded"}}'  ; 
	
	awk 'NR==FNR {a[$1];next} !(($NF) in a)' tmp/"$f".excluded $f > tmp/$f-tBLASTn.included ;
	
	blastdbcmd -db ~/EDL933/blastn/refseq_genomic/refseq_genomic -entry_batch  tmp/$f.batch -outfmt %a%f | 
	transeq -filter | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' | 
	awk '{if ($1 !~ ">") print $0}' | paste -d '\t' tmp/$f-tBLASTn.included - > tmp/bdb-batch-$f &

done ;	

wait ;


cat tmp/bdb-batch-"$bed"_plus*[0-9][0-9] > blastdb-batch-"$bed"_plus ;

cat tmp/bdb-batch-"$bed"_minus*[0-9][0-9] > blastdb-batch-"$bed"_minus ; 


	
#######
###Finding Similarity 
#######

#PLUS 

awk 'NR==FNR{a[$1];next} ($1 in a)' "$bed"_plus.tBLASTn.txt query-list-"$bed"_plus.txt > query-list-"$bed"_plus2.txt ;

split -dn l/14 query-list-"$bed"_plus2.txt query-list-"$bed"_plus.txt  ;

for f in query-list-"$bed"_plus.txt[0-9][0-9] ;    
do 
	cat $f | while read -r query etc ;  
do  
	awk '/'$query'/{flag=1; print; next}/>/{flag=0} flag' "$bed"_all.fasta > fastas/$query.fasta ; 
	awk -F "\t" '{if($1=="'$query'" && $NF!~"*") print ">"$(NF-1)"#" "\n" $(NF)}' blastdb-batch-$bed > tmp/$query-seq ;   
	needleall -auto -stdout -aformat3 simple fastas/$query.fasta tmp/$query-seq > tmp/$query-needle ;  

	cat tmp/$query-needle | awk -F "\t" '{if ($1~"# 2:" || $1~"# Similarity:") print $0}' |
	awk 'NR%2{printf "%s ",$0;next;}1' | awk '{print $3 "\t" $7}' > tmp/$query-plus-needle2 ;
	
done & 
done ;

wait ; 


#MINUS 

awk 'NR==FNR{a[$1];next} ($1 in a)' "$bed"_minus.tBLASTn.txt query-list-"$bed"_minus.txt > query-list-"$bed"_minus2.txt ;

split -dn l/14 query-list-"$bed"_minus2.txt query-list-"$bed"_minus.txt

for f in query-list-"$bed"_minus.txt[0-9][0-9] ;    
do 
	cat $f | while read -r query etc ;  
do  
	awk '/'$query'/{flag=1; print; next}/>/{flag=0} flag' "$bed"_all.fasta > fastas/$query.fasta ; 
# Note that the ' $NF!~"*" ' criterion is removed when considering intergenic sequences rather than ORFs  
	awk -F "\t" '{if($1=="'$query'" && $NF!~"*") print ">"$(NF-1)"#" "\n" $(NF)}' blastdb-batch-$bed > tmp/$query-seq ;   
	needleall -auto -stdout -aformat3 simple fastas/$query.fasta tmp/$query-seq > tmp/$query-needle ;   
	
	cat tmp/$query-needle | awk -F "\t" '{if ($1~"# 2:" || $1~"# Similarity:") print $0}' |
	awk 'NR%2{printf "%s ",$0;next;}1' | awk '{print $3 "\t" $7}' > tmp/$query-minus-needle2 ;
	
done & 
done ;

wait ;

cat tmp/*-plus-needle2 | sed -e 's|#||g' > needleall_plus-$bed ;
cat tmp/*-minus-needle2 | sed -e 's|#||g' > needleall_minus-$bed ;


#######
###Match 
#######


awk 'NR==FNR{a[$1];next} (($(NF-1)) in a)' needleall_plus-$bed blastdb-batch-"$bed"_plus > blastdb-batch-"$bed"_plus-filtered ;

awk 'NR==FNR{a[$1];next} (($(NF-1)) in a)' needleall_minus-$bed blastdb-batch-"$bed"_minus > blastdb-batch-"$bed"_minus-filtered ;



paste -d '\t' blastdb-batch-"$bed"_plus-filtered <(sort -n needleall_plus-$bed ) > concat_plus-$bed ;

paste -d '\t' blastdb-batch-"$bed"_minus-filtered <(sort -n needleall_minus-$bed ) > concat_minus-$bed ;

cat concat_plus-$bed concat_minus-$bed > concat-$bed ;



grep -o "'.*'" $tree > tree-names.txt ;
cat tree-names.txt | awk -F "\t" '{split ($1,a,",| "); if (a[3]!="uncultured") print a[3] "\t" substr(a[1],2) "\t" $0 }' > tree-names-genera.txt ;


awk -F "\t" '{split($(NF-4),a," "); print a[1]}' concat-$bed | sort | uniq > concat-$bed-genera ;
 

 
genera=$(awk 'NR==FNR{a[$1];next} (($(NF-1)) in a)' <(awk '{print $1}' tree-names-genera.txt | sort | uniq ) concat-$bed-genera | wc -l ) ;
matches=$(wc -l concat-$bed-genera | awk '{print "'$genera'"/$1 *100}') ;
echo APPROX "$matches"% OF GENERA MAP TO TREE ;
 
 

# Average distance from origin (e.g. E coli) for each genus - average of random sample of up to 10 species for each genus  

	# add genus and percentage at end of line  
cat concat-$bed | 
awk -F "\t" '{split ($(NF-4),a," "); print $0 "\t" a[1] "\t" substr($NF,2, length($NF)-2) }'  > concat2-$bed ;

	# list of genera, restricted to those present in tree, to avoid error message 
cat concat2-$bed | awk -F "\t" '{print $(NF-1)}' | sort | uniq  |
awk 'NR==FNR{a[$1];next} ($1 in a)' tree-names-genera.txt - > genera-$bed ;



	# distances 
split -dn l/$cores genera-$bed genera-$bed ; 

for f in genera-$bed[0-9][0-9] ; 

do cat $f | while read -r genus ; 
do 
	awk '{if ($1=="'$genus'") print $2}' tree-names-genera.txt | shuf -n 10  | while read -r acc ;
do 
	nw_distance -n -m l $tree2 $acc $origin | 
	paste -d "\t"  - -  | awk '{print $2 + $4}' ;
done | 
	awk '{ SUM += $1} END { if (NR>0) print "'$genus'" "\t" SUM/NR}' ;
done > distances-$f & 
done ;

wait ; 
cat distances-genera-$bed[0-9][0-9] > all-LCA-distances.txt ;
 
 
 
# All genera hit for each query ORF, with average similarity for hits in genus    
 
 
 
 
# Distance to each genus by matching with earlier step    
 
 
 
 
# Plotted!  
 

gnuplot 

 
plot "noblast-plot-test.txt" using 1:2:3 with points pointtype 7 pointsize 1.3 palette notitle 
set title "Putative Genes with BLASTp hits" font ", 20"
unset xtics 
set xlabel "ORFs (250)" font ", 15"
set ylabel "Evolutionary Distance (from 16S rRNA tree)" font ", 15"
set colorbox vertical user size 0.1, 0.4
replot






#################################################################

#################################################################

#################################################################


 
 


Find NCBI taxonomy IDs with blastdbcmd 
Map to SILVA tax IDs with Megan file (or file straight from SILVA DB?)

https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/taxonomy/tax_slv_ssu_128.acc_taxid   > maps silva 

 


 
 
 

 
# Add Average Similarity to Gene, for all Subject Sequences for Each Species (e.g. all tTBLASTn hits in E. coli) 
 

for f in out/*distances ; 
do 
	query=$(echo $f | sed -e 's|species-||g' | sed -e 's|-distances||g') ; 
 
	cat int_with_spp-$query | sed -e 's|(||g' | sed -e 's|%||g' | sed -e 's|)||g' > tmp-$query ; 
 
	awk -F "\t" '{if (NF==6) print "'$f'" "\t" $3 "\t" $4+$6}' $f | 
	sed -e 's|species-||g' | sed -e 's|-distances||g' | 
	while read -r a b c ; 
do
	awk '{ if ("'$b'"~$NF) sum += $(NF-2); n++} {if ("'$b'"~$NF) count++; n++}
	END { if (n > 0) print "'$a'" "\t" "'$b'" "\t" "'$c'" "\t" sum / count; }' tmp-$query ; 
done &
done > dist-with-sim.txt ;





cat query-list.txt | grep blast > query-list-blast.txt ;
cat query-list.txt | grep -v blast > query-list-noblast.txt ;

i=0 ;
cat query-list-noblast.txt | while read -r a ; do 
i=$(echo "$i+1" | bc)
awk '{if ($1=="'$a'") print "'$i'" "\t" $(NF-1) "\t" $NF }' dist-with-sim.txt ; 
done > noblast-plot.txt ;

j=0 ;
cat query-list-blast.txt | while read -r a ; do 
j=$(echo "$j+1" | bc)
awk '{if ($1=="'$a'") print "'$j'" "\t" $(NF-1) "\t" $NF }' dist-with-sim.txt ; 
done > blast-plot.txt ;









 

 
 











# get one genome per NCBI species, use this to find NCBI tax ID of species
# map this NCBI tax ID to SILVA tax ID and SILVA species name [using SSU ... synonyms map file and taxmap...128.txt file ]
# map to SILVA LTP names list 



#list of representative species  :     ###### GENERA          
 

cat concat-$bed | awk -F "\t" '{split($(NF-4),a," "); {if (a[3]~"subsp.") print $0 "\t" a[1] "_" a[2] "_" a[3] "_" a[4] "__" ; 
else if (a[2]=="sp.") print $0 "\t" a[1] "_" a[2] "_" a[3] "__" ;
else if (a[3]!~"sp.") print $0 "\t" a[1] "_" a[2] "__" }}' | sort | uniq > concat-species-$bed ; 


cat concat-species-$bed | awk -F "\t" '{print $NF}' | sort | uniq > species-list-$bed ;

cat species-list-$bed | while read -r species; 
do 
	awk -F "\t" '{ if ($NF=="'$species'") print $2}' concat-species-$bed | awk -F "\t" '{split($1,a,"|"); print a[4]}' | head -1 ;
done > genome-reps-list-$bed ;


-db ~/EDL933/blastn/refseq_genomic/refseq_genomic -entry_batch  genome-reps-list-$bed -outfmt %T > taxon-list-$bed

paste -d '\t' species-list-$bed taxon-list-$bed > taxon-list-species-$bed ;


# exact matches against taxmap file: 












# Get SILVA names list  
 
#cat LTPs128_SSU_tree.txt | tail -n +14 | 
#sed -e 's|(||g' | sed -e 's|)||g' | 
#awk -F ":" '{if ($1 !="" && $1 != ";" && $1~"_") print $1}' | 
#awk -F "_" '{if ($3~"subsp.") {print $1 "_" $2 "_" $3 "_" $4 "__" "\t" $0}
# else {print $1 "_" $2 "__" "\t" $0}}' > LTP-names.txt ; 











