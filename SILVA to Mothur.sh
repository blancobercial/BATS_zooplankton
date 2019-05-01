
#to transform the taxonomy file from SILVA to mothur. Worked in SILVA 128. 
awk -F $'\t' 'BEGIN { FS="\t" }; {gsub(/ /,"_")};{print $1"."$2"."$3"\t"$4""$5";"}' taxmap_slv_ssu_ref_128.txt > taxmap_slv_ssu_ref_128_almostREady.txt
#removing the first, empty, line
sed 1d taxmap_slv_ssu_ref_128_almostREady.txt > taxmap_slv_ssu_ref_128_REady.txt


 # to modify the alignment for mothur
 cat YOURSILVAFASTAALIGNMENT.txt | sed 's/ /\t100\t/' > b.txt  # substitute the first space in silva for a tab then, and adds a fake align_ident_slv value (mothur is expecting it; but I would not trust the species ID anyway - blast your reads of interest) and a second tab
 cat b.txt | sed -e 's/ /_/g' > c.txt ## substitutes spaces by underscores in the taxonomy string
 cat c.txt | awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' > READYtoGOfasta.txt
 
#Once with these outfiles, you can select your region using Mothur
  
