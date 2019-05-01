#make contigs with the file created according to mothur instructions; only accept bases unique to a strand if Q>30; trim the overlaps (useful in V9)
make.contigs(file=BZall_Clean2.txt, insert=30, trimoverlap=T, processors=60)
summary.seqs(fasta=current)
#remove short reads, or those with any ambiguity
screen.seqs(fasta=current, group=current, qfile=current, minlength=115, maxambig=0)
summary.seqs(fasta=current)
unique.seqs(fasta=current)
summary.seqs(fasta=current, name=current)
count.seqs(name=current, group=current)
summary.seqs(fasta=current, name=current)
chimera.uchime(fasta=current, name=current, group=current, dereplicate=t)
remove.seqs(fasta=current, accnos=current, count=current, dups=f)
summary.seqs(fasta=current, name=current)
#aling to the fasta created with all the SILVA sequences
align.seqs(fasta=current, template=Silva128_ShortAll.fas) 
#if the alignment is not of 692 sequences, you might need to restart here (pcr.seqs will fail). This will allow you to feed your files here
get.current()
summary.seqs(fasta=current)
#remove all extra bases
pcr.seqs(fasta=current, taxonomy=taxmap_slv_ssu_ref_128_ShortAll.txt, keepdots=false, start=0, end=692, count=current)
#the next two steps will leave only perfect reads, that cover the whole region
#remove all reads that do not reach the end
screen.seqs(fasta=current, count=current, end=692)
#remove all reads that do not start at the beggining
screen.seqs(fasta=current, count=current, start=1)
degap.seqs(fasta=current)
unique.seqs(fasta=current, count=current)
chimera.uchime(fasta=current, count=current, dereplicate=t)
remove.seqs(fasta=current, accnos=current, count=current, dups=f)
#pairwise - this is very time consuming, I would not recommend do this anymore. There are better options in the new mothur versions
pairwise.seqs(fasta=current, countends=F, calc=nogaps)
get.current()
#clustering at 0.016 will cluster sequences that are up to 2 bp different - the centroid and 1 bp error around
cluster(column=BZall_Clean2.trim.contigs.good.unique.pick.pcr.good.good.ng.unique.pick.dist, count=BZall_Clean2.trim.contigs.good.unique.pick.pcr.good.good.ng.denovo.uchime.pick.pick.count_table, method=opti, cutoff=0.016 precision=1000)
make.shared(list=BZall_Clean2.trim.contigs.good.unique.pick.pcr.good.good.ng.unique.pick.opti_mcc16.list, count=BZall_Clean2.trim.contigs.good.unique.pick.pcr.good.good.ng.denovo.uchime.pick.pick.count_table, label=0.016)
classify.seqs(fasta=BZall_Clean2.trim.contigs.good.unique.pick.pcr.good.good.ng.unique.pick.fasta, template=Silva128_ShortAll.ng.fasta, taxonomy=taxmap_slv_ssu_ref_128_ShortAll.txt, count=current)
classify.otu(taxonomy=current, count=current, list=current, label=0.016)