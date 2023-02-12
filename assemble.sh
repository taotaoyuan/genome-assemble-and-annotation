
#BUSCO completeness assessments are performed for each assembly generated in the pipeline
#run_BUSCO.py -i sinella.contigs.fasta -c 48 -o Tetranychus_cinnabarinus -m geno -l /home/host/install/busco-3.0.2/lineage/arthropoda_odb9/ 

#quality control for illumina raw data
clumpify.sh in1=illumina.1.fq.gz in2=illumina.2.fq.gz out1=1.clumped.fq.gz out2=2.clumped.fq.gz pigz dedupe 1>>bbtool.log 2>&1
bbduk.sh in1=1.clumped.fq.gz in2=2.clumped.fq.gz out1=1.trim.fq.gz out2=2.trim.fq.gz pigz ref=/home/server/install/bbtools-38.67/resources/adapters.fa qtrim=rl trimq=15 minlen=15 ecco=t maxns=5 trimpolya=10 trimpolyg=10 trimpolyc=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo 1>>bbtool.log 2>&1
mv 1.trim.fq.gz 1.fq.gz && mv 2.trim.fq.gz 2.fq.gz
khist.sh in1=1.fq.gz in2=2.fq.gz histcol=2 khist=khist.txt k=21 peaks=peaks.txt 1>>bbtool.log 2>&1


#genome size estimation with genomescope
#histogram files from khist.txt (bbnorm from illumina) or canu/trimming/0-mercounts/sinella.ms22.histogram (pacbio)
cd 1-survey/
cat ../raw/illumina/khist.txt | sed '1d' | sed "s/	/ /g" > khist.txt
Rscript ~/install/genomescope2-1.0.0/genomescope.R -i khist.txt -o out10000 -k 21 -p 2 -m 10000


#merge and convert the bam file into fasta
samtools merge -@ 111 pacbio.subreads.bam *.bam
samtools fasta -@ 111 pacbio.subreads.bam > pacbio.fa
#all two bam files: 15919614 reads, 235,787,168,300 bp, N/L50 4084402/21.193 KB, mean 14811.11 bp



#assembly with flye
flye --pacbio-raw raw/pacbio/pacbio.fa -g 2000m -o 2-flye -t 111 -i 1 -m 3000 --asm-coverage 50
cat 40-polishing/filtered_contigs.fasta | seqkit seq -m 4999 > contigs.filter5k.fasta


#assemble with wtdbg2
wtdbg2.pl -x sq -t 110 -o wtdbg -g 2g ../raw/pacbio/pacbio.fa
#wtdbg2.pl -x sq -t 110 -o wtdbg -g 2g -p 17 -a -RAS1 ../raw/pacbio/pacbio.fa
#polish using short reads with nextpolish
minimap2 -ax sr wtdbg.cns.fa ../raw/illumina/1.fq.gz ../raw/illumina/2.fq.gz -t 110 | samtools view -bS -@ 110 | samtools sort -@ 110 > map1.illumina.bam
samtools index -@ 110 map1.illumina.bam
samtools faidx wtdbg.cns.fa
python ~/install/NextPolish-1.0.5/lib/nextPolish.py -g wtdbg.cns.fa -t 1 -s map1.illumina.bam -p 110 > wtdbg.sgspolish1.fa

minimap2 -ax sr wtdbg.sgspolish1.fa ../raw/illumina/1.fq.gz ../raw/illumina/2.fq.gz -t 110 | samtools view -bS -@ 110 | samtools sort -@ 110 > map2.illumina.bam
samtools index -@ 110 map2.illumina.bam
samtools faidx wtdbg.sgspolish1.fa
python ~/install/NextPolish-1.0.5/lib/nextPolish.py -g wtdbg.sgspolish1.fa -t 2 -s map2.illumina.bam -p 110 > wtdbg.sgspolish2.fa
#wtdbg2 assembly has larger genome size and higher BUSCO duplicates


#reduce heterozygous regions with 3 rounds of Purge_Dups
cd 3-purge_dups/
minimap2 -xmap-pb -t 110 ../2-flye/contigs.filter5k.fasta ../0-raw/pacbio/pacbio.fa | pigz -c -p 110 > pacbio.paf.gz
pbcstat pacbio.paf.gz
calcuts PB.stat > cutoffs 2>calcults.log

split_fa ../2-flye/contigs.filter5k.fasta > asm.split
minimap2 -xasm5 -DP -t 110 asm.split asm.split | pigz -c -p 110 > asm.split.self.paf.gz
purge_dups -2 -a 50 -f 0.9 -l 5000 -T cutoffs -c PB.base.cov asm.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -p flye dups.bed ../2-flye/contigs.filter5k.fasta


#prepare a mapping file and index both bam and fasta inputs
cd 3-purge_haplotigs/
samtools faidx ../flye/assembly.fasta

minimap2 -ax map-pb -t 110 ../flye/assembly.fasta  ../raw/pacbio/pacbio.fa | samtools view -bS -@ 110 | samtools sort -@ 110 > pacbio.bam
samtools index pacbio.bam -@ 110
#samtools flagstat pacbio.bam
purge_haplotigs hist -b pacbio.bam -g ../flye/assembly.fasta -t 110
purge_haplotigs cov -i pacbio.bam.gencov -o coverage_stats.csv -l 10 -m 42 -h 190 -j 60 -s 60
#l, m and h are read-depth cutoffs for low read depth, the midpoint between the 'haploid' and 'diploid' peaks, and high read-depth cutoff chosen from the histogram .png file; -j and -s can be lower to 60 if high heterozygosity although lower values have weak impacts
purge_haplotigs purge -g ../flye/assembly.fasta -c coverage_stats.csv -b pacbio.bam -t 110  -a 50 -d  #'-r repeats.bed' can be added when repeats marked by windowmasker; -a can be lower to 50
#purge-haplotigs result is better than purge_dups 

#Nextpolish polishing with illumina reads (2 rounds)
cd 4-nextpolish
minimap2 -ax sr ../purge_haplotigs/curated.fasta ../raw/illumina/1.fq.gz ../raw/illumina/2.fq.gz -t 110 | samtools view -bS -@ 110 | samtools sort -@ 110 > map1.illumina.bam
samtools index -@ 110 map1.illumina.bam
samtools faidx ../purge_haplotigs/curated.fasta
python ~/install/NextPolish-1.0.5/lib/nextPolish.py -g ../purge_haplotigs/curated.fasta -t 1 -s map1.illumina.bam -p 110 > purge.sgspolish1.fa

minimap2 -ax sr purge.sgspolish1.fa ../raw/illumina/1.fq.gz ../raw/illumina/2.fq.gz -t 110 | samtools view -bS -@ 110 | samtools sort -@ 110 > map2.illumina.bam
samtools index -@ 110 map2.illumina.bam
samtools faidx purge.sgspolish1.fa
python ~/install/NextPolish-1.0.5/lib/nextPolish.py -g purge.sgspolish1.fa -t 2 -s map2.illumina.bam -p 110 > purge.sgspolish2.fa


##scaffolding with Hi-C data using juicer-3d-dna pipeline
cd 5-chromosome/juicer
mkdir fastq 1-genome 2-restriction_sites
cp FOLDER/*_R*fastq.gz fastq/

cd 1-genome
cat ../../../4-nextpolish/purge.sgspolish2.fa | seqkit seq -s -w 0 | sed -r 's/[Nn]+/\n/g' | cat -n | seqkit tab2fx | seqkit replace -p '.+' -r "Contig{nr}" --nr-width 5 | seqkit seq -m 5000 > genome.fa
bwa index genome.fa

cd 2-restriction_sites
python2 ~/install/juicer-1.6.2/misc/generate_site_positions.py MboI hic ../1-genome/genome.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' hic_MboI.txt > hic.chrom.sizes

cd juicer
bash /home/server/install/juicer-1.6.2/CPU/juicer.sh -d /home/server/diskA/genomes/Trichonephila_antipodiana/5-chromosome/juicer -D /home/server/install/juicer-1.6.2/CPU  -y /home/server/diskA/genomes/Trichonephila_antipodiana/5-chromosome/juicer/2-restriction_sites/hic_MboI.txt -z /home/server/diskA/genomes/Trichonephila_antipodiana/5-chromosome/juicer/1-genome/genome.fa -p /home/server/diskA/genomes/Trichonephila_antipodiana/5-chromosome/juicer/2-restriction_sites/hic.chrom.sizes -s MboI -t 110

#assemble with 3d-dna
cd 3ddna
~/install/3d-dna-180922/run-asm-pipeline.sh -g 100 -m haploid ../juicer/1-genome/genome.fa ../juicer/aligned/merged_nodups.txt
#bad results so we use genome.2 to construct the corrected contigs
~/install/3d-dna-180922/run-asm-pipeline-post-review.sh -g 100 -r genome.2.assembly ../../juicer/1-genome/genome.fa ../../juicer/align
ed/merged_nodups.txt
#too many contigs and thus using salsa to break the possible assembly misjoins

#Hi-C scaffolding with SALSA
#we mainly use its misjoin correction function 
bamToBed -i ../juicer/aligned/alignable.bam > alignment.bed
sort -k 4 alignment.bed > tmp && mv tmp alignment.bed

cp ../juicer/1-genome/genome.fa ./
samtools faidx genome.fa
python ~/install/SALSA-2.2/run_pipeline.py -a genome.fa -l genome.fa.fai -b alignment.bed -e GATC -o scaffolds -m yes -i 1
#the run can be stopped once the corrected contigs generated (assembly.cleaned.fasta or input_breaks)

#continue allhic scaffolding
cp ../salsa/scaffolds/assembly.cleaned.fasta contigs.fa
samtools fastq -@ 100 -1 hic.1.fq -2 hic.2.fq ../juicer/aligned/alignable.bam 
minimap2 -ax sr contigs.fa hic.1.fq hic.2.fq -t 110 | samtools view -hb -@ 110 > alignable.bam

allhic extract alignable.bam contigs.fa --RE GATC
allhic partition alignable.counts_GATC.txt alignable.pairs.txt 13 --minREs 10
for num in $(seq 1 13); do allhic optimize alignable.counts_GATC.13g"$num".txt alignable.clm; done
allhic build alignable.counts_GATC.13g*.tour contigs.fa allhic.chr.fasta
allhic plot alignable.bam *tour  #
#get unplaced contigs
cat allhic.chr.agp | grep "Contig" | cut -f 6 > chr.contigs.list
cat contigs.fa | seqkit grep -v -f chr.contigs.list > unplaced.fa
cat allhic.chr.fasta  unplaced.fa > allhic.fa


#check possible contamminants by blastn

#filter human_mouse contaminants
for FILE in *.fa
  do
    makeblastdb -in $FILE -input_type fasta -dbtype nucl
    windowmasker -in $FILE -infmt blastdb -mk_counts -out $FILE.counts
    windowmasker -in $FILE.counts -sformat obinary -out $FILE.counts.obinary -convert
    hs-blastn index $FILE
    hs-blastn align -db $FILE -window_masker_db $FILE.counts.obinary -query /home/server/diskB/database/NCBI20190601/hm.fa -out $FILE.blast -outfmt 6 -num_threads 110 -evalue 1e-5
    cat $FILE.blast | awk '{print $2}' | sort -n | uniq > seqid_hm.list 
    for seq in $(cat seqid_hm.list); do cat $FILE.blast | grep "$seq" | sort -t $'\t' -k2n -k4nr | head -n 1 >> temp.blast; done
    cat temp.blast | sort -t $'\t' -k4nr > conta_hm.blast
    rm -rf temp.blast
  done

#check the possible contaminants with length shorter than 10000 bp
for id in $(cat seqid_hm.list); do cat input.fa | seqkit grep -r -p "$id" | seqkit stats >> conta_hm.length; done
#no sequence shorter tha 10000 bp

#extract the assembly with contaminants excluded
#seqkit grep -v -f conta1.list pilon.fasta | seqkit seq -u -w 0 > non_conta1.fa
#no obvious hm contaminants

#filter the vector contaminants
#makeblastdb -in UniVec -dbtype nucl -parse_seqids -out univec
blastn -db /home/server/diskB/database/univec/univec -query allhic.fa -out univec.blast -outfmt 6 -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000
#no obvious vectors


##filter the microbe contaminants
for FILE in non_conta2.fa
  do
    makeblastdb -in $FILE -input_type fasta -dbtype nucl
    windowmasker -in $FILE -infmt blastdb -mk_counts -out $FILE.counts
    windowmasker -in $FILE.counts -sformat obinary -out $FILE.counts.obinary -convert
    hs-blastn index $FILE
    hs-blastn align -db $FILE -window_masker_db $FILE.counts.obinary -query /home/server/diskB/database/NCBI20190601/microbes/microbes.fa -out $FILE.blast -outfmt 6 -num_threads 110 -evalue 1e-5
    cat $FILE.blast | awk '{print $2}' | sort -n | uniq > seqid_microbes.list 
    for seq in $(cat seqid_microbes.list); do cat $FILE.blast | grep "$seq" | sort -t $'\t' -k2n -k4nr | head -n 1 >> temp.blast; done
    cat temp.blast | sort -t $'\t' -k4nr > conta_microbes.blast
    rm -rf temp.blast
  done

#check the sequences with very high identity (>90%)
for id in $(cat seqid_microbes.list); do cat non_conta2.fa | seqkit grep -r -p "$id" | seqkit stats >> conta_microbes.length; done

cd microbes/
for id in $(cat ../seqid_microbes.list); do cat ../input.fa | seqkit grep -r -p "$id" > $id.fasta; done
for conta in g8.fasta
do
  cat $conta | seqkit subseq -r 1:26583123 | sed '1d' > t1.fa
  cat $conta | seqkit subseq -r 26611442:198974281 | sed '1d' > t2.fa
  cat $conta | head -n 1 > head.fa
  cat head.fa t*.fa > $conta
  rm head.fa t*.fa
done

seqkit grep -v -f ../seqid_microbes_final.list ../non_conta2.fa | seqkit seq -u -w 0 > temp.fa
cat temp.fa *fasta > ../non_conta3.fa
#some obvious microbes


#get the final assembly

sortbyname.sh in=non_conta3.fa out=temp1.fa length descending
cat temp1.fa | awk '/^>/{print ">Nlau_Chr" ++i; next}{print}' | seqkit seq -u -w 0 | seqkit replace -p .+ -r "Nlau_Chr{nr}" --nr-width 3 > genome.Trichonephila_antipodiana.fa


#map illumina raw data to assembly   97.23%
#map RNA-seq to assembly 96.78%%
#map pacbio long reads to assembly  97.61%



#mask repeats
#build a de novo species specific repeat library
cd 7-repeats/
BuildDatabase -name mydb -engine ncbi ../genome*.fa
RepeatModeler -database mydb -pa 110 -LTRStruct  #-pa=threads-1
#Combine de novo custom library and RepBase
~/install/RepeatMasker-4.0.9/util/queryRepeatDatabase.pl -species all > repeatmasker.taxon.fa
cat repeatmasker.taxon.fa mydb-families.fa > taxon.repeatlib.fa
rm -rf RM* mydb.* repeatmasker.taxon.fa *gz

#mask repeats using RepeatMasker
#this step can be skipped because it will be integrated into MAKER analyses
RepeatMasker -e ncbi -pa 110 -xsmall -lib taxon.repeatlib.fa genome*.fa
#-s as slow search and -q as quick search; -e hmmer can be more accurate but slower
#ouput of soft-masked assembly can be used for downstream BRAKER and MAKER analyses
~/install/RepeatMasker-4.0.9/util/rmOutToGFF3.pl genome.Trichonephila_antipodiana.fa.out > repeats.Trichonephila_antipodiana.gff
~/install/RepeatMasker-4.0.9/util/buildSummary.pl genome.Trichonephila_antipodiana.fa.out > repeats.Trichonephila_antipodiana.summary
#extract complex repeats which will be passed to maker
cat repeats.Trichonephila_antipodiana.gff | grep -v -e "Satellite" -e ")n" -e "-rich" | perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' > repeats.Trichonephila_antipodiana.reformat.gff


#genome-guided aeesmbly of RNA-seq data
#RNA-seq data were downloaded from NCBI (SRR3189772.1, SRR3204354.1, SRR3204356.1, SRR3204357.1)
cd 8-transcripts/
bbduk.sh in1=RNA.1.fq.gz in2=RNA.2.fq.gz out1=1.trim.fq.gz out2=2.trim.fq.gz qtrim=rl trimq=15 minlen=20 ecco=t maxns=5 trimpolya=10 trimpolyg=10 trimpolyc=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo  ref=~/install/bbtools-38.67/resources/adapters.fa  1>>bbtool.log 2>&1

hisat2-build -p 32 ../genome*fa genome
hisat2 -q --dta -x genome -p 32 -1 1.trim.fq.gz -2 2.trim.fq.gz | samtools view -bS -@ 32 | samtools sort -@ 32 > map.rna.bam
samtools  flagstat map.rna.bam > map.rna.stat

#assemble transcripts
stringtie map.rna.bam -o ./trans.stringtie.gtf -p 110 -l Tant_trans
TransBorrow -r trans.stringtie.gtf -g ../genome*fa -b map.rna.bam -s unstranded -t 110
#gffread -w transcripts.stringtie.fa -g ../genome*fa trans.gtf
#stringtie has more unigenes (20430 vs 17344) with comparable completeness value

#reduce duplicated transcripts
#cd-hit-est -i ../NCBI/GCF_004193835.1_ASM419383v1_rna.fna -o  transcripts.NCBI.cdhit.fa -c 0.99 -n 9 -M 28000 -T 32
~/install/redundans-0.13c/redundans.py -f transcripts.stringtie.fa -t 32 --log log.txt --identity 0.9 --overlap 0.8



#Ab Initio Gene Prediction using BRAKER2 (GeneMark+Augustus)
#the aim of this step is to get Genemark and Augustus gene models or ab-initio prediction gff
cd 9-braker
 ~/install/BRAKER-2.1.5/scripts/braker.pl --species=Trichonephila_antipodiana --genome=../7-repeats/genome.Trichonephila_antipodiana.fa.masked --bam=../8-transcripts/map.rna.bam --prot_seq=/home/server/diskB/database/ortho10v1/odb10v1.arthropoda.fa --etpmode --cores 96 --softmasking
#transform gtf annotation file into gff format
gtf2gff.pl < braker.gtf --out=braker.gff3 --gff3


#genome annotation using MAKER
#download  protein sequences from the NCBI
cd 10-gene/maker/
#remove redundant isoforms from proteins
for file in $(ls raw-proteins/); do cd-hit -i raw-proteins/$file -o $file -c 0.99 -n 5 -M 28000 -d 0 -T 16; done
rm *clstr
#create a set of configuration files
cp ../../7-repeats/genome*fa.masked ../../7-repeats/repeats*reformat.gff ../../8-transcripts/transcripts.stringtie.fa ../../9-braker/braker/braker.gff3 ./

maker -CTL
#modify the paramters in maker_opts.ctl and maker_exe.ctl
#genome=genome.Trichonephila_antipodiana.fa.masked
#est=transcripts.stringtie.fa
#protein=Drosophila_melanogaster.fa,Ixodes_scapularis.fa,Stegodyphus_mimosarum.fa,Trichonephila_clavipes.fa,Parasteatoda_tepidariorum.fa,Strigamia_maritima.fa,Daphnia_pulex.fa 
#model_org=
#rm_gff=repeats.Trichonephila_antipodiana.reformat.gff
#pred_gff=braker.gff3
#est2genome=1

#be cautious! The "mpiexec" executed here is from MPICH rather than the default OpenMpi, otherwise Warning will arise.

~/install/maker-2.31.10/exe/mpich2/bin/mpiexec -n 32 maker


cd genome*.maker.output/
fasta_merge -d genome*_master_datastore_index.log
gff3_merge -d genome*_master_datastore_index.log
#maker 2> stderr.log 
#check the posible errors when maker failed "XXX FAILED" in genome*log file
run_BUSCO.py -i *all.maker.proteins.fasta -m prot -l /home/zf/install/busco-3.0.2/lineage/insecta_odb9/ -o Trichonephila_antipodiana_prot -c 16 -sp fly

grep -n "##FASTA" genome*.fa.all.gff
head -n 8518363 genome*.fa.all.gff > genome.maker.without_fasta.gff
#rm -rf *datastore *blastdb tmp* genome*.fa.all.gff
#cd .. && rm -rf *gff *fa *masked *Inline 

#remove redundant isoforms and keep the longest isoforms
mkdir maker1
cat ../genome.*.all.maker.proteins.fasta | seqkit seq -n | cut -d " " -f 1 | sort | grep -v ".t1" | cut -d "." -f 1 | sort | uniq > redundant.gene.abinitio.list
cat ../genome.*.all.maker.proteins.fasta | seqkit seq -n | cut -d " " -f 1 | grep "est2genome-gene" | grep -v "mRNA-1" | sed "s/-mRNA-[0-9]//g" | sort | uniq > redundant.gene.mRNA.list
for loci in $(cat redundant.gene.abinitio.list); do cat ../genome.*.all.maker.proteins.fasta | seqkit grep -r -p "$loci.t" | seqkit seq -n | sort -k9r -t '|' | cut -d " " -f 1 | sed '1d' >> redundant.seq.list; done
for loci in $(cat redundant.gene.mRNA.list); do cat ../genome.*.all.maker.proteins.fasta | seqkit grep -r -p "$loci" | seqkit seq -n | sort -k9r -t '|' | cut -d " " -f 1 | sed '1d' >> redundant.seq.list; done

cat ../genome.maker*gff | grep "     maker   " > temp.gff
cat ../genome.maker*gff | grep "maker     mRNA" | grep -f redundant.seq.list | awk '{print $9}' | cut -d ";" -f 1 | sed "s/ID=//g" > redundant.mRNA.list
#check the possible errors "two mRNAs extracted for single redundant seq"
#for seq in $(cat redundant.seq.list); do num=$(cat ../genome.maker*gff | grep "maker      mRNA" | grep -c "$seq"); echo "$seq     $num" >> temp; done 
for mrna in $(cat redundant.mRNA.list); do sed -i -e "/ID="$mrna"/d" -e "s/,"$mrna"//g" temp.gff; done


#delete proteins of length smaller than 50
cat ../genome.*.all.maker.proteins.fasta | seqkit grep -v -f redundant.seq.list | seqkit seq -n -M 49 | cut -d " " -f 1 > short.list
cat redundant.seq.list short.list > remove.list
cat ../genome.*.all.maker.proteins.fasta | seqkit grep -v -f remove.list > maker.proteins.fasta
cat ../genome.*.all.maker.transcripts.fasta | seqkit grep -v -f remove.list > maker.transcripts.fasta

cat temp.gff | grep "maker	mRNA" | grep -f short.list | awk '{print $9}' | cut -d ";" -f 1 | sed -e "s/ID=//g" -e "s/-mRNA-[0-9]//g" > short.ID.list
cat temp.gff | grep -v -f short.ID.list > maker.gff

#Post Processing of Annotations
cd 9-gene/function/
#revise the names in sequence
maker_map_ids --prefix Trichonephila_antipodiana_ --justify 8 ../maker/genome.*maker.output/maker1/maker.gff > maker.map
cat ../maker/genome.*maker.output/maker1/maker.proteins.fasta | seqkit seq -n | cut -d " " -f 1 > proteins.id.list
for protID in $(cat proteins.id.list); do   mrnaID=$(cat ../maker/genome.*maker.output/maker1/maker.gff | grep "$protID;" | awk '{print $9}' | cut -d ";" -f 1 | sed "s/ID=//g");  cat maker.map | grep "$mrnaID" | sed "s/$mrnaID;/$protID;/g" >> maker.seq.map; done

cp ../maker/genome.*maker.output/maker1/maker.gff maker.renamed.gff
cp ../maker/genome.*maker.output/maker1/maker.proteins.fasta maker.proteins.renamed.fasta
cp ../maker/genome.*maker.output/maker1/maker.transcripts.fasta maker.transcripts.renamed.fasta
map_gff_ids maker.map maker.renamed.gff
for protID in $(cat proteins.id.list); do   mrnaID=$(cat maker.seq.map | grep "$protID	" | cut -d "	" -f 2); sed -i "s/$protID;/$mrnaID;/g" maker.renamed.gff; done
map_fasta_ids maker.seq.map maker.proteins.renamed.fasta
map_fasta_ids maker.seq.map maker.transcripts.renamed.fasta

#annotate gene functions with eggnog
cd function/eggnog/
~/install/eggnog-mapper-2.0.1/emapper.py -i ../maker.proteins.renamed.fasta -m diamond -o eggnog --cpu 110
#egnog results particularly predicted gene names are much fewer than Uniprot-blastps, thus not added to annotations

#assign homology-based gene functions by blastp or diamond against UniProtKB (SwissProt+TrEMBL) databases
cd function/blastp
#pigz -d -p 16 -c ~/database/Uniprot-20190527/uniprot_sprot.fasta.gz > uniprot_sprot.fasta
#pigz -d -p 16 -c ~/database/Uniprot-20190527/uniprot_trembl.fasta.gz > uniprot_trembl.fasta
#cat uniprot*fasta > uniprot.fasta
#diamond makedb --in uniprot.fasta -d uniprot
#diamond blastp -d uniprot -q ../../genome.maker.*output/genome.*all.maker.proteins.fasta -o output.blastp -f 6 --sensitive -p 16 -e 1e-5 --max-target-seqs 1
#diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot
#diamond makedb --in uniprot_trembl.fasta -d uniprot_trembl

diamond blastp -d ~/diskB/database/Uniprot-20200316/uniprot_sprot -q ../maker.proteins.renamed.fasta -o output_sprot.blastp -f 6 --more-sensitive -p 110 -e 1e-5 --max-target-seqs 1
#12666 aligned
cat output_sprot.blastp | awk '{print $1}' | uniq > sprot.list
seqkit grep -v -f sprot.list ../maker.proteins.renamed.fasta > proteins.trembl.fasta

diamond blastp -d ~/diskB/database/Uniprot-20200316/uniprot_trembl -q proteins.trembl.fasta -o output_trembl.blastp -f 6 --more-sensitive -p 110 -e 1e-5 --max-target-seqs 1
#5637
cat output*blastp > output.blastp
rm proteins.trembl.fasta
#makeblastdb -in uniprot.fasta -dbtype prot -out uniprot
#blastp -query maker.all.maker.proteins.fasta -db uniprot -evalue 1e-5 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out output.blastp -num_threads 48


#identify protein domains by InterProscan
#edit maxnumber.of.embedded.workers=(number of threads -1) in 'interproscan.properties' file
cd ../interproscan/
interproscan.sh -dp -b genome.iprscan -f TSV,GFF3 -goterms -iprlookup -pa -t p -cpu 110 -i ../maker.proteins.renamed.fasta -appl Pfam,Panther,Gene3D,Superfamily,CDD


#Post Processing of Annotations
cd maker/function/

#add putative gene functions using UniProt outputs
maker_functional_gff ~/diskB/database/Uniprot-20200316/uniprot.fasta blastp/output.blastp maker.renamed.gff > maker.uniprot.gff #only with UniProt results
#cat maker.proteins.uniprot.fasta | seqkit seq -n | sed "s/ protein Name:\"/\t/g" | sed "s/\"/\;\t/g" | cut -f 1,2 | sed "s/Similar to/Note=Similar to/g" | sed "s/-RA\t/\t/g" > maker.uniprot.map
#for gene in $(cat gene.id); do DESC=$(cat maker.uniprot.map | grep "$gene" | awk 'BEGIN{FS="\t"}{print $2}' | sed "s/\// or /g"); grep -n -E "ID="$gene";|ID="$gene"-RA;" maker.renamed.gff | cut -d ":" -f 1 > line.list; for line in $(cat line.list); do sed -i "$line s/$/Note=Similar to $DESC;/" maker.uniprot.gff; done; done
#sed -i "s/ or /\//g" maker.uniprot.gff
maker_functional_fasta ~/diskB/database/Uniprot-20200316/uniprot.fasta blastp/output.blastp maker.proteins.renamed.fasta | seqkit seq -u -w 0 > maker.proteins.uniprot.fasta
maker_functional_fasta ~/diskB/database/Uniprot-20200316/uniprot.fasta blastp/output.blastp maker.transcripts.renamed.fasta | seqkit seq -u -w 0 > maker.transcripts.uniprot.fasta

#add interproscan protein domain annotation to gffs
#iprscan2gff3 interproscan/genome.iprscan.tsv maker.renamed.gff > maker.iprscan.gff #only with interproscan resuts
ipr_update_gff maker.uniprot.gff interproscan/genome.iprscan.tsv > maker.uniprot_ipr.gff #Uniprot+Interpro gff file

#add ggNOG annotations (KEGG_ko)
cat eggnog/eggnog.emapper.annotations | sed '1,4d' | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$9}' | grep -e "ko:" > eggnog/ko.list
cat eggnog/ko.list | awk '{print $1}' | sed "s/-RA//g" > eggnog/ko.gene.list
cp maker.uniprot_ipr.gff maker.final.gff
for gene in $(cat eggnog/ko.gene.list); do KO=$(cat eggnog/ko.list | grep "$gene" | awk '{print $2}'); grep -n -e "$gene;" maker.final.gff | cut -d ":" -f 1 > line.list; for line in $(cat line.list); do sed -i "$line s/$/KEGG_ko="$KO";/" maker.final.gff; done; done

#give statistics for the gene annotation
cd statistics/

#get cds sequens (protein-->nucleotide)
gffread ../maker.final.gff -g ../../../genome.Trichonephila_antipodiana.fa -x genome.maker.cds.fa -o genome.maker.genes.bed #CDs 19,001, sum 26,729,514 bp, mean 1,406.74 bp

#get exons (cds)
bedtools getfasta -fi ../../../genome.Trichonephila_antipodiana.fa -bed genome.maker.genes.bed -name > test.fa
cat test.fa | seqkit grep -r -p "exon" > exon.fa  #137,611, sum 34,052,800 bp, mean 247.46 bp
cat test.fa | seqkit grep -r -p "CDS" > CDS.fa #135,384, sum 26,729,514 bp, mean 197.43 bp

#get complete gene sequences
grep 'maker     gene' ../maker.final.gff  > gene.gff
sed -i "s/maker gene/maker      exon/g" gene.gff
gffread gene.gff -g ../../../genome.Trichonephila_antipodiana.fa -w gene.fa  #19,001, sum 468,821,229 bp, mean 24673.50 bp

#get intron size
perl get_intron_size.pl ../maker.final.gff | sed "s/0//g" | wc -w  #116,383, sum 434,768,429 bp, mean 3735.67 bp

#get evidence from interproscan
cat ../interproscan/genome.iprscan.tsv | awk '{print $1}' | uniq | wc -l  #14,705
grep "GO:" ../interproscan/genome.iprscan.tsv | awk '{print $1}' | uniq | wc -l #9,329
#grep "KEGG:" ../interproscan/genome.iprscan.tsv | awk '{print $1}' | uniq | wc -l #732
#grep "MetaCyc:" ../interproscan/genome.iprscan.tsv | awk '{print $1}' | uniq | wc -l #622
grep "Reactome:" ../interproscan/genome.iprscan.tsv | awk '{print $1}' | uniq | wc -l # 3,102

#get evidence from Uniprot
cat ../blastp/output.blastp | wc -l  #18,303 blstp, 3,570 "Uncharacterized protein", 698 "unknown function"

#get evidence from eggNOG
cat eggnog.emapper.annotations | sed '1,4d' | awk 'BEGIN{FS="\t";OFS="\t"}{print $6}' | awk NF | wc -l  # 5,508 predicted gene names(sed删除了1-4行的注释行)
cat eggnog.emapper.annotations | grep -c "GO:" #8,862 GOs
cat eggnog.emapper.annotations | sed '1,4d' | awk 'BEGIN{FS="\t";OFS="\t"}{print $8}' | awk NF | wc -l #3,183 EC numbers
-f 9 #9,465 KEGG_ko
-f 10 #5,788 KEGG pathways
-f 21 #14,325 COG Functional Category

#12,226 sequences have GOs (InterPro + eggnog)
cat ../eggnog/eggnog.emapper.annotations ../interproscan/genome.iprscan.tsv | grep "GO:" | awk '{print $1}' | sort | uniq | wc -l

#identigy non-coding RNAs (ncRNAs) with infernal and tRNAscan
mkdir -p 11-ncRNAs/infernal 11-ncRNAs/tRNAscan
#annotate other RNAs with infernal
#download the newest version of Rfam database (Rfam.cm.gz and Rfam.clanin from ftp://ftp.ebi.ac.uk/pub/databases/Rfam)
cd infernal
cmpress Rfam.cm
cmscan --rfam --cut_ga --nohmmonly --tblout rna.tblout --fmt 2 --clanin Rfam.clanin -o rna.cmscan --cpu 16 Rfam.cm ../../genome.*fa

cat family.types | awk 'BEGIN {FS=OFS="\t"}{split($3,x,";");class=x[2];print $1,$2,$3,class}' > rfam_anno_class.txt
awk 'BEGIN{OFS="\t";} {if(FNR==1) print "target_name\taccession\tquery_name\tquery_start\tquery_end\tstrand\tscore\tEvalue"; if(FNR>2 && $20!="=" && $0!~/^#/) print $2,$3,$4,$10,$11,$12,$17,$18; }' rna.tblout > rna.tblout.final
##统计数量
awk 'BEGIN{OFS=FS="\t"} ARGIND==1{a[$2]=$4;} ARGIND==2{type=a[$1]; if(type=="") type="Others"; count[type]+=1;}
 END{for(type in count) print type, count[type];}' rfam_anno_class.txt rna.tblout.final > infernal.result
##统计长度
awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$2]=$4;}ARGIND==2{type=a[$1]; if(type=="") type="Others"; if($6=="-")len[type]+=$4-$5+1;if($6=="+")len[type]+=$5-$4+1}END{for(type in len) print type, len[type];}' rfam_anno_class.txt rna.tblout.final >infernal_length.result

#find tRNAs with tRNAscan
cd tRNAscan
tRNAscan-SE -b trnascan.bed -o trna.out -f trna.structure -s trna.isospecific -m trna.stats -l log.txt -H -y --detail ../../genome*.fa
~/install/tRNAscan-SE-2.0.5/bin/EukHighConfidenceFilter -i trna.out -s trna.isospecific -o . -p filter
cat filter.out | grep "high confidence set" > trna.final.out

#merge infernal and tRNAscan results
cd 11-ncRNAs
cat infernal/rna.tblout.final | grep -v 'tRNA' | sed '1d' | sed -e "s/^/Infernal\t/g" -e "s/^/.\t/g" -e "s/     RF/\tNote=Rfam_accession_RF/g" | awk '{print $5,$2,$3,$6,$7,$9,$8,$1,$4}' | sed "s/ /\t/g"  > infernal.gff3
cat tRNAscan/trna.final.out | awk '{print $1, $2}' | sed "s/^/.tRNA - /g" | awk '{print $3, $1, $4, $2}' | sed "s/ //g" > trna.id
for id in $(cat trna.id); do strand=$(grep "$id" tRNAscan/trnascan.bed | awk '{print $6}'); echo "$id	$strand" | grep "Tant" >> trna.strand; done
paste tRNAscan/trna.final.out trna.strand | sed -e "s/^/tRNAscan\t/g" -e "s/^/.\t/g" | awk '{print $3,$2,$7,$5,$6,$11,$21,$1}' | sed "s/ /\t/g" | sed "s/tRNAscan\t/tRNAscan\ttRNA-/g" > t1
cat tRNAscan/trna.final.out | awk '{print $6}' | sed "s/^/Note=Anticodon_/g" > t2
paste t1 t2 > tRNAscan.gff3
cat *gff3 | sort -k1 -k3 | sed '1i ##the score numbers of column 6 are infernal scores' | sed '1i ##non-coding RNAs in Trichonephila_antipodiana' | sed '1i ##gff-version 3' > ncRNAs.gff3


##infer orthogroups
cd 12-orthofinder/
#simplify the sequence name
ls proteins-raw/ | sed "s/.fa//g" > species.list
for species in $(cat species.list);
  do
    cat proteins-raw/$species.fa | seqkit seq -n | awk '{print $1}' | sed "s/^/$species-/g" > t1 ##-是为了后面串联的时候，保持id一致，把-后面的删掉，前面种名保留
    cat proteins-raw/$species.fa | seqkit seq -s -w 0 > t2
    paste t1 t2 | seqkit tab2fx | seqkit seq -w 0 > $species.fa
    rm t1 t2
  done
#remove th redundant protein sequences
for species in $(cat species.list); do cd-hit -i $species.fa -o proteins/$species.fa -c 0.99 -n 5 -M 280000 -d 0 -T 16; done
sed -i "s/*//g" proteins/*fa && sed -i "s/ //g" proteins/*fa
rm *fa proteins/*clstr

python ~/install/OrthoFinder-2.3.8/orthofinder.py -f proteins -t 16 -a 16 -S diamond -og
mv proteins/OrthoFinder/Results* ./ && rm -rf proteins/OrthoFinder

##extract number of orthogroups used for bar figures
cd Results*/statistics
#the universal orthogroups present in all species
cat ../Orthogroups/Orthogroups.GeneCount.tsv | tail -n +2 | grep -v "        0       " | cut -f 1 > universal_orthogroups.list
#universal single-copy orthogroups/genes
ls ../Single_Copy_Orthologue_Sequences/* | cut -d "/" -f 3 | sed "s/.fa//g" > single-copy.list
for num in $(seq 11); do echo '236' >> 1-single-copy; done
#universal multi-copy orthogroups
cat universal_orthogroups.list | grep -v -f single-copy.list > multi-copy.list
#multi-copy genes/orthologs for each species
for family in $(cat multi-copy.list); do cat ../Orthogroup_Sequences/$family.fa >> temp; done
for species in $(cat ../../species.list); do grep -c "$species" temp >> 2-multi-copy; done
rm temp
##orthogroups only present in all Araneae (spiders) 
cat ../Orthogroups/Orthogroups.GeneCount.tsv | tail -n +2 | cut -d "     " -f 1,6,7,10,11 | grep -v " 0        " | cut -d "    " -f 1 > temp.list
for family in $(cat temp.list); do a=$(grep "$family" ../Orthogroups/Orthogroups.GeneCount.tsv | awk '{print $2}'); b=$(grep "$family" ../Orthogroups/Orthogroups.GeneCount.tsv | awk '{print $3}'); c=$(grep "$family" ../Orthogroups/Orthogroups.GeneCount.tsv | awk '{print $4}'); d=$(grep "$family" ../Orthogroups/Orthogroups.GeneCount.tsv | awk '{print $5}'); e=$(grep "$family" ../Orthogroups/Orthogroups.GeneCount.tsv | awk '{print $8}'); f=$(grep "$family" ../Orthogroups/Orthogroups.GeneCount.tsv | awk '{print $9}'); g=$(grep "$family" ../Orthogroups/Orthogroups.GeneCount.tsv | awk '{print $12}'); test "$a" = 0 -a "$b" = 0 -a "$c" = 0 -a "$d" = 0 -a "$e" = 0 -a "$f" = 0 -a "$g" = 0 && echo $family >> Araneae.orthogroups.list; done
#（test 测试 -a并且 -o或者）对于6.7.10.11	中出现的OG,对应的OG赋值给其他物种，test测试对于的其他物种OG都是否为0，若成了则打印OG.

cat Araneae.orthgroups.list | wc -l  #662
for family in $(cat Araneae.orthogroups.list); do cat ../Orthogroup_Sequences/$family.fa >> temp; done
for species in $(cat ../../species.list); do grep -c "$species" temp >> 3-Araneae; done
rm temp*
#Number of genes in species-specific orthogroups
cat ../Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sed -n '10p' | sed 's/Number of genes in species-specific orthogroups\t//g' | sed "s/\t/\n/g" | col -b > 4-species-specific
#number of genes in the remaining orthogroups
cat ../Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sed -n '3p' | sed "s/Number of genes in orthogroups\t//g" | sed "s/\t/\n/g" | col -b > num1
paste 1-single-copy 2-multi-copy 3-Acaneae 4-species-specific | sed "s/\t/\+/g" | bc > num2
paste num1 num2 | sed "s/\t/-/g" | bc > 5-remaining
rm num*
#Number of unassigned genes
cat ../Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv | sed -n '4p' | sed 's/Number of unassigned genes\t//g' | sed "s/\t/\n/g" | col -b > 6-unassigned_genes

paste ../../species.list 1* 2* 3* 4* 5* 6* | sed  "1i\Species\tSingle-copy orthologs\tMulti-copy orthologs\tCommon genes unique to Araneae\tSpecies-specific genes\tOther orthologs\tUnassigned genes" > statistics.txt


#extract 236 single copy orthologs
cd tree/
ls ../Results*/Single_Copy_Orthologue_Sequences/ | cut -d "/" -f 3 | sed "s/.fa//g" > loci.list

cd tree/0-raw/
for loci in $(cat ../loci.list)
  do
    cat ../../Results*/Single_Copy_Orthologue_Sequences/$loci.fa | seqkit seq -n | cut -d "-" -f 1 > t1 ##把-后面的本身序列名字删除，只要前面的种名好串联
    cat ../../Results*/Single_Copy_Orthologue_Sequences/$loci.fa | seqkit seq -s -w 0 > t2
    paste t1 t2 | seqkit tab2fx > $loci.fa
    rm t1 t2
  done

#phylogenetic tree
cd tree/
mkdir 1-mafft 2-trim 3-symtest 4-ML
#align the sequences
for gene in $(sed -n '1,300p' loci.list); do linsi --thread 5  0-raw/$gene.fa > 1-mafft/$gene.mafft.fa; done

#trim the alignments
cd 2-trim
for gene in $(cat ../loci.list); do java -jar ~/install/BMGE-1.12/BMGE.jar -i ../1-mafft/$gene.mafft.fa -t AA -of $gene.fas >> ../log.trim; done
perl ~/install/FASconCAT-G/FASconCAT-G*pl -s -l

#symtest (IQTREE) to exclude partitions that violate models (i.e. that may have evolved under non-SRH conditions) prior to tree reconstruction 
cd 3-symtest
mv ../2-trim/FcC* ./
iqtree -s FcC_supermatrix.fas -p FcC_supermatrix_partition.txt --symtest-remove-bad --symtest-pval 0.05 #stop the computation once the symtest finished, p-value can be 0.01-0.1
cat FcC_supermatrix_partition.txt.bad.nex | grep "charset" | cut -d " " -f 4 | sed "s/.aa//g" > symtest.loci.list
#27 loci were removed and 209 were left
rm *gz

#reconstruct the partitioned ML tree
cd 4-ML/loci
for loci in $(cat ../../3-symtest/symtest.loci.list); do cp ../../2-trim/$loci.fas ./; done
perl ~/install/FASconCAT-G/FASconCAT-G*pl -s -l
mv FcC* ../ && cd ..

#reroot the outgroup as the first taxon in the fas file
    echo "Tachypleus_tridentatus" > list
    cat *fas | seqkit grep -f list | seqkit seq -u -w 0 > outgroup.fa
    cat *fas | seqkit grep -v -f list | seqkit seq -u -w 0 > ingroup.fa
    cat outgroup.fa ingroup.fa > FcC_supermatrix.fas
    rm ingroup* outgroup* list

iqtree -s FcC_supermatrix.fas -p FcC_supermatrix_partition.txt -m MFP --mset LG --msub nuclear --rclusterf 10 -B 1000 --alrt 1000 -T AUTO  ##主要LG是蛋白模型，核酸可以选GTR

#estimate divergence time using mcmctree
cd 5-mcmctree
#reroot the tree using the outgroup
#transform sequence format into PAML-phylip using readal in trimal package
for loci in $(cat ../3-symtest/symtest.loci.list); do readal -in ../2-trim/$loci.fas -out 0-phylip/$loci.phy -phylip_paml; done (-phylip_paml在最新的版本中是指的-phylip3.2格式)
cat 0-phylip/*phy > input.phy
#modify the calibration in the tree, check the possible nodes/calibrations errors using FigTree
#note: The time unit is 100 Million years (Myr)
#for the large data, use the approximation likelihood calculation instead of "usedata = 1"

#run using "usedata = 3" to construct an approximation to the likelihood function
cd 1-BV
mcmctree mcmctree.ctl

#for the protein sequences
rm out.BV rst*
cp ~/install/paml-4.9j/dat/lg.dat ./
for file in tmp*ctl; do sed -i "s/model = 0/model = 2/g" $file; sed -i "s/aaRatefile = /aaRatefile = lg.dat/g" $file; sed -i "s/method = 0/method = 1/g" $file; done
for file in tmp*ctl; do codeml $file; cat rst2 >> in.BV; done

#re-run using "usedata = 2" to finish the final step
cd 2-divergence
mcmctree mcmctree.ctl
#repeat the run at least twice
#check the convergence "mcmc.txt" using TRACER
#final results in FigTree.tre

#estimate divergence using r8s
cd 5-r8s
python ~/install/CAFE-4.2.1/scripts/python_scripts/cafetutorial_prep_r8s.py -i input.tre -o r8s_ctl_file.txt -s 75868 -p 'Centruroides_sculpturatus,Dermatophagoides_pteronyssinus,Galendromus_occidentalis,Ixodes_scapularis,Parasteatoda_tepidariorum,Stegodyphus_mimosarum,Tachypleus_tridentatus,Tetranychus_urticae,Trichonephila_antipodiana,Trichonephila_clavipes,Varroa_destructor' -c '541'

## seqkit stat MultipleSequenceAlignments/SpeciesTreeAlignment.fa(用作-s的参数)
## seqkit stat 2-trim/*| awk '{sum += $8}; END {print sum}' ##1081641

#modify the ctl file
r8s -b -f r8s_ctl_file.txt > r8s_tmp.txt
tail -n 1 r8s_tmp.txt | cut -c 16- > r8s_ultrametric.txt



#gene family evolution using CAFE4 and 5
#prepare a binary, rooted, ultrametric tree
cd 13-cafe
cp ../12-orthofinder/tree/5-mcmctree/2-divergence/run2/FigTree.tre ./  #revise the tree format and the time scale (X100)
#prepare gene count data file
cat ../12-orthofinder/Results*/Orthogroups/Orthogroups.GeneCount.tsv | sed "s/^/\t/g" | sed "1s/\tOrthogroup/Description\tID/" | cut -d " " -f 1-12 > data.txt
cafexp -i data.txt -t FigTree.tre -p -k 3 -o gamma
cafexp -i data.txt -t FigTree.tre -p -o base
#get the lambda value for CAFE4 analyses

#cafe4
##no "_" in species name
cafe
load -i data.rename.txt -p 0.01 -t 16 -l log.txt -r 1000
tree (((((Galendromus:146.7892,Varroa:146.7892):221.4606,Ixodes:368.2497):84.4784,(Dermatophagoides:330.0525,Tetranychus:330.0525):122.6757):30.8185,((((Tantipodiana:17.8274,Tclavipes:17.8274):101.7749,Parasteatoda:119.6024):24.6833,Stegodyphus:144.2857):300.3262,Centruroides:444.6119):38.9347):11.1995,Tachypleus:494.7461)
lambda -s -t (((((1,1)1,1)1,(1,1)1)1,((((1,1)1,1)1,1)1,1)1)1,1) #lambdamu can be used to separate birth and death parameters
#lambdamu -s -t (((((1,1)1,(((1,1)1,(((((1,1)1,1)1,((1,1)1,1)1)1,1)1,1)1)1,1)1)1,1)1,1)1,1) 
lambda -l 0.001265590477607
#lambdamu -l 0.00021776549904 -m 0.00138349331557
report report.lambda

python /home/zf/install/CAFE-4.2.1/python_scripts/cafetutorial_report_analysis.py -i report.lambda.cafe -o summary_lambda
python /home/zf/install/CAFE-4.2.1/scripts/python_scripts/cafetutorial_draw_tree.py -i summary_lambda_node.txt -t '(((((Galendromus:146.7892,Varroa:146.7892):221.4606,Ixodes:368.2497):84.4784,(Dermatophagoides:330.0525,Tetranychus:330.0525):122.6757):30.8185,((((Tantipodiana:17.8274,Tclavipes:17.8274):101.7749,Parasteatoda:119.6024):24.6833,Stegodyphus:144.2857):300.3262,Centruroides:444.6119):38.9347):11.1995,Tachypleus:494.7461)' -d '(((((Heliconius<0>,Danaus<2>)<1>,((((((Heliothis<4>,Helicoverpa<6>)<5>,Spodoptera<8>)<7>,Trichoplusia<10>)<9>,Operophtera<12>)<11>,(Bombyx<14>,Manduca<16>)<15>)<13>,(Ostrinia<18>,Amyelois<20>)<19>)<17>)<3>,Plutella<22>)<21>,Stenopsyche<24>)<23>,Drosophila<26>)<25>' -o tree_rapid.pdf -y Rapid

cd expansion
cat ../cafe4/summary_lambda_fams.txt | grep "Tantipodiana<" | sed "s/,/\n/g" | grep "+" | cut -d "[" -f 1 > expansion.list #revise the first line
#113 significantly expanded families




##自己编的程序
cat ../cafe4/summary_lambda_fams.txt | grep "Tantipodiana<" | sed "s/,/\n/g" | grep "-" | cut -d "[" -f 1 > contruction.list
##自己编的程序
for id in  $(cat expansion.list)
 do
  cat report.lambda.cafe| grep "$id"  >> expansion.result
 done
 ###自己编的程序
 cat expansion.result | cut -d "	" -f 4 > evalue.tree
 cat evalue.tree | sed 's/(//g' |  sed 's/)//g' |  sed 's/,/	/g' | cut - f 11 > evalue
将 evalue 值进行排序，查看其小于0.01的值，即为显著的
 一共找到了
 
 

 
 ###



for id in $(cat expansion.list)
 do
  cat ../../12-orthofinder/Results*/Orthogroups/Orthogroups.txt | grep "$id" | sed "s/ /\n/g" | grep "Trichonephila_antipodiana" > $id.list
  echo "$id	$(cat $id.list | wc -l)" >> gene.number
 done
sed -i "s/Trichonephila_antipodiana-//g" OG*list

cat ../../10-gene/function/maker.final.gff | grep "maker	gene" > gene.gff3
for id in $(cat expansion.list)
 do
   for gene in $(cat $id.list | sed "s/-RA//g"); do cat gene.gff3 | grep "$gene" | cut -d ";" -f 4 | sed "s/Note=Similar to//g" >> $id.function; done
   for gene in $(cat $id.list); do cat genome.iprscan.tsv | grep "$gene" | grep "Pfam" | cut -d "	" -f 6 >> $id.pfam; done
done


#Gene Enrichment analyses for significantly expanded gene families
cd 14-enrichment
#GO enrichment
cd GO-expansion
#obtain the background file (the whole genome)
cat ../eggnog.emapper.annotations ../genome.iprscan.tsv | grep "GO:" | awk '{print $1}' | sort | uniq > all_gene_with_go.list
for gene in $(cat all_gene_with_go.list)
  do
    cat ../eggnog.emapper.annotations | grep "$gene" | grep "GO:" | cut -d "        " -f 7 | sed "s/,/\n/g" > temp.list
    cat ../genome.iprscan.tsv | grep "$gene" | grep "GO:" | cut -d "   " -f 14 | sed "s/|/\n/g" >> temp.list
    cat temp.list | sort | uniq | sed "s/$/\t$gene/g" >> term2gene.txt
    rm temp.list
  done
##最新的eggnog的结果中GOs在第10行



#get description from the TBtools GO name Parse function or extract from go-basic.obo (website of GO)
#cat go-basic.obo | sed -n "/^id: GO:/,/^namespace:/p"| awk '{if(NR%3!=0)ORS=" ";else ORS="\n"}1' | sed "s/id: //g" | sed "s/ name: /\t/g" | sed "s/ namespace: /\t/g" > GO.descriptions
cat ~/diskB/database/GO/GO.descriptions | cut -d "  " -f 1,2 > term2name.txt

#get the list of candidate (target) gene set
for family in $(cat ../../13-cafe/expansion/expansion.list); do cat ../../12-orthofinder/Results*/Orthogroup_Sequences/$family.fa | grep "Trichonephila_antipodiana" | sed "s/>Trichonephila_antipodiana-//g"  >> temp1.list; done
cat all_gene_with_go.list | grep -f temp1.list | tr "\n" "," | sed "s/,$//g" > gene.txt
rm temp1.list

#R package 'clusterProfiler'
gene <- read_delim("C:/Users/nj/Desktop/GO/raw/gene.txt", ",", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
term2gene <- read_delim("C:/Users/nj/Desktop/GO/raw/term2gene.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
term2name <- read_delim("C:/Users/nj/Desktop/GO/raw/term2name.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
x <- enricher(gene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
barplot(x, showCategory=20)
dotplot(x, showCategory=20)

#KEGG enrichment
#obtain the background file (the whole genome)
cd KEGG-expansion
cat ../eggnog.emapper.annotations | tail -n +5 | grep ",map" | cut -d "  " -f 1 > all_gene_with_ko.list
###yt
大空格：ctrl + v + i   不能复制，会变格式（4个空格）
###
for gene in $(cat all_gene_with_ko.list); do cat ../eggnog.emapper.annotations | grep "$gene" | cut -d " " -f 10 | sed "s/,map/\t/g" | cut -d "   " -f 1 | sed "s/,/\n/g" | sed "s/$/\t$gene/g" >> term2gene.txt; done
###
for gene in $(cat all_gene_with_ko.list); do cat ../eggnog.emapper.annotations | grep "$gene" |cut -d "   " -f13| sed "s/,map/\t/g"|cut -d "      " -f1|sed 's/,/\n/g'| sed "s/$/\t$gene/g" >term2gene;done
---
for gene in $(cat all_gene_with_ko.list); do cat ../eggnog.emapper.annotations | grep "$gene"|cut -d "    " -f13| sed "s/,map/\t/g" | cut -d "    " -f 1| sed "s/,/\n/g"|sed '/-/d'|sed "s/$/\t$gene/g" >> term2gene.txt; done
###

以上kegg文件准备可以用R做
```
##1.先把需要的行读到eggnogmapper中
eggnogmapper<- read_delim(
  file = 'eggnog.emapper.annotations',
  "\t",
  escape_double =FALSE,
  col_names=FALSE,
  comment = "#",trim_ws = TRUE)%>%
  dplyr::select(GID=X1,
                COG=X7
                SYMBOL=X9,
                GO=X10,
                KO=X12,
                Pathway=X13,
                OG=X7,
                Gene_Name=X21)
##2.把Pathway和GID提取出来
temp = dplyr::select(eggnogmapper,Pathway,GID) ##此时为一个GID对应多个Pathway，我们希望1对1，多行。
此时需要dplyr包中的separate_rows（）来对行进行拆分
emp = dplyr::select(eggnogmapper,Pathway,GID) %>% separate_rows(Pathway, sep = ',',convert =F) 对Pathway列以','为分隔符进行分割，convert =F指不自动将类似ko03008换成数字。

同时去掉map开头的行以及没有对应上的行
filter(str_detect(Pathway,'ko')) ##检测Pathway行是不是包含ko字符，TURE则保留，

temp = dplyr::select(eggnogmapper,Pathway,GID) %>% separate_rows(Pathway, sep = ',',convert =F) %>% filter(str_detect(Pathway,'ko'))
最后去掉ko字符
mutate(Pathway=str_remove(Pathway,'ko'))
汇总：

Pathway2gene = dplyr::select(eggnogmapper,Pathway,GID) %>% separate_rows(Pathway, sep = ',',convert =F) %>% filter(str_detect(Pathway,'ko'))%>%mutate(Pathway=str_remove(Pathway,'ko'))
```







#get the list of candidate (target) gene set
for family in $(cat ../../13-cafe/expansion/expansion.list); do cat ../../12-orthofinder/Results*/Orthogroup_Sequences/$family.fa | grep "Trichonephila_antipodiana" | sed "s/>Trichonephila_antipodiana-//g"  >> temp.list; done
cat all_gene_with_ko.list | grep -f temp.list | tr "\n" "," | sed "s/,$//g" > gene.txt
rm temp.list

#get KO orthology from https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir=
cat ko00001.keg | grep '\[PATH:' | sed "s/^.\{11\}//g" | sed "s/]//g" | sed "s/ \[PATH:/\t/g" | cut -d "   " -f 2 > c1
###yt
cat ko00001.keg | grep '\[PATH:' | sed "s/^.\{11\}//g" |cut -d ':' -f2|sed 's/]//g'>c1
###
cat ko00001.keg | grep '\[PATH:' | sed "s/^.\{11\}//g" | sed "s/]//g" | sed "s/ \[PATH:/\t/g" | cut -d "   " -f 1 > c2
###yt
cat ko00001.keg | grep '\[PATH:' | sed "s/^.\{11\}//g" |cut -d '[' -f1 >c2
###
paste c1 c2 > term2name.txt
rm c1 c2

#R package 'clusterProfiler', same as above

#GO enrichment for species-specific families
cd GO-species-specific
cp ../GO-expansion/term* ../GO-expansion/all* ./
cat ../../12-orthofinder/Results*/Orthogroups/Orthogroups.GeneCount.tsv | tail -n +2 | cut -d "     " -f 1 > temp.list
for family in $(cat temp.list); do a=$(grep "$family" ../../12-orthofinder/Results*/Orthogroups/Orthogroups.GeneCount.tsv | cut -f10); b=$(grep "$family" ../../12-orthofinder/Results*/Orthogroups/Orthogroups.GeneCount.tsv | cut -f13 | col -b); test "$a" = "$b" && echo $family >> species-specific.orthogroups.list; done
for family in $(cat species-specific.orthogroups.list); do cat ../../12-orthofinder/Results*/Orthogroup_Sequences/$family.fa >> temp; done
cat temp | seqkit seq -n | sed "s/Trichonephila_antipodiana-//g" | grep -f all_gene_with_go.list | tr "\n" "," | sed "s/,$//g" > gene.txt
rm temp*

#KEGG enrichment for species-specific families
cd KEGG-species-specific
cp ../KEGG-expansion/term* ../KEGG-expansion/all* ./
for family in $(cat ../GO-species-specific/species-specific.orthogroups.list); do cat ../../12-orthofinder/Results*/Orthogroup_Sequences/$family.fa >> temp; done
cat temp | seqkit seq -n | sed "s/Trichonephila_antipodiana-//g" | grep -f all_gene_with_ko.list | tr "\n" "," | sed "s/,$//g" > gene.txt


##synteny analyses using Tbtools
#Tant vs Ttri
cd 15-synteny


##WGD (whole genome duplication) analyses

#ks-based wgd
wgd dmd genome.maker.cds.fa
wgd ksd -n 16 wgd_dmd/genome.maker.cds.fa.mcl genome.maker.cds.fa
wgd kde -b 100 -r 0 3 wgd_ksd/genome.maker.cds.fa.ks.tsv
wgd mix --method bgmm -b 50 wgd_ksd/genome.maker.cds.fa.ks.tsv > log.mix_bgmm
#view the figures manually
bokeh serve &
wgd viz -i -ks wgd_ksd/genome.maker.cds.fa.ks.tsv


##extract gene families from genome and proteins
cd 17-gene-families

#download reference protein sequences from NCBI, such GST famiy

#blastp-like search against proteinsDB
##mmseqs easy-search Bmori_GST.fasta genome.fa blastp.out tmp -s 7.5 --alignment-mode 3 --num-iterations 4 -e 0.001 --format-mode 2

#erect query and target DB
mmseqs createdb genome.fa genomeDB
mmseqs createdb proteins.fa proteinsDB

#erect index file for the targetDB
mmseqs createindex genomeDB tmp
mmseqs createindex proteinsDB tmp

cd GST
mmseqs createdb ../Bmori_GST.fasta queryDB
cd GST/AA
mmseqs search ../queryDB ../../proteinsDB result_AADB tmp -s 7.5 --alignment-mode 3 --num-iterations 4 -e 0.001       #(-e   10     1e-5) , --min-seq-id 0.5 
mmseqs convertalis ../queryDB ../../proteinsDB result_AADB blastp.out --format-mode 2 
cat blastp.out | awk '{print $2}' | sort | uniq > blastp.AA.list
cat ../../proteins.fa | seqkit grep -f blastp.AA.list | seqkit seq -w 0 > blastp.AA.fa
cd-hit -i blastp.AA.fa -o blastp.AA.cdhit95.fa -c 0.95 -n 5 -M 12000 -d 0 -T 8  #0.99 for CSP and very conserved genes

#check GSTs by HMMER3 search (cutoff E-value = 0.001) using the Pfam database to confirm conserved domains (online or local interproscan)
HMMER-Pfam results were further checked by blastp in the nonredundant GenBank protein database.
#cat hmmer-pfam.results | grep "GST_" | awk '{print $1}' | sed "s/>//g" | sort | uniq > hmmer-pfam.list
#cat blastp.AA.fa | seqkit grep -f hmmer-pfam.list | seqkit seq -u -w 0 > hmmer-pfam.fa
#cd-hit -i hmmer-pfam.fa -o hmmer-pfam.cdhit99.fa -c 0.99 -n 5 -M 12000 -T 16
#check sequences too short or too long, 40 remaining


#tblatn-like search against genomeDB using mmseq2
cd GST/genome
mmseqs search ../queryDB ../../genomeDB result_genomeDB tmp -s 7.5 --alignment-mode 3 --num-iterations 4 -e 10
mmseqs convertalis ../queryDB ../../genomeDB result_genomeDB tblastn.out --format-output 'query,target,qstart,qend,qlen,tstart,tend,tlen,evalue,alnlen,pident,qaln,taln'
#mmseqs easy-search template.fasta genome.picauta_chinensis.fa tblastn_aln.cce tmp -s 7.5 --alignment-mode 3 --num-iterations 4 -e 0.001 --format-output 'query,target,evalue,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,mismatch,bits,qcov,tcov,qframe,tframe,qaln,taln'

#sort by target, qstart and length using excel or csvtk
cat tblastn.out | csvtk sort -H -t -k2 -k6:n -k9:n | nl -b a -n rz | sed "s/^/>/g" > tblastn.temp1
#reduce the multi-hits to the same gene region
cat tblastn.temp1 | awk '{print $1,$14}' | sed "s/-//g" | sed "s/ /\n/g" > tblastn.temp2.fa
cd-hit -i tblastn.temp2.fa -o tblastn.temp2.cdhit99.fa -c 0.99 -n 5 -M 12000 -T 8
cat tblastn.temp2.cdhit99.fa | grep ">" | sed "s/>//g" | sort > tblastn.temp3.list
for list in $(cat tblastn.temp3.list); do cat tblastn.temp1 | grep ">$list" >> tblastn.temp4;done #delete invalid matches
#filter the identical regions matching GSTs prodicted from annotated proteins 
mmseqs easy-search tblastn.temp2.cdhit99.fa ../AA/blastp.AA.cdhit95.fa tblastn.blastp.out tmp -s 7.5 --alignment-mode 3 -e 0.001 --format-output 'query,target,qstart,qend,qlen,tstart,tend,tlen,evalue,alnlen,pident,qaln,taln' --min-seq-id 0.9   #0.99 for CSP
cat tblastn.blastp.out | awk '{print $1}' | sort | uniq > matched.list
for list in $(cat matched.list); do sed -i "s/>"$list"/>"$list"***/g" tblastn.temp4; done

#manually merge exons
cat ../genome.fa | seqkit grep -r -p "Sexi_chr31" | seqkit subseq -r 7359364:7359604 | seqkit seq -w 0 > test.fasta
#accurately determine Reading frames and intron/exon boundaries
#align candidates to annotated proteins or use tools GeneWise, exonerate or GenomeThreader
gth -xmlout -o gth.xml -genomic test.fasta -protein ../Bmori_GST.fasta -cdna transcripts.fa
cat gth.xml | grep "<predicted_protein_sequence>" | sed "s/<predicted_protein_sequence>//g" | sed "s/<\/predicted_protein_sequence>//g" >> test.aa.fasta 

#compare newly preidicted genes to AA
mmseqs easy-search unmatched.fa ../AA/blastp.AA.cdhit95.fa unmatched.blastp.out tmp -s 7.5 --alignment-mode 3 -e 0.001 --format-output 'query,target,qstart,qend,qlen,tstart,tend,tlen,evalue,alnlen,pident,qaln,taln'
cd ..
cat genome/unmatched.fa AA/blastp.AA.cdhit95.fa > GSTs.fa   #44 totally
#cd-hit -i P450s.fa -o P450s.cdhit95.fa -c 0.95 -n 5 -M 12000 -T 8

#check their positions on chromosomes
cd tree
#cat ../CCEs.fa | sed "s/ protein /_protein_/g" | sed "s/_protein_Name/ /g" > CCEs.fa
linsi --thread 8 ../ABCs.fa > ABC.mafft.fas
trimal -in ABC.mafft.fas -out ABC.trim.fas -automated1
iqtree -s input.fas -mset LG -bb 1000 -nt AUTO

mmseqs easy-search ../ABCs.fa ../../genome.fa blastp.out tmp -s 7.5 --alignment-mode 3 -e 0.001 --format-output 'query,target,qstart,qend,qlen,tstart,tend,tlen,evalue,alnlen,pident,qaln,taln' --min-seq-id 0.9


##circos figures for genome characteristics using Tbtools

#get the Chr length
cat ../genome.Trichonephila_antipodiana.fa | seqkit grep -r -p "_Chr" | sed "s/Ofur_//g" > genome.fa
cat genome.fa | seqkit fx2tab -n -l | csvtk sort -k1 | sed "s/\t\t\t/\t/g" > Chr.length

#get the GC content
cat genome.fa | seqkit sliding -s 1000 -W 1000 | seqkit fx2tab -n -g | sed "s/_sliding:/\t/g" | sed "s/-/\t/g" | sed "s/\t\t\t/\t/g" > GC.tab

#get the gene information
cat maker.final.gff | grep "maker  gene" | grep 'Ofur_Chr' | cut -f 1-5 | sed "s/maker\tgene/1/g" | awk '{OFS="\t"}{print $1,$3,$4,$2}' > gene.txt

