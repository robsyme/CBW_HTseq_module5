#set up


docker run --privileged -v /tmp:/tmp --network host -it -w $PWD -v $HOME:$HOME -v /media:/media --user $UID:$GROUPS -v /etc/group:/etc/group -v /etc/passwd:/etc/passwd c3genomics/genpipes:0.8
export WORK_DIR_M5=$HOME/workspace/HTseq/Module5/
export REF=$HOME/workspace/HTseq/Module5/reference
mkdir -p $WORK_DIR_M5
cd $WORK_DIR_M5
ln -s $HOME/CourseData/HT_data/Module5/* .

module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/samtools/1.9 mugqic/GenomeAnalysisTK/4.1.0.0 mugqic/R_Bioconductor/3.5.0_3.7 mugqic/bcftools/1.6 mugqic/tabix/0.2.6 mugqic/Delly/0.7.8



samtools view bam/NA12878/NA12878_S1.chr20.20X.pairs.readSorted.bam \
     | python scripts/pairend_distro.py \
     -r 101 \
     -X 4 \
     -N 10000 \
     -o SVvariants/NA12878_S1.chr20.20X.pairs.histo \
 > SVvariants/NA12878_S1.chr20.20X.pairs.params


#NA12891 
samtools view bam/NA12891/NA12891_S1.chr20.20X.pairs.readSorted.bam \
     | python scripts/pairend_distro.py \
     -r 101 \
     -X 4 \
     -N 10000 \
     -o SVvariants/NA12891_S1.chr20.20X.pairs.histo \
 > SVvariants/NA12891_S1.chr20.20X.pairs.params


#NA12892. 
samtools view bam/NA12892/NA12892_S1.chr20.20X.pairs.readSorted.bam \
     | python scripts/pairend_distro.py \
     -r 101 \
     -X 4 \
     -N 10000 \
     -o SVvariants/NA12892_S1.chr20.20X.pairs.histo \
 > SVvariants/NA12892_S1.chr20.20X.pairs.params
 

head -n 10 SVvariants/NA12878_S1.chr20.20X.pairs.histo

R
size_dist <- read.table('SVvariants/NA12878_S1.chr20.20X.pairs.histo')
pdf(file = "SVvariants/fragment.hist.pdf") 
layout(matrix(1:3))
plot(size_dist[,1], size_dist[,2], type='h', main="NA12878 insert size") 
size_dist <- read.table('SVvariants/NA12891_S1.chr20.20X.pairs.histo') 
plot(size_dist[,1], size_dist[,2], type='h', main="NA12891 insert size") 
size_dist <- read.table('SVvariants/NA12892_S1.chr20.20X.pairs.histo') 
plot(size_dist[,1], size_dist[,2], type='h', main="NA12892 insert size") 
dev.off()
quit("no")

#NA12878
delly call -g $REF/hg19.fa -o SVvariants/NA12878.bcf -x $REF/hg19.excl bam/NA12878/NA12878_S1.chr20.20X.pairs.posSorted.bam

#NA12891
delly call -g $REF/hg19.fa -o SVvariants/NA12891.bcf -x $REF/hg19.excl bam/NA12891/NA12891_S1.chr20.20X.pairs.posSorted.bam

#NA12892
delly call -g $REF/hg19.fa -o SVvariants/NA12892.bcf -x $REF/hg19.excl bam/NA12892/NA12892_S1.chr20.20X.pairs.posSorted.bam
  
#bcftools view SVvariants/NA12878.bcf | less -S

delly merge -m 500 -n 1000000 -o SVvariants/sv.bcf -b 500 -r 0.5 SVvariants/NA12878.bcf SVvariants/NA12891.bcf SVvariants/NA12892.bcf

#bcftools view SVvariants/sv.bcf | less -S

#NA12878
delly call -g $REF/hg19.fa -v SVvariants/sv.bcf -o SVvariants/NA12878.geno.bcf -x $REF/hg19.excl bam/NA12878/NA12878_S1.chr20.20X.pairs.posSorted.bam

#NA12891
delly call -g $REF/hg19.fa -v SVvariants/sv.bcf -o SVvariants/NA12891.geno.bcf -x $REF/hg19.excl bam/NA12891/NA12891_S1.chr20.20X.pairs.posSorted.bam

#NA12892
delly call -g $REF/hg19.fa -v SVvariants/sv.bcf -o SVvariants/NA12892.geno.bcf -x $REF/hg19.excl bam/NA12892/NA12892_S1.chr20.20X.pairs.posSorted.bam

#bcftools view SVvariants/NA12878.geno.bcf | less -S

bcftools merge -O b -o SVvariants/merged.bcf SVvariants/NA12878.geno.bcf SVvariants/NA12891.geno.bcf SVvariants/NA12892.geno.bcf
bcftools index SVvariants/merged.bcf
bcftools view SVvariants/merged.bcf > SVvariants/merged.vcf
bgzip -c SVvariants/merged.vcf > SVvariants/merged.vcf.gz
tabix -p vcf SVvariants/merged.vcf.gz


