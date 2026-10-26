sample=PTD858-1_L1
index=/NAS/home/yszhou/data/ensemble-104
hisat2 -p 5   -q -x   $index/Homo_sapiens.GRCh38.dna.primary_assembly.fa  -1   ./Rawdata/$sample.new\_1.fq.gz  -2   ./Rawdata/$sample.new\_2.fq.gz  -S  $sample.sam
samtools view -b  -q 30  -Su  $sample.sam | samtools sort  -@ 5 - > $sample.sorted.bam
python2 /NAS/home/yszhou/software/REDItools-1.2.1/REDItools-install/bin/REDItoolDenovo.py  -o  $sample.REDItools   -i  $sample.sorted.bam  -f   $index/Homo_sapiens.GRCh38.dna.primary_assembly.fa   -t  5  -e -d -l -U [AG,TC,CT,GA] -p -u -m 60 -T 5-5 -W  -n  0.0
