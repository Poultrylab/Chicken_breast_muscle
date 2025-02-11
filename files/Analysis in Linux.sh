

#1. RNA-Seq alignment and sorting using HISAT2 and Samtools.
#hisat2-2.2.1 samtools-1.11

cd ${cleanPath}
for file in `cat file`
do

$hisat/hisat2 -p 20 -x $ref -1 ${file}_1.fq.gz -2 ${file}_2.fq.gz |$samtools/samtools sort -@ 8 -o $out/${file}.sorted.bam
done


# High-Throughput Transcript Quantification Pipeline Using StringTie v2.1.5: Automated Multi-Threaded Assembly and Expression Analysis
module load StringTie/2.1.5

cd ${cleanPath}
for file in `cat file`
do

stringtie -p 20 -e -G $gtf  -o $out/${file}.gtf  -A $out/${file}.tsv $input/${file}.sorted.bam
done
