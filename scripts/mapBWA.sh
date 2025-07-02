module load bwa
module load samtools

bwa index $3

bwa aln -l 1024 -n 0.02 -o 2 -t 10 $3 $1 | bwa samse $3 - $1 | samtools view -F 4 -q 25 -@ 10 -uS - | samtools sort -@ 10 -o $2
