#/bin/bash

bam=$1
txt=eQTL_variants.txt

# Index bam file if it does not exists
if ! test -f ${bam}.bai; then
  echo "Index $bam file"
  samtools index $bam
fi

echo "sample,variant,T,C,A,G" > ${bam%.bam}.csv

for var in $(cat $txt)
do
  position=$(echo $var | cut -d"_" -f 1)
  pos=$(echo $position | cut -d":" -f 2)
  position=${position}-${pos}

  echo $position

  # get all nucleotides
  nucleotides=$(samtools mpileup -r $position $bam | cut -f 5)

  # make all nucleotides capital
  nucleotides=${nucleotides^^}

  # count nucleotides
  T=$(echo $nucleotides | grep -o "T" | wc -l)
  C=$(echo $nucleotides | grep -o "C" | wc -l)
  A=$(echo $nucleotides | grep -o "A" | wc -l)
  G=$(echo $nucleotides | grep -o "G" | wc -l)

  # save nucleotide count as CSV file

  echo "$bam,$var,$T,$C,$A,$G" >> ${bam%.bam}.csv
done