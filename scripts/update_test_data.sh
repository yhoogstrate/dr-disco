#!/bin/bash

corrected_alignment="tests/fix-chimeric/test_terg_01.filtered.bam" &&

# Header
headerfile="/tmp/header.sam" &&
samtools view -H "$corrected_alignment" > "$headerfile" &&


# ---------------------------- fix chimeric -------------------------- #


# fix-chimeric, test 02
target_alignment="tests/fix-chimeric/test_terg_02.filtered.bam" &&
samtools view "$target_alignment" | cut -f 1 | sort | uniq > "/tmp/readnames.txt" &&
sed -i.bak -e ':a' -e 'N' -e '$!ba' -e 's/\n/|/g' "/tmp/readnames.txt" &&
samtools view "$corrected_alignment" | grep -P $(cat "/tmp/readnames.txt") > "/tmp/body.sam" &&
cat "$headerfile" "/tmp/body.sam" > "/tmp/alignment_new.sam" &&
samtools view -bS "/tmp/alignment_new.sam" > "/tmp/alignment_new.unsorted.bam" &&
samtools sort -o "/tmp/alignment_new.sorted.bam" "/tmp/alignment_new.unsorted.bam" &&
rm "$target_alignment".bai &&
mv "/tmp/alignment_new.sorted.bam" "$target_alignment" &&
samtools index "$target_alignment"


# ------------------------- intron detection ------------------------- #

corrected_alignment="tests/detect-intronic/test_terg_01.filtered.fixed.bam" &&

rm -f "tmp/*" &&
python setup.py install --user ; python tests/test_fix_chimeric_alignment.py &&
rm "$corrected_alignment".bai
mv "tmp/test_terg_01.filtered.fixed.bam" "$corrected_alignment"
samtools index "$corrected_alignment"

# Header
headerfile="/tmp/header.sam" &&
samtools view -H "$corrected_alignment" > "$headerfile" &&


# intron, test 01
target_alignment="tests/detect-intronic/test_terg_01.sub_01.filtered.fixed.bam" &&
samtools view "$target_alignment" | cut -f 1 | sort | uniq > "/tmp/readnames.txt" &&
sed -i.bak -e ':a' -e 'N' -e '$!ba' -e 's/\n/|/g' "/tmp/readnames.txt" &&
samtools view "$corrected_alignment" | grep -P $(cat "/tmp/readnames.txt") > "/tmp/body.sam" &&
cat "$headerfile" "/tmp/body.sam" > "/tmp/alignment_new.sam" &&
samtools view -bS "/tmp/alignment_new.sam" > "/tmp/alignment_new.unsorted.bam" &&
samtools sort -o "/tmp/alignment_new.sorted.bam" "/tmp/alignment_new.unsorted.bam" &&
rm "$target_alignment".bai &&
mv "/tmp/alignment_new.sorted.bam" "$target_alignment" &&
samtools index "$target_alignment" &&


# intron, test 02
target_alignment="tests/detect-intronic/test_terg_01.sub_02.filtered.fixed.bam" &&
samtools view "$target_alignment" | cut -f 1 | sort | uniq > "/tmp/readnames.txt" &&
sed -i.bak -e ':a' -e 'N' -e '$!ba' -e 's/\n/|/g' "/tmp/readnames.txt" &&
samtools view "$corrected_alignment" | grep -P $(cat "/tmp/readnames.txt") > "/tmp/body.sam" &&
cat "$headerfile" "/tmp/body.sam" > "/tmp/alignment_new.sam" &&
samtools view -bS "/tmp/alignment_new.sam" > "/tmp/alignment_new.unsorted.bam" &&
samtools sort -o "/tmp/alignment_new.sorted.bam" "/tmp/alignment_new.unsorted.bam" &&
rm "$target_alignment".bai &&
mv "/tmp/alignment_new.sorted.bam" "$target_alignment" &&
samtools index "$target_alignment" &&


# intron, test 03
target_alignment="tests/detect-intronic/test_terg_01.sub_03.filtered.fixed.bam" &&
samtools view "$target_alignment" | cut -f 1 | sort | uniq > "/tmp/readnames.txt" &&
sed -i.bak -e ':a' -e 'N' -e '$!ba' -e 's/\n/|/g' "/tmp/readnames.txt" &&
samtools view "$corrected_alignment" | grep -P $(cat "/tmp/readnames.txt") > "/tmp/body.sam" &&
cat "$headerfile" "/tmp/body.sam" > "/tmp/alignment_new.sam" &&
samtools view -bS "/tmp/alignment_new.sam" > "/tmp/alignment_new.unsorted.bam" &&
samtools sort -o "/tmp/alignment_new.sorted.bam" "/tmp/alignment_new.unsorted.bam" &&
rm "$target_alignment".bai &&
mv "/tmp/alignment_new.sorted.bam" "$target_alignment" &&
samtools index "$target_alignment" &&


# intron, test 04
target_alignment="tests/detect-intronic/test_terg_01.sub_04.filtered.fixed.bam" &&
samtools view "$target_alignment" | cut -f 1 | sort | uniq > "/tmp/readnames.txt" &&
sed -i.bak -e ':a' -e 'N' -e '$!ba' -e 's/\n/|/g' "/tmp/readnames.txt" &&
samtools view "$corrected_alignment" | grep -P $(cat "/tmp/readnames.txt") > "/tmp/body.sam" &&
cat "$headerfile" "tests/detect-intronic/test_terg_01.sub_04.filtered.splice-offset.sam" "/tmp/body.sam" > "/tmp/alignment_new.sam" &&
samtools view -bS "/tmp/alignment_new.sam" > "/tmp/alignment_new.unsorted.bam" &&
samtools sort -o "/tmp/alignment_new.sorted.bam" "/tmp/alignment_new.unsorted.bam" &&
rm "$target_alignment".bai &&
mv "/tmp/alignment_new.sorted.bam" "$target_alignment" &&
samtools index "$target_alignment" &&


# intron, test 05
target_alignment="tests/detect-intronic/test_terg_01.sub_05.filtered.fixed.bam" &&
samtools view "$target_alignment" | cut -f 1 | sort | uniq > "/tmp/readnames.txt" &&
sed -i.bak -e ':a' -e 'N' -e '$!ba' -e 's/\n/|/g' "/tmp/readnames.txt" &&
samtools view "$corrected_alignment" | grep -P $(cat "/tmp/readnames.txt") > "/tmp/body.sam" &&
cat "$headerfile" "tests/detect-intronic/test_terg_01.sub_05.filtered.small-insert-size.sam" "/tmp/body.sam" > "/tmp/alignment_new.sam" &&
samtools view -bS "/tmp/alignment_new.sam" > "/tmp/alignment_new.unsorted.bam" &&
samtools sort -o "/tmp/alignment_new.sorted.bam" "/tmp/alignment_new.unsorted.bam" &&
rm "$target_alignment".bai &&
mv "/tmp/alignment_new.sorted.bam" "$target_alignment" &&
samtools index "$target_alignment" &


echo ":)"

