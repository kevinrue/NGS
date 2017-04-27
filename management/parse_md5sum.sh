
# With `path_to_md5sum` a file with records:
# c4f9d24f1052da9e224247d527a61c98  WTCHG_310996_249_1.fastq.gz

# Extract full list of *.fastq.gz files
grep \.fastq\.gz$ "$inputfile" | cut -d' ' -f3

# Extract the list of unique basenames
grep \.fastq\.gz$ "$inputfile" | cut -d' ' -f3 | sed 's/_..fastq.gz//g' | uniq


