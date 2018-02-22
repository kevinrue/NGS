# https://www.gnu.org/software/parallel/parallel_tutorial.html

# echo list of files in current directory
ls . | parallel echo :::: - 
