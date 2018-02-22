#!/bin/bash

# # root
# ftp://ftp.ncbi.nlm.nih.gov/geo/
# # range subdirectory: nnn = last three 
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE12nnn/

cd "/gfs/work/kralbrecht/04_lincRNAs/data/benoist"
pwd

# GSE to download
gse_list=(GSE12345 GSE67890) # dummy examples

for gse_id in ${gse_list[@]}
do
        echo "gse_id: ${gse_id}"
        range_sub=${gse_id:0:$((${#gse_id}-3))}
        echo "range_sub: ${range_sub}"
        file_basename="${gse_id}_RAW.tar"
        ftp_link="ftp://ftp.ncbi.nlm.nih.gov/geo/series/${range_sub}nnn/${gse_id}/suppl/${file_basename}"
        echo "ftp_link: ${ftp_link}"
        cmd="wget -nv ${ftp_link}"
        echo "cmd: ${cmd}"
        #eval $cmd                                                                                                                                                                             
        cmd="mkdir -v ${gse_id}"
        echo "cmd: ${cmd}"
        eval $cmd
        cmd="tar -xvC ${gse_id} -f ${file_basename}"
        echo "cmd: ${cmd}"
        eval $cmd
done

echo "Complete"
