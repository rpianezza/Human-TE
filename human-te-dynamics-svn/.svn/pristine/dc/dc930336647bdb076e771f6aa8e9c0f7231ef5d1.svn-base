#!/bin/zsh
# read parameters
if [ $# -lt 2 ]
   then
   echo "Usage $0 infile outputdir"
   # infile: file with a list of all cram files; each file on new line
   # outputdir: absolute path of outputdir
exit 2
fi

####  INPUT  ####
infile=$1
outputdir=$2
################


### FIXED parameters ####
minmapq=10 # minimum mapping quality is set to 10, remove ambiguous
threads=10
samtoolsmaxmem="3G"
cramref="/Volumes/Temp2/human_TEs/pipeline_test/GRCh38_full_analysis_set_plus_decoy_hla.fa" # reference genome used for generating the cram
ref_fasta="/Volumes/Temp2/human_TEs/pipeline_test/human-te-dynamics/refg/reflibrary_humans_v6.2.fasta"
########################


#### SOFTWARE ####
md5="/sbin/md5"                         			# version BSD; June 6, 2004 (could not find more info)
curl="/usr/local/anaconda3/bin/curl"    			# version curl 7.71.1 (x86_64-apple-darwin13.4.0) libcurl/7.71.1 OpenSSL/1.1.1i zlib/1.2.11 libssh2/1.9.0
ascp="/Users/fschwarz/Applications/AsperaCLI/bin/ascp"
# preparatory work
mkdir -p $outputdir

while read -r line  # read each line of input, store in variable cramid
do
# extract info from file
sampleid=`echo $line|awk '{print $1}'`
crampath=`echo $line |awk '{print $2}'`
md5expect=`echo $line |awk '{print $3}'`

# temp cram-file
cramfile="$outputdir/$sampleid.cram"
tmp=`echo $crampath |perl -pe 's{ftp://ftp.1000genomes.ebi.ac.uk/}{}'`
ascpfile="fasp-g1k@fasp.1000genomes.ebi.ac.uk:$tmp"

command="$ascp --file-checksum=md5 -i /Users/fschwarz/Applications/AsperaCLI/etc/asperaweb_id_dsa.openssh -QTr -l 100m -k1 -P33001 -L- $ascpfile $cramfile"
echo $command
eval $command 

done < $infile






