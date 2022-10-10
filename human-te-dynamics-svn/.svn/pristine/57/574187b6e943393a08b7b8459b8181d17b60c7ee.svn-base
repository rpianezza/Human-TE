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
samtoolsmaxmem="10G"
cramref="/Volumes/Temp2/human_TEs/pipeline_test/GRCh38_full_analysis_set_plus_decoy_hla.fa" # reference genome used for generating the cram
ref_fasta="/Volumes/Temp2/human_TEs/refseqs/reflibrary_humans_v5.fasta"
########################



#### SOFTWARE ####
md5="/sbin/md5"                         			# version BSD; June 6, 2004 (could not find more info)
curl="/usr/local/anaconda3/bin/curl"    			# version curl 7.71.1 (x86_64-apple-darwin13.4.0) libcurl/7.71.1 OpenSSL/1.1.1i zlib/1.2.11 libssh2/1.9.0
samtools="/usr/local/bin/samtools"				# 1.10 (using htslib 1.10.2)				
bwa="/usr/local/bin/bwa"          	     			# version 0.7.17-r1188
mosdepth="/usr/local/anaconda3/bin/mosdepth"			# version 0.3.1
python="/usr/local/bin/python3"					# version 3.8.5
mapstat="/Volumes/Temp2/human_TEs/pipeline_test/human-te-dynamics/humante-mapstat.py"

# preparatory work
mkdir -p $outputdir
bwa index $ref_fasta
samtools faidx $ref_fasta
ref_fai="$ref_fasta.fai"

while read -r line  # read each line of input, store in variable cramid
do
# extract info from file
sampleid=`echo $line|cut -f1`
crampath=`echo $line |cut -f2`
md5expect=`echo $line |cut -f3`

# temp cram-file
cramfile="$outputdir/$sampleid.cram"

# CURL and md5
command="$curl -# $crampath --output $cramfile"
eval $command       # i exchange echo and eval for debugging

md5observe=`md5 -q $cramfile`    # -q quite mode, to only print the md5 checksum
echo "md5check $sampleid $md5observe $md5expect"  # allows to check whether the md5 failed for some

# 
sortbam="$outputdir/$sampleid.sort.bam"
prefix="$outputdir/$sampleid"

# samtools fastq -@10 -N -o - HGDP00001-Brahui.cram |head -1000000| bwa mem -t 10 /Volumes/Temp2/human_TEs/refseqs/reflibrary_humans_v5.fasta - | samtools view -@10 -b | samtools sort -m 4G -O bam -o test.sort.bam -@10 -T quatsch 
command="$samtools fastq -@$threads -N -o - $cramfile | $bwa mem -t $threads $ref_fasta - | $samtools view -F 0x04 -@$threads -b | samtools sort -@$threads -O bam -o $sortbam -m $samtoolsmaxmem -T $prefix"
eval $command
$samtools index $sortbam

# Mosdepth
mdoutput="$outputdir/$sampleid"
command="$mosdepth $mdoutput $sortbam"
eval $command

# Mapstat
# samtools view HGDP00001-Brahui.sort.bam |python ../human-te-dynamics/humante-mapstat.py --fai /Volumes/Temp2/human_TEs/refseqs/reflibrary_humans_v5.fasta.fai --sam - --min-mq 10
mapstatoutput="$outputdir/$sampleid.mapstat"
command="$samtools view $sortbam | $python $mapstat --sam - --fai $ref_fai --min-mq $minmapq > $mapstatoutput"
eval $command

# cleanup
# rm $cramflie

done < $infile






