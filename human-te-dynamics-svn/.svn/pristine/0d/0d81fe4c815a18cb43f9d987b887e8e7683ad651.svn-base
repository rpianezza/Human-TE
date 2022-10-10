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
threads=20
samtoolsmaxmem="3G"
cramref="/Volumes/Temp2/human_TEs/pipeline_test/GRCh38_full_analysis_set_plus_decoy_hla.fa" # reference genome used for generating the cram
ref_fasta="/Volumes/Temp2/human_TEs/pipeline_test/human-te-dynamics/refg/reflibrary_humans_v6.2.fasta"
########################


#### SOFTWARE ####
md5="/sbin/md5"                         			# version BSD; June 6, 2004 (could not find more info)
curl="/usr/local/anaconda3/bin/curl"    			# version curl 7.71.1 (x86_64-apple-darwin13.4.0) libcurl/7.71.1 OpenSSL/1.1.1i zlib/1.2.11 libssh2/1.9.0
samtools="/usr/local/bin/samtools"				# 1.10 (using htslib 1.10.2)				
bwa="/usr/local/bin/bwa"          	     			# version 0.7.17-r1188
mosdepth="/usr/local/anaconda3/bin/mosdepth"			# version 0.3.1
python="/usr/local/anaconda3/bin/python3"					# version 3.8.5
java="/usr/bin/java"						# java version "1.8.0_192"
gzip="/usr/bin/gzip"						# Apple gzip 251
mptsync="/Volumes/Temp2/human_TEs/pipeline_test/human-te-dynamics/scripts/mpileup2sync.jar"
mapstat="/Volumes/Temp2/human_TEs/pipeline_test/human-te-dynamics/scripts/humante-mapstat-weight.py"

# preparatory work
mkdir -p $outputdir
bwa index $ref_fasta
samtools faidx $ref_fasta
ref_fai="$ref_fasta.fai"

while read -r line  # read each line of input, store in variable cramid
do
# extract info from file
sampleid=`echo $line|awk '{print $1}'`
crampath=`echo $line |awk '{print $2}'`
md5expect=`echo $line |awk '{print $3}'`

# temp cram-file
cramfile="$outputdir/$sampleid.cram"


# 
sortbam="$outputdir/$sampleid.sort.bam"
prefix="$outputdir/$sampleid"

# samtools fastq -@2 -N -o - HGDP00001-Brahui.cram |head -1000000| bwa mem -t 10 /Volumes/Temp2/human_TEs/refseqs/reflibrary_humans_v5.fasta - | samtools view -@10 -b | samtools sort -m 4G -O bam -o test.sort.bam -@10 -T quatsch 
# bwa mem -t 10 /Volumes/Temp2/human_TEs/refseqs/reflibrary_humans_v5.1.fasta test.fastq.gz | samtools view -F 0x04 -@10 -b |samtools sort -@10 -m10G | samtools view -b -o test2.bam --write-index
command="$samtools fastq -@2 -N -o - $cramfile | $bwa mem -t $threads $ref_fasta - | $samtools view -F 0x04 -@2 -b | $samtools sort -@2 -m $samtoolsmaxmem -T $prefix | $samtools view -b -o $sortbam --write-index"
eval $command

# Mosdepth
mdoutput="$outputdir/$sampleid"
command="$mosdepth $mdoutput $sortbam"
eval $command

# Mapstat
# samtools view HGDP00001-Brahui.sort.bam |python ../human-te-dynamics/humante-mapstat.py --fai /Volumes/Temp2/human_TEs/refseqs/reflibrary_humans_v5.fasta.fai --sam - --min-mq 10
mapstatoutput="$outputdir/$sampleid.mq${minmapq}.mapstat"
command="$samtools view $sortbam | $python $mapstat --sam - --fai $ref_fai --min-mq $minmapq > $mapstatoutput &"
eval $command

mapstatoutputmq="$outputdir/$sampleid.mq0.mapstat"
command2="$samtools view $sortbam | $python $mapstat --sam - --fai $ref_fai --min-mq 0 > $mapstatoutputmq &"
eval $command2

# sync
# samtools mpileup map_cov35_rl150_uniform_error0.01_to_reflib_v6.sort.bam -B --min-BQ 0 -A -d 100000 -f reflibrary_humans_v6.fasta | java -jar ../human-te-dynamics/mpileup2sync.jar --input /dev/stdin --output /dev/stdout | gzip -c > test3.sync.gz
syncfile="$outputdir/$sampleid.mq${minmapq}.sync.gz"
command="$samtools mpileup -q $minmapq -B --min-BQ 0 -d 100000 -A -f $ref_fasta $sortbam | $java -jar $mptsync --threads $threads --input /dev/stdin --output /dev/stdout | $gzip -c > $syncfile &"
eval $command

# cleanup
# rm $cramflie

done < $infile






