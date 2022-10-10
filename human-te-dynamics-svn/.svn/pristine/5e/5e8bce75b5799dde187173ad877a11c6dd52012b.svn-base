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
ascp="/Users/fschwarz/Applications/AsperaCLI/bin/ascp"
samtools="/usr/local/bin/samtools"				# 1.10 (using htslib 1.10.2)				
bwa="/usr/local/bin/bwa"          	     			# version 0.7.17-r1188
mosdepth="/usr/local/anaconda3/bin/mosdepth"			# version 0.3.1
python="/usr/local/anaconda3/bin/python3"					# version 3.8.5
java="/usr/bin/java"						# java version "1.8.0_192"
gzip="/usr/bin/gzip"						# Apple gzip 251
mptsync="/Volumes/Temp2/human_TEs/pipeline_test/human-te-dynamics/scripts/mpileup2sync.jar"
mapstat="/Volumes/Temp2/human_TEs/pipeline_test/human-te-dynamics/humante-mapstat-weight.py"

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
tmp=`echo $crampath |perl -pe 's{ftp://ftp.1000genomes.ebi.ac.uk/}{}'`
ascpfile="fasp-g1k@fasp.1000genomes.ebi.ac.uk:$tmp"

command="$ascp --file-checksum=md5 -i /Users/fschwarz/Applications/AsperaCLI/etc/asperaweb_id_dsa.openssh -QTr -l 100m -k1 -P33001 -L- $ascpfile $cramfile"
eval $command

# map and sort
sortbam="$outputdir/$sampleid.sort.bam"
prefix="$outputdir/$sampleid"

command="$samtools fastq -@$threads -N -o - $cramfile | $bwa mem -t $threads $ref_fasta - | $samtools view -F 0x04 -@$threads -b | $samtools sort -@$threads -m $samtoolsmaxmem -T $prefix | $samtools view -b -o $sortbam --write-index"
eval $command

# Mosdepth
mdoutput="$outputdir/$sampleid"
command="$mosdepth $mdoutput $sortbam"
eval $command

# Mapstat
mapstatoutput="$outputdir/$sampleid.mq${minmapq}.mapstat"
command="$samtools view $sortbam | $python $mapstat --sam - --fai $ref_fai --min-mq $minmapq > $mapstatoutput &"
eval $command

mapstatoutputmq="$outputdir/$sampleid.mq0.mapstat"
command2="$samtools view $sortbam | $python $mapstat --sam - --fai $ref_fai --min-mq 0 > $mapstatoutputmq &"
eval $command2

# sync
syncfile="$outputdir/$sampleid.mq${minmapq}.sync.gz"
command="$samtools mpileup -q $minmapq -B --min-BQ 0 -d 100000 -A -f $ref_fasta $sortbam | $java -jar $mptsync --threads $threads --input /dev/stdin --output /dev/stdout | $gzip -c > $syncfile &"
eval $command

# cleanup
# rm $cramflie
done < $infile






