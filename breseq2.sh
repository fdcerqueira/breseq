#!/bin/bash
###all fastq files have to have a name like: Inv/Res-whatever-S(same number)_L00*_R*_001.fastq.gz
##in case breseq does not run, due to an error in lib.. from R version, change line 160 (choose R version to install)
#the name of the bacterial references both .fasta and .gbk files need to have the same name.
###directory of the bash script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
###function to show the options in case the number of arguments put is wrong, or the options selected are wrong
function show_help () {
    
    printf "The script it will activate conda and run fastp, bbplist and breseq\n"
    printf "How to use: $0 parameters\n"
    printf "\n"
    printf "OPTIONS:\n"
    printf "for  breseq\n"
    printf " -b <clonal or polymorphism mode in breseq>\n"
    printf " -r <select the genome reference number shown above (1,2, etc)>\n"
    printf " -n <select the number of cores to be used>\n"
    printf " -h <print help>\n"
    printf "\n"
    printf " e.g.  bash breseq.sh -b clonal -r 1 -n 12 \n"
    printf "\n"
    
return 0
}
###list all the annotated genomes in the folder
files=$(ls *.gbk)
i=1
###loop to create options to be chosen
printf " \n"
echo "the list of annotated genomes:"
for j in $files
    do echo "$i.$j" #print option number i and correspondent file (concatenation)
    file[i]=$j # [i] parameter to store value i which will be equal to variable j
    i=$((i+1)) # $(()) arithmetic expansion to evaluate an arithmetic expression and the substitution of the result
done

printf "\n"

printf "\n"
###while loop with case to allow to choose parameters in for the downstream programs, print the options selected, or in case of -h to give the show_help function and esit program. the *) means in case of wrong parameter written it will exit the program
while [ ! -z "$1" ]; do
    case "$1" in
    -h)
    show_help
    exit
    ;;
    -b)
    shift 1
    ##get variable mode, to insert below in the breseq commands
    MODE="$1"
    echo "the breseq mode selected was: ${MODE}"
    printf "\n"
    ;;
    -r)
    shift 1
    ##get variable reference, to insert below in the breseq commands
    ### extract the filename which later allows the breseq to run according the reference chosen
    REFERENCE="$1"
    filename=$(basename -- "${file[$REFERENCE]}")
    filename1="${filename%%.*}"
    echo "the reference selected selected was: ${filename1}"
    printf "\n"
    ;;
    -n)
    shift 1
    CORES="$1"
    echo "number of cores selected was: ${CORES}"
    printf "\n"
    ;;
    -mode) shift 2; single=true;;
    *)
    echo "ERROR: Incorrect input provided, run the program again"
    printf "\n"
    show_help
    exit
    ;;
esac
shift
done

###error messages if the annotated genomes (.gbk), reference fasta files (.fa) are not present
if [[ -z "$(ls -1 *.gbk 2>/dev/null )" ]]
 then
    echo "There are no reference "*.gbk" files"
    exit
fi

if [[ -z "$(ls -1 *.fa 2>/dev/null )" ]]
 then
    echo 'There are no files for bbsplit "*.fa" files'
    exit
fi
###error in case fastq files  are not present according the reference chosen
if [[ -z "$(ls -1 Res*.fastq.gz 2>/dev/null )" ]] && [[ $REFERENCE == 2 ]]
then
    echo "There are no Resident E.coli fastq files to process"
    exit
fi
###error in case fastq files  are not present according the reference chosen
if [[ -z "$(ls -1 Inv*.fastq.gz 2>/dev/null )" ]] && [[ $REFERENCE == 1 ]]
then
    echo "There are no Invader E.coli fastq files to process"
    exit
fi


###Error in case variable for breseq mode, reference or numbers of cores are empty
if [ -z $MODE ]
    then
        echo "ERROR: no input for breseq mode"
        printf "\n"
        show_help
        exit
elif [ -z $REFERENCE ]
    then
        echo "ERROR: no input for reference"
        printf "\n"
        show_help
        exit
elif [ -z $CORES ]
    then
        echo "ERROR: no input for number of cores"
        printf "\n"
        show_help
        exit
fi

###error in case the number of cores selected for breseq are not appropriated
if [[ $CORES -gt 12 ]] || [[ $CORES -eq 0 ]]
then
 echo "ERROR:number of cores too high or 0"
   exit
fi

printf "\n"
printf "\n"
###automated workflow to check if the bbplsit, fastp and breseq are not installed in a conda environment
echo "Do you have fastp, bbsplit and breseq installed in a conda env? (y/n)?"
read CONT
printf "\n"
if [[ $CONT == y ]]
then
  echo "Proceeding with fastp";
else
###choose to create conda env
### prompt for conda env name
printf "\n"
  echo "Creating new conda environment, choose name"
  read input
  echo "Name ${input} was chosen";
    #list name of packages
        echo "installing base packages"
        conda create --name $input
#search for conda to be able to activate it and allow to install the required programs
        eval "$(conda shell.bash hook)"
        conda activate $input
        conda install -c bioconda/label/cf201901 bbmap
        conda install -c bioconda/label/cf201901 fastp
        conda install -c bioconda/label/cf201901 breseq
#if there is an error with R in conda, try to change the version installed (more recent)
        conda install -c conda-forge r-base=4.1.0
  echo "conda env ${input} created, and the required packages as well"
fi

###activate the conda environment and select the environment you have fastp, bbsplit, breseq
eval "$(conda shell.bash hook)"
echo "the list of your conda environments"
conda env list
echo "Write the one you wanted"
read conda
conda activate $conda
echo "Conda environment ${conda} activated"


PASTA=${SCRIPT_DIR}/good_ones/
####in case directory already exists, overwrite it
if [ ! -d "$PASTA" ]
 then
 	mkdir ${PASTA}
 	mkdir ${PASTA}/FASTP_REPORTS
 else
    	echo "Directory ${PASTA} exist. Do you want overwrite it? (y/n)"
    	read sim
        	if [ $sim = "n" ]
            	then
            	exit
    	else
        	rm -rf $PASTA
        	echo "Creating/overwriting ${PASTA} directory"
        	mkdir ${PASTA}
        	mkdir ${PASTA}/FASTP_REPORTS
        fi

fi
#loop for quality control with fastp 
#if else statements in case of single end and pair end

if [ "$single" = "true" ]
     then
        for f1 in *merged*.fastq.gz
            do 
            fastp -q 20 \
            -u 50 \
            --length_required 100 \
            -p 3 -5 -M 20 -W 4 \
            -i $f1 \
            -o ./good_ones/t-"$f1" \
            -h ${PASTA}/FASTP_REPORTS/${f1}.html
        done
    else
        for f1 in *_R1_001.fastq.gz
            do
            f2=${f1%%_R1_001.fastq.gz}_R2_001.fastq.gz
            fastp -q 20 \
            -u 50 \
            --length_required 100 \
            --detect_adapter_for_pe \
            -p 3 -5 -M 20 -W 4 -c \
            -i $f1 \
            -I $f2 \
            -o ./good_ones/t-"$f1" -O ./good_ones/t-"$f2" \
            -j ${PASTA}/FASTP_REPORTS/${f1}.json \
            -h ${PASTA}/FASTP_REPORTS/${f1}.html
        done
fi

cd good_ones
#for loop to merge fastq from the same sample that are from different lanes.
if [ "$single" != "true" ]
    then

        for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)
            do echo "Merging R1 $i"
                cat "$i"_L00*_R1_001.fastq.gz > M_"$i"_L001_R1_001.fastq.gz
                echo "Merging R2 $i"
                cat "$i"_L00*_R2_001.fastq.gz > M_"$i"_L001_R2_001.fastq.gz
        done
fi

#remove original files not merged in the good_ones folder, copy the fasta files for decontamination, copy the reference input and make a directory "breseq"
cp ../*.fa ./
cp ../${file[$REFERENCE]} ./
mkdir breseq
if [ "$single" != "true" ]
	then
		rm t-*.fastq.gz
fi
##loop to run bbpsplit
if [ "$single" != "true" ]
    then 
        if [ "$REFERENCE" = "1" ]
            then

                for prefix in $(ls *.fastq.gz | sed -E "s/_R[12]_001[.]fastq.gz//" | uniq) #  sed program to separate the extension and then get the unique names of samples ;Extended regexps, makes use of less /
                    do
                    bbsplit.sh in1=${prefix}_R1_001.fastq.gz \
                    in2=${prefix}_R2_001.fastq.gz \
                    ambig2=best\
                    ref=1_invaderBacterialGenome.fa,2_plasmid1_genome_gtggataagt.fa,3_plasmid2_genome_gtgaatcaa.fa,4_residentBacterialGenome.fa \
                    basename=./breseq/${prefix%.fastq.gz}%_#.fastq.gz
                done
        else
                for prefix in $(ls *.fastq.gz | sed -E "s/_R[12]_001[.]fastq.gz//" | uniq) #  sed program to separate the extension and then get the unique names of samples ;Extended regexps, makes use of less /
                    do
                    bbsplit.sh in1=${prefix}_R1_001.fastq.gz \
                    in2=${prefix}_R2_001.fastq.gz \
                    ambig2=best\
                    ref=4_residentBacterialGenome.fa,2_plasmid1_genome_gtggataagt.fa,3_plasmid2_genome_gtgaatcaa.fa,1_invaderBacterialGenome.fa \
                    basename=./breseq/${prefix%.fastq.gz}%_#.fastq.gz
                done
        fi
fi

if [ "$single" = "true" ]
    then
        if [ "$REFERENCE" = "1" ]
            then

                for prefix in $(ls *.fastq.gz | sed -E "s/_R1_001[.]fastq.gz//" | uniq) #  sed program to separate the extension and then get the unique names of samples ;Extended regexps, makes use of less /
                    do
                    bbsplit.sh in=${prefix}_R1_001.fastq.gz \
                    ambig2=best\
                    ref=1_invaderBacterialGenome.fa,2_plasmid1_genome_gtggataagt.fa,3_plasmid2_genome_gtgaatcaa.fa,4_residentBacterialGenome.fa \
                    basename=./breseq/${prefix%.fastq.gz}%_#.fastq.gz
                done
        else
                for prefix in $(ls *.fastq.gz | sed -E "s/_R1_001[.]fastq.gz//" | uniq) #  sed program to separate the extension and then get the unique names of samples ;Extended regexps, makes use of less /
                    do
                    bbsplit.sh in=${prefix}_R1_001.fastq.gz \
                    ambig2=best\
                    ref=4_residentBacterialGenome.fa,2_plasmid1_genome_gtggataagt.fa,3_plasmid2_genome_gtgaatcaa.fa,1_invaderBacterialGenome.fa \
                    basename=./breseq/${prefix%.fastq.gz}%_#.fastq.gz
                done
        fi
fi


cd breseq
cp ../${file[$REFERENCE]} ./

#loop accoring the reference and the files: inv fastq files to invader reference and res fastq files to resident reference

if [ "$single" != "true" ]
    then
        if [ $filename1 = "invaderBacterialGenome" ] && [ $MODE = "polymorphism" ]
            then
                for i in $(ls *Inv*$filename1*.fastq.gz | sed -E "s/_[12][.]fastq.gz//" | uniq); do
                    breseq -j ${CORES} \
                    -r ${file[$REFERENCE]} \
                    -b 30 -m 20 -p \
                    --polymorphism-frequency-cutoff 0 \
                    --polymorphism-minimum-variant-coverage 5 \
                    --polymorphism-minimum-variant-coverage-each-strand 2 \
                    --per-position-file \
                    -o ${i}_breseq \
                    -n ${i}_breseq \
                    ${i}_1.fastq.gz \
                    ${i}_2.fastq.gz
                done
        elif [ $filename1 = "residentBacterialGenome" ] && [ $MODE = "polymorphism" ]
            then
                for j in $(ls *Res*$filename1*.fastq.gz | sed -E "s/_[12][.]fastq.gz//" | uniq); do
                    breseq -j ${CORES} \
                    -r ${file[$REFERENCE]} \
                    -b 30 -m 20 -p \
                    --polymorphism-frequency-cutoff 0 \
                    --polymorphism-minimum-variant-coverage 5 \
                    --polymorphism-minimum-variant-coverage-each-strand 2 \
                    --per-position-file \
                    -o ${j}_breseq \
                    -n ${j}_breseq \
                    ${j}_1.fastq.gz \
                    ${j}_2.fastq.gz
                done
        elif [ $filename1 = "invaderBacterialGenome" ] && [ $MODE = "clonal" ]
            then
                for i in $(ls *Inv*$filename1*.fastq.gz | sed -E "s/_[12][.]fastq.gz//" | uniq); do
                    breseq -j ${CORES} \
                    -r ${file[$REFERENCE]} \
                    -b 30 -m 20 \
                    --polymorphism-frequency-cutoff 0 \
                    --polymorphism-minimum-variant-coverage 5 \
                    --polymorphism-minimum-variant-coverage-each-strand 2 \
                    --per-position-file \
                    -o ${i}_breseq \
                    -n ${i}_breseq \
                    ${i}_1.fastq.gz \
                    ${i}_2.fastq.gz
                done
        elif [ $filename1 = "residentBacterialGenome" ] && [ $MODE = "clonal" ]
            then
                for j in $(ls *Res*$filename1*.fastq.gz | sed -E "s/_[12][.]fastq.gz//" | uniq); do
                    breseq -j ${CORES} \
                    -r ${file[$REFERENCE]} \
                    -b 30 -m 20 \
                    --polymorphism-frequency-cutoff 0 \
                    --polymorphism-minimum-variant-coverage 5 \
                    --polymorphism-minimum-variant-coverage-each-strand 2 \
                    --per-position-file \
                    -o ${j}_breseq \
                    -n ${j}_breseq \
                    ${j}_1.fastq.gz \
                    ${j}_2.fastq.gz
                done
        fi
fi

if [ "$single" = "true" ]
    then
        if [ $filename1 = "invaderBacterialGenome" ] && [ $MODE = "polymorphism" ]
            then
                for i in $(ls *Inv*$filename1*.fastq.gz | sed -E "s/_[1][.]fastq.gz//" | uniq); do
                    breseq -j ${CORES} \
                    -r ${file[$REFERENCE]} \
                    -b 30 -m 20 -p \
                    --polymorphism-frequency-cutoff 0 \
                    --polymorphism-minimum-variant-coverage 5 \
                    --per-position-file \
                    -o ${i}_breseq \
                    -n ${i}_breseq \
                    ${i}_1.fastq.gz 
                done
        elif [ $filename1 = "residentBacterialGenome" ] && [ $MODE = "polymorphism" ]
            then
                for j in $(ls *Res*$filename1*.fastq.gz | sed -E "s/_[1][.]fastq.gz//" | uniq); do
                    breseq -j ${CORES} \
                    -r ${file[$REFERENCE]} \
                    -b 30 -m 20 -p \
                    --polymorphism-frequency-cutoff 0 \
                    --polymorphism-minimum-variant-coverage 5 \
                    --per-position-file \
                    -o ${j}_breseq \
                    -n ${j}_breseq \
                    ${j}_1.fastq.gz 
                done
        elif [ $filename1 = "invaderBacterialGenome" ] && [ $MODE = "clonal" ]
            then
                for i in $(ls *Inv*$filename1*.fastq.gz | sed -E "s/_[1][.]fastq.gz//" | uniq); do
                    breseq -j ${CORES} \
                    -r ${file[$REFERENCE]} \
                    -b 30 -m 20 \
                    --polymorphism-frequency-cutoff 0 \
                    --polymorphism-minimum-variant-coverage 5 \
                    --per-position-file \
                    -o ${i}_breseq \
                    -n ${i}_breseq \
                    ${i}_1.fastq.gz 
                done
        elif [ $filename1 = "residentBacterialGenome" ] && [ $MODE = "clonal" ]
            then
                for j in $(ls *Res*$filename1*.fastq.gz | sed -E "s/_[1][.]fastq.gz//" | uniq); do
                    breseq -j ${CORES} \
                    -r ${file[$REFERENCE]} \
                    -b 30 -m 20 \
                    --polymorphism-frequency-cutoff 0 \
                    --polymorphism-minimum-variant-coverage 5 \
                    --per-position-file \
                    -o ${j}_breseq \
                    -n ${j}_breseq \
                    ${j}_1.fastq.gz 
                done
        fi
fi

#--polymorphism-minimum-variant-coverage-each-strand 2 \
