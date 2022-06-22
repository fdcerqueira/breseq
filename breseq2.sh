#!/bin/bash

#all fastq files have to have a name like: Inv/Res-whatever-S(same number)_L00*_R*_001.fastq.gz
#in case breseq does not run, due to an error in lib.. from R version, change line 160 (choose R version to install)
#the name of the bacterial references both .fasta and .gbk files need to have the same name.
#directory of the bash script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

#references directory
references=$SCRIPT_DIR/references

#genones parder
download_genomes=$SCRIPT_DIR/midas_parse_taxa.py

#function to show the options in case the number of arguments put is wrong, or the options selected are wrong
function show_help () {
    
    echo "The script it will activate conda and run fastp,flash, bbplist and breseq"
    echo "How to use: $0 parameters"
    echo ""
    echo "OPTIONS:"
    echo ""
    echo " -input <folder with the raw fastq files >"
    echo " -output <location to output the results>"
    echo " -raw_data < name of the folder to put fastq file to change their names>"
    echo " -b <clonal or polymorphism mode in breseq>"
    echo " -r <select the genome reference number shown above (1,2, etc)>"
    echo " -n <select the number of cores to be used>"
    echo " -mode <write (y/n) if you want to run pipeline with merged reads (single-end) or not (pair-end)>"
    echo " -contaminant <write this option to remove contaminant genomes"
    echo " -h <print help>"
    echo ""
    echo " e.g. ./breseq.sh -input /path/to/raw_fastq -output /path/to/ouput_folder -b clonal -r 1 -n 10 -single y"
    echo "e.g. ./breseq.sh -input /path/to/raw_fastq -output /path/to/ouput_folder -b polymorphism -r 1 -n 10 -single n -contaminant"
    echo ""
    
return 0
}

#list all the annotated genomes in the folder
files=$(ls $references/*.gbk)

#loop to create options to be chosen
echo ""
echo "the list of annotated genomes:"
i=1
for j in $files
    do echo "$i.$j" #print option number i and correspondent file (concatenation)
    file[i]=$j # [i] parameter to store value i which will be equal to variable j
    i=$((i+1)) # $(()) arithmetic expansion to evaluate an arithmetic expression and the substitution of the result
done

echo ""
echo ""

#while loop with case to allow to choose parameters in for the downstream programs, print the options selected, or in case of -h to give the show_help function and esit program. the *) means in case of wrong parameter written it will  the program
while [ ! -z "$1" ]; do
    case "$1" in
    -h) show_help ;exit ;;
    -input) shift 1; input=$1;;
    -output) shift 1; output=$1;;
    -raw_data) shift 1; raw_data=$1;; 
    -b) shift 1; MODE=$1;  echo "the breseq mode selected was: ${MODE}"; echo "";;#get variable mode, to insert below in the breseq commands
    -r) shift 1; REFERENCE=$1; filename=$(basename -- "${file[$REFERENCE]}"); filename1="${filename%%.*}"; echo "the reference selected selected was: ${filename1}";printf "\n";; #get variable reference, to insert below in the breseq commands and extract the filename which later allows the breseq to run according the reference chosen
    -n) shift 1; CORES=$1; echo "number of cores selected was: ${CORES}" ;echo "";;
    -mode) shift 1; single=$1;;
    -contaminant) shift 2; contaminant=true;;
    *) echo "ERROR: Incorrect input provided, run the program again"; echo "";show_help; exit;;
    esac
    shift
done

#error messages if the annotated genomes (.gbk), reference fasta files (.fa) are not present
if [ -z "$(ls -1 $references/*.gbk 2>/dev/null )" ]
    then
        echo "There are no reference "*.gbk" files"
        exit
fi

if [ -z "$(ls -1 $references/*.fa 2>/dev/null )" ]
    then
        echo 'There are no files for  "*.fa" files'
        exit
fi

#Error in case variable for breseq mode, reference or numbers of cores are empty
if [ -z $MODE ]
    then
        echo "ERROR: no input for breseq mode"
        echo ""
        show_help
        exit 1
elif [ -z $REFERENCE ]
    then
        echo "ERROR: no input for reference"
        echo ""
        show_help
        exit 1
elif [ -z $CORES ]
    then
        echo "ERROR: no input for number of cores"
        echo ""
        show_help
        exit 1
elif [ -z "$raw_data" ] 
	then
	    echo "ERROR: no input for raw_data folder"	
	    show_help
	    exit 1
elif [ -z "$output" ]
	then
	    echo "ERROR: no output folder defined"
	    show_help
	    exit 1
elif [ -z "$input" ]
	then
	    echo "ERROR: no input folder defined"
	    show_help
	    exit 1
fi

###error in case the number of cores selected for breseq are not appropriated
if [ $CORES -gt 12 ] || [ $CORES -eq 0 ]
then
 echo "ERROR:number of cores too high or 0"
   exit
fi

echo ""
echo ""
###automated workflow to check if the bbplsit, fastp and breseq are not installed in a conda environment
echo "Do you have fastp, bbsplit and breseq installed in a conda env? (y/n)?"
read CONT
echo ""
if [ "$CONT" = "y" ]
then
  echo "Proceeding with pipeline"
else
    #choose to create conda env, prompt for conda env name
    echo ""
    echo "Creating new conda environment, choose name"
    read input
    echo "Name ${input} was chosen"
    
    #list name of packages
    echo "installing base packages"
    conda create -y -n $input
    
    #search for conda to be able to activate it and allow to install the required programs
    eval "$(conda shell.bash hook)"
    conda activate $input
    conda install -y -c bioconda bbmap
    conda install -y -c bioconda fastp=0.20.0
    conda install -y -c bioconda breseq=0.35.1
    conda install -y -c bioconda flash
    pip install biopython
    
    #if there is an error with R in conda, try to change the version installed (more recent)
    conda install -y -c conda-forge r-base=4.1.0
    
    # in case is necessary to install FLASh download from their website: https://ccb.jhu.edu/software/FLASH/
    # then on the folder the tar files is present run to extract and uncompress:
    # tar -xvzf FLASH-x.x.x-Linux-x86_64echo "running breseq"

    #to add to .bashrc path to avoid point to prograam everytime to run this. Put the path were the folder with the program is:
    # echo "export PATH=/path/to/FLASH-x.x.x-Linux-x86_64/:\$PATH" >> ~/.bashrc
    # source ~/.bashrc        
    echo "conda env ${input} created, and the required packages as well"
fi

#activate the conda environment and select the environment you have fastp, bbsplit, breseq
eval "$(conda shell.bash hook)"
echo "the list of your conda environments"
conda env list
echo "Write the one you wanted"
read conda
conda activate $conda

if [ $? -ne 0 ]
    then
        echo "Exiting now"
        exit 1 
    else 
        echo "Conda environment ${conda} activated"
fi

#directories
fastp=$output/fastp
processed_fastq=$output/processed_fastq
merged_reads=$output/merged_reads
breseq=$output/breseq
raw_data=$output/$raw_data
fastp_reports=$fastp/FASTP_REPORTS
bbsplit=$output/bbsplit
merged_reads=$output/merged_reads
contami=$SCRIPT_DIR/contaminant

#references bbsplit
inv=$references/1_invaderBacterialGenome.fa
res=$references/4_residentBacterialGenome.fa
plas1=$references/3_plasmid2_genome_gtgaatcaa.fa
plas2=$references/2_plasmid1_genome_gtggataagt.fa

#function bbpsplit
function bbsplitz () {
    if [ "$contaminant" = "true" ]
        then
            echo "contminanat"  
            cont=$(ls $contami/*.fa* | tr "\n" "," | sed "s/,$//")
            echo "$cont"
            if [ "$REFERENCE" = "1" ] 
                then
                    refs=(${inv},${plas1},${plas2},${res},${cont}) 
                    echo "$refs"
                elif [ "$reference" = "2" ] 
                    then
                        refs=(${res},${plas1},${plas2},${inv},${cont}) 
                        echo "$refs"
            fi 
    elif [ -z "$contaminant" ]
    then
        echo "no contamninant"
            if [ "$REFERENCE" = "1" ]
            then
                echo "running breseq"
                    refs=(${inv},${plas1},${plas2},${res}) 
                    echo "$refs"
            elif [ "$REFERENCE" = "2" ] 
                then
                    refs=(${res},${plas1},${plas2},${inv}) 
                    echo "$refs"
            fi
    fi 

    if [ "$single" = "y" ]
        then
            echo "bbsplit single"
            bbsplit.sh in=$1 ambig2=best ref=$refs basename=$2
    elif [ "$single" = "n" ]
        then
            echo "bbsplit pair end"
            bbsplit.sh in=$1 in2=$2 ref=$refs basename=$3 
    fi
}

#function breseq
function breseqz () {
    
    if [ "$single" = "y" ]
        then
            if [ "$MODE" = "polymorphism" ]
                then
                    breseq -j $1 -r $2 -b 30 -m 20 -p \
                    --polymorphism-frequency-cutoff 0 \
                    --polymorphism-minimum-variant-coverage 5 \
                    --per-position-file \
                    -o $3 \
                    -n $4 $5 
            elif [ "$MODE" = "clonal" ]
                then
                    breseq -j $1 -r $2 -b 30 -m 20 -o $3 -n $4 $5 
            fi
    elif [ "$single" = "n" ]
        then
            if [ "$MODE" = "polymorphism" ]
                then
                    breseq -j $1 -r $2 -b 30 -m 20 -p \
                    --polymorphism-frequency-cutoff 0 \
                    --polymorphism-minimum-variant-coverage 5 \
                    --per-position-file \
                    -o $3 \
                    -n $4 $5 $6 
            elif [ "$MODE" = "clonal" ]
                then
                    breseq -j $1 -r $2 -b 30 -m 20 -o $3 -n $4 $5 $6 
            fi
    fi
}

#in case directory already exists, overwrite it
if [ ! -d "$output" ]
    then
 	    mkdir -p $output
 	    mkdir -p $fastp_reports
 else
        echo "Directory ${output} exist. Do you want overwrite it? (y/n)"
    	read sim
        	if [ $sim = "y" ]
            	then
                    rm -rf $output 
        	        rm -rf $fastp_reports
        	        echo "Creating/overwriting ${output} directory"
        	        mkdir ${output}

    	    else
    	        exit 
        	fi

fi

#in case contaminant folder already exists
if [ "$contaminant" = "true" ] && [ ! -d "$contami" ]
    then 
        mkdir -p $contami
        python $download_genomes $contami
elif [ "$contaminant" = "true" ] && [ -d "$contami" ]
    then
        rm -rf $contami
        echo "Creating/overwriting ${contami} directory"
        mkdir -p $contami
        python $download_genomes $contami
fi

#create folders
mkdir -p $fastp
mkdir -p $processed_fastq
mkdir -p $merged_reads
mkdir -p $fastp_reports
mkdir -p $breseq
mkdir -p $raw_data
mkdir -p $bbsplit
mkdir -p $merged_reads

#check the names of the files begin with Inv or Res
if [ "$REFERENCE" = "1" ]
  then
     for i in $input/*.fastq.gz 
        do
        file=$(basename $i)
  
		    if [[ $file == "Inv"* ]]
		        then
			        echo "copying $file"
			        cp $i $raw_data
		    else
		  	    echo "copying $file, add prefix Inv_"
			    cp $i $raw_data/Inv_${file}
		fi
     done
else 
    for i in $input/*.fastq.gz
       do
		file=$(basename $i)

	        if [[ $file == "Res"* ]]
		        then
		            echo "copying $file"
		            cp $i $raw_data
		    else
		        echo "copying $file"
		        cp $i $raw_data/Res_${file}
		fi
	done		
fi

#check if file names end with R*_001.fastq.gz
for f in $raw_data/*[R_]1.fastq.gz;
    do 
        mv "$f" "$(echo "$f" | sed s/[R_]1/_R1_001/)" 2>/dev/null 
done	 
	 
for f in $raw_data/*[R_]2.fastq.gz;
    do 
        mv "$f" "$(echo "$f" | sed s/[R_]2/_R2_001/)"  2>/dev/null
done

#merge flash output for breseq
#loop for quality control with fastp 
#if else statements in case of single end and pair end
for f1 in $raw_data/*_R1_001.fastq.gz
    do

        file=$(basename $f1)
	    file1=${filename%%_R1_001.fastq.gz}
        f2=${f1%%_R1_001.fastq.gz}_R2_001.fastq.gz
        file2=$(basename $f2)
           
        fastp -q 20 \
        -u 50 \
        --length_required 100 \
        --detect_adapter_for_pe \
        -p 3 -5 -M 20 -W 4 -c \
        -i $f1 \
        -I $f2 \
        -o $fastp/${file} \
        -O $fastp/${file2} \
        -j $fastp_reports/${file1}.json \
        -h $fastp_reports/${file1}.html
done

for i in $(find $fastp -type f -name "*.fastq.gz" \
    | while read F; do basename $F \
    | rev \
    | cut -c 22- \
    | rev; done \
    | sort  \
    | uniq)
        do
          file=${i%%_L00*}

          echo "Merging ${file} R1"
          cat $fastp/${file}_L00*_R1_001.fastq.gz > $merged_reads/${i}_R1_001.fastq.gz
          echo "Merging ${file} R2"
          cat $fastp/${file}_L00*_R2_001.fastq.gz > $merged_reads/${i}_R2_001.fastq.gz
done
	
if [ "$single" = "y" ]
    then
          
        flash=$output/flash
	    mkdir -p $flash
        
        for i in $merged_reads/*_R1_001.fastq.gz
            do
	            file=$(basename $i)	
		        f2=${file%%_R1_001.fastq.gz}_R2_001.fastq.gz
            	prefix=${file%%.fastq.gz}
		
		        flash $merged_reads/${file} $merged_reads/${f2} \
		        --max-overlap 100 \
		        --max-mismatch-density 0 \
		        --output-directory $flash \
		        --compress-prog=gzip \
		        --suffix=gz \
		        -t $CORES \
		        -o $prefix		 
        done

	    for i in $(find $flash -type f -name "*.fastq.gz" \
    	    | while read F; do basename $F \
    	    | rev \
    	    | cut -c 22- \
    	    | rev; done \
    	    | sort  \
    	    | uniq)
            do
                
                file=${i%%.[en]*}
                echo "Merging FLASH ${file} "
                cat $flash/${file}*.fastq.gz > $processed_fastq/${file}_R1_001.fastq.gz	
	done
else
    cp $merged_reads/*fastq.gz $processed_fastq
fi

#run bbsplit
for prefix in $(ls $processed_fastq/*.fastq.gz | sed -E "s/_R[12]_001[.]fastq.gz//" | uniq) #  sed program to separate the extension and then get the unique names of samples ;Extended regexps, makes use of less /
    do
                    
        file=$(basename $prefix)

        if [ "$single" = "y" ]
            then
                bbsplitz ${prefix}_R1_001.fastq.gz $bbsplit/${file%.fastq.gz}%_#.fastq.gz
        elif [ "$single" = "n" ]
            then
                bbsplitz ${prefix}_R1_001.fastq.gz ${prefix}_R2_001.fastq.gz $bbsplit/${file%.fastq.gz}%_#.fastq.gz
            
        fi
done

#run breseq
for i in $(ls $bbsplit/*$filename1*.fastq.gz | sed -E "s/_[12][.]fastq.gz//" | uniq)
    do
        file=$(basename $i)

        echo "inside loop"
        echo "$filename1"
        echo "$file"
        echo "$REFERENCE"
        echo "${file[$REFERENCE]} "

        if [ "$single" = "y" ]
            then 
                breseqz $CORES ${file[$REFERENCE]} $breseq/$file $breseq/${file} ${i}_1.fastq.gz 
        elif [ "$single" = "n" ]
            then
                breseqz $CORES ${file[$REFERENCE]} $breseq/$file $breseq/$file ${i}_1.fastq.gz ${i}_2.fastq.gz 
        fi
done

echo "#########################################################################"
echo "#                                                                       #"
echo "#                     Script finished successfully                      #"
echo "#                                                                       #"
echo "#########################################################################"


