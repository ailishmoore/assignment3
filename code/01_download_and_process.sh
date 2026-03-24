#download SRA accession codes
while read run_id; do prefetch $run_id; done < Assignment_3_SRA_Accession_Codes.txt

for sra_file in */*.sra; do
  fastq-dump --split-files --gzip --outdir "$(dirname "$sra_file")" "$sra_file"
done

#trimmomatic
for i in */ ; 
do 	
	
	java -jar PATH_TO_TRIMMOMATIC/Trimmomatic-0.36/trimmomatic-0.36.jar  PE -threads 22 -phred33 -trimlog trimLogFile  $i/*_1.fastq.gz $i/*_2.fastq.gz $i/${i%/}output_forward_paired.fq.gz $i/${i%/}output_forward_unpaired.fq.gz $i/${i%/}output_reverse_paired.fq.gz $i/${i%/}output_reverse_unpaired.fq.gz ILLUMINACLIP:PATH_TO_TRIMMOMATIC/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50; 
	
	rm $i/*_unpaired.fq.gz;
	
done
#Copy paste from here -------------------------------------------------------
# Set input and output directories
BASE_DIR="/media/sorbaralab2/Extra_storage/Ailish/Assignment_3"
OUTPUT_DIR="/media/sorbaralab2/Extra_storage/Ailish/Assignment_3/kraken2_output"
DB_DIR="/media/sorbaralab2/Extra_storage/kraken2/k2_standard_20251015"

mkdir -p "$OUTPUT_DIR"

# Loop over all sample directories
for sample_dir in "$BASE_DIR"/*/; do
    sample_name=$(basename "$sample_dir")

    fwd=$(find "$sample_dir" -maxdepth 1 -type f -iname "*output_forward_paired.fq.gz")
    rev=$(find "$sample_dir" -maxdepth 1 -type f -iname "*output_reverse_paired.fq.gz")

    report_file="${OUTPUT_DIR}/${sample_name}.kraken.report"
    output_file="${OUTPUT_DIR}/${sample_name}.kraken.out"

    # Skip if already processed
    if [[ -f "$report_file" ]]; then
        echo "Skipping $sample_name (already processed)"
        continue
    fi

    if [[ -f "$fwd" && -f "$rev" ]]; then
        echo "Running Kraken2 on $sample_name"

        kraken2 \
          --db "$DB_DIR" \
          --paired "$fwd" "$rev" \
          --use-names \
          --report "$report_file" \
          --output "$output_file"

        echo "Finished $sample_name"
    else
        echo "Missing paired FASTQ files for $sample_name — skipping"
    fi
done

echo "All samples processed with Kraken2."

#To here!-----------------------------------------------------------------------------




