conda —version

CONDA_SUBDIR=osx-64 conda create -n bracken_env -c bioconda -c conda-forge bracken kraken2

y

conda activate bracken_env
conda config --env --set subdir osx-64

for file in /Users/ailishm/Desktop/Kraken2_Ailish/*.kraken.report; do
  sample=$(basename "$file" .kraken.report)
  echo "Processing $sample..."
  bracken \
    -d /Users/ailishm/Desktop/Kraken2_Ailish \
    -i "$file" \
    -o /Users/ailishm/Desktop/Kraken2_Ailish/${sample}.bracken \
    -w /Users/ailishm/Desktop/Kraken2_Ailish/${sample}.bracken.report \
    -r 150 \
    -l S
  echo "Done with $sample!"
done
echo "All samples complete!"

