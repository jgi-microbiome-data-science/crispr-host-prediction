# crispr-host-prediction

Below is a step-by-step walkthrough of using the CRISPR database to assign host lineages to viral genome sequences

### Clone the repo
```bash
wget https://github.com/jgi-microbiome-data-science/crispr-host-prediction
cd crispr-host-prediction
```

### Download the database
```bash
mkdir global-crispr-db
cd global-crispr-db
wget https://portal.nersc.gov/cfs/m342/global_crispr_db/crispr_spacers_filtered_clustered.tsv
```

### Prepare the database files
```bash
cut -f1,14 crispr_spacers_filtered_clustered.tsv | sed 1d | sed 's/^/>/' | tr '\t' '\n' > spacers.fna
cut -f1,7 crispr_spacers_filtered_clustered.tsv > spacer_lineage.tsv
makeblastdb -in spacers.fna -out spacers.fna -dbtype nucl
cd ..
```

### BLAST viruses to database
```bash
blastn -query test_viruses.fna -db global-crispr-db/spacers.fna -dust no -word_size 8 -max_target_seqs 1000 -outfmt '6 std qlen slen' -num_threads 64 > blastn.tsv
```

### Assign GTDB host lineages
```bash
assign_host.py -i blastn.tsv -d global-crispr-db -o host_predictions.tsv
```


