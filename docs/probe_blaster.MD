# Probe Blaster

Probe blaster accepts a `.fasta` file containing genome probes, and runs them via `tblastn` against the 
provided genome database.

It will then provide a `.tsv` file on a per-fasta-record basis, output to the directory provided 
(defaulting to `<current working directory>/OUTPUT`)

## Usage

#### Arguments

| Argument | Verbose | Description | Type | Required | Default |
|:--------:|:-------:| ----------- |:----:|:--------:|:-------:|
| -f       | --file  | Source fasta file containing the sample probe records to run through tblastn | Filepath | True | N/A |
| -db      | --database_path | Path to the genome database against which the probe records are to be BLASTed | str |True | N/A |
| -o       | --output_dir | Location to which to write the result files | str | False | _current_working_directory_/OUTPUT |
| -a       | --alignment_length_threshold | Minimum length threshold that BLAST result alignment sequence lengths should exceed | int | False | 400 | 

#### Example

``
python src/probe_blaster.py -f data/fasta_file.fasta -db dbs/genome_db -o results/probe_blaster_output
``

#### Output

_NB: Column headers are not present in the data output, they are provided here for illustrative purposes only_

| Accession ID | Scaffold | Scaffold length | Start | End | E-value | Alignment length | Acc. Seq | Scaffold Alignment | Frame |
| ------------ | -------- | --------------- | ----- | --- | ------- | ---------------- | -------- | ------------------ | ----- |
| XXX11111.1 | Scaf_00000111 | 12345678 | 992016 | 994340 | 1.63E-114 | 797 | MVLDLKDCFFSIPLAEQDRE... | LILDLKDCFFTIPLHPDD... | 1 |
| ---------- | ------------- | -------- | ------ | ------ | --------- | --- | ----------------------- | --------------------- | - |
| ---------- | ------------- | -------- | ------ | ------ | --------- | --- | ----------------------- | --------------------- | - |
| ---------- | ------------- | -------- | ------ | ------ | --------- | --- | ----------------------- | --------------------- | - |