<link href="formatting/style.css" rel="stylesheet"></link>

# Probe Blaster

Probe blaster accepts a `.fasta` file containing genome probes, and runs them via `tblastn` against the 
provided genome database.

It will then provide a `.tsv` file on a per-fasta-record basis, output to the directory provided 
(defaulting to `<current working directory>/OUTPUT`)

## Usage

#### Arguments

| Argument | Verbose | Description | Type | Required | Default |
|:--------:|:-------:| ----------- |:----:|:--------:|:-------:|
| `-f`       | `--file`  | Source fasta file containing the sample probe records to run through tblastn | `Filepath` | True | N/A |
| `-db`      | `--database_path` | Path to the genome database against which the probe records are to be BLASTed | `str` |True | N/A |
| `-o`       | `--output_dir` | Location to which to write the result files | `str` | False | `<current_working_directory>/OUTPUT` |
| `-a`       | `--alignment_length_threshold` | Minimum length threshold that BLAST result alignment sequence lengths should exceed | `int` | False | `400` | 

#### Example

``
python src/probe_blaster.py -f data/fasta_file.fasta -db dbs/genome_db -o results/probe_blaster_output
``

#### Output

_NB: Column headers are not present in the data output, they are provided here for illustrative purposes only_

<div class="data-table-wrapper> markdown="block">

<table class="data-table">
<tr>
<th>Accession ID</th> 
<th>Scaffold</th> 
<th>Scaffold length</th> 
<th>Start</th> 
<th>End</th> 
<th>E-value</th> 
<th>Alignment length</th> 
<th>Acc. Seq</th> 
<th>Scaffold Alignment</th> 
<th>Frame</th> 
</tr>
<tr>
<td><code>XXX11111.1</code></td> 
<td><code>Scaf_00000111</code></td> 
<td><code>12345678</code></td> 
<td><code>992016</code></td> 
<td><code>994340</code></td> 
<td><code>1.63E-114</code></td> 
<td><code>797</code></td> 
<td><code>MVLDLKDCFFSIPLAEQDRE...</code></td> 
<td><code>LILDLKDCFFTIPLHPDD...</code></td> 
<td><code>1</code></td>
</tr>
<tr>
<td><code>----------</code></td> 
<td><code>-------------</code></td> 
<td><code>--------</code></td> 
<td><code>------</code></td> 
<td><code>------</code></td> 
<td><code>---------</code></td> 
<td><code>---</code></td> 
<td><code>-----------------------</code></td> 
<td><code>---------------------</code></td> 
<td><code>-</code></td> 
</tr>
<tr>
<td><code>----------</code></td> 
<td><code>-------------</code></td> 
<td><code>--------</code></td> 
<td><code>------</code></td> 
<td><code>------</code></td> 
<td><code>---------</code></td> 
<td><code>---</code></td> 
<td><code>-----------------------</code></td> 
<td><code>---------------------</code></td> 
<td><code>-</code></td> 
</tr>
<tr>
<td><code>----------</code></td> 
<td><code>-------------</code></td> 
<td><code>--------</code></td> 
<td><code>------</code></td> 
<td><code>------</code></td> 
<td><code>---------</code></td> 
<td><code>---</code></td> 
<td><code>-----------------------</code></td> 
<td><code>---------------------</code></td> 
<td><code>-</code></td> 
</tr>
</table>

</div>