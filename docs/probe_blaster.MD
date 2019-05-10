<link href="formatting/style.css" rel="stylesheet"></link>

# Probe Blaster

Probe blaster accepts a `.fasta` file containing genome probes, and runs them via `tblastn` against the 
provided genome database.

It will then provide a `.tsv` file on a per-fasta-record basis, output to the directory provided 
(defaulting to `<current working directory>/OUTPUT`)

## Usage

#### Arguments

<div class="data-table-wrapper">

<table class="data-table">

<tr>
<th>Argument</th>
<th>Verbose</th>
<th>Description</th>
<th>Type</th>
<th>Required</th>
<th>Default</th>
</tr>
<tr>
<td><code>-f</code></td>
<td><code>--file</code></td>
<td>Source fasta file containing the sample probe records to run through tblastn</td>
<td><code>Filepath</code></td>
<td>True</td>
<td></td>
<tr>
<tr>
<td><code>-db</code></td>
<td><code>--database_path</code></td>
<td>Path to the genome database against which the probe records are to be BLASTed</td>
<td><code>str</code></td>
<td>True</td>
<td></td>
<tr>
<tr>
<td><code>-o</code></td>
<td><code>--output_dir</code></td>
<td>Location to which to write the result files</td>
<td><code>str</code></td>
<td>False</td>
<td><code>&lt;current_working_directory&gt;/OUTPUT</code></td>
<tr>
<tr>
<td><code>-a</code></td>
<td><code>--file</code></td>
<td>Minimum length threshold that BLAST result alignment sequence lengths should exceed</td>
<td><code>int</code></td>
<td>False</td>
<td><code>400</code></td>
<tr>

</table>
</div>

| Argument | Verbose | Description | Type | Required | Default |
|:--------:|:-------:| ----------- |:----:|:--------:|:-------:|
| `-f`       | `--file`  | Source fasta file containing the sample probe records to run through tblastn | `Filepath` | True | N/A |
| `-db`      | `--database_path` | Path to the genome database against which the probe records are to be BLASTed | `str` |True | N/A |
| `-o`       | `--output_dir` | Location to which to write the result files | `str` | False | `<current_working_directory>/OUTPUT` |
| `-a`       | `--alignment_length_threshold` | Minimum length threshold that BLAST result alignment sequence lengths should exceed | `int` | False | `400` | 

</div>

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
<td class="data-table-cell"><code>XXX11111.1</code></td> 
<td class="data-table-cell"><code>Scaf_00000111</code></td> 
<td class="data-table-cell"><code>12345678</code></td> 
<td class="data-table-cell"><code>992016</code></td> 
<td class="data-table-cell"><code>994340</code></td> 
<td class="data-table-cell"><code>1.63E-114</code></td> 
<td class="data-table-cell"><code>797</code></td> 
<td class="data-table-cell"><code>MVLDLKDCFFSIPLAEQDRE...</code></td> 
<td class="data-table-cell"><code>LILDLKDCFFTIPLHPDD...</code></td> 
<td class="data-table-cell"><code>1</code></td>
</tr>
<tr>
<td class="data-table-cell"><code>----------</code></td> 
<td class="data-table-cell"><code>-------------</code></td> 
<td class="data-table-cell"><code>--------</code></td> 
<td class="data-table-cell"><code>------</code></td> 
<td class="data-table-cell"><code>------</code></td> 
<td class="data-table-cell"><code>---------</code></td> 
<td class="data-table-cell"><code>---</code></td> 
<td class="data-table-cell"><code>-----------------------</code></td> 
<td class="data-table-cell"><code>---------------------</code></td> 
<td class="data-table-cell"><code>-</code></td> 
</tr>
<tr>
<td class="data-table-cell"><code>----------</code></td> 
<td class="data-table-cell"><code>-------------</code></td> 
<td class="data-table-cell"><code>--------</code></td> 
<td class="data-table-cell"><code>------</code></td> 
<td class="data-table-cell"><code>------</code></td> 
<td class="data-table-cell"><code>---------</code></td> 
<td class="data-table-cell"><code>---</code></td> 
<td class="data-table-cell"><code>-----------------------</code></td> 
<td class="data-table-cell"><code>---------------------</code></td> 
<td class="data-table-cell"><code>-</code></td> 
</tr>
<tr>
<td class="data-table-cell"><code>----------</code></td> 
<td class="data-table-cell"><code>-------------</code></td> 
<td class="data-table-cell"><code>--------</code></td> 
<td class="data-table-cell"><code>------</code></td> 
<td class="data-table-cell"><code>------</code></td> 
<td class="data-table-cell"><code>---------</code></td> 
<td class="data-table-cell"><code>---</code></td> 
<td class="data-table-cell"><code>-----------------------</code></td> 
<td class="data-table-cell"><code>---------------------</code></td> 
<td class="data-table-cell"><code>-</code></td> 
</tr>
</table>

</div>