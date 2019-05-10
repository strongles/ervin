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
<td class="data-table-cell"><code>-f</code></td>
<td class="data-table-cell"><code>--file</code></td>
<td class="data-table-cell">Source fasta file containing the sample probe records to run through tblastn</td>
<td class="data-table-cell"><code>Filepath</code></td>
<td class="data-table-cell">True</td>
<td class="data-table-cell"></td>
<tr>
<tr>
<td class="data-table-cell"><code>-db</code></td>
<td class="data-table-cell"><code>--database_path</code></td>
<td class="data-table-cell">Path to the genome database against which the probe records are to be BLASTed</td>
<td class="data-table-cell"><code>str</code></td>
<td class="data-table-cell">True</td>
<td class="data-table-cell"></td>
<tr>
<tr>
<td class="data-table-cell"><code>-o</code></td>
<td class="data-table-cell"><code>--output_dir</code></td>
<td class="data-table-cell">Location to which to write the result files</td>
<td class="data-table-cell"><code>str</code></td>
<td class="data-table-cell">False</td>
<td class="data-table-cell"><code>&lt;current_working_directory&gt;/OUTPUT</code></td>
<tr>
<tr>
<td class="data-table-cell"><code>-a</code></td>
<td class="data-table-cell"><code>--alignment_len_threshold</code></td>
<td class="data-table-cell">Minimum length threshold that BLAST result alignment sequence lengths should exceed</td>
<td class="data-table-cell"><code>int</code></td>
<td class="data-table-cell">False</td>
<td class="data-table-cell"><code>400</code></td>
<tr>
<tr>
<td class="data-table-cell"><code>-e</code></td>
<td class="data-table-cell"><code>--e_value</code></td>
<td class="data-table-cell">Maximum e-value threshold that BLAST result e-values should exceed</td>
<td class="data-table-cell"><code>float</code></td>
<td class="data-table-cell">False</td>
<td class="data-table-cell"><code>0.009</code></td>
<tr>

</table>
</div>


#### Examples
<div class="data-table-wrapper">

<code>python src/probe_blaster.py -f data/fasta_file.fasta -db dbs/genome_db</code>


<code>python src/probe_blaster.py -f data/fasta_file.fasta -db dbs/genome_db -o results/probe_blaster_output</code>


<code>python src/probe_blaster.py -f data/fasta_file.fasta -db dbs/genome_db -o results/probe_blaster_output -a 500</code>


<code>python src/probe_blaster.py -f data/fasta_file.fasta -db dbs/genome_db -o results/probe_blaster_output -e 0.0008</code>


<code>python src/probe_blaster.py -f data/fasta_file.fasta -db dbs/genome_db -o results/probe_blaster_output -a 800 -e 0.01</code>
</div>

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