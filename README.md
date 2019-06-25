# ERVin

This is a tool to allow for the detection of ERVs in genome segments

This has been designed primarily with a view to be used on OSX, cross-compatibility with other UNIX-based architectures may exist, but it almost certainly will not run on Microsoft Windows systems

### Installation

`pip install ervin`

### Requirements
- Python 3.6+ ([Download](https://www.python.org/downloads/))
- NCBI BLAST suite must be installed locally ([Download](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/))
- Local genome db to be queried
    - This can be located in a directory of your choosing, but must be named in a `config.json` file
        - There is a `config.json.templ` file which will be used to create a `config.json` file from with the contained defaults at first run if you do not provide your own
 
    
### Current functionality
ERViN Currently:
- When provided with a `.fasta` file of probe sequences
    - Runs local `tblastn` against the specified genome database, filtering the results based on alignment length and e-value (optional arguments which result in default values of >400 and <0.009 respectively when omitted)
    - Parses and merges filtered results where appropriate 
    - Runs resultant fasta records against a local Viruses refseq database (a copy will be downloaded if not user provided, and will be kept up-to-date) using `tblastn`, grouping the records in a final set of output files based on their top hit
    
### Usage

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
<td class="data-table-cell"><code>-gdb</code></td>
<td class="data-table-cell"><code>--genome_database</code></td>
<td class="data-table-cell">Name of the genome database against which the probe records are to be BLASTed (located in the genome db store specified in the config file</td>
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

<code>ervin -f data/fasta_file.fasta -gdb genome_db</code>


<code>ervin -f data/fasta_file.fasta -gdb genome_db -o results/probe_blaster_output</code>


<code>ervin -f data/fasta_file.fasta -gdb genome_db -o results/probe_blaster_output -a 500</code>


<code>ervin -f data/fasta_file.fasta -gdb genome_db -o results/probe_blaster_output -e 0.0008</code>


<code>ervin -f data/fasta_file.fasta -gdb genome_db -o results/probe_blaster_output -a 800 -e 0.01</code>
</div>