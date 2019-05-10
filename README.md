# ERVin

This is a tool to allow for the detection of ERVs in genome segments

This has been designed primarily with a view to be used on OSX, cross-compatibility with other UNIX-based architectures may exist, but it almost certainly will not run on Microsoft Windows systems

### Requirements
- Python 3.6+ ([Download](https://www.python.org/downloads/))
- NCBI BLAST suite must be installed locally ([Download](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/))
- Local genome db to be queried
    - This should be located in a directory named `data` and specified either using the command line argument or by defining in the `config/config.json` file
    
### Current functionality
- [probe_blaster](docs/probe_blaster.MD)