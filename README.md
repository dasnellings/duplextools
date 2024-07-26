# duplextools

<!-- badges: start -->

<!-- badges: end -->

<img src="https://github.com/user-attachments/assets/895780b3-ea8e-48f8-8cd5-a5b9ca3cc07d" width=400 align="right">

duplextools is a collection of tools for processing and analyzing data generated from a range of duplex sequencing methods.

duplextools can be used to:
1) Extract and trim duplex barcodes from fastq files
2) Pair reads by molecule of origin (a read family)
3) Call variants from read families
4) Filter variant calls
5) Compute genome-wide mutation burden

### Installation
1. [Install go.](https://go.dev/doc/install)
2. Run `go install github.com/dasnellings/duplextools/...@latest`

`duplextools` binary will be present in `~/go/bin`

### Supported Methods
The included tools are designed to be flexible for use with any duplex sequencing method that either 
1) uses duplex barcodes inline with read1 and read2
2) the user can provide a fastq containing the duplex barcodes

duplextools has been tested (and has presets) for the following methods
* [META-CS](https://pubmed.ncbi.nlm.nih.gov/33593904/)
