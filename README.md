<foreignObject> <style>
body{color:#edf2f4ff !important;}
H1{color:#80ed99ff !important;}
H2{color:#57cc99ff !important;}
p{color:#5bc0beff !important;}
span{color:#fdc500ff !important;}
p, ol{color:#6fffe9ff !important;}
li{color:#38a3a5ff !important;}
code {
  font-family: Consolas,"courier new";
  color: #6fffe9ff;
  background-color: #3a506bff;
  padding: 2px;
  font-size: 105%;
}
</style>
<link rel="stylesheet"
      href="//cdnjs.cloudflare.com/ajax/libs/highlight.js/11.0.1/styles/default.min.css">
<script src="//cdnjs.cloudflare.com/ajax/libs/highlight.js/11.0.1/highlight.min.js"></script></foreignObject> 
<body>

# JDA Batch Summary
<p>
The samples were sequenced using Illumina in NOVOGENE genomics. </br>
The total output of data on the sequencer: <span>Raw data 95.6 GB.</span> </br>
The detail statistics for the quality of sequencing data are shown in Table 1.</br>



| Sample      | Library Flowcell Lane           | Raw reads      | Raw Data | Effective(%) | Error(%) | Q20(%) | Q30(%)| GC(%) |
| :---------: | :-----------------------------: | :------------: | :------: | :----------: | :------: | :----: | :---: | :---: |
| JDA5        | FDSW210086500-1r_H33GCDSX2_L4   | 220638234      |   33.1   |    99.33     |  0.03    |  97.81 | 93.65 | 35.77 |
| JDA6        | FDSW210086501-1r_H33GCDSX2_L4   | 216705178      |   32.5   |    99.16     |  0.03    |  97.57 | 93.11 | 36.04 |
| JDA25       | FDSW210086502-1r_H33GCDSX2_L4   | 200234308      |   30.0   |    99.24     |  0.03    |  97.76 | 93.53 | 35.89 |

</br>
<p style="color:#6fffe9ff;">
<span>Sample:</span> sample name</br>
<span>Library_Flowcell_Lane:</span> Library ID_Flowcell ID_lane ID, for raw data file naming.</br>
<span>Raw reads:</span> Total amount of reads of raw data, each four lines taken as one unit. For paired-end sequencing, it equals the amount of read1 and read2, otherwise it equals the amount of read1 for single-end sequencing.</br>
<span>Raw data:</span> (Raw reads) * (sequence length), calculating in G. For paired-end sequencing like PE150, sequencing length equals 150, otherwise it equals 50 for sequencing like SE50.</br>
<span>Effective:</span> (Clean reads/Raw reads)*100%</br>
<span>Error:</span> base error rate</br>
<span>Q20, Q30:</span> (Base count of Phred value > 20 or 30) / (Total base count)</br>
<span>GC:</span> (G & C base count) / (Total base count)
</p>
</p>

## Sequencing Data Format
<p>
The original raw data from Illumina platform are transformed to Sequenced Reads, known as Raw Data or RAW Reads, by base calling. Raw data are recorded in a FASTQ file, which contains sequencing reads and corresponding sequencing quality. Every read in FASTQ format is stored in four lines, as indicated below <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847217/">(Cock P.J.A. et al. 2010):</a>
</p>
<span>
@HWI-ST1276:71:C1162ACXX:1:1101:1208:2458 1:N:0:CGATGTNAAGAACACGTTCGGTCACCTCAGCACACTTGTGAATGTCATGGGATCCAT
+#55???BBBBB?BA@DEEFFCFFHHFFCFFHHHHHHHFAE0ECFFD/AEHH
</span>
</br>

<p>Line 1 begins with a <span>'@'</span> character and is followed by the Illumina Sequence Identifiers and an optional description. The details of Illumina sequence identifier are listed in the table below.
</p>

|Identifier            |Meaning                                                       |
| :------------------  | :----------------------------------------------------------  |
|HWI-ST1276 Instrument | unique identifier of the sequencer                           
|71 Run number -       | Run number on instrument                                     
|C1162ACXX Flow Cell ID| - ID of flow cell
|1                     | Lane Number - positive integer
|1101                  | Tile Number - positive integer
|1208                  | X - x coordinate of the spot. Integer which can be negative
|2458                  | Y - y coordinate of the spot. Integer which can be negative
|1                     | Read number. 1 can be single read or Read 2 of paired-end.
|N                     | Y if the read is filtered (did not pass), N otherwise.
|0                     | Control number - 0 when none of the control bits are on, otherwise it is an even number
|C                     | GATGTIllumina index sequences
</br>

- Line 2 is the raw sequence of the read.

- Line 3 begins with a '+' character and is optionally followed by the same sequence identifiers and descriptions as in Line 1.

- Line 4 encodes the quality values for the bases in Line 2 and contains the same number of characters as the bases in the read <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847217/">(Cock P.J.A. et al. 2010)</a>.

## Explanation of Sequencing Data Related
<ol>
<li>The data delivered is a compressed file in format of <code>.fq.gz</code>. Before data delivery, we will calculate the md5 value of each compressed file and please check it when you get the data. There are two ways to check the md5 value. In Linux environment, you can use <code>md5sum -c <*md5.txt></code> command under the data directory. In Windows environment, you can use a calibration tool e.g. hashmyfiles. If the md5 value of compressed file doesn't match with the one we provide in md5 file in data directory, the file may have been damaged during the transmitting procedure.</li>
<li> For paired-end (PE) sequencing, every sample should have 2 data flies (read1 file and read2 file). These 2 files have the same line number, you could use <code>wc -l</code> command to check the line number in Linux environment. The line number divide by 4 is the number of reads.</li>
<li> The data size is the space occupied by the data in the hard disk. It's related to the format of the disk and compression ratio. And it has no influence on the quantity of sequenced bases. Thus, the size of read1 file may be unequal to the size of read2 file.</li>
<li> When customer's samples need large amount of data e.g. whole genome sequencing data, we would use separate-lane sequencing strategy to make sure the quality of data, so it's possible that one sample has several parts sequencing data. For example, if sample 1 has two read1 files, <span>sample1_L1_1.fq.gz</span> and <span>sample1_L2_1.fq.gz</span>, that means this sample was sequenced on different lanes.</li>
<li> If we agree to deliver the clean data before the project starts, we will filter the data strictly according to the standard to obtain high quality clean data which can be used for further research and paper writing. We will discard the paired reads in the following situation: when either one read contains adapter contamination; when either one read contains more than 10 percent uncertain nucleotides; when either one read contains more than 50 percent low quality nucleotides (base quality less than 5). The data analysis results based on the clean data that is filtered by this standard can be approved by high level magazines <a href="https://www.nature.com/articles/nsmb.2660?page=9">(Yan L.Y. et al . 2013)</a>. If you want to get more information, please refer to the official website of Novogene <a href="[www.novogene.com](https://en.novogene.com/)">novogene</a>.</li>
<li> The Index is normally in the middle of the adapter during the process of experiment and sequencing except the special library. We can get the Read1 sequence and Read2 sequence according to the Index read. They are all the sequence of samples so that it's no necessary to dispose the beginning and end of reads in the downstream analysis (e.g. mapping).</li>
<li> 30 days after the data delivery, we will delete outdated data. So please keep your data properly. If you have any questions or concerns, please contact us as soon as possible.</li>

## Sequences of adapter

  * 5' Adapter:

 <span> 5'-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT-3'</span>

  * 3' Adapter:

  <span>5'-GATCGGAAGAGCACACGTCTGAACTCCAGTCAC-3'</span>

# Steps to process the FASTQ files into BAM (Binary Alignment Mapping) files
<img src="./fig01.png"></img>
<p><strong>Quality  control of raw reads</strong> </p>
<p>Quality  control of RNA-seq raw reads consists of analysis of sequence quality, GC  content, adaptor content, overrepresented k-mers, and duplicated reads,  dedicated to detecting sequencing errors, contaminations, and PCR artifacts.  Read quality decreases towards the 3' end of reads, bases with low quality,  therefore, they should be removed to improve mappability. In addition to the  quality of raw data, quality control of raw reads also includes the analysis of  read alignment (read uniformity and GC content), quantification (3' bias,  biotypes, and low-counts), and reproducibility (correlation, principal  component analysis, and batch effects).</p>
<p align="center">Table 1. The tools for quality control of  RNA-seq raw reads.</p>
<div class=" table-responsive">
<table border="0" cellspacing="0" cellpadding="0" class=" table-bordered tablecontentshow ke-zeroborder">
<tbody>
<tr>
<td width="177"><strong>Tools</strong> </td>
<td width="482"><strong>Applications</strong> </td>
</tr>
<tr>
<td width="177">FastQC</td>
<td width="482">Quality control of raw reads generated by any platforms.</td>
</tr>
<tr>
<td width="177">Trimmomatic</td>
<td width="482">Trimming the adapter content present in reads.</td>
</tr>
</tbody>
</table>
</div>
<p><strong>Read alignment</strong> </p>
<p>There are generally three strategies for read alignment,  genome mapping, transcriptome mapping, and <em>de novo</em> assembly. Regardless of whether a  genome or transcriptome reference is available, reads may map uniquely or be  assigned to multiple positions in the reference, which are referred to as  multi-mapped reads or multireads. Genomic multireads are generally due to  repetitive sequences or shared domains of paralogous genes. Transcriptome multi-mapping  arises more often due to gene isoforms. Therefore, transcript identification  and quantification are important challenges for alternatively expressed genes. When a reference is not available, RNA-seq reads are assembled <em>de novo</em> using  packages such as SOAPdenovo-Trans, Oases, Trans-ABySS, or Trinity. PE  strand-specific and long-length reads are preferred since they are more  informative. Emerging long-read technologies, such as PacBio SMRT sequencing and Nanopore sequencing, can generate full-length transcripts for most genes.</p>
<p style="align:center;"><img border="0" src="./fig22.png"> <br />Figure  2. Three basic strategies for RNA-seq read mapping (Conesa <em>et al</em>. 2016).  Abbreviations: GFF, General Feature Format; GTF, gene transfer format; RSEM,  RNA-seq by Expectation Maximization.</p>
<p align="center">Table  2. The comparison of genome-based and <em>de novo</em> assembly strategies for RNA-seq  analysis.</p>
<div class=" table-responsive">
<table border="0" cellspacing="0" cellpadding="0">
<tbody>
<tr>
<td width="151"><strong> </strong></td>
<td width="273"><strong>Genome-based</strong> </td>
<td width="288"><strong><em>De novo</em></strong><strong> assembly</strong></td>
</tr>
<tr>
<td width="151">Method</td>
<td width="273">Alignment to a reference genome</td>
<td width="288">Not using a reference genome</td>
</tr>
<tr>
<td width="151">Advantages</td>
<td width="273">
<ul class=" ullist">
<li>Efficient computing</li>
<li>Eliminates contaminating reads</li>
<li>Very sensitive and can assemble transcripts of low abundance</li>
<li>Can discover novel transcripts without annotation</li>
</ul>
</td>
<td width="288">
<ul class=" ullist">
<li>Reference genome is not required</li>
<li>Correct alignment of reads to known splice site is not required</li>
<li>Trans-spliced transcripts can be assembled</li>
</ul>
</td>
</tr>
<tr>
<td width="151">Disadvantages</td>
<td width="273">Requires high-quality reference genome</td>
<td width="288" valign="top">
<ul>
<li>More computational intense</li>
<li>Sensitive to sequencing error</li>
</ul>
</td>
</tr>
<tr>
<td width="151">Recommended depth</td>
<td width="273">Approximately 10x</td>
<td width="288">Beyond 30x</td>
</tr>
</tbody>
</table>
</div>
<p align="center">Table 3. The public sources of RNA-seq data.</p>
<div class="table-responsive">
<table width="712" border="0" cellspacing="0" cellpadding="0">
<tbody>
<tr>
<td width="237" valign="top"><strong>Transcriptomic Database</strong></td>
<td width="237"><strong>Data Type</strong></td>
<td width="237"><strong>Website</strong></td>
</tr>
<tr>
<td width="237">Gene Expression Omnibus (GEO)</td>
<td width="237">Both microarray and sequencing data</td>
<td width="237">https://www.ncbi.nlm.nih.gov/geo/</td>
</tr>
<tr>
<td width="237">ArrayExpress</td>
<td width="237">Both microarray and sequencing data</td>
<td width="237">https://www.ebi.ac.uk/arrayexpress/</td>
</tr>
<tr>
<td width="237">ENCODE: Encyclopedia of DNA Elements</td>
<td width="237">Public ENCODE Consortium data</td>
<td width="237">https://www.encodeproject.org/</td>
</tr>
<tr>
<td width="237">Sequence Read Archive (SRA)</td>
<td width="237">Sequencing data</td>
<td width="237">https://www.ncbi.nlm.nih.gov/sra</td>
</tr>
<tr>
<td width="237">European Nucleotide Archive (ENA)</td>
<td width="237">Sequencing data</td>
<td width="237">https://www.ebi.ac.uk/ena</td>
</tr>
<tr>
<td width="237">DDBJ Sequence Read Archive (DRA)</td>
<td width="237">Sequencing data</td>
<td width="237">https://www.ddbj.nig.ac.jp/dra</td>
</tr>
</tbody>
</table>
</div>
<p align="left"><strong>Transcript  quantification</strong> </p>
<p align="left">	Transcript quantification can be used to estimate gene and transcript expression levels.</p>

<img border="0" src="./fig3.png"> <br/><center>Figure 3. The tools for isoform expression quantification.</center>
<p><strong>Differential  expression testing</strong> </p>
<p>Differential expression testing is used to evaluate if  one gene is differentially expressed in one condition compared to the other(s).  Normalizing methods need to be adopted before comparing different samples. RPKM  and TPM normalize away the most important factor, sequencing depth. TMM, DESeq,  and UpperQuartile can ignore highly variable and/or highly expressed features.  Other factors that interfere with intra-sample comparisons involve transcript  length, positional biases in coverage, average fragment size, and GC content,  which can be normalized by tools, such as DESeq, edgeR, baySeq, and NOISeq.  Batch effects may still be present after normalization, which can be minimized  by appropriate experimental design, or removed by methods such as COMBAT or  ARSyN.</p>
<p align="center">Table 5. The normalization tools for  differential expression testing.</p>
<div class="table-responsive">
<table border="0" cellspacing="0" cellpadding="0" width="663">
<tbody>
<tr>
<td width="79" valign="top"><strong>Package</strong></td>
<td width="225"><strong>Read count distribution assumptions</strong> </td>
<td width="140"><strong>Input</strong></td>
<td width="96"><strong>Replicates</strong> </td>
<td width="123"><strong>Normalization</strong> </td>
</tr>
<tr>
<td width="79">DESeq</td>
<td width="225">Negative binomial distribution</td>
<td width="140">Raw counts</td>
<td width="96">No</td>
<td width="123">Library size</td>
</tr>
<tr>
<td width="79">edgeR</td>
<td width="225">Bayesian methods for negative binomial  distribution</td>
<td width="140">Raw counts</td>
<td width="96">Yes</td>
<td width="123">Library size<br />TMM<br />RLE<br />Upperquartile</td>
</tr>
<tr>
<td width="79">baySeq</td>
<td width="225">Bayesian methods for negative binomial distribution</td>
<td width="140">Raw counts</td>
<td width="96">Yes</td>
<td width="123">Library size<br />Quantile<br />TMM</td>
</tr>
<tr>
<td width="79">NOISeq</td>
<td width="225">Non-parametric</td>
<td width="140">Raw or normalized counts</td>
<td width="96">No</td>
<td width="123">Library size<br />RPKM<br />TMM<br />Upperquartile</td>
</tr>
</tbody>
</table>
</div>

# Steps involve in this process to convert BAM files

## 1. Downloading the Raw data

The data is present on the cloud storage bucket [plantik-prod-data](https://console.cloud.google.com/storage/browser/plantik-prod-data/data_type%3Dsra/JDA/raw_data?authuser=1&cloudshell=false&orgonly=true&project=plantik-prod&supportedpurview=organizationId&pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false) in the file [path](plantik-prod-data/data_type=sra/JDA/raw_data) linked here. You can download the data to your cloud space by using the command <code>gsutil -m cp -r gs://plantik-prod-data/data_type=sra/JDA/raw_data/ /path/you/wanto/save/the/data/</code> (in this command the <span>gsutil</span> serves as a python application that lets you access Cloud storage from the command line and <span>-m</span> refers to the multi-threaded which means it uses multiple-processing at a time to parallelize the process, <span>cp</span> stands for copying the files from one location to other location and <span>-r</span> enables directories, buckets, and bucket subdirectories to be copied recursively. If you don't use this option it will skip the directories).  

## 2. Preprocessing the Raw Data

The data downloaded is raw type can contain noise, unwanted text or some kind of useless data in it so, we need to clean the data at first hand to proceed for further proccesses. To determine this we actually do the FASTQC check on the raw data <span>e.g.two read1 files, sample1_L1_1.fq.gz</span> and <span>sample1_L2_1.fq.gz</span>

Forgot to tell you few points :thinking: here to setup this we need to have tools as follows:
* FASTQC
* Trimmomatic
* HISAT2
* SAMTOOLS
* Seqtk

And links to download and install them is below as follows :point_down:

* [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip)</br>
Actually installing `FastQC` is as simple as unzipping the zip file it comes in into a suitable location. That's it. Once unzipped it's ready to go.
```bash
~$ wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
~$ unzip fastqc_v0.11.9.zip
~$ chmod 755 fastqc
#but once you have done that you can run it directly
~$ ./fastqc
#or place a link in /usr/local/bin to be able to run the program from any location:
~$ sudo ln -s /path/to/FastQC/fastqc /usr/local/bin/fastqc
```

* [Trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip)</br>
* Actually installing `Trimmomatic` is as also same `FASTQC` by unzipping the zip file it comes into a suitable location. That's it. Once unzipped it's ready to go.
```bash
~$ wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
~$ unzip Trimmomatic-0.39.zip
~$ cd Trimmonatic-0.39/
~$ java -jar trimmomatic-0.39.jar --help # usgae commands will be displayed
```

* [HISAT2](https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download)</br>
It's one more simple process just unzip the file to use the commands: 
```bash
~$ wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download -O hisat2.zip
~$ unzip hisat2
~$ cd hisat2/
~$ ./hisat2 --help #help usage will be printed
```
* [Samtools](https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2)</br>

 See `INSTALL` in each of the source directories for further details.

The executable programs will be installed to a `bin` subdirectory under your specified prefix, so you may wish to add this directory to your `$PATH`: 
```bash
~$ cd samtools-1.x    # and similarly for bcftools and htslib
~$ ./configure --prefix=/where/to/install
~$ make
~$ make install
~$ export PATH=/where/to/install/bin:$PATH    # for sh or bash users
~$ setenv PATH /where/to/install/bin:$PATH    # for csh users
```

* [Seqtk](https://github.com/lh3/seqtk.git)</br>
The only library dependency is zlib should be present for this `seqtk` tool
```bash
~$ git clone https://github.com/lh3/seqtk.git;
~$ cd seqtk; make # 
~$ ./seqtk #usage details
```

:tada: Well Done! :clap: Now you are all set to go! :tada:

## 3. Quality control check on the raw reads using <span>FASTQC</span>
Modern high throughput sequencers can generate hundreds of millions of sequences in a single run. Before analysing this sequence to draw biological conclusions you should always perform some simple quality control checks to ensure that the raw data looks good and there are no problems or biases in your data which may affect how you can usefully use it.

Most sequencers will generate a QC report as part of their analysis pipeline, but this is usually only focused on identifying problems which were generated by the sequencer itself. FastQC aims to provide a QC report which can spot problems which originate either in the sequencer or in the starting library material.

FastQC can be run in one of two modes. It can either run as a stand alone `interactive application` for the immediate analysis of small numbers of FastQ files, or it can be run in a `non-interactive mode` where it would be suitable for integrating into a larger analysis pipeline for the systematic processing of large numbers of files.

```bash
~$ ./fastqc /path/to/fastq/gz/files/*.gz -o /path/to/output-dir/qc-files/
# output files will be in html format so download them at first to check the quality of raw reads
```
The ouput files consists of list of items which are important to check them:
* Basic Statistics
* Per base sequence quality
* Per tile sequence quality
* Per base sequence content
* Per sequence GC content
* Per base N content
* Sequence Length Distributuion
* Overrepresented sequences
* Adapter content
<ol type="I">
<strong><li>Basic Sequence Statistics:-</strong>
The Basic Statistics module generates some simple composition statistics for the file analysed.

* Filename: The original filename of the file which was analysed
* File type: Says whether the file appeared to contain actual base calls or colorspace data which had to be converted to base calls
* Encoding: Says which ASCII encoding of quality values was found in this file.
Total Sequences: A count of the total number of sequences processed. There are two values reported, actual and estimated. At the moment these will always be the same. In the future it may be possible to analyse just a subset of sequences and estimate the total number, to speed up the analysis, but since we have found that problematic sequences are not evenly distributed through a file we have disabled this for now.
* Filtered Sequences: If running in Casava mode sequences flagged to be filtered will be removed from all analyses. The number of such sequences removed will be reported here. The total sequences count above will not include these filtered sequences and will the number of sequences actually used for the rest of the analysis.
* Sequence Length: Provides the length of the shortest and longest sequence in the set. If all sequences are the same length only one value is reported.
* %GC: The overall %GC of all bases in all sequences
</li>
<strong><li>Per base sequence quality:-</strong></br>
This view shows an overview of the range of quality values across all bases at each position in the FastQ file.</br>

<img src="./boxplot.png"></br>

> For each position a BoxWhisker type plot is drawn. The elements of the plot are as follows:

* The central red line is the median value
* The yellow box represents the inter-quartile range (25-75%)
* The upper and lower whiskers represent the 10% and 90% points
* The blue line represents the mean quality

> The y-axis on the graph shows the quality scores. The higher the score the better the base call. The background of the graph divides the y axis into very good quality calls (green), calls of reasonable quality (orange), and calls of poor quality (red). The quality of calls on most platforms will degrade as the run progresses, so it is common to see base calls falling into the orange area towards the end of a read.

> It should be mentioned that there are number of different ways to encode a quality score in a FastQ file. FastQC attempts to automatically determine which encoding method was used, but in some very limited datasets it is possible that it will guess this incorrectly (ironically only when your data is universally very good!). The title of the graph will describe the encoding FastQC thinks your file used.

#### Warning
A warning will be issued if the `lower quartile` for any base is `less than 10`, or if the `median` for any base is `less than 25`.

#### Failure
This module will raise a failure if the `lower quartile` for any base is `less than 5` or if the `median` for any base is `less than 20`.
</br>

<center><strong style="color:#c7f9ccff;">The Good Illumina Data will be shown in below figure:</strong></center></br>
<img src="./goodboxplot.png"></br>
</br>

> Here for the plot it shows you that the reads in your raw data are in good quality. But we need to trim the adapters from the raw reads. We use `Trimmomatic` for further process.
</li>

<strong><li>Overrepresented sequences:-</strong></br>

A normal high-throughput library will contain a diverse set of sequences, with no individual sequence making up a tiny fraction of the whole. Finding that a single sequence is very overrepresented in the set either means that it is highly biologically significant, or indicates that the library is contaminated, or not as diverse as you expected.

This module lists all of the sequence which make up more than 0.1% of the total. To conserve memory only sequences which appear in the first 100,000 sequences are tracked to the end of the file. It is therefore possible that a sequence which is overrepresented but doesn't appear at the start of the file for some reason could be missed by this module.

For each overrepresented sequence the program will look for matches in a database of common contaminants and will report the best hit it finds. Hits must be at least 20bp in length and have no more than 1 mismatch. Finding a hit doesn't necessarily mean that this is the source of the contamination, but may point you in the right direction. It's also worth pointing out that many adapter sequences are very similar to each other so you may get a hit reported which isn't technically correct, but which has very similar sequence to the actual match.

Because the duplication detection requires an exact sequence match over the whole length of the sequence any reads over 75bp in length are truncated to 50bp for the purposes of this analysis. Even so, longer reads are more likely to contain sequencing errors which will artificially increase the observed diversity and will tend to underrepresent highly duplicated sequences.
#### Warning
This module will issue a warning if any sequence is found to represent `more than 0.1%` of the total.

#### Failure
This module will issue an error if any sequence is found to represent `more than 1%` of the total.
>And for the Overrepresented sequences we can consider till 2% of representation. You can ignore if they are `less than 2%`, but if they exceed the `2%` you need to trim the sequences.</br>

> These three are the most important features to check whenever the Fquality check is implemented to consider to step forward with the process.

## 4. Trimming the raw data using <span>Trimmomatic</span>

Trimmomatic is a fast, multithreaded command line tool that can be used to trim and crop
Illumina (FASTQ) data as well as to remove adapters. These adapters can pose a real problem
depending on the library preparation and downstream application. </br>
There are two major modes of the program: `Paired end mode` and `Single end mode`. The
paired end mode will maintain correspondence of read pairs and also use the additional
information contained in paired reads to better find adapter or PCR primer fragments
introduced by the library preparation process.</br>
Trimmomatic works with FASTQ files (using `phred + 33` or `phred + 64` quality scores,
depending on the Illumina pipeline used). Files compressed using either `„gzip‟ or „bzip2‟` are
supported, and are identified by use of `„.gz‟ or „.bz2‟` file extensions. </br>
## Implemented trimming steps (Quick reference)
Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single
ended data. The selection of trimming steps and their associated parameters are supplied on
the command line.</br>
The current trimming steps are:</br>
* ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
* SLIDINGWINDOW: Performs a sliding window trimming approach. It starts
scanning at the 5‟ end and clips the read once the average quality within the window
falls below a threshold.
* MAXINFO: An adaptive quality trimmer which balances read length and error rate to
maximise the value of each read
* LEADING: Cut bases off the start of a read, if below a threshold quality
* TRAILING: Cut bases off the end of a read, if below a threshold quality
* CROP: Cut the read to a specified length by removing bases from the end
* HEADCROP: Cut the specified number of bases from the start of the read
* MINLEN: Drop the read if it is below a specified length
* AVGQUAL: Drop the read if the average quality is below the specified level
* TOPHRED33: Convert quality scores to Phred-33
* TOPHRED64: Convert quality scores to Phred-64 </br>

## Running Trimmomatic

### Processing Order
The different processing steps occur in the order in which the steps are specified on the
command line. It is recommended in most cases that adapter clipping, if required, is done as
early as possible, since correctly identifying adapters using partial matches is more difficult.
#### Single End Mode
For single-ended data, one input and one output file are specified. The required processing
steps (trimming, cropping, adapter clipping etc.) are specified as additional arguments after
the input/output files.
```bash
~$ java -jar <path to trimmomatic jar> SE [-threads <threads>] [-phred33 | -phred64] [-trimlog
<logFile>] <input> <output> <step 1>
```
or
```bash
~$ java -classpath <path to trimmomatic jar> org.usadellab.trimmomatic.TrimmomaticSE [-
threads <threads>] [-phred33 | -phred64] [-trimlog <logFile>] <input> <output> <step 1>
``` 

> `-phred33` or `-phred64` specifies the base quality encoding. If no quality encoding is specified,
it will be determined automatically (since version 0.32). The prior default was <span>-phred64</span>.</br>

Specifying a trimlog file creates a log of all read trimmings, indicating the following details:</br>
* the read name
* the surviving sequence length
* the location of the first surviving base, aka. the amount trimmed from the start
* the location of the last surviving base in the original read
* the amount trimmed from the end
> Multiple steps can be specified as required, by using additional arguments at the end as
described in the section processing steps.

#### Paired End Mode
For paired-end data, two input files, and 4 output files are specified, 2 for the 'paired' output
where both reads survived the processing, and 2 for corresponding 'unpaired' output where a
read survived, but the partner read did not.

<img src="./TrimmPE.png" /></br><center><em>Figure shows the flow of reads in Trimmomatic Paired End mode</em></center>
```bash
java -jar <path to trimmomatic.jar> PE [-threads <threads] [-phred33 | -phred64] [-trimlog
<logFile>] >] [-basein <inputBase> | <input 1> <input 2>] [-baseout <outputBase> |
<unpaired output 1> <paired output 2> <unpaired output 2> <step 1>
```
or
```bash
java -classpath <path to trimmomatic jar> org.usadellab.trimmomatic.TrimmomaticPE [-
threads <threads>] [-phred33 | -phred64] [-trimlog <logFile>] [-basein <inputBase> | <input
1> <input 2>] [-baseout <outputBase> | <paired output 1> <unpaired output 1> <paired
output 2> <unpaired output 2> <step 1>
```

> `-phred33` or `-phred64` specifies the base quality encoding. If no quality encoding is specified,
it will be determined automatically (since version 0.32). The prior default was <span>-phred64</span>.</br>
> <span>-threads</span> indicates the number of threads to use, which improves performance on multi-core computers. If not specified, it will be chosen automatically. Specifying a trimlog file creates a log of all read trimmings, indicating the following details:</br>
* the read name
* the surviving sequence length
* the location of the first surviving base, aka. the amount trimmed from the start
* the location of the last surviving base in the original read
* the amount trimmed from the end
> Multiple steps can be specified as required, by using additional arguments at the end as
described in the section processing steps. 

*** 
### Examples

#### Paired-End
```bash
~$ java -jar trimmomatic-0.30.jar PE s_1_1_sequence.txt.gz s_1_2_sequence.txt.gz
lane1_forward_paired.fq.gz lane1_forward_unpaired.fq.gz lane1_reverse_paired.fq.gz
lane1_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3
TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
> This will perform the following <b>in this order</b></br>

* Remove Illumina adapters provided in the TruSeq3-PE.fa file (provided). Initially
Trimmomatic will look for seed matches (16 bases) allowing maximally <b>2</b>
mismatches. These seeds will be extended and clipped if in the case of paired end
reads a score of <b>30</b> is reached (about 50 bases), or in the case of single ended reads a
score of <b>10</b>, (about 17 bases).
* Remove leading low quality or N bases (below quality <b>3</b>)
* Remove trailing low quality or N bases (below quality <b>3</b>)
* Scan the read with a 4-base wide sliding window, cutting when the average quality per
base drops below <b>15</b>
* Drop reads which are less than <b>36</b> bases long after these steps</br>

#### Single End
```bash
~$ java -jar trimmomatic-0.30.jar SE s_1_1_sequence.txt.gz lane1_forward.fq.gz
ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
MINLEN:36
```

> This will perform the same steps, using the single-ended adapter file. (Of course the <b>:30:</b>
parameter has no effect, but a value has to be specified nevertheless)

> For more reference read the linked document [here!](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

### Output files produced by trimmomatic are as follows:

Paired-end mode requires 2 input files (for forward and reverse reads) and 4 output files (for
forward paired, forward unpaired, reverse paired and reverse unpaired reads).</br>
Since these files often have similar names, the user has the option to provide either the
individual file names, or just one name from which the file names can be derived.</br>
> For input files, either of the following can be used:</br>
* Explicitly naming the 2 input files
* Naming the forward file using the -basein flag, where the reverse file can be
determined automatically. The second file is determined by looking for common
patterns of file naming, and changing the appropriate character to reference the reverse
file.</br> Examples which should be correctly handled include:
  * `Sample_Name_R1_001.fq.gh` -> <span>Sample_Name_R2_001.fq.gz</span>
  * `Sample_Name.f.fastq` -> <span>Sample_Name.r.fastq</span>
  * `Sample_Name.1.sequence.txt` -> <span>Sample_Name.2.sequence.txt</span>
> For output files, either of the following can be used:
* Explicity naming the 4 output files
* Providing a base file name using the -baseout flag, from which the 4 output files can
be derived. </br>If the name "mySampleFiltered.fq.gz" is provided, the following 4 file
names will be used:
  * `mySampleFiltered_1P.fq.gz` - <span>for paired forward reads</span>
  * `mySampleFiltered_1U.fq.gz` - <span>for unpaired forward reads</span>
  * `mySampleFiltered_2P.fq.gz` - <span>for paired reverse reads</span>
  * `mySampleFiltered_2U.fq.gz` - <span>for unpaired reverse reads</span>
> For input and output files adding .gz to an extension tells Trimmomatic that the file is
provided in gzipped format or that Trimmomatic should gzip the file, respectively. This
extension can be used with both explicitly named and template-based file naming.

After Trimming sequenceswe need to perform the `FastQC` again on the files named with above mentioned output files with `paired forward`, `unpaired forward`, `paired reverse` and `unpaired reverse` reads.

```bash
~$ ./fastqc /path/of/the/trimmed/files/*.gz -o /path/to/save/output/files/
```
Check all the files with `.html` format and if the Reads have good quality then we can move the next step i.e Mapping reads.

## 5. Mapping the Reads and performing alignment with reference genome (<span>cs10</span>) using <span>HISAT2</span> &amp; <span>Seqtk</span>

## I. HISAT2:
> HISAT2 is a fast and sensitive alignment program for mapping next-generation sequencing reads (both DNA and RNA) to a population of human genomes as well as to a single reference genome.


## II. Seqtk:

>Seqtk is a fast and lightweight tool for processing sequences in the FASTA or FASTQ format. It seamlessly parses both FASTA and FASTQ files which can also be optionally compressed by gzip.

</body>
