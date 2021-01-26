---
title: "Biopython Tutorial: from beginner to advanced."
layout: post
date: 2021-01-05 11:38
image: /assets/images/biopython_logo.svg
headerImage: true
tag:
- Python
- Biopython
- Bioinformatics
star: false
category: blog
hidden: false
author: davidboo
description: Biopython tutorial and quick tips
---
Biopython is a Python library for reading and writing many common biological data formats.

Biopython is a Python library that allows us to perform bioinformatics computations on many common biological data formats. It is a powerful, useful library which provides a wide range of functions from reading large files with biological data.


#### 1. First Steps
   * What is Biopython
   * Installing and upgrading

#### 2. Working with Sequences
   * Creating a sequence
   * Parsing a sequence file
   * Counting sequence length & number of occurrences of a nucleotide
   * Calculating GC-content
   * Calculating molecular weight 
   * Get the reverse complement of a sequence, transcription & translation
   * Slicing a sequence
   * Concatenating sequences
   * Finding the starting index of a subsequence
   * Writing sequences to a file
   * Converting a FASTQ file to FASTA file & other formats
   * Separate sequences by ids from a list of ids

#### 3. NCBI Entrez databases
   * General Guidelines
   * Accessing Pubmed, Medline, Genbank & others
   
#### 4. BLAST
   * Running a web BLAST
   * Parsing a BLAST output
   * Other sequence search tools: SearchIO (QueryResult, Hit...)
   
#### 5. Multiple Sequence Alignment
   * Reading a MSA
   * Creating an alignment using different algorithms: ClustalW, MUSCLE...
   * Pairwise sequence alignment
   
#### 6. Phylogenetics
   * Constructing a phylogenetic tree
   * Modifying an existing tree
   
#### 7. Sequence motif analysis

#### 8. PDB: 3D structure protein analysis
   * Count atoms in a PDB structure
   
---

#### 1. First Steps

Our first step is installing Biopython. The easiest way is to Biopython is using ```pip``` package such as follows. (Please note that Biopython also requires NumPy package, which will be installed automatically if you install Biopython with ```pip```).

If you already have Biopython installed you can upgrade it using ```pip``` aswell. 

```python
pip install biopython
pip install --upgrade biopython
```

You can test if Biopython is properly installed by executing:

```python
import Bio
```

If no error messages appear, then everything is working and you're all set.

---

#### 2. Working with sequences

**Creating a sequence**

The easiest way to create a sequence is using the Biopython ```Seq``` object. The bio.seq module is really useful and we will make use of it constantly. You can check the [Wiki page](https://biopython.org/wiki/Seq) for more information and advice. Here is a small example.

```python
from Bio.Seq import Seq
seq_example = Seq("AGTACACTGGT")
print("Sequence is", seq_example)
>>> Sequence is AGTACACTGGT
```

**Parsing a sequence file**

[Biopython’s SeqIO](https://biopython.org/wiki/SeqIO) (Sequence Input/Output) can be used to read sequence files. The main function is ```Bio.SeqIO.parse()```, which requires a filename and format and returns a SeqRecord iterator. This is an example of parsing a FASTA file.
```python
from Bio import SeqIO
for record in SeqIO.parse("example_sequence.fasta", "fasta"):
    print(record.seq)
```
```record.id``` will return the identifier of the sequence. ```record.seq``` will return the sequence itself. ```record.description``` will return the sequence description.

**Counting sequence length & number of occurrences of a nucleotide**

Counting sequence length is really easy using ```len()```. Using our previous example:

```python
from Bio.Seq import Seq
seq_example = Seq("AGTACACTGGT")
print("Sequence is", seq_example)
>>> Sequence is AGTACACTGGT
print("The length of the sequence is", len(seq_example))
>>> The length of the sequence is 11
```
You can count numer of occurrences of a certain nucleotide as follows: (Make sure to create/load sequence as we've done earlier)

```python
print(seq_example.count("C"))
>>> 2
```

**Calculating GC-content**

[GC content](https://en.wikipedia.org/wiki/GC-content) is a measure of the percentage of guanine (G) or cytosine (C) related to the total nitrogenous bases (A+T on DNA or A+U on RNA). GC content is important because it indicates a higher melting temperature and also key factor in bacteria taxonomy. Biopython has a ```Bio.SeqUtils``` [package](https://biopython.org/docs/1.75/api/Bio.SeqUtils.html) that simplifies this task. Feel free to check wiki for more information.

```python
from Bio.SeqUtils import GC
from Bio.Seq import Seq
seq_example = Seq("AGTACACTGGT")
print(GC(seq_example))
>>> 45.45
```

**Calculating molecular weight**

Again, molecular weight can be easily calculated using ```Bio.SeqUtils``` as follows:

```python
from Bio.SeqUtils import molecular_weight
from Bio.Seq import Seq
seq_example = Seq("AGTACACTGGT")
print(molecular_weight(seq_example))
>>> 3436.19
```

**Get the reverse complement of a sequence, transcription & translation**

Using Biopython these three processes ar really straightforward. To get the reverse complement use the single function ```reverse_complement()```.

```python
print("Reverse complement:", seq_example.reverse_complement())
>>> The reverse complement if the sequence is ACCAGTGTACT
```

Transcription is also a single function, ```transcribe()```

```python
print("Transcription:", seq_example.transcribe())
>>> Transcription: AGUACACUGGU
```

Translation is again another function, ```translate()```. Please note that if your sequence is not suitable for correct translation you will recieve a warning from Biopython (our example gets this warning). This happens if your sequence is not a multiple of three, which makes partial codons on translation (we'll come back to this later). You can have more information about this biological process [here](https://www.genome.gov/genetics-glossary/Translation)

```python
print("Translation:", seq_example.translate())
>>> Translation: STL
```

**Slicing a sequence**

With Biopython you can slice effectively sequences and extract any part of it.

```python
from Bio.Seq import Seq
seq_example = Seq("AGTACACTGGT")
print(seq_example[4:9])
>>> CACTG
```

It follows the usual indexing conventions for Python strings, with the first element of the sequence numbered 0. When you do a slice the first item is included (i.e. 4) and the last is excluded (9 in this case). We can get the third codon positions of this DNA sequence:

```python
from Bio.Seq import Seq
seq_example = Seq("AGTACACTGGT")
print(seq_example[0::3])
>>> AACG
```

**Concatenating sequences**

You can add any two ```Seq``` objects together

```python
from Bio.Seq import Seq
seq_example = Seq("AGTACACTGGT")
seq_add = Seq("TTTT")
print(seq_example + seq_add)
>>> AGTACACTGGTTTTT
```

Biopython ```Seq``` also has a ```.join``` method which can be really useful

```python
from Bio.Seq import Seq
fragments = [Seq("AGTA"), Seq("CACT"), Seq("GGT")]
filler = Seq("G"*10)
print(filler.join(fragments))
>>> AGTAGGGGGGGGGGCACTGGGGGGGGGGGGT
```

**Find the starting index of a subsequence**

Easiest way to find starting index of a subsequence is using ```find()```

```python
from Bio.Seq import Seq
seq_example = Seq("AGTACACTGGT")
print("CAC index:", seq_example.find("CAC"))
>>> CAC index: 4
```

**Writing sequences to a file**

SeqIO (Sequence Input/Output) package can be used to write sequences to files. Our main function is ```SeqIO.write()``` used as:

```python
SeqIO.write(records, file, format)
```

Where records will be a list of sequence records you wish to save, file is a file open for writing and format should be ‘fasta’.

ME QUEDO AQUÍ DAVID

```python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
sequences = ["AAACGTGG", "TGAACCG", "GGTGCA", "CCAATGCG"]
records = (SeqRecord(Seq(seq, generic_dna), str(index)) for index,seq in enumerate(sequences))
with open("example.fasta", "w") as output_handle:
    SeqIO.write(records, output_handle, "fasta")
# This code will result in a FASTA file with sequence ids starting from 0. If you want to give a custom id and a description you can create the records as follows.
sequences = ["AAACGTGG", "TGAACCG", "GGTGCA", "CCAATGCG"]
new_sequences = []
i=1
for sequence in sequences:
    record = SeqRecord(sequence, id="Seq_"+str(i), name="Seq_"+str(i), description="<custom description>")
    new_sequences.append(record)
with open("example.fasta", "w") as output_handle:
    SeqIO.write(new_sequences, output_handle, "fasta")
# The SeqIO.write() function will return the number of sequences written.
```

**Converting a FASTQ file to FASTA file & other formats**

We need to convert DNA data file formats in certain applications. For example, we can do file format conversions from FASTQ to FASTA as follows.
```python
from Bio import SeqIO
with open("path/to/fastq/file.fastq", "r") as input_handle, open("path/to/fasta/file.fasta", "w") as output_handle:
    sequences = SeqIO.parse(input_handle, "fastq")        
    count = SeqIO.write(sequences, output_handle, "fasta")        
print("Converted %i records" % count)
```
If you want to convert a GenBank file to FASTA format,
```python
from Bio import SeqIO
with open("path/to/genbank/file.gb", "rU") as input_handle, open("path/to/fasta/file.fasta", "w") as output_handle:
    sequences = SeqIO.parse(input_handle, "genbank")
    count = SeqIO.write(sequences, output_handle, "fasta")
print("Converted %i records" % count)
```

**Separate sequences by ids from a list of ids**

Assume that you have a list of sequence identifiers in a file named list.lst where you want to separate the corresponding sequences from a FASTA file. You can run the following and write those sequences to a file.
```python
from Bio import SeqIO
ids = set(x[:-1] for x in open(path+"list.lst"))
with open(path+'list.fq', mode='a') as my_output:
    
    for seq in SeqIO.parse(path+"list_sequences.fq", "fastq"):
        
        if seq.id in ids: 
            my_output.write(seq.format("fastq"))
```

#Final Thoughts
#Hope you got an idea of how to use Seq, SeqRecord and SeqIO Biopython functions and will be useful for your research work.
#Thank you for reading. I would love to hear your thoughts. Stay tuned for the next part of this article with more usages and Biopython functions.
#Cheers, and stay safe!

---

#### 3. NCBI Entrez databases

**General Guidelines**

[Entrez](https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html) is a data retrieval system that provides users access to NCBI’s databases such as PubMed, GenBank, GEO, and many others. Using Biopython's module ```Bio.Entrez``` allows you to access Entrez and search & download records from within a Python script. You must provide an email to access Entrez - NCBI can block anonymous requests. Also,  ```Bio.Entrez``` might not be the best option to download huge amounts of data as it can clog the servers. 

```python
from Bio import Entrez
Entrez.email = "davidboo@example.com"
```

**Accessing PubMed, Medline, Genbank & others**

To search any of these databases, we use ```Bio.Entrez.esearch()``` module, which requires some parameters:

```python
handle = Entrez.esearch(db="value", term="keywords", retmax=100)
```

Where db are databases such as Pubmed, nucleotide, protein, snp, omim, unigene... keywords are what you are interested in looking at (organisms, Human, Biopython, urea...) and retmax are the number of identifies returned. 

Let's do an example accessing PubMed. We will be searching in PubMed for 20 publications that include breast cancer in their title.

```python
from Bio import Entrez
Entrez.email = "davidboo@example.com"
handle = Entrez.esearch(db="pubmed", term="breast[title] AND cancer[title]", retmax=20)
record = Entrez.read(handle)
identifiers = records['IdList']
```
handle = Entrez.esearch(db="value", term="keywords", retmax=100)


---


#### 4. BLAST

**Running Web BLAST**
Using Biopython, you can align sequences with Web BLAST which is the online version of BLAST. For this, we will be using the qblast() function in the Bio.Blast.NCBIWWW module.
You can check the documentation as follows.
```python
from Bio.Blast import NCBIWWW
help(NCBIWWW.qblast)
```

If you have a FASTA file, you can search it against the nucleotide database (nt) using Nucleotide BLAST (BLASTN) as follows.
```python
from Bio import SeqIO
record = SeqIO.read("sample.fasta", format="fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))
Now we can save the result to a file.
# with open("my_blast_result.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```

**Parsing a BLAST output**

We can load our saved BLAST result as follows.

```python
result_handle = open("my_blast_result.xml")
```
Let’s check the returned hits. Since we have obtained the BLAST output in XML format, we can parse the result using NCBIXML. Here we have used one query sequence and hence we get only one record.
```python
from Bio.Blast import NCBIXML
blast_record = NCBIXML.read(result_handle)
```

You can get the alignments from blast_record.alignments.
```python
print("Alignments for sequence", blast_record.query)
for alignment in blast_record.alignments:
    print("Accession number:", alignment.accession)
    print("Sequence:", alignment.title)
    print("Length:", alignment.length)
    print()
```
You can list the various attributes of each object using the dir() function.
If you have multiple query sequences, you can parse the result as follows.

```python
blast_records = NCBIXML.parse(result_handle)
# You can use a for loop to access the records as follows.
for blast_record in blast_records:
    print("BLAST result for sequence:", blast_record.query)
    print("Number of alignments:", len(blast_record.alignments))
    print()
```
You can read more about how to use Biopython with BLAST from the Biopython Tutorial and Cookbook.

**Other sequence search tools: SearchIO (QueryResult, Hit...)**

a

---


#### 5. Multiple Sequence Alignment

**Reading a MSA**

aa

**Creating an alignment using different algorithms: ClustalW, MUSCLE...**

We can perform multiple sequence alignment (MSA) where we compare only more than two sequences. If you want to know more about MSA, you can read my article.

Biopython provides command-line wrappers for MSA tools such as Clustal Omega, T-Coffee, ClustalW and DIALIGN.

For example, if you want to use Clustal Omega, then first, you have to download its precompiled binaries. Then you can get an executable command as follows.
```python
from Bio.Align.Applications import ClustalOmegaCommandline
in_file = "sample.fsa"
out_file = "aligned.fasta"
clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
print(clustalomega_cline)
```
clustalomega_cline will be the command which you have to run. You can simply copy-paste it on your terminal.

**Pairwise sequence alignment**

We will be using the Bio.pairwise2 module for PSA.

```python
from Bio import pairwise2

# We have two sequences in two files. Let’s read the sequences from the files.
from Bio import SeqIO
seq1 = next(SeqIO.parse("seq1.fasta", "fasta"))
seq2 = next(SeqIO.parse("seq2.fasta", "fasta"))

# Now we can align the two sequences.

alignments = pairwise2.align.globalxx(seq1, seq2)

# Let’s print out the alignment.

from Bio.pairwise2 import format_alignment
for a in alignments:
    print(format_alignment(*a))
```
You can do a local alignment of the two sequences as follows.
```python
alignments = pairwise2.align.localxx(seq1, seq2)
```
You can do a global alignment and change the scoring scheme (assign custom values for matches, mismatches and gaps). For example, matching characters are given 2 points, 1 point is deducted for each mismatching character. 0.5 points are deducted when opening a gap, and 0.1 points are deducted when extending it.
```python
alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -0.5, -0.1)
```

---


#### 6. Phylogenetics

**Constructing a phylogenetic tree**

Phylogenetic trees represent evolutionary relationships between organisms or genes. We can use the Bio.Phylo module for this.
As an example, we will consider the sequences in a file named as msa.phy that can be found in the official biopython test material for tree construction. Make sure to download the msa.phy file.

```python
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO

# Let’s read the sequences from msa.phy and align them.
aln = AlignIO.read('msa.phy', 'phylip')
# We have to get the distance matrix.
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)
# Now we can construct the phylogenetic tree.
constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)
# We can draw the phylogenetic tree on the terminal.
Phylo.draw_ascii(tree)
```

**Modifying an existing tree**

aaa
   
---


#### 7. Sequence motif analysis

aaa

---


#### 8. PDB: 3D structure protein analysis

**Count atoms in a PDB structure**

aaa

---

