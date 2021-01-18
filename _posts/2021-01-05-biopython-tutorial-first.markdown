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
   * Finding codon in a sequence
   * Identifying open reading frames
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

The latest version available when I’m writing this article is biopython-1.77 released in May 2020.
You can install Biopython using pip
Biopython requires NumPy which will be installed automatically if you install Biopython with pip (see below for compiling Biopython yourself).

```python
pip install biopython
pip install --upgrade biopython
```

You can test whether Biopython is properly installed by executing the following line in the python interpreter.
```python
import Bio
```

If you get an error such as ImportError: No module named Bio then you haven’t installed Biopython properly in your working environment. If no error messages appear, we are good to go.


---


#### 2. Working with sequences

**Creating a sequence**
To create your own sequence, you can use the Biopython Seq object. Here is an example.

```python
>>> from Bio.Seq import Seq
>>> my_sequence = Seq("ATGACGTTGCATG")
>>> print("The sequence is", my_sequence)
The sequence is ATGACGTTGCATG
>>> print("The length of the sequence is", len(my_sequence))
The length of the sequence is 13
```

**Get the reverse complement of a sequence**
You can easily get the reverse complement of a sequence using a single function call reverse_complement().
```python
>>> print("The reverse complement if the sequence is", my_sequence.reverse_complement())
The reverse complement if the sequence is CATGCAACGTCAT
```

**Count the number of occurrences of a nucleotide**
You can get the number of occurrence of a particular nucleotide using the count() function.
```python
>>> print("The number of As in the sequence", my_sequence.count("A"))
The number of As in the sequence 3
```

**Find the starting index of a subsequence**
You can find the starting index of a subsequence using the find() function.
```python
>>> print("Found TTG in the sequence at index", my_sequence.find("TTG"))
Found TTG in the sequence at index 6
```

**Reading a sequence**
Biopython’s SeqIO (Sequence Input/Output) interface can be used to read sequence files. The parse() function takes a file (with a file handle and format) and returns a SeqRecord iterator. Following is an example of how to read a FASTA file.
```python
from Bio import SeqIO
for record in SeqIO.parse("example.fasta", "fasta"):
    print(record.id)
record.id will return the identifier of the sequence. record.seq will return the sequence itself. record.description will return the sequence description.
```

**Writing sequences to a file**
Biopython’s SeqIO (Sequence Input/Output) interface can be used to write sequences to files. Following is an example where a list of sequences are written to a FASTA file.
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

**Convert a FASTQ file to FASTA file**
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


1. Running Web BLAST
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

2. Reading a BLAST result
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

3. Pairwise sequence alignment

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

4. Multiple sequence alignment
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

5. Construct a phylogenetic tree
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

