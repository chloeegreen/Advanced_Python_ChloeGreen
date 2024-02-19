# Advanced_Python_ChloeGreen
This is a final portfolio with my Advanced Python information.

##Sequence Objects (videos 1-4)
```python
from Bio.Seq import Seq
```


```python
# this can be used to name a sequence that wants to be evaluated
my_seq = Seq("GATCG")
```


```python
# a number can be assigned to every letter in the sequence
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# We can also print the length of each sequence
print(len(my_seq))
```

    5



```python
# we can recall the specific letter that matches the number in the sequence (following two steps)
print(my_seq[0])
```

    G



```python
print(my_seq[4])
```

    G



```python
print(my_seq[2])
```

    T



```python
# we can count how many times a specific segment occurs in a sequence, but this does not include overelapping
Seq("AAA").count("AA")
```




    1




```python
my_seq = Seq("GATCGATGGGCCTATAGGATCGAAAATCGC")
```


```python
# this gives the number of letters within the sequence 
len(my_seq)
```




    30




```python
# can count the number of a specific letter 
my_seq.count("G")
```




    9




```python
# this can be used to find a specific count of letters
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    50.0




```python
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GATYCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
# this can be used to slice sequences to give a fraction of the letters we are looking at
gc_fraction(my_seq)
```




    0.46875




```python
my_seq[4:12]
```




    Seq('CGATGGGC')




```python
# every third nucleotide can be identified or any value that needs to be evaluated
my_seq[0::3]
```




    Seq('GYAGCTAAGAC')




```python
my_seq[1::3]
```




    Seq('ACTGTAGTAAG')




```python
my_seq[2:3]
```




    Seq('T')




```python
# by using a negative, the sequence can be reprinted backwards
my_seq[::-1]
```




    Seq('CGCTAAAAGCTAGGATATATCCGGGTAGCYTAG')




```python
# this command turns it into its original string
str(my_seq)
```




    'GATYCGATGGGCCTATATAGGATCGAAAATCGC'




```python
# a name can be given to tell what the fasta string was
fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
print(fasta_format_string)
```

    >Name
    GATYCGATGGGCCTATATAGGATCGAAAATCGC
    



```python
# multiple sequences can be identified at a time
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
# the identified sequences can be added together
seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
seq2 + seq1
```




    Seq('AACCGGACGT')




```python
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
```


```python
# N could not be found within the sequence and can serve as a spacer
spacer = Seq("N" *10)
```


```python
# this will take the spacer object and joins all of the contigs together
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')




```python
dna_seq = Seq("acgtACGT")
```


```python
dna_seq
```




    Seq('acgtACGT')




```python
# this changes everything in the sequence to upper case
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
# this changes everything in the sequence to lower case
dna_seq.lower()
```




    Seq('acgtacgt')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
# this proves that searching is case sensitive
"gtac" in dna_seq
```




    False




```python
"GTAC" in dna_seq
```




    False




```python
dna_seq
```




    Seq('acgtACGT')




```python
# this assigns the dna_seq name to the sequence
dna_seq = dna_seq.upper()
```


```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GATCGATGGGGCCTATATAGGATCGAAAATCGC")
```


```python
# this gives the complement strand of the seqyence
my_seq.complement()
```




    Seq('CTAGCTACCCCGGATATATCCTAGCTTTTAGCG')




```python
# this does the reverse while using the negative form together
my_seq.reverse_complement()
```




    Seq('GCGATTTTCGATCCTATATAGGCCCCATCGATC')




```python
# this assigns both a sequence to the protein and its complement
protein_seq = Seq("EVYRAA")
protein_seq.complement()
```




    Seq('EBRYTT')




```python
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
# this is the reverse complement of the template DNA strand
template_dna = coding_dna.reverse_complement()
```


```python
template_dna
```




    Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')




```python
# this converts the coding dna from t's to u's to make it messenger 
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# this creates the same affect as what is shown on the messenger rna
template_dna.reverse_complement().transcribe()
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# this is reverse transcription that turns it into the original coding strand
messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# the astriks are stop codons, but this makes the previous strand into single letter codons
messenger_rna.translate()
```




    Seq('MAIVMGR*KGAR*')




```python
# this takes out the premature stop codon represented by the atrik above
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('MAIVMGRWKGAR*')




```python
# the table referenced before also may have a specific number to give the same result
coding_dna.translate(table = 2)
```




    Seq('MAIVMGRWKGAR*')




```python
# this causes it to stop where the stop codon was marked by the astrik previously
coding_dna.translate(to_stop = True)
```




    Seq('MAIVMGR')




```python
# the astrik is taken away because of the to_stop
coding_dna.translate(table = 2, to_stop=True)
```




    Seq('MAIVMGRWKGAR')




```python
# the symbol for stop can be changed if directed to
coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('MAIVMGRWKGAR!')




```python
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


```python
# this translate the table to a bacterial one like we did with the Vertebrate mitochondrial before
gene.translate(table = "Bacterial")
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HR*')




```python
# this publishes it till it reaches the bacterial stop codon 
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
# this tells Biopython that it is a complete gene and it should start with methionine
gene.translate(table = "Bacterial", cds = True)
```




    Seq('MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
from Bio.Data import CodonTable
```


```python
# we have now named the table as standard
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
# for reference, this is table 2 and that could be used instead of Vertebrate Mitochondrial. This is a function from the codon table 
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
# the table can be visualized
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
# now we can see the mitochondrial table
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
# you can ask for the specific stop codons for the specific table you want to read
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
# you can do the same with start
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
# now we start sequence comparison
seq = Seq("ACGT")
```


```python
"ACGT" == seq1
```




    True




```python
# the name of the sequence is equivalent to the content of the sequence 
seq1 == "ACGT"
```


```python
# this defines the sequence as 10 base pairs long
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
# this answers how long the sequence is; helpful in comparing gene sizes
len(unknown_seq)
```




    10




```python
seq = Seq({117512683: "TTGAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


```python
# this tells us that we have a length of 20, but don't know what those 20 are
seq[1000:1020]
```




    Seq(None, length=20)




```python
# this gives us an answer since we previously defined it 
seq[117512690:117512700]
```




    Seq('CTGAATGTGA')




```python
# this recalls the sequence through the end while also giving the length; shows you can have partially defined sequences
seq[117512670:]
```




    Seq({13: 'TTGAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)




```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length = 10)
```


```python
# this tells us the sequence information we do know along with information concerning what we don't know to add it together
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```


```python
# importing this allows for a position in a sequence to be changed 
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# we are saying that we want position 5 in our defined sequence to become a C
mutable_seq[5] = "C"
```


```python
# this displays the changed sequence
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# this reverses the defined sequence
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
# this changes the sequence back to being protected and unable to be changed
new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
```


```python
# this gives the reversal an the complement of the sequence
reverse_complement(my_string)
```




    'CTAACCAGCAGCACGACCACCCTTCCAACGACCCATAACAGC'




```python
# this turns it into an RNA strand
transcribe(my_string)
```




    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUUAG'




```python
# this reverts it to its original before the transcription
back_transcribe(my_string)
```




    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG'




```python
# this turns it into the amino acid version
translate(my_string)
```




    'AVMGRWKGGRAAG*'

