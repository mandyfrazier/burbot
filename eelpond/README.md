# Burbot RNAseq analysis using eelpond

We are going to use [eelpond](https://dib-lab.github.io/eelpond/) to build an annotated transcriptome from Burbot RNAseq data, and then map our reads back to that transcriptome for quantification and analysis.
Building a transcriptome from all 40 samples would take a great deal of time and compute power. Further, after a certain point, adding more data simply introduces more error and makes the assembly more difficult rather than better. As such, we're going to do this analysis in two parts. First, we'll build and annotate a transcriptome using a subset of our samples, then we'll run a pipeline to align all the samples to it.

This is the overall workflow for this project:

1. Trim reads and assemble a subset of four individuals using trimmomatic, trinity
2 Annotate using dammit (if you use some number other than 4 individuals, change line 35 of 1_burbot_assemble.yaml to the right number)
3. Trim all reads to prep for Salmon via eelpond
4. Quantify expression levels using Salmon via eelpond
5. Determine differentially expressed genes using DESeq2 via eelpond

This should only require running two overall workflows though, one for the subset of data, and one for the overall data.

## Detailed steps and notes

### 1 Trim, assemble, annotate

This requires two files: my_transcriptome_samples.tsv and 1_burbot_assemble.yaml 
 - my_transcriptome_samples.tsv:  a tab separated file that contains all the relevant information for the subset of samples that will be used in the assembly
 - 1_burbot_assemble.yaml: a yaml file detailing the programs to be used and the specific parameters for each program
 
#### To run

You should probably run this in a screen session so it can run while you're away:

`screen`

`cd eelpond/`

`conda activate eelpond`

`./run_eelpond 1_burbot_assemble.yaml assemble annotate --threads 46`

then *briefly* hit: 'control' and 'r' at the same time, then the letter 'd' to get out of your screen session.

When you come back to check on it, you re-open your screen session by typing 

`screen -r`

### 2 Trim and prep all reads

This requires two files: my_transcriptome_samples.tsv and 1_burbot_assemble.yaml 
 - all_my_samples.tsv:  a tab separated file that contains all the relevant information for all of samples that will be used in the assembly
 - 2_burbot_quantify.yaml: a yaml file detailing the programs to be used and the specific parameters for each program

#### To run
./run_eelpond 2_burbot_quantify.yaml assemblyinput quantify diffexp


## Run notes and important details for paper


In this protocol we are running trimmomatic with *just* Illuminaclip, because we're starting with reads that were already trimmed and Rcorrected, but the pipeline will break if we skip the trimming entirely. So, we're doing a trim step to generate the intermediate files, but the trim step won't actually do any trimming. The trim parameters used for the real trimming step were: sliding window 4:5, leading:5, trailing:5, minlen:25

For the assembly, we're using four individuals from the same family (Ma) that are all non-cannibals. I chose to use the non-cannibals because this is the more "normal" feeding strategy and a lake-adapted family because that's more representative of all the samples we have. I used the Ma family arbitrarily, it was the first family in the list. Reducing the genetic variability and using four individuals instead of 40 will help reduce the number of total transcripts and will make the assembly better because there should be fewer mutations. 

After assembly, I compared the 40-sample transcriptome to the 4-sample transcriptome using ./TrinityStats, sourmash and BUSCO. The number of total trinity transcripts decreased from 793,202 to 387,175. 


### Fix for 'real' runthrough

- no thread calls in the yaml, do it in the origianl eelpond call
- just `experiment1` not `_experiment1`

## Sourmash: 

Build a signature for the large transcriptome: 

`sourmash compute --scaled 10000 -k 31 ~/data/3_Trinity/burbot_trinity.fasta -o burbot_large.sig`

Build a signature for the small transcriptome: 

`sourmash compute --scaled 10000 -k 31 /home/ubuntu/data/burbot_assemble_out__experiment1/assembly/burbot_assemble_trinity.fasta -o burbot_small.sig`

Compare transcriptomes:

`sourmash search -k 31 burbot_large.sig burbot_small.sig --containment`

Result was 54.8% similarity (54.8% of the 40-sample transcriptome is contained in the 4-sample transcriptome). This is only comparing the k-mers, it doesn't tell you which transcriptome is more complete. 

## BUSCO: 

40-sample transcriptome: 

Download eukaryota dataset
`wget https://busco.ezlab.org/datasets/eukaryota_odb9.tar.gz
gunzip eukaryota_odb9.tar.gz
tar -xvf eukaryota_odb9.tar`

Run BUSCO
`run_BUSCO.py -i Trinity.fasta -o burbot_busco_eukaryota \
-l eukaryota_odb9 -m transcriptome --cpu 4`

Results: C:97.4%[S:24.8%,D:72.6%],F:2.6%,M:0.0%,n:303 
Running time: <1 hr 

Download fish dataset
`wget https://busco.ezlab.org/datasets/actinopterygii_odb9.tar.gz
gunzip actinopterygii_odb9.tar.gz
tar -xvf actinopterygii_odb9.tar`

Run BUSCO
`run_BUSCO.py -i Trinity.fasta -o burbot_busco_fish -l actinopterygii_odb9 -m transcriptome --cpu 4`

Results: C:81.0%[S:27.9%,D:53.1%],F:13.6%,M:5.4%,n:4584
Running time: ~3 days 



4-sample transcriptome: 

`run_BUSCO.py -i /home/ubuntu/data/burbot_assemble_out__experiment1/assembly/burbot_assemble_trinity.fasta -o burbot2_busco_eukaryota -l eukaryota_odb9 -m transcriptome --cpu 6`

Output: C:99.0%[S:36.0%,D:63.0%],F:0.0%,M:1.0%,n:303   
Time to compute: 20 minutes 

`run_BUSCO.py -i /home/ubuntu/data/burbot_assemble_out__experiment1/assembly/burbot_assemble_trinity.fasta -o burbot2_busco_fish -l actinopterygii_odb9 -m transcriptome --cpu 7`

Results: C:91.9%[S:36.9%,D:55.0%],F:4.8%,M:3.3%,n:4584 
Time to compute: ~6 hours 


So, the BUSCO score of the 40-sample transcriptome against actinopterygii was 81%, but it was 92% for the 4-sample transcriptome. This tells us that the loss of transcripts in Trinity for the 4-sample was mostly error-filled transcripts! 


## Code to parse dammit output 

The first time we ran the differential expression, we were using trinity contig numbers. But this included all of the different isoforms as separate contigs, so the count was super big (no collapsing). The second time, we ran it using the trinity-identified "genes", which collapsed based on what Trinity thinks are likely genes. However, what we really want to do is collapse genes based on the dammit annotation from the Ensembl Atlantic cod transcriptome. The intitial GFF file from dammit includes all of the different gene "options" for each transcript that it found, but we need to assign it what we think is the most likely annoation. There are different metrics you can use to choose which annotation we assign to each transcript, we chose based on the longest match length. This code is how we splice the output from dammit to generate one annotation per transcript. 

Are the reads them collapsed in tximport? 

```
gff_file = "~/data/burbot_assemble_out_run2/annotation/burbot_assemble_trinity.fasta.dammit/burbot_assemble_trinity.fasta.dammit.gff3"
annotations = GFF3Parser(filename=gff_file).read()
# Keeps track of long the annotation is 
annotations["length"] = annotations["end"].subtract(annotations["start"], fill_value=0)
# make new table for each with seqid, Name, start, end, length
annotations = annotations.loc[annotations['database'] == "Edit_Gadus_morhua.gadMor1.pep.all.fa"]
annotations = annotations.sort_values(by=['seqid','length'],ascending=False).drop_duplicates(subset='seqid')[['seqid', 'Name','start','end','length']]
annotations = annotations.rename(columns = {'Name':'Ensembl'})
print('ensembl annotations',annotations.shape)
new_file = annotations.dropna(axis=0,how='all')
new_file.head()

# to rename in transmap file: 
conversion = pd.read_csv("~/data/burbot_assemble_out_run2/annotation/burbot_assemble_trinity.fasta.dammit/burbot_assemble_trinity.fasta.dammit.namemap.csv")
conversion['contig'], conversion['info'] = conversion['original'].str.split(' ', 1).str
conversion['seqid'] = conversion['renamed']
conversion = conversion[['contig','seqid']]
contigs_w_expression_conversion = pd.merge(conversion, new_file, on="seqid")
contigs_w_expression_conversion = contigs_w_expression_conversion[['Ensembl','seqid']]
contigs_w_expression_conversion.to_csv("~/data/burbot_assemble_out_run2/annotation/burbot_assemble_trinity.fasta.dammit/burbot_assemble_Ensembl.fasta.dammit.namemap.csv", index=False)
```






