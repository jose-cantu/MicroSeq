---
layout: page 
title: Gui Walkthrough # <title> element shown in browser tab 
permalink: /gui/
--- 

## Launch 

* Note in order to launch the GUI you need to have conda activated and from there you type in the  terminal `microseq-gui`. 

* The window will show the defaults:
Control	Default	Meaning
DB	gg2	Reference database
ID %	97	BLAST identity cut-off
Q-cov %	80	Query coverage cut-off
Threads	4	CPU threads
Mode	Fast (megablast)	Alignment algorithm


## Select Input 
* First you click Browse...and pick in this case a folder with ab1 files. In order to access folder selection it's impertative you click cancel first then choose folder will pop up for it. Then you can go down into the nested folders until you find the desired folder which directly in that file houses the .ab1 files you want to blast against. You can also type in the path in the search bar as well. 

## Configuring The Run 
* In a future update I will have all the dials adjustable for now they are the default settings so keep that in mind. 
* If you want the BIOM file you will need to tick the BIOM box when clicking full pipeline run note you will need to supply a metadata table (GUI will warn if missing). 
* Choose Fast(megablast) for routine Sanger reads; Comprehensive (blastn) for divergent amplicon reads. Note I do have it set up where if the coverage or percent ID doesn't reach 90% blastn will rerun instead given the nature of megablast algorithm.(Fast = NCBI Megablast—large 28-bp word, tuned for ≥95 % identity hits between nearly identical Sanger traces; Comprehensive = standard blastn—11-bp word, slower but ∼5× more sensitive to <90 % identity or divergent/gapped matches; MicroSeq auto-falls back to blastn whenever the initial Megablast hit covers <90 % of the read or shows <90 % identity). 
> Current limitation (to be patched): > The Comprehensive (blastn) path is still being refactored; a rare QPainter crash can occur after the run finishes. For production runs keep the default Fast (megablast), which is stable. An updated patch is under review in a separate branch that I am testing for stability and will be shippedin the next release. 

## Execution What do the buttons do??

* Button Stages Performed 
* **Run QC** Ab1 -> Fastq -> trim -> pass/fail 
* **Run Blast** BLAST only existing fasta/contigs 
* **Full Pipeline** (BIOM unchecked) QC + BLAST + Taxonomy (Note:Assembly for full length will be added in an update so I may just include a checkbox for partial 16s thenclicking FullPipeline as a workaround.)
* **Post-BLAST** Collapse hits to best one + merge metadata to BIOM file that is outputted 

During execution the status bar shown the current stage and ETA; log pane streams MicroSeq's standard logging I have embedded in it. 

## Example of a Pilot Run Using Data From Quintara 
I unzipped this file here:
`1038530_195270`. 

I then walked through the file system using MicroSeq GUI until I chose that folder specifically that has the ab1 files I need directly inside it and no other folders.I exectued in clicking Full Pipeline it will go through the stages described above.
It will output the same name with 'microseq' appended to it so you know MicroSeq was being used for the run and the results are in the directory. 

I will go through how the tree structure looks and the run...

```
1038530_195270_microseq/
├── raw_ab1/              # original chromatograms
├── raw_fastq/            # AB1→FASTQ before trimming
├── passed_qc_fastq/
├── failed_qc_fastq/
├── qc/
│   ├── *_avg_qual.txt    # Individual output metric avg Q for each sequence 
│   ├── trim_summary.tsv  # Summary output metrics for avg Q for each file combined
│   └── trimmed.fasta     # Placed here for providence - same file as reads.fasta
├── reads.fasta           # merged QC-passed reads
├── hits.tsv              # BLAST results
└── hits_tax.tsv          # taxonomy-joined results
```
> Note: duplicate FASTA files > `qc/trimmed.fasta` and `reads.fasta` contain the exact same set of QC-passed reads. The first is emitted by the Trim step and kept inside qc/ for provenance; the second is copied to the top-level so downstream tools (assembly, BLAST, etc.) can find a single “canonical” FASTA without digging into sub-folders.

* raw_ab1/ inputs retained for provenance (copies of original ab1 files) will include option to symlink later. 
* raw_fastq/ FASTQ files produced from the ab1 files before trimming 
* qc/ per-read stats and trimmed FASTA quality‑control summary files:
  - *_avg_qual.txt is the per‑read length and Phred averages.
  - trim_summary.tsv is the combined metrics for all files.
  - trimmed.fasta is the final FASTA used for assembly or BLAST that you can manually hand off to NCBI website for example.
* passed_qc_fastq/ quality filter outcome with default at (≥ Q20); these are individual FASTQ files that passed that metric. 
* failed_qc_fastq/ quality filter outcome with default at (≥ Q20); these are individual files that failed that metric whose average quality was below the threshold. 
* reads.fasta the merged FASTA of all QC‑passed reads (made from passed_qc_fastq/). Created during the pipeline when FASTQs are converted to a single FASTA that is used for blasting by MicroSeq. 
* hits.tsv `microseq blast` output 
* hits_tax.tsv This file has the taxonomy needed to interpret the blast results by adding the taxonomy column using `microseq add_taxonomy` 



