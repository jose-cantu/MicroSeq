# Project 

Here are each of the test driven scripts for this project along with the data repository of 2 examples of ab1 files {One featuring reverse only and the other being a paired example of F/R ab1 files}.

Meant to be for demonstrative use paired with the doc where I go into more detail how to use them.

Features and Data directories are used as post hoc examples for how to run and set up the post-blast run to get the BIOM file.

tests/data/ – holds a tiny synthetic FASTA plus metadata that let the post‑BLAST pipeline run in tests without needing large datasets. The post‑BLAST integration test explicitly wires demo.fasta and meta.tsv into a full microseq blast → add_taxonomy → postblast round‑trip to assert BIOM/CSV outputs are produced. Look at `test/test_pipeline.py` for details on the code test there for this. 

tests/fixtures/ – contains a handful of real Sanger AB1 traces (and supporting TSVs) used to smoke‑test the AB1 workflow. For example, the “single AB1” full‑pipeline test reads fixtures/39764260.ab1, stubs out BLAST/ taxonomy, and verifies that the pipeline still emits a non‑empty reads.fasta, confirming the trimming/QC path works on genuine chromatograms.

Look at `tests/test_full_pipeline_single_ab1.py` for details on this. 


