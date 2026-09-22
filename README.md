## Proposed analysis pipeline
Phase 1: dataset identification and metadata <source sourceId="302-309" title="" sourceType="library" />

1. Search curated human microbiome resources, SRA and ENA for tonsil tissue, saliva, oral rinse, plaque, buccal, tongue and oropharyngeal shotgun metagenomes.

2. Record study accession, sample ID, body site, disease state, age, sex, smoking, antibiotic exposure, sequencing platform, read type, read length, depth and host-depletion method.

3. Exclude datasets without adequate metadata or with unsuitable sequencing depth for the proposed analysis.

4. Separate tonsil tissue from other oral niches rather than pooling them immediately.

Phase 2: read processing and profiling <source sourceId="311-317" title="" sourceType="library" />

1. Run FastQC and adapter or quality trimming.

2. Remove human reads using an appropriate host-reference procedure and document the reference version.

3. Profile taxonomy using MetaPhlAn 4.

4. Record sample-level read counts, retained bases, taxonomic profiles and batch variables.

5. Apply prevalence and coverage filters before comparing candidate taxa or genes.

Phase 3: assembly and functional annotation <source sourceId="319-326" title="" sourceType="library" />

1. Assemble reads with MEGAHIT, initially using individual assemblies for straightforward interpretation.

2. Consider co-assembly only for related samples with sufficient coverage and compatible biological design.

3. Predict genes and proteins from contigs.

4. Annotate with eggNOG-mapper, KofamKOALA and targeted enzyme databases.

5. Search specifically for Lcp, RoxA, RoxB, oxidoreductase, esterase and hydrolase candidates using homology and profile-HMM approaches.

6. Track whether each candidate is assembled from a contig, located in a MAG, supported by coverage and assigned to a taxon.

Phase 4: MAG reconstruction <source sourceId="328-336" title="" sourceType="library" />

1. Map reads back to contigs to estimate coverage.

2. Bin contigs with MetaBAT 2 and consider a second binner for comparison.

3. Dereplicate bins with dRep or an equivalent genome-comparison method.

4. Assess completeness and contamination with CheckM.

5. Classify acceptable MAGs with GTDB-Tk.

6. Annotate each MAG and link candidate enzyme genes to MAG, taxonomy and sample prevalence.

7. Report MAG quality according to MIMAG-style fields.

Phase 5: candidate enzyme ranking <source sourceId="338-348" title="" sourceType="library" />

1. Build the positive reference set from experimentally characterised enzymes.

2. Keep Lcp, RoxA and RoxB as separate enzyme families rather than assuming they are interchangeable.

3. Generate candidate sequences from all oral and tonsil metagenomes and MAGs.

4. Remove fragments, duplicates, host proteins and sequences with inadequate quality.

5. Run homology, profile-HMM and domain searches before ML.

6. Generate sequence embeddings and train simple baseline classifiers before more complex models.

7. Use strict separation of training and test data to reduce homolog leakage.

8. Use masked-reference recovery, cross-validation and negative controls to test whether the model can recover known enzymes.

9. Rank candidates using model score, uncertainty, sequence evidence, active-site evidence, structural plausibility, taxonomic origin, MAG quality and sample prevalence.
