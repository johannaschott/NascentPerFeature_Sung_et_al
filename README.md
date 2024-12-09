## Feature-wise estimation of nascent RNA from T-to-C transitions

This repository provides access to the custom code used in "Stress-induced nuclear speckle reorganization is linked to activation of immediate early gene splicing" for the analysis of splicing efficiency within nascent and pre-existing RNA during ribotoxic stress.<sup>1</sup>

### Overview:
Nascent RNA labelled with the modified nucleotide 4-Thiouridine (4sU) can be detected by converting 4sU chemically into cytidine.<sup>2</sup> Because only a minor percentage (typically 2 - 12%) of uridines is replaced by 4sU, not all nascent fragments will display a T-to-C transition, although the conversion is nearly complete. The size of this "invisible" nascent fraction depends on the incorporation rate and the number of Us per fragment. For an incorporation rate of 2% and 10 Us per fragment (p = 0.02, n = 10), ~ 82% of fragments will not show any T-to-C transitions. For an incorporation rate of 10% and 40 Us (p = 0.1, n = 40), only ~1.5% of reads will not show any transitions.

![binom1](https://user-images.githubusercontent.com/37538623/235510082-b82756c5-270f-4349-b904-18f66d959d61.png)

In addition, non-nascent fragments will show a minor rate of T-to-C transitions due to errors during reverse transcription or sequencing ($p_{bg}$), which also contribute to the observed distribution of T-to-C transition counts. Therefore, the proportion of nascent reads ($\pi$) can be estimated from a binomial mixture model<sup>3</sup>:

![binom2](https://user-images.githubusercontent.com/37538623/235511542-efab4876-92d1-4f09-a68f-bdc90050d99f.png)

In "Stress-induced nuclear speckle reorganization is linked to activation of immediate early gene splicing", we estimated these parameters separately for intronic and spliced fragments as well as for regulatory groups of genes in order to draw conclusions about splicing efficiency within nascent and pre-existing RNA during ribotoxic stress.<sup>1</sup>  

### Steps
* Alignment to the genome using STAR<sup>4</sup>
* Identification of SNPs (from an external set of sequences; theoretically, this can be achieved from the same data, because SNPs should lead to a much higher T-to-C transition rate than 4sU incorporation)
* Removal of reads that overlap putative SNPs
* Identification and annotation of intronic and exon-exon junction reads with featureCounts<sup>5</sup>
* Feature-wise counting of T-to-C transitions (i.e. at the gene-level)
* Estimation of parameters (transition probability and proportion of nascent reads, background transition rate within non-nascent reads) for groups of genes using non-linear regression in R 

### Tools
In "Stress-induced nuclear speckle reorganization is linked to activation of immediate early gene splicing", we used the following tools:
* STAR v2.5.3a
* samtools v1.7
* featureCounts v1.5.2
* R v4.0.5
  
### Preparation
If you want to recapitulate the full analysis presented in Sung et al. 2023, you have to download the following additional files:
* Fastq files from GEO (GSE231520) (place them into the directory ./raw_data)
* Genome Reference Consortium Human Build 38 patch release 10 (GRCh38.p10), as a fasta file in ./genome_files/GRCh38.p10.genome.fa
In addition, the gtf files provided in ./genome_files have to be de-compressed.
  
### References
1. Sung HM, Schott J, Boss P, Lehmann JA, Hardt MR, Lindner D, Messens J, Bogeski I, Ohler U, Stoecklin G. Stress-induced nuclear speckle reorganization is linked to activation of immediate early gene splicing. J Cell Biol. 2023 Dec 4;222(12):e202111151.
2. Schott J, Reitter S, Lindner D, Grosser J, Bruer M, Shenoy A, Geiger T, Mathes A, Dobreva G, and Stoecklin G. 2021. 'Nascent Ribo-Seq measures ribosomal loading time and reveals kinetic impact on ribosome density', Nat Methods, 18: 1068-74.
3. Jurges C, Dolken L, and Erhard F. 2018. 'Dissecting newly transcribed and old RNA using GRAND-SLAM', Bioinformatics, 34: i218-i26.
4. Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, and Gingeras TR. 2013. 'STAR: ultrafast universal RNA-seq aligner', Bioinformatics, 29: 15-21.
5. Liao Y, Smyth GK, and Shi W. 2014. 'featureCounts: an efficient general purpose program for assigning sequence reads to genomic features', Bioinformatics, 30: 923-30.
