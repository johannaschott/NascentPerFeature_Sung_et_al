## Feature-wise estimation of nascent RNA from T-to-C transitions

### Overview:
Nascent RNA labelled with the modified nucleotide 4-Thiouridine (4sU) can be detected by converting 4sU chemically into cytidine.[^1] Because only a minor percentage (typically 2 - 12%) of uridines is replaced by 4sU, not all nascent fragments will display a T-to-C transition, although the conversion is nearly complete. The size of this "invisible" nascent fraction depends on the incorporation rate and the number of Us per fragment. For an incorporation rate of 2% and 10 Us per fragment (p = 0.02, n = 10), ~ 82% of fragments will not show any T-to-C transitions. For an incorporation rate of 10% and 40 Us (p = 0.1, n = 40), only ~1.5% of reads will not show any transitions.

![binom1](https://user-images.githubusercontent.com/37538623/235510082-b82756c5-270f-4349-b904-18f66d959d61.png)

In addition, non-nascent fragments will show a minor rate of T-to-C transitions due to errors during reverse transcription or sequencing ($p_{bg}$), which also contribute to the observed distribution of T-to-C transition counts. Therefore, the proportion of nascent reads ($\pi$) can be estimated from a binomial mixture model[^2]:

![binom2](https://user-images.githubusercontent.com/37538623/235511542-efab4876-92d1-4f09-a68f-bdc90050d99f.png)

In <Title of mansucript>, we estimated these parameters separately for intronic and spliced fragments as well as for regulatory groups of genes.[^3]  

### Steps
* Alignment to the genome using STAR[^4]
* Identification of SNPs (from an external set of sequences; theoretically, this can be achieved from the same data, because SNPs should lead to a much higher T-to-C transition rate than 4sU incorporation)
* Removal of reads that overlap putative SNPs
* Identification and annotation of intronic and exon-exon junction reads with featureCounts[^5]
* Feature-wise counting of T-to-C transitions (i.e. at the gene-level)
* Estimation of parameters (transition probability and proportion of nascent reads, background transition rate within non-nascent reads) for groups of genes using non-linear regression in R 

### Tools
In <Title of mansucript>, we used the following tools:
* STAR v2.5.3a
* samtools v???
* featureCounts v???
* R v4.0.5
  
### References
