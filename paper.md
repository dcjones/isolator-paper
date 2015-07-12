
---
title: 'Isolator: accurate and consistent analysis of RNA-Seq experiments'
author:
  - Daniel C. Jones
  - Kavitha T. Kuppusamy
  - Xinxia Peng
  - Michael G. Katze
  - Hannele Ruohola-Baker
  - Walter L. Ruzzo
abstract: |
    While RNA-Seq has enabled great progress towards the goal of wide-scale
    isoform-level mRNA quantification, short reads have limitations when
    resolving complex and similar sets of isoforms.  As a result, estimates of
    isoform abundance, compared to those at the gene level, carry far more
    uncertainty. When confronted with this uncertainty, commonly used methods
    produce estimates that are unstable; small perturbations in the data often
    produce dramatically different results. We introduce a new method, Isolator,
    which analyzes all samples in an experiment in unision using a simple
    Bayesian hierarchical model. Combined with agressive bias correction, it
    produces estimates that are simultaneously accurate and stable.  In a
    comprehensive comparison of accuracy and stability, we show that this
    property is unique to Isolator. We further demonstrate that the approach of
    modeling an entire experiment enables new analyses by examining splicing
    monotonicity across several time points in the devolpment of human
    cardiomyocyte cells.
...


<!--
    Keep in mind
    Abstract: ~150 words.
    Paper: 2500-3500 words.
           <= 6 figures or tables

    Structure
      * Introduction
      * Results
      * Discussion
      * Online Methods
        Use subheadings here, go into details, etc.


What are our six figures?

  * 1. seqc consistency
  * 2. taqmac
  * 3. primepcr
  * 4. spike-in
  * 5. a. batch-effects
  * 5. b. batch-effects, bias-correction
  * 6. rhemac chart

  What about qPCR validation?

-->

# Introduction

Since the introduction of RNA-Seq it has rapidly become a preferred method of
studying gene expression, having proved to be reliable, reproducible, and
increasingly affordable. In principle, RNA-Seq enables measurements of
expression at a higher resolution than has previously been possible or practical
at the same scale. Rather than an often vaguely defined notion of a gene,
expresison of specific isoform or exons can be measured. In practice this
promise is difficult to realize. There are no standardized names or symbols for
isoforms, annotations are often incomplete or in disagreement, and methods
designed to detect changes in splicing rarely concur (TODO: cite).

Because a short read is often compatible with many isoforms, the signal produced
from sequencing must be deconvoluted, requiring a certain level of algorithmic
sophistication. This problem was tackled early on and with great success by
Cufflinks, and indeed many of the methods developed in the interim have followed
from the same basic model of RNA-Seq in which transcripts are represented as
particular distributions over possible reads, and the problem is then to infer
mixing coefficients (relative expression values) inducing a mixture model that
best explains the data.

We built on this work, implementing a new method called Isolator. Though it
retains the same underlying model of RNA-Seq, we approach the problem from the
perspective that isoform-level estimates possess a degree uncertainty that
renders maximum likelihood point-estimates alone a fundamentally inadequate
solution to the problem. As the purpose of gene expression studies is to compare
expression between multiple samples, our goal, beyond maximizing a narrow
definition of accuracy, is to provide a tool that effectively accounts for this
uncertainty in a coherent and reliable manner.

Rather than estimating expression of individual samples, Isolator uses a
hierarchical model of an entire experiment. In addition to infering expression
values for the individual samples, we introduce parameters representing
condition-wise and experiment-wise expression and splicing, as well as
per-transcript variance parameters that are shared across conditions. As
expression studies often use only a small number of replicates, shared variance
parameters allows us to make more accurate estimates of biological variability.

An efficient MCMC algorithm is used to generate samples over the parameters of
this model, which are saved in an HDF5 based format. The output from the sampler
can then be processed to produce point estimates, credible intervals, posterior
probabilities, and diagnostics. Though designed to be run on multiple samples
concurrently, Isolator can be run as a conventional tool on individual samples.

The broad effect of this model is to encode the common assumption that genes
tend to be similarly expressed and similarly spliced between replicates, and to
a lesser extent between conditions. In the absense of sufficient evidence to the
contrary, this informative prior shrinks estimates towards a baseline of no
change, producing more conservative estimates of effect size than one would
obtain by using, for example, a flat Dirichlet prior. This feature of the model
is particularly important for robust treatment of low-expression genes.  With
few reads there is little power to determine splicing patterns. Considering the
samples in isolation will often produce very different estimates.  Statistical
tests based on these estimates, without considering the full posterior
distribution, make approximations by either either disregarding the relative or
compositional nature of the data and assuming independence, or by relying on
asymptotic properties that may not hold when a gene is not deeply sequenced.

In addition to  approach to we attempt to aggressively correct for technical
effects. The naive model of RNA-Seq assumes fragments are sampled uniformly at
random from transcripts in proportion to their abundance. Several have observed
that this uniformity assumption is far from an accurate description of real
sequencing data (TODO: cite). Though severe at the nucleotide resolution, gene
level expression estimates are only moderately affected by this bias, as it
tends to average out over long transcripts. When one considers expression at the
isoform level, attempting to correct for this nonuniformity becomes vital, as
the difference between two isoforms is often only a few nucleotides.

We previously developed an efficient and versatile model to account for sequence
specific bias that commonly occur at the ends of fragments (TODO: cite). In
Isolator, we supplemented this with methods accounting for 3' bias, fragment
GC-content bias, and positional bias caused by end-effects of fragmentation.

Isolator is available under a permissive open source license at:
https://github.com/dcjones/isolator

# Results

The convention for evaluating the accuracy of RNA-Seq analysis methods has
consisted chiefly of comparing point estimate accuracy in simulations,
supplemented with comparisons to qPCR or other technologies. Tuning methods
exclusively to maximize accuracy on simulated data creates the risk of solving
only an idealized mathematical model of RNA-Seq while disregarding noise, bias,
and nonuniformity. As no gold standard for whole-transcriptome isoform-level
estimation exists, we adopt multi-faceted approach, with a focus measuring the
stability or consistency of estimates, along with the accuracy.

Comparisons of gene or transcript expression is sensitive to how one chooses to
measure similarity or distance. Because gene expression typically varies across
several orders of magnitude, measurements like the Pearson correlation, $L_1$ or
$L_2$ distance, or root-mean-square error are often dominated by a small subset
of the most highly expressed genes. Recognizing this, a popular alternative is
the Spearman rank correlation, which conversely can overemphasize the majority
of transcripts which have no, or very low expression. These issues are discussed
in depth by Lovell et al [@Lovell:2015il], who propose measurements of
proportionality as a more meaningful metric. We adopt their "proportionality
correlation", defined as $\text{pcor}(x, y) = 2\text{cov}(\log(x), \log(y)) /
(\text{var}(\log(x)) + \text{var}(\log(y)))$. Similar to other measures of
correlation, proportionality correlation varies from -1 to 1, with 1 indicating
perfectly proportionality, and -1 perfect reciprocality. Zeros are accounted for
by adding 0.1 to all TPM values.

Many methods have been proposed to aid in the analysis of RNA-Seq experiments.
We limit our discussion specifically to those that seek to estimate the relative
abundance a set of known transcripts. These include Cufflinks, RSEM, eXpress,
BitSeq, Sailfish, Salmon, BitSeq, and Kallisto. RSEM was used to produce both
maximum likelihood and posterior mean estimates, which we label “RSEM/ML” and
“RSEM/PM” respectively. This is not an exhaustive account of such methods, but
does represent a wide variety of popular approaches.


## Agreement with qPCR and spike-in controls

Spike-in controls and qPCR are far more limited than RNA-Seq, measuring only the
abundance of specific loci in the case of qPCR, or known proportions of simple
artificial transcripts with spike-in controls. However, both are considered
reliable enough to serve as proxy gold standards for gene-level expression
estimates.

We use data generated by the SEQC project, consisting of four reference samples
(labeled A, B, C, D).  The RNA-seq methods were run using approximately 25
million reads per sample, obtaining estimates that were then compared to qPCR in
Tables TODO and TODO, and the known ERCC spike-in mixtures in Table TODO. In
each of these we find that Isolator produces the highest correlation of the
methods compared, though sometimes by negligible margin. Despite two different
technologies, and two different sets of genes, the qPCR benchmarks agree
closely, ranking the methods nearly identically. Cufflinks, eXpress, Salmon, and
RSEM/ML typically perform similarly while BitSeq, Kallisto, Sailfish, and
RSEM/PM trail slightly. Results are similar with ERCC spike-in controls, the
notable exceptions being Kallisto and Sailfish, which show higher accuracy, and
eXpress, which shows relatively lower accuracy.


Method             A          B          C          D
---------  ---------  ---------  ---------  ---------
Isolator   **0.903**  **0.908**  **0.890**  **0.884**
eXpress        0.901      0.902      0.884      0.874
Cufflinks      0.900      0.901      0.865      0.872
Salmon         0.899      0.897      0.882      0.871
RSEM/ML        0.899      0.896      0.881      0.870
Kallisto       0.893      0.885      0.874      0.860
BitSeq         0.875      0.865      0.863      0.848
Sailfish       0.879      0.857      0.854      0.834
RSEM/PM        0.878      0.866      0.864      0.848

Table: Proportionality correlation between gene level quantification of 806 genes
using TaqMan qPCR and RNA-Seq quantification.


Method             A          B          C          D
---------  ---------  ---------  ---------  ---------
Isolator   **0.878**  **0.866**  **0.839**  **0.852**
Cufflinks      0.870      0.856      0.799      0.841
eXpress        0.870      0.855      0.829      0.840
Salmon         0.866      0.852      0.826      0.836
RSEM/ML        0.865      0.851      0.825      0.835
BitSeq         0.840      0.821      0.802      0.813
Kallisto       0.858      0.840      0.817      0.826
Sailfish       0.844      0.814      0.797      0.802
RSEM/PM        0.840      0.822      0.803      0.811

Table: Proportionality correlation between gene level quantification of 18353 genes
using PrimePCR qPCR and RNA-Seq quantification.


Method             A         B         C         D
---------  --------- --------- --------- ---------
Isolator   **0.979** **0.978** **0.981** **0.982**
Salmon         0.976     0.975     0.978     0.979
Kallisto       0.972     0.972     0.973     0.976
Sailfish       0.970     0.969     0.969     0.972
Cufflinks      0.967     0.969     0.972     0.974
RSEM/PM        0.943     0.949     0.944     0.949
RSEM/ML        0.941     0.948     0.945     0.951
BitSeq         0.940     0.949     0.943     0.946
eXpress        0.931     0.939     0.935     0.942

Table: Proportionality correlation between known proportions of 92 ERCC spike-in
controls and RNA-Seq quantification.


## Estimation consistency

The four SEQC samples consist of two commercially available reference RNA
samples, labeled $A$ and $B$, and two synthetic samples formed by mixing $A$ and
$B$ in specific proportions. Sample $C$ is composed of 75% $A$ and 25% $B$,
while sample $D$ is 25% $A$ and 75% $B$. This design allows us to directly
measure expression consistency. If we estimate transcript abundance for each
sample producing estimates $a, b, c, d$, then for these estimates to be in
agreement we should expect $c$ to be approximated by $0.25a + 0.75b$, and
similarly $d$ by $0.75a + 0.25b$. In this manner we measured consistency for
each method using proportionality correlation (Table TODO).

Results from RSEM are particularly informative. We generated maximum likeilood
and posterior mean estimates using the same aligned reads. Merely switching the
estimate from the former to the latter increased the correlation from 0.922 to
0.968. This stands in contrast to comparisons to qPCR which suggest reduced
accuracy of posterior mean estimates, demonstrating that consistency and
accuracy, while related, are two separate axis of comparison.

Estimates made at the gene level (Supplementary Table TODO) agree with these
results, but are of uniformly higher correlation and show a much smaller gap
between posterior mean and maximum likelihood estimate. Posterior mean methods
(Isolator, BitSeq, RSEM/PM) again show the highest correlation and are in very
close agreement.

Method      $c$ vs $0.75a + 0.25b$   $d$ vs $0.25a + 0.75b$
--------- ------------------------  -----------------------
Isolator                 **0.975**                **0.975**
BitSeq                       0.967                    0.967
RSEM/PM                      0.968                    0.967
Sailfish                     0.932                    0.925
RSEM/ML                      0.922                    0.919
Salmon                       0.916                    0.914
Kallisto                     0.907                    0.902
eXpress                      0.903                    0.899
Cufflinks                    0.870                    0.916

Table: Proportionality correlation between transcript-level estimates for the
mixed samples C and D and weighted averages of estimates for A and B,
corresponding to the mixture proportions for C and D.


Method      $c$ vs $0.75a + 0.25b$   $d$ vs $0.25a + 0.75b$
--------- ------------------------  -----------------------
Isolator                 **0.995**                **0.995**
BitSeq                       0.994                    0.994
RSEM/PM                      0.994                    0.993
RSEM/ML                      0.987                    0.987
Salmon                       0.986                    0.986
eXpress                      0.983                    0.983
Sailfish                     0.978                    0.976
Kallisto                     0.972                    0.970
Cufflinks                    0.941                    0.982

Table: Proportionality correlation between gene-level estimates for the
mixed samples C and D and weighted averages of estimates for A and B,
corresponding to the mixture proportions for C and D. TODO: move this table to
the supplement.


## Batch effects

 <!--![](batch-effect-figure.pdf)-->
\begin{figure*}
\includegraphics[width=\textwidth]{batch-effect-figure.pdf}
\caption{
\textbf{a} A heatmap showing pairwise proportionality correlation between the sample
sampled sequenced on two flowcells each at five sites. Flowcells are numbered 1
or 2 and sequencing sites are abbreviated with three letter codes: Australian
Genome Research Facility (AGR), Beijing Genome Institute (BGI), Cornell
University (CNL), Mayo Clinic (MAY), and Novartis (NVS). \textbf{b} The absolute
change in correlation induced by enabling bias correction on methods that
support it. For clarity this plot excludes points for BitSeq estimates of "MAY
2", as bias correction has an extreme detremental effect on these.
}
\end{figure*}

We further examined the question of consistency by comparing the same samples
sequenced on different flow cells and at different sites, again using data from
the SEQC. We used 10 samples from the larger SEQC dataset, each constituting a single
flowcell lane from two separate flowcells sequenced at five different sites, all
on Illumina HiSeq 2000 instruments, and compared pairwise agreement (Figure TODO).

These correlations are somewhat smaller than those in Table TODO, largely
because of shallower sequencing: there five flowcell lanes were combined,
totaling ~25 millions reads per sample, here only one lane or ~5 million reads
are used. Nevertheless, we again posterior mean methods with showing much higher
agreement between pairs.

When accounting for batch effects, bias correction methods can have a
significant impact. To examine the efficacy of bias correction, we disabled bias
correction functionality on those methods that support it, repeated the
experiment, and measured the change in pairwise correlation (Figure TODO). We
see that bias correction is largely beneficial for each method,

Isolator, Salmon, and Cufflinks show similar improvements with bias correction,
though in a few cases Cufflinks slightly decreases agreement. Kallisto show a
consistent, but slight improvement. BitSeq's bias correction sometimes has a
very positive effect, but other times negative or disastrous effect. Estimates
from the "MAY 2" samples in particular had far worse agreement with other
samples with bias correction enabled.

While this benchmark is informative, it should be considered a lower bound on
the batch effects and bias: these samples were identically prepared and
sequenced.  The fact that technical effects is indicative that RNA-Seq is not so
consistent that batch effects can be under any circumstanced ignored.


## Accuracy in simulated data

To demonstrate Isolator's capacity to produce more accurate estimates by
modeling entire experiments, we simulated a simple RNA-Seq experiment consisting
of two conditions each with three replicated. Expression values were generated
from a two-component log-normal model with parameters fit to Cufflinks estimates
of data from Kuppusamy et al (TODO: cite). We generated simulated RNA-Seq reads
using rlsim (TODO: cite).

In this benchmark, Isolator significantly improves on existing methods. It
has the unique advantage of sharing information between replicates. With the
assumption that genes tend to be similarly spliced between replicates, we are
able to more effectively resolve transcript expression in complex loci.
Kallisto, Salmon, RSEM/ML, and Cufflinks all perform very similarly. While
comparisons to qPCR showed eXpress and Cufflinks with a significant advantage,
effective bias correction is less important in simulated data, plausibly
explaining the difference. Although rlsim models some forms of bias, technical
effects in RNA-Seq are not entirely understood, so these effects are likely
understated.


Method       Correlation
---------    -----------
Isolator           0.919
Kallisto           0.887
Salmon             0.886
RSEM/ML            0.881
Cufflinks          0.881
eXpress            0.825
Sailfish           0.816
RSEM/PM            0.806
BitSeq             0.796

Table: Proportionality correlation between ground truth and estimates produced
by each method on simulated RNA-Seq.


## Agreement between two sequencing technologies

RNA-Seq reads become more informative with length and quantity. We compared data
from the same sample sequenced twice, once with 300pb paired-end reads on a
MiSeq sequencer and again with 100bp paired-end reads using a HiSeq 2000,
yielding approximately 4.5 million and 11 million reads, respectively. To
elucidate the effect of deeper sequencing, we sampled without replacement
smaller subsets of the 11 million HiSeq 2000 reads and measured correlation
between estimates these subsamples with estimates from the MiSeq data (Figure
TODO).

Isolator show the highest correlation, but only when sequencing depth exceeds
several million reads, where it overtakes eXpress by a small margin. The
outliers in the this test are BitSeq and RSEM/PM. The former shows extremely high
correlation with small numbers of reads, and contrary to all other methods,
correlation decreases with deeper sequencing after 1 million reads. When
posterior mean estimates are generated from a model using a uninformative prior,
as is the case with BitSeq, transcripts with few or no reads may produce
estimates that are highly influenced by transcript length, which plausibly
explains this phenomenon. Yet, RSEM/PM does increase in correlation
monotonically with sequencing depth, though the correlation remains very low.
Even using all 11 million reads, RSEM/PM had a correlation of only 0.567.

\begin{figure}
\includegraphics[width=0.5\textwidth]{rhemac-figure.pdf}
\end{figure}


## Finding monotonic differential splicing

By modeling the entire RNA-Seq experiment and saving samples generated during
MCMC sampling, Isolator is uniquely able to compute posterior probabilities
corresponding to arbitrarily complex questions, within the confines of the
model. As a case study, we previously used Isolator to analyze splicing dynamics
in human cardiomyocyte cells during maturation.

The experiment consisted of a number of conditions at different time points so
did not easily lend itself to analysis using pairwise hypothesis tests. We
instead looked for genes showed a change in splicing that was consistent between
mature and immature cells by computing a probability of "monotonic splicing", or
specifically that an observed change in splicing occurs in a consistent
direction

TODO: talk about validation.


# Discussion

To assess the suitability of Isolator we have compared it to large number of
alternative methods using a wide variety of benchmarks, with a particular focus
on estimation consistency or stability. Data generated by the SEQC is ideally
suited to examine the question consistency. Our data agrees with the results
published by the SEQC, who compared Cufflinks and BitSeq and found BitSeq to
produce more consistent estimates. Here we control for alignment method by using
the same alignments for every method (excluding the alignment-free methods
Salmon, Sailfish, and Kallisto), and compare wider variety methods. In this
broader context it becomes clear that this separation is a consequence of the
choice in estimators. Unmoderated maximum likelihood methods, though often
accurate, frequently demonstrate lower consistency than methods that use MCMC
sampling and report posterior means or medians. Conversely, most methods that
report posterior mean estimates show seemingly poor accuracy.

Much of the cause of the reduced accuracy of the MCMC based models is driven by
low-expression transcripts. A transcript with no reads has a maximum likelihood
expression of zero, but likelihood is positive over a range of non-zero
values, so the posterior mean must necessarily also be positive. With
a flat, or "non-informative" prior expression values for these low-expression
transcripts are often quite high and highly determined by transcript length.
Isolator avoids this problem by using a model with informative priors: namely
encoding assumptions that most transcripts have low expression and that
expression and splicing tend to be similar between samples. In this way, it is
able to avoid what seems like a tradeoff and produce estimates that are both
consistent and accurate.

Fully Bayesian methods like the one presented here possess a degree of
subjectivity that sometimes give researchers pause. Although we have found the
results from Isolator to be insensitive to precise values given to
hyperparameters, they are chosen in advance and with a degree of arbitrariness.
What is often ignored is that the alternative, methods based on unregularized
maximum likelihood estimates, in practice often necessitate a more insidious
form of subjectivity: contriving cutoffs to filter out unreliable point
estimates of low expression genes or transcripts. While the priors in Bayesian
methods are purely explicit, the means by which data was filtered during an
analysis is often a form of off-the-books subjectivity--critical to the results,
yet unmentioned or only alluded to in manuscripts.

# Online Methods


## Hierarchical Modeling of RNA-Seq experiments

To model RNA-Seq experiments we use a hierarchical model consisting of three
levels: samples or replicates, conditions, and the experiment. Any number of
conditions, and any number of replicates within those conditions. We model
condition-wise and experiment-wise expression for each transcript using a Gamma
distribution, parameterized by mean and shape. Shape, which capture biological
variability, are shared across conditions.

We further parameterize the model to capture splicing. Despite the compasitional
nature, for effeciency purposes we use independent Normal distributions the
model splicing rates at the condition and experiment level, where a splicing
rate is an isoforms proportion relative to total expression of the its gene.

## Efficient Sampling

The isoform quantification or deconvolution problem is most easily thought of as
an additive mixture model in which the component distributions (the isoforms)
are known, but the mixture coefficients (isoform abundances) must be inferred.
The standard approach to this problem is ether expectation maximization (EM) if
the maximum likelihood solution is desired, or Gibbs sampling if, as we are,
estimating the posterior distribution. We have taken a different approach and
instead rely on slice sampling.

The trade-offs between slice sampling and Gibbs sampling for the isoform
quantification problem are subtle. Subsequent samples drawn from a slice sampler
are generally less autocorrelated than those drawn from a Gibbs sampler. As a
result, fewer samples need to be generated to adequately explore a distribution.
However, each of these samples is more computationally expensive to compute, so
one might compensate for autocorrelation in a Gibbs sampler by drawing more
samples in less time. It is not obvious then that one approach is inherently
more efficient than the other.

Because the likelihood function in this problem is simple, and we
sample from the posterior probability directly, without introducing latent
variables, we are afforded some optimization opportunities not available in a
Gibbs sampler. Specifically, we implement the likelihood function using SIMD
(single instruction, multiple data) instructions, which are available nearly all
CPUs from the last decade. This, in combination with some numerical
approximations, allow us to compute the likelihood function over an order of
magnitude faster on conventional hardware, than a naive C implementation.

TODO: point to timing table in supplement

## Correcting for Bias

Isolator attempts to model and correct for multiple factors that can conflate
the transcript expression estimates.

Perhaps the most prominant source of nonuniformity is sequence specific bias
present at fragment ends. A probable source of this bias is cDNA synthesis by
random priming, which is part of popular RNA-Seq protocols such as Illumina's
TruSeq kit. We have previously published a model that sucessfully accounts for
effect by training a sparse graphical model (TODO: cite).

Beyond this random priming bias, there is risidual fragment GC-bias, which may
be an artifact of PCR amplification. We observe that fragments with very high or
very low GC-content are less likely to be observed than expected, even after
accounting for priming bias. We model this by binning fragments by GC content,
computing expected and observed frequencies, and adjusting by their ratio.

Selection of polyadenylated transcripts is a common step mRNA-Seq used to avoid
sequencing introns, partially degraded transcripts, and ribosomal RNA. This can
cause enrichment of reads at the 3' end of transcripts if only partial
transcripts are captured. To model the effect we fit a one parameter model in
which a transcript is truncated at any position with probability $p$. The
probability of observing a fragment at a position $k$ of length $\ell$, counting
from the 3' end, is then proportional to the probability that the transcript was
*not* truncated before $k + \ell$, which is $(1-p)^{k + \ell}$. When fit, $p$ is
typically quite small, on the order of $10^{-5}$, in which case this correction
has little effect on transcripts shorter that several kilobases.

Lastly, subtle implications arise from random fragmentation steps in many
protocals. Existing statisticaly models (see [@Pachter:2011wm] for a review),
assume fragments are sampled uniformly at random from a transcript. However,
this does not exactly match the implications of random fragmentation in which,
under an ideal model, *breakpoints* rather than fragments are introduced
uniformly at random. For a fragment to be observed it must pass size selection
and have fallen between two breakpoints, or one breakpoint and the end of the
transcript.  Since the ends of a transcript act as fixed breakpoints, the result
is some enrichment of fragments at either end of the transcript. We compute this
enrichment exactly for various transcript length, and use interpolation to
approximate th effect.


