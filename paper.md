
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
    which uses a simple Bayesian model combined with agressive bias correction
    to produce estimates that are simultaneously accurate and stable.  In a
    comprehensive comparison of accuracy and stability, we show that this
    property is unique to Isolator.

    This model is then extended to encompass an entire RNA-Seq experiment.
    Analyzing all samples in unison allows us to share information, yielding
    more accurate estimates and enabling more sophisticated analysis beyond
    pairwise hypothesis tests. We demonstrate this by analyzing splicing
    monotonicity across several time points in the devolopment of human
    cardiomyocyte cells.
...


<!--
    Keep in mind
    Abstract: ~150 words.
    Paper: 2500-3500 words.
           <= 6 figures or table

    Structure
      * Introduction
      * Results
      * Discussion
      * Online Methods
        Use subheadings here, go into details, etc.
-->



# Introduction


# Results


# Discussion


# Online Methods

## Hierarchical Modeling of RNA-Seq experiments

TODO

## Efficient Sampling

TODO

## Correcting for Bias

TODO








<!-- -->



# Methods

## Moderated effect-size estimates

Anyone setting out to analyze RNA-Seq data today will find a curious bifurcation
in the literature and tools that have emerged around it. One branch focuses
first on accurately estimating transcript abundance, before extending their
analysis to consider differential expression or splicing. In this group we may
include for example Cufflinks, RSEM, IsoEM, BitSeq, and Sailfish. The second
branch instead focuses first on differential expression using count data, and
have worked to refine these methods to consider isoform or exon level events
(edgeR, DESeq, BaySeq, DEXSeq).

Neither approach is a panacea. Methods that consider count data only are less
able to account for the nastier details in RNA-Seq data: multiple alignments,
positional bias, sequence bias, and fragment length distributions, to name
several. These details are considered by transcript abundance methods, but under
this regime, differential expression is typically not informed by the full
likelihood model, but rather treated as a secondary step to be performed after
abundance estimation.

<!-- I need to back this up with specific critism. -->

The shortcomings of count-based differential expression approaches has been well
articulated by others. As Trapnell et al demonstrates, count-based approaches
can perform poorly when it comes to overlapping isoforms in which reads can not
be easily assigned. Even if analysis is done on the gene level, separate
regulation of constituent isoforms of the gene can lead to incorrect
conclusions.

Approaches 


Our method was developed from the perspective that the compartmentalization of
transcript abundance and differential expression is artificial and harms
analysis. Differential expression analysis should be fully informed by the
uncertainty in estimates of transcript abundance, and abundance estimation can
and should be informed by the structure of th experiment.



~~~~~~~~~~~~~~~~~~~~~~~


Studies of gene or transcript expression typcially seek to understand the
dynamics of expression: what changes in response to a perturbation. 

<!--What's a catchy way I could describe the ignorance to effect size.-->

Moderated or regularized estimated are not a new concept in gene expression
analysis. Recently 



A recent paper by Love et al. (TODO: cite) acknowledges this deficiency and takes
steps to address it, presenting a sequel to the popular differential expression
package DESeq (TODO: cite). They present an empirical Bayes method in which
estimates of fold-change are shrunk towards 


Emperical Bayes approaches, like the one developed by Love et al, are an
epproximation to a full hierarchical Bayesian 


~~~~~~~~~~~~~~~~~~~~~~~~~


## A hierarchical model for RNA-Seq experiments

## Correcting for Bias

<!-- Sequence bias, GC bias, 3' bias, fragmentation bias. -->

The standard model adopted for most RNA-Seq analysis (with a few notable
exceptions such as TODO) assumes short fragments are sampled uniformly at random
from RNA transcripts (TODO: cite pachter's review paper). This model was
expanded to account for sources of positional bias by TODO cite Roberts paper.
Yet modeling specific sources of bias is often

Despite the severity of the bias, Hansen and others than followed (TODO: bias
correction papers) understated the bias. Because common library preparation
protocols like Illumina's TruSeq kit perform both first- and second-strand cDNA
synthesis using random priming, fragments are biased at both ends. At the
fragment level, kkkkkkkkkkkkkkkkkkkkkkkk


Modeling and correcting for bias is central to doing accurace quantification at
the transcript level.

TODO: Point out that positional bias can have a huge effect on relative
transcript abundance.







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

This approach has significant room to scale. SIMD units have been continually
grown in subsequence CPU generations, while current GPUs are able to operate on
far more data at once than CPUs. 

# Results

The most obvious and frequently employed means of evaluating a method for
RNA-Seq quantification is to compare point estimates made by the program to
either ground truth, if the data is simulated, or if the data is real, to some
proxy for ground truth, such as qPCR.

Studies of gene and transcript expression typically seek to understand the
dynamics of expression.

TODO: I don't know what I intended to get at here.



Comparisons of gene or transcript expression is sensitive to how one chooses to
measure simularity or distance. Because gene expression typically varies across
several orders of magnitude, measurements like the Pearson correlation, $L_1$ or
$L_2$ distance, or root-mean-square error are nearly entirely dominated by a
small subset of the most highly expressed genes. Recognizing this, a popular
alternative is the Spearman rank correlation. It is generally more robust than
Pearson correlation, and more considerate of low expression genes, but comes
with its own problems. Since Spearman correlation considers only rank, there is
no penalty whatsoever if a method produces estimates that are highly distored,
so long as that distortion is monotonic. Furthermore, genes can swap ranks with
arbitrarilly small perturbations to expression estimates, which can make
estimates appear inconsistent that are in fact in very tight agreement.

With these problems in mind, we have adopted the approach of measuring Pearson
correlation on log-transformed values. To account for the fact that most methods
report estimates of zero for transcripts with no compatible reads, we add a
small constant $c$ to compute $log(x_i + c)$ for each expression estimate $x_i$.
In the microarray literature this is commonly described as a "started-log"
transform. (TODO: citation) The value of $c$ is chosen as the smallest nonzero
$x_i$ (i.e. $c = \min_{i, x_i \neq 0} x_i$).

The justification for this is to
maintain a significant separation between zero and non-zero estimates, 



Is there *really* a clear justification for this. Really what I want is to
choose $c_1, c_2$, such that min(x_i) is colinear with max(y_i).





For consistency we applied the started-log transformation to all estimates
regardless of the presence or absence of zeros.

Though we consider this a more reasonable, consistent, and interpretable method
for comparing RNA-Seq gene or transcript expression estimates, Spearman
correlation is more conventional and so we repeated these analysis using
Spearman correlation and included these results in the supplement (TODO:
supplementary section). The two measures agree qualitatively in the majority
of cases, but the magnitude of the difference between methods can vary.

## Agreement with qPCR

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


## Agreement with Spike-in Controls

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


## Estimation Consistency

Method      $C$ vs $0.75A + 0.25B$   $D$ vs $0.25A + 0.75B$
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

Method      $C$ vs $0.75A + 0.25B$   $D$ vs $0.25A + 0.75B$
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
corresponding to the mixture proportions for C and D.


## Batch effects

![](analysis/seqc/batch-effect-comparison.pdf)


![](analysis/seqc/batch-effect-bias-correction.pdf)

The SEQC data gives a unique oppourtunity to examine how quantification methods
cope with batch effects. The analysis performed by the SEQC Consortium evaluated
three pipelines for transcript-level quantification: BitSeq using alignments
from SHRiMP2 (TODO: cite), and Cufflinks using alignments from TopHat2 run with
and without existing gene models. Here we present a broader anlysis on a subset
of the SEQC data, in which choice of alignment tool is consistent and thus
controlled for, and a wide variety of popular tools are compared.

We used 10 samples from the larger SEQC dataset, each consituting a single
flowcell lane from two separate flowcells sequenced at five different sites, all
on Illumina HiSeq 2000 instruments.

The SEQC Consortium showed much higher consistency in estimates from BitSeq as
compared to Cufflinks, but was not in a position to explain why. Here, In the
context of other methods it becomes obvious that the defining separation is by
choice of estimator: maximum likelihood (Cufflinks, eXpress, RSEM/ML, Sailfish,
Salmon, Kallipso) or posterior mean (Isolator, BitSeq, RSEM/PM). Beyond this
separation, the methods differ greatly in their approach to bias correction.

To 




## Accuracy of isoform-level estimates in simulations

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
by each method.


## Aggressive bias correction significantly improves accuracy



## Posterior mean estimates are more consistent that maximum likelihood estimates



# Discussion




The convention for evaluating the accuracy of RNA-Seq analysis methods has until
now consisted chiefly of comparing point estimate accuracy in simulations, and
occasionally accompanying these results with comparisons to qPCR or other
technoligies. This severely limited evaluation regime has incentivized methods
tuned to solve a theororical model of RNA-Seq that bears only a resemblence to
real data. Though RNA-Seq data has been known to suffer from rather severe
sources of bias (TODO: cite), newer methods 




<!-- Probably cut this stuff. -->
Data from previous studies (TODO: cite), as well as the data presented here,
have consistently demonstrated that bias correction elicits a moderate
improvement in agreement with qPCR. 

additionally presented data from a more nuanced simulation of RNA-Seq that
includes commonly understood sources of bias.

A more
nuanced simulation of RNA-Seq that includes commonly understood

It should not be surprising then that this improvement is greatly amplified when
working at the isoform level, as the simulation results suggest. At the gene
level, position-level biases may average out, dampening their overall effect.
Yet the difference between two isoforms is often slight, sometimes just several
nucleotides. Understanding the technical effects that may cause non-uniform
sampling surrounding these nucleotides is thus vital to accurately inferring
the degree to which jjjjjjjjjjj.




Furthermore, in focusing soley on accurate point estimates in individual
isolated samples, the role of RNA-Seq as a means of comparing multiple
biological conditions has been occluded. Our data builds on the important work
done by the SEQC Consortium in understanding consistency and reproducability. We
show that MCMC and posterior mean estimates can provide a more practical
estimate of expression without sacrificing accuracy, simply by being more
predictable is the presense of limited data.



Our method, based on read alignment and MCMC sampling, though generally faster
than other methods than rely of full read alignment, is not currently
competitive in performance with the newer alignment-free methods Sailfish,
Salmon, and Kallisto.
<!--TODO: 1. sampling is inherently much more expensive maximum likelihood, but
it's probably worth it. 2. we've demonstrated a basic methodology in which a
very efficient sampler can be implemented for this problem (slice samplers with
SIMD).
-->


Among RNA-Seq gene expression studies, it's aknowledged, though sometimes
implicitly, that effect size matters. The SEQC Consortium recommends filtering
differential expression calls based on fold change and overall expression in
order to achieve high consistency.
(TODO: cite)


I want to make the case that filtering RNA-Seq is common:
Examples:

 * A comparison of methods for differential expression analysis of RNA-seq data
    Filter out genes for which the total count across samples was less than 10.

 * 



It is surpsising then to see so little attention payed to producing reliable
estimates of effect size. The recent work by Love et al on DESeq2, which
developes a regularized logarithm transformation, effectively moderating effect
size estimates, is an important and welcomed step towards taking effect size
seriously. Yet DESeq2 is not able to do full isoform deconvolution and
transcript expression estimates, as Isolator and the other tools evaluated here
do.


Fully Bayesian methods like the one presented here possess a degree of
subjectivity that sometimes give researchers pause. Although we have found the
results from Isolator to be insensitive to precise values given to these
hyperparameters, they are chosen in advance and with a degree of arbitrariness.
What is often ignored is that the alternative, methods based on unregularized
maximum likelihood estimates, in practice often necessitate a more insidious
form of subjectivity: contriving cutoffs to filter out unreliable point
estimates of low expression genes or transcripts. While the priors in Bayesian
methods are purely explicit, the means by which data was filtered during an
analysis is often a form of off-the-books subjectivity--critical to the results,
yet unmentioned or only alluded to in manuscripts.







### Future Work

Not all detectable or statistically significant changes in gene expression or
splicing are of scientific interest, as the effect size is often minute. If the
goal is to uncover interesting changes in RNA expression, explicitly
considering effect size moves us far closer asking the right question.

Though we consider the approach described here an major improvement, it
simplistically defines a minimally interesting effect size cutoff uniformly
across genes and transcripts. This ``one effect size fits all'' model is clearly
wrong. Wh


<!--

    Effect sizes should not be uniform, but we lack an appropriate method to

-->






