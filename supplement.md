
# Measuring agreement between expression estimates

Most of the evaluations of the methods in this paper consist of comparing two
estimates of transcript expression, wether it be multiple runs from an RNA-Seq
experiment, qPCR, or ground truth in a simulation. While this task seems
trivial, many of the obvious approaches one might take, and indeed often *are*
taken, have subtle faults that can lead to unreliable or unintuitive results.
The approach used in this paper has not been previously applied to our
knowledge, so here we describe its motivation.

There is not a "correct" method to measure similarity or disimilarity. Any one
of a variety of measures may be suitable depending on one's preception of the
important features of the data. Our chief concern has been choosing a method
that is both consistent and interpretable.



## 

We adopt the use of mean relative difference (MRD) in this paper. It was
inspired by the use by [@Bray:2015uj] of median relative difference, but
corrects for several several shortcomings of that metric.

$$mrd(u, v) = \frac{1}{n} \sum_{i=1}^{n} \left| \frac{2(u_i - v_i)}{u_i + v_i} \right| $$

## Normalization

MRD has the desirable property that it, along with each constituent term of the
sum, varies between 0 and 1. It would seem then that no small subset of the
transcripts can dominate the MRD, but unfortunately this is not the case in
general. 

Because RNA-Seq is only capable of measuring relative, or compositional
expression, estimates are given in units such as fragments per exonic kilobase
per million fragments (FPKM), reads per million (RPM), or transcripts per
million (TPM). In TPM, for example, estimates must always sum to $10^6$. For
this property to be preserved, if the relative expression of one transcript goes
up, the relative expression of every other transcript must go down to
compensate. In this way, small fluctuations in the estimates of extremely
abundant transcripts can have a sizeable effect, shifting the estimates of all
other transcripts. The end result is that MRD, despite its intent, suffers from
similar problems to Pearson correlation, RMSE, and Euclidean distance. The
measure can by heavily affected by a small number of transcripts.

This notion has been well explored in the context of detecting differential
expression in RNA-Seq by [@Bullard:2010go], among others. The concept is no
different here and the same principles apply: reliable comparison necessitates
normalization beyond FPKM or TPM. To account for this, we adopt the simple
upper-quartile normalization procedure recommended by [@Bullard:2010go]. Before
computing the MRD, we scale the vectors $u$ and $v$ so that their respective
third quartiles are equal.

## Handling Zeros



# Bias Correction

Isolator attempts to model and correct for multiple factors that can conflate
the transcript expression estimates.


## 3' bias

Selection of polyadenylated transcripts is a common step mRNA-Seq used to avoid
sequencing introns, partially degraded transcripts, and ribosomal RNA.

To model the effect we fit a one parameter model in which a transcript is
truncated at any position with probability $p$. The probability of observing a
fragment at a position $k$ of length $\ell$, counting from the 3' end, is then
proportional to the probability that the transcript was *not* truncated before
$k + \ell$, which is $(1-p)^{k + \ell}$. When fit, $p$ is typically quite small,
on the order of $10^{-5}$, in which case this correction has little effect on
transcripts shorter that several kilobases.

## Fragmentation effects

Popular RNA-Seq protocols such as Illumina's TruSeq typically involve a random
fragmentation step. Subtle implications arise from this process.  Existing
statisticaly models (see [@Pachter:2011wm] for a review), assume fragments are
sampled uniformly at random from a transcript. This does not exactly match the
implications of random fragmentation in which, under an ideal model,
*breakpoints* rather than fragments are introduced uniformly at random.  For a
fragment to be observed it must pass size selection and have fallen between two
breakpoints, or one breakpoint and the end of the transcript. Since the ends of
a transcript act as fixed breakpoints, the result is some enrichment of
fragments at either end of the transcript.

## Fragment GC content



