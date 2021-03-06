\name{allezPlot}
\alias{allezPlot}

\title{
Display an image of gene scores by functional sets
}
\description{
Plots functional sets against gene scores.  Sets are ordered by highest
sum of gene scores not in preceding sets, and genes are ordered by scores
within set. Higher gene scores are indicated by darker color.
}
\usage{
allezPlot(allez.out, n.low = 5, n.upp = 500, n.cell=2, nominal.alpha=0.01, gmax=20, maxSets=10,
glab=c("none","gene_id","symbol"), slab=c("none","z.score","set.means"),
...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allez.out}{list output from \code{allez} function, or
    concatenated output from \code{allezC}.}
  \item{n.low}{numeric, smallest number of genes in a set to consider
    for ranking}
  \item{n.upp}{numeric, largest number of genes in a set to consider for
    ranking}
    \item{n.cell}{numeric, minimum number of genes in a set with score>0,
  mainly used in conjunction with binary scores}
 \item{nominal.alpha}{numeric, keep sets with z-scores exceeding Bonferroni level at this value}
  \item{gmax}{numeric, maximum number of genes for x-axis labelling}
  \item{maxSets}{numeric, maximum number of sets in plot}
  \item{glab}{character string, adds gene labels to the x-axis: "none"
    indicates no gene labels, "gene_id" indicates entrez id, and
    "symbol" is gene symbol. Useful mainly for binary scores.}
  \item{slab}{character string, adds set statistics to right axis:
    "none" indicates no labels, "z.score" indicates z.score, and
    "set.means" is set means (from \code{setscores} element)}
  \item{\dots}{any extra arguments for \code{text} or \code{axis}, such
    as \code{cex} or \code{cex.axis}}
}
\details{
\code{allezPlot} creates an indicator matrix of the highest ranking sets
using \code{allezMat}, adds gene score information, and orders the genes
by highest ranking set and gene score.

If x-axis labelling is not chosen (\code{glab} = "none"), the plot width
must be longer than the width of the set labels (the text appears in the
plot region); the user may need to adjust the width of the figure in
order to encompass the text.

If x-axis labelling is chosen (\code{glab} = "gene_id" or "symbol"), the
plot width will be predefined, based on character size; the user may
need to adjust the width of the figure in order to encompass the text.
}
\author{
Aimee Teo Broman and Michael Newton
}
\seealso{
  \code{allezMat}
  \code{allez}
}
\examples{
data(npc)
scores  <- (1/2)*sqrt(28)*log((1-npc)/(1+npc))

npc.kegg <- allez( scores=scores, lib="hgu133plus2", sets="KEGG")

allezPlot(npc.kegg)

\dontrun{
npc.go <- allez( scores=scores, lib="hgu133plus2", sets="GO")
npc.all <- allezC(npc.go, npc.kegg)

allezPlot(npc.all)
}
}
