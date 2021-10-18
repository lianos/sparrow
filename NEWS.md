# sparrow 1.0 (2021-10-XX)

## Enhancements

* Released to Bioconductor
* Adds support for use of BiocSet as a means by which users can bring their
  genesets to -- or take them from -- sparrow.
* Improvements to the corplot() functionality contributed by by Arkadiusz Gladki
  (@gladki). Users can specify the size of the text reported in the bottom half
  of the pair plot, and spurious/annoying warnings that were produced after a
  a totally valid call are no longer produced.

## Breaking Changes from Pre-release

* First two parameters in ora() function have been swapped so that the
  first parameter (`x`) is the object (data.frame) to run an over
  representation analysis against, and the second parameter is the GeneSetDb.
* scoreSingleSamples no longer drops features in `y` that are not found
  in the GeneSetDb used for scoring. This was changed so that gsva and ssGSEA
  scores match the scores produced by a normal GSVA::gsva call. You can set
  the `drop.unconformed = TRUE` to retain the older behavior.
