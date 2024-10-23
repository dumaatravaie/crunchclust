
Crunchclust is an efficient clustering algorithm that is capable of handling the most common Roche's 454 sequencing error ( Homopolymers ). It uses Levenshtein distance for sequence comparison during clustering. It is also used successfully for the clustering of Illumina Miseq sequences. The software had been developed with the supervision of late Prof. Richard Christen.

Paper: http://www.nature.com/ismej/journal/vaop/ncurrent/full/ismej201284a.html


Introduction

Crunchclust is an efficient clustering algorithm that is capable of handling the most common Roche's 454 sequencing error ( Homopolymers ). It uses Levenshtein distance for sequence comparison during clustering. It can also handle Illumina Miseq sequences.

Details

Crunchclust is an efficient clustering ( Search version is under construction and will be released soon ), alignment algorithm that is capable of handling the most common Roche's 454 sequencing error ( Errors due to the presence of Homopolymers in the sequence ). It is also effective in clustering illumina miseq sequences. It uses Levenshtein distance for sequence comparison during clustering. It has the additional capability of performing strict dereplication on the raw sequence dataset and sorting of the raw sequences according to their abundance before doing actual clustering. In this way it reduces the dimension of the datasets, increases the clustering accuracy and reduces the computational complexity. It also provides the user the option of running it for different distance thresholds between two intervals ( kmin and kmax ) with a single command. Parallel version of Crunchclust is under construction to meet the special needs of scientific community.

Levenshtein distance is a metric for measuring the amount of differences between two sequences (ie an edit distance). As the 454 sequences tend to start at the same place, Levenshtein distance is the most appropriate measure for calculating distances between them.
Crunchclust has given the user the freedom of choice whether or not to count the end gaps between the sequences as distance. This is useful as the 454 sequences may or may not end at the same place.
