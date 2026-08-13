# Syncmer Deduplication

## Installation

###  Requirements
- Requires CMake 3.25+
- Requires zlib 1.2.3+
- Will fetch [zstr](https://github.com/mateidavid/zstrr)

### Installation
```
cmake --workflow release
cmake --install build
```

## Description

This program is intended to provide a way of deduplicating reads which combines the best of string deduplication and mapping based deduplication while also handling substitution errors during sequencing.

Mapping based deduplication requires ... mapping, but also can by stymied by differences at fragment termini. For example if one performs quality trimming of the reads, then the termini between duplicates may no be different. If instead one skips the quality trimming, a sequencing error leading to soft clipped bases on one duplicates termini will also escape deduplication (this can also happen with quality trimming). Additionally mapping based deduplication ignores the actual sequence. This is much less often a problem, but one could conceive of a situation where both termini of a molecule are the same and map to the same place, but the interior ends of the reads are different. In one case they might be soft clipped, and not in the other and they would be marked as duplicates regardless.

Enter a read sketch based method that assesses if two reads have the same sequence.

The concept is to take a sketch of each read where a sketch is the set of bounded syncmers and their positions in the original read. Overlapping reads will share some syncmers, but those syncmers will not be at the same position within the read. After indexing sketch elements across all reads, candidate matches for a read can be quickly identified as those which share an above threshold number of sketch elements. The optimal values for thresholds and sketching parameters should be estimable from the data. 

This is the core concept, however UMI's solve this problem completely. So it probably isn't worth pursuing.

# Dev Notes

Testing requires that Catch2 v3+ is installed and findable by CMake
