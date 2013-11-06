ADSP - Algorithms and Data Structures for Protein Sequence Analysis

A series of tools developed since 1986 for the analysis of protein and DNA
sequences.

Not everything is here yet, more will be added as I add tests and migrate
to Linux. I and starting with a suite of database scanning tools and extending
from there.

Some people might be interested in exploring the linked list implementation
of the sequence model, which along with the raw memory allocation provides a
fast but hair-raising ride for the programmer - not how I would do this today!

These programs are mainly of historical interest for me. They were originally
written to run on VAX/VMS, then Solaris 2.2 and later releases. Multithreading
using Solaris threads was added to take advantage of multi cpu servers in the mid 1990s.

Now, the programs are being updated to run under Linux, using the pthreads library, and
OSX (in particular SOMAP, the screen editor for multiple alignments).

Methods include multiple motif analysis using a position based weight matrix approach,
ROC analysis, multiple sequence alignment (interactive, colour terminal based) and support
for the Prints database of protein fingerprints.

Scanning of databases is multithreaded with memory mapped files for optimisation
of efficiency.

David J Parry-Smith
2013
