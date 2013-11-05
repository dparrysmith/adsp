#include <stdio.h>
#include "rr.h"
#include "yheader.h"

int TRANSLATION = TRUE;

main(int argc, char **argv)
{
void init_gencode();
void test_gencode();
rootres *translate(rootres *, int);
rootres *read_file( int, int);
int write_fasta(rootres *);
rootres *first;

if (argc < 2)
{
	fprintf(stderr,"usage: translate file\n");
	exit(1);
}
init_gencode();
if ((data = fopen(argv[1],"r")) == NULL)
{
	perror(argv[1]);
	exit(1);
}
while (feof(data) == 0)
{
	first = read_file(1, SCAN_FA);
	translate(first, SIX_FRAME);
	write_fasta(first->next);
	purge_all(first);
}
exit(0);
}
