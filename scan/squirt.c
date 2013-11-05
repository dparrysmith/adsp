#include <stdio.h>
#include "rr.h"
#include "yheader.h"

#define MIN 30

int TRANSLATION = FALSE;

main()
{
rootres *read_file( int, int);
rootres *first;
int pmax, pmin, ns;
int total, nmin;
char tcode[50];



data = stdin;
pmax = 0;
pmin = 1000;
total = 0;
ns = 0;
nmin = 0;
first = read_file(1, SCAN_FA);

do
{
	total += first->numres;
	if (first->numres > pmax)
	{
		pmax = first->numres;
		strcpy(tcode, first->code);
	}
	else if (first->numres < pmin) pmin = first->numres;
	if (first->numres < MIN) 
	{
		printf("Code: %s Residues: %d !\n", first->code, first->numres);
		nmin++;	
	}
	ns++;
	clear_gaps(first);
	write_fasta(first);
	purge_seq(first);
	first = read_file(1, SCAN_FA);
} while (first != NULL);

printf("Number of sequences: %d\n", ns);
printf("Max sequence length: %d\n", pmax);
printf("Min sequence length: %d\n", pmin);
printf("Ave sequence length: %d\n", total/ns);
printf("Sequences less than %d residues: %d\n", MIN, nmin);
printf("Longest sequence was: %s\n", tcode);
}
