#include <stdio.h>
#include "rr.h"
#include "yheader.h"
#include <string.h>

int TRANSLATION = FALSE;

int pairs[23][23];
char template[] = "ABCDEFGHIKLMNPQRSTVWXYZ";

main(int argc, char **argv)
{
rootres *read_file( int, int);
rootres *first;
int i;
int db_size;

if (argc < 2)
{
	fprintf(stderr,"usage: pair fasta_database\n");
	exit(1);
}
if ((data = fopen(argv[1],"r")) == NULL)
{
	perror(argv[1]);
	exit(1);
}
/* Calculate throughout the database the number of occurrences of each
pair of residues */

/* initialise the pairs matrix */
for (i = 0; i < 23; i++)
{
	int j;
	for (j = 0; j < 23; j++)
	{
		pairs[i][j] = 0;
	}
}

db_size = 0;

while ((first = read_file(1, SCAN_FA)) != NULL)
{
	char *sp1;
	char *sp2;

	db_size += first->numres;
	sp1 = &first->seq[1];
	sp2 = sp1 + 1;
	for (i = 1; i < first->numres; i++)
	{
		pairs[index(sp1)][index(sp2)]++;
		sp1++, sp2++;
	}
	purge_all(first);
}/* End of while */

printf("Residue pairs program.\n");
printf("Database: %s contains %d residues.\n", argv[1], db_size);
printf("The format for the matrix is %s\n\n", template);

printf("Raw data:\n");

for (i = 0; i < 23; i++)
{
	int j;
	for (j = 0; j < 23; j++)
	{
		printf("%d ", pairs[i][j]);
	}
	printf("\n");
}

} /* End of main */


int index(char *sp)
{
char *ptr;

ptr = strchr(template, (int) *sp);
if (ptr == (char *)NULL)
{
	fprintf(stderr, "--- null pointer reference in index.\n");
	abort();
}
return ptr - template;

}/* End of index */
