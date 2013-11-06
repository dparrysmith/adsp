/*

ROOTRES database and sequence handling data structures.

Author:		D J Parry-Smith
Date :		July 1989
Updated:	February 1994 (translations)
Version:	2.0	
Implemenation:	Unix System V and VAX/VMS

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

/* #include <malloc.h>*/ /* for SUNOS */
#include "rr.h"
#include "yheader.h"

#define SLOW 1

#define MAX_LINE BUFSIZ

static char gb_template[]={"ARNDCQEGHILKMFPSTWYVBZX-"};
static char gcode_template[]={"TCAGN"};
static char gcode_template6[]={"AGTCN"};
int nseqs;
static int modified = 1;
static char *place = (char *) NULL;
char file_nm[]="rr.seqs";
extern char *src;
extern FILE *data;
extern int TRANSLATION;
FILE *gcode;

/* for the translator */
char *gencode[5][5];

rootres *read_file(int num, int mask)
{
/* num = number of sequences to read at once */
/* mask = value to indicate which sort of sequences to read */
/* current values are 0 (all) or 1 (only >P1 sequences) */
/* Read sequences from the sequence file */
rootres *create_root_residue(rootres *);
int blank();
int validate(char);
rootres *root;
rootres *root_one;
char line[MAX_LINE];
char tempcode[MAX_LINE];
char res;
int cgap;
register int cres;
register int numres;
int start_seq;
int counter;
char modifier[20];
char *pos;
char *p;

if (mask == SCAN_ALL)
{
	strcpy(modifier,">%*[^;];%s");
}
else if (mask == SCAN_P1)
{
	strcpy(modifier,">P1;%s");
}
/* 6/4/93 add mask for Fasta format files */
else if (mask == SCAN_FA)
{
	strcpy(modifier,">%s");
}
else
{
	fprintf(stderr,"** modifier mask not recognised <%d> !\n", mask);
	exit(0);
}
counter = 0; /* Make sure don't read in more than num sequences */
while (fgets(line, MAX_LINE, data) != NULL)
{
	if ((pos = strchr(line,'\n')) != 0) *pos = '\0';
	if (blank(line) == 1) continue;
	if (sscanf(line,modifier,tempcode) != 1) continue;
	/* New sequence found so establish a new structure for it */
	root = create_root_residue((rootres *)NULL);
	if (root == NULL)
	{
		fprintf(stderr,"Problem allocating root!\n");
		perror("root");
		exit(1);
	}
	if (counter == 0) root_one = root;
	root->code = (char *)malloc(strlen(tempcode)+1); /* For '\0' */
	if ((p = strpbrk(tempcode, ",;.!?£$&^*")) != NULL) *p = '\0';
	if (root->code != NULL) strcpy(root->code, tempcode);
	else 
	{
		printf("Malloc failed in read_file tempcode !");
		exit(0);
	}
	if (mask != SCAN_FA) /* 6/4/93 */
	{
	/* Pir format has a separate title line */
	if (fgets(line, MAX_LINE, data) == NULL)
	{
		printf("Premature end of file - expecting comment !");
		exit(0);
	}
	}
	else
	{
		/* SCAN_FA */
		/* Ignore code part of title */
		/* Find the first space */
		if (( pos = strchr(line,' ')) == NULL)
		{
			/* no title for this code */
			strcpy(line, "<no title found>");
		}
		else
		{
			/* Make pos the beginning of the line */
			strcpy(line, pos); /* This may cause problems */
		}
	}
	line[strlen(line)] = '\0';
	strcpy(tempcode,line);
	root->comment = (char *)malloc(strlen(tempcode)+1);
	if (root->comment != NULL)
	{
		strncpy(root->comment, tempcode, strlen(tempcode));
		root->comment[strlen(tempcode)] = '\0';
	}
	else 
	{
		printf("Malloc failed in read_file comment !\n");
		exit(0);
	}
	/* Now find how many residues there are in this sequence */
	start_seq = ftell(data);
	cgap = 0;
	cres = 0;
	if (mask != SCAN_FA)
	{
	while((res = getc(data)) != '*')
	{
		if (validate(res) == 0)
		{
			if (res == '-') cgap++;
			else cres++;
		}
	}
	}
	else
	{
		/* Fasta format */
		while((res = getc(data)) != '>')
		{
			if (res == EOF) break;
			else if (validate(res) == 0)
			{
				if (res == '-') cgap++;
				else cres++;
			}
/*
			if (feof(data) != 0) break;
*/
		}
		if (res == '>') ungetc(res, data); /* Push > back */
	}
	/*printf("creating space for %d residues (%d gaps)..\n",cres,cgap);*/
	root->seq = (char *)calloc(cres+3,sizeof(char));/* +2 = 0th residue + \n */
	if (root->seq == NULL)
	{
		printf("*** calloc failed root->seq\n");
		exit(0);
	}
	root->gap = (int *)calloc(cres+3,sizeof(int));/* same length as seq */
	if (root->gap == NULL)
	{
		printf("*** calloc failed root->gap\n");
		exit(0);
	}
	root->fdata = (float *)calloc(cres+3,sizeof(float));
	if (root->fdata == NULL)
	{
		printf("*** calloc failed root->fdata\n");
		exit(0);
	}
	if ((fseek(data,start_seq,0)) == EOF)
	{
		printf("*** improper use of fseek() !\n");
		exit(0);
	}
	root->assoc = (pres *)NULL;
	numres = 0;
	root->numres = 0;
	root->numgap = 0;	
	/* The 0th residue is always special character */
	root->seq[0] = ':';
	root->gap[0] = 0;
	numres++;
	if (mask != SCAN_FA)
	{
	while((res = getc(data)) != '*')
	{
		if (validate(res) == 0)
		{
			if (res == '-')
			{
				root->gap[numres]++;
				root->numgap++;
			}
			else
			{
			if (TRANSLATION == TRUE)
			{
			if (strchr(gcode_template, res) == NULL) res = 'N';
			}
			root->seq[numres] = res;
			numres++;
			}
		}
	}
	}
	else
	{
	while((res=getc(data)) != '>')
        {
		if (res == EOF) break;
                else if (validate(res) == 0)
                {
                        if (res == '-')
                        {
                                root->gap[numres]++;
                                root->numgap++;
                        }
                        else
                        {
			if (TRANSLATION == TRUE)
			{
			if (strchr(gcode_template, res) == NULL) res = 'N';
			}
                        root->seq[numres] = res;
                        numres++;
                        }
                }
/*
		if (feof(data) != 0) break;
*/
        }
	if (res == '>') ungetc(res, data);
	}
	root->seq[numres] = '\0';
	numres -= 1; /* make numres reflect the actual number of residues*/
	root->numres = numres;
	counter++;
	if (counter >= num) return(root_one);
}/* End while loop */
if (counter > 0) return(root_one);
else return((rootres *)NULL); /* Have run out of file space */
}/* End of read_file */

#ifndef SLOW

rootres *read_file_quick(int fdin, int num, int mask)
{
char modifier[20];
rootres *create_root_residue(rootres *);
rootres *root, *root_one;

if (mask == SCAN_ALL)
{
	strcpy(modifier,">%*[^;];%s");
}
else if (mask == SCAN_P1)
{
	strcpy(modifier,">P1;%s");
}
/* 6/4/93 add mask for Fasta format files */
else if (mask == SCAN_FA)
{
	strcpy(modifier,">%s");
}
else
{
	fprintf(stderr,"** modifier mask not recognised <%d> !\n", mask);
	exit(0);
}
/* read only num sequences */
counter = 0;
/* set place to beginning of file */
if (place == (char *) NULL)
{
	place = src;
	newp = strchr(place, '>');
	if (newp == (char *) NULL)
	{
		fprintf(stderr,"*** can't find valid FASTA format code\n");
		exit(0);
	}
}
else newp = place; /* should be a > */
root = create_root_residue((rootres *)NULL);
if (root == NULL)
{
	fprintf(stderr,"Problem allocating root!\n");
	perror("root");
	exit(1);
}
if (counter == 0) root_one = root;
if (sscanf(newp, modifier, tempcode) != 1)
{
	fprintf(stderr,"** bad database format %s did not match\n",modifier) 
	exit(1);
}
root->code = (char *)malloc(strlen(tempcode) + 1);
if ((p = strpbrk(tempcode, ".;,!?£$&^*")) != NULL) *p = '\0';
if (root->code != NULL) strcpy(root->code, tempcode);
else
{
	fprintf(stderr,"Malloc failed in read_file_quick tempcode\n");
	exit(0);
}
if (mask != SCAN_FA)
{
	/* PIR format has a separate line for the title */
	if ((p1 = strchr(newp,'\n') != NULL) *p1 = '\0';

	newp++; /* beginning of the next line */
	strcpy(tempcode, newp);
	tempcode[80] = '\0';
}
else
{
	/* SCAN_FA */
	if ((p1 = strchr(newp,'\n') != NULL) *p1 = '\0';
	if ((pos = strchr(newp,' ')) == NULL)
	{
		/* no title for this code */
		strcpy(tempcode,"<No title found>");
	}
	else
	{
		strcpy(tempcode, pos);
	}
	/* terminators ? */
	tempcode[80] = '\0';
	newp = p1;
	newp++;
}
root->comment = (char *)malloc(strlen(tempcode) +1);
if (root->comment != NULL)
{
	strncpy(root->comment, tempcode, strlen(tempcode));
}
else
{
	fprintf(stderr,"Malloc failed in read_file comment !\n");
	exit(0);
}
/* Now find the next > */

}/* End of read_file_quick */
#endif


int blank(char *line)
{
/* Function determines whether the input string is all white space */
int len;
int i;

len = strlen(line);
if (len == 0) return(0);
for (i=0; i<=len; i++)
{
	if (isspace(line[i]) == 0) return(0);
}
return(1);
}/* End of blank */


rootres *create_root_residue(rootres *mres)
{
/* Make space for one more root residue structure */
rootres *current;
rootres *last;

current = (rootres *)malloc(sizeof(rootres));
if (current == NULL)
{
	printf("*** No space returned by malloc in create_root_residue !");
	exit(0);
}
/* Step along the chain of residues until find the end and add new structure */
if (mres == NULL)
{
/* this is the first structure */
current->previous = NULL; /* There is no predecessor of this one */
current->next = NULL; /* there is no successor of this one */
current->assoc = (pres *)NULL;
return(current);
}
/* Add it to the end of the list */
last = mres;
while (last->next != NULL) last = last->next;
last->next = current;
current->previous = last;
current->next = NULL; /* No successor yet */
current->assoc = (pres *)NULL;
return(current);
}/* end of create_root_residue */


void dump(rootres *mres)
{
/* Debug procedure */
/* List the contents of structures allocated */
rootres *current;
int count;
int inc;

current = mres;
count = 1;
while (current != NULL)
{
	printf("Sequence: %d at %p\n",count,current);
	printf("code = %s\n",current->code);
	printf("comment = %s\n",current->comment);
	if (current->assoc != NULL)
	{
		printf("assoc = %p\n", current->assoc);
		printf("assoc->assoc = %p, %s\n", current->assoc->assoc, current->assoc->assoc);
		printf("assoc->next = %p\n", current->assoc->next);
	}
	printf("Purged sequence :\n%s\n",current->seq);
	printf("\nAligned sequence :\n");
	for (inc=0; inc<=current->numres; inc++)
	{
		if (current->gap[inc] != 0)
		{
		int loop;
			for (loop = 0; loop <= current->gap[inc]; loop++)
			{
				putchar('-');
			}
		}
		putchar(current->seq[inc]);
	}
	putchar('\n');
	current = current->next;
	count++;
}
return;
}/* End of dump */


int validate(char res)
{
/* this procedure validates the given residue, if it is not in the valid set of
residues return -1 to indicate that this sequence should be skipped */
char *strchr();
char *la;

la = strchr(gb_template,res);
if (la == NULL)
{
	return(-1);
}
else
{
	return(0);
}
}/* End of check_residue */

int numbers(rootres *mres)
{
/* print out some statistics about residues and so forth */
/* also set the global variable nseqs for use elsewhere */
rootres *current;
char msg[80];
int tot_seqs;
int tot_res;
int tot_gap;

current = mres;
tot_seqs = 0;
tot_gap = 0;
tot_res = 0;
while (current != NULL)
{
	tot_seqs++;
	tot_res += current->numres;
	tot_gap += current->numgap;
	current = current->next;
}
sprintf(msg,"Total number of sequences : %d  (%d residues, %d gaps)\n",
	tot_seqs,tot_res,tot_gap);
puts(msg);
return(tot_seqs);
}/* End of numbers */

int residues(rootres *mres)
{
	return(mres->numres);
}/* End of residues */

int sequences(rootres *mres)
{
int tot_seqs;
rootres *current;

tot_seqs = 0;
current = mres;

while (current != NULL)
{
	tot_seqs++;
	current = current->next;
}
return(tot_seqs);
}/* End of sequences */

void gen_pos(int start_at, rootres *seq_struct)
{
/*
find out where abouts in the seq array to start
*/
int count;

if ((start_at == 0) || (start_at <= seq_struct->gap[0])) 
{
	seq_struct->marker[0] = 0;
	seq_struct->marker[1] = seq_struct->gap[0];
	return;
}
if (start_at > (seq_struct->numres + seq_struct->numgap))
{
/* Out of range */
	seq_struct->marker[0] = -1;
	seq_struct->marker[1] = 0;
	return;
}
count = 0;
seq_struct->marker[0] = 0;
seq_struct->marker[1] = 0;
while (count <= start_at)
{
	count += seq_struct->gap[seq_struct->marker[0]];
	seq_struct->marker[0]++;
	count++;
}
if (count-start_at > 0) 
{
	seq_struct->marker[1]=count-start_at -1;
	seq_struct->marker[0]--;
}
else if (count-start_at == 0 )
{
	seq_struct->marker[1]=0;
	seq_struct->marker[0]--;
}
else if (count-start_at < 0)
{
	seq_struct->marker[0]++;
}
return;
}/* End of gen_pos */


void gen_seq(rootres *seq_struct, char *string, int number, int position)
{
/*
Create a gapped sequence string of length number from the sequence seq_struct.
Read the sequence header to find out where to start.
*/
int i;
int length;
int strp;
int point;
int kk;
int maxl;

if (seq_struct->marker[0] == -1)
{
/* return a blank padded string */
for (i=0; i<number; i++)
{
	string[i] = ' ';
}
string[number] = '\0';
return;
}/* That was for past the end of sequence */
maxl = seq_struct->numres + seq_struct->numgap;
point = seq_struct->marker[0];
strp = 0;
length = number;
kk = (seq_struct->marker[1] > length) ? length : seq_struct->marker[1];
if (kk > 0)
{
for (i=0; i < kk; i++)
{
	string[i] = '-';
	length--;
	strp++; /* Points to the next element of string to be filled */
}
}/* End if kk>0 */
else
{
i = 0;
}
if (i > number)
{
	string[number] = '\0';
	return;
}
string[strp++] = seq_struct->seq[point++];
/* Do it fast if not past sequence end else include check */
if ((maxl - position) < number)
{
/* We will run out of sequence before end of string */
/* How many residues are left ? */
int left;
left = seq_struct->numres - seq_struct->marker[0] +1;
for (i = 0; i<=left; i++)
{
int j;
	for (j=0; j<seq_struct->gap[point]; j++)
	{
		string[strp++] = '-';
	}
string[strp++] = seq_struct->seq[point++];
}/* End for i */
/* Now fill up with blanks */
for (i = strp; i < length; i++)
{
	string[i] = ' ';
}
string[number] = '\0';
return;
}
else
{
/* Will not run out of sequence */
for (i=0; i<length; i++)
{
int j;
int index;
	index = (i + seq_struct->gap[point] > length) ? length - i : seq_struct->gap[point];
	for (j=0; j < index; j++,i++)
	{
		string[strp++] = '-';
	}
string[strp++] = seq_struct->seq[point++];
}/* End for i */
string[number] = '\0';
return;
}
}/* End of gen_seq */


write_als_file()
{
int ret_val;
char msg[60];

modified=1; /* Set modified flag to allow output even if buffer unmodified */
printf("File name :");
scanf("%s",file_nm);
ret_val=write_als();
sprintf(msg,"%d lines written to %s",ret_val,file_nm);
puts(msg);
return(0);
}

write_als(rootres *mres)
{
int j,lc;
double tlen,iter;
double ceil();
char seq[100];
rootres *current;

if (!modified) return(0);
data = fopen(file_nm,"w");
modified =0;
lc =0;
current = mres;
while (current != NULL)
{
	tlen = current->numres + current->numgap;
	tlen /= 80.0;
	iter = ceil(tlen);
	fprintf(data,">P1;%s\n",current->code);
	lc++;
	fprintf(data,"%s\n",current->comment);
	lc++;
	for(j=1;j<(int)iter*80;j += 80)
	{
	gen_pos(j,current);
	gen_seq(current,seq,80,j);
	fprintf(data,"%-.80s\n",seq);
	lc++;
	}
fprintf(data,"*\n");
lc++;
current = current->next;
}
fclose(data);
return(lc);
}


write_fasta(rootres *start_with)
{
int j,lc;
double tlen,iter;
double ceil();
char seq[100];
rootres *current;

lc =0;
current = start_with;
while (current != NULL)
{
	tlen = current->numres + current->numgap;
	tlen /= 80.0;
	iter = ceil(tlen);
	fprintf(stdout,">%s " ,current->code);
	lc++;
	fprintf(stdout,"%s\n",current->comment);
	lc++;
	for(j=1;j<(int)iter*80;j += 80)
	{
		gen_pos(j,current);
		gen_seq(current,seq,80,j);
		fprintf(stdout,"%-.80s\n",seq);
		lc++;
	}
lc++;
current = current->next;
}
return(lc);
}


void clear_gaps(rootres *seqstr)
{
/* zero the gaps array */
int i;

for (i=0; i<=seqstr->numres; i++)
{
	seqstr->gap[i] = 0;
}
seqstr->numgap = 0;
return;
}/* End of clear_gaps */


void purge_seq(rootres *seqstr)
{
/* free memory associated with seqstr */
/* First preserve link pointers */
/* Is there a previous residue ? */
if (seqstr->previous != NULL)
{
	seqstr->previous->next = seqstr->next;
}
if (seqstr->next != NULL)
{
	seqstr->next->previous = seqstr->previous;
}
/* Now deallocate arrays */
if (seqstr->seq != NULL) cfree(seqstr->seq);
if (seqstr->code != NULL) free(seqstr->code);
if (seqstr->comment != NULL) free(seqstr->comment);
if (seqstr->gap != NULL) cfree(seqstr->gap);
if (seqstr->fdata != NULL) cfree(seqstr->fdata);
free(seqstr);
return;
}


void purge_all(rootres *seqstr)
{
/* starting with one seqstr, delete space allocated by all descendents */
/* If a list, move to the end of the list */
rootres *p0, *p1;
pres *p2;
char *code_freed;
char *comment_freed;

p0 = seqstr;
code_freed = (char *)NULL;
comment_freed = (char *)NULL;
while (p0 != NULL)
{
	/* if there is a valid next, save its address */
	p1 = p0->next;
	/* Now deallocate arrays */
	if (p0->seq != NULL) cfree(p0->seq);
	if (p0->code != NULL)
	{
		/* Check to see if it was already freed */

		if (p0->code != code_freed)
		{
			code_freed = p0->code;
			cfree(p0->code);
		}
	}
/* debug
	printf("free (%d): code:%p comment:%p root:%p\n", n, p0->code, p0->comment, p0);
*/
	if (p0->comment != NULL)
	{
		/* Check to see if it was already freed */

		if (p0->comment != comment_freed)
		{
			comment_freed = p0->comment;
			cfree(p0->comment);
		}
	}
	if (p0->gap != NULL) cfree(p0->gap);
	if (p0->fdata != NULL) cfree(p0->fdata);
/*
	if (p0->assoc != NULL)
	{
	pres *p3;
		p2 = p0->assoc;
		while (p2 != NULL)
		{
			p3 = p2->next;
			cfree(p2->assoc);
			p2 = p3;
		}
		cfree(p0->assoc);
	}
*/
	free(p0);
	p0 = p1;
}
return;
}

void translate(rootres *dna, int nframe)
/*
rootres * dna : pointer to the sequence structure to use
int nfram : three or six frame translation, 0 three, 1 six
*/
{
/*
This is a 6-frame translator for ADSP scanning routines.
*/
int i;
int ii;
int size_trans;
rootres *temp;
char frame[50];
char *dna_ptr;
char *prot_ptr;
int ind[3];
char *ap;
char *it;
int a;
int last;

/* Calculate the length of the translation */
size_trans = (dna->numres/3) + 5; /* Defensive programming ? */
/* Initialise pointer to rootres */
temp = dna;
/* i keeps a check on the frame we are in, i+1 = actual frame */
for (i = 0; i < 3 ; i++)
{
	/* allocate a rootres structure for each translation */
	temp->next = (rootres *)malloc(sizeof(rootres));
	if (temp->next == (rootres *)NULL)
	{
		perror("rootres");
		exit(0);
	}
	/* temp is dna at the moment, make it the new structure */
	temp = temp->next;
	/* terminate the list here */
	temp->next = (rootres *)NULL;
	/* allocate seq array */
	temp->seq = (char *)calloc(size_trans+3,sizeof(char));
	if (temp->seq == NULL)
	{
		printf("*** calloc failed temp->seq\n");
		exit(0);
	}
	temp->gap = (int *)calloc(size_trans+3,sizeof(int));/* same length as seq */
	if (temp->gap == NULL)
	{
		printf("*** calloc failed temp->gap\n");
		exit(0);
	}
	temp->fdata = (float *)calloc(size_trans+3,sizeof(float));
	if (temp->fdata == NULL)
	{
		printf("*** calloc failed root->fdata\n");
		exit(0);
	}
	temp->assoc = (pres *)NULL;
	temp->comment = dna->comment; /* db comment is the same */
	temp->code = dna->code; /* db code is the same */
	/* Annotate the Frame */
	temp->type = DNA;
	temp->frame = i + 1;
	/* Initialise numres etc */
	dna_ptr = &dna->seq[i + 1]; /* Get in the right frame */
	temp->seq[0] = ':';
	prot_ptr = &temp->seq[1];
	for (ii = i + 1; ii < dna->numres; ii += 3)
	{
		for (a = 0; a < 3; a++)
		{
			ap = strchr(gcode_template, *dna_ptr);
			if (ap != (char *)NULL)
			{
				ind[a] = (char *)ap - (char *)gcode_template; 
			}
			else
			{
				printf("Null pointer will be referenced.\n");
				printf("%s, frame: %d\n", temp->code, i + 1);
				printf("'%c'\n", *dna_ptr);
				dump(dna);
				exit(1);
			}
			dna_ptr++; /* Next base of the codon */
		}
		/* Now do the look-up in gencode */
		it = gencode[ind[0]][ind[1]];
		*prot_ptr = it[ind[2]];
		prot_ptr++;
	}
	*prot_ptr = '\0';
	temp->numres = strlen(temp->seq);
	temp->numgap = 0;
}
/* Only execute this if 6 frame translation required */
if (nframe == THREE_FRAME) return;
/* Now the three reverse frames */
for (i = 0; i < 3 ; i++)
{
	/* allocate a rootres structure for each translation */
	temp->next = (rootres *)malloc(sizeof(rootres));
	if (temp->next == (rootres *)NULL)
	{
		perror("rootres");
		exit(0);
	}
	/* temp is dna at the moment, make it the new structure */
	temp = temp->next;
	/* terminate the list here */
	temp->next = (rootres *)NULL;
	/* allocate seq array */
	temp->seq = (char *)calloc(size_trans+3,sizeof(char));
	if (temp->seq == NULL)
	{
		printf("*** calloc failed temp->seq\n");
		exit(0);
	}
	temp->gap = (int *)calloc(size_trans+3,sizeof(int));/* same length as seq */
	if (temp->gap == NULL)
	{
		printf("*** calloc failed temp->gap\n");
		exit(0);
	}
	temp->fdata = (float *)calloc(size_trans+3,sizeof(float));
	if (temp->fdata == NULL)
	{
		printf("*** calloc failed root->fdata\n");
		exit(0);
	}
	temp->assoc = (pres *)NULL;
	temp->comment = dna->comment; /* db comment is the same */
	temp->code = dna->code; /* db code is the same */
	/* Annotate the Frame */
	/* For the reverse frames also store the number of bases here */
	temp->type = DNA;
	temp->frame = (i + 1) * -1;
	temp->bases = dna->numres;
	/* Initialise numres etc */
	dna_ptr = &dna->seq[i + 1]; /* Get in the right frame */
	temp->seq[0] = ':';
	/* Fill array from the end to the beginning */
	last = (dna->numres - i)/3;
	prot_ptr = &temp->seq[last];
	temp->seq[last + 1] = '\0';
	for (ii = i + 1; ii < dna->numres - 1; ii += 3)
	{
		for (a = 0; a < 3; a++)
		{
			ap = strchr(gcode_template6, *dna_ptr);
			if (ap != (char *)NULL)
			{
				ind[a] = (char *)ap - (char *)gcode_template6; 
			}
			else
			{
				printf("Null pointer will be referenced.\n");
				printf("%s, frame: %d\n", temp->code, i + 1);
				printf("'%c'\n", *dna_ptr);
				dump(dna);
				exit(1);
			}
			dna_ptr++; /* Next base of the codon */
		}
		/* Now do the look-up in gencode */
		it = gencode[ind[2]][ind[1]];
		*prot_ptr = it[ind[0]];
		prot_ptr--;
	}
	temp->numres = strlen(temp->seq);
	temp->numgap = 0;
}
return;
}/* End of translate */

void init_gencode()
{
int i, j;
char line[MAX_LINE];
char *p;
char *a;

if ((a = getenv("ADSP_GENETIC_CODE")) == (char *)NULL)
{
	fprintf(stderr,"- please set ADSP_GENETIC CODE environment variable.\n");
	exit(1);
}
if ((gcode = fopen(a,"r")) == NULL)
{
	perror(a);
	exit(1);
}
/* We can only guess the residue if the N is in the last base of a codon */
for (i = 0; i < 5; i++)
{
	for (j = 0; j < 5; j++)
	{
		fgets(line, MAX_LINE, gcode);
		if ((p = strchr(line,'\n')) != NULL) *p = '\0';
		gencode[i][j] = (char *) calloc(strlen(line) + 2, sizeof(char));
		if (gencode[i][j] == (char *)NULL)
		{
			perror("malloc gencode");
			exit(0);
		}
		strncpy(gencode[i][j], line, strlen(line));
	}
}
fclose(gcode);
return;
}/* End of init gencode */

void test_gencode()
{
int i, j;

for (i = 0; i < 4; i++)
{
	for (j = 0; j < 4; j++)
	{
		printf("%s\n", gencode[i][j]);
	}
}
return;
}/* End of test gencode */
