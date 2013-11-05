/*
Title           : SCAN_COMMAND
Module          : SCAN_COMMAND.C
Version         : 2.0
Author          : D.J.Parry-Smith
Date            : 6-JUNE-1989
History         : Bugs fixed to 28-FEB-1995
		  Incorporates plot as well as scan.
		  Segment input is ok.
		  Incorporates translation of DNA
		  Support for searching PRINTS format databases added
Description

SCAN_COMMAND contains the main controlling loop for the ADSP
database SCAN function.

It can read single or multiple table files or segment files.

*/

#ifdef VMS
#include stdio
#include string
#include unixio
#else
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <thread.h>
#include <synch.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#endif
#include "yheader.h"
#include "rr.h"

#define ADSP_VERSION "v3.1"
#define DCL_STRING 256
#define NCPU 8
#define MAX_FEATURES 2300
#define MAX_DATABASES 50
#define MAX_THREADS 30
#define MAX_SEG_LEN 50   /* The maximum length of each individual feature */
#define MAX_PRINTS_CODE_LEN 30
#define MAX_DB_CHUNK 500
#define SEGMENT_DATA 0
#define TABLE_DATA 1
#define PRINTS_DATA 2
#define SINGLE_SCAN 0    /* These are the scan styles */
#define PAIR_SCAN 1
#define NSINGLE_SCAN 2
#define NPAIR_SCAN 3
#define SINGLE 0
#define PAIR 1
#define SCAN_ALL 0
#define SCAN_P1 1
#define SCAN_FA 2
#define MAX_LINE BUFSIZ
#define DB_ARG 1  /* db list arg */
#define MT_ARG 2  /* method arg */
#define IN_ARG 3  /* in list arg */
#define FM_ARG 4  /* data format arg */
#define MD_ARG 5  /* scan modifier arg */
#define PL_ARG 6  /* plot determination */
#define TR_ARG 7  /* translation determination */
#define CT_ARG 8
#define NH_ARG 9  /* number of hits arg */
#define NC_ARG 10  /* count arg */
#define RT_ARG 11 /* results template arg */
/* if you want to use the old slow reading mechanism uncomment this
#define SLOW 1
*/
/* Data structures for the feature representation */

typedef struct fth fth;
typedef struct hlist hlist;
typedef struct link_residue link_residue;
typedef struct link_down link_down;

struct fth
{
/* Feature header structure */
	int len_feature;
	int max_score_feature;
	char file_spec[DCL_STRING];
	char header[80];
	char code[MAX_PRINTS_CODE_LEN];
	link_residue *start; 
	/* The address of the beginning of the feature structure */
};

struct hlist
{
	char code[20];            /*Pointer to sequence code */
	char comment[100];         /*Pointer to sequence comment */
	char segment[MAX_SEG_LEN];/*Pointer to the sequence segment */
	int start;                /*The sequence position of the start of the segment*/
	float score;              /*The percent score attained */
	int frame;		  /*Frame for translation mode */
	int dna_size;		  /*Size of DNA sequence for translation */
	hlist *previous;
	hlist *next;
};

struct link_residue
{
	int position;
	int pairs;
	link_residue *across;
	link_down *next;
};

struct link_down
{
	char res1;
	char res2;
	int sep;
	int score_good;
	link_down *next;
};

/* Global variables */

FILE *in_data;
FILE *in_data_1;
FILE *hit_data;
FILE *db_list;
FILE *plot_file;

fth *feature[MAX_FEATURES]; /* Can scan with up to MAX_FEATURES features */
hlist *hits[MAX_FEATURES];  /* start of each hitlist */
hlist *hite[MAX_FEATURES];   /* end of each hilist - used in merge_hits */
char *db_names[MAX_DATABASES]; /* maximum number of databases scanned */
int PLOT = FALSE;
int TRANSLATION = FALSE;
void scan_single();
void scan_pair();
void *scan_nsingle();
void scan_npair();
void max_score(int,int);
/* feature index, scan style */
void purge_seq();
void init_hitlist(int, int);
/* feature index, number of hits */
rootres *read_file(int, int);
rootres *read_file_quick(int, int, off_t);
/* number od sequences to be read, scan_modifier */
void merge_hits(rootres *,int);
/* feature index */
void type_hits(char *,char *,char *, int, int, int);
/* command line, feature list name, total sequences, total residues, method */
void read_segments(FILE *, char *, int, int);
/* file pointer, file spec, mode, feature index, data type */
void read_table(FILE *, char *, int, int);
/* file pointer, file spec, mode, feature index */
link_residue *allocate_across(int);
/* feature index */
void read_position_single(FILE *,int);
/* file pointer,feature index */
void read_position_pair(FILE *,int);
/* file pointer,feature index */
void purge_seq(rootres *);
/* sequence pointer */
void clear_gaps(rootres *);
/* sequence pointer */
void min_fix(int);
/* feature index */
int read_positions(FILE *, int);
/* returns -1 on end */
/* file pointer, feature index */
int read_header(FILE *, int);
/* Returns number of comments ignored */
/* file pointer, feature index */
int residues(rootres *);
int sequences(rootres *);
void type_prints_hits(char *, char *, int , int , int , int );
void translate(rootres *, int);
void purge_all(rootres *);
char *hitext;
char line[MAX_LINE];  
int total_features;
int num_hits; /* Now global because used in merge_hits */
int PRINTS = 1;
int cutoff = 0;
int gb_threshold;  /* How many contibuting matches before modification of score */
int scan_modifier;/* SCAN_ALL or SCAN_P1 only */
int error;
int total_residues = 0; /* total number of residues processed */
int seq_count = 0;
int frames = SIX_FRAME; /* three or six frame, THREE_FRAME, SIX_FRAME */
int ccount = 0;
int count;      /* Increment to indicate how many sequences have been processed */
int input_data_type;/* segments or table files */
/* Threads stuff */
mutex_t read_file_lock;
mutex_t seq_count_lock;
mutex_t merge_lock;
mutex_t translate_lock;
thread_t thread_id[MAX_THREADS];
/* mmap stuff */
char *src;
struct stat statbuf; /* for the fast db reading */


main(int argc, char **argv)
{
int scan_proc;  /* Scanning procedure */
int fc;/* Feature index counter */
int mode;/* SINGLE or PAIR mode */
char fn[DCL_STRING];/* database name */
char command[DCL_STRING]; /* Command line */
int i;/* total number of sequences processed */
int prev_i;  /* number of sequences processed from previous database*/
int prev_res;/* number of residues processed from previous database*/
/* FastA format databases now allowed, SCAN_FA */
int db_num; /* Number of databases scanned */
int th;
int din; /* file descriptor for fast db reading */

printf("ADSP Scan facility %s\n\n", ADSP_VERSION);
printf("Features: multithreaded, memory mapped databases.\n");
if (argc == 1)
{
	fprintf(stderr, "Usage:\n");
	fprintf(stderr, "scan dblist -method list_names prints|segments|tables all|p|fasta plot|scan six|three cutoff hits count results\n");
	fprintf(stderr, "dblist list of databases, one perline\n");
	fprintf(stderr, "method is one of s,p,n or a\n");
	fprintf(stderr, "(single, pair, nsingle, npair)\n");
	fprintf(stderr, "list of segment/table file names\n");
	fprintf(stderr, "keyword : segments or tables to indicate data format\n");
	fprintf(stderr, "          or use PRINTS database\n");
	fprintf(stderr, "keyword : all or just >P1 (p) sequences\n");
	fprintf(stderr, "keyword : write results to plotfile or merge\n");
	fprintf(stderr, "keyword : six|three frame translation of input\n");
	fprintf(stderr, "integer : only show hits above cutoff (PRINTS only)\n");
	fprintf(stderr, "number of hits to be listed\n");
	fprintf(stderr, "verification count of sequences processed\n");
	fprintf(stderr, "file extension for the created results files\n");
	fprintf(stderr, "- or the name of the plotfile to create\n");
	fprintf(stderr, "\n\n");
	exit(0);
}
else if (argc > 12)
{
	fprintf(stderr,"** Excess command line arguments ignored.\n");
}

if ((sscanf(argv[NH_ARG], "%d", &num_hits)) != 1)
{
	fprintf(stderr,"** Bad argument (number of hits) : `%s' not acceptable.\n",argv[NH_ARG]);
	exit(0);
}
if ((sscanf(argv[NC_ARG], "%d", &count)) != 1)
{
	fprintf(stderr,"** Bad argument (count) : `%s' not acceptable.\n",argv[NC_ARG]);
	exit(0);
}
if (strcmp(argv[MT_ARG], "-s") == 0)
{
	scan_proc = SINGLE_SCAN;
	mode = SINGLE;
	printf("Singlet scan chosen\n");
}
else if (strcmp(argv[MT_ARG], "-p") == 0)
{
	scan_proc = PAIR_SCAN;
	mode = PAIR;
	printf("Pairwise scan chosen\n");
}
else if (strncmp(argv[MT_ARG], "-n",2) == 0)
{
	scan_proc = NSINGLE_SCAN;
	mode = SINGLE;
	if (sscanf(argv[MT_ARG],"-n%d", &gb_threshold) != 1) gb_threshold = 0;
	printf("\nNsingle scan chosen with threshold of %d\n",gb_threshold);
}
else if (strncmp(argv[MT_ARG],"-a",2) == 0)
{
	scan_proc = NPAIR_SCAN;
	mode = PAIR;
	if (sscanf(argv[MT_ARG],"-a%d", &gb_threshold) != 1) gb_threshold = 0;
	printf("\nNpair scan chosen with threshold of %d\n",gb_threshold);
}
else 
{
	fprintf(stderr,"\n** Invalid option `%s' - unknown scanning method !\n",argv[MT_ARG]);
	exit(0);
}
/* Read segments or tables ? */
if (strcmp(argv[FM_ARG], "segments") == 0)
{
	input_data_type = SEGMENT_DATA;
}
else if (strcmp(argv[FM_ARG], "tables") == 0)
{
	input_data_type = TABLE_DATA;
}
else if (strcmp(argv[FM_ARG], "prints") == 0)
{
	input_data_type = PRINTS_DATA;
	PRINTS = 0;
}
else
{
	fprintf(stderr,"** Unrecognised input data format specified - `%s' !\n",argv[FM_ARG]);
	exit(0);
}
/* Scan all sequences or only >P1 sequences */
if (strcmp(argv[MD_ARG], "all") == 0)
{
	scan_modifier = SCAN_ALL;
	printf("Scanning all sequences\n");
}
else if (strcmp(argv[MD_ARG], "fasta") == 0)
{
	scan_modifier = SCAN_FA;
	printf("Scanning sequences in FASTA format\n");
}
else
{
	scan_modifier = SCAN_P1;
	printf("Scanning only >P1 sequences\n");
}

/* Determine whether plotting results or merging for database scan */
if (strcmp(argv[PL_ARG], "plot") == 0)
{
	PLOT = TRUE;
	printf("Unmerged results go to plotfile '%s'\n", argv[RT_ARG]);
	if ((plot_file  = fopen(argv[RT_ARG], "w")) == NULL)
	{
		fprintf(stderr, "Error opeining file '%s'\n", argv[RT_ARG]);
		perror(argv[RT_ARG]);
		exit(0);
	}
}
else PLOT = FALSE;
/* New: 21 February 1995 */
/* Determine whether translation is required */
if (strcmp(argv[TR_ARG],"six") == 0)
{
	TRANSLATION = TRUE;
	frames = SIX_FRAME;
	init_gencode();
	printf("6 frame translation will be performed on nucleotide data\n");
}
else if (strcmp(argv[TR_ARG],"three") == 0)
{
	TRANSLATION = TRUE;
	frames = THREE_FRAME;
	init_gencode();
	printf("3 frame translation will be performed on nucleotide data\n");
}
else TRANSLATION = FALSE;
/* Get the cutoff value for PRINTS */
if (sscanf(argv[CT_ARG],"%d", &cutoff) != 1)
{
	fprintf(stderr,"- expecting integer for cutoff value, '%s'.",argv[CT_ARG]);
	exit(1);
}

/* The file extension to use for the hit files */
hitext = argv[RT_ARG]; /* The file extension to use for the hit files */
/* Open the database */
/* First check for user defined database */
/* The file pointed to is a list of databases to be searched */
strcpy(fn,argv[DB_ARG]);
if ((db_list = fopen(fn,"r")) == NULL)
{
	fprintf(stderr, "** Unable to open database list <%s>\n",fn);
	exit(0);
}

/* Open the list of input file names */
if (input_data_type == PRINTS_DATA)
{
char *a;
	/* Open the PRINTS database */
	if ((a = getenv("ADSP_PRINTS_FILE")) == NULL)
	{
		fprintf(stderr,"- please setenv ADSP_PRINTS_FILE.\n");
		exit(1);
	}
	if ((in_data = fopen(a,"r")) == NULL)
	{
		perror(a);
		exit(1);
	}
}
else if ((in_data = fopen(argv[IN_ARG],"r")) == NULL)
{
	fprintf(stderr,"** Unable to open input file `%s' !",argv[IN_ARG]);
	exit(0);
}
/* Now read all the segment/table files */
if (input_data_type == PRINTS_DATA)
{
	printf("Reading PRINTS");
	fflush(stdout);
	fc = read_prints(mode, scan_proc);
	printf("..%d fingerprint elements\n", fc);
	fclose(in_data);
}
else
{
fc = 0;
while (fgets(line, MAX_LINE, in_data) != NULL)
{
char *p;
	if (line[0] == '%') continue; /* Ignore comment lines */
	/* Get rid of the newline from fgets */
	p = line;
	while (*p != '\0')
	{
		if (*p == '\n') *p = '\0';
		p++;
	}
	if ((in_data_1 = fopen(line,"r")) == NULL)
	{
		fprintf(stderr,"** Unable to open input file `%s' !",line);
		exit(0);
	}
	fc++;
	if (input_data_type == SEGMENT_DATA)
	{
		read_segments(in_data_1, line, mode, fc);
	}
	else
	{
		/* Must be table file format */
		read_table(in_data_1, line, mode, fc);
	}
	fclose(in_data_1);
	/* The max score achievable depends on whether modified method is used */
	max_score(fc, scan_proc);
	/* Now initialise a hitlist for the number of hits required */
	/* if (PLOT == FALSE) */ init_hitlist(fc, num_hits);
	/* Get the next file name */
}
printf("\n%d input files\n",fc);
}
total_features = fc;
total_residues = 0;
ccount = 0;
i = 0;
prev_i = 0;
prev_res = 0;
db_num = 0;
printf("\nInitialisation complete. Scanning..\n");
fflush(stdout); /* Try to force a write on the output file, may not work on VMS*/
while (fgets(line, MAX_LINE, db_list) != NULL)
{
char *p;

	p = line;
	while (*p != '\0')
	{
		if (*p == '\n') *p = '\0';
		p++;
	}
	strcpy(fn,line);
	if (scan_modifier != SCAN_FA) strcat(fn,".seq");
#ifdef SLOW
	if ((data = fopen(fn,"r")) == NULL)
#else
	if ((din = open(fn, O_RDONLY)) < 0)
#endif
	{
		fprintf(stderr, "\n** Cannot access %s as database\n",fn);
		continue; /* Next database ? */
	}
#ifndef SLOW
	if (fstat(din, &statbuf) < 0)
	{
		fprintf(stderr,"fstat error on database\n");
		exit(1);
	}
	/* open the file mmap'd */
	if ((src = mmap(0, statbuf.st_size, PROT_READ, MAP_SHARED, din, 0)) == (caddr_t) -1)
	{
		fprintf(stderr, "mmap error for input database\n");
		exit(1);
	}
#endif
	printf("database : %s\n",fn);
	db_names[db_num] = (char *)malloc(strlen(fn)+1);
	if (db_names[db_num] == NULL)
	{
		fprintf(stderr,"\n** memory corruption error in main (db_names) !\n");
		exit(0);
	}
	strcpy(db_names[db_num], fn);
	db_num++;
	/* thr_setconcurrency(NCPU); */
	for (th = 0; th < MAX_THREADS; th++)
	{
		error = thr_create((void *)NULL, (size_t)NULL, scan_nsingle, (void *)NULL, THR_NEW_LWP, &thread_id[th]) ;
		if (error != 0)
		{
			fprintf(stderr,"thr_create returned non-zero (%d)\n", error);
			perror("thr_create");
			exit(error);
		}
	}
	printf("Concurrency set to: %d\n", thr_getconcurrency());
/* Now catch the threads as they terminate */
for (th = 0; th < MAX_THREADS; th++)
{
	thr_join((thread_t)0, NULL, NULL);
}
	fclose(data);
	printf("\n%d sequences processed, %d residues\n",seq_count-prev_i,total_residues-prev_res);
	fflush(stdout);
	prev_i = seq_count;
	prev_res = total_residues;
}/* End while list of databases */
printf("\n%d sequences processed, %d residues in total.\n",seq_count,total_residues);
printf("%d databases scanned.\n",db_num);
fflush(stdout);
if (PLOT == FALSE)
{
printf("Writing out hitlists..\n");
sprintf(command,"%s %s %s %s %s %s %s %s %s",
argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8]);
if (input_data_type == PRINTS_DATA)
{
	type_prints_hits(command, argv[DB_ARG],cutoff, seq_count, total_residues, scan_proc);
}
else
{
	type_hits(command,argv[DB_ARG],argv[IN_ARG],seq_count ,total_residues,scan_proc);
}
}
if (PLOT == TRUE)
{
	fclose(plot_file);
	printf("\nEnd of PLOT\n");
}
else printf("\nEnd of SCAN\n");
}/* End of main */



void max_score(int f_ind, int style)
{
/* Calculate the maximum score that can be achieved by a feature */
int dumb;
link_residue *current;
link_down *down;
int last_score;

last_score = 0; 
current = feature[f_ind]->start;
dumb = 0;
while (current != NULL)
{
	down = current->next;
	while (down != NULL)
	{
		if (down->score_good > last_score) last_score = down->score_good;
		down = down->next;
	}
	dumb += last_score;
	current = current->across;
	last_score = 0;
}
if (style > PAIR_SCAN) dumb *= feature[f_ind]->len_feature;
feature[f_ind]->max_score_feature = dumb;
/*
printf("\nMaximum feature score is %d\n",feature[f_ind]->max_score_feature);
*/
return;
}/* End of max_score */


void read_segments(FILE * fp, char * fsp, int d_format, int f_ind)
{
/* Read segments from a segment file into ADSP data structures
 
Arguments:

file pointer to the stream opened in another function
file specification for storage
data format (single or pairwise)
feature index
*/
char *p;

/* Make a new feature header structure and add to the list */
feature[f_ind] = (fth *) malloc(sizeof(fth));
demand(feature[f_ind], malloc failed in read_segments);
/* New for PRINTS version */
feature[f_ind]->code[0] = '\0';
strcpy(feature[f_ind]->file_spec,fsp); /* fsp is really a pointer to line ! */
fgets(line, MAX_LINE, fp);
while (line[0] == '%') fgets(line, MAX_LINE, fp);
sscanf(line, "%d", &feature[f_ind]->len_feature);
/* Next line is the formal title line */
fgets(line, MAX_LINE, fp);
/* Remove newline generated by fgets */
p = line;
while (*p != '\0')
{
	if (*p == '\n') *p = '\0';
	p++;
}
line[80] = '\0'; /* Make certain string is terminated */
strcpy(feature[f_ind]->header,line);
printf("\n%s\n",feature[f_ind]->file_spec);
printf("%s\n",feature[f_ind]->header);
feature[f_ind]->start = allocate_across(f_ind);
/* f_ind tells allocate across which feature to use */
if (d_format == SINGLE) read_position_single(fp,f_ind);
else read_position_pair(fp,f_ind);
return;
}/*End of read_segments */


void read_table(FILE *fp, char *fsp, int d_format, int f_ind)
{
/*
Read data from table file into ADSP data structures

Arguments:

file pointer to the stream opened in another function
file specification for storage
data format (single or pairwise)
feature index
*/

printf("Reading table file\n");
/* Make a new feature header structure and add to the list */
feature[f_ind] = (fth *) malloc(sizeof(fth));
demand(feature[f_ind], malloc failed in read_segments);
strcpy(feature[f_ind]->file_spec,fsp); /* fsp is really a pointer to line ! */
read_header(fp, f_ind);
printf("\n%s\n",feature[f_ind]->file_spec);
printf("%s\n",feature[f_ind]->header);
feature[f_ind]->start = allocate_across(f_ind);
for (;;)
{
	if (read_positions(fp, f_ind) == -1) break;
}
min_fix(f_ind);
return;
}/* End of read table file */


int read_header(FILE *fp, int f_ind)
{
/* Arguments -  file pointer, feature index */
/* Read the standard motif file header */
/* At present this just consists of the title of the motif */
/* comments and overall length of the motif */

int comments;
char *fgets();
char *p;

fgets(line,MAX_LINE,fp);
p = line;
while (*p != '\0')
{
	if (*p == '\n') *p = '\0';
	p++;
}
line[80] = '\0'; /* Make certain string is terminated */
strcpy(feature[f_ind]->header,line);
comments = 0;
fgets(line,MAX_LINE,fp);
while(line[0] == '%')
{
	fgets(line,LINE_LEN,fp);
	comments++;
}
sscanf(line,"%d",&feature[f_ind]->len_feature);
printf("Feature length is %d.\n",feature[f_ind]->len_feature);
return(comments);
}


link_residue *allocate_across(int f_ind)
{
/* Allocate dynamic memory for across structure */
link_residue *result;
link_residue *previous;
link_residue *addr_across;
int i;

result = (link_residue *) malloc(sizeof(link_residue));
demand(result, malloc failed in allocate_across);
result->position = 0;
result->next = (link_down *)NULL;
result->pairs = 0;
previous = result;
addr_across  = result; /* Don't loose the start of the list */
for (i=1; i<= feature[f_ind]->len_feature; i++)
{
	result = (link_residue *) malloc(sizeof(link_residue));
	demand(result, malloc failed in allocate_across);
	result->position = i;
	previous->across = result;
	result->next = (link_down *)NULL;
	result->pairs = 0; /* ->pairs keeps total pairs to be scanned in this sector */
	previous = result;
}
/* Put null pointer in the last structure */
result->across = (link_residue *)NULL;
return(addr_across);
}/* End of allocate_across */


void read_position_single( FILE *fp, int f_ind)
{
/* Reads the motif entry for each position which has a weight in the motif */
/* new structures of type link_down are created for each pair */
/* The next record in the file must contain the position number, starting at 0 */
int i;
int j;
int num_motifs;
int num_ignored;
short found_it;
link_residue *link;
link_down *current;
link_down *res;
link_down *list_end;

num_motifs = 0;
num_ignored = 0;
printf("read_position_single\n");
while (fgets(line, LINE_LEN, fp) != NULL)
{
link = feature[f_ind]->start; /* The addr of the first link */
if (strspn(line,"TCAGNSPFLYHQVKDEIWRMBXZ-") != 	feature[f_ind]->len_feature)
{
	printf("%-.*s ignored\n",feature[f_ind]->len_feature,line);
	num_ignored++;
	continue;
}
	num_motifs++;
	for (i=0; i<feature[f_ind]->len_feature; i++)
	{
		if (line[i] != '-')
		{
		for (j=i; j<= i; j++) /* Only one iteration necessary */
		{
			found_it = 0;
			if (line[j] == '-') continue;
			/* If there is already a node for this sep with these residues, then use it */
			if (link->next != NULL)
			{
				current = link->next;
				while (current != NULL)
				{
					if (current->sep != j-i)
					{
						current = current->next;
					}
					else
					{
						if ((current->res2 == line[j]) && (current->res1 == line[i]))
						{
						current->score_good++;
						found_it = 1;
						break;
						}
						else
						{
						current = current->next;
						}
					}
				}
			}
			if (found_it == 0)
			{
			res = (link_down *) malloc(sizeof(link_down));
			demand(res,malloc failed in read_position_single);
			/* DJPS 6 April 1994 */
			res->score_good = 0;
			res->next = NULL;
			if (link->pairs == 0)
			{
				link->next = res;
			}
			else
			{
			/* Hop to the end of the list */
			list_end = link->next;
			while (list_end->next != NULL) list_end = list_end->next;
			list_end->next = res;
			}
			res->res1 = line[i];
			res->res2 = line[j];
			res->sep = j-i;
			res->score_good++;
			link->pairs++;
			} /* End an if */
		}/* End for j */
		}/* End if*/
		link = link->across;
	}
/* now read the next motif */
} /* End while fgets */
printf("\n%d motifs were read, %d accepted, %d ignored\n",num_motifs+num_ignored,
	num_motifs,num_ignored);
return;
}/* read_position_single */


void read_position_pair(FILE *fp, int f_ind)
{
/* Reads the motif entry for each position which has a weight in the motif */
/* new structures of type link_down are created for each pair */
/* The next record in the file must contain the position number, starting at 0 */


int i;
int j;
int num_motifs;
int num_ignored;
short found_it;
link_residue *link;
link_down *current;
link_down *res;
link_down *list_end;

num_motifs = 0;
num_ignored = 0;
printf("read_position_pair\n");
while (fgets(line, LINE_LEN, fp) != NULL)
{
link = feature[f_ind]->start; /* addr of first link */
if (strspn(line,"TCAGNSPFLYHQVKDEIWRMBXZ-") != 	feature[f_ind]->len_feature)
{
	printf("%-.*s ignored\n",feature[f_ind]->len_feature,line);
	num_ignored++;
	continue;
}
	num_motifs++;
for (i=0; i<feature[f_ind]->len_feature; i++)
	{
		if (line[i] != '-')
		{
		for (j=i+1; j< feature[f_ind]->len_feature; j++)/*Don't uses 0 separation */
		{
			found_it = 0;
			if (line[j] == '-') continue;
			/* If there is already a node for this sep with these residues, then use it */
			if (link->next != NULL)
			{
				current = link->next;
				while (current != NULL)
				{
					if (current->sep != j-i)
					{
						current = current->next;
					}
					else
					{
						if ((current->res2 == line[j]) && (current->res1 == line[i]))
						{
						current->score_good++;
						found_it = 1;
						break;
						}
						else
						{
						current = current->next;
						}
					}
				}
			}
			if (found_it == 0)
			{
			res = (link_down *) malloc(sizeof(link_down));
			demand(res,malloc failed in read_position_pair);
			if (link->pairs == 0)
			{
				link->next = res;
			}
			else
			{
			/* Hop to the end of the list */
			list_end = link->next;
			while (list_end->next != NULL) list_end = list_end->next;
			list_end->next = res;
			}
			res->res1 = line[i];
			res->res2 = line[j];
			res->sep = j-i;
			res->score_good++;
			link->pairs++;
			} /* End an if */
		}/* End for j */
		}/* End if*/
		link = link->across;
	}
/* now read the next motif */
} /* End while fgets */
printf("\n%d motifs were read, %d accepted, %d ignored\n",num_motifs+num_ignored,
	num_motifs,num_ignored);
return;
}/* read_position_pair */


int read_positions(FILE *fp, int f_ind)
{
/* Arguments file pointer, feature index */
/* Reads the motif entry for each position which has a weight in the motif */
/* new structures of type link_down are created for each pair */

/* The next record in the file must contain the position number, starting at 0 */


int start_at;
int i;
link_residue *link;
link_down *res;
link_down *previous;
char *fgets();

if (fgets(line,MAX_LINE,fp) == NULL)
{
	return(-1);
}
else
{
/* Which starting position are we talking about */
sscanf(line,"%d",&start_at);
/* Now generate the first link_down structure */
link = feature[f_ind]->start;
if (start_at != 0)
{
/* must hop across links to find required structure */
for (i=0; i< start_at; i++)
{
	link = link->across;
}
if (link->position != start_at)
{
	fprintf(stderr,"** internal data structure not compatible with input data !\n");
	fprintf(stderr,"** expecting link->position = %d but found %d !\n",
		link->position,start_at);
}
}/* End if */
link->pairs = 0;
/* The next few records contain the data for structure link_down */
/* starting at start_at */
fgets(line,MAX_LINE,fp);
while(line[0] != '!')
{
	res = (link_down *) malloc(sizeof(link_down));
	demand(res,malloc failed in read_positions);
	if (link->pairs == 0)
	{
		link->next = res;
		previous = res;
	}
	else
	{
		previous->next = res;
		previous = res;
	}
	sscanf(line,"%c %c %d %d",&res->res1,&res->res2,
		&res->sep,&res->score_good);
	fgets(line,LINE_LEN,fp);
	link->pairs++;
}/* End while */
res->next = 0; /* End of list */
return(0);
}/* End main if */
/*return(0);*/
}/* End of read_positions */


void min_fix(int f_ind)
{
/* Argument - feature index */
/* 
ADSP does not allow negative scores in table files, this function
adds a constant value to make all scores start at zero.
*/
link_residue *current;
link_down *down;
int last_score;

/* find the largest negative value in the motif */
last_score = 0;
current = feature[f_ind]->start;
while(current != NULL)
{
	down = current->next;
	while (down != NULL)
	{
		if (down->score_good < last_score) last_score = down->score_good;
		down = down->next;
	}
        current = current->across;
}
if (last_score >= 0) return; /* No fix to apply in this case */
last_score = abs(last_score);
printf("Correction is %d\n",last_score);
current = feature[f_ind]->start;
while (current != NULL)
{
	down = current->next;
	while (down != NULL)
	{
		down->score_good += last_score;
		down = down->next;
	}
	current = current->across;
}
printf("Feature scores corrected - negative scores not allowed.\n");
return;
}/* End of min_fix */



/* Functions for dealing with hitlists */

void init_hitlist(int f_ind, int n_hits)
{
/* Generate a new hitlist structure of the required size */
hlist *current;
int i;

/*
printf("Creating hitlist number %d\n",f_ind);
*/
hits[f_ind] = (hlist *)malloc(sizeof(hlist));
if (hits[f_ind] == NULL)
{
	printf("\n*** malloc failed in init_hitlist\n");
	exit(0);
}
hits[f_ind]->code[0] = '\0'; /* We test for empty code later on */
hits[f_ind]->previous = NULL;
current = hits[f_ind];
current->score = 0.0;
current->frame = 0;
current->dna_size = 0;
for (i=1; i< n_hits; i++)
{
	current->next = (hlist *)malloc(sizeof(hlist));
	if (current->next == NULL)
	{
		printf("\n*** malloc failed in init_hitlist\n");
		exit(0);
	}
	current->next->previous = current;
	current = current->next;
	current->code[0] = '\0';
	current->score = 0.0;
	current->frame = 0;
	current->dna_size = 0;
}
hite[f_ind] = current; /* Used in merge_hits */
current->next = NULL;
return;
}/* End of init_hitlist */


#define HITSCORE mres->fdata[val]

void merge_hits(rootres *mres, int f_ind)
{
/* Merge a list of hits from space into the hitlist */
int max_hit(rootres *);
int opt_num;
int val;
hlist *current;
hlist *inserted;
int ok;
void type_hits(char *,char *,char *, int, int, int);
/* command line, feature list name, total sequences, total residues, method */
int loop;

opt_num = mres->numres;
ok = 0;
loop = 0;
current = hits[f_ind]; /* Start at the top of the list */
/*
while ((val = max_hit()) != 0 )
*/
val = max_hit(mres);
/*
for (val = 1; val <= mres->numres; val++)
*/
if (val != 0)
{
	/* Some sensitivity checks */
	if (HITSCORE <= hite[f_ind]->score) return;
	if ((int)HITSCORE < cutoff) return;

	current = hits[f_ind];
	while(current != NULL)
	{
	if (current->score >= HITSCORE) current = current->next;
	else 
	{
	/* Score of current is less than HITSCORE so insert hitscore here */
	/* Get a spare node from the end of the hitlist */
	inserted = hite[f_ind];
	if (current != hite[f_ind])
	{
	hite[f_ind] = hite[f_ind]->previous;
	hite[f_ind]->next = NULL;
	inserted->previous = current->previous;
	if (current->previous == NULL) 
	{
		hits[f_ind] = inserted; /* This is root node */
		inserted->next = current;
		current->previous = inserted;
	}
	else
	{
		inserted->next = current;
		current->previous->next = inserted;
		current->previous = inserted;
	}
	if (current->previous == current->next)
	{
		printf("\nDud node : current\n");
		type_hits("n/s","n/s","n/s",0,0,0);
		exit(0);
	}
	if (inserted->previous == inserted->next)
	{
		printf("\nDud node : inserted\n");
		type_hits("n/s","n/s","n/s",0,0,0);
		exit(0);
	}
	}
	/* Node has been inserted -- now fill it */
	strcpy(inserted->code,mres->code);
 	strcpy(inserted->comment,mres->comment);
	inserted->comment[80] = '\0';
	strncpy(inserted->segment,&mres->seq[val],feature[f_ind]->len_feature);
	inserted->segment[feature[f_ind]->len_feature] = '\0';
	inserted->score = HITSCORE;
	inserted->start = val;
	if (TRANSLATION == TRUE)
	{
		inserted->frame = mres->frame;
		inserted->dna_size = mres->bases;
	}
	ok = 1;
	}
	if (ok == 1) break;
	} /* End while */
ok = 0;
}
return;
}/* End of merge_hits */

void type_prints_hits(char *command, char *db, int cut, int num_seqs,
			int num_res, int method)
{

hlist *current;
int f_ind;
int i;
int db_inc;
int dna_start;
int frame;
char strmethod[10];

if ((hit_data = fopen("prints.hits","w")) == NULL)
{
	fprintf(stderr,"- could not open file 'prints.hits'.");
	perror("prints.hits");
	exit(1);
}

if (method == SINGLE_SCAN) strcpy(strmethod,"Single");
else if (method == PAIR_SCAN) strcpy(strmethod,"Pair");
else if (method == NSINGLE_SCAN) strcpy(strmethod,"Nsingle");
else if (method == NPAIR_SCAN) strcpy(strmethod,"Npair");
fprintf(hit_data,"ADSP Scan (%s option) %s\n\n",strmethod, ADSP_VERSION);
fprintf(hit_data,"Command line was     : %s\n",command);
db_inc = 0;
fprintf(hit_data,"Database(s) scanned  : %s\n",db_names[db_inc++]);
while (db_names[db_inc] != NULL)
{
	fprintf(hit_data,"                       %s\n",db_names[db_inc++]);
}
fprintf(hit_data,"Number of sequences  : %d\n",num_seqs);
fprintf(hit_data,"Number of residues   : %d\n",num_res);
fprintf(hit_data,"PRINTS cutoff        : %d\n",cut);


f_ind = 1;
while(f_ind < MAX_FEATURES && feature[f_ind] != NULL)
{
int prints_code;

	prints_code = 0;
	current = hits[f_ind];
	i = 1;
	if (TRANSLATION == FALSE)
	{
	while (current != NULL)
	{
		if (current->code[0] == '\0') break;
		if ((int)current->score < cut) break;
		else
		{
			/* put out prints code if not already done */
			if (prints_code != 1)
			{
				fprintf(hit_data,"\n\n%s\n",feature[f_ind]->code);
				prints_code = 1;
			}
		}
		fprintf(hit_data,"%4d) %6.2f    %-14.12s%4d - %4d   %s\n",i++,current->score,
		current->code,current->start, 
		current->start+feature[f_ind]->len_feature-1,current->segment);
		current = current->next;
	}
	}
	else /* TRANSLATION must be true */
	{
	while (current != NULL)
	{
		if (current->code[0] == '\0') break;
		if ((int)current->score < cut) break;
		else
		{
			/* put out prints code if not already done */
			if (prints_code != 1)
			{
				fprintf(hit_data,"\n\n%s\n",feature[f_ind]->code);
				prints_code = 1;
			}
		}
		/* calculate the starting point in the DNA */
		/* Get the frame first */
		frame = current->frame;
		/* deal with +1, +2 and +3 first */
		if (frame > 0)
		{
		fprintf(hit_data,"%4d) %6.2f    %-16.14s%5d - %5d %2d %s\n",
			i++,
			current->score,
			current->code,
			(current->start * 3) - (3 - frame),
			((current->start * 3) - (3 -frame)) + ((feature[f_ind]->len_feature-1) * 3),
			current->frame,
			current->segment);
			current = current->next;
		}
		else /* Frame < 0 */
		{
		/* Extract dna length from current */
		dna_start = current->dna_size;
		fprintf(hit_data,"%4d) %6.2f    %-16.14s%5d - %5d %2d %s\n",
			i++,
			current->score,
			current->code,
			dna_start +1 - (current->start * 3) - (3 - frame),
			(dna_start +1 - (current->start * 3) -(3 -frame)) - ((feature[f_ind]->len_feature-1) * 3),
			current->frame,
			current->segment);
			current = current->next;
		}
		
	}
	}
	/* Now for the titles (do we want title for this ? */
	f_ind++;
}
return;
}/* End of type prints hits */

void type_hits(char *command, char *db, char *motif_file,
		int num_seqs, int num_res, int method)
{
/* print out the hitlists from beginning to end */
/* Caller passes the complete command line that invoked scan and the name of the
motif file separately
*/
hlist *current;
int i;
char line[512];
char u_temp[DCL_STRING];
char strmethod[10];
char hitfile[DCL_STRING];
int f_ind;
char *p;
char *q;
char *slash;
int db_inc;
int dna_start;
int frame;

f_ind = 1;
while(f_ind < MAX_FEATURES && feature[f_ind] != NULL)
{
	/* Generate the new hitfile name */

	/* The following is necessary because of the awkward way VMS files*/
	/* are named. The hitlist is called by the name of the segment */
	/* or table file used to generate it. Its device, directory and */
	/* file type are taken from a template file specification supplied */
	/* on the command line */

	/* Copy the device, or logical name, if there is one */
	p = hitext;
	q = hitfile;
#ifdef VMS
	if (strchr(hitext, ':') != NULL)
	{
		while (*p != ':')
		{
			*q++ = *p++;
		}
		*q++ = *p++; /* Copy the ':' as well */
	}
	/* Copy the directory from hitext */
	if (*p == '[') /* There is a directory */
	{
		while (*p != ']')
		{
			*q++ = *p++;
		}
		*q++ = *p++; /* Copy the ']' as well */
	}
	/* Ignore any file name - that is supplied from the input file name */
	/* Now copy the file name from input file */
	p = strchr(feature[f_ind]->file_spec,']'); /* Is there a directory */
	if (p == NULL)
	{
		p = strchr(feature[f_ind]->file_spec,':');
	}
	if (p == NULL) p = feature[f_ind]->file_spec;
	else p++; /* Ignore : or ] */
	/* p now points to the beginning of the file name */
	while (*p != '.')
	{
		*q++ = *p++;
	}
	p = strrchr(hitext,'.');
	if (p == NULL)
	{
		printf("No file type supplied\n");
		*q++ = '.';
		*q = '\0';
	}
	else
	{
		while (*p != '\0')
		{
			if (*p == ';') break;
			*q++ = *p++;
		}
		*q = '\0';
	}
#endif
#ifdef unix
	/* Deal with a unix file name */
	/* Find the last `/' character*/
	slash = strrchr(hitext,'/');
	if (slash != NULL)
	{
		while (p != slash)
		{
			*q++ = *p++;
		}
		*q++ = *p++; /* Copy the '/' as well */
	}
	/* Now get the input file name */
	p = strrchr(feature[f_ind]->file_spec,'/');
	if (p == NULL)
	{
		p = feature[f_ind]->file_spec;
	}
	else p++; /* Ignore the '/'*/
	/* p now points to the beginning of the file name */
	while (*p != '.')
	{
		if (*p == '\0') break;
		*q++ = *p++;
	}
	p = strrchr(hitext,'.');
	if (p == NULL)
	{
		printf("No file type supplied\n");
		*q = '\0';
	}
	else
	{
		while (*p != '\0')
		{
			*q++ = *p++;
		}
		*q = '\0';
	}	
#endif
	if ((hit_data = fopen(hitfile,"w")) == NULL)
	{
		printf("\n\nCould not open output file <%s>!\n",hitfile);
		printf("Will attempt to write output to unique filename.\n");
		strcpy(u_temp,"scnXXXXXX");
		if (mktemp(u_temp) == NULL)
		{
			fprintf(stderr,"Unable to create unique filename !\n");
			fprintf(stderr,"Exiting\n");
		}
		else
		{
			if ((hit_data = fopen(u_temp,"w")) == NULL)
			{
				fprintf(stderr,"Could not open output file <%s>!\n",u_temp);
				fprintf(stderr,"Exiting\n");
			}
			else
			{
				printf("%s\n",u_temp);
			}
		}
	}
	else
	{
		printf("%s\n",hitfile);
	}
	if (method == SINGLE_SCAN) strcpy(strmethod,"Single");
	else if (method == PAIR_SCAN) strcpy(strmethod,"Pair");
	else if (method == NSINGLE_SCAN) strcpy(strmethod,"Nsingle");
	else if (method == NPAIR_SCAN) strcpy(strmethod,"Npair");
	fprintf(hit_data,"ADSP Scan (%s option) %s\n\n",strmethod, ADSP_VERSION);
	fprintf(hit_data,"Command line was     : %s\n",command);
	db_inc = 0;
	fprintf(hit_data,"Database(s) scanned  : %s\n",db_names[db_inc++]);
	while (db_names[db_inc] != NULL)
	{
	fprintf(hit_data,"                       %s\n",db_names[db_inc++]);
	}
	fprintf(hit_data,"Number of sequences  : %d\n",num_seqs);
	fprintf(hit_data,"Number of residues   : %d\n",num_res);
	fprintf(hit_data,"Feature file         : %s\n",feature[f_ind]->file_spec);
	fprintf(hit_data,"Feature file title   : %-.80s\n",feature[f_ind]->header);
	fprintf(hit_data,"Part                 : %d of %d\n\n",f_ind,total_features);
	fprintf(hit_data,"Summary listing of hits follows :\n\n");
if (TRANSLATION == FALSE)
{
	fprintf(hit_data,"      %%SCORE    NAME           FROM    TO   SEQUENCE\n");
	fprintf(hit_data,"      ------    ----           ----    --   --------\n\n");
}
else
{
	fprintf(hit_data,"      %%SCORE    NAME             FROM      TO FM SEQUENCE\n");
	fprintf(hit_data,"      ------    ----             ----      -- -- --------\n\n");
}
	/* That was to be compatible with Dix's Scanning program */
	current = hits[f_ind];
	i = 1;
	if (TRANSLATION == FALSE)
	{
	while (current != NULL)
	{
		if (current->code[0] == '\0') break;
		fprintf(hit_data,"%4d) %6.2f    %-14.12s%4d - %4d   %s\n",i++,current->score,
		current->code,current->start, 
		current->start+feature[f_ind]->len_feature-1,current->segment);
		current = current->next;
	}
	}
	else /* TRANSLATION must be true */
	{
	while (current != NULL)
	{
		if (current->code[0] == '\0') break;
		/* calculate the starting point in the DNA */
		/* Get the frame first */
		frame = current->frame;
		/* deal with +1, +2 and +3 first */
		if (frame > 0)
		{
		fprintf(hit_data,"%4d) %6.2f    %-16.14s%5d - %5d %2d %s\n",
			i++,
			current->score,
			current->code,
			(current->start * 3) - (3 - frame),
			((current->start * 3) - (3 -frame)) + ((feature[f_ind]->len_feature-1) * 3),
			current->frame,
			current->segment);
			current = current->next;
		}
		else /* Frame < 0 */
		{
		/* Extract dna length from current */
		dna_start = current->dna_size;
		fprintf(hit_data,"%4d) %6.2f    %-16.14s%5d - %5d %2d %s\n",
			i++,
			current->score,
			current->code,
			dna_start +1 - (current->start * 3) - (3 - frame),
			(dna_start +1 - (current->start * 3) -(3 -frame)) - ((feature[f_ind]->len_feature-1) * 3),
			current->frame,
			current->segment);
			current = current->next;
		}
		
	}
	}
	/* Now for the titles */
	fprintf(hit_data,"\n\n\n");
	current = hits[f_ind];
	i = 1;
	while (current != NULL)
	{
		if (current->code[0] == '\0') break;
		fprintf(hit_data,"%5d  %-14.12s  %.80s\n",i++,current->code,current->comment);
		current = current->next;
	}
	fclose(hit_data);
	f_ind++;
}/* End while f_ind */
return;
} /* End of type hits */



/* Scanning functions go here */

void *scan_nsingle()
{
/* Scan_nsingle - modified singlet scan 

Arg: threshold indicates the number of contributing matches that must 
be exceeded before score modification takes place. Normally 0.

The hitlist is modified by multiplying
the score achieved for each segment by the number of matches with the matrix
for that segment 
*/

int num_matches;
int num_scans;
link_residue *current_residue;
link_down *current_down;
rootres *mine;
rootres *mine_one;
int i;
char *anchor;
int *anchor_score;
char *seq_anch;
int found;
int f_ind;
int ret_flag;
int block;

ret_flag = 0;

for ( ; ; )
{
block = 1;
#ifdef SLOW
mutex_lock(&read_file_lock);
if ((mine = read_file(block, scan_modifier)) == NULL) ret_flag = 1;
mutex_unlock(&read_file_lock);
#else
mutex_lock(&read_file_lock);
if ((mine = read_file_quick(block, scan_modifier, statbuf.st_size)) == NULL) ret_flag = 1;
mutex_unlock(&read_file_lock);
#endif
if (ret_flag == 1)
{
	return;
}
/* If translation is required then do it */
mine_one = mine;
if (TRANSLATION == TRUE)
{
	translate(mine, frames);
	/* Dont search the DNA, search the translation */
	if (mine->next != (rootres *)NULL) mine = mine->next;
}
/* We have to perform one scan with each feature */
while(mine != NULL)
{
clear_gaps(mine);
f_ind = 1;
while (f_ind < MAX_FEATURES && feature[f_ind] != NULL)
{
	found = 0;
	num_scans = mine->numres - feature[f_ind]->len_feature;
	anchor = mine->seq;
	anchor_score = mine->gap;
	for (i=0; i<=num_scans; i++)
	{
	/* Anchor to the first position of motif and check position in pattern */
		num_matches = 0;
       		seq_anch = anchor;
		current_residue = feature[f_ind]->start;
		while (current_residue != NULL)
		{
			current_down = current_residue->next;
			while (found == 0)
			{
				if (current_down == NULL ) break;
				if (current_down->res1 == *seq_anch)
				{
		  			*anchor_score += current_down->score_good;
					(current_down->score_good > gb_threshold) ? num_matches++ : num_matches;
					found = 1; /* Go on for the next position */
				}
				else
				{
					current_down = current_down->next;
				}
			}
			current_residue = current_residue->across;
			seq_anch++;
			found = 0;
		}
	*anchor_score *= num_matches; /* Special scan_novel modification */
	anchor++;
	anchor_score++; /* Move on to next segment of sequence */
	}
	/* Now work out the percentage scores */
	for (i=0; i<=num_scans; i++)
	{
		mine->fdata[i] = (100.0 * (float)mine->gap[i])/(float)feature[f_ind]->max_score_feature;
	}
	if (PLOT == TRUE) /* then write data to plotfile */
	{
		/* Write a title for this plot */
		fprintf(plot_file, "\n# %s\n", feature[f_ind]->header);
		for (i = 0; i <= num_scans; i++)
		{
			fprintf(plot_file, "%7.2f\n", mine->fdata[i]);
		}
		fprintf(plot_file,"\n");
	}
	else
	{
		mutex_lock(&merge_lock);
		merge_hits(mine, f_ind);
		mutex_unlock(&merge_lock);
	}
	for (i=0; i<=num_scans; i++)
	{
		mine->fdata[i] = 0;
		mine->gap[i] = 0;
	}
	f_ind++;
}/* End while */
mine = mine->next; /* Get the next sequence in the list */
} /* End while mine != NULL */
mutex_lock(&seq_count_lock);
seq_count += 1;
total_residues += residues(mine_one);
ccount += 1;
if (ccount >= count)
{
	printf("%d\n", seq_count);
	ccount = 0;
}
mutex_unlock(&seq_count_lock);
purge_all(mine_one);
mine = NULL;
} /* End for ever */
}/* End of scan_nsingle */



int max_hit(rootres *mres)
{
/* Return the seq array index of the next best hit */
register int i;
int max;
int index;

max=0;
index=0;
for (i=0; i<mres->numres; i++)
{
	if (mres->gap[i] > max) 
	{
		max = mres->gap[i];
		index = i;
	}
}
/*printf("finished max (%d) index (%d)\n",max,index);*/
return(index);
}/* End of max_hit */


int read_prints(int d_format, int scan_proc)
{
char *p;
int f_ind;
void read_prints_single(int);

f_ind = 0;

while (feof(in_data) == 0 )
{
	while (fgets(line, MAX_LINE, in_data) != NULL)
	{
		if (strncmp(line, "fc;", 3) == 0) break;
	}
	if (feof(in_data) != 0) break;
	if (++f_ind == MAX_FEATURES)
	{
		fprintf(stderr,"- sorry too many fingerprints.\n");
		fprintf(stderr,"- internal maximum %d\n", MAX_FEATURES);
		exit(1);
	}
	/* Make a new feature header structure and add to the list */
	feature[f_ind] = (fth *) malloc(sizeof(fth));
	demand(feature[f_ind], malloc failed in read_segments);
	strcpy(feature[f_ind]->file_spec, "PRINTS");

	/* Remove newline generated by fgets */

	p = line;
	if ((p = strchr(line, '\n')) != NULL) *p = '\0';
	line[80] = '\0'; /* Make certain string is terminated */
	strncpy(feature[f_ind]->code, &line[4], MAX_PRINTS_CODE_LEN);
	fgets(line, MAX_LINE, in_data);
	if (strncmp(line, "fl;", 3) != 0)
	{
		fprintf(stderr, "PRINTS format corrupt: fl specifier not found.\n");
		fprintf(stderr, "code: %s\n", feature[f_ind]->code);
		exit(1);
	}
	sscanf(&line[4], "%d", &feature[f_ind]->len_feature);
	fgets(line, MAX_LINE, in_data);
	if (strncmp(line, "ft;", 3) != 0)
	{
		fprintf(stderr, "PRINTS format corrupt: ft specifier not found.\n");
		fprintf(stderr, "code: %s\n", feature[f_ind]->code);
		exit(1);
	}
	p = line;
	if ((p = strchr(line, '\n')) != NULL) *p = '\0';
	strcpy(feature[f_ind]->header,&line[4]);
	feature[f_ind]->start = allocate_across(f_ind);
	/* f_ind tells allocate across which feature to use */
	if (d_format == SINGLE) read_prints_single(f_ind);
	else 
	{
	/*	read_position_pair(f_ind); */
		fprintf(stderr,"Method not yet implemented (pairwise).\n");
		exit(1);
	}
	max_score(f_ind, scan_proc);
	init_hitlist(f_ind, num_hits);
/*
	printf("Code: %s (%d)\n", feature[f_ind]->code, f_ind);
*/
}/* end while */
return(f_ind);
	
}/* End of read_prints */



void read_prints_single(int f_ind)
{
/* Reads the motif entry for each position which has a weight in the motif */
/* new structures of type link_down are created for each pair */
/* The next record in the file must contain the position number, starting at 0 */
int i;
int j;
int num_motifs;
int num_ignored;
short found_it;
link_residue *link;
link_down *current;
link_down *res;
link_down *list_end;
char *p;

num_motifs = 0;
num_ignored = 0;
while (fgets(line, LINE_LEN, in_data) != NULL)
{
/* bb signals the end of a section in prints */
if (strncmp(line, "bb;", 3) == 0) break;
link = feature[f_ind]->start; /* The addr of the first link */
p = &line[4]; /* Miss out the code portion of prints */
if (strspn(p,"TCAGNSPFLYHQVKDEIWRMBXZ-") != feature[f_ind]->len_feature)
{
	printf("%s %-.*s ignored\n",feature[f_ind]->code, feature[f_ind]->len_feature,p);
	num_ignored++;
	continue;
}
	num_motifs++;
	for (i=0; i<feature[f_ind]->len_feature; i++)
	{
		if (p[i] != '-')
		{
		for (j=i; j<= i; j++) /* Only one iteration necessary */
		{
			found_it = 0;
			if (p[j] == '-') continue;
			/* If there is already a node for this sep with these residues, then use it */
			if (link->next != NULL)
			{
				current = link->next;
				while (current != NULL)
				{
					if (current->sep != j-i)
					{
						current = current->next;
					}
					else
					{
						if ((current->res2 == p[j]) && (current->res1 == p[i]))
						{
						current->score_good++;
						found_it = 1;
						break;
						}
						else
						{
						current = current->next;
						}
					}
				}
			}
			if (found_it == 0)
			{
			res = (link_down *) malloc(sizeof(link_down));
			demand(res,malloc failed in read_position_single);
			/* DJPS 6 April 1994 */
			res->score_good = 0;
			res->next = NULL;
			if (link->pairs == 0)
			{
				link->next = res;
			}
			else
			{
			/* Hop to the end of the list */
			list_end = link->next;
			while (list_end->next != NULL) list_end = list_end->next;
			list_end->next = res;
			}
			res->res1 = p[i];
			res->res2 = p[j];
			res->sep = j-i;
			res->score_good++;
			link->pairs++;
			} /* End an if */
		}/* End for j */
		}/* End if*/
		link = link->across;
	}
/* now read the next motif */
} /* End while fgets */
/*
printf("\n%d motifs were read, %d accepted, %d ignored\n",num_motifs+num_ignored,
	num_motifs,num_ignored);
*/
return;
}/* read_prints_single */
