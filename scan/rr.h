#define THREE_FRAME 0
#define SIX_FRAME 1
#define SCAN_ALL 0
#define SCAN_P1 1
#define SCAN_FA 2
#define DNA 0
#define PROTEIN 1

typedef struct pres pres;

struct pres
{
	pres *next;
	char *assoc;
};

typedef struct rootres rootres;

struct rootres
{
	rootres *next;      /* Next rootres structure in list */
	rootres *previous;  /* Previous rootres structure in list */
	rootres *link_s;    /* Next link partner (NULL if not set */
	rootres *link_p;    /* Previous link partner */
	char *seq;          /* Pointer to  the sequence array */
	pres *assoc;        /* Pointer to first block of per residue info */
	char *code;         /* Pointer to the sequence code */
	char *comment;      /* Pointer to the sequence comment */
	int *gap;           /* Pointer to the gap array (per residue) */
	float *fdata;       /* Pointer to the floating point data array */
	int type;           /* DNA or PROTEIN */
	int frame;          /* -3, -2, -1, 1, 2 or 3 */
	int bases;          /* Used for calulating reverse frame stats */
	int numres;         /* Total number of residues in this sequence */
	int numgap;         /* Total number of gaps in this sequence */
	int marker[2];      /* Marker[0] number of gaps before the first residue */
                            /* Marker[1] the number of the first residue */
};

extern rootres *first;
FILE *data;
