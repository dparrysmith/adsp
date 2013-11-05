/* Structure definitions for modules that link with y.c */

typedef struct link_residue link_residue;
typedef struct link_down link_down;

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
	int  sep;
	int  score_good;
	link_down *next;
};

extern link_residue *original_residue;
extern int len_motif;
extern int max_score_motif;
