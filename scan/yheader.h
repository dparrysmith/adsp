/*
Title       : Yheader
Module      : Yheader.h
Author      : D.J.Parry-Smith
Date        : 5-JUNE-1989
History     :
Description : Header information for Y
*/

#define		YMVERSION	"V1.0"
#define		demand(fact, remark) {\
		if (!(fact)) {\
		fprintf(stderr,"demand not met: fact\n");\
		fprintf(stderr,"remark\n");\
		exit(0);\
		}\
		}
typedef		int boolean;
#define		FALSE		0
#define		TRUE		!FALSE
#define		LINE_LEN	512

