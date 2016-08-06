 
/* mystr.h for ANSI C */

#ifndef MYSTR_H
#define MYSTR_H
 

void mystr_Delete (char S[], size_t index, size_t len);



void mystr_Insert (char Res[], char Source[], size_t Pos);



void mystr_ItemS (char R[], char S[], const char T[], unsigned int N);



int mystr_Match (char Source[], char Pattern[]);



void mystr_Slice (char R[], char S[], size_t P, size_t L);



void mystr_Subst (char Source[], char OldPattern[], char NewPattern[]);



void mystr_Position (char Substring[], char Source[], size_t at,
	size_t * pos, size_t * found);

 

#endif
 
