#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define MAX_DIM 100
#define LINE_LEN 500
#define HIR(A) cout << "HIER" << (A) << endl;
#define MAX_DOUBLE 10000000000.0
#define FORALL(A,B) for((A)=0; (A) < (B);A++)
#define MIN(A,B)  (((A) < (B)) ?  (A) : (B))
typedef char* charPtr;

// functions
double sphere(int,double*);


int main(int argc, char **argv)
{
  char *line, *ok;              
  FILE *datei;                 // File-Descriptor
  char **lines;                // File-Buffer    
  int count = 0;               // Line-Counter for input-request
  double x[MAX_DIM];           // Objaect variables read from file
  double res;                  // Result of Function
  char *inputFileName;         // in-File
  char *outputFileName;        // res-File
  
  if (argc != 3) 
    {
      printf("function <in-File> <res-File>\n"); exit(0);
    }  
  inputFileName = argv[1];      
  outputFileName = argv[2];
  

  // read x values from file ***********************
  lines = (char**)malloc(MAX_DIM*sizeof(charPtr));	
  datei =  fopen(inputFileName, "r");
  if (!datei) {printf("FEHLER: 'in'-File nicht gefunden\n");exit(0);}
  
  do 
    {
      line = (char*)malloc(sizeof(char)*LINE_LEN);
      ok = fgets(line, LINE_LEN-1, datei);
      if (ok) {lines[count] = line; count++;}
    } while (ok != NULL && count < MAX_DIM);
  
  fclose(datei);
  for(int i = 1; i < count; i++) {
    x[i-1] = atof(lines[i]);
    
  }
  
  //******************* function call *********************
  res = sphere(atoi(lines[0]), x);  // dimension length = count - 1
  //*******************************************************

  // write back result ****************************
  datei =  fopen(outputFileName, "w");
  fprintf(datei,"%.38f\n",res);
  fclose(datei);
}

double sphere(int n,double *x)
{
  // Problem 2.17 of Schwefel
	int i;
	double	res = 0.0;

	for (i=0;i<n;i++){
	  res += x[i]*x[i];
  }
  return(res);
}


