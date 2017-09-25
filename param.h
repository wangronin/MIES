#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <time.h>
#include <math.h>

#ifndef _PARAM_TIP
# define _PARAM_TIP

#define LINE_LEN 500
#define MAX_FILENAME_LENGTH 100
#define MAX_DIM 200
#define MAX_LAMBDA 100
#define MAX_DOUBLE 10000000000.0
#define FORALL(A,B) for((A)=0; (A) < (B);A++)
#define PRINTOUT 0

#define DEFAULT_LOWER      -10.0
#define DEFAULT_UPPER       10.0
#define DEFAULT_ITERATION  1000
#define DEFAULT_SEED        42
#define DEFAULT_PROBLEM_NAME "a"
#define DEFAULT_VERBOSE_MODE 2
#define DEFAULT_VERBOSE_ITERATION 1
#define DEFAULT_STEP 0.1
#define MAX_PARALLEL_NUMBER 20
typedef char* charPtr;

enum class TYPE
{
  CONTINUOUS = 0,
  INTEGER = 1,
  NOMINAL_DISCRETE = 2,

  MIN = CONTINUOUS,
  MAX = NOMINAL_DISCRETE,

  DEFAULT = CONTINUOUS
};

/**
 * methods
 */
bool   valid_type(int type);
void   printFatalMessage(void);
void   print(int, double, double, double*);
char*  getTarget(void);
const char*  getProblemName(void);
char*  getParameter(char*);
char*  getBenchmarkHome();
int    isNumber(char*);
int    readParameters(int, char**);
int    checkTermination(int, double);
int    getDimension(void);
int    getIterations(void);
int    getVerboseMode(void);
int    getVerboseIterations(void);
int    getConfigFileSET(void);
int    getRandomSeed();
double getUpperBound(int);
double getLowerBound(int);
TYPE   getType(int);
double getXInit(int);
double getStep(int);
double gauss(void);
double getResult(char *x1, char *x2, char *x3);
int    setResult(char *x1, double f);
int    getXvector(char *x1,double* x);
int    setXvector(char* x1, double* x2, int x3);
double getResultX();
int    getParallelNumber(void);
int    setXvectors(double** x1, char* x2, int x3, int x4);
double* getResults(char* x1, char* x2, char* x3, int x4);
int    getMu();
int    getLambda();
int    getKappa();

/*
 * global optimization parameters
 */
int             dimension;  /* dimension of the problem                           */
int            randomSeed;  /* random seed number                                 */
int            iterations;  /* number of iterations algorithms should make        */
int      verboseIteration;  /* output every verboseIteration                      */
int           verboseMode;  /* verbose mode (0 <=> only f / 1 <=> f, f_best, x's) */
int         configFileSET;  /* indicator. shows if config file name is provided   */
int          thresholdSET;  /* indicator whether threshold f-value is set         */
int          iterationSET;  /* indicator whether iteration number is provided     */
bool                muSET;
bool            lambdaSET;
bool             kappaSET;
int                 mu_in;
int             lambda_in;
int              kappa_in;
char*              target;  /* target function to be optimised                    */
const char*   problemName;  /* name of the problem                                */
char*      configFileName;  /* name of the config file for strategy parameters    */
char           inFileName[MAX_FILENAME_LENGTH];
char          resFileName[MAX_FILENAME_LENGTH];  
int          configNumber;  /* number of config parameters given in file          */
double     xInit[MAX_DIM * MAX_PARALLEL_NUMBER];  /* initial design site          */
double     steps[MAX_DIM * MAX_PARALLEL_NUMBER];  /* step length                  */
double     upper[MAX_DIM * MAX_PARALLEL_NUMBER];  /* upper bounds                 */ 
double     lower[MAX_DIM * MAX_PARALLEL_NUMBER];  /* lower bounds                 */
TYPE       types[MAX_DIM * MAX_PARALLEL_NUMBER];  /* variable types               */
double         fThreshold;  /* threshold value for target function value          */
char  configName[MAX_DIM][LINE_LEN];  /*  name of config parameter                */
char configValue[MAX_DIM][LINE_LEN];  /*  value of config parameter               */
int parallelNumber;

/**
 * Parameter showing whether arguments are already read and initialised ?
 */
int INITIALISED;

#endif
