/**
 * param.c
 * Description: Parameter read file for optimizers
 *              This file contains routines useful for parsing arguments for optimizers.
                set- and get-methods for file interface are included.
                This c-file is included by TOP-Optimizers.
		(TOP = Tool pour Optimization Parametrique)
 * Authors:     Michael Emmerich, Mutlu Oezdemir
 * Date:        2002
 */

#include "param.h"


int readParameters(int argc, char** argv) {
  int   found;
  int   foundInt, foundInt2;
  int   error;
  double foundDouble;
  char* foundString;
  int i, j, k, count, malloced;
  int reserved[argc];
  int boundsProvided = 1;
  double u;
  char *line, *ok;              
  FILE *datei;                 
  char **lines;  
  char name[MAX_FILENAME_LENGTH];
  muSET = false;
  lambdaSET = false;
  kappaSET = false;
  
  INITIALISED = 0;
  malloced    = 0;
  /*********************************************************************
   *                        NOTHING GIVEN !                            * 
   ********************************************************************/
  if(argc == 1) {
    printf("missing arguments\n");
    printf("Try `a.out --help' for more information.\n");
    exit(0);
  }

  /*********************************************************************
   *                           HELP !                                  * 
   ********************************************************************/
  if((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help")  == 0))
    {
      printf("Parameter configuration options for optimizers:\n\n");
      printf("-d, --dimension         dimension of the problem (MUST be provided !) \n");
      printf("-t, --target            target function to be optimized (MUST be provided !)\n");
      printf("-p, --problemName       name of the problem. (DEFAULT = %s)\n", DEFAULT_PROBLEM_NAME);
      printf("-m, --maxIt             maximum number of iterations for algorithm (DEFAULT = %i)\n", DEFAULT_ITERATION );
      printf("-f, --functValue        threshold function value (DEFAULT: No threshold function value)\n");
      printf("-V, --VERBOSE           verbose mode every iteration. 0 <=> print nothing, 1 <=> print iteration, f\n");
      printf("                        2 <=> print it, f, f_best, 3 <=> printf it, f, f_best, x-values\n");
      printf("                        4 <=> print iteration number at the end (DEFAULT = %i)\n", DEFAULT_VERBOSE_MODE);
      printf("-v, --verbose           verbose mode every x iterations, verbose mode first, x second parameter.\n");
      printf("                        (DEFAULT verbose iterations = %i)\n", DEFAULT_VERBOSE_ITERATION);
      printf("-r, --random            random seed (DEFAULT = %i)\n", DEFAULT_SEED);
      printf("-B, --BOUNDS            global bounds, upper and lower to be provided \n");
      printf("                        (DEFAULT LOWER BOUND= %2.2f, DEFAULT UPPER BOUND = %2.2f)\n", DEFAULT_LOWER, DEFAULT_UPPER);
      printf("-b, --bounds            individual bounds, to be provided for each variable, 2xDIM values \n");
      printf("-Y, --TYPES             global types to be provided. 0 = continuous, 1 = integer, 2 = nominal discrete\n");
      printf("                        (DEFAULT TYPE= continuous\n");
      printf("-y, --types             individual type, to be provided for each variable, DIM values \n");
      printf("-I, --INITGLOBAL        global initial point, must be within bounds. \n");
      printf("                        (DEFAULT for initial points is: Set randomly within bounds)\n");
      printf("-i, --initIndividual    individual initial points. DIM values to be provided within bounds\n");
      printf("-i, --initByFile\n");
      printf("-I, --INITBYFILE        name of file containing initial files\n");
      printf("-S, --STEPSGLOBAL       global step size for parameters. Number between 0.0 and 1.0, related to whole space.\n"); 
      printf("                        (DEFAULT for step sizes: %f\n", DEFAULT_STEP);
      printf("-s, --stepsIndividual   individual step sizes for parameters.\n");
      printf("-c, --configFile        name of the configuration file for strategy-specific parameters\n");
      printf("                        structure of a line of this config file should be:\n");
      printf("                        <name of parameter>   <value of parameter> \n");
      printf("                        EXAMPLE:    lambda    7\n");
      printf("                        (DEFAULT: default strategy parameters)\n");
      printf("--mu                    Size of the parent population\n");
      printf("--lambda                Size of the offspring population\n");
      printf("--kappa                 Maximum number of generations an individual can survive\n");
      printf("-h, --help              display this help and exit\n\n");
      exit(0);
    }
  

  /**
   * initiate reserved
   */
  for(i = 0; i < argc; i++) {
    reserved[i] = 0;
  }

  
    
  /*********************************************************************
   *                           DIMENSION                               * 
   ********************************************************************/
  foundInt = -1;
  for(i = 0; i < argc - 1; i++) {
    if(isNumber(argv[i+1]) == 0
       &&
       reserved[i] == 0     
       &&
       reserved[i+1] == 0     
       &&
       (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--dimension")  == 0))
      {
	foundInt = atoi(argv[i+1]);
	reserved[i]   = 1;
	reserved[i+1] = 1; 
	break;
      }    
  }
  
  if(foundInt < 1) {
    printf("You provided no dimension or dimension < 1 ! \n");
    printf("Dimension is provided by the -d or --dimension parameter. \n");
    printf("Example: -d 10\n");
    exit(0);
  }
  else
    dimension = foundInt;

  /*********************************************************************
   *                          TARGET FUNCTION                          *
   ********************************************************************/
  found = -1;
  for(i = 0; i < argc - 1; i++) {
    if(reserved[i] == 0     
       && 
       reserved[i+1] == 0
       &&
       (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--target") == 0))
      {
	found         = 1;
	foundString   = argv[i+1];
	reserved[i]   = 1;
	reserved[i+1] = 1; 
	break;
      }    
  }
  
  if(found == -1) {
    printf("You provided no target function !\n");
    printf("Target function is provided by the -t or --target parameter.\n");
    printf("Example: -t <target-function-file>\n");
    exit(0);
  }
  else
    target = foundString;


  /*********************************************************************
   *                     RANDOM NUMBER (optional)                      *
   ********************************************************************/
  found    = -1;
 
  for(i = 0; i < argc - 1; i++) {
    if(reserved[i] == 0     
       && 
       (reserved[i+1] == 0   && isNumber(argv[i+1]) == 0)
       &&
       (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--random") == 0))
      {
	found         = 1;
	foundInt      = atoi(argv[i+1]);
	reserved[i]   = 1;
	reserved[i+1] = 1; 
	break;
      }    
  }

  if(found != -1) 
    randomSeed = foundInt;
 else 
    randomSeed = 42;



  /*********************************************************************
   *                   THRESHOLD F-VALUE (optional)                    *
   ********************************************************************/
  found    = -1;
  thresholdSET = 0;
  for(i = 0; i < argc - 1; i++) {
    if(reserved[i] == 0     
       && 
       (reserved[i+1] == 0   && isNumber(argv[i+1]) == 0)
       &&
       (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--functValue") == 0))
      {
	found         = 1;
	foundDouble   = atof(argv[i+1]);
	reserved[i]   = 1;
	reserved[i+1] = 1; 
	break;
      }    
  }

  if(found != -1) {
    thresholdSET = 1;
    fThreshold   = foundDouble;    
  }
  //else 
    //randomSeed = 42;
  

  /*********************************************************************
   *                   ITERATION NUMBER (optional)                     *
   ********************************************************************/

  found        = -1;
  iterationSET =  0;
 
  for(i = 0; i < argc - 1; i++) {
    if(reserved[i] == 0     
       && 
       (reserved[i+1] == 0   && isNumber(argv[i+1]) == 0)
       &&
       (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--maxIt") == 0))
      {
	found         = 1;
	foundInt      = atoi(argv[i+1]);
	reserved[i]   = 1;
	reserved[i+1] = 1; 
	break;
      }    
  }
  
  if(found != -1 && foundInt > 0) {
    iterations = foundInt;
    iterationSET = 1;
  }
  else iterations = DEFAULT_ITERATION;
  

  /*********************************************************************
   *                   VERBOSE MODE (optional)                         *
   ********************************************************************/
  verboseMode = DEFAULT_VERBOSE_MODE;
  verboseIteration = DEFAULT_VERBOSE_ITERATION;
  
  found    = -1;
  
  for(i = 0; i < argc - 1; i++) {
    if(reserved[i] == 0     
       && 
       (reserved[i+1] == 0   && isNumber(argv[i+1]) == 0)
       &&
       (strcmp(argv[i], "-V") == 0 || strcmp(argv[i], "-VERBOSE") == 0))
      {
	found         = 1;
	foundInt      = atoi(argv[i+1]);
	reserved[i]   = 1;
	reserved[i+1] = 1; 
	break;
      }    
  }
  
  if(found != -1) {  
    if(foundInt >= 0 && foundInt < 5) {
      verboseMode      = foundInt;
      verboseIteration = 1; 
    }
    else {
      printf("Verbose Mode should be 0, 1 or 2! \n");
      printf("0: print f    1: print f, f_best  2: print f, f_best, x\n");
      printf("Example: -V 0\n");
      exit(0);
    }
  }  
  else {
    foundInt  = 0;
    foundInt2 = 0;
    for(i = 0; i < argc - 2; i++) {
      if(reserved[i] == 0     
	 && 
	 (reserved[i+1] == 0   && isNumber(argv[i+1]) == 0)
	 &&
	 (reserved[i+2] == 0   && isNumber(argv[i+2]) == 0)
	 &&
	 (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0))
	{
	  found         = 1;
	  foundInt      = atoi(argv[i+1]);
	  foundInt2     = atoi(argv[i+2]);
	  reserved[i]   = 1;
	  reserved[i+1] = 1; 
	  reserved[i+2] = 1;
	  break;
	}    
    }
    
    if(found == 1) {
      if(foundInt >= 0 && foundInt < 5) {
	if(foundInt2 > 0) {
	  verboseMode      = foundInt;
	  verboseIteration = foundInt2; 
	}
	else {
	  printf("Verbose Iteration number should be greater 0 ! \n");
	  printf("Example: -v 0 10\n");
	  printf("That means: Every ten iteration, f is printed out\n");
	  exit(0);
	}
      }
      else {
	printf("Verbose Mode should be 0, 1 or 2! \n");
	printf("0: print f  1: print f, f_best   2: print f, f_best, x\n");
	printf("Example: -v 0 10\n");
	exit(0);
      }
    }
    
  }
  

  /**********************************************************************
   *                PROBLEM NAME (optional)                             *
   *********************************************************************/
  found    = -1;
 
  for(i = 0; i < argc - 1; i++) {
    if(reserved[i] == 0     
       && 
       (reserved[i+1] == 0   && isNumber(argv[i+1]) != 0)
       &&
       (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--problemName") == 0))
      {
	found         = 1;
	foundString   = argv[i+1];
	reserved[i]   = 1;
	reserved[i+1] = 1; 
	break;
      }    
  }
  
  if(found != -1) {
    problemName = foundString;
  }
  else 
    problemName = DEFAULT_PROBLEM_NAME;

  sprintf(inFileName, "%s.in",  problemName);
  sprintf(resFileName,"%s.res", problemName);
 

  /**
   * initialise history file and tag-file
   */
  sprintf(name, "%s.history", problemName);  
  datei =  fopen(name, "w");
  fclose(datei);
  sprintf(name, "%s.tag", problemName);
  datei =  fopen(name, "w");
  fclose(datei);

  /**********************************************************************
   *                CONFIG FILE (optional)                              *
   *********************************************************************/
  found         = -1;
  configFileSET =  0;

  for(i = 0; i < argc - 1; i++) {
   
    if(reserved[i] == 0     
       && 
       (reserved[i+1] == 0)
       &&
       (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--configFile") == 0))
      {
	found         = 1;
	foundString   = argv[i+1];
	reserved[i]   = 1;
	reserved[i+1] = 1; 
	break;
      }    
  }
  
  if(found != -1) {  
    configFileName = foundString;
    configFileSET  = 1;  
    malloced       = 1;
    
    count = 0;
    lines = (char**)malloc(MAX_DIM*sizeof(charPtr));	
    datei =  fopen(configFileName, "r");
    if (!datei) {
      printf("Config file could not be found !\n");
      printf("Please provide a config file with following structure: \n");
      printf("<name>  <value>\n");
      printf("<name>  <value>\n");
      printf(" ...      ...  \n"); 
      printf("<name>  <value>\n");
      exit(0);
    }
    
    do 
      {
	line = (char*)malloc(sizeof(char)*LINE_LEN);
	ok = fgets(line, LINE_LEN-1, datei);
	if (ok) {lines[count] = line; count++;}
      } while (ok != NULL && count < MAX_DIM);
    
    fclose(datei);

    
    for(i = 0; i < count; i++) {
      sscanf(lines[i], "%s %s", configName[i], configValue[i]);
    }
  }
  configNumber = count;
  
  /*********************************************************************
   *                           MU (optional)                           * 
   ********************************************************************/
  foundInt = -1;
  for(i = 0; i < argc - 1; i++) {
    if(isNumber(argv[i+1]) == 0
       &&
       reserved[i] == 0     
       &&
       reserved[i+1] == 0     
       &&
       (strcmp(argv[i], "--mu")  == 0))
    {
      foundInt = atoi(argv[i+1]);
      reserved[i]   = 1;
      reserved[i+1] = 1; 
      break;
    }    
  }

  if(foundInt < 1) {
    printf("You provided no mu or mu < 1 ! \n");
    printf("Mu is provided by the --mu parameter. \n");
    printf("Example: -mu 10\n");
    exit(0);
  }
  else
  {
    mu_in = foundInt;
    muSET = true;
  }
  
  /*********************************************************************
   *                           LAMBDA (optional)                       * 
   ********************************************************************/
  foundInt = -1;
  for(i = 0; i < argc - 1; i++) {
    if(isNumber(argv[i+1]) == 0
       &&
       reserved[i] == 0     
       &&
       reserved[i+1] == 0     
       &&
       (strcmp(argv[i], "--lambda")  == 0))
    {
      foundInt = atoi(argv[i+1]);
      reserved[i]   = 1;
      reserved[i+1] = 1; 
      break;
    }    
  }

  if(foundInt < 1) {
    printf("You provided no lambda or lambda < 1 ! \n");
    printf("Lambda is provided by the --lambda parameter. \n");
    printf("Example: -lambda 10\n");
    exit(0);
  }
  else
  {
    lambda_in = foundInt;
    lambdaSET = true;
  }
   
  /*********************************************************************
   *                           KAPPA (optional)                        * 
   ********************************************************************/
  foundInt = -1;
  for(i = 0; i < argc - 1; i++) {
    if(isNumber(argv[i+1]) == 0
       &&
       reserved[i] == 0     
       &&
       reserved[i+1] == 0     
       &&
       (strcmp(argv[i], "--kappa")  == 0))
    {
      foundInt = atoi(argv[i+1]);
      reserved[i]   = 1;
      reserved[i+1] = 1; 
      break;
    }    
  }

  if(foundInt < 1) {
    printf("You provided no kappa or kappa < 1 ! \n");
    printf("Kappa is provided by the --kappa parameter. \n");
    printf("Example: -kappa 10\n");
    exit(0);
  }
  else
  {
    kappa_in = foundInt;
    kappaSET = true;
  }

  /***********************************************************************
   *                           BOUNDS (optional)                         *
   ***********************************************************************/
  
  found    = -1;
  
  for(i = 0; i < argc - 2; i++) {
    if(reserved[i] == 0     
       && 
       (reserved[i+1] == 0   && isNumber(argv[i+1]) == 0)
       &&
       (reserved[i+2] == 0   && isNumber(argv[i+2]) == 0)
       &&
       (strcmp(argv[i], "-B") == 0 || strcmp(argv[i], "--BOUNDS") == 0))
      {
	found         = 1;
	reserved[i]   = 1;
	reserved[i+1] = 1; 
	reserved[i+2] = 1;
	lower[0]      = atof(argv[i+1]);
	upper[0]      = atof(argv[i+2]);
	break;
      }    
  }
  
  if(found == 1) {
    if(lower[0] >= upper[0]) {
      printf("You provided an upper bound ");
      printf("<= lower bound ! \n");
      exit(0);
    }
    for(i = 1; i < dimension; i++) {
      upper[i] = upper[0];
      lower[i] = lower[0];
    }    
  }
  
  else {
    error = 0;
    for(i = 0; i < argc - 1; i++) {
      if(reserved[i] == 0     
	 &&
	 (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--bounds") == 0))
	{
	  found = 1;
	  for(j = i + 1; j < i + 1 + 2 * dimension; j += 2) {
       
	    if(reserved[j] == 0     && reserved[j+1] == 0
	       && isNumber(argv[j])==0 && isNumber(argv[j+1])==0
	       && atof(argv[j+1])    > atof(argv[j]))
	      {
		lower[(j-i-1)/2] = atof(argv[j]);
		upper[(j-i-1)/2] = atof(argv[j+1]);
	      }
	    else {
	      error = 1;
	      break;
	    }
	    
	  }
	}    
    }
  }
  if(error == 1)
    {
      printf("You did not set the individual bounds correctly \n");
      printf("or their number is not enough. Please check !\n");
      printf("Example: Dimension is 2. Then individual parameters are ");
      printf("set by: \n -b <lower1> <upper1> <lower2> <upper2>.\n");
      printf("Consider that upper > lower !!\n");
      exit(0);
    }
  
  if(found == -1) {
    /**
     * no information about bounds found..
     * so we set them loosely..
     */
    for(i = 0; i < dimension; i++) {
      lower[i] = DEFAULT_LOWER;
      upper[i] = DEFAULT_UPPER;
    }    
    boundsProvided = -1;
  }
  
  



  /***********************************************************************
   *                           TYPES (optional)                         *
   ***********************************************************************/
  
  found    = -1;
  int temp_type = -1;

  for(i = 0; i < argc - 1; i++)
  {
    if(reserved[i] == 0     
        && 
        (reserved[i+1] == 0 && isNumber(argv[i+1]) == 0)
        &&
        (strcmp(argv[i], "-Y") == 0 || strcmp(argv[i], "--TYPES") == 0))
    {
      found         = 1;
      reserved[i]   = 1;
      reserved[i+1] = 1;
      temp_type     = atoi(argv[i+1]);
      if(valid_type(temp_type))
      {
        types[0]    = static_cast<TYPE>(temp_type);
      }
      break;
    }    
  }

  if(found == 1)
  {
    for(i = 1; i < dimension; i++)
    {
      types[i] = types[0];
    }    
  }
  else
  {
    error = 0;
    for(i = 0; i < argc; i++)
    {
      if(reserved[i] == 0     
          &&
          (strcmp(argv[i], "-y") == 0 || strcmp(argv[i], "--types") == 0))
      {
        if(argc < i + dimension)
        {
          error = 1;
          break;
        }
        found = 1;
        for(j = i + 1; j < i + 1 + dimension; ++j)
        {
          if(reserved[j] == 0 && isNumber(argv[j]) == 0)
          {
            temp_type = atoi(argv[j]);
            if(valid_type(temp_type))
            {
              types[j - i - 1] = static_cast<TYPE>(temp_type);
            }
          }
          else
          {
            error = 1;
            break;
          }

        }
      }
    }
  }

  if(error == 1)
  {
    printf("You did not set the individual types correctly \n");
    printf("or their number is not enough. Please check !\n");
    printf("Example: Dimension is 2. Then individual parameters are ");
    printf("set by: \n -y <type1> <type2>.\n");
    printf("Valid types are: 0 = continuous, 1 = integer, 2 = nominal discrete\n");
    exit(0);
  }

  if(found == -1)
  {
    /**
     * no information about types found..
     * so we set them all to continuous..
     */
    for(i = 0; i < dimension; i++)
    {
      types[i] = TYPE::DEFAULT;
    }    
  }

  // Print set types
//  for(i = 0; i < dimension; ++i)
//  {
//    printf("types[%d]: %d\n", i, static_cast<int>(types[i]));
//  }

  /************************************************************************
   *                  INITIAL POINTS (optional)                           *
   ***********************************************************************/
  found    = -1;
  
  for(i = 0; i < argc - 1; i++) {
    if(reserved[i] == 0     
       && 
       (reserved[i+1] == 0   && isNumber(argv[i+1]) == 0)
       &&
       (strcmp(argv[i], "-I") == 0 || strcmp(argv[i], "--INITGLOBAL") == 0))
      {
	found         = 1;
	reserved[i]   = 1;
	reserved[i+1] = 1; 
	xInit[0]      = atof(argv[i+1]);
	break;
      }    
  }
  
  if(found == 1) {
    for(i = 1; i < dimension; i++) {
      xInit[i] = xInit[0];
    }    
  }
  
  else {
    error = 0;
    for(i = 0; i < argc - 1; i++) {
      if(reserved[i] == 0     
	 &&
	 (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--initIndividual") == 0))
	{
	  if(isNumber(argv[i+1]) == 0 
	     &&
	     argc < i + 1 + dimension)
	    {
	      printf("You did not provide enough initial design sites !!\n");
	      exit(0);
	    }
	    
	  if(isNumber(argv[i+1]) == 0) {
	    found = 1;
	    
	    for(j = i + 1; j < i + 1 + dimension; j += 1) {
	      
	      if(reserved[j] == 0 && isNumber(argv[j]) == 0)
		{
		  xInit[j - i - 1] = atof(argv[j]);
		}
	      else {
		error = 1;
		break;
	      }
	      
	    }
	  }    
	}
    }
  }
  if(error == 1)
    {
      printf("You did not set the initial points correctly \n");
      printf("or their number is not enough. Please check !\n");
      printf("Example: Dimension is 2. Then initial points maybe ");
      printf("set by: \n -i x1 x2\n");
      printf("Consider x1 and x2 should be within bounds if provided !!\n");
      exit(0);
    }

  
  /**
   * are the initial files given by a file ?
   */
  error = 0;
  if(found == -1) {
    for(i = 0; i < argc - 1; i++) {
    if(reserved[i] == 0     
       && 
       reserved[i+1] == 0   && !(isNumber(argv[i+1]) == 0)
       &&
       (strcmp(argv[i], "-I") == 0 || strcmp(argv[i], "--INITBYFILE") == 0 ||  
	strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--initByFile") == 0))
      {

	if(PRINTOUT) printf("READING x-Values from file !!\n");
	/* read x values from file *********************** */
	count    = 0;
	malloced = 1;
	lines = (char**)malloc(MAX_DIM*sizeof(charPtr));	
	datei =  fopen(argv[i+1], "r");
	if (!datei) {
	  printf("Initial file for initial design sites could not be found !\n");
	  exit(0);
	}
	
	do 
	  {
	    line = (char*)malloc(sizeof(char)*LINE_LEN);
	    ok = fgets(line, LINE_LEN-1, datei);
	    if (ok) {lines[count] = line; count++;}
	  } while (ok != NULL && count < MAX_DIM);
	
	fclose(datei);
	if(count != dimension+1) {
	  printf("The number of variables in in-File is not correct !\n");
	  exit(0);
	}

	
	for(k = 1; k < count; k++) {
	  xInit[k-1] = atof(lines[k]);
	}

	found         = 1;
	reserved[i]   = 1;
	reserved[i+1] = 1; 
	break;
      }    
    }

  }
  

  if(found == -1) {
    /**
     * No information about initial points found.
     * Setting them randomly within bounds.
     */
    for(i = 0; i < dimension; i++) {
      u = ((double)(rand() % 10000))/10000.0;
      if(boundsProvided == -1) 
	xInit[i] = DEFAULT_LOWER + (DEFAULT_UPPER - DEFAULT_LOWER) * u;
      else {
	xInit[i] = lower[i] + (upper[i] - lower[i]) * u;
      }

    }
  }

  /**
   * Check whether given initial point is within bounds
   * if bounds were provided
   */
  for(i = 0; i < dimension; i++) {
    if(!(xInit[i] <= upper[i] && xInit[i] >= lower[i])) {
      printf("Given initial point is not within bounds !!\n");
      printf("Please take the bounds into consideration !\n");
      printf("x[%i] = %f\n",i,xInit[i]);
      exit(0);
    }
  }
  


  /************************************************************************
   *                      STEP LENGTH (optional)                          *
   ***********************************************************************/
  
  found    = -1;
  
  for(i = 0; i < argc - 1; i++) {
    if(reserved[i] == 0     
       && 
       (reserved[i+1] == 0   && isNumber(argv[i+1]) == 0)
       &&
       (strcmp(argv[i], "-S") == 0 || strcmp(argv[i], "--STEPSGLOBAL") == 0))
      {
	found         = 1;
	reserved[i]   = 1;
	reserved[i+1] = 1; 
	if(atof(argv[i+1]) > 0.0 && atof(argv[i+1]) <= 1.0)
	  steps[0]      = atof(argv[i+1]);
	else {
	  printf("You provided a step length out of the valid range [0.0, 1.0]\n");
	  printf("This is not allowed, please check !\n");
	  exit(0);
	}
	break;
      }    
  }
  
  if(found == 1) {
    for(i = 1; i < dimension; i++) {
      steps[i] = steps[0];
    }    
  }
  
  else {
    error = 0;
    for(i = 0; i < argc - 1; i++) {
      if(reserved[i] == 0     
	 &&
	 (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--stepsIndividual") == 0))
	{
	  if(isNumber(argv[i+1]) == 0 
	     &&
	     argc < i + 1 + dimension)
	    {
	      printf("You did not provide enough initial step lengths !!\n");
	      exit(0);
	    }
	    
	  found = 1;
	    
	    for(j = i + 1; j < i + 1 + dimension; j += 1) {
	      
	      if(reserved[j] == 0 && isNumber(argv[j]) == 0)
		{
		  if(atof(argv[j]) > 0.0 && atof(argv[j]) <= 1.0) 
		    steps[j - i - 1] = atof(argv[j]);
		  else {
		    printf("You provided a step length out of the valid range [0.0, 1.0]\n");
		    printf("This is not allowed, please check !\n");
		    exit(0);
		  }
		}
	      else {
		error = 1;
		break;
	      }
	      
	    }
	}
    }
  }
  if(error == 1)
    {
      printf("You did not set the init. indiv. step length correctly !! \n");
      exit(0);
    }
  
  if(found == -1) {
    /**
     * set the step lengths randomly as they are not provided by user ...
     */
    for(i = 0; i < dimension; i++) {
      u = ((double)(rand() % 10000))/10000.0;
      steps[i] = u;
    }
  }
  if(PRINTOUT) {
    printf("dimension        = %i\n", dimension);
    printf("target           = %s\n", target);
    printf("problemName      = %s\n", problemName);
    printf("iterations       = %i\n", iterations);
    printf("random seed      = %i\n", randomSeed);
    printf("verboseIteration = %i\n", verboseIteration);
    printf("verboseMode      = %i\n", verboseMode);
    printf("configFileSet    = %i\n", configFileSET);
        
    for(i = 0; i < dimension; i++) {
      printf("x[%i]=%f   lb[%i]=%f   ub[%i]=%f s[%i] = %f\n", i, xInit[i], i, lower[i], i, upper[i], i, steps[i]);
    }
  }

  INITIALISED = 1;

  /**
   * some garbage collection ...
   */
  if(malloced) {
    //for(i = 0; i < MAX_DIM; i++) 
    //  free(lines[i]);
    free(lines);
  }

  return 1;
}
    
// Return true if the type is valid, exit the program with a warning if it is
// not
bool valid_type(int type)
{
  if(type < static_cast<int>(TYPE::MIN) || type > static_cast<int>(TYPE::MAX))
  {
    printf("You provided an invalid type. Options are: 0 = continuous, ");
    printf("1 = integer, 2 = nominal discrete \n");
    exit(0);
  }
  else
  {
    return true;
  }
}

/**
 * get the parent population size
 */
int getMu() {
  if(muSET)
    return mu_in;
  else
  {
    return 0;
  }
}


/**
 * get the offspring population size
 */
int getLambda() {
  if(lambdaSET)
    return lambda_in;
  else
  {
    return 0;
  }
}


/**
 * get the maximal lifetime of an individual
 */
int getKappa() {
  if(kappaSET)
    return kappa_in;
  else
  {
    return 0;
  }
}

/**
 * get the dimension of the problem 
 */
int getDimension() {
  if(INITIALISED)
    return dimension;
  else
    printFatalMessage();
}

/**
 * get the target function
 */
char* getTarget() {
  if(INITIALISED)
    return target;
  else
    printFatalMessage();
}

/**
 * get name of the problem
 */
const char* getProblemName() {
  if(INITIALISED)
    return problemName;
  else
    printFatalMessage();
}

/**
 * get iteration number
 */
int getIterations() {
  if(INITIALISED)
    return iterations;
  else
    printFatalMessage();
}

/**
 * get verbose mode
 */
int getVerboseMode() {
  if(INITIALISED)
    return verboseMode;
  else
    printFatalMessage();
}

/**
 * get iteration number for output
 */
int getVerboseIterations() {
  if(INITIALISED)
    return verboseIteration;
  else
    printFatalMessage();
}

/**
 * check if config file name is provided by user
 */
int getConfigFileSET() {
  if(INITIALISED)
    return configFileSET;
  else
    printFatalMessage();
}


/**
 * get i'th entry of upper bound
 */
double getUpperBound(int index) {
  if(INITIALISED) {
    if(index >= dimension) {
      printf("given index to access upper bounds >= dimension ! \n");
      printf("please consider !\n");
      exit(0);
    }
    else {
      return(upper[index]);
    }
  }
  else
    printFatalMessage();
  
}

/**
 * get i'th entry of upper bound
 */
double getLowerBound(int index) {
  if(INITIALISED) {
    if(index >= dimension) {
      printf("given index to access lower bounds >= dimension ! \n");
      printf("please consider !\n");
      exit(0);
    }
    else {
      return(lower[index]);
    }
  }
  else
    printFatalMessage();
}

/**
 * get i'th entry of type
 */
TYPE getType(int index) {
  if(INITIALISED) {
    if(index >= dimension) {
      printf("given index to access type >= dimension ! \n");
      printf("please consider !\n");
      exit(0);
    }
    else {
      return(types[index]);
    }
  }
  else
    printFatalMessage();
}

/**
 * get i'th entry of initial design site
 */
double getXInit(int index) {
  if(INITIALISED) {
    if(index >= dimension) {
      printf("given index to access initial design site >= dimension ! \n");
      printf("please consider !\n");
      exit(0);
    }
    else {
      return(xInit[index]);
    }
  }
  else
    printFatalMessage();
}


/**
 * get i'th entry of steps, absolute value, not relative !
 */
double getStep(int index) {
 if(INITIALISED) {
    if(index >= dimension) {
      printf("given index to access initial step sizes >= dimension ! \n");
      printf("please consider !\n");
      exit(0);
    }
    else {
      return(steps[index] * (getUpperBound(index) - getLowerBound(index)));
    }
  }
  else
    printFatalMessage();

}

/**
 * get random seed number
 */
int getRandomSeed() {
 if(INITIALISED)
    return randomSeed;
  else
    printFatalMessage();
}


/**
 * Get value of parameter with given name
 */
const char* getParameter(const char* name) {
  int i;
 
  if(INITIALISED) {
    if(configFileSET) {
      for(i = 0; i < configNumber; i++) {
	if(strcmp(configName[i], name) == 0) 
	  return configValue[i];
      }
      return "empty";
    }
    else {
      printf("Config file is not given ! \n");
      printf("Please first check if it is given, then use this method !!");
      exit(0);      
    }
        
  }
  else
    printFatalMessage();
}


/**
 * get benchmark home directory
 */
char* getBenchmarkHome() {
  if(getenv("BENCH_HOME") != NULL && (strcmp(getenv("BENCH_HOME"),"") != 0))
    return getenv("BENCH_HOME");
  else {
    printf("Environment variable BENCH_HOME not set or is empty !\n");
    printf("Please set this variable, this is absolutely necessary !\n");
    printf("Try to set it by 'export BENCH_HOME=...' or 'set BENCH_HOME=..'\n");
    printf("The set command depends on the shell You use.\n");
    exit(0);
  }
  
}

/**
 * 
 */
void printFatalMessage() {
  printf("Please first initialise getParam-method with arguments !\n");
  printf("After that parameters will be available !\n");
  printf("HAVE A NICE DAY ! \n");
  exit(0);
}

/**
 * print out depending of VERBOSE_MODE
 */
void print(int iteration, double fValue, double fBest, double *x) {
  int p;
  FILE *datei;
  char name[MAX_FILENAME_LENGTH];

  sprintf(name, "%s.history", getProblemName());
 
  
  datei =  fopen(name, "a");
  fprintf(datei, "%i  %f", iteration, fValue);
  
  fprintf(datei," %f", fBest);
  for(p = 0; p < getDimension(); p++) { 
    fprintf(datei, " %f", x[p]);
  }
  fprintf(datei,"\n");
  fclose(datei);

  if(getVerboseMode() != 4 && getVerboseMode() > 0) 
    {
      if(((iteration + 1) % getVerboseIterations()) == 0
	 ||
	 iteration == 0
	 ||
	 checkTermination(iteration, fBest)) 
	{
	  printf("%i  %f", iteration, fValue);
	  if(getVerboseMode() > 1) 
	    {
	      printf(" %f", fBest);
	      if(getVerboseMode() == 3) 
		for(p = 0; p < getDimension(); p++) { 
		  printf(" %f", x[p]);
		}
	    }
	  printf("\n");
	  
	}
    }
  else {
    if(checkTermination(iteration, fBest) == 1 && getVerboseMode() != 0)
      {
	printf("%i\n", iteration);
      }

  }

}

/**
 * Return 1 if termination criterion is fulfilled, 
 * 0 else
 * remove tag-file before returning 1
 */
int checkTermination(int iteration, double fBest) {
  FILE* datei;
  char cmd[1000];
  char name[MAX_FILENAME_LENGTH];

  /**
   * exit, if there is no tag-file !!
   */
  sprintf(name, "%s.tag", getProblemName());
  datei =  fopen(name, "r");
  if(!datei)
    exit(0);
  else
    fclose(datei);

  
  if(iterationSET == 0 && thresholdSET == 1) 
    {
      if(fBest <= fThreshold) {
	sprintf(cmd,"rm %s.tag", getProblemName());
	system(cmd);	
	return 1;
      }
      else
	return 0;
    }      
  
  if(thresholdSET == 0) {
    if(iteration >= (getIterations() - 1)) {
      sprintf(cmd,"rm %s.tag", getProblemName());
      system(cmd);
      return 1;
    }
    else
      return 0;
  }

  if(thresholdSET == 1 && iterationSET == 1) {
    if((iteration >= (getIterations() - 1)) || (fBest <= fThreshold)) {
      sprintf(cmd,"rm %s.tag", getProblemName());
      system(cmd);
      return 1;
    }
    else
      return 0;
  }
}


/**
 * returns 0 if the string argument is a number
 * else 1
 */
int isNumber(char* string) {
  int intValue  = atoi(string);
  
  int number = 1;
  
  if(intValue != 0 
     || 
     (intValue == 0 && 
      ((strncmp(string, "0", 1) == 0) || (strncmp(string, "-0", 2) == 0))))
    number = 0;

  return number;
}


/**
 * Method needed by evolutionary strategies
 */
double gauss() {
  double nv[11];  /* Gauss-distribution table */
  int flag = 0;
  double t;      /* Approximately normal distributed value */
  double u;      /* Uniform Distributed value */          
 
  int ub,lb;

  nv[0]=0.00;          /* 11 Values of Gauss-Curve */
  nv[1]=0.14; nv[2]=0.26; nv[3]=0.39; nv[4]=0.53;
  nv[5]=0.68; nv[6]=0.86; nv[7]=1.04; nv[8]=1.30;
  nv[9]=1.65; nv[10]=4.0;
      
  u =((double)(rand() % 10000))/10000.0;
  
  flag = 0;
  
  if (u < 0.5) { u = 1.0 - u; flag = 1;}
  u = (u-0.5)*20.0;
  lb = (int)u;
  ub = lb + 1;
  
  t = nv[lb] + (u - (double)lb) / ((double)ub -(double)lb) 
               * (nv[ub] - nv[lb]);
  if(flag == 1) return -t; else return t;
  
}


/**
 * 8 Methods(4 extended, 4 simple) for file interface following ...
 *
 */  

double getResultX() {
  return(getResult(getTarget(), inFileName, resFileName));
}

int setResult(char* outputFileName, double result) {
  FILE *datei;
  datei =  fopen(outputFileName, "w");
  if (!datei) {printf("FEHLER: output-File nicht gefunden\n");exit(0);}

  fprintf(datei, "%f\n", result); 
  fclose(datei);
  return 0;
}

double getResult(char* function, char* inputFileName,char* outputFileName) {
  char cmd[1000];
  double f;
  FILE *datei;                
  char *line;

  sprintf(cmd,"%s %s %s",function,inputFileName,outputFileName);
  system(cmd);
  datei =  fopen(outputFileName, "r");
  line = (char*)malloc(sizeof(char)*LINE_LEN);
  fgets(line, LINE_LEN-1, datei);
  f = atof(line);
  fclose(datei);
  return f;
}

int setXvector(char* inputFileName, double* xE, int dim) {
  FILE *datei;
  int i;

  datei =  fopen(inputFileName, "w");
  if (!datei) {printf("FEHLER: 'in'-File nicht gefunden\n");exit(0);}

  fprintf(datei, "%i\n", dim);
  for(i = 0; i < dim; i++)
    fprintf(datei,"%f\n", xE[i]);
  fclose(datei);

  return 0;
}

int getXvector(char* inputFileName, double* xE) {
  int i,count = 0, dimension;
  char ** lines;
  FILE *datei;
  char *line, *ok; 
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
  for(i = 1; i < count; i++) {
    xE[i-1] = atof(lines[i]);
  }
  dimension = atoi(lines[0]);

  for(i = 0; i < MAX_DIM; i++) 
    free(lines[i]);
  free(lines);
  return dimension;  
}

