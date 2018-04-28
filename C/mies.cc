/******************************************************/
/*                Mkl strategy                        */
/*      -  with Modified Interval-Bounds Treatment    */ 
/*    (c) by Michael Emmerich, Jing Zhou              */
/******************************************************/

#include "param.cc"

#define DEFAULT_LAMBDA 28
#define DEFAULT_MU 4
// For both RECO and RECO_S:
// 0 = always take parent 1
// 1 = intermediate recombination (take the arithmetic mean of parents 1 and 2)
// Greater than 1 = dominant recombination (select parent uniformly at random)
#define DEFAULT_RECO 2
#define DEFAULT_RECO_S 1
#define DEFAULT_KAPPA 1

enum class STEP_SIZE_MODE
{
  SINGLE,
  MULTI
};

double getStepMIES(TYPE type, double lb, double ub, int dim_d);
double uniform();
int get_dim_r(const TYPE (&type)[MAX_DIM]);
int get_dim_i(const TYPE (&type)[MAX_DIM]);
int get_dim_d(const TYPE (&type)[MAX_DIM]);
double get_tau_global(int dim);
double get_tau_local(int dim);
double mutate_step_r(double val, double tau, double tau_prime, double u,
                     STEP_SIZE_MODE mode);
double mutate_r(double val, double step);
double mutate_step_i(double val, double tau, double tau_prime, double u,
                     STEP_SIZE_MODE mode);
double mutate_i(double val, double step, int dim);
double mutate_step_d(double val, double tau, double tau_prime, double u,
                     STEP_SIZE_MODE mode, int dim_d);
double mutate_d(double val, double step, double lb, double ub);

int main(int argc, char** argv)
{
  char* function;
  char inFile[MAX_FILENAME_LENGTH];
  char resFile[MAX_FILENAME_LENGTH];
  char bestInFile[MAX_FILENAME_LENGTH];
  char bestResFile[MAX_FILENAME_LENGTH];

  int i, j, k, l, c1, c2, cnt;
  double a, b, x, u;

  int gen = 0;

  // Variables for Sort
  int min;
  double t;
  int tt;

  double fn_old, fn_best_all_times, fn_best;
  fn_best_all_times = MAX_DOUBLE;
  int dim, lambda, mu, reco, reco_s, kappa, sel, zfa, n_better=1;
//  double tau_1, tau_2;
  int dim_r; // Number of real parameters
  int dim_i; // Number of integer parameters
  int dim_d; // Number of discrete parameters
  double tau_r; // Global learning rates for real parameters
  double tau_i; // Global learning rates for integer parameters
  double tau_d; // Global learning rates for discrete parameters
  double tau_prime_r; // Local learning rates for real parameters
  double tau_prime_i; // Local learning rates for integer parameters
  double tau_prime_d; // Local learning rates for discrete parameters
  int first_r = -1; // Index to the first continuous variable
  int first_i = -1; // Index to the first integer variable
  int first_d = -1; // Index to the first discrete variable
  STEP_SIZE_MODE step_size_mode_r = STEP_SIZE_MODE::MULTI;
  STEP_SIZE_MODE step_size_mode_i = STEP_SIZE_MODE::MULTI;
  STEP_SIZE_MODE step_size_mode_d = STEP_SIZE_MODE::SINGLE;
  //  int parameter_value[5];
  //   parameter_value[0] = DEFAULT_LAMBDA;
  //   parameter_value[1] = DEFAULT_MU;
  //   parameter_value[2] = DEFAULT_RECO;
  //   parameter_value[3] = DEFAULT_RECO_S;
  //   parameter_value[4] = DEFAULT_KAPPA;

  //   char *parameter_name[] = {"lambda", "mu", "reco", "reco_s", "kappa"};

  double xe_ini[MAX_DIM], xe[MAX_LAMBDA][MAX_DIM], xn[MAX_LAMBDA][MAX_DIM];
  double s, se[MAX_LAMBDA][MAX_DIM],  sn[MAX_LAMBDA][MAX_DIM]; 
  double lb[MAX_DIM]; // Lower bound for object vars
  double ub[MAX_DIM]; // Upper bound for object vars
  int age[MAX_LAMBDA];
  double fe[MAX_LAMBDA], fn[MAX_LAMBDA];
  TYPE type[MAX_DIM]; // 0 = continuous, 1 = integer, 2 = categorial


  /**
   * read parameters from arguments 
   */
  readParameters(argc, argv);

  srand((unsigned int)getRandomSeed());         

  dim = getDimension();
  mu = getMu();
  lambda = getLambda();
  kappa = getKappa();

  //   FORALL(i,5)
  //     {
  //       if(getConfigFileSET() && strcmp(getParameter("parameter_name[i]"), "empty") != 0) 
  // 	parameter_value[i] = atoi(getParameter("parameter_name[i]"));
  //     }
  //   lambda = parameter_value[0];
  //   mu = parameter_value[1];
  //   reco = parameter_value[2];
  //   reco_s = parameter_value[3];
  //   kappa = parameter_value[4];
  if(getConfigFileSET() 
      &&
      strcmp(getParameter("lambda"), "empty") != 0) {
    lambda = atoi(getParameter("lambda"));
  } 
  else if(lambda == 0)
    lambda = DEFAULT_LAMBDA;
  if(getConfigFileSET() 
      &&
      strcmp(getParameter("mu"), "empty") != 0) {
    mu = atoi(getParameter("mu"));
  } 
  else if(mu == 0)
    mu = DEFAULT_MU;
  if(getConfigFileSET() 
      &&
      strcmp(getParameter("reco"), "empty") != 0) {
    reco = atoi(getParameter("reco"));
  } 
  else 
    reco = DEFAULT_RECO;
  if(getConfigFileSET() 
      &&
      strcmp(getParameter("reco_s"), "empty") != 0) {
    reco_s = atoi(getParameter("reco_s"));
  } 
  else 
    reco_s = DEFAULT_RECO_S;
  if(getConfigFileSET() 
      &&
      strcmp(getParameter("kappa"), "empty") != 0) {
    kappa = atoi(getParameter("kappa"));
  } 
  else if(kappa == 0)
    kappa = DEFAULT_KAPPA;

  zfa      = getIterations();
  function = getTarget();
  sprintf(inFile, "%s.in",getProblemName());
  sprintf(resFile,"%s.res", getProblemName());
  sprintf(bestInFile,"%s.best.in", getProblemName());
  sprintf(bestResFile,"%s.best.res", getProblemName());

  cnt = 0;
  FORALL(i, dim)
  {
    lb[i]   = getLowerBound(i);    
    ub[i]   = getUpperBound(i);
    xe_ini[i]   = getXInit(i);
  }

  s = getStep(0); 

  // Retrieve types
  FORALL(i,dim)
  {
    type[i] = getType(i);
  }

  dim_r = get_dim_r(type);
  dim_i = get_dim_i(type);
  dim_d = get_dim_d(type);
  tau_r = get_tau_global(dim_r);
  tau_prime_r = get_tau_local(dim_r);
  tau_i = get_tau_global(dim_i);
  tau_prime_i = get_tau_local(dim_i);
  tau_d = get_tau_global(dim_d);
  tau_prime_d = get_tau_local(dim_d);

  //tau_1 = 1 / sqrt(2*sqrt(dim));
  //tau_2 = 1 / sqrt(2*dim);
// TODO: xe_ini is currently a wasted iteration, fix this
// TODO: More importantly: xe_ini might currently use continuous values even
//       for integer and discrete variables!!
  //setXvector(inFile, xe_ini, dim);

  //fn_best = getResultX();
  //fn_best_all_times = fn_best;

  //print(0, fn_best, fn_best, xe_ini);

  // Initialise the initial population of mu parent individuals
  FORALL(i,dim)
  {
//    type[i] = getType(i);

    if(type[i] == TYPE::CONTINUOUS && first_r == -1)
    {
      first_r = i;
    }
    else if(type[i] == TYPE::INTEGER && first_i == -1)
    {
      first_i = i;
    }
    else if(type[i] == TYPE::NOMINAL_DISCRETE && first_d == -1)
    {
      first_d = i;
    }

    FORALL(j,mu)
    {
      u = uniform();// ((double) (rand() % 10000)) / 10000.0;

      // random startvalue
      if(type[i] == TYPE::CONTINUOUS)
      {
        xe[j][i] = lb[i] + u * (ub[i]-lb[i]);

        if(step_size_mode_r == STEP_SIZE_MODE::SINGLE && i > first_r)
        {
          // Single step size, copy the first set step size of this type
          se[j][i] = se[j][first_r];
        }
        else
        {
          se[j][i] = getStepMIES(type[i], lb[i], ub[i], dim_d);
        }
      }
      else if(type[i] == TYPE::INTEGER)
      {
        xe[j][i] = ceil((lb[i] - 1) + u * (ub[i]-(lb[i] - 1)));

        if(step_size_mode_i == STEP_SIZE_MODE::SINGLE && i > first_i)
        {
          // Single step size, copy the first set step size of this type
          se[j][i] = se[j][first_i];
        }
        else
        {
          se[j][i] = getStepMIES(type[i], lb[i], ub[i], dim_d);
        }
      }
      else if(type[i] == TYPE::NOMINAL_DISCRETE)
      {
        xe[j][i] = ceil((lb[i] - 1) + u * (ub[i]-(lb[i] - 1)));

        if(step_size_mode_d == STEP_SIZE_MODE::SINGLE && i > first_d)
        {
          // Single step size, copy the first set step size of this type
          se[j][i] = se[j][first_d];
        }
        else
        {
          se[j][i] = getStepMIES(type[i], lb[i], ub[i], dim_d);
        }
      }

      age[j] = 1;

      //se[j][i] = getStepMIES(type[i], lb[i], ub[i]);//s;
      //se[j][i] = s*(ub[i]-lb[i]);
    }
  }

// printf("Generation %d:\n", gen);
  fn_best = MAX_DOUBLE;
  FORALL(j,mu) // Evaluate initial population
  {
    setXvector(inFile, xe[j], dim);
    fe[j] = getResultX();

    if(fe[j] < fn_best)
    {
      fn_best = fe[j];
      fn_best_all_times = fn_best;
      setXvector(bestInFile, xe[j], dim);
      setResult(bestResFile, fe[j]);
    }
    print(j, fe[j], fn_best_all_times, xe[j]);

    if(checkTermination(cnt, fn_best_all_times) == 1)
    {
      setResult(resFile, fn_best);
      exit(0);
    }
    ++cnt;
  }

  // printf("Initial parents:\n");
  // FORALL(j,mu)
  // {
  //   printf("f: %f; x: ", fe[j]);
  //   FORALL(k,dim - 1)
  //   {
  //     printf("%f, ", xe[j][k]);
  //   }
  //   printf("%f; ", xe[j][dim - 1]);

  //   printf("age: %d\n", age[j]);
  // }

//printf("Initial parent step sizes:\n");
//FORALL(j,lambda)
//{
//  FORALL(k,dim)
//  {
//    printf(" %f", se[j][k]);
//  }
//printf("\n");
//}

  FORALL(i,(zfa - mu)/lambda)//+1)
  {
++gen;
    FORALL(j,lambda)
    {
      c1 = (int)(/*((double) (rand() % 10000)) / 10000.0*/uniform() * ((double)mu));
      c2 = (int)(/*((double) (rand() % 10000)) / 10000.0*/uniform() * ((double)mu));

      /**********************************************
       *     GENERATE LAMBDA OFFSPRING              *
       **********************************************/

      FORALL(k, dim) 
      {
        if(reco == 1)
          xn[j][k] = (xe[c1][k] + xe[c2][k]) / 2.0;
        if(reco == 0)
          xn[j][k] = xe[c1][k];
        if(reco > 1)
        {
          u = uniform();// ((double) (rand() % 10000)) / 10000.0;
          if(u > 0.5)
            xn[j][k] = xe[c1][k];
          else
            xn[j][k] = xe[c2][k];
        }
      }

      FORALL(k, dim) 
      {
        if(reco_s == 1)
          sn[j][k] = (se[c1][k] + se[c2][k]) / 2.0;
        if(reco_s == 0)
          sn[j][k] = se[c1][k];
        if(reco_s > 1)
        {
          u = uniform();// ((double) (rand() % 10000)) / 10000.0;
          if(u > 0.5)
            sn[j][k] = se[c1][k];
          else
            sn[j][k] = se[c2][k];
        }
      }

      // Mutation 
      u = gauss();

      FORALL(k,dim) 
      {
        if(type[k] == TYPE::CONTINUOUS)
        {
          sn[j][k] = mutate_step_r(sn[j][k], tau_r, tau_prime_r, u,
                                   step_size_mode_r);
          xn[j][k] = mutate_r(xn[j][k], sn[j][k]);
        }
        else if(type[k] == TYPE::INTEGER)
        {
          sn[j][k] = mutate_step_i(sn[j][k], tau_i, tau_prime_i, u,
                                   step_size_mode_i);
          xn[j][k] = mutate_i(xn[j][k], sn[j][k], dim_i);
        }
        else if(type[k] == TYPE::NOMINAL_DISCRETE)
        {
//          int step_idx = k;
//          if(step_size_mode_d == STEP_SIZE_MODE::SINGLE)
//          {
//            step_idx = first_d;
//          }
          sn[j][k] = mutate_step_d(sn[j][k], tau_d, tau_prime_d, u,
                                   step_size_mode_d, dim_d);
          xn[j][k] = mutate_d(xn[j][k], sn[j][k], lb[k], ub[k]);
        }
        //sn[j][k] *= exp(u * tau_2 + gauss() * tau_1);
        //xn[j][k] += gauss() * sn[j][k];

        if(type[k] == TYPE::CONTINUOUS || type[k] == TYPE::INTEGER)
        {
          /**
           * Check bound (Modified Interval Bounds Treatment)
           */
          if(xn[j][k] < lb[k] || xn[j][k] > ub[k]) 
          {
            x = xn[j][k];
            a = lb[k];
            b = ub[k];
            x = (x - a) / (b - a);
            if((abs((int)x) % 2) == 0) x = fabs( x - (int)x); 
            else x = fabs(1-fabs(x - (int)x));
            xn[j][k] = a + (b - a) * x;
          }
        }
      }
    }

//printf("Mutated step sizes:\n");
//FORALL(j,lambda)
//{
//  FORALL(k,dim)
//  {
//    printf(" %f", sn[j][k]);
//  }
//printf("\n");
//}

    /* Mutation and Recobination ended, now Selection... */  

    /********************************************************************
     *   Precise Evaluation and the Produce for needed various File      *
     ********************************************************************/

    fn_old = fn_best; // the best value of previous generation among with mu plus lambda individuals 
    fn_best = MAX_DOUBLE; // for the finding of best value in new lambda individuals
    n_better = 0; // How many individuals are better than old-best value
// printf("Generation %d:\n", gen);
    FORALL(j,lambda)
    {
      setXvector(inFile, xn[j], dim);
      fn[j] = getResultX();

      if (fn[j] < fn_best) {
        fn_best=fn[j]; 
        sel=j;
      }
      if (fn[j] < fn_best_all_times) {
        fn_best_all_times=fn[j]; 
        sel=j;
        setXvector(bestInFile, xn[j], dim);
        setResult(bestResFile, fn[j]);
      }

      if (fn[j] < fn_old) n_better++;	

      // print(cnt, fn[j], fn_best_all_times, xn[j]);
      

      if(checkTermination(cnt, fn_best_all_times) == 1)
      {
        setResult(resFile, fn_best);
        exit(0);
      }
      ++cnt;
    }
    // printf("iteration %i, %.20f\n", gen, log10(fn_best_all_times));
    /********************************************************************
     *   Replace mu parent individuals by lambda best offspring         *
     ********************************************************************/

    for(j = 0; j < lambda-1; j++)
    {
      min = j;
      for(k = j+1; k < lambda; k++)
      {
        if(fn[k] < fn[min])
          min = k;
      }

      t = fn[j];
      fn[j] = fn[min];
      fn[min] = t;

      for(l = 0; l < dim; l++)
      {
        t = xn[j][l];
        xn[j][l] = xn[min][l]; 
        xn[min][l] = t;

        t = sn[j][l];
        sn[j][l] = sn[min][l]; 
        sn[min][l] = t;
      }
    }
    for(j = 0; j < mu-1; j++)
    {
      min = j;
      for(k = j+1; k < mu; k++)
      {
        if(fe[k] < fe[min])// || age[j] >= kappa)
          min = k;
      }

      t = fe[j];
      fe[j] = fe[min];
      fe[min] = t;
      tt = age[j];
      age[j] = age[min];
      age[min] = tt;

      for(l = 0; l < dim; l++)
      {
        t = xe[j][l];
        xe[j][l] = xe[min][l]; 
        xe[min][l] = t;

        t = se[j][l];
        se[j][l] = se[min][l]; 
        se[min][l] = t;
      }
    }

    /************************************************
     * MERGE SORTED POPULATIONS                     *
     ************************************************/

    double f_buf[MAX_LAMBDA];
    int age_buf[MAX_LAMBDA];
    double x_buf[MAX_LAMBDA][MAX_DIM];
    double s_buf[MAX_LAMBDA][MAX_DIM];

    int ie = 0,in = 0;
    j = 0;

    while(j < mu)
    {
//      if(fn[in] <= fe[ie] || age[ie] >= kappa)
//      {
//        if (fn[in] < fn_best) {
//          fn_best = fn[in];
//
//        }
//        f_buf[j] =  fn[in];
//        age_buf[j] = 1;
//        FORALL(l, MAX_LAMBDA)
//        {
//          x_buf[j][l] = xn[in][l];
//          s_buf[j][l] = sn[in][l];
//        }
//        in++;
//      }
//      else
//      {
//        f_buf[j] =  fe[ie];
//        age_buf[j] = age[ie] + 1;
//        FORALL(l, MAX_LAMBDA)
//        {
//          x_buf[j][l] = xe[ie][l];
//          s_buf[j][l] = se[ie][l];
//        }
//        ie++;
//      }
//      j++;
      if(fn[in] <= fe[ie] || ie >= mu)
      {
        if (fn[in] < fn_best)
        {
          fn_best = fn[in];
        }
        f_buf[j] =  fn[in];
        age_buf[j] = 1;
        FORALL(l, MAX_LAMBDA)
        {
          x_buf[j][l] = xn[in][l];
          s_buf[j][l] = sn[in][l];
        }
        ++in;
        ++j;
      }
      else
      {
        if(age[ie] < kappa)
        {
          f_buf[j] =  fe[ie];
          age_buf[j] = age[ie] + 1;
          FORALL(l, MAX_LAMBDA)
          {
            x_buf[j][l] = xe[ie][l];
            s_buf[j][l] = se[ie][l];
          }
          ++j;
        }
        ++ie;
      }
    }

    FORALL(j,mu)
    {
      fe[j] = f_buf[j];
      age[j] = age_buf[j];
      FORALL(k,dim)
      {
        xe[j][k] = x_buf[j][k];
        se[j][k] = s_buf[j][k];
      }

    }

  // printf("New parents:\n");
  // FORALL(j,mu)
  // {
  //   printf("f: %.10e; x: ", fe[j]);
  //   FORALL(k,dim - 1)
  //   {
  //     printf("%e, ", xe[j][k]);
  //   }
  //   printf("%e; ", xe[j][dim - 1]);

  //   FORALL(k,dim)
  //   {
  //     printf("%f, ", se[j][k]);
  //   }

  //   printf("age: %d\n", age[j]);
  // }

  }

  // Print final population
  // printf("Final population:\n");
  // FORALL(j,mu)
  // {
  //   printf("f: %e; x: ", fe[j]);
  //   FORALL(k,dim - 1)
  //   {
  //     printf("%f, ", xe[j][k]);
  //   }
  //   printf("%f\n", xe[j][dim - 1]);
  // }
  printf("%e", fn_best_all_times);
  return(0);
}

double getStepMIES(TYPE type, double lb, double ub, int dim_d)
{
  if(type == TYPE::CONTINUOUS || type == TYPE::INTEGER)
  {
    return 0.05 * (ub - lb);
  }
  else //(type == TYPE::NOMINAL_DISCRETE)
  {
    return 1.0 / static_cast<double>(dim_d);
  }
}

double uniform()
{
  double x = 0.0;
  while (x == 0.0)
    x = ((double)(rand() % 10000)) / 10000.0;
  return x;
}

int get_dim_r(const TYPE (&type)[MAX_DIM])
{
  int count = 0;

  for(int i = 0; i < getDimension(); ++i)
  {
    if(type[i] == TYPE::CONTINUOUS)
    {
      ++count;
    }
  }

  return count;
}

int get_dim_i(const TYPE (&type)[MAX_DIM])
{
  int count = 0;

  for(int i = 0; i < getDimension(); ++i)
  {
    if(type[i] == TYPE::INTEGER)
    {
      ++count;
    }
  }

  return count;
}

int get_dim_d(const TYPE (&type)[MAX_DIM])
{
  int count = 0;

  for(int i = 0; i < getDimension(); ++i)
  {
    if(type[i] == TYPE::NOMINAL_DISCRETE)
    {
      ++count;
    }
  }

  return count;
}

double get_tau_global(int dim)
{
  return 1 / sqrt(2 * dim);
}

double get_tau_local(int dim)
{
  return 1 / sqrt(2 * sqrt(dim));
}

// Mutate step sizes of continuous variables
double mutate_step_r(double val, double tau, double tau_prime, double u,
                     STEP_SIZE_MODE mode)
{
  if(mode == STEP_SIZE_MODE::SINGLE)
  {
    return val * exp(tau * u);
  }
  else //(mode == STEP_SIZE_MODE::MULTI)
  {
    return val * exp(u * tau + gauss() * tau_prime);
  }
}

// Mutate continuous variables
double mutate_r(double val, double step)
{
  return val + gauss() * step;
}

// Mutate step size of integer variables
double mutate_step_i(double val, double tau, double tau_prime, double u,
                     STEP_SIZE_MODE mode)
{
  double out;

  if(mode == STEP_SIZE_MODE::SINGLE)
  {
    out = val * exp(tau * u);
  }
  else //(mode == STEP_SIZE_MODE::MULTI)
  {
    out = val * exp(u * tau + gauss() * tau_prime);
  }

  if(out < 1)
  {
    out = 1;
  }

  return out;
}

// Mutate integer variables
double mutate_i(double val, double step, int dim)
{
  double u1 = uniform();
  double u2 = uniform();
  double psi = 1 - (step / dim) / (1 + sqrt(1 + pow(step / dim, 2)));
  int G1 = floor(log(1 - u1) / log(1 - psi));
  int G2 = floor(log(1 - u2) / log(1 - psi));

  return val + G1 - G2;
}

// Mutate step sizes of nominal discrete variables
double mutate_step_d(double val, double tau, double tau_prime, double u,
                     STEP_SIZE_MODE mode, int dim_d)
{
  double out;
  double lb;
  double ub = 0.5;

  if(mode == STEP_SIZE_MODE::SINGLE)
  {
    out = 1 / (1 + ((1 - val) / val) * exp((-tau) * u));
    lb = 1.0 / static_cast<double>(dim_d);
    if(lb > (1.0 / 3.0)) { lb = 1.0 / 3.0; }
  }
  else //(mode == STEP_SIZE_MODE::MULTI)
  {
    out = 1 / (1 + ((1 - val) / val) * exp((-tau) * u - tau_prime * gauss()));
    lb = 1.0 / (3.0 * static_cast<double>(dim_d));
  }

  // Modified boundary treatment to keep variable in (lb,ub)
  if(out < lb || out > ub) 
  {
    double x = out;
    double a = lb;
    double b = ub;
    x = (x - a) / (b - a);
    if((abs((int)x) % 2) == 0) x = fabs( x - (int)x); 
    else x = fabs(1-fabs(x - (int)x));
    out = a + (b - a) * x;
  }

  return out;
}

// Mutate nominal discrete variables
double mutate_d(double val, double step, double lb, double ub)
{
  double out = val;

  if(uniform() < step)
  {
    while(out == val)
    {
      out = ceil((lb - 1) + uniform() * (ub-(lb - 1)));
    }
  }

  return out;
}
