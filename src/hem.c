/**************************************************************************/
/*                                                                        */
/*                                                                        */
/*          Error Model for Anlaysis of Microarray Data                   */
/*          (C function to perform MCMC implementation)                   */
/*                                                                        */
/*           by  HyungJun Cho, PhD, and Jae K. Lee, PhD,                  */
/*         at the University of Virginia School of Medicine               */
/*                                                                        */
/*                                                                        */
/**************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "random4f.h"

static double ****obs;
static int ****mis;

static int nchip, ngene, ngroup, nrep, total_group; 
static int Brep, nquantile;
static int burnin, nsample;
static int *group;
static int **rep;
static int ***nn;

static int method_array;
static int method_bio;
static int method_total;
static int iii;

static double mu;
static double *gene;
static double *cell;
static double **inter;
static double ***expr;
static double **expr_est; 
static double **sigma2_bio;
static double max_sigma2_bio;
static double ***sigma2_array;

static double mean_mu;
static double *mean_gene;
static double *mean_cell;
static double **mean_inter;
static double ***mean_expr;
static double **mean_expr_est; 
static double **mean_sigma2_bio;
static double ***mean_sigma2_array;

static double sigma2g;
static double sigma2c;
static double sigma2r;

static double alpha_bio;
static double beta_bio;
static double **lpe_bio;

static double alpha_array;
static double beta_array;
static double **lpe_array;
static double ***boot_array;
static double **quantile_array;

#define SEP ','	       
#define MAX_ENTRY 20   


/* ****************************************************************** */
/*                                                                    */
/* READ DATA                                                          */
/*                                                                    */
/* fp: data file (seperator ,)                                        */
/* fpar: parameter file (seperator ,)                                 */
/*                                                                    */
/* ngene: number of genes                                             */
/* nchip: number of chips                                             */
/* ngroup: number of groups (or conditions)                           */
/* group[ngroup]: group size                                          */
/* nrep: number of replicates (if 0, various replicates)              */
/* rep : number of replicates                                         */
/* nn : same as rep (but -1 per a missing value)                      */
/*                                                                    */
/* burnin: number of burn-ins                                         */
/* nsample: number of MCMC samples after bur-ins                      */
/* printout1,printout2: if 1, print out                               */
/* sigma2g, sigma2c, sigma2r: priors of gene,condition,interaction    */
/* alpha_bio, beta_bio: Beta priors of biological error               */
/* alpha_array, beta_array: Beta priors of array error                */
/* max_sigma2_bio: max limit of sigma2_bio                            */
/*                                                                    */
/* obs: data vector (size ngene*nchip)                                */
/* mis: missing vector (1=missing, 0=exist)                           */
/*                                                                    */
/* ****************************************************************** */

void read_data_file(double *dat, int *grp, int *nrp)  
{

	int	i, j, k, l, n;
 
	/* Read parameters (group size) */
        total_group=0;
 	group = (int *) calloc(ngroup, sizeof(int));        
	for (j = 0; j < ngroup; j++)
	{
             group[j] = grp[j]; 
	     total_group += group[j];
	}


	/* Read parameters (number of replicates: rep=nn for each gene) */
 	rep = (int **) calloc(ngroup, sizeof(int));
        for(j = 0; j < ngroup; j++) rep[j] = (int *) calloc(group[j], sizeof(int)); 

 	nn = (int ***) calloc(ngene, sizeof(int));
        for(i = 0; i < ngene; i++) nn[i] = (int **) calloc(ngroup, sizeof(int)); 
        for(i = 0; i < ngene; i++) 
        {
            for(j = 0; j < ngroup; j++) 
            {
                nn[i][j] = (int *) calloc(group[j], sizeof(int));
            }
        }

	n = 0;
        for (j = 0; j < ngroup; j++)
	{
             for (k = 0; k < group[j]; k++)
	     {                      
                  rep[j][k] = nrp[n];
 	          for(i = 0; i < ngene; i++) nn[i][j][k] = nrp[n];
                  n += 1;
             }
        }

	/* Allocate space */
	obs = (double ****) calloc(ngene, sizeof(double));
	mis = (int ****) calloc(ngene, sizeof(int));
	gene = (double *) calloc(ngene, sizeof(double));
	cell = (double *) calloc(ngroup, sizeof(double));
	expr = (double ***) calloc(ngene, sizeof(double));
	inter = (double **) calloc(ngene, sizeof(double));
	sigma2_bio = (double **) calloc(ngene, sizeof(double));
	sigma2_array = (double ***) calloc(ngene, sizeof(double));
	
        for(i = 0; i < ngene; i++) 
        {
            obs[i] = (double ***) calloc(ngroup, sizeof(double)); 
            mis[i] = (int ***) calloc(ngroup, sizeof(int)); 
            expr[i] = (double **) calloc(ngroup, sizeof(double)); 
            inter[i] = (double *) calloc(ngroup, sizeof(double)); 
            sigma2_bio[i] = (double *) calloc(ngroup, sizeof(double)); 
            sigma2_array[i] = (double **) calloc(ngroup, sizeof(double)); 
	   	   
       }

        for(i = 0; i < ngene; i++) 
        {
            for(j = 0; j < ngroup; j++) 
            {  
                obs[i][j] = (double **) calloc(group[j], sizeof(double));
                mis[i][j] = (int **) calloc(group[j], sizeof(int));
                expr[i][j] = (double *) calloc(group[j], sizeof(double));
                sigma2_array[i][j] = (double *) calloc(group[j], sizeof(double)); 
           }
        }

        for(i = 0; i < ngene; i++) 
        {
            for(j = 0; j < ngroup; j++) 
            {  
                for(k = 0; k < group[j]; k++) 
                {  
		    obs[i][j][k] = (double *) calloc(nn[i][j][k], sizeof(double));
		    mis[i][j][k] = (int *) calloc(nn[i][j][k], sizeof(int));
                }
            }
        }


        /*Read data*/
        n = 0;
	for (i = 0; i < ngene; i++)
	{
		for (j = 0; j < ngroup; j++)
  		{
		for (k = 0; k < group[j]; k++)
		{
		for (l = 0; l < rep[j][k]; l++)
		{
		       mis[i][j][k][l] = 0;
		       obs[i][j][k][l] = dat[n];
                       if(obs[i][j][k][l] <= -9) {
                          obs[i][j][k][l] = -9;
         		  mis[i][j][k][l] = 1;
                       }
	               n += 1;		
		      /*printf(" %lf ",obs[i][j][k][l]);*/
		}}}
		/*printf("==> %d \n ",i);*/
	}

}


/* ****************************************************************** */
/*                                                                    */
/* READ DATA                                                          */
/* (One Layer HEM)                                                    */
/*                                                                    */
/* ****************************************************************** */

void read_data_file_onelayer(double *dat, int *grp)  
{
	int	i, j, k, l, n;
 
	/* Read parameters (group size) */
        total_group=0;
 	group = (int *) calloc(ngroup, sizeof(int));        
	for (j = 0; j < ngroup; j++)
	{
             group[j] = grp[j]; 
	     total_group += group[j];
	}


	/* Read parameters (number of replicates: rep=nn for each gene) */
        nrep = 1;
 	rep = (int **) calloc(ngroup, sizeof(int));
        for(j = 0; j < ngroup; j++) rep[j] = (int *) calloc(group[j], sizeof(int)); 

 	nn = (int ***) calloc(ngene, sizeof(int));
        for(i = 0; i < ngene; i++) nn[i] = (int **) calloc(ngroup, sizeof(int)); 
        for(i = 0; i < ngene; i++) 
        {
            for(j = 0; j < ngroup; j++) 
            {
                nn[i][j] = (int *) calloc(group[j], sizeof(int));
            }
        }

        
            for (j = 0; j < ngroup; j++)
	    {
                 for (k = 0; k < group[j]; k++)
	        {
                      rep[j][k]=nrep;
 	              for(i = 0; i < ngene; i++) nn[i][j][k] = nrep;
                }
            }

	/* Allocate space */
	obs = (double ****) calloc(ngene, sizeof(double));
	mis = (int ****) calloc(ngene, sizeof(int));
	gene = (double *) calloc(ngene, sizeof(double));
	cell = (double *) calloc(ngroup, sizeof(double));
	expr = (double ***) calloc(ngene, sizeof(double));
	inter = (double **) calloc(ngene, sizeof(double));
	sigma2_bio = (double **) calloc(ngene, sizeof(double));
	
        for(i = 0; i < ngene; i++) 
        {
            obs[i] = (double ***) calloc(ngroup, sizeof(double)); 
            mis[i] = (int ***) calloc(ngroup, sizeof(int)); 
            expr[i] = (double **) calloc(ngroup, sizeof(double)); 
            inter[i] = (double *) calloc(ngroup, sizeof(double)); 
            sigma2_bio[i] = (double *) calloc(ngroup, sizeof(double)); 
	   	   
        }

        for(i = 0; i < ngene; i++) 
        {
            for(j = 0; j < ngroup; j++) 
            {  
                obs[i][j] = (double **) calloc(group[j], sizeof(double));
                mis[i][j] = (int **) calloc(group[j], sizeof(int));
                expr[i][j] = (double *) calloc(group[j], sizeof(double));
           }
        }

        for(i = 0; i < ngene; i++) 
        {
            for(j = 0; j < ngroup; j++) 
            {  
                for(k = 0; k < group[j]; k++) 
                {  
		    obs[i][j][k] = (double *) calloc(nn[i][j][k], sizeof(double));
		    mis[i][j][k] = (int *) calloc(nn[i][j][k], sizeof(int));
                }
            }
        }

 
        /*Read data*/
        n = 0;
	for (i = 0; i < ngene; i++)
	{
		for (j = 0; j < ngroup; j++)
  		{
		for (k = 0; k < group[j]; k++)
		{
		for (l = 0; l < rep[j][k]; l++)
		{
		       mis[i][j][k][l] = 0;
		       obs[i][j][k][l] = dat[n];
                       if(obs[i][j][k][l] <= 0.0) {
                          obs[i][j][k][l] = 0.0;
         		  mis[i][j][k][l] = 1;
                       }
	               n = n +1;		
		       /*printf(" %lf ",obs[i][j][k][l]);*/
		}}}
		/*printf("==> %d \n ",i);*/
	}
}


/* ****************************************************************** */
/*                                                                    */
/* READ PRIORS                                                        */
/* (Two layer EM)                                                     */
/*                                                                    */
/* ****************************************************************** */
void read_prior_file(int *Bsize, double *vare, double *varb)  
{
	int	i, j, k, l, n;
        /*double  mean_array, var_array;*/

 
    /*Read bootstrap samples for exp error*/
    if(method_array ==3) {
	boot_array = (double ***) calloc(ngroup, sizeof(double));
	quantile_array = (double **) calloc(ngroup, sizeof(double));

        Brep = Bsize[0];
        nquantile = Bsize[1];  

        for(j = 0; j < ngroup; j++) 
        {
            boot_array[j] = (double **) calloc(Brep, sizeof(double));
	    quantile_array[j] = (double *) calloc(nquantile, sizeof(double));
        }

        for(j = 0; j < ngroup; j++) 
        {
            for(k = 0; k < Brep; k++) 
            {
                boot_array[j][k] = (double *) calloc(nquantile, sizeof(double));
	    }
	}	


        n = 0;
	for (j = 0; j < ngroup; j++)
	{
		for (k = 0; k < Brep; k++)
  		{
      		        for (l = 0; l < nquantile; l++)
 	        	{
                             boot_array[j][k][l] = vare[n];	     
                             n +=1;
               	             /*printf(" %lf  ", boot_array[j][k][l]);*/
                        }
			/*printf("==> array error %d \n ",j);*/
 		}
        	for (l = 0; l < nquantile; l++)
	        {
                     quantile_array[j][l] = vare[n];
                     n +=1;
		     /*printf("  %lf ", quantile_array[j][l]);*/
                }
		/*printf("==> quantile %d \n ",j);*/
	}
    }
  

    /*Read priors for exp error (Parametric empirical Bayes HEM)*/
    if(method_array ==2) {
	lpe_array = (double **) calloc(nquantile, sizeof(double));
	quantile_array = (double **) calloc(nquantile, sizeof(double));

        for(i = 0; i < ngene; i++) 
        {
            lpe_array[i] = (double *) calloc(ngroup, sizeof(double));
            quantile_array[i] = (double *) calloc(ngroup, sizeof(double));
        }
	

        n = 0;
	for (j = 0; j < ngroup; j++)
	{

      		for (l = 0; l < nquantile; l++)
 	        {
                     lpe_array[j][l] = vare[n];	     
                     n +=1;
                }

        	for (l = 0; l < nquantile; l++)
	        {
                     quantile_array[j][l] = vare[n];
                     n +=1;
                }
	}
    }


    /*Read priors for biological error (Parametric empirical Bayes EM)*/
    if(method_bio ==2) {
	lpe_bio = (double **) calloc(ngene, sizeof(double));
        for(i = 0; i < ngene; i++) 
        {
            lpe_bio[i] = (double *) calloc(ngroup, sizeof(double));
        }

        n =0;	
	for (i = 0; i < ngene; i++)
	{
		for (j = 0; j < ngroup; j++)
  		{
                     lpe_bio[i][j] = varb[n];
                     n += 1;     
 		}
	}
    }

}


/* ****************************************************************** */
/*                                                                    */
/* READ PRIORS                                                        */
/* (One layer HEM)                                                    */
/*                                                                    */
/* ****************************************************************** */
void read_prior_file_onelayer(int *Bsize, double *vart)  
{
	int	i, j, k, l, n;
        /*double  mean_array, var_array;*/

 
    /*Read bootstrap samples for exp error*/
    if(method_total ==3) {
	boot_array = (double ***) calloc(ngroup, sizeof(double));
	quantile_array = (double **) calloc(ngroup, sizeof(double));

        Brep = Bsize[0];
        nquantile = Bsize[1];  

        for(j = 0; j < ngroup; j++) 
        {
            boot_array[j] = (double **) calloc(Brep, sizeof(double));
	    quantile_array[j] = (double *) calloc(nquantile, sizeof(double));
        }

        for(j = 0; j < ngroup; j++) 
        {
            for(k = 0; k < Brep; k++) 
            {
                boot_array[j][k] = (double *) calloc(nquantile, sizeof(double));
	    }
	}	


        n = 0;
	for (j = 0; j < ngroup; j++)
	{

		for (k = 0; k < Brep; k++)
  		{
      		        for (l = 0; l < nquantile; l++)
 	        	{
                             boot_array[j][k][l] = vart[n];	     
                             n +=1;
                       }
 		}
        	for (l = 0; l < nquantile; l++)
	        {
                     quantile_array[j][l] = vart[n];
                     n +=1;
                }
	}


    }


   /*Read priors for biological error (Parametric empirical Bayes HEM)*/
    if(method_total ==2) {
	lpe_bio = (double **) calloc(ngene, sizeof(double));
        for(i = 0; i < ngene; i++) 
        {
            lpe_bio[i] = (double *) calloc(ngroup, sizeof(double));
        }

        n =0;	
	for (i = 0; i < ngene; i++)
	{
		for (j = 0; j < ngroup; j++)
  		{
                     lpe_bio[i][j] = vart[n];
                     n += 1;     
 		}
	}
    }


}


/* ****************************************************************** */
/*                                                                    */
/* INITIALIZE THE PARAMETER VALUES                                    */
/*                                                                    */
/* ****************************************************************** */

void initialize()
{
	int i,j,k,l, n, n1;
	/*double tmp, alpha;*/

	for (i=0 ; i < ngene; i++)
	{
	     for (j=0; j < ngroup; j++)
	     {
		  for (k=0; k < group[j]; k++)
		  {
          	       for (l=0; l < rep[j][k] ; l++) expr[i][j][k] += obs[i][j][k][l];
 		       inter[i][j] += expr[i][j][k]/group[j];
 	          }
	     }
	}

	for (j = 0; j < ngroup; j++) {
             for (i = 0; i < ngene; i++) {
		  gene[i] += inter[i][j] / ngroup;
                  cell[j] += inter[i][j] / ngene;
	     }
        }
	
        if(method_array==4) { /*constant array error variance*/
	     sigma2_array[0][0][0] = 1./ Gamma(alpha_array, beta_array);
             for (i=0 ; i < ngene; i++)
	     {
	         for (j=0; j < ngroup; j++)
	         {
        		for (k=0; k < group[j]; k++)
	        	{   
	                     sigma2_array[i][j][k] = sigma2_array[0][0][0];
			}
	         }
	    }
	}

        if(method_array==1) { 
             for (i=0 ; i < ngene; i++)
	     {
	         for (j=0; j < ngroup; j++)
	         {
        		for (k=0; k < group[j]; k++)
	        	{   
	                     sigma2_array[i][j][k] = 1./ Gamma(alpha_array, beta_array);
			}
	         }
	    }
	}

        if(method_array==2) {
    	      for (i=0; i < ngene; i++)
	      {
		   for (j=0; j < ngroup; j++)
		   {
        		for (k=0; k < group[j]; k++)
	        	{   
			         n=nquantile-1;
                     		 for (l=0; l < nquantile-1; l++)
	                	 {
				      if (expr[i][j][k] <= quantile_array[l][j]) {
					  n = l;
                                          break;
                                      } 
                                 }
          	                 sigma2_array[i][j][k] =  lpe_array[n][j];
 			}  
		    }
	      }
	}

        if(method_array==3) {
    	      for (i=0; i < ngene; i++)
	      {
		   for (j=0; j < ngroup; j++)
		   {
        		for (k=0; k < group[j]; k++)
	        	{   
			         n=nquantile-1;
                     		 for (l=0; l < nquantile-1; l++)
	                	 {
				      if (expr[i][j][k] <= quantile_array[j][l]) {
					  n = l;
                                          break;
                                      } 
                                 }

                                 n1 = RandomUniform()*(Brep-1);
          	                 sigma2_array[i][j][k] =  boot_array[j][n1][n];
 			}  
		   }
	      }
	}

        if(method_bio==1) {
             for (i=0 ; i < ngene; i++)
	     {
	         for (j=0; j < ngroup; j++)
	         {
	              sigma2_bio[i][j] = 1./ Gamma(alpha_bio, beta_bio);
	         }
	    }
	}


        if(method_bio==2) {
    	      for (i=0; i < ngene; i++)
	      {
		   for (j=0; j < ngroup; j++)
		   {
          	        sigma2_bio[i][j] =  lpe_bio[i][j];
		   }
	      }
	}


        if(method_total==1) {
             for (i=0 ; i < ngene; i++)
	     {
	         for (j=0; j < ngroup; j++)
	         {
	              sigma2_bio[i][j] = 1./ Gamma(alpha_bio, beta_bio);
	         }
	    }
	}


        if(method_total==2) {
    	      for (i=0; i < ngene; i++)
	      {
		   for (j=0; j < ngroup; j++)
		   {
          	        sigma2_bio[i][j] =  lpe_bio[i][j];
 		   }
	      }
	}


        if(method_total==3) {
       	      for (i=0; i < ngene; i++)
	      {
		   for (j=0; j < ngroup; j++)
		   {
			         n=nquantile-1;
                     		 for (l=0; l < nquantile-1; l++)
	                	 {
				      if (inter[i][j] <= quantile_array[j][l]) {
					  n = l;
                                          break;
                                      } 
                                 }
                                 n1 = RandomUniform()*(Brep-1);
          	                 sigma2_bio[i][j] =  boot_array[j][n1][n];  
		   }
	      }
	}


        n =  0;
	for (i=0; i < ngene; i++){
             for (j=0; j < ngroup; j++){
                  for (k=0; k < group[j]; k++){
                       for (l=0; l < nn[i][j][k]; l++) {
                            mu += obs[i][j][k][l];
                       }
                       n  +=  nn[i][j][k];
                  }
              }
        }
	        mu = mu / n;

	for (j = 0; j < ngroup; j++) cell[j] = cell[j] - mu;
	for (i = 0; i < ngene; i++) gene[i] = gene[i] - mu;
	for (i = 0; i < ngene; i++) {
             for (j = 0; j < ngroup; j++) inter[i][j] = inter[i][j] - mu;
        }
}

/* ****************************************************************** */
/*                                                                    */
/* UPDATE PARAMETERS AND MISSING DATA BY GIBBS SAMLPING               */
/*                                                                    */
/* ****************************************************************** */

double update_mu()
{
	int i, j, k;
	double mean, std, var;

	mean = 0.;
	var  = 0.;

	for (i=0; i < ngene; i++)
	{
		for (j=0; j < ngroup; j++)
		{
			var +=  group[j] / sigma2_bio[i][j];
		}
	}
	var = 1. / var;

	for (i=0; i < ngene; i++)
	{
		for (j=0; j < ngroup; j++)
		{
			for (k=0; k < group[j]; k++)
			{
				mean += var * ( expr[i][j][k] -  gene[i] -
					cell[j] - inter[i][j] ) / sigma2_bio[i][j];
			}
		}
	}

	std = sqrt( var );
	mu = std * StdNormal() + mean;
	
	return var;
}


void update_gene()
{
	int i,j,k;
	double mean, std, var;

	for (i=0 ; i < ngene; i++)
	{
		var  = 1. / sigma2g;
		for (j=0; j < ngroup; j++)
		{
			var +=  group[j] / sigma2_bio[i][j];
		}
		var = 1. / var;

		mean = 0.;
		for (j=0; j < ngroup; j++)
		{
			for (k=0; k < group[j]; k++)
			{
				mean += var * ( expr[i][j][k] - mu -  cell[j] -
					inter[i][j] ) / sigma2_bio[i][j];
			}
		}
		std = sqrt( var );
		gene[i] = std * StdNormal() + mean;
	}
}


void update_cell()
{
	int i,j,k;
	double mean, std, var;


	for (j=0 ; j < ngroup; j++)
	{
		var  = 1. / sigma2c;
		for (i=0; i < ngene; i++)
		{
			var +=  group[j] / sigma2_bio[i][j];
		}
		var = 1. / var;

		mean = 0.;
		for (i=0; i < ngene; i++)
		{
			for (k=0; k < group[j]; k++)
			{
				mean += var * (expr[i][j][k] - mu -  gene[i] -
					inter[i][j]) / sigma2_bio[i][j];
			}
		}

		std = sqrt( var );
		cell[j] = std * StdNormal() + mean;
	}

}


void update_inter()
{
	int i,j,k;
	double mean, std;

	for (i=0 ; i < ngene; i++)
	{
		for (j=0; j < ngroup; j++)
		{
			mean = 0.;
			for (k=0; k < group[j]; k++)
			{
				mean += sigma2r / ( sigma2_bio[i][j] / group[j] + sigma2r )*
					( expr[i][j][k] - mu - gene[i] - cell[j] ) /
					group[j];
			}
			std = sqrt( 1. / ( group[j] / sigma2_bio[i][j] + 1. / sigma2r ) );
			inter[i][j] = std * StdNormal() + mean;
		}
	}

}



void update_expr()
{
	int i,j,k,l;
	double mean, std;

	for (i=0 ; i < ngene; i++)
	{
		for (j=0; j < ngroup; j++)
		{
			for (k=0; k < group[j]; k++)
			{
			        if (nn[i][j][k] > 0) /*exist*/
				{
        			        mean =0.;
                			for (l=0; l < rep[j][k]; l++) /*rep is the same for each gene*/
	        		        {
					     mean += obs[i][j][k][l];
 					}
        				mean = sigma2_bio[i][j]  *  mean / nn[i][j][k];
					mean += sigma2_array[i][j][k] / nn[i][j][k] *
  	                                        (mu + gene[i] + cell[j] + inter[i][j]);
                                        mean = mean / (sigma2_bio[i][j] + sigma2_array[i][j][k] / nn[i][j][k]);
	         			std = sqrt( 1. / ( 1. / sigma2_bio[i][j] + nn[i][j][k] / sigma2_array[i][j][k] ) );
				} else /*missing*/
				{
			        	mean =  mu + gene[i] + cell[j] + inter[i][j];
				        std = sqrt( sigma2_bio[i][j] );
				}
				expr[i][j][k] = std * StdNormal() + mean;
                               /*if(i==52 && j==0) printf("expr= %lf <= %lf %lf \n", expr[i][j][k], mean, std);*/
			}
		}
	}

}



           

void update_sigma2_bio()
{
	int i,j,k;
	double beta, tmp, alpha;


	for (i=0; i < ngene; i++)
	{
		for (j=0; j < ngroup; j++)
		{
			alpha = group[j] / 2. + alpha_bio;
                        if(method_bio ==2) beta_bio = lpe_bio[i][j] * (alpha_bio-1); /*parametric EB*/
			beta  = beta_bio;


			for (k=0; k < group[j]; k++)
			{
			     tmp = expr[i][j][k] - (mu + gene[i] + cell[j] + inter[i][j]);
			     beta += tmp*tmp/2.;
 			}

			sigma2_bio[i][j] = 1. / Gamma(alpha, beta);
			if(sigma2_bio[i][j] > max_sigma2_bio) {
                           sigma2_bio[i][j] = (max_sigma2_bio + beta_bio)/2.;
                        }
		}
	}
}


void update_sigma2_array_const()
{
	int i, j, k, l;
	double alpha, beta, tmp;

	alpha= alpha_array;
	beta = beta_array;

	for (i=0; i < ngene; i++)
	{
		for (j=0; j < ngroup; j++)
		{
			for (k=0; k < group[j]; k++)
			{
         			for (l=0; l < rep[j][k]; l++)
			        {
				     tmp = (obs[i][j][k][l] - expr[i][j][k])*(1 - mis[i][j][k][l]) ;
          			     beta += tmp*tmp/2;
                                }                         
     	                        alpha += nn[i][j][k]/2.0;
			}
		}
	}


        sigma2_array[0][0][0] = 1. / Gamma(alpha, beta);  
        for (i=0 ; i < ngene; i++)
	{
	     for (j=0; j < ngroup; j++)
	     {
        	  for (k=0; k < group[j]; k++)
	          {   
	               sigma2_array[i][j][k] = sigma2_array[0][0][0];
		  }
	     }
        }


}


void update_sigma2_array()
{
	int i,j,k,l,n;
	double alpha, beta, tmp;

	for (i=0; i < ngene; i++)
	{
		for (j=0; j < ngroup; j++)
		{
			for (k=0; k < group[j]; k++)
			{
                                if(method_array ==2) { /*parametric EB*/
			             n=nquantile-1;
                     		     for (l=0; l < nquantile-1; l++)
	                	     {
				          if (expr[i][j][k] <= quantile_array[l][j]) {
					      n = l;
                                              break;
                                          } 
                                     }
                                     beta_array = lpe_array[n][j] * (alpha_array-1);
                                }

      		                beta =0.;
         			for (l=0; l < rep[j][k]; l++)
			        {
				     tmp = (obs[i][j][k][l] - expr[i][j][k])*(1 - mis[i][j][k][l]) ;
          			     beta += tmp*tmp;
                                }                         
 	                        alpha = nn[i][j][k] / 2. + alpha_array;
	                        beta = beta / 2. + beta_array;
        	                sigma2_array[i][j][k] = 1. / Gamma(alpha, beta);  
 			}
		}
	}
}



void update_sigma2_array_nonpar()
{
	int i, j, k, l, n, n1;
	double temp, temp2, MHratio, U, old, new, raccept;

        raccept=0.0;
	for (i=0; i < ngene; i++)
	{
		for (j=0; j < ngroup; j++)
		{
        		for (k=0; k < group[j]; k++)
	        	{

			         n=nquantile-1;
                     		 for (l=0; l < nquantile-1; l++)
	                	 {
				      if (expr[i][j][k] <= quantile_array[j][l]) {
					  n = l;
                                          break;
                                      } 
                                 }
                
                                 n1 = RandomUniform()*(Brep-1);
                                 new = boot_array[j][n1][n];
                                 old = sigma2_array[i][j][k];                        			
  
                                 MHratio=1;
	        	         for (l=0; l < nn[i][j][k]; l++) MHratio *= sqrt(old/new);

                                 temp2=0.0;
         		         for (l=0; l < rep[j][k]; l++)
			         {
			              temp = (obs[i][j][k][l] - expr[i][j][k])*(1 - mis[i][j][k][l]) ;
          		              temp2 += temp*temp/2;
                                 } 

                                 MHratio *= exp(temp2*(1.0/old-1.0/new));
                                 if(MHratio < 1) {
                                      U = RandomUniform();
		                      if(U <= MHratio) sigma2_array[i][j][k] = new;
                                      else raccept ++;
				 } else
				 {
                                      sigma2_array[i][j][k] = new;
				 }                
 			}
		}
	}
	/*printf("sigma2_array= %lf \n", sigma2_array[1][1][1]);*/
	/*printf("Rejection rate= %10.5f (#accepted= %5.0f/%d) \n", raccept/(ngene*60),ngene*60.0-raccept,ngene*60 );*/
}


void update_expr_total()
{
	int i,j,k;

	for (i=0 ; i < ngene; i++)
	{
		for (j=0; j < ngroup; j++)
		{
			for (k=0; k < group[j]; k++)
			{
			      
				expr[i][j][k] = obs[i][j][k][0];
			}
		}
	}

}


void update_sigma2_total()
{
	int i,j,k;
	double beta, tmp, alpha;

	for (i=0; i < ngene; i++)
	{
		for (j=0; j < ngroup; j++)
		{
			alpha= group[j] / 2. + alpha_bio;
                        if(method_total ==2) beta_bio = lpe_bio[i][j] * (alpha_bio-1); /*param EB*/
			beta  = beta_bio;

			for (k=0; k < group[j]; k++)
			{
			     tmp = expr[i][j][k] - mu -  gene[i] - cell[j] - inter[i][j];
			     beta += tmp*tmp/2.;
 			}

			sigma2_bio[i][j] = 1. / Gamma(alpha, beta);
			if(sigma2_bio[i][j] > max_sigma2_bio) sigma2_bio[i][j] = (max_sigma2_bio + beta_bio)/2.;
		}
	}
}



void update_sigma2_total_nonpar()
{
	int i, j, k, l, n, n1;
	double temp, temp2, MHratio, U, old, new, raccept;

        raccept=0.0;
	for (i=0; i < ngene; i++)
	{
		for (j=0; j < ngroup; j++)
		{
        	
			         n=nquantile-1;
                     		 for (l=0; l < nquantile-1; l++)
	                	 {
				      if (expr_est[i][j] <= quantile_array[j][l]) {
					  n = l;
                                          break;
                                      } 
                                 }
                
                                 n1 = RandomUniform()*(Brep-1);
                                 new = boot_array[j][n1][n];
                                 old = sigma2_bio[i][j];                        			
  
                                 MHratio=1;
	        	         for (k=0; k < group[j]; k++) MHratio *= sqrt(old/new);

                                 temp2=0.0;
         		         for (k=0; k < group[j]; k++)
			         {
			             temp = (expr[i][j][k]-expr_est[i][j]) ;
          		             temp2 += temp*temp/2;
                                 } 

                                 MHratio *= exp(temp2*(1.0/old-1.0/new));
                                 if(MHratio < 1) {
                                    U = RandomUniform();
		                    if(U <= MHratio) sigma2_bio[i][j] = new;
                                    else raccept ++;
				 } else
				 {
                                    sigma2_bio[i][j] = new;
				 }                
		}
	}
	/*printf("Rejection rate= %lf (#accepted= %5.0f/%d) \n", raccept/(ngene*60),ngene*60.0-raccept,ngene*60 );*/
}


/* ****************************************************************** */
/*                                                                    */
/*    MAIN FUNCTION                                                   */
/*                                                                    */
/* ****************************************************************** */

void twolayerhem(double *dat, int *opt, int *ndat, int *grp, int *nrp, int *nMCMC, double *par, int *Bsize, double *vare, double *varb,  double *fstat, double *mexprest, double *mexpr, double *msigma2b, double *msigma2e, double *MCMCsamp) 
{
	int i, j, k, m, n;
        double mean, var;

        method_array = opt[0];
        method_bio   = opt[1];
        method_total = opt[2];

        ngene  = ndat[0];
        nchip  = ndat[1];
        ngroup = ndat[2];

        burnin  = nMCMC[0];
        nsample = nMCMC[1];

        sigma2g     = par[0];
        sigma2c     = par[1];
        sigma2r     = par[2];
        alpha_bio   = par[3];
        beta_bio    = par[4];
        alpha_array = par[5];
        beta_array  = par[6];
        max_sigma2_bio = 10;



	/* Read data */
	read_data_file(dat, grp, nrp);

 	/* Empirical Bayes */
        if(method_bio==2) read_prior_file(Bsize, vare, varb);  


	/* Allocate memory */
	mean_gene = (double *) calloc(ngene, sizeof(double));
	mean_cell = (double *) calloc(ngroup, sizeof(double));
	mean_inter = (double **) calloc(ngene, sizeof(double));
	mean_sigma2_bio = (double **) calloc(ngene, sizeof(double));
	mean_sigma2_array = (double ***) calloc(ngene, sizeof(double));
	mean_expr = (double ***) calloc(ngene, sizeof(double));
	expr_est = (double **) calloc(ngene, sizeof(double));
	mean_expr_est = (double **) calloc(ngene, sizeof(double));

        for(i = 0; i < ngene; i++) 
        {
            mean_inter[i] = (double *) calloc(ngroup, sizeof(double)); 
            mean_sigma2_bio[i] = (double *) calloc(ngroup, sizeof(double)); 
            mean_sigma2_array[i] = (double **) calloc(ngroup, sizeof(double)); 
            mean_expr[i] = (double **) calloc(ngroup, sizeof(double)); 
            expr_est[i] = (double *) calloc(ngroup, sizeof(double)); 
            mean_expr_est[i] = (double *) calloc(ngroup, sizeof(double)); 
        }
        for(i = 0; i < ngene; i++) 
        {
            for(j = 0; j < ngroup; j++) 
            {  
               mean_expr[i][j] = (double *) calloc(group[j], sizeof(double));
               mean_sigma2_array[i][j] = (double *) calloc(group[j], sizeof(double)); 
            }
        }

	/* Set zero for output variables */
	mean_mu = 0.;
	for (j = 0; j < ngroup; j++) mean_cell[j] = 0.;

	for (i = 0; i < ngene; i++)
	{
		mean_gene[i] = 0.;
		for (j = 0; j < ngroup; j++)
		{
			mean_inter[i][j] = 0.;
			mean_sigma2_bio[i][j] = 0.;
  			mean_expr_est[i][j] = 0.;

			for (k = 0; k < group[j]; k++)
			{
				mean_expr[i][j][k] = 0.;
                        	mean_sigma2_array[i][j][k] = 0.;
			}
		}
	}


	/* Set and print initial values*/ 
	initialize();

   
	/* Save MCMC samples*/
        m=0;i=0;j=0;k=0;	     
        MCMCsamp[m] = expr[i][j][k]; m = m + 1;
        MCMCsamp[m] = sigma2_array[i][j][k]; m = m + 1;
        MCMCsamp[m] = mu; m = m + 1;
        MCMCsamp[m] = gene[i]; m = m + 1;
        MCMCsamp[m] = cell[j]; m = m + 1;
        MCMCsamp[m] = inter[i][j]; m = m + 1;
        MCMCsamp[m] = sigma2_bio[i][j]; m = m + 1;
      
	
        /*** Burn-ins ***/
        k = 0 ;
	for (n = 0; n < burnin; n++)
	{

             update_expr();
             update_inter();
             update_gene();
             update_cell();
             var = update_mu();
             update_sigma2_bio();

             switch(method_array)
             {
             case 3:
                  update_sigma2_array_nonpar();
                  break;
             case 4:
                  update_sigma2_array_const();
                  break;
             default:
                  update_sigma2_array();
             }

	    
             /* Save MCMC samples*/
             i=0;j=0;k=0;	     
             MCMCsamp[m] = expr[i][j][k]; m = m + 1;
             MCMCsamp[m] = sigma2_array[i][j][k]; m = m + 1;
             MCMCsamp[m] = mu; m = m + 1;
             MCMCsamp[m] = gene[i]; m = m + 1;
             MCMCsamp[m] = cell[j]; m = m + 1;
             MCMCsamp[m] = inter[i][j]; m = m + 1;
             MCMCsamp[m] = sigma2_bio[i][j]; m = m + 1;
	    
             k +=1; 
             if(k == 50) {k = 0; printf(".");}
             /*printf("iter= %d \n",n);*/

        }


        /*** Keep MCMC samples ***/
	for (n = 0; n < nsample; n++)
	{
             update_expr();
             update_inter();
             update_gene();
             update_cell();
             var = update_mu();
             update_sigma2_bio();

             switch(method_array)
             {
             case 3:
                  update_sigma2_array_nonpar();
                  break;
             case 4:
                  update_sigma2_array_const();
                  break;
             default:
                  update_sigma2_array();
             }


	    
             /* Save MCMC samples*/
             i=0;j=0;k=0;	     
             MCMCsamp[m] = expr[i][j][k]; m = m + 1;
             MCMCsamp[m] = sigma2_array[i][j][k]; m = m + 1;
             MCMCsamp[m] = mu; m = m + 1;
             MCMCsamp[m] = gene[i]; m = m + 1;
             MCMCsamp[m] = cell[j]; m = m + 1;
             MCMCsamp[m] = inter[i][j]; m = m + 1;
             MCMCsamp[m] = sigma2_bio[i][j]; m = m + 1;
	    
             /*
             k +=1; 
             if(k == 50) { k = 0; printf(".");}
             printf("iter= %d \n",n);
             */
 

  
	     mean_mu += mu / nsample;
     	     for (j = 0; j < ngroup; j++) mean_cell[j] += cell[j] / nsample;

	     for (i = 0; i < ngene; i++)
             {
        	    mean_gene[i] += gene[i] / nsample;
		    for (j = 0; j < ngroup; j++)
		    {
			 mean_inter[i][j] += inter[i][j] / nsample;
			 mean_sigma2_bio[i][j] += sigma2_bio[i][j] / nsample;
                         expr_est[i][j] = mu + gene[i] + cell[j] + inter[i][j];
                         mean_expr_est[i][j] += expr_est[i][j] / nsample;

			 for (k = 0; k < group[j]; k++)
			 {
			      mean_expr[i][j][k] += expr[i][j][k] / nsample;
                      	      mean_sigma2_array[i][j][k] += sigma2_array[i][j][k] / nsample;
			 }
		   }
	       }
        } 
        /*printf("\n");*/




        /*Save output*/
        n = 0;
        m = 0;
	for (i = 0; i < ngene; i++)
	{
		for (j = 0; j < ngroup; j++)
  		{
		       mexprest[n] = mean_expr_est[i][j];
		       msigma2b[n] = mean_sigma2_bio[i][j];
	               n = n + 1;		
		       for (k = 0; k < group[j]; k++)
		       {
    		            mexpr[m] = mean_expr[i][j][k];
		            msigma2e[m] = mean_sigma2_array[i][j][k];
	                    m = m + 1;				
		       }
                }
		/*printf("==> %d \n ",i);*/
	}




        /*Compute F-scores*/
        m = 0;
	for (j = 0; j < ngroup; j++) m += group[j];

	for (i = 0; i < ngene; i++)
	{
		mean =0.0;
		for (j = 0; j < ngroup; j++) mean += mean_expr_est[i][j]/ngroup;    

	        fstat[i] = 0.0;
		for (j = 0; j < ngroup; j++)
		{ 
		     var = mean_sigma2_bio[i][j];
                     if(method_array==4) var += mean_sigma2_array[i][j][0];
                     else for (k = 0; k < group[j]; k++) var += mean_sigma2_array[i][j][k]/group[j];
                     fstat[i] += (group[j] * (mean_expr_est[i][j]-mean) * (mean_expr_est[i][j]-mean) / var);
                }
                fstat[i] = fstat[i]/m;
	}    


}



void onelayerhem(double *dat,  int *opt, int *ndat, int *grp, int *nMCMC, double *par,  int *Bsize, double *vart, double *fstat,  double *mexprest, double *msigma2, double *MCMCsamp) {

	int i, j, k, m, n;
        double mean, var;

        method_array = opt[0];
        method_bio   = opt[1];
        method_total = opt[2];

        ngene  = ndat[0];
        nchip  = ndat[1];
        ngroup = ndat[2];

        burnin  = nMCMC[0];
        nsample = nMCMC[1];

        sigma2g     = par[0];
        sigma2c     = par[1];
        sigma2r     = par[2];
        alpha_bio   = par[3];
        beta_bio    = par[4];
        max_sigma2_bio = 10;



	/* Read data */
	read_data_file_onelayer(dat, grp);

 
 	/* Empirical Bayes */
        if(method_total > 1) read_prior_file_onelayer(Bsize, vart);  


	/* Allocate memory */
	mean_gene = (double *) calloc(ngene, sizeof(double));
	mean_cell = (double *) calloc(ngroup, sizeof(double));
	mean_inter = (double **) calloc(ngene, sizeof(double));
	mean_sigma2_bio = (double **) calloc(ngene, sizeof(double));
	expr_est = (double **) calloc(ngene, sizeof(double));
	mean_expr_est = (double **) calloc(ngene, sizeof(double));

        for(i = 0; i < ngene; i++) 
        {
            mean_inter[i] = (double *) calloc(ngroup, sizeof(double)); 
            mean_sigma2_bio[i] = (double *) calloc(ngroup, sizeof(double)); 
            expr_est[i] = (double *) calloc(ngroup, sizeof(double)); 
            mean_expr_est[i] = (double *) calloc(ngroup, sizeof(double)); 
        }

	/* Set zero for output variables */
	mean_mu = 0.;
	for (j = 0; j < ngroup; j++) mean_cell[j] = 0.;

	for (i = 0; i < ngene; i++)
	{
		mean_gene[i] = 0.;
		for (j = 0; j < ngroup; j++)
		{
			mean_inter[i][j] = 0.;
			mean_sigma2_bio[i][j] = 0.;
  			mean_expr_est[i][j] = 0.;
		}
	}

	/* Set and print initial values*/
	initialize();

         /* Save MCMC samples*/
         m=0;i=1;j=1;
         MCMCsamp[m] = mu; m = m + 1;
         MCMCsamp[m] = gene[i]; m = m + 1;
         MCMCsamp[m] = cell[j]; m = m + 1;
         MCMCsamp[m] = inter[i][j]; m = m + 1;
         MCMCsamp[m] = sigma2_bio[i][j]; m = m + 1;

	
        /*** Burn-ins ***/
        k = 0;
	for (n = 0; n < burnin; n++)
	{
   	     update_expr_total();
 	     update_inter();
	     update_gene();
 	     update_cell();
 	     var = update_mu();
	     for (i = 0; i < ngene; i++){
		  for (j = 0; j < ngroup; j++){
                       expr_est[i][j] = mu + gene[i] + cell[j] + inter[i][j];
		  }
             }
 	     if(method_total==3) update_sigma2_total_nonpar();
 	     else update_sigma2_total();


             /* Save MCMC samples*/
             i=1;j=1;	     
             MCMCsamp[m] = mu; m = m + 1;
             MCMCsamp[m] = gene[i]; m = m + 1;
             MCMCsamp[m] = cell[j]; m = m + 1;
             MCMCsamp[m] = inter[i][j]; m = m + 1;
             MCMCsamp[m] = sigma2_bio[i][j]; m = m + 1;

             /*
             k +=1; 
             if(k == 50) {k = 0; printf(".");}
             */
        }


        /*** Keep MCMC samples ***/
	for (n = 0; n < nsample; n++)
	{
   	     update_expr_total();
 	     update_inter();
	     update_gene();
 	     update_cell();
 	     var = update_mu();

	     for (i = 0; i < ngene; i++){
		  for (j = 0; j < ngroup; j++){
                       expr_est[i][j] = mu + gene[i] + cell[j] + inter[i][j];
		  }
             }

 	     if(method_total==3) update_sigma2_total_nonpar();
 	     else update_sigma2_total();
                        
	     mean_mu += mu / nsample;
     	     for (j = 0; j < ngroup; j++) mean_cell[j] += cell[j] / nsample;

	     for (i = 0; i < ngene; i++)
	     {
		  mean_gene[i] += gene[i] / nsample;
		  for (j = 0; j < ngroup; j++)
		  {
		       mean_inter[i][j] += inter[i][j] / nsample;
		       mean_sigma2_bio[i][j] += sigma2_bio[i][j] / nsample;
                       mean_expr_est[i][j] += expr_est[i][j] / nsample;
		  }
	     }


             /* Save MCMC samples*/
             i=1;j=1;	     
             MCMCsamp[m] = mu; m = m + 1;
             MCMCsamp[m] = gene[i]; m = m + 1;
             MCMCsamp[m] = cell[j]; m = m + 1;
             MCMCsamp[m] = inter[i][j]; m = m + 1;
             MCMCsamp[m] = sigma2_bio[i][j]; m = m + 1;

             /*
             k +=1; 
             if(k == 50) {k = 0; printf(".");}
             */

        } 
        /*printf("\n");*/




        /*Save output*/
        n = 0;
        m = 0;
	for (i = 0; i < ngene; i++)
	{
		for (j = 0; j < ngroup; j++)
  		{
		       mexprest[n] = mean_expr_est[i][j];
		       msigma2[n] = mean_sigma2_bio[i][j];
	               n = n + 1;		
                }
		/*printf("==> %d \n ",i);*/
	}



        /*Compute F scores*/
        m = 0;
	for (j = 0; j < ngroup; j++) m += group[j];

	for (i = 0; i < ngene; i++)
	{
		mean =0.0;
		for (j = 0; j < ngroup; j++) mean += mean_expr_est[i][j]/ngroup;    

	        fstat[i] = 0.0;
		for (j = 0; j < ngroup; j++)
		{ 
		     var = mean_sigma2_bio[i][j];
                     fstat[i] += (group[j] * (mean_expr_est[i][j]-mean) * (mean_expr_est[i][j]-mean) / var);
                }
                fstat[i] = fstat[i]/m;
	}    
 
}

/* ****************************************************************** */
/*                                                                    */
/*    PRINTING FUNCTION                                               */
/*                                                                    */
/* ****************************************************************** */

void chk_prt()
{
  int i, j, k, l;
  FILE *fout;
  fout = fopen("tmp.out","w");

  i=1;
  j=1;
  k=1;
  l=1;

  printf("\n");

  
  /*  
  printf("obs         = %lf ",obs[i][j][k][l]);
  printf("expr        = %lf ",expr[i][j][k]);
  fprintf(fout, "inter       = %lf ",inter[i][j]);
  fprintf(fout, "gene        = %lf ",gene[i]);
  fprintf(fout, "cell        = %lf ",cell[j]);
  fprintf(fout, "mu          = %lf ",mu);
  fprintf(fout, "sigma2_array= %lf ",sigma2_array[i][j][k]);
  fprintf(fout, "sigma2_bio  = %lf ",sigma2_bio[i][j]);
  */

  fprintf(fout, "%lf, ",inter[i][j]);
  fprintf(fout, "%lf, ",gene[i]);
  fprintf(fout, "%lf, ",cell[j]);
  fprintf(fout, "%lf, ",mu);
  fprintf(fout, "%lf, ",sigma2_array[i][j][k]);
  fprintf(fout, "%lf ",sigma2_bio[i][j]);

  if(iii >10) fclose(fout);

}



/******** END OF PROGRAM ***********************************************************************/
