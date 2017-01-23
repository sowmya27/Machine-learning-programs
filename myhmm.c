//
//  main.c
//  myhmm
//
//  Created by Annapurna Annadatha on 9/13/15.
//  Copyright (c) 2015 Annapurna Annadatha. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N 2
#define M 27
#define ALPHA_LOWER {"abcdefghijklmnopqrstuvwxyz "}
#define ALPHA_UPPER {"ABCDEFGHIJKLMNOPQRSTUVWXYZ "}
#define MAX_CHARS 500


struct hmmStruct
{
    int obs;
    double c;
    double alpha[N];
    double beta[N];
    double gamma[N];
    double diGamma[N][N];
};


void cal_alpha(struct hmmStruct *pass_struct,double pi[],double A[N][N],double B[N][M],int T)
{
    int i,j,t;
    double temp;
    
    // compute alpha 0
    temp = 0.0;
    for(i = 0; i < N; ++i)
    {
        pass_struct[0].alpha[i] = pi[i] * B[i][pass_struct[0].obs];
        temp += pass_struct[0].alpha[i];
    }
    
    //Scale alpha 0
    pass_struct[0].c = 1.0 / temp;
    for(i = 0; i < N; ++i)
    {
        pass_struct[0].alpha[i] /= temp;
    }
    
    // compute alpha t
    for(t = 1; t < T; ++t)
    {
        temp = 0.0; //c0=0
        for(i = 0; i < N; ++i)
        {
            pass_struct[t].alpha[i] = 0.0;
            for(j = 0; j < N; ++j)
            {
                pass_struct[t].alpha[i] += pass_struct[t - 1].alpha[j] * A[j][i];
            }
            pass_struct[t].alpha[i] *= B[i][pass_struct[t].obs];
            temp += pass_struct[t].alpha[i];
        }
       
        // scale alpha t
        pass_struct[t].c = 1.0 / temp; // c0= 1/c0
        for(i = 0; i < N; ++i)
        {
            pass_struct[t].alpha[i] /= temp;
        }
        
    }
    
}


void cal_beta(struct hmmStruct *pass_struct,double pi[],double A[N][N],double B[N][M],int T)
{
    int i,
    j,
    t;
    
    // Let Beta T-1 (i) = 1 scaled by C(T-1)
    for(i = 0; i < N; ++i)
    {
        pass_struct[T - 1].beta[i] = pass_struct[T - 1].c;
    }
    
    // beta pass
    for(t = T - 2; t >= 0; --t)
    {
        for(i = 0; i < N; ++i)
        {
            pass_struct[t].beta[i] = 0.0;
            for(j = 0; j < N; ++j)
            {
                pass_struct[t].beta[i] += A[i][j] * B[j][pass_struct[t + 1].obs] * pass_struct[t + 1].beta[j];
            }
            
            // scale beta
            pass_struct[t].beta[i] *= pass_struct[t].c;
        }
        
    }
    
}


void cal_gamma(struct hmmStruct *pass_struct,double pi[],double A[][N],double B[][M],int T)
{
    int i,
    j,
    t;
    
    double denom;

    
    // compute gamma and digamma
    for(t = 0; t < T - 1; ++t)
    {
        denom = 0.0;
        for(i = 0; i < N; ++i)
        {
            for(j = 0; j < N; ++j)
            {
                denom += pass_struct[t].alpha[i] * A[i][j] * B[j][pass_struct[t + 1].obs] * pass_struct[t + 1].beta[j];
            }
        }

        
        for(i = 0; i < N; ++i)
        {
            pass_struct[t].gamma[i] = 0.0;
            for(j = 0; j < N; ++j)
            {
                pass_struct[t].diGamma[i][j] = (pass_struct[t].alpha[i] * A[i][j] * B[j][pass_struct[t + 1].obs] * pass_struct[t + 1].beta[j])/ denom;
                pass_struct[t].gamma[i] += pass_struct[t].diGamma[i][j];
            }
        }
        
    }
    
}


void cal_reestimatePi(struct hmmStruct *pass_struct,double pi_new[])
{
    int i;
    
    for(i = 0; i < N; ++i)
    {
        pi_new[i] = pass_struct[0].gamma[i];
    }
    
}


void cal_reestimateA(struct hmmStruct *pass_struct,double A_new[][N],int T)
{
    int i,
    j,
    t;
    
    double numer,
    denom;
    
    for(i = 0; i < N; ++i)
    {
        for(j = 0; j < N; ++j)
        {
            numer = denom = 0.0;
            
            for(t = 0; t < T - 1; ++t)
            {
                numer += pass_struct[t].diGamma[i][j];
                denom += pass_struct[t].gamma[i];
                
            }
            
            A_new[i][j] = numer / denom;
            
        }
        
    }
    
}



void cal_reestimateB(struct hmmStruct *pass_struct,double B_new[][M],int T)
{
    int i,
    j,
    t;
    
    double numer,
    denom;
    
    
    for(i = 0; i < N; ++i)
    {
        for(j = 0; j < M; ++j)
        {
            numer = denom = 0.0;
            for(t = 0; t < T - 1; ++t)
            {
                if(pass_struct[t].obs == j)
                {
                    numer += pass_struct[t].gamma[i];
                }
                denom += pass_struct[t].gamma[i];
                
            }
            
            B_new[i][j] = numer / denom;
            
        }
        
    }
    
}


void initialization(double pi[],double A[][N],double B[][M],int seed)
{
    int i,
    j;
    
    double prob,
    temp,
    temp2;
    
    // initialize pseudo-random number generator
    srand(seed);
    
    // initialize pi
    prob = 1.0 / (double)N;
    temp = prob / 10.0;
    temp2 = 0.0;
    for(i = 0; i < N; ++i)
    {
        if((rand() & 0x1) == 0)
        {
            pi[i] = prob + (double)(rand() & 0x7) / 8.0 * temp;
        }
        else
        {
            pi[i] = prob - (double)(rand() & 0x7) / 8.0 * temp;
        }
        temp2 += pi[i];
        
    }
    
    for(i = 0; i < N; ++i)
    {
        pi[i] /= temp2;
    }
    
    // initialize A[][]
    prob = 1.0 / (double)N;
    temp = prob / 10.0;
    for(i = 0; i < N; ++i)
    {
        temp2 = 0.0;
        for(j = 0; j < N; ++j)
        {
            if((rand() & 0x1) == 0)
            {
                A[i][j] = prob + (double)(rand() & 0x7) / 8.0 * temp;
            }
            else
            {
                A[i][j] = prob - (double)(rand() & 0x7) / 8.0 * temp;
            }
            temp2 += A[i][j];
            
        }// next j
        
        for(j = 0; j < N; ++j)
        {
            A[i][j] /= temp2;
        }
        
    }
    
    
    // initialize B[][]
    prob = 1.0 / (double)M;
    temp = prob / 10.0;
    for(i = 0; i < N; ++i)
    {
        temp2 = 0.0;
        for(j = 0; j < M; ++j)
        {
            if((rand() & 0x1) == 0)
            {
                B[i][j] = prob + (double)(rand() & 0x7) / 8.0 * temp;
            }
            else
            {
                B[i][j] = prob - (double)(rand() & 0x7) / 8.0 * temp;
            }
            temp2 += B[i][j];
            
        }
        
        for(j = 0; j < M; ++j)
        {
            B[i][j] /= temp2;
        }
        
    }
    
}



int Cal_T(char fname[],int startPos,int startChar,int maxChars)
{
    FILE *in;
    
    int i,
    j,
    len,
    thisStartPos,
    totalNum,
    num;
    
    char temp[MAX_CHARS + 1];
    
    char space[1] = {" "};
    
    char alphabetUpper[M] = ALPHA_UPPER;
    char alphabetLower[M] = ALPHA_LOWER;
    
    in = fopen(fname, "r");
    if(in == NULL)
    {
        fprintf(stderr, "\nError opening file %s\n\n", fname);
        exit(0);
    }
    

    // count
    totalNum = num = 0;
    while(fgets(temp, MAX_CHARS, in) != NULL)
    {
        len = strlen(temp);
        
        // each line should end with a single space
        while((strncmp(&temp[len - 1], space, 1) == 0) && (len > 0))
        {
            --len;
        }
        
        strncpy(&temp[len], space, 1);
        
        thisStartPos = startPos;
       
        while((strncmp(&temp[thisStartPos], space, 1) == 0) && (thisStartPos < len))
        {
            ++thisStartPos;
        }
        
        for(i = thisStartPos; i <= len; ++i)
        {
            
            for(j = 0; j < M; ++j)
            {
                if((strncmp(&temp[i], &alphabetLower[j], 1) == 0)
                   || (strncmp(&temp[i], &alphabetUpper[j], 1) == 0))
                {
                    ++totalNum;
                    if(totalNum >= startChar)
                    {
                        ++num;
                        if((maxChars > 0) && (num >= maxChars))
                        {
                            return(num);
                        }
                        
                    }
                    
                    break;
                    
                }
                
            }
            
        }
        
    }
    
    fclose(in);
    
    return(num);
    
}


int Cal_Observations(char fname[],struct hmmStruct *pass_struct,int T,int startPos,int startChar,int maxChars)
{
    FILE *in;
    
    int i,j,len,num,thisStartPos,totalNum;
    
    char temp[MAX_CHARS + 2];
    
    char space[1] = {" "};
    
    char alphabetUpper[M] = ALPHA_UPPER;
    char alphabetLower[M] = ALPHA_LOWER;
    
    in = fopen(fname, "r");
    if(in == NULL)
    {
        fprintf(stderr, "\nError opening file %s\n\n", fname);
        exit(0);
    }

    
    totalNum = num = 0;
    while(fgets(temp, MAX_CHARS, in) != NULL)
    {
        len = strlen(temp);
        
     
        while((strncmp(&temp[len - 1], space, 1) == 0) && (len > 0))
        {
            --len;
        }
        strncpy(&temp[len], space, 1);
        
        thisStartPos = startPos;
        
        while((strncmp(&temp[thisStartPos], space, 1) == 0) && (thisStartPos < len))
        {
            ++thisStartPos;
        }
        
        for(i = thisStartPos; i <= len; ++i)
        {
          
            for(j = 0; j < M; ++j)
            {
                if((strncmp(&temp[i], &alphabetLower[j], 1) == 0)
                   || (strncmp(&temp[i], &alphabetUpper[j], 1) == 0))
                    
                {
                    ++totalNum;
                    if(totalNum >= startChar)
                    {
                        pass_struct[num].obs = j;

                        ++num;
                        if(num > T)
                        {
                            printf("\nError in Cal_Observations()\n\n");
                            exit(0);
                        }
                        if((maxChars > 0) && (num >= maxChars))
                        {
                            return(num);
                        }
                        
                    }
                    
                    break;
                    
                }
            }
            
        }
        
    }
    
    fclose(in);
    
    return(num);
    
}


void printPi(double pi[])
{
    int i;
    double temp=0.0;
    for(i = 0; i < N; ++i)
    {
        printf("%8.5f ", pi[i]);
        temp += pi[i];
    }
   // printf(",  sum = %f\n", temp);
    
}


void printA(double A[][N])
{
    int i,j;
    double temp;
    
    for(i = 0; i < N; ++i)
    {
        temp = 0.0;
        for(j = 0; j < N; ++j)
        {
            printf("%8.5f ", A[i][j]);
            temp += A[i][j];
        }
      printf("\n");
        
    }
    
}


void printBT(double B[][M])
{
    int i,j;
    double temp;
    char alphabet[M] = ALPHA_LOWER;
    
    for(i = 0; i < M; ++i)
    {
        printf("%c ", alphabet[i]);
        for(j = 0; j < N; ++j)
        {
            printf("%8.5f ", B[j][i]);
        }
        printf("\n");
    }
    for(i = 0; i < N; ++i)
    {
        temp = 0.0;
        for(j = 0; j < M; ++j)
        {
            temp += B[i][j];
        }
       // printf("sum[%d] = %f ", i, temp);
    }
    printf("\n");
    
}


int main()
{
    int startPos,startChar,maxChars,maxIters,i,j,T,iters,seed;
    double oldLogProb,newLogProb;
    double pi[N],pi_new[N],A[N][N],A_new[N][N],B[N][M],B_new[N][M];
    
    char fname[50];
    
    struct hmmStruct *pass_struct;
    
    printf("enter file name:");
    gets(fname);
    printf("enter iterations:");
    scanf("%d",&maxIters);
    //printf("enter seed value: ");
    //scanf("%d",&seed);
    
    startPos = 15;
    startChar = 0;
    maxChars = 50000;
    seed = 1241;
    
    // number of observations
    
    fflush(stdout);
    T = Cal_T(fname,startPos,startChar,maxChars);
    
    // memory
    fflush(stdout);
    if((pass_struct = calloc(T + 1, sizeof(struct hmmStruct))) == NULL)
    {
        fprintf(stderr, "\nUnable to allocate alpha\n\n");
        exit(0);
    }
       // read in the observations from file
    fflush(stdout);
    T = Cal_Observations(fname,pass_struct,T,startPos,startChar,maxChars);
    
    srand(seed);
    
    initialization(pi, A, B, seed);
    
    printf("\nN = %d \nM = %d \nT = %d\n", N, M, T);
    printf("\nInitial Pi Matrix =\n");
    printPi(pi);
    printf("\nInitial A  Matrix=\n");
    printA(A);
    printf("\nInitial B Matrix Transpose =\n");
    printBT(B);
 
    iters = 0;
    oldLogProb = -1.0;
    newLogProb = 0.0;
    
    while((iters < maxIters) && (newLogProb > oldLogProb))
    {
        printf("\n Iteration = %d\n", iters);
        
        oldLogProb = newLogProb;
     
        fflush(stdout);
        cal_alpha(pass_struct, pi, A, B, T);
       
        fflush(stdout);
        cal_beta(pass_struct, pi, A, B, T);
       
        fflush(stdout);
        cal_gamma(pass_struct, pi, A, B, T);
        
        fflush(stdout);
        cal_reestimatePi(pass_struct, pi_new);
        
        fflush(stdout);
        cal_reestimateA(pass_struct, A_new, T);
        
        fflush(stdout);
        cal_reestimateB(pass_struct, B_new, T);
        

        printf("\npi_new  =\n");
        printPi(pi_new);
        printf("\nA_new =\n");
        printA(A_new);
        printf("\nB_new^T = \n");
        printBT(B_new);
        
        
        for(i = 0; i < N; ++i)
        {
            pi[i] = pi_new[i];
            
            for(j = 0; j < N; ++j)
            {
                A[i][j] = A_new[i][j];
            }
            
            for(j = 0; j < M; ++j)
            {
                B[i][j] = B_new[i][j];
            }
            
        }
        
        // compute log
        newLogProb = 0.0;
        for(i = 0; i < T; ++i)
        {
            newLogProb += log(pass_struct[i].c);
        }
        newLogProb = -newLogProb;
        
        
        if(iters == 0)
        {
            oldLogProb = newLogProb - 1.0;
        }
        
        printf("Iteration = %d, \nLog Probability = %f\n",iters, newLogProb);
        
        ++iters;
        
    }
    
    
    printf("Final Pi Matrix =\n");
    printPi(pi);
    printf("\nFinal A Matrix =\n");
    printA(A);
    printf("\nFinal B Matrix Transpose =\n");
    printBT(B);
    printf("\nT = %d \nN = %d \nM = %d \nNum of iterations = %d", T, N, M, iters);
    printf("\nLog Probability = %f\n\n", newLogProb);
    
}



