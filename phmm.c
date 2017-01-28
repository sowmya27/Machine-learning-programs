//
//  PHMM code phmm.c
//  prob5a
//
//  Created by Annapurna Annadatha on 9/22/15.
//  Copyright (c) 2015 Annapurna Annadatha. All rights reserved.
//

#include <stdio.h>
#include <math.h>

double MX[4][4]= {
    {0.0,0.0,0.0,0.0},
    {0.0,0.625,0.125,0.125},
    {0.0,0.111,0.111,0.333},
    {0.0,0.125,0.125,0.625}
};


double IX[4][4] = {
    {0.0,0.0,0.0,0.0},
    {0.0,0.25,0.25,0.25},
    {0.0,0.333,0.333,0.222},
    {0.0,0.25,0.25,0.25}
};

double MM[4] = {0.625,0.714,0.25,0.833};

double MI[4] = {0.125,0.143,0.5,0.167};

double MD[4] = {0.25,0.143,0.25,0.0};

double IM[4] = {0.333,0.333,0.5,0.5};

double II[4] = {0.333,0.333,0.375,0.5};

double ID[4] = {0.333,0.333,0.125,0.0};

double DM[4] = {0.0,0.5,0.333,0.66};

double DI[4] = {0,0.25,0.333,0.333};

double DD[4] = {0.0,0.25,0.333,0.0};



double fun_M(int i, int j);
double fun_I(int i, int j);
double fun_D(int i, int j);


//FjM(i)
double fun_M(int i, int j)
{
    double result =0.0;

    //Base case F0M(0)
    if (i==0 && j==0 )
        result =0.0;

    // excluding undefined terms
    else if (j==0 && i>0)
    {
        result =  log( (MI[j] * exp(fun_M(i-1,j))) + (II[j] * exp(fun_I(i-1,j))));
    }

    // excluding undefined terms
    else if (i==1 || j==1 || (i==0 && j>0) )
    {
        result = log(MX[j][i]/0.333);
    }

    //recursive relation
    else
    {
    result = log(MX[j][i]/0.333) + log( (MM[j-1] * exp(fun_M(i-1,j-1))) + (IM[j-1] * exp(fun_I(i-1,j-1))) + (DM[j-1] * exp(fun_D(i-1,j-1))));

    }

    printf("\n FM_%d(%d) = %f",i,j,result);
    return result;
}



//FjI(i)
double fun_I(int i, int j)
{
    double result =0.0;


    if (i==0 && j==0 )
        result =0.0;

    // excluding undefined terms
    else if (j==0 && i>0)
    {
        result =  log( (MI[j] * exp(fun_M(i-1,j))) + (II[j] * exp(fun_I(i-1,j))));
    }

    // excluding undefined terms
    else if (i==1 || j==1 || (i==0 && j>0) || (j==0 && i>0) )
    {
        result = log(IX[j][i]/0.333);
    }

    //recursive relation
    else
    {
        result = log(IX[j][i]/0.333) + log( (MI[j] * exp(fun_M(i-1,j))) + (II[j] * exp(fun_I(i-1,j))) + (DI[j] * exp(fun_D(i-1,j))));
    }

    printf("\n FI_%d(%d) = %f",i,j,result);

    return result;
}


//FjDi
double fun_D(int i, int j)
{
    double result =0.0;


    if (i==0 && j==0 )
        result =0.0;

    // excluding undefined terms
    else if ( (j==1)&& (i>=0) )
    {
        result = log( (MD[j-1] * exp(fun_M(i,j-1))) + (ID[j-1] * exp(fun_I(i,j-1))));
    }

    // excluding undefined terms
    else if ( (i==0 && j>0) || (j==0 && i>0) )
    {
        result = 0.0;
    }

    // excluding undefined terms
    else
    {
        result = log( (MD[j-1] * exp(fun_M(i,j-1))) + (ID[j-1] * exp(fun_I(i,j-1))) + (DD[j] * exp(fun_D(i,j-1))));
    }


    printf("\n FD_%d(%d) = %f",i,j,result);

    return result;
}


int main(int argc, const char * argv[]) {

    double score =0.0;

    //computing score
    score = log( (MM[3] * exp(fun_M(3,3))) + (IM[3] * exp(fun_I(3,3))) + (DM[3] * exp(fun_D(3,3))) );

    printf("\n\n\n PHMM Score = %f\n\n", score);


    return 0;
}
