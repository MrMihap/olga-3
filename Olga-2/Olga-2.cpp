#include "stdafx.h"

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#define MATH_PI 3.14159265358979323846
// ******************
// MAIN
int n = 16;
double d1 = 1, d2 = 0, d3 = 1, d4 = 0, d5 = 0;
//int NOT_VOID_ELEMENTS = 7 * n * n + 6 * n + 1;
int FuncType = 1;
double Eps = 1.e-4;
double AlphaAngle = 0;
// ********************

const int MAX_WAIT_RANGE = 150;
double h;

double Fxy(double x, double y)
{
	//SafeMod

	switch(FuncType)
	{
	case 1:
		return 1;
	case 2:
		return d1 * x + d2 * y;
	case 3:
		return d3 * x * x + d4 * y * y + d5 * x * y;
	}
	return 0;
}

double FF(void)

{


	switch(FuncType)
	{
                    
        case 1:
		return sin(AlphaAngle);
		
		
	case 2:
		return d1 * d1 * (sin(AlphaAngle) / 3.0 
        + cos(AlphaAngle) * sin(AlphaAngle) / 2.0
        + cos(AlphaAngle) * cos(AlphaAngle) * sin(AlphaAngle) / 3.0 ) 
        + d1 * d2 * (sin(AlphaAngle) * sin(AlphaAngle) / 2.0 
        + 2.0 * cos(AlphaAngle) * sin(AlphaAngle) * sin(AlphaAngle) / 3.0)
        + d2 * d2 * sin(AlphaAngle) * sin(AlphaAngle) * sin(AlphaAngle) / 3.0;
                
		
	case 3:
		return d3 * d3 * sin(AlphaAngle) / 5.0 
        + ( 2.0 * d3 * cos(AlphaAngle) / sin(AlphaAngle) + d4 ) 
        * ( 2.0 * d3 * cos(AlphaAngle) / sin(AlphaAngle) + d4 )
        * sin(AlphaAngle) * sin(AlphaAngle) * sin(AlphaAngle) / 9.0
        + ( d3 * cos(AlphaAngle) * cos(AlphaAngle) / ( sin(AlphaAngle) * sin(AlphaAngle ) ) 
        + d4 * cos(AlphaAngle) / sin(AlphaAngle) + d5 )
        * ( d3 * cos(AlphaAngle) * cos(AlphaAngle) / ( sin(AlphaAngle) * sin(AlphaAngle ) ) 
        + d4 * cos(AlphaAngle) / sin(AlphaAngle) + d5 )
        * sin(AlphaAngle) * sin(AlphaAngle) * sin(AlphaAngle) * sin(AlphaAngle) * sin(AlphaAngle) / 5.0
        + d3 * ( 2.0 * d3 * cos(AlphaAngle) / sin(AlphaAngle) 
	+ d4 ) * sin(AlphaAngle) * sin(AlphaAngle) / 4.0
        + 2.0 * d3 *( d3 * cos(AlphaAngle) * cos(AlphaAngle) / ( sin(AlphaAngle) * sin(AlphaAngle) ) 
        + d4 * cos(AlphaAngle) / sin(AlphaAngle) + d5 )
        * sin(AlphaAngle) * sin(AlphaAngle) * sin(AlphaAngle) / 9.0
        + ( 2.0 * d3 * cos(AlphaAngle) / sin(AlphaAngle) + d4 )
        * ( d3 * cos(AlphaAngle) * cos(AlphaAngle) / ( sin(AlphaAngle) * sin(AlphaAngle) )
        + d4 * cos(AlphaAngle) / sin(AlphaAngle) + d5 )
        * sin(AlphaAngle) * sin(AlphaAngle) * sin(AlphaAngle) * sin(AlphaAngle) / 4.0;   
                
		

	}
	return 0;
}
/*{
	double Ax = 0.5;
	double Bx = 1.0 + 0.5 * cos(AlphaAngle);
	double Cx = 0.5 + 1.0 * cos(AlphaAngle);
	double Dx = 0.0 + 0.5 * cos(AlphaAngle);
	double Ex = 0.5 + 0.5 * cos(AlphaAngle);


	double Ay = 0.0;
	double By = 0.5 * sin(AlphaAngle);
	double Cy = 1.0 * sin(AlphaAngle);
	double Dy = 0.5 * sin(AlphaAngle);
	double Ey = 0.5 * sin(AlphaAngle);

	double result = 0;

	result +=pow(Fxy(Ax, Ay), 2) * 1;
	result +=pow(Fxy(Bx, By), 2) * 1;
	result +=pow(Fxy(Cx, Cy), 2) * 1;
	result +=pow(Fxy(Dx, Dy), 2) * 1;
	result +=pow(Fxy(Ex, Ey), 2) * 2;

	result *= sin(AlphaAngle) / 6;
	return result;
}*/

double CalcB(double x, double y, int type)
{
	double result = 0;
	double Angle = AlphaAngle;
	switch(type)
	{
	case 9:
		result += Fxy(x + h * cos(Angle)		/ 2, y + h * sin(Angle) / 2); //A
		result += Fxy(x - h * (1 + cos(Angle))	/ 2, y + h * sin(Angle) / 2); //B
		result += Fxy(x - (h)                   / 2, y                     ); //C
		result += Fxy(x - h * cos(Angle)		/ 2, y + h * sin(Angle) / 2); //D
		result += Fxy(x + h * (1 + cos(Angle))	/ 2, y - h * sin(Angle) / 2); //E
		result += Fxy(x + h              		/ 2, y                     ); //F
		result *= h * h * sin(Angle) / 6;
		break;
	case 8:
		result += Fxy(x + h * cos(Angle)		/ 2, y + h * sin(Angle) / 2) * 2; //A
		result += Fxy(x - h*(1 + cos(Angle))	/ 2, y + h * sin(Angle) / 2) * 2; //B
		result += Fxy(x - (h)                   / 2, y                     ) * 1; //C
		result += Fxy(x + h              		/ 2, y                     ) * 1; //F
		result *= h * h * sin(Angle) / 12;
		break;
	case 7:
		result += Fxy(x + h * cos(Angle)		/ 2, y + h * sin(Angle) / 2) * 1; //A
		result += Fxy(x - h * (1 + cos(Angle))	/ 2, y + h * sin(Angle) / 2) * 2; //B
		result += Fxy(x - (h)                   / 2, y                     ) * 1; //C
		result *= h * h * sin(Angle) / 12;
		break;
	case 6:
		result += Fxy(x + h * cos(Angle)		/ 2, y + h * sin(Angle) / 2); //A
		result += 2 * Fxy(x - h * (1 + cos(Angle))	/ 2, y + h * sin(Angle) / 2); //B
		result += 2 *Fxy(x - (h)                   / 2, y                     ); //C
		result += Fxy(x - h * cos(Angle)		/ 2, y + h * sin(Angle) / 2); //D
		result *= h * h * sin(Angle) / 12;
		break;
	case 5:
		result += Fxy(x - (h)                   / 2, y                     ); //C
		result += Fxy(x - h * cos(Angle)		/ 2, y + h * sin(Angle) / 2); //D
		result *= h * h * sin(Angle) / 12;
		break;
	case 4:
		result += Fxy(x - (h)                   / 2, y                     ) * 1; //C
		result += Fxy(x - h * cos(Angle)		/ 2, y + h * sin(Angle) / 2) * 2; //D
		result += Fxy(x + h * (1 + cos(Angle))	/ 2, y - h * sin(Angle) / 2) * 2; //E
		result += Fxy(x + h              		/ 2, y                     ) * 1; //F
		result *= h * h * sin(Angle) / 12;
		break;
	case 3:
		result += Fxy(x - h * cos(Angle)		/ 2, y + h * sin(Angle) / 2) * 1; //D
		result += Fxy(x + h * (1 + cos(Angle))	/ 2, y - h * sin(Angle) / 2) * 2; //E
		result += Fxy(x + h              		/ 2, y                     ) * 1; //F
		result *= h * h * sin(Angle) / 12;

		break;
	case 2:
		result += Fxy(x + h * cos(Angle)		/ 2, y + h * sin(Angle) / 2) * 1; //A
		result += Fxy(x - h * cos(Angle)		/ 2, y + h * sin(Angle) / 2) * 1; //D
		result += Fxy(x + h * (1 + cos(Angle))	/ 2, y - h * sin(Angle) / 2) * 2; //E
		result += Fxy(x + h              		/ 2, y                     ) * 2; //F
		result *= h * h * sin(Angle) / 12;
		break;
	case 1:
		result += Fxy(x + h * cos(Angle)		/ 2, y + h * sin(Angle) / 2) * 1; //A
		result += Fxy(x + h              		/ 2, y                     ) * 1; //F
		result *= h * h * sin(Angle) / 12;
		break;
	}
	return result;
}

double Scalar(double *element_1, double *element_2, int Size)
{
	int i;
	double result = 0;
	for(i = 0; i < Size; i++)
	{
		result += element_1[i] * element_2[i];
	}
	return result;
}

double NormaVL2(double *element, int Size)
{
	int i;
	double result = 0;
	for(i = 0; i < Size; i++)
	{
		result += element[i] * element[i];
	}
	return sqrt(result);
}

double* Summ(double *a, double *b, int Size, double c1, double c2)
{
	int i;
	double *result = (double *)malloc(sizeof(double) * Size);
	for(i = 0; i < Size; i++)
	{
		result[i] = a[i] * c1 + b[i] * c2;
	}
	return result;
}

double* AlterMatrixVector(double* Matrix, int *IndexArray, int IndexSize, double *Vector, int VectorSize)
{
	int k;
	int NE = IndexSize;
	double *target = (double *) malloc (VectorSize * sizeof(double));
	for (k = 0; k < VectorSize; k++)
	{
		target[k] = 0.0;
	}

	for(k = 0; k < NE; k++)
	{
		target[IndexArray[k + IndexSize]] += Matrix[k] * Vector[IndexArray[k]];
	}
	return target;
}


int main()
{
	int i,j, k = 0;
	int N = (n + 1) * (n + 1); //Full Matrix Size
	int NOT_VOID_ELEMENTS = ( 7 * n * n + 6 * n + 1);
	double *c[MAX_WAIT_RANGE]; //Vector Array
	double *r[MAX_WAIT_RANGE]; //vector Array
	double *p[MAX_WAIT_RANGE]; //Vector Array
	double *V[MAX_WAIT_RANGE]; //vector Array
	double *Swap =			(double *)malloc(N * sizeof(double)); // Vector
	double *w =			(double *)malloc(N * sizeof(double)); // vector
	double *B =			(double *)malloc(N * sizeof(double)); // vector
	double *Betha = 	(double *)malloc(MAX_WAIT_RANGE * sizeof(double));  // unit
	double *Etha =		(double *)malloc(MAX_WAIT_RANGE * sizeof(double));  // unit
	double *Alpha = 	(double *)malloc(MAX_WAIT_RANGE * sizeof(double));  // unit
	double *lambda =	(double *)malloc(MAX_WAIT_RANGE * sizeof(double));  // unit
	double *Ksi =	  	(double *)malloc(MAX_WAIT_RANGE * sizeof(double));  // unit
	double Delta = 0;

	int *INDEX = (int *) malloc(NOT_VOID_ELEMENTS * 2 * sizeof(int));
	double *VALUE = (double *) malloc(NOT_VOID_ELEMENTS * sizeof(double));

	//Allocate memory for Vector Arrays
	for(i = 0; i < MAX_WAIT_RANGE; i++)
	{
		c[i] = (double *)malloc(N * sizeof(double));
		r[i] = (double *)malloc(N * sizeof(double));
		p[i] = (double *)malloc(N * sizeof(double));
		V[i] = (double *)malloc(N * sizeof(double));
	}
	k = 0;
	//prepearing B vector array
	h = 1/(float)n;

	// Enter alpha value in degrees
	printf("Enter alpha value in degrees \n");
	printf("\n");
	scanf("%lf", &AlphaAngle);
	AlphaAngle *= MATH_PI / 180;

	//Prepearing A Matrix as Arrays of indexes and values  of non zero elements
	for (i = 0; i < n + 1; i++)
	{
		for (j = 0; j < n + 1; j++)
		{
			int ii, jj;

			if( i == j && i >= 1 && i < n + 1 - 1)
			{
				// M matrix Patter
				for(ii = 0; ii < n + 1; ii++)
				{
					for(jj = 0; jj < n + 1; jj++)
					{
						if(ii == jj + 1 || ii == jj - 1)
						{
							INDEX[k] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 2;
							k++;
						}
						if( ii == jj && ii > 0 && ii < n + 1 - 1)
						{
							INDEX[k] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 12;
							k++;
						}
						if( ii == jj && ii == 0)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 6;
							k++;
						}
						if( ii == jj && ii == n)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 6;
							k++;
						}
					}
				}
			}
			if(i == j - 1 && j >= 1 && j < n + 1)
			{
				//L matrix pattern
				for(ii = 0; ii < n + 1; ii++)
				{
					for(jj = 0; jj < n + 1; jj++)
					{
						if(ii == jj + 1)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 2;
							k++;
						}
						if( ii == jj && ii > 0 && ii < n + 1 - 1)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 2;
							k++;
						}
						if( ii == jj && ii == 0)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 1;
							k++;
						}
						if( ii == jj && ii == n)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 1;
							k++;
						}
					}
				}

			}
			if(i == j + 1 && j >= 0 && j < n + 1 - 1)
			{
				//lt matrix pattern
				for(ii = 0; ii < n + 1; ii++)
				{
					for(jj = 0; jj < n + 1; jj++)
					{
						if(ii == jj - 1)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 2;
							k++;
						}
						if( ii == jj && ii > 0 && ii < n + 1 - 1)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 2;
							k++;
						}
						if( ii == jj && ii == 0)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 1;
							k++;
						}
						if( ii == jj && ii == n)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 1;
							k++;
						}
					}
				}
			}
			if(i == j && i == 0)
			{
				// K matrix Pattern
				for(ii = 0; ii < n + 1; ii++)
				{
					for(jj = 0; jj < n + 1; jj++)
					{
						if(ii == jj + 1 || ii == jj - 1)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 1;
							k++;
						}
						if( ii == jj && ii > 0 && ii < n + 1 - 1)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 6;
							k++;
						}
						if( ii == jj && ii == 0)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 2;
							k++;
						}
						if( ii == jj && ii == n)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 4;
							k++;
						}
					}
				}
			}
			if(i == j && i == n)
			{
				// N matrix Pattern
				for(ii = 0; ii < n + 1; ii++)
				{
					for(jj = 0; jj < n + 1; jj++)
					{
						if(ii == jj + 1 || ii == jj - 1)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 1;
							k++;
						}
						if( ii == jj && ii > 0 && ii < n + 1 - 1)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 6;
							k++;
						}
						if( ii == jj && ii == 0)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 4;
							k++;
						}
						if( ii == jj && ii == n)
						{
							INDEX[k * 1] = i * (n + 1) + ii;
							INDEX[k + NOT_VOID_ELEMENTS] = j * (n + 1) + jj;
							VALUE[k] = 2;
							k++;
						}
					}
				}
			}
		}
	}
	for(k = 0; k < NOT_VOID_ELEMENTS; k++)
	{
		//printf("Value[%d][%d] = %lf \n", INDEX[k * 1] + 1, INDEX[k + NOT_VOID_ELEMENTS] + 1, VALUE[k]);
		VALUE[k] *= h * h * sin(AlphaAngle)/24;
	}
	//printf("Matrix succsesful calced \n");
	//system("pause");
	for(FuncType = 1; FuncType <= 3; FuncType++)
	{

		for(i = 0; i < n + 1; i++)
		{
			for(j = 0; j < n + 1; j++)

			{
				if(i == j && j == 0)
					B[(n + 1) * i + j] = CalcB(i * h + j * h * cos(AlphaAngle), j * h * sin (AlphaAngle), 1);

				if(j > 0 && j < n + 1 - 1 && i == 0)
					B[(n + 1) * i + j] = CalcB(i * h + j * h * cos(AlphaAngle), j * h * sin (AlphaAngle), 2);

				if(j == n + 1 - 1 && i == 0)
					B[(n + 1) * i + j] = CalcB(i * h + j * h * cos(AlphaAngle), j * h * sin (AlphaAngle), 3);

				if(j == n + 1 - 1 && i > 0 && i < n + 1 - 1)
					B[(n + 1) * i + j] = CalcB(i * h + j * h * cos(AlphaAngle), j * h * sin (AlphaAngle), 4);

				if(j == n + 1 - 1 &&  i == n + 1 - 1)
					B[(n + 1) * i + j] = CalcB(i * h + j * h * cos(AlphaAngle), j * h * sin (AlphaAngle), 5);

				if(j > 0 && j < n + 1 - 1 && i == n + 1 - 1)
					B[(n + 1) * i + j] = CalcB(i * h + j * h * cos(AlphaAngle), j * h * sin (AlphaAngle), 6);

				if(j == 0 && i == n + 1 - 1)
					B[(n + 1) * i + j] = CalcB(i * h + j * h * cos(AlphaAngle), j * h * sin (AlphaAngle), 7);

				if(j == 0 && i > 0 && i < n + 1 - 1)
					B[(n + 1) * i + j] = CalcB(i * h + j * h * cos(AlphaAngle), j * h * sin (AlphaAngle), 8);

				if(j > 0 && j < n + 1 - 1 && i > 0 && i < n + 1 - 1)
					B[(n + 1) * i + j] = CalcB(i * h + j * h * cos(AlphaAngle), j * h * sin (AlphaAngle), 9);
			}
		}
		for(i = 0; i < N; i++)
		{
			//printf("B[%d] = %lf\n", i, B[i]);

		}
		//printf("B succsesful calced \n");
		Eps = 1.e-4;
		// обнулим массивы и массивы векторов
		for(i = 0; i < MAX_WAIT_RANGE; i++)
		{
			c[i] = (double *) malloc(N * sizeof(double));
			r[i] = (double *) malloc(N * sizeof(double));
			p[i] = (double *) malloc(N * sizeof(double));
			V[i] = (double *) malloc(N * sizeof(double));
		}
		//printf("CRPV succsesful Created \n");
		for(j = 0; j < MAX_WAIT_RANGE; j++)
		{
			for(i = 0; i < N; i++)
			{
				c[j][i] = 0;
				r[j][i] = 0;
				p[j][i] = 0;
				V[j][i] = 0;
			}
		}
		for(i = 0; i < MAX_WAIT_RANGE; i++)
		{
			lambda[i] = 0;
			Betha[i] = 0;
			Etha[i] = 0;
			Alpha[i] = 0;
			Ksi[i] = 0;
			w[i] = 0;
		}
		//printf("All arrays succsesful became as zero\n");
		//for( i = 0; i < (n + 1) * (n + 1); i++)
		//{
		//	//printf(" %lf ", b[i]);
		//}

		//****************************************
		// Главный цикл
		//****************************************

		// подготовка

		//Ro = b - ACo
		Swap = AlterMatrixVector(VALUE, INDEX, NOT_VOID_ELEMENTS , c[0], N );
		r[0] = Summ (B, Swap, N, 1, -1);
		free(Swap);
		/*for( j = 0; j < N; j++)
		{
			printf("r[0][%d] =  %lf \n", j, r[0][j]);
		}*/
		Ksi[1] = Betha[0] = NormaVL2(r[0], N);
		//printf("Ro succsesful calced \n");
		//		printf("\n\nV1:");
		for( j = 0; j < N; j++)
		{
			V[1][j] = r[0][j]/Betha[0];
			//printf("V[1][%d] =  %lf \n", j, V[1][j]);
		}
		//system("pause");
		//printf("V[1][j] succsesful calced \n");
		//		printf("\nRo:");
		//for( j = 0; j < (n + 1) * (n + 1); j++)
		//{
		//	//printf(" %lf ", r[0][j]);
		//}
		//		printf("\n");

		lambda[1] = 0;
		Betha[1] = 0;

		// тело цикла
		for( j = 1; j < MAX_WAIT_RANGE - 1; j++)
		{
			////1:  w = AVj - BETHAj * Vj-1
			Swap = AlterMatrixVector(VALUE,INDEX, NOT_VOID_ELEMENTS, V[j], N);
			w = Summ(Swap, V[j - 1], N, 1, -1 * Betha[j]);
			free(Swap);

			//2: ALPHAj = (w, Vj)
			Alpha[j] = Scalar(w, V[j], N);
			//system("pause");

			if(j > 1)
			{
				lambda[j] = Betha[j] / Etha[j - 1];
				Ksi[j] = -lambda[j] * Ksi[j - 1];
			}
			//5: ETHAj = ALPHAj - LAMBDAj * BETHAj
			Etha[j] = Alpha[j] - lambda[j] * Betha[j];

			//6: Pj = (Vj - BETHAj * Pj-1) / ETHAj
			p[j] = Summ(V[j], p[j - 1], N, 1/Etha[j], -Betha[j]/Etha[j]);

			//7: Cj = Cj-1 + KSIj * Pj
			c[j] = Summ(c[j - 1] , p[j], N, 1, Ksi[j]);

			//8: w = w - ALPHAj * Vj
			w = Summ(w, V[j], N, 1, - Alpha[j]);

			//9: BETHAj+1 = ||w||
			Betha[j + 1] = NormaVL2(w, N);

			//10: Vj+1 = w/BETHAj+1
			for(i = 0; i < N; i++)
			{
				V[j + 1][i] = w[i]/Betha[j + 1];
			}
			//printf("\n");

			//11 Ri = B - Ac[j]
			Swap = AlterMatrixVector(VALUE, INDEX, NOT_VOID_ELEMENTS, c[j], N);
			r[j] = Summ(B, Swap, N, 1, -1);
			//printf(" R[%d]/R[0] = %.3e \n", j, NormaVL2(r[j], N)/ NormaVL2(r[0], N) );
			for(i = 0; i < N; i++)
			{
				//printf("r[%d][%d] = %.3e \n", j, i, r[j][i]);
			}
			//system("pause");
			free(Swap);


			if(NormaVL2(r[j], N)/ NormaVL2(r[0], N) < Eps)
			{
				printf("\n\n\n");
				printf("FuncType = %d", FuncType);
				printf(" Step Num = %d ", j);

				printf(" R[%d]/R[0] = %.3e \n", j - 1, NormaVL2(r[j - 1], N)/ NormaVL2(r[0], N) );//?
				printf(" R[%d]/R[0] = %.3e \n", j, NormaVL2(r[j], N)/ NormaVL2(r[0], N) );//?

				Delta = fabs( FF() - 2 * Scalar(B, c[j], N) + Scalar(AlterMatrixVector(VALUE, INDEX, NOT_VOID_ELEMENTS, c[j], N), c[j], N));

				//BUG, JUST TEST!
				//Delta = fabs(Scalar(B, c[j], N) - Scalar(AlterMatrixVector(VALUE, INDEX, NOT_VOID_ELEMENTS, c[j], N), c[j], N));

				printf("DELTA =  %.3e (Eps = %.1e) \n", sqrt(Delta), Eps);
				if(Eps < 1.e-5)
				{
					printf("\n");
					break;
				}
				else
				{
					Eps = 1.e-6;
					if(NormaVL2(r[j], N)/ NormaVL2(r[0], N) < Eps)
					{

						printf("FuncType = %d", FuncType);
						printf(" Step Num = %d ", j);
						printf(" R[%d]/R[0] = %.3e \n", j - 1, NormaVL2(r[j - 1], N)/ NormaVL2(r[0], N) );//?
						printf(" R[%d]/R[0] = %.3e \n", j, NormaVL2(r[j], N)/ NormaVL2(r[0], N) );//?
						Delta = FF() - 2 * Scalar(B, c[j], N) + Scalar(AlterMatrixVector(VALUE, INDEX, NOT_VOID_ELEMENTS, c[j], N), c[j], N);
						if (fabs(Delta) < 1e-15)
						{Delta = 0.0;}
						//printf("DELTA1 =  %.3e (Eps = %.1e) \n", Delta, Eps);
						printf("DELTA =  %.3e (Eps = %.1e) \n", sqrt(Delta), Eps);
						printf("\n");
						break;

					}
					continue;
				}
			}

			//printf("\n r[%d]", j);
			for(i = 0; i < N; i++)
			{
				V[j + 1][i] = w[i]/Betha[j + 1];
				//printf(" %lf ", r[j][i]);
			}
		}
	}
	//system("pause");
	return 0;
}



