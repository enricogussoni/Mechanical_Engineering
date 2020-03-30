#include <stdio.h>

#define dim 4
int main()
{
    float A[dim][dim],At[dim][dim];
    int x,y;

                                        /*Inserimento matrici*/
 for(x=0;x<dim;x++)
    {
     for(y=0;y<dim;y++)
         {
          printf("Inserire l'elemento (%d,%d):\n",x,y);
          scanf("%f",&A[x][y]);
         }
    }

                                                /*Inverto A*/
 for(y=0;y<dim;y++)
    {
     for(x=0;x<4;x++)
        {
         At[y][x]=A[x][y];
        }
    }

                /*Stampa (per ogni riga y tutte le colonne x) */
 for(y=0;y<dim;y++)
    {
     for(x=0;x<dim;x++)
         printf("%f",At[y][x]);
     printf("\n");
    }

 return 0;
}
