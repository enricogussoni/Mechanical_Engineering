#include <stdio.h>

#define r 10  /*righe*/
#define c 10  /*colonne*/
              /*in questo caso matrice quadrata*/
int main ()
{
 int campo[r][c],appoggio[r][c];
 int n,x,y,z,p,q;                        /*Conteggi*/
 int k,intorno=0;  /*k cicli e somma dei "vicini"*/
 int ins1,ins2;/*Utile per termionare inserimento
                         o rilanciare programma*/
 int N=(r*c); /*Numero celle*/

                                        /*INTRO*/
printf("CONWAY'S GAME OF LIFE\n");
printf("versione del gioco con matrice %d x %d\n",r,c);
printf("\nInserimento dei valori:"
       "\ninserire le coordinate x,y delle celle vive:\n"
       "(inserire -1,-1 per concludere l'inserimento)\n");

for(x=0;x<r;x++)
   {
    for(y=0;y<c;y++)
       {
        campo[x][y]=0;
       }
	}

for(x=0;x<r;x++)
   {
    for(y=0;y<c;y++)
       {
        printf("%d ", campo[x][y]);
       }
		printf ("\n");
	}

	for(x=0;x<r;x++)
   {
    for(y=0;y<c;y++)
       {
        appoggio[x][y]=0;
       }
	}

for(x=0;x<r;x++)
   {
    for(y=0;y<c;y++)
       {
        printf("%d ", appoggio[x][y]);
       }
		printf ("\n");
	}

                     /*Inserimento celle "vive"*/
do{
    printf("ascissa:\n");
    scanf("%d",&ins1);
    printf("ordinata:\n");
    scanf("%d",&ins2);
    if(ins1!=-1 && ins2!= -1)
        campo[ins1][ins2]=1;
		x++;
	}while(x<N && (ins1!=-1 && ins2!=-1));

	                     /*Stampa campo iniziale*/
for(x=0;x<r;x++)
   {
    for(y=0;y<c;y++)
       {
        printf("%d ", campo[x][y]);
       }
		printf ("\n");
	}

/*Inizio computazione*/
printf("Inserire numero di cicli:\n");
scanf("%d",&k);

for(n=0;n<k;n++)
   {
    for(x=1;x<r-1;x++)                 /*Scorri righe matrice r-1,c-1*/
       {
        for(y=1;y<c-1;y++)          /*Scorri colonne matrice r-1,c-1**/
           {                        /*Somma delle caselle dell'intorno di x,y*/
            intorno=(intorno+
                    campo[x-1][y-1]+
                    campo[x][y-1]+
                    campo[x+1][y-1]+
                    campo[x-1][y]+
                    campo[x+1][y]+
                    campo[x-1][y+1]+
                    campo[x][y+1]+
                    campo[x+1][y+1]);

                                        /*Creazione nuovo campo in "appoggio"*/
            if((campo[x][y]==0)&&(intorno==3))
               appoggio[x][y]=1;
            else
                if((campo[x][y]==1)&&((intorno==2)||(intorno==3)))
                   appoggio[x][y]==1;
                else
                   appoggio[x][y]==0;
            }
            intorno=0;
        }

      /*riporto a 0 i casi limite*/

    for(x=0;x<c;x++)
       {appoggio[0][x]=0;
        appoggio[x][0]=0;
        appoggio[r][x]=0;
        appoggio[x][c]=0;}

    /*Copio "appoggio" in "campo"*
    for(x=0;x<r;x++)
       {
        for(y=0;y<c;y++)
            {campo[x][y]=appoggio[x][y];}
       }

    /*Stampo stato*/
    for(x=1;x=r-1;x++)
       {
        for(y=1;y<c-1;y++)
            {
             printf("%d",appoggio[x][y]);
            }
        printf("\n");
       }
    }

return 0;
}
