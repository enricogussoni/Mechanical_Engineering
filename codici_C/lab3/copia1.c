#include <stdio.h>

#define max 100
int main ()
{
                                             /*Copia stringa senza strcopy*/
char stringa1[max],stringa2[max];
int lung1,lung2,i;
                                                   /*Acquisizione stringa1*/
printf("Quanto � lunga stringa1?\n");
scanf("%d",lung1);

do
  {
   printf("La stringa deve contenere meno di 100 elementi, reinserire un valore<100?\n");
   scanf("%d",lung1);
  }while((lung1<0)||(lung1>100));

printf("S'inseriscano i caratteri di stringa1:\n");
scanf("%s",stringa1);
printf("CONTROLLO: stringa1 contiene:%s\n",stringa1);

                                                    /*Acquisizione stringa2*/
printf("Quanto � lunga stringa2?\n");
scanf("%d",lung2);

do
  {
   printf("La stringa deve contenere meno di 100 elementi, reinserire un valore<100?\n");
   scanf("%d",lung2);
  }while((lung2<0)||(lung2>100));

if(lung1>lung2)
   printf("Impossibile eseguire copia: la stinga 1 � pi� lunga della stringa 2\n");
else
    for(i=0;i=lung1;i++)
        {
         stringa2[i]=stringa1[i];
        }
printf("Ora stringa2 contiene: %s\n",stringa2);

return 0;
}
