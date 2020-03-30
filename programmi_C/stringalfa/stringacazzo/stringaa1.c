#include <stdio.h>
#include <string.h>

#define dim 20
int main()
{
int c=0;
char stringa1[dim],stringa2[dim];

printf("Inserire la stringa 1:\n");
scanf("%s",stringa1);
printf("Inserire stringa2:\n");
scanf("%s",stringa2);

do{
   printf("CONTROLLO: c=%d, stringa1[%d]=%c e stringa2[%d]=%c\n",c,c,stringa1[c],c,stringa2[c]);
   c++;
  }while(stringa1[c]==stringa2[c]);

if(stringa1[c]>stringa2[c])
   printf("%s precede %s in ordine alfabetico\n",stringa2,stringa1);
else
   printf("%s precede %s in ordine alfabetico\n",stringa1,stringa2);

return 0;
}
