#include <stdio.h>
#include <string.h>

#define dim 10
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
    }while(stringa1[c]!=stringa2[c]);

if(stringa1[c]>stringa2[c])
   printf("2 precede 1\n");
else
   printf("1 precede 2\n");

return 0;
}
