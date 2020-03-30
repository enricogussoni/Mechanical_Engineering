#include <stdio.h>
#include <string.h>

#define dim 10

int main() {

    int c,flag=0;
    char stringa1[dim],stringa2[dim];

printf("Inserire la stringa 1:\n");
scanf("%s",stringa1);
printf("Inserire stringa2:\n");
scanf("%s",stringa2);

for(c=0;(c<dim)||(flag=1);c++)
    {
     if (stringa1[c]!= stringa2[c])
        flag=1;
    }
printf("CONTROLLO: c=%d\n",c);

if(flag==0)
   printf("Stringhe ugulai\n");
else
   {
    if(stringa1[c]>stringa2[c])
       printf("2 precede 1\n");
    else
        printf("1 precede 2\n");
   }

return 0;
}
