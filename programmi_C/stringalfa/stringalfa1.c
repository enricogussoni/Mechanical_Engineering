#include <stdio.h>

int main()
{
int main()
{ int c,lungmax,flag=0;
char stringa1[10],stringa2[10];

printf("Inserire stringa1:\n");
scanf("%s",stringa1);
printf("Inserire stringa2:\n");
scanf("%s",stringa2);

if(strlen(stringa1)>strlen(stringa2))
    lungmax=strlen(stringa1);
else lungmax=strlen(stringa2);

for(c=0;(c<lungmax)||(flag=1);c++)
    {
     if (stringa1[c]!= stringa2[c]) flag=1; }
printf("CONTROLLO: c=%d\n",c);

if(stringa1[c]>stringa2[c])
    printf("2 precede 1\n");
else
    printf("1 precede 2\n");

return 0;
}
