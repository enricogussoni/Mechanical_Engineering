#include <stdio.h>

int main ()
{
int n,x,y;
char operazione;

printf("S'inserisca x\n");
scanf("%d",&x);

printf("S'inserisca y\n");
scanf("%d",&y);

printf("Scegliere l'operazione: +,-,*,/");
scanf("%c ",&operazione);

if(operazione='+')
   printf("d",x+y);
else
   if(operazione='-')
      printf("d",x-y);
   else
       if(operazione='*')
          printf("d",x*y);
       else
           if(operazione='/')
              printf("d",x/y);
           else
              printf("Operazione non riconosciuta");


return 0;
}
