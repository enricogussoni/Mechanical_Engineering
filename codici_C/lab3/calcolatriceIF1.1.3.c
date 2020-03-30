#include <stdio.h>

int main ()
{
int n,x,y;
int operazione;

printf("S'inserisca x\n");
scanf("%d",&x);

printf("S'inserisca y\n");
scanf("%d",&y);

do{
   printf("Scegliere l'operazione: 1)+,\n2)-,3)\n*,4)\n/\n");
   scanf("%d",&operazione);
  }while((n<1)||(n>4));

if(operazione==1)
   printf("d",x+y);
else
   if(operazione==2)
      printf("d",x-y);
   else
       if(operazione==2)
          printf("d",x*y);
       else
          printf("d",x/y);



return 0;
}
