#include <stdio.h>

#define max 10
int main()
{
int a1[max],a2[max-1];
int i,n;
int lung1,lung2;
int contenuto,partenza;

{
printf("Quanto � lunga a1?\n");
scanf("%d",&lung1);

if(lung1>max)
   do{
	   printf("a1 non pu� avere pi� di %d elementi\nPrego, inserire la lunghezza di a1\n",max);
	   scanf("%d",&lung1);
	 }while((lung1<0)||(lung1>max));

for(i=0;i<lung1;i++)
    {
	  do{
		  printf("S'inserisca il numero in a1[%d]:\n",i);
	      scanf("%d",&a1[i]);
	    }while((a1[i]<0)||(a1[i]>max));
	}

printf("\nCONTROLLO:La stinga a1 � composta dagli elementi:\n");
for(i=0;i<lung1;i++)
    printf("%d; ",a1[i]);

}

{
printf("\n\nQuanto � lunga a2?\n");
scanf("%d",&lung2);

if(lung1>(lung1-1))
   do{
	   printf("a2 non pu� avere pi� di %d elementi\nPrego, inserire la lunghezza di a2\n",lung1-1);
	   scanf("%d",&lung2);
	 }while((lung2<0)||(lung2>(lung1-1)));

for(i=0;i<lung2;i++)
    {
	  do{
		  printf("S'inserisca il numero in a2[%d]:\n",i);
	      scanf("%d",&a2[i]);
	    }while((a2[i]<0)||(a2[i]>max));
	}

printf("\nCONTROLLO:La stinga a2 � composta dagli elementi:\n");
for(i=0;i<lung2;i++)
    printf("%d; \n",a2[i]);
}

for(n=0;n<=lung1;i++)
    {
     if(a1[n]=a2[0])
        partenza=a1[n];
    }

for(i=0,n=a1[n];i<=lung2,n<))

return 0;
}
