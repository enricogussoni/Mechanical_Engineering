n=input('Dimensione matrice:\n')
t=zeros(n);
for x=1:n;
	for y=1:n;
		t(x,y)=t(1,x)*t(y,1);
		end
		
end
t