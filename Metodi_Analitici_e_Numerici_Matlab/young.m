function E=young(lambda,eps)

% inserire controllo su lunghezz vettori, positività ecc...

E=zeros(length(lambda),1);

for c=1:length(lambda)-1
    lambda_media=(lambda(c)+lambda(c+1))/2;
    eps_media=(eps(c)+eps(c+1))/2;
    E(c)=lambda_media/eps_media;
end

% calcolo il modulo di Young medio per l'ultima coppia di valori
E(length(lambda))=((lambda(length(lambda))+lambda(length(lambda)-1))/2)/((eps(length(lambda))+eps(length(lambda)-1))/2);