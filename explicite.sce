function [y]=fT(X,t)
    y=2*t+X*(1-X);
endfunction
function[h,X]=maillage(a,b,N)
    h=(b-a)/(N+1);
    for i=1:N+1;
        X(i)=i*h;
    end
 endfunction

a=0;
b=1;
c=1;
Ti=0;
Tf=1;
alfaa=0;
betaa=1;
gamma=0
N=10;
M=201;
[h,X]=maillage(a,b,N);
[k,t]=maillage(Ti,Tf,M);
U=zeros(N+1,M+1);  // U est une matrice a N+1 lgs et M+1 cls
  // entree des conditions aux limites
  for j=1:M+1
      U(1,j)=alfaa;
      U(N+1,j)=betaa
      
  end
  //conditions initiales
  for i=1:N+1
      U(i,1)=1;
  end
  for j=1:M
      for i=2:N
          U(i,j+1)=(1-2*c*k/h^2)*U(i,j)+(c*k/h^2)*U(i+1,j)+(c*k/h^2)*U(i-1,j)+k*fT(X(i),t(j));
      end
  end
  //representation  de la solution
  for j=1:M+1
  plot(U(:,j));
  xpause(1);                   // pause sur laffichage
end                         //on fixe j et on parcours i // on voit la solotion                               n'est pas stable i=1..2..3
xlabel("positon(x)");
ylabel("temps(t)");
title("visualisation propagation de la chaleur dans un mur d''epaisseur l")
 
  // pour faire une pause sur laffichage des courbes on 
  
   //ffor j=1:M+1
  //plot(U(:,j));
 // xpause(1);
  //end          
  
