function[y]=U0(x)
    y=1;
 endfunction
 function[h,x]=maillage(a,b,N)
    h=(b-a)/(N+1);
    for i=1:N+1;
        x(i)=i*h;
    end
 endfunction

function[y]=f(x,t)
    y=2*t+x*(1-x);
 endfunction
 a=0;
    b=1;
    Ti=0;
    Tf=1;
    alfa =0;
    betta=0;
    N=10;
    M=5;
    [h,x]=maillage(a,b,N);
    [k,t]=maillage(Ti,Tf,M);
     U=zeros(N+1,M+1);
     for j=1:M+1 
        U(1,j)=alfa;
         U(N+1,j)=betta;
    end
     for i=1:N+1
        U(i,1)=U0(x(i));
     end
     A=zeros(N-1,N-1);
     b=zeros(N-1,1);
     for i=1:N-1
         A(i,i)=(1/k+2/h^2);
     end
     for i=1:N-2
         A(i,i+1)=-1/(h^2);
         A(i+1,i)=-1/(h^2);
     end
     for j=2:M+1
         b(1)=f(x(2),t(j))+(1/k)*U(2,j-1)+(1/h^2)*U(1,j);
          b(N-1)=f(x(N),t(j))+(1/k)*U(N,j-1)+(1/h^2)*U(N+1,j);
          for i=2:N-2
               b(i)=f(x(i+1),t(j))+(1/k)*U(i+1,j-1);
     end
     // Resolution
        U(2:N,j)=inv(A)*b;
        for j=1:M+1
            plot(U(:,j));
            xlabel("position(x)");
            ylabel("temps(t)");
            title("visualisation propagation de la chaleur dans un mur d''epaisseur");
            
            
            end
end

