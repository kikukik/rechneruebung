function [ xopt,B,message, iter, Zielfktnswert] = SimplexBland( A,b,c,Binit,xB)
%function [ xopt,B ,message, iter] = primalSimplex( A,b,c,Binit,xB )
%
% Primales Simplexverfahren
%
% Input:  A, b, c   - Daten für LP in primaler Standardform
%                      min c'x s.t. Ax=b, x>=0
%         Binit, xB - primal zulaessige Basis, zugehörige Basislösung
%                     (optional)
% Output: xopt      - optimale Lsg
%         B         - zugehörige Basis
%         message   - Information über Optimallösung oder Unbeschraenktheit
%         iter      - Anzahl der Iterationen
%
% Patrick Nowak, Yannick Gläser, Tim Rauch, Ben , DATUM

% Toleranz Definieren!(siehe Blatt)
    tol=1e-6;

% Eingabefehler abfangen
    m=size(A,1);
    n=size(A,2);
    if length(b)~=m || length(c)~=n
        error('A,b,c nicht Dimensionsverträglich')
    end
    if rank(A)~=m
        error('A hat nicht vollen Zeilenrang')
    end
    A_Binit=A(:,Binit);
    if rank(A_Binit)~=m
        error('Binit ist keine Basis')
    elseif any(A_Binit\b<0)
        error('Binit is nicht primal zulaessig')
    elseif nargin==5 && any(A_Binit\b~=xB)
            error('xB stimmen nicht')
    end

% Initialisierung
    B=Binit;
    N=setdiff((1:n),B);
    if nargin<5
        xB=A_Binit/b;
    end
    xopt=zeros(n,1);
    for i=1:m
        xopt(B(i))=xB(i);
    end
    Zielfktnswert=c'*xopt
    
for iter=1:1000
    iter
% Einzelnen Schritte des Algorithmus:

% (1) BTRAN:  
    %loese y'A_B=c_B'
    y=(c(B).'/A(:,B)).';

% (2) Pricing:
    %berechne z_N=c_N-A_N'*y
    z_N=c(N)-A(:,N).'*y;
    if all(z_N>=0)
        message='LP hat optimallsg';
        return
    else
        k=1;
        while z_N(k)>=0
            k=k+1;
        end
        j=N(k);
    end
% (3) FTRAN:
    %loese A_B w=A.k
    w=A(:,B)\A(:,j);
% (4) Ratiotest:
    if any(w<=tol) %all or any ??
        message='LP ist unbeschraenkt';
        Zielfktnswert=NaN;
        xopt=xopt*NaN;
        return
    else
        M=NaN*zeros(1,m);
        for i=1:m
            if(w(i)>0)
                M(i)=(xopt(B(i)))/w(i);
            end
        end
        [gamma,i]=min(M);
    end 
        
% (5) Update:
   % xopt=xB-gamma*w;
    for l=1:m
        xopt(B(l))=xopt(B(l))-gamma*w(l);
    end
    N(k)=B(i);
    B(i)=j;
    xopt(j)=gamma;    
    Zielfktnswert=c'*xopt;
end
