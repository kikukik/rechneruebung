% Patrick Nowak, Yannick Gläßer, Tim Rauch, Ben Meyer, 17.01.21 
Out=string(zeros(6,4));
Out(1,1)='';
Out(1,2)='Zielwert Dantzig';
Out(1,3)='Zielwert Bland';
Out(1,4)='Matlab linprog Zielwert';
for i=1:5
    text=append('Beispiel',int2str(i),'.mat');
    load (text);
    Out(i+1,1)=text;
    
    bool=0;
    try [op,Basis,status,it,wert]=SimplexDantzig(A,b,c,Binit,xB);
    catch
        warning('error in init occured')
        outD='error in init';
        bool(1)=1;
    end
    if bool==0
        outD=int2str(wert);
    end
    Out(i+1,2)=outD;
    
    bool=0;
    try [op,Basis,status,it,wert]=SimplexBland(A,b,c,Binit,xB);
    catch
        warning('error in init occured')
        outB='error in init';
        bool=1;
    end
    if bool==0
        outB=int2str(wert);
    end
    Out(i+1,3)=outB;
    
    [x,fval]=linprog(c,A,b,A,b,zeros(1,size(A,2)));
    Out(i+1,4)=int2str(fval);
end
Out

