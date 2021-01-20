% Patrick Nowak, Yannick Gläßer, Tim Rauch, Ben Meyer, 20.01.2021

% Initialisierung mit Werten der Aufgabenstellung (Form Ax=b, x>=0)
A=[-2 -9 1 9 1 0 ; 1/3 1 -1/3 -2 0 1];
b=[0;0];
Binit=[5 6]; % Matrix mit Schlupfvariablen hat vollen Rang
c=[-2;-3;1;12;0;0];

% Berechnung der Lösung mit SimplexDantzig & Ausgabe der Iterationen der Basis -> Kreiseln
[xopt_SD,B_SD,message_SD, iter_SD, Zielfktnswert_SD] = SimplexDantzig(A,b,c,Binit)

% Alternative: Verwenden des SimplexBand
[xopt_SB,B_SB,message_SB, iter_SB, Zielfktnswert_SB] = SimplexBland(A,b,c,Binit);

