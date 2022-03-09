clear; close all; clc;

%TODO
% - check correctness of fluid flow field

% Esercitazione dal corso di Fluidodinamica Numerica
% 22 novembre 2021
% Prof. G. Coppola
%  
% Codice di esempio per la discretizzazione delle equazioni di
% Navier-Stokes 2D in variabili primitive u,v,p.
% La procedura di integrazione segue il classico metodo di proiezione.
% La discretizzazione spaziale si basa sullo staggering alla Harlow-Welch, 
% nel quale le variabili componenti di velocità e pressione sono 
% localizzate sul mesh secondo lo schema:
% 
% 
%                                V_(i,j+1)
%                            o-----^-----o
%                            |           |
%                            |   p_(i,j) |
%                    U_(i,j) >     x     > U_(i+1,j)
%                            |           |
%                            |           |
%                            o-----^-----o
%                                V_(i,j)
% 
%  con i simboli:
%     x centro-cella (pressione)
%     o vertici di cella
%     > facce verticali (componenti u)
%     ^ facce orizzontali (componenti v)
%
% Il numero di incognite per U e V è dato da:
%
%                   U  --> Nx*(Ny-1)     
%                   V  --> Ny*(Nx-1)     
%                   p  --> (Nx-1)*(Nx-1) 
% 
% La integrazione temporale viene effettuata mediate il classico schema 
% RK4, nel quale la proiezione del campo di velocità sul sottospazio a
% divergenza nulla viene fatta ad ogni stage.
% 

global Nx Ny Re h hq Lapg            % Dichiarazione globale delle varabili
global UNord USud UWest UEst VNord VSud VWest VEst
global x y

%data matfile filename
filename="dataN60.mat";

warning off
Lx = 1;    Ly = 1;                   % Dimensioni del dominio
Nx = 60;   Ny = 60;   N = Nx;        % Numero di nodi lungo ogni lato
x  = linspace(0,Lx,Nx);              % Mesh (uniforme) lungo x
y  = linspace(0,Ly,Ny);              % Mesh (uniforme) lungo y
hx = x(2) - x(1);                    % Passo spaziale lungo x
hy = y(2) - y(1);                    % Passo spaziale lungo y
h  = hx;           hq = h*h;        
Re = 1000;                           % Numero di Reynolds
T  = 30;                             % Tempo finale della simulazione
% Preallocazione delle variabili. 
% Le variabili U e V sono preallocate in array che contengono anche le
% variabili di bordo che cadono sul boundary. Per la variabile P si
% considerano invece solo i valori interni al dominio.
% La istruzione di preallocazione introduce anche le condizioni iniziali.
U  = zeros(Nx,Ny-1);     V = zeros(Nx-1,Ny);    P = zeros(Nx-1,Ny-1);

% Condizioni al contorno
UWest = 0; UEst = 0; VNord = 0; VSud = 0;
UNord = 1; USud = 0; VWest = 0; VEst = 0;
%Flag for periodic lid velocity
UFLAG=1;

Uref  = 1;                           % Velocità di riferimento
C     = 1.2;       beta = 0.8;       % Parametri della discretizzazione numerica
Dt    = min([C*h/Uref,beta*hq*Re]);  % Dt per la stabiità
Nt    = round(T/Dt);                 % Numero di step temporali
% Inserimento delle BCs negli array U e V. Questi valori di bordo non
% saranno aggiornati durante la integrazione.
U(1,:)  = UWest;     U(Nx,:) = UEst; 
V(:,Ny) = VNord;     V(:,1)  = VSud;
% Calcolo dell'operatore di Laplace per l'ellittica di pressione
G = numgrid('S',N+1);    Lap = -delsq(G)/hq;
% Implementazione delle condizioni al contorno alla Neumann sull'operatore
% di Laplace discreto
for i = 2:Nx
    j = 2;    k = G(i,j);   Lap(k,k) = Lap(k,k) + 1/hq;
    j = Ny;   k = G(i,j);   Lap(k,k) = Lap(k,k) + 1/hq;
end
for j = 2:Ny
    i = 2;    k = G(i,j);   Lap(k,k) = Lap(k,k) + 1/hq;
    i = Nx;   k = G(i,j);   Lap(k,k) = Lap(k,k) + 1/hq;
end
% Calcolo dell'operatore di Laplace con condizioni al contorno alla 
% Dirichlet su un mesh collocato, con nodi ai vertici del nostro mesh.
% Questo operatore verrà usato per calcolare la Psi da U e V ai fini della
% grafica.
G   = numgrid('S',Nx);     Lapg = -delsq(G);     Lapg = Lapg/hq;


% condizione iniziale particle tracking
Np=5;
[xp,yp]=initialPoints(Np,Lx,Ly,"yline",0.5);

Xp=nan(Nt,Np);
Yp=nan(Nt,Np);
Up=nan(Nt,Np);
Vp=nan(Nt,Np);
% Avanzamento temporale
time=linspace(0,T,Nt);


% create mat file 
matObj=matfile(filename,"Writable",true);

matObj.time=time;
matObj.Dt=Dt;
matObj.x=x;
matObj.y=y;
matObj.Np=Np;
%velocity variables
matObj.U=zeros(Nx,Ny-1,Nt);
matObj.V=zeros(Nx-1,Ny,Nt);
%store Boundary conditions
matObj.USud=USud; matObj.VSud=VSud;
matObj.UEst=UEst; matObj.VEst=VEst;
matObj.UWest=UWest; matObj.VWest=VWest;
matObj.UNord=UNord; matObj.VNord=VNord;

matObj.UFLAG=UFLAG;
if UFLAG
% set lid boundary velocity
u_nord=setLidVelocity(time,0.05);
matObj.u_nord=u_nord;
end

%Unsteady loop
for it = 1:Nt

    Xp(it,:)=xp;
    Yp(it,:)=yp;

    %update lid driven velocity
    if UFLAG
        UNord=u_nord(it);
    end
%%%%%%%%%%%%%%%%
% fractional step approach: pressure is accurate at first order
% STAGE 1 
    U1 = U;         V1 = V;          [Fu1,Fv1] = RHS_HW(U1,V1);
 %stage 1 path lines
    xp1=xp; yp1=yp;
    [up1,vp1]=getV(U1,V1,xp1,yp1);

%%%%%%%%%%%%%%%%
% STAGE 2
% Calcolo del campo asteriscato
    U2  = U + Dt*0.5*Fu1;             V2 = V + Dt*0.5*Fv1;
% Soluzione della ellittica di pressione
    DIV = DivCalc(U2,V2);             P = Lap\DIV(:);  
% Calcolo del gradiente di pressione    
    P   = reshape(P,Nx-1,Ny-1);      [Px,Py]   = GradCalc(P);
% Step di proiezione
    U2  = U2 - Px;   V2 = V2 - Py;   [Fu2,Fv2] = RHS_HW(U2,V2);
  
  % stage 2 pathlines
  xp2=xp + Dt*0.5*up1;  yp2=yp+ Dt*0.5*vp1;
  [up2,vp2]=getV(U2,V2,xp2,yp2);

%%%%%%%%%%%%%%%%    
% STAGE 3
    U3  = U + Dt*0.5*Fu2;             V3 = V + Dt*0.5*Fv2;
    DIV = DivCalc(U3,V3);             P  = Lap\DIV(:);   
    P   = reshape(P,Nx-1,Ny-1);      [Px,Py]   = GradCalc(P);
    U3  = U3 - Px;   V3 = V3 - Py;   [Fu3,Fv3] = RHS_HW(U3,V3);

  % stage 3 pathlines
  xp3=xp + Dt*0.5*up2;  yp3=yp+ Dt*0.5*vp2;
  [up3,vp3]=getV(U3,V3,xp3,yp3);  

%%%%%%%%%%%%%%%%    
% STAGE 4
    U4  = U + Dt*Fu3;                 V4 = V + Dt*Fv3;
    DIV = DivCalc(U4,V4);             P  = Lap\DIV(:);   
    P   = reshape(P,Nx-1,Ny-1);      [Px,Py] = GradCalc(P);
    U4  = U4 - Px;   V4 = V4 - Py;   [Fu4,Fv4] = RHS_HW(U4,V4);

  % stage 4 pathlines
  xp4=xp + Dt*1*up3;  yp4=yp+ Dt*1*vp3;
  [up4,vp4]=getV(U4,V4,xp4,yp4);

%%%%%%%%%%%%%%%%    
%   Step Finale
    U   = U + Dt*((1/6)*(Fu1 + Fu4) + (1/3)*(Fu2 + Fu3));    
    V   = V + Dt*((1/6)*(Fv1 + Fv4) + (1/3)*(Fv2 + Fv3));
    DIV = DivCalc(U,V);               P  = Lap\DIV(:);
    P   = reshape(P,Nx-1,Ny-1);      [Px,Py] = GradCalc(P);
    U   = U - Px;    V = V - Py;
    
    % stage finale pathlines
    xp= xp + Dt*((1/6)*(up1 + up4) + (1/3)*(up2 + up3));
    yp= yp + Dt*((1/6)*(vp1 + vp4) + (1/3)*(vp2 + vp3));

    %store u and v velocities
    [Up(it,:),Vp(it,:)]=getV(U,V,xp,yp);

matObj.U(:,:,it)=U;
matObj.V(:,:,it)=V;

%check divergence
if mod(it,100)==0
     MaxDiv = max(max(abs(DivCalc(U,V))));
     disp(['Time ',num2str(it*Dt)])
     C=max(max(U,[],"all"),max(V,[],"all"))*Dt/h;
     disp(['Max Courant ', num2str(C)])
     disp(['Massima divergenza sul campo = ',num2str(MaxDiv)])
end

end

%save Lagrangian trajectories and velocities
matObj.Xp=Xp;
matObj.Yp=Yp;
matObj.Up=Up;
matObj.Vp=Vp;

