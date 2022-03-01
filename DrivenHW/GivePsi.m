function PSI = GivePsi(U,V)

% Function per il calcolo della funzione di corrente Psi per un dato campo
% di velocita' staggerato alla Harlow-Welch
% La funzione di corente e' ottenuta calcolando preliminarmente il campo di
% vorticità e risolvendo quindi una equazione di Poisson
% 

global h Lapg
[Nx,~] = size(U);           [~,Ny] = size(V);
Uy     = nan(Nx,Ny);        Vx     = nan(Nx,Ny);
% Calcolo delle derivate di velocita'
i = 2:Nx-1;  j = 2:Ny-1;
Uy(i,j) = (U(i,j)-U(i,j-1))/h;
Vx(i,j) = (V(i,j)-V(i-1,j))/h;
% Calcolo del campo di Zita dentro al dominio
ZITA = Uy - Vx;             ZITA = ZITA(2:end-1,2:end-1);
% Risoluzione dela equazione di Poisson
zita = ZITA(:);
Psi  = Lapg\zita;
% Agginta delle BCs
PSI  = reshape(Psi,Nx-2,Ny-2);
PSI  = [zeros(1,Nx-2);PSI;zeros(1,Nx-2)];
PSI  = [zeros(Ny,1),PSI,zeros(Ny,1)];







end

