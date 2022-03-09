%Script for Driven_HW post processing

clear;clc;close all;
%data loading
% load "data.mat" %for large matfiles this isn't optimal

data=matfile("data.mat");
Dt=data.Dt;
time=data.time;
x=data.x;
y=data.y;
USud=data.USud; UNord=data.UNord; VSud=data.VSud; VNord=data.VNord;
UEst=data.UEst; UWest=data.UWest; VEst=data.VEst; VWest=data.VWest;

UFLAG=data.UFLAG;
if UFLAG
    u_nord=data.u_nord;
end
Np=data.Np;
Xp=data.Xp; Yp=data.Yp; Up=data.Yp; Vp=data.Vp;
Nx=length(x);Ny=length(y);
Lx=x(end);Ly=y(end);

%Flow visualization
for it=1:length(time)
    
% Output
     if mod(it,30) == 0

         if UFLAG
             UNord=u_nord(it);
         end
         t = it*Dt;
     %Verifica sulla divergenza alla fine dello step
         %MaxDiv = max(max(abs(DivCalc(U,V))));
         %disp(['Massima divergenza sul campo = ',num2str(MaxDiv)])

% Interpolazione ai nodi
        i = 1:Nx;     j = 2:Ny-1;   Iu(i,j)  = (data.U(i,j,it) + data.U(i,j-1,it))/2; 
        Iu(i,1) = USud;             Iu(i,Ny) = UNord;
        i = 2:Nx-1;   j = 1:Ny;     Iv(i,j)  = (data.V(i,j,it) + data.V(i-1,j,it))/2; 
        Iv(1,j) = VWest;            Iv(Nx,j) = VEst;
% Mesh 2D
        X = repmat(x,Ny,1);         Y = repmat(y',1,Nx);
% Grafica
        figure(1); clf; 
        %pressure
        % pcolor(X(1:end-1,1:end-1)+h/2,Y(1:end-1,1:end-1)+h/2,P); shading interp 
        % colormap jet;
        
        %velocity magnitude node points
        pcolor(X,Y,sqrt(Iu.^2 + Iv.^2)');
        colormap jet; shading flat;
        colorbar; caxis([0,1])
        hold on

        %path lines
        for i=1:Np
            plot(Xp(1:it,i),Yp(1:it,i),"-","LineWidth",2);
            hold on; 
            % plot(xp(i),yp(i),"o","MarkerSize",8); %plot current position
            % hold on
        end
        %stream lines
        l=streamslice(X,Y,Iu',Iv');
        set(l,'Color','w')
        axis([0 Lx 0 Ly ]);
        axis square;
        title(['Harlow-Welch. Driven cavity. t = ',num2str(t)]);
        hold off;
        drawnow;

%         figure(2); clf;
%         PSI  = GivePsi(U(:,:,it),V(:,:,it));                 % Calcolo della Psi da U e V
%         vneg = linspace(min(min(PSI)),0,10);
%         vpos = linspace(0,max(max(PSI)),10);
%         contour(x,y,PSI',vneg,'k'); axis square; hold on;
%         contour(x,y,PSI',vpos,'r'); 
%         title({'Harlow-Welch. Driven cavity. Streamlines da \psi'})
%         drawnow; hold off
     end
end


%plot velocities and position of p
figure(3)
subplot(2,2,1);
plot(time,Xp)
subtitle("p y positon"); xlabel("t [s]");
grid on
subplot(2,2,2)
plot(time,Yp);
%legend(["x","y"]);
subtitle("p y positon"); xlabel("t [s]");
grid on
subplot(2,2,3);
plot(time,Up)
subtitle("p u velocity"); xlabel("t [s]");
grid on
subplot(2,2,4)
plot(time,Vp); 
%legend(["u","v"]);
subtitle("p v velocity"); xlabel("t [s]");
grid on

% plot of complete lagrangian trajectory
figure(4)
for i=1:Np
    plot(Xp(:,i),Yp(:,i),"-","LineWidth",2);
    hold on; 
end
ylim([0,Ly]);xlim([0,Lx]); axis square
hold off  
grid on