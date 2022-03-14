%Script for Driven_HW post processing

clear;clc;close all;
%data loading
% load "data.mat" %for large matfiles this isn't optimal

data=matfile("dataN60Re1000.mat");
%setup
Dt=data.Dt;
time=data.time;
x=data.x;
y=data.y;
h=x(2)-x(1);
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

count=1;
reIn_timeStep=data.reIn_timeStep;
REINJECTION_FLAG=data.REINJECTION_FLAG;
%Flow visualization
for it=1:length(time)

    %reInjection 
    if REINJECTION_FLAG && mod(it,reIn_timeStep)==0
        count=count+reIn_timeStep;
    end

% Output
     if mod(it,30) == 0
 
%check divergence
     MaxDiv = max(max(abs(DivCalc(data.U(:,:,it),data.V(:,:,it)))));
     disp(['Time ',num2str(it*Dt)])
     C=max(max(data.U(:,:,it),[],"all"),max(data.V(:,:,it),[],"all"))*Dt/h;
     disp(['Max Courant ', num2str(C)])
     disp(['Massima divergenza sul campo = ',num2str(MaxDiv)])

         if UFLAG
             UNord=u_nord(it);
         end

         %current time
         t = it*Dt;

% Interpolazione ai nodi
        i = 1:Nx;     j = 2:Ny-1;   Iu(i,j)  = (data.U(i,j,it) + data.U(i,j-1,it))/2; 
        Iu(i,1) = USud;             Iu(i,Ny) = UNord;
        i = 2:Nx-1;   j = 1:Ny;     Iv(i,j)  = (data.V(i,j,it) + data.V(i-1,j,it))/2; 
        Iv(1,j) = VWest;            Iv(Nx,j) = VEst;
% Mesh 2D
        X = repmat(x,Ny,1);         Y = repmat(y',1,Nx);
% Grafica
        figure(1); clf; 
        
        %velocity magnitude node points
        pcolor(X,Y,sqrt(Iu.^2 + Iv.^2)');
        colormap jet; shading flat;
        colorbar; caxis([0,1])
        hold on

        %path lines
        for i=1:Np
            plot(Xp(count:it,i),Yp(count:it,i),"-","LineWidth",2);
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

     end
end


%plot velocities and position of p
figure(3)
subplot(2,2,1);
plot(time,Xp)
subtitle("y positon"); xlabel("t");
grid on
subplot(2,2,2)
plot(time,Yp);
%legend(["x","y"]);
subtitle("y positon"); xlabel("t");
grid on
subplot(2,2,3);
plot(time,Up)
subtitle("u velocity"); xlabel("t");
grid on
subplot(2,2,4)
plot(time,Vp); 
%legend(["u","v"]);
subtitle("v velocity"); xlabel("t");
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