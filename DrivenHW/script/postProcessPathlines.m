%Script for Driven_HW post processing and streakline visualization

clear;clc;close all;
%data loading
% load "data.mat" %for large matfiles this isn't optimal

data=matfile("dataN60Re1000.mat");
Dt=data.Dt;
time_data=data.time;
x=data.x;
y=data.y;
h=x(2)-x(1);
USud=data.USud; UNord=data.UNord; VSud=data.VSud; VNord=data.VNord;
UEst=data.UEst; UWest=data.UWest; VEst=data.VEst; VWest=data.VWest;

%read periodc lid velocity
UFLAG=data.UFLAG;
if UFLAG
    u_nord=data.u_nord;
end
Np=data.Np;
Nx=length(x);Ny=length(y);
Lx=x(end);Ly=y(end);

%pathline
ns=5;
pathDt=ns*Dt; %this ha to be a multiple of Dt
Np=data.Np;
Lx=data.x(1,end); Ly=data.y(1,end);

%[x0,y0]=initialPoints(Np,Lx,Ly,"yline",0.5);
%[x0,y0]=initialPoints(3,Lx,Ly,"xline",0.8);
[x0,y0]=initialPoints(3,Lx,Ly,"single",0.5,0.5);

T0=0;Tfp=30;
[Xp,Yp,Up,Vp,time]=calcPathlines(data,length(x0),x0,y0,pathDt,T0,Tfp);


%simulation time
Tf=Tfp;
%Flow visualization
n0=round(T0/Dt);
for it_path=1:length(time)

    if time(it_path)> Tf
        break
    end
    it=ns*it_path+n0;
% Output
     if mod(it,20) == 0
 
%check divergence
     MaxDiv = max(max(abs(DivCalc(data.U(:,:,it),data.V(:,:,it)))));
     disp(['Time ',num2str(it*Dt)])
     C=max(max(data.U(:,:,it),[],"all"),max(data.V(:,:,it),[],"all"))*Dt/h;
     disp(['Max Courant ', num2str(C)])
     disp(['Massima divergenza sul campo = ',num2str(MaxDiv)])

         if UFLAG
             UNord=u_nord(it);
         end
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
%         pcolor(X,Y,sqrt(Iu.^2 + Iv.^2)');
%         colormap jet; shading flat;
%         colorbar; caxis([0,1])
%         hold on

        %path lines
        if time_data(it)>=T0 && time_data(it)<=Tfp
        for j=1:size(Xp,2)
            plot(Xp(1:it_path,j),Yp(1:it_path,j),"-","LineWidth",2);
            hold on; 
        end
        end
        %stream lines
        l=streamslice(X,Y,Iu',Iv');
        set(l,'Color','k')
        axis([0 Lx 0 Ly ]);
        axis square;
        title(['Harlow-Welch. Driven cavity. t = ',num2str(t)]);
        hold off;
        drawnow;

     end
end

figure(2)
%plot velocities and position of particles tracked
subplot(2,2,1);
plot(time,Xp)
subtitle("x positon"); xlabel("t");
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