%Script for Driven_HW post processing and streakline visualization

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

%for periodic lid velocity
UFLAG=data.UFLAG;
if UFLAG
    u_nord=data.u_nord;
end
Np=data.Np;
Nx=length(x);Ny=length(y);
Lx=x(end);Ly=y(end);

%streakline start
t0_streakline=0;
%Flow visualization

%final time
Tf=20;

%unsteady loop
for it=1:length(time)

    if time(it)>Tf
        break
    end

    % Output
    if mod(it,50) == 0

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

        %velocity magnitude node points contour
        %         pcolor(X,Y,sqrt(Iu.^2 + Iv.^2)');
        %         colormap jet; shading flat;
        %         colorbar; caxis([0,1])
        %         hold on

        %streakline
        if time(it)>t0_streakline
            [Xp,Yp]=calcStreakline(data,0.5,0.5,0.1,Dt*it,t0_streakline,"spline");
            plot(Xp,Yp,LineWidth=2,Color="red");
            hold on
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

