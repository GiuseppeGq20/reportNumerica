function [Xp,Yp,Up,Vp,time]=calcPathlines(data,Np,x0,y0,dt,T0,Tf)
%CALCPATHLINES  integrate pathlines from known velocity field
%the function implement the Heun method to avoid time interpolation in the
%velocity field
%
%INPUT:
% data: matfile object 
% Np: number of particles
% x0: initial x position of the particles
% y0: initial y position of the particles
% dt: approximate timestep size for particle traking integration
% Varargin{1} (optional): time for particles reinjection
%
%OUTPUT:
% position, velocity ant time arrarys of tracked points

% globals needed for 'getV()' function
global x y h
global UNord USud UWest UEst VNord VSud VWest VEst

%setup
Dt=data.Dt;
time=data.time;
Nt=length(time);
Nt=round(Tf/Dt);
x=data.x;
y=data.y;
h=x(2)-x(1);
USud=data.USud; UNord=data.UNord; VSud=data.VSud; VNord=data.VNord;
UEst=data.UEst; UWest=data.UWest; VEst=data.VEst; VWest=data.VWest;

Xp=nan(Nt,Np);
Yp=nan(Nt,Np);
Up=nan(Nt,Np);
Vp=nan(Nt,Np);
xp=x0;yp=y0;


% particle traking timestep
if dt<Dt
    error("dt must be greater than flow data time step size")
end
ns=round(dt/Dt); %this is to enlarge the timestep size

%set initial time step of integration
if T0==0
    n0=1;
else
n0=round(T0/Dt);
end

Np_times=n0:ns:Nt-ns;

%unsteady cycle
for it=Np_times

    Xp(it,:)=xp;
    Yp(it,:)=yp;
    

    % Heun method
    %first stage
    xp1=xp; yp1=yp;
    [up1,vp1]=getV(data.U(:,:,it),data.V(:,:,it),xp1,yp1);
    %second stage
    xp2=xp + Dt*ns*up1;  yp2=yp+ Dt*ns*vp1;
    [up2,vp2]=getV(data.U(:,:,it+ns),data.V(:,:,it+ns),xp2,yp2);
    %final stage
    xp=xp + (Dt/2)*ns*(up1+up2);
    yp=yp + (Dt/2)*ns*(vp1+vp2);
    %store u and v velocities
    [Up(it,:),Vp(it,:)]=getV(data.U(:,:,it+ns),data.V(:,:,it+ns),xp,yp);

end

Xp=Xp(Np_times,:); Yp=Yp(Np_times,:);
Up=Up(Np_times,:); Vp=Vp(Np_times,:);
time=time(Np_times);

end
