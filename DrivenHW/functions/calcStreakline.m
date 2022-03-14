function [Xp,Yp,varargout]=calcStreakline(data,x0,y0,dt,Tf,T0,varargin)
%CALCSTREAKLINE function to evaluate a streakline given source point and
%flow data
%INPUT:
% data: datastructure holding necessary flow data
% x0,y0: source coordinate position
% dt: timestep;
% Tf: final time of streak line integration
%varargin{1}: interpolation type, for types of interpolation see help of
% the function interp1

%OUTPUT:
% Xp,Yp: streak line coordinate array, for the "moving" case each row represent a streak line


%global for the getV function
global x y h
global UNord USud UWest UEst VNord VSud VWest VEst
%setup
Dt=data.Dt;
x=data.x;
y=data.y;
h=x(2)-x(1);
USud=data.USud; UNord=data.UNord; VSud=data.VSud; VNord=data.VNord;
UEst=data.UEst; UWest=data.UWest; VEst=data.VEst; VWest=data.VWest;

%step of timestep
ns=round(dt/Dt);

%Number of timestep

[~,Nt]=min(abs(data.time- Tf));
       
        %initial position
        Xp=x0*ones(1,Nt-ns);
        Yp=y0*ones(1,Nt-ns);
    
        %time indeces
        if T0 ==0
            n0=1;
        else
        n0=round(T0/Dt);
        end
        Nt_times=n0:ns:Nt-ns;
        
        %unsteady cycle
        for i=Nt_times

            %incrementally update particle position
            %     Xp(i)=x0; Yp(i)=y0;
            [Xp(1:i),Yp(1:i)]=updateHeun(data,Dt,i,ns,Xp(1:i),Yp(1:i));
        end

        Xp=Xp(Nt_times);Yp=Yp(Nt_times);
        
        time=data.time;
        if ~isempty(varargin{1})
            Xp=interp1(time(Nt_times),Xp,time(Nt_times(1):Nt_times(end)),varargin{1});
            Yp=interp1(time(Nt_times),Yp,time(Nt_times(1):Nt_times(end)),varargin{1});
        end

end



function [xp,yp]=updateHeun(data,Dt,it,ns,xp_i,yp_i)

% global for getV function
global x y h
global UNord USud UWest UEst VNord VSud VWest VEst

% Heun method
%first stage
xp1=xp_i; yp1=yp_i;
[up1,vp1]=getV(data.U(:,:,it),data.V(:,:,it),xp1,yp1);
%second stage
xp2=xp_i + Dt*ns*up1;  yp2=yp_i+ Dt*ns*vp1;
[up2,vp2]=getV(data.U(:,:,it+ns),data.V(:,:,it+ns),xp2,yp2);
%final stage
xp=xp_i + (Dt/2)*ns*(up1+up2);
yp=yp_i + (Dt/2)*ns*(vp1+vp2);

end