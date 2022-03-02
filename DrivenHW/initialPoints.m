function [xp,yp]=initialPoints(Np,Lx,Ly,name,varargin)
% function to set up different initial point distributions


switch name

    case "yline"
        xp=varargin{1}*ones(1,Np)*Lx;
        yp=linspace(0.1,0.9,Np)*Ly;
    case "xline"
        xp=linspace(0.1,0.9,Np)*Lx;
        yp=varargin{1}*ones(1,Np)*Ly;
    case "grid"
        xp=linspace(0.1,0.9,Np)*Lx;
        yp=linspace(0.1,0.9,Np)*Ly;
    case "centreNormal"
        xp=(0.5+0.1*rand(1,Np))*Lx;
        yp=(0.5+0.1*rand(1,Np))*Ly;
    case "singleP"
        xp=varargin{1}*Lx;
        yp=varargin{2}*Ly;
end