function [xp,yp]=initialPoints(Np,Lx,Ly,name,varargin)
%INITIALPOINTS function to set up different initial point distributions
%INPUTS:
% Np: number of particles
% Lx: length of the domanin in the x direction
% Ly: length of the domanin in the y direction
% name: initial condition considered
% varargin{1}: x positon percentage
% varargin{2}: y positon percentage
%OUTPUT:
% xp: x coordinates position of the particles
% yp: y coordinates position of the particles

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