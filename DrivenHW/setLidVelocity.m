function  [u_nord,v_nord] = setLidVelocity(t,freq,varargin)
%SETLIDVELOCITY set moving lid boundary condition
%INPUT:
% t: time [s]
% freq: frequency of the u and v velocity component [Hz]
% varargin{1}: x stencil of the lid
% varargin{2}: Lx, length of the lid
%OUTPUT:
% u_nord: u velocity component on the lid
% v_nord (optional): v velocity component on the lid

%u component
u_nord=sin(2*pi*freq*t);

%optional v component
if nargout==2 && nargin==4
    x=varargin{1}; Lx=varargin{2};
    v_nord=sin(4*pi*x/Lx)*sin(2*pi*t*freq);
end

end