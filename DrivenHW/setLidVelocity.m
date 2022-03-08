function  [u_nord,v_nord] = setLidVelocity(t,freq,varargin)

u_nord=sin(2*pi*freq*t);

if nargout==2 && nargin==4
    x=varargin{1}; Lx=varargin{2};
    v_nord=sin(4*pi*x/Lx)*sin(2*pi*t*freq);
end
