function [u,v]=getV(U,V,xp,yp)
% GETV function to query velocity at location of lagrangian particles
%INPUT
% U: u velocity field given in the classical Harlow Welch grid
% V: v velocity field given in the classical Harlow Welch grid
% xp: x position vector of the lagrangian particles
% yp: y postition vector of the lagrangian particles
%OUTPUT
% u: u velocity vector of the particles
% v: v velocity vector of the particles

global x y h UNord 

%TODO compare two approaces: boundary included or if
y_interp=[0,y(1:end-1)+h/2,y(end)];
x_interp=[0,x(1:end-1)+h/2,x(end)];

U_interp=[zeros(1,size(U,1));U';UNord.*ones(1,size(U,1))];
V_interp=[zeros(size(V,2),1), V',zeros(size(V,2),1)];

u=interp2(x,y_interp,U_interp,xp,yp,"linear");
v=interp2(x_interp,y,V_interp,xp,yp,"linear");
%u=interp2(x,y(1:end-1)+h/2,U',xp,yp,"linear");
%v=interp2(x(1:end-1)+h/2,y,V',xp,yp,"linear");



%wall collision detection
% for i=1:length(xp)
    
%     % Harlow welch u component
%     if (xp(i) > x(end))|| (xp(i) < x(1)) || (yp(i) < y(1)+h/2)
%         u(i)=0;
%     end
% 
%     %top boundary condition
%     if yp(i)>y(end)-h/2
%         u(i)= UNord;
%     end

    % HW v component
%     if (xp(i) > x(end)-h/2)|| (xp(i) < x(1)+h/2) || (yp(i) < y(1))
%         v(i)=0;
%     end
% 
%     %top boundary condition
%     if yp(i) > y(end)
%         v(i)= VNord;
%     end
% end

end