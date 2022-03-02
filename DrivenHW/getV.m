function [u,v]=getV(U,V,xp,yp)
global x y h UNord VNord

u=interp2(y(1:end-1)+h/2,x,U,xp,yp,"spline");
v=interp2(y,x(1:end-1)+h/2,V,xp,yp,"spline");

%

%wall collision detection

for i=1:length(xp)
    % for Harlow welch u grid
    if (xp(i) > x(end))|| (xp(i) < x(1)) || (yp(i) < y(1)+h/2)
        u=0;

    end

    %top boundary condition
    if yp(i)>y(end)-h/2
        u= UNord;
    end

    % HW v component
    if (xp(i) > x(end)-h/2)|| (xp(i) < x(1)+h/2) || (yp(i) < y(1))
        v=0;
    end

    %top boundary condition
    if yp(i) > y(end)
        v= VNord;
    end
end

end