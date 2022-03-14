%plot quantity from foamLog output

meanCourant=readmatrix("logs/CourantMean_0");
maxCourant=readmatrix("logs/CourantMax_0");
figure(1)
plot(meanCourant(:,1),meanCourant(:,2), maxCourant(:,1),maxCourant(:,2))
title("Courant Number")
xlabel("t")
legend("mean","max")
grid on

pRes=readmatrix("logs/pFinalRes_0");
figure(2)
semilogy(pRes(:,1),pRes(:,2));
title("pressure residuals")
xlabel("t")
grid on

uxRes=readmatrix("logs/Ux_0");
uyRes=readmatrix("logs/Uy_0");
figure(3)
semilogy(uxRes(:,1),uxRes(:,2),uyRes(:,1),uyRes(:,2));
title("velocity component residuals")
legend("Ux", "Uy");
xlabel("t");
grid on