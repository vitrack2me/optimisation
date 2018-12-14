
clear all;
close all;
load acq;

N = 1024;
precision = 1e-25;


A=1.004929580510462e+03;
C=6.367934463142592;
phi=0.049829472652009;
wb=0.518658841275729;
error = (C + A*sin(wb*(1:1024) + phi)) - data;
X=[C A phi wb]';
F = [sum(C+A*sin(wb*(1:N)+phi));
     sum(sin(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)));
     sum(A*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)));
     sum((1:N)*A.*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)))];
J=[N sum(sin(wb*(1:N)+phi)) sum(A*cos(wb*(1:N)+phi)) sum((1:N)*A.*cos(wb*(1:N)+phi));
   sum(sin(wb*(1:N)+phi)) sum(sin(wb*(1:N)+phi).^2) sum(cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)*A.*cos(wb*(1:N)+phi))) sum((1:N).*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)+phi).*(1:N)*A.*cos(wb*(1:N)+phi));
   sum(A*cos(wb*(1:N)+phi)) sum(cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)*A.*cos(wb*(1:N)+phi))) sum(-A*sin(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+A^2*(cos(wb*(1:N)+phi).^2)) sum(-(1:N)*A.*sin(wb*(1:N)+phi).*(C+A*sin(1:N)-data)+(1:N)*A^2.*(cos(wb*(1:N)+phi).^2));
   sum((1:N)*A.*cos(wb*(1:N)+phi)) sum((1:N).*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)+phi).*(1:N)*A.*cos(wb*(1:N)+phi)) sum(-(1:N)*A.*sin(wb*(1:N)+phi).*(C+A*sin(1:N)-data)+(1:N)*A^2.*(cos(wb*(1:N)+phi).^2)) sum(-((1:N).^2)*A.*sin(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+((1:N).^2)*A^2.*(cos(wb*(1:N)+phi).^2))];
    
while (error > precision)
    X = X-inv(J)*F;
    error = (C + A*sin(wb*(1:1024) + phi)) - data;
    J=[N sum(sin(wb*(1:N)+phi)) sum(A*cos(wb*(1:N)+phi)) sum((1:N)*A.*cos(wb*(1:N)+phi));
       sum(sin(wb*(1:N)+phi)) sum(sin(wb*(1:N)+phi).^2) sum(cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)*A.*cos(wb*(1:N)+phi))) sum((1:N).*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)+phi).*(1:N)*A.*cos(wb*(1:N)+phi));
       sum(A*cos(wb*(1:N)+phi)) sum(cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)*A.*cos(wb*(1:N)+phi))) sum(-A*sin(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+A^2*(cos(wb*(1:N)+phi).^2)) sum(-(1:N)*A.*sin(wb*(1:N)+phi).*(C+A*sin(1:N)-data)+(1:N)*A^2.*(cos(wb*(1:N)+phi).^2));
       sum((1:N)*A.*cos(wb*(1:N)+phi)) sum((1:N).*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)+phi).*(1:N)*A.*cos(wb*(1:N)+phi)) sum(-(1:N)*A.*sin(wb*(1:N)+phi).*(C+A*sin(1:N)-data)+(1:N)*A^2.*(cos(wb*(1:N)+phi).^2)) sum(-((1:N).^2)*A.*sin(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+((1:N).^2)*A^2.*(cos(wb*(1:N)+phi).^2))];
    F=[sum(C+A*sin(wb*(1:N)+phi));
         sum(sin(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)));
         sum(A*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)));
         sum((1:N)*A.*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)))];
end