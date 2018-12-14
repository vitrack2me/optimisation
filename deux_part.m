
clear all;
close all;
load acq;

N = 1024;
% precision = 1e-25;


A=1.004929580510462e+03;
C=6.367934463142592;
phi=0.049829472652009;
wb=0.518658841275729;

error = (C + A*sin(wb*(1:1024) + phi)) - data;

X=[C A wb phi]';

F1_1 = sum(C+A*sin(wb*(1:N)+phi)-data);
F2_1 = sum(sin(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data));
F3_1 = sum((1:N)*A.*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data));
F4_1 = sum(A*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data));


F=[F1_1;
   F2_1;
   F3_1;
   F4_1];

% F = [sum(C+A*sin(wb*(1:N)+phi));
%      sum(sin(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)));
%      sum(A*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)));
%      sum((1:N)*A.*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)))];
J1_1= N;
J1_2 = sum(sin(wb*(1:N)+phi));
J1_3 = sum((1:N)*A.*cos(wb*(1:N)+phi));
J1_4 = sum(A*cos(wb*(1:N)+phi));


J2_1 = sum(sin(wb*(1:N)+phi));
J2_2 = sum(sin(wb*(1:N)+phi).^2);
J2_3 = sum((1:N).*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)+phi).*(1:N)*A.*cos(wb*(1:N)+phi));
J2_4 = sum(cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)*A.*cos(wb*(1:N)+phi)));


J3_1 = sum((1:N)*A.*cos(wb*(1:N)+phi));
J3_2 = sum((1:N).*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)+phi).*(1:N)*A.*cos(wb*(1:N)+phi));
J3_3 = sum(-((1:N).^2)*A.*sin(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+((1:N).^2)*A^2.*(cos(wb*(1:N)+phi).^2));
J3_4 = sum(-(1:N)*A.*sin(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+(1:N)*A^2.*(cos(wb*(1:N)+phi).^2));


J4_1 = sum(A*cos(wb*(1:N)+phi));
J4_2 = sum(cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)*A.*cos(wb*(1:N)+phi)));
J4_3 = sum(-(1:N)*A.*sin(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+(1:N)*A^2.*(cos(wb*(1:N)+phi).^2));
J4_4 = sum(-A*sin(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+A^2*(cos(wb*(1:N)+phi).^2));


J=[J1_1 J1_2 J1_3 J1_4;
   J2_1 J2_2 J2_3 J2_4;
   J3_1 J3_2 J3_3 J3_3;
   J4_1 J4_2 J4_3 J4_4];


% J=[N sum(sin(wb*(1:N)+phi)) sum(A*cos(wb*(1:N)+phi)) sum((1:N)*A.*cos(wb*(1:N)+phi));
%    sum(sin(wb*(1:N)+phi)) sum(sin(wb*(1:N)+phi).^2) sum(cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)*A.*cos(wb*(1:N)+phi))) sum((1:N).*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)+phi).*(1:N)*A.*cos(wb*(1:N)+phi));
%    sum(A*cos(wb*(1:N)+phi)) sum(cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)*A.*cos(wb*(1:N)+phi))) sum(-A*sin(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+A^2*(cos(wb*(1:N)+phi).^2)) sum(-(1:N)*A.*sin(wb*(1:N)+phi).*(C+A*sin(1:N)-data)+(1:N)*A^2.*(cos(wb*(1:N)+phi).^2));
%    sum((1:N)*A.*cos(wb*(1:N)+phi)) sum((1:N).*cos(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+sin(wb*(1:N)+phi).*(1:N)*A.*cos(wb*(1:N)+phi)) sum(-(1:N)*A.*sin(wb*(1:N)+phi).*(C+A*sin(1:N)-data)+(1:N)*A^2.*(cos(wb*(1:N)+phi).^2)) sum(-((1:N).^2)*A.*sin(wb*(1:N)+phi).*(C+A*sin(wb*(1:N)+phi)-data)+((1:N).^2)*A^2.*(cos(wb*(1:N)+phi).^2))];

i = 0;
while (i < 10)
    X = X-inv(J)*F;
%     error = (C + A*sin(wb*(1:1024) + phi)) - data;
    J=[J1_1 J1_2 J1_3 J1_4;
       J2_1 J2_2 J2_3 J2_4;
       J3_1 J3_2 J3_3 J3_3;
       J4_1 J4_2 J4_3 J4_4];
    
    F=[F1_1;
       F2_1;
       F3_1;
       F4_1];
    
    i = i + 1;
end

error = (X(1) + X(2)*sin(X(3)*(1:1024) + X(4))) - data;
v=var(error);

figure(1);
plot(1:N,data);
hold on;
plot(1:N,(X(1) + X(2)*sin(X(3)*(1:1024) + X(4))));
