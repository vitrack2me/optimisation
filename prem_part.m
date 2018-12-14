
clear all;
close all;
load acq;

N = 1024;

wb = 2*pi*fin/fs;
wbn = wb*(1:1024);
sin_wbn = sin(wbn);
cos_wbn = cos(wbn);
sin2_wbn = sin(wbn).^2;
cos2_wbn = cos(wbn).^2;
sincos_wbn = cos(wbn).*sin(wbn);

M = [N sum(sin_wbn) sum(cos_wbn); sum(sin_wbn) sum(sin2_wbn) sum(sincos_wbn); sum(cos_wbn) sum(sincos_wbn) sum(cos2_wbn)];
D = [sum(data);sum(data.*sin_wbn);sum(data.*cos_wbn)];

res = inv(M)*D;

a0 = res(1);
as = res(2);
ac = res(3);

C = a0;
A = sqrt(as.^2 + ac.^2);
phi = acos(as/A);

error = (C + A*sin(wb*(1:1024) + phi)) - data;
v = var(error);

figure(1);
plot(1:N,error)


%B-I-5
w0 = wb;
pas = w0*(4*10 ^-4)/200;
W=(1:200)*pas + w0*(1-2*10^-4);


C_mat=(1:200);
A_mat=(1:200);
phi_mat=(1:200);
error_mat=zeros(200,1024);


for (i=1:1:200)
    N = 1024;
    wb = W(i);
    wbn = wb*(1:1024);
    sin_wbn = sin(wbn);
    cos_wbn = cos(wbn);
    sin2_wbn = sin(wbn).^2;
    cos2_wbn = cos(wbn).^2;
    sincos_wbn = cos(wbn).*sin(wbn);
    
    M = [N sum(sin_wbn) sum(cos_wbn); sum(sin_wbn) sum(sin2_wbn) sum(sincos_wbn); sum(cos_wbn) sum(sincos_wbn) sum(cos2_wbn)];
    D = [sum(data);sum(data.*sin_wbn);sum(data.*cos_wbn)];
    
    res = inv(M)*D;
    
    a0 = res(1);
    as = res(2);
    ac = res(3);
    C_mat(i)= a0;
    A_mat(i) = sqrt(as.^2 + ac.^2);
    phi_mat(i)= acos(as/A_mat(i));
    error_mat(i,:)=(C_mat(i) + A_mat(i)*sin(wb*(1:1024) + phi_mat(i))) - data;
    Var = var(error_mat');    
end
figure(2);
plot(W,C_mat);
figure(3);
plot(W,A_mat);
figure(4);
plot(W,phi_mat);
figure(5);
plot(W,Var);





