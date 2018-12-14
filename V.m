function [ v_n ] = V( n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    C = a0;
    A = sqrt(as.^2 + ac.^2);
    phi = acos(as/A);
    v_n = C + A*sin(wb*n + phi);
end