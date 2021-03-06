clear all
clc 

%% Input data
d.N = 10; %Number of divisions (panels)
d.L = 1; %Chord
d.p = 0.4; % x coordinate of maximum chamber
d.m = 0.02; % maximum chamber
d.a = 0; % AoA
d.U_inf = 1; %Upstream velocity

%% Coordinates structure

[c.X c.Xc c.Xp c.Nc] = discretization(d.N,d.L,d.a,d.m,d.p); %In order to create the matrices

%% Problem solution

[s.Gamma s.cL gamma_adim]  = solver(c.Xc,c.Nc,c.Xp,d.a,d.U_inf,d.N,d.L);

