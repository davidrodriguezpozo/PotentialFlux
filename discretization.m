function [X Xc Xp Nc] = discretization(N,L,a,m,p)

Delta_x = L/N;

%X: division points for the CVs, N+1
X(:,1) = linspace(0,L,N+1);
%Xp: center of pressure of each CV (at 1/4)
Xp = zeros(N,2);
%Xc: control point of each CV (at 3/4)
Xc = zeros(N,2);
%Nc: normal vectors to the CV line
Nc = zeros(N,2);


for i=1:N
    X(i,2) = chamber(X(i,1),p,m);
    Xp(i,1) = X(i,1)+1/4*Delta_x;
    Xc(i,1) = X(i,1)+3/4*Delta_x;
    Nc1 = -slope(Xc(i,1),p,m);
    Nc2 = 1;
    Nc(i,1) = Nc1 / sqrt(Nc1^2 + Nc2^2);
    Nc(i,2) = Nc2 / sqrt(Nc1^2 + Nc2^2);
    Xp(i,2) = chamber(Xp(i,1),p,m);
    Xc(i,2) = chamber(Xc(i,1),p,m);
end

