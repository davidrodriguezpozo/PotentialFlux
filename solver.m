function [gamma c_L gamma_adim] = solver(Xc,Nc,Xp,a,U_inf,N,L)

M = N-1;


for i=1:M
    for j=1:M
        r2 = (Xc(i,1) - Xp(j,1))^2 + (Xc(i,2) - Xp(j,2))^2;
        u_ind = (1/(2*pi*r2))*((Xc(i,2) - Xp(j,2)));
        w_ind = ((-1)/(2*pi*r2))*((Xc(i,1) - Xp(j,1)));
        A(i,j) = u_ind*Nc(i,1) + w_ind*Nc(i,2);
    end
    b(i) = -U_inf*(cos(deg2rad(a))*Nc(i,1) + sin(deg2rad(a))*Nc(i,2));
end

%-- Circulation --
gamma = inv(A)*b';
gamma_adim = gamma/(U_inf*L);
%-- Lift coefficient --
c_L = (2/(U_inf*L))*sum(gamma);

end
