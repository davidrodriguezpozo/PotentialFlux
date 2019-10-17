function [gamma c_L gamma_adim] = solver(Xc,Nc,Xp,a,U_inf,N,L)

M = N-1;


for i=1:M
    for j=1:M
        x_ij = Xc(i,1) - Xp(j,1);
        z_ij = (Xc(i,2) - Xp(j,2));
        r2 = x_ij^2 + z_ij^2; %r^2
        u_ind = z_ij/(2*pi*r2);
        %u_ind = (1/(2*pi*r2))*((Xc(i,2) - Xp(j,2)));
        w_ind = x_ij/(2*pi*r2);
        %w_ind = ((-1)/(2*pi*r2))*((Xc(i,1) - Xp(j,1)));
        A(i,j) = u_ind*Nc(i,1) - w_ind*Nc(i,2);
    end
    %b(i) = -U_inf*(cos(deg2rad(a))*Nc(i,1) + sin(deg2rad(a))*Nc(i,2));
    U = [U_inf*cos(deg2rad(a)) U_inf*sin(deg2rad(a))];
    b_2(i) = -U*Nc(i,:)';
end

%-- Circulation --
gamma = inv(A)*b';
gamma_adim = gamma/(U_inf*L);
%-- Lift coefficient --
c_L = (2/(U_inf*L))*sum(gamma);

end
