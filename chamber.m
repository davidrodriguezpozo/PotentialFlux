function y = chamber(x,p,m)

if x < p
    y = m/p^2*(2*p*x-x^2);
else
    y = m/(1-p)^2*((1-2*p)+2*p*x-x^2);
end

end