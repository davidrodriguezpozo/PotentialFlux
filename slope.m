function dzdx = slope(x,p,m)

if x<p
    dzdx=(2*m*(p-x))/p^2;
else
    dzdx=(2*m*(p-x))/(1-p^2);
end

end