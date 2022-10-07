function D = phase_struc_func(k, D, L0, Cn2, p)
    a = 11/3;
    k0 = (2*pi)/L0;
    d=2*Cn2*(k^2)*D*gamma(a-1)*sin((a-3)*(pi/2));
    D = d* ( (k0^(2-a)/(a-2)) - ((a*2^(-a/2))/gamma(1+a/2))*((p/k0)^((a/2)-1))* besselk(1-(a/2),k0*p));
end
    