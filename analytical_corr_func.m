function B=analytical_corr_func(k, D, L0, Cn2, p)
    a = (11/3);
    k0 = (2*pi)/L0;
    b = (a*2^(-a/2)*gamma(a-1))/gamma(1+(a/2));
    B = abs(b*sin((a-3)*(pi/2)) * k^2 * D * Cn2 * ((p/k0)^((a/2)-1)) * besselk(1-(a/2),k0*p));
    

    