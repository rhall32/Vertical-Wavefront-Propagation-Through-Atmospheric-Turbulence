function var=analytical_var(D, k, Cn2, L0)
    k0 = 2*pi/L0;
    var = 0.782*D*(k^2)*Cn2*(k0^(-5/3));
end