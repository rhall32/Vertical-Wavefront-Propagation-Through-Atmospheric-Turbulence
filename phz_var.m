function phase_variance = phz_var(R, k, Cn2, L0)
    phase_variance = 0.782*R*(k^2)*Cn2*((2*pi/L0)^(-5/3));
end