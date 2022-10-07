function vari=variance(A)
    vari = var(A(:))/(mean(A(:))^2);
end