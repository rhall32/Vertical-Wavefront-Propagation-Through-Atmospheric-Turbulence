function st_func = phz_st_func(corr)
    %cov = corr ./ mean(corr(:)).^2;
    %st_func = (2.*(vrnc-cov));%+100;
    N=size(corr,1);
    st_func = 2.*(corr(N/2,N/2)-corr(N/2,N/2+1:N));
end