function c = corr2_ft_m(u1, u2, mask, delta)
% function c = corr2_ft(u1, u2, mask, delta)
    N = size(u1, 1);
    c = zeros(N);
    %mask = ones(N);
    delta_f = 1/(N*delta); % frequency grid spacing [m]

    U1 = ft2(u1 .* mask, delta); % DFTs of signals
    %U1 = fftshift(fft2(fftshift(u1)))*delta^2;
    U2 = ft2(u2 .* mask, delta);
    %U2 = fftshift(fft2(fftshift(u2)))*(N*delta_f)^2;
    U12corr = ift2(conj(U1) .* U2, delta_f);
    %U12corr = ifftshift(ifft2(ifftshift(U2.*conj(U1))))*(N*delta_f)^2 ;
    
    maskcorr = ift2(abs(ft2(mask, delta)).^2, delta_f);% .* delta^2;
    %mc = abs(fftshift(fft2(fftshift(mask))))*delta^2;
    %maskcorr =(ifftshift(ifft2(ifftshift(mc))))*(N*delta_f)^2;
    idx = logical(maskcorr);
    c(idx) = U12corr(idx) ./ maskcorr(idx) .* mask(idx);
    c = real(c);
