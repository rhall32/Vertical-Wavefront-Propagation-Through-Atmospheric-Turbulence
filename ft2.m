function g = ft2(G, delta_f)
% function g = ft2(G, delta_f)
    N = size(G, 1);
    g = fftshift(fft2(fftshift(G))) * (N * delta_f).^2;