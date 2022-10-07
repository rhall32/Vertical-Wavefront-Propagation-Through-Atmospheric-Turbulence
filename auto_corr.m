function corr = auto_corr(U)
    %deltaf = (1/size(U,1))^2;
    %real(fftshift(ifft(fft(a).*conj(fft(b)))));
    corr = real(ifftshift(ifft2( ifftshift(fftshift(fft2(fftshift(U))).*conj(fftshift(fft2(fftshift(U))))) )));%.*deltaf);
end