function corr=correlate(a, b)
    corr = real(ifftshift(ifft(ft2(a, 1).*conj(ft2(b, 1)))));
end