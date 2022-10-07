function f=convolve(a,b,delta)
    f = real(ift2(ft2(a,delta).*ft2(b,delta),delta));
end