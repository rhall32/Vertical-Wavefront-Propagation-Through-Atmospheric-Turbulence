function c = circ(r, N)
    [xx, yy] = meshgrid(1:N,1:N);
    c = zeros(size(xx));
    c(((xx-N/2+1).^2+(yy-N/2+1).^2)<r^2)=1;
