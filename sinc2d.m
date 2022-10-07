function z=sinc2d(x,y)
    x_sinc = sinc(x);
    y_sinc = sinc(y);
    z = x_sinc.*y_sinc;