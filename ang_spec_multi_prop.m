function [xn, yn, Uout] = ang_spec_multi_prop (Uin, wvl, delta1, deltan, z, t)
    N = size(Uin, 1);
    [nx, ny] = meshgrid((-N/2 : 1 : N/2 - 1)); k = 2*pi/wvl; % optical wavevector
    % super-Gaussian absorbing boundary
    %nsq = nx.^2 + ny.^2;
    %w = 0.47*N;
    %sg = exp(-nsq.^8/w^16); clear('nsq', 'w');
    %z = [0 z]; % propagation plane locations 
    n = length(z);
    % propagation distances
    Delta_z = z(2:n) - z(1:n-1);
    % grid spacings
    alpha = z / z(n);
    delta = (1-alpha) * delta1 + alpha * deltan;
    m = delta(2:n) ./ delta(1:n-1);
    %x1 = nx * delta(1);
    %y1 = ny * delta(1);
    %r1sq = x1.^2 + y1.^2;
    %Q1 = exp(1i*k/2*(1-m(1))/Delta_z(1)*r1sq);
    %Uin = Uin .* Q1 .* t(:,:,1);

    % spatial frequencies (of i^th plane) 
    
    for idx=1:(n-1)
    % compute the propagated field
        deltaf = 1 / (N*delta(idx));
        fX = nx * deltaf;
        fY = ny * deltaf;
        fsq = fX.^2 + fY.^2;
        deltaz = Delta_z(idx);
        Q2 = exp(-1i*pi^2*2*deltaz/m(idx)/k*fsq);
        
        deltax = delta(idx);
        Ufreq = fftshift(fft2(fftshift(Uin/m(idx))))*deltax^2;
        Uprop = Q2.*Ufreq;
        Ux = ifftshift(ifft2(ifftshift(Uprop)))*(N*deltaf)^2;
        Uin = t(:,:,idx+1).*Ux;
        
        %{
        uin = ft2(Uin, delta(idx));
        uin = uin.*Q2;
        Uin = ift2(uin, deltaf);
        Uin = t(:,:,idx+1) .* Uin;
        %Uin = Uin/((mean(Uin(:)).^2));
        %}
        %
        %max(Uin(:))
        
    end
    %Uin = normalize(Uin);
    % observation-plane coordinates
    xn = nx * delta(n);
    yn = ny * delta(n);
    rnsq = xn.^2 + yn.^2;
    Q3 = exp(1i*k/2*(m(n-1)-1)/(m(n-1)*deltaz)*rnsq);
    Uout = Q3 .* Uin;