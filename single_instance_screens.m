n_lam = 21;
frames =100;
npix=480;
lam = linspace(400,1000,n_lam);
lam=fliplr(lam);
lam = lam.*1e-9;
k = (2*pi)./lam;

delta=(3.60/(npix/2));
min_r0=10; %planewave
l0 = 1e-3;
L0 = 100;
phzs = zeros(npix,npix,n_lam,frames);
r0=zeros(1,n_lam);
d=zeros(1,n_lam);
for l = 1:n_lam
    r0(l) = min_r0*((lam(n_lam)/lam(l))^(-6/5));
    d(l) = (round((npix/2)*(lam(1)/lam(l))));
end
%r0=fliplr(r0);
r0=r0./100;
%d=fliplr(d);
for f = 1:frames
    seed=round(10*rand(1));
    for l =1:n_lam
        phzs(:,:,l,f) = instance_phase_screen(r0(l), npix, delta, L0, l0,seed);
    end
end
%{
fitswrite(phzs,'phases.fits')
figure()
subplot(2,2,1)
image(phzs(:,:,1,1),'CDataMapping','scaled');
colorbar;
subplot(2,2,2)
image(phzs(:,:,8,1),'CDataMapping','scaled');
colorbar;
subplot(2,2,3)
image(phzs(:,:,15,1),'CDataMapping','scaled');
colorbar;
subplot(2,2,4)
image(phzs(:,:,21,1),'CDataMapping','scaled');
colorbar;
%}