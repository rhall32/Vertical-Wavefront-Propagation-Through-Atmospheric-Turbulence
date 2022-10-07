%%%%%%%%%%determine geometry%%%%%%%%%%%
n_lam=1;
%cn2;
n_trial=10; % number of data frames

lam = 671;
lam = lam.*1e-9;

%%%%%%%%%Physical properties%%%%%%%%%%
%lam =500e-9;
k = (2*pi)./lam;
D2 = 0.4; % diameter of the observation aperture [m]
Dz = 2080; % propagation distance [m]
elevation = 0;%18e3; % observation elevation
N = 2048;
N_final = 2048;
z = pi/2; % angular distance from zenith [radians]

%%%%%%%%%screen properties%%%%%%%%%%
n = 5; % number of screens
l0 = 1e-3;
L0 = 3;
%l0= linspace(1e-3,3e-2,n);
%L0 = linspace(3,300,n);
%rytov = 0.1;

%%%%%%%%%grid spacings%%%%%%%%%
delta1= 0.83e-3; %source-plane grid spacing [m]
deltan= 0.83e-3; %observation-plane grid spacing [m]
dz=((1:n-1)*Dz/(n-1))+elevation;
dz=[elevation dz];
alpha = dz / (dz(n));
delta = ((1-alpha)*delta1)+(alpha*deltan);
del_z = dz(n)-dz(n-1);

%DROI=5; %diam of obs-plane region of interest [m]
%D1 = lam*Dz / DROI;
arg = D2/(lam*Dz);
R = Dz; % wavefront radius of curvature [m]
x = (-N/2 : N/2 - 1).*delta1;
[x1, y1] = meshgrid(x);
[theta1, r1] = cart2pol(x1, y1);
Cn2 = ones(1,n);
Cn2 = Cn2.*1.0 * (10.0^(-14));
%Cn2(:) = mean(Cn2);

amp_U=zeros(N_final,N_final,n_lam,n_trial);
phz = zeros(N,N,n,n_lam);
unwrap=zeros(N_final,N_final,n_lam,n_trial);

all_st_func = zeros(N_final/2,n_lam,n_trial);
all_corr_func = zeros(N_final/2, n_lam,n_trial);
all_var = zeros(1,n_lam,n_trial);
all_U_var = zeros(1,n_lam,n_trial);
intensity_var = zeros(1,n_lam,n_trial);
intensity_var_U = zeros(1,n_lam,n_trial);
st_func = zeros(N_final/2,n_lam,n_trial);
an_corr = zeros(N_final/2,n_lam,n_trial);
r0i = zeros(n, n_lam);
MCF2 = zeros(N);
MCF2_col = zeros(N);

%trial_seed= abs(randperm(1000,n_trial));
for frame = 1:n_trial
    frame
    for l = 1:n_lam

    %%%%%%%%%atmospheric properties%%%%%%%%%

        r0pw = (0.423 * (k(l).^2) * Cn2 * del_z).^(-3.0/5); %planewave

        for idx = 1 : n
            phz(:,:,idx,l) = phase_screen(r0pw(idx), N, delta(idx), L0, l0);
        end
        
%{
        w = round(0.1524/delta1);
        sg=rect(w,N);
        a = 0;
        b = 2*pi;
        sqr_phz = (b-a).*rand(N) + a;
        U1 = sg.*exp(1i.*sqr_phz);
        t = ones(N,N,n);
        t = cat(3,U1, t);
        [~, ~, U2] = ang_spec_multi_prop(U1, lam(l), delta1, deltan, dz, t);
        
        w2=(round((2*D2)/delta1));
        U2_prime = rect(w2,N).*U2;
        t(:,:,1) = U2_prime;
        [~, ~, U1_prime] = ang_spec_multi_prop(U2_prime, lam(l), delta1, deltan, -dz, t);
        %}
        x = linspace(-N/2,N/2,N);
        [X, Y] = meshgrid(x);
        X = X./(lam*Dz/(2*D2)/deltan); Y = Y./(lam*Dz/(2*D2)/deltan);
        [~, r1] = cart2pol(x1, y1);
        A = lam * Dz / deltan;
        pt = A * exp(-i*k/(2*Dz/deltan) * r1.^2) * (arg/deltan)^2 .* sinc(X) .* sinc(Y) .*exp(-(r1/(4*D2)).^2);
        
        t = ones(N,N,n);
        t = cat(3,pt, t);
        [~, ~, U2] = ang_spec_multi_prop(pt, lam(l), delta1, deltan, dz, t);
        
        sg = exp(-(x1/(0.47*N)).^16) .* exp(-(y1/(0.47*N)).^16);
        t = sg.*exp(1i.*phz(:,:,:,l));
        %U1_prime = sg.*exp(1i.*ones(N));
        %t = cat(3,pt, t);%cat(3,sg.*U1_prime, t);
        [xn, yn, U] = ang_spec_multi_prop(pt, lam(l), delta1, deltan, dz, t);
        U_col = U.*exp(-i*pi/(lam*Dz)*(xn.^2+yn.^2));
        ap=circ(round(D2/2/deltan),N);
        U_ap = U.*ap;
        U_col_ap = U_col.*ap;
        amp_U(:,:,l,frame) = abs(U((N/2)-(N_final/2)+1:(N/2)+(N_final/2),(N/2)-(N_final/2)+1:(N/2)+(N_final/2)));
        intensity = amp_U.^2;
        intensity_var_U(l) = variance(intensity);

        holder = GoldsteinUnwrap(U);
        a="here2";
        holder(isnan(holder))=0;
        unwrap(:,:,l,frame)=holder((N/2)-(N_final/2)+1:(N/2)+(N_final/2),(N/2)-(N_final/2)+1:(N/2)+(N_final/2));
   
        %%%%%%%simulated%%%%%%%
        corr_func = corr2_ft(unwrap(:,:,l,frame), unwrap(:,:,l,frame), delta1);
        t_st_func = phz_st_func(corr_func);
        MCF2 = MCF2 + corr2_ft_m(U_col, U_col, ap, deltan);
        MCF2_col = MCF2_col + corr2_ft_m(U_col_ap, U_col_ap, ap,deltan);

        all_st_func(:,l,frame) = t_st_func; 
        all_corr_func(:,l,frame) = corr_func(N_final/2,N_final/2+1:N_final);
        all_var(:,l,frame) = variance(unwrap(:,:,l,frame));
        all_U_var(:,l,frame) = var(U(:));

        %%%%%%analytical%%%%%%
        i=1;
        for a=delta1:delta1:N_final*delta1/2
            an_corr(i,l) = analytical_corr_func(k(l), Dz, L0(1), sum(Cn2), a);
            st_func(i,l) = phase_struc_func(k(l), Dz, L0(1), sum(Cn2), a);
            i=i+1;
        end
        clear i

        %clear U unwrap 

    end
end


avg_st_func = mean(all_st_func,3);
avg_corr_func = mean(all_corr_func,3);
avg_var = mean(all_var(:));
MCDOC2 = abs(MCF2) / (MCF2(N/2+1,N/2+1));
MCDOC2_col = abs(MCF2_col) / (MCF2_col(N/2+1,N/2+1));

p = linspace(0,Dz,n);
rytov = zeros(1,n_lam);
iso_theta = zeros(1,n_lam);
for l = 1:n_lam
    rytov(l) =  0.563 * k(l)^(7/6) * Dz* sum( Cn2 .* (1-p./Dz).^(5/6) .* (p(2)-p(1)));
    iso_theta(l) =  (2.91 * k(l)^(2) * sum( Cn2 .* (p).^(5/3) .* (p(2)-p(1))))^(-3/5);
end
iso_theta = rad2deg(iso_theta).*3600; %arc-seconds
%r0 = (0.423 *(cos(z)^-1) * (k.^2) * sum(Cn2(1:n)*del_z)).^(-3.0/5);
r0 = (0.423 * (k.^2) * sum(Cn2(1:n)*del_z)).^(-3.0/5);
r=0:delta1:(delta1*N/2-delta1);
mu = exp(-1.46* k.^(2) * r.^(5/3.0) * sum(Cn2(1:n))*del_z);
%{
figure
title('Object Plane')
subplot(2,2,1)
imagesc(abs(pt)), colormap(gray), colorbar, axis square, title('Obj Amp');
subplot(2,2,2),imagesc(angle(pt)), colormap(gray), colorbar, axis square, title('Obj Phz');
subplot(2,2,3),imagesc(abs(U2)), colormap(gray), colorbar, axis square, title('Receiver (Vacuum) Amp');
subplot(2,2,4),imagesc(angle(U2)), colormap(jet), colorbar, axis square, title('Receiver (Vacuum) Phz');

figure
title('Receiver Plane')
subplot(2,2,1)
imagesc(abs(U2)), colormap(gray), colorbar, axis square, title('Original Amp');
subplot(2,2,2),imagesc(angle(U2)), colormap(gray), colorbar, axis square, title('Original Phz');
subplot(2,2,3),imagesc(abs(U2_prime)), colormap(gray), colorbar, axis square, title('Windowed Amp');
subplot(2,2,4),imagesc(angle(U2_prime)), colormap(jet), colorbar, axis square, title('Windowed Phz');

figure
title('Recovered')
subplot(2,2,1)
imagesc(abs(U_col_ap)), colormap(gray), colorbar, axis square, title('Post-Turbulence Amp');
subplot(2,2,2),imagesc(angle(U_col_ap)), colormap(gray), colorbar, axis square, title('Post-Turbulence Phz');
subplot(2,2,3),imagesc(abs(U_ap)), colormap(gray), colorbar, axis square, title('Post-Turbulence Amp (Aperture)');
subplot(2,2,4),imagesc(angle(U_ap)), colormap(jet), colorbar, axis square, title('Post-Turbulence Phz (Aperture)');

x = linspace(-1024,1024,8192);
y = sinc(x*(lam*Dz/0.1524/deltan));
y = y(8192/2+1:8192);
x1=0:1:1023;
x2=linspace(0,1024,4096);
fig=figure
hax=axes;
hold on
plot(x1, MCDOC2(N_final/2+1,N_final/2+1:N_final))
plot(x2, y);
SP=(lam*Dz/0.1524)/deltan; %your point goes here 
line([SP SP],get(hax,'YLim'),'Color',[1 0 0])
axis([0 20 0 inf])
title('Coherence Factor')
grid

x = linspace(-1024,1024,2048);
y = sinc(x/(lam*Dz/(2*D2)/deltan));
y = y(N/2+1:N);
x= x(N/2+1:N);

fig=figure
hax=axes;
hold on
plot(r/r0, MCDOC2_col(N_final/2+1,N_final/2+1:N_final))
plot(r/r0,mu)
%SP=(lam*Dz/(2*D2)); %your point goes here 
%line([SP SP],get(hax,'YLim'),'Color',[1 0 0])
axis([0 1.5 0 inf])
title('Coherence Factor')
grid
%}

