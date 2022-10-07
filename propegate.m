%%%%%%%%%%determine geometry%%%%%%%%%%%
n_lam=40;
%n_lam=3;
%cn2;
n_trial=11; % number of data frames

%%%%%%%%Polychromatic%%%%%%%%%%
lam = linspace(400,1000,n_lam);
lam=fliplr(lam);
lam = lam.*1e-9;

%%%%%%%%Monochromatic%%%%%%%%%%
%lam = 637e-9;

%%%%%%%%%Physical properties%%%%%%%%%%
k = (2*pi)./lam;
D2 = 4.3; % diameter of the observation aperture [m]
Dz = 30e3; % propagation distance [m]
elevation = 2400;% observation elevation
N = 2048;
N_final = 1024;
z = 0; % angular distance from zenith [radians]

%%%%%%%%%screen properties%%%%%%%%%%
n = 10; % number of screens
l0 = linspace(3e-3,3e-2,n);
L0 = linspace(10,2000,n);
%l0= linspace(1e-3,3e-2,n);
%L0 = linspace(3,300,n);
%rytov = 0.1;

%%%%%%%%%grid spacings%%%%%%%%%
delta1= 4*D2/N; %3.52e-3; %source-plane grid spacing [m]
deltan= 4*D2/N; %3.52e-3; %observation-plane grid spacing [m]
dz=((1 : n-1) * Dz / (n-1));
dzt=[0 dz];
alpha = dzt / (dzt(n));
delta = ((1-alpha)*delta1)+(alpha*deltan);
del_z = dz(n-1)-dz(n-2);


DROI= 4*D2; %diam of obs-plane region of interest [m]
%D1 = lam*Dz / DROI;
R = Dz; % wavefront radius of curvature [m]
x = (-N/2 : N/2 - 1).*delta1;
[x1, y1] = meshgrid(x);
[theta1, r1] = cart2pol(x1, y1);
  
%sg = exp(-(x1/(0.47*N*d1)).^16) .* exp(-(y1/(0.47*N*d1)).^16);
sg = exp(-(x1/(0.47*N)).^16) .* exp(-(y1/(0.47*N)).^16);
sg = repmat(sg, [1 1 n]);

Cn2 = fliplr(HV_Cn2(dzt+elevation));
%Cn2(:) = mean(Cn2);

amp_U=zeros(N_final,N_final,n_lam,n_trial);
phz = zeros(N,N,n,n_lam);
unwrap=zeros(N_final,N_final,n_lam,n_trial);

%{
all_st_func = zeros(N_final/2,n_lam,n_trial);
all_corr_func = zeros(N_final/2, n_lam,n_trial);
all_var = zeros(1,n_lam,n_trial);
all_U_var = zeros(1,n_lam,n_trial);
intensity_var = zeros(1,n_lam,n_trial);
intensity_var_U = zeros(1,n_lam,n_trial);
st_func = zeros(N_final/2,n_lam,n_trial);
an_corr = zeros(N_final/2,n_lam,n_trial);
%}
r0i = zeros(n, n_lam);

trial_seeds= abs(randperm(n*n_trial,n*n_trial));

for frame = 1:n_trial
    frame
    for l = 1:n_lam

    %%%%%%%%%atmospheric properties%%%%%%%%%
    %r0 = 0.1; %[m]
    %Cn2 = ((r0)^(-5.0/3))/( 0.423*(k^2)*Dz );
        

        %Cn2 = rytov/ ( 0.312*(k(l).^(7.0/6))*(Dz.^(11.0/6)) );
    %r0sw = (0.423 * (k.^2) * Cn2 * (3.0/8) * Dz)**(-3.0/5) %spherical wave

        r0pw = (0.423 * (k(l).^2) * (cos(z)^-1) * Cn2 * del_z).^(-3.0/5); %planewave
        r0i(:,l) = r0pw;
        %r0pw
        %clear U sg t corr_func t_st_func phz xn yn
        U = ones(N);
        %rng(trial_seed(frame))
        %rng(frame);
        %a="here";
        seeds = trial_seeds(((frame-1)*n)+1:(frame*n));
        for idx = 1 : n
            %seed=abs(round(100*randn(1)));
            seed=seeds(idx);
            phz(:,:,idx,l) = instance_phase_screen(r0pw(idx), N, delta(idx), L0(idx), l0(idx),seed);
        end
        arg = D2/(lam(l)*Dz);
        %x = linspace(-N/2,N/2,N);
        %[X, Y] = meshgrid(x);
        %X = X./(lam(l)*Dz/(2*D2)/deltan); Y = Y./(lam(l)*Dz/(2*D2)/deltan);
        %[~, r1] = cart2pol(x1, y1);
        A = lam(l) * Dz / deltan;
        
        D1 = lam(l)*Dz / DROI;
        pt = exp(-1i*k(l)/(2*Dz) * r1.^2) / D1^2 .* sinc(x1/D1) .* sinc(y1/D1) .* exp(-(r1/(4*D1)).^2);

        %t = sg.*exp(1i.*phz(:,:,:,l));
        %t(:,:,1) = ones(N);
        %t = cat(3,ones(N), t);
        [xn, yn, U] = ang_spec_multi_prop(U, lam(l), delta1, deltan, dz, sg.*exp(1i.*phz(:,:,:,l)));
        %U = U .* exp(-1i*pi/(lam(l)*R)*(xn.^2+yn.^2)); %collimate beam
        amp_U(:,:,l,frame) = abs(U((N/2)-(N_final/2)+1:(N/2)+(N_final/2),(N/2)-(N_final/2)+1:(N/2)+(N_final/2)));
        intensity = amp_U.^2;
        %intensity_var_U(l) = variance(intensity);
        %
        %warning('off', 'Images:initSize:adjustingMag');
        %intensity = mat2gray(intensity);
        %imshow(intensity);
        %image(intensity,'CDataMapping','scaled');
        %colorbar;
        if N ~= N_final
            mask = rect(N,N_final);
            holder = GoldsteinUnwrap(U,mask);
        else
            holder = GoldsteinUnwrap(U);
        end
        a="here2";
        holder(isnan(holder))=0;
        unwrap(:,:,l,frame)=holder((N/2)-(N_final/2)+1:(N/2)+(N_final/2),(N/2)-(N_final/2)+1:(N/2)+(N_final/2));
        %clear holder
        %imagesc(unwrap), colormap(jet), colorbar, axis square, axis off, title('GS Unwrapped phase');
        %unwrap = normalize(unwrap);
        %unwrapped_intensity = abs(unwrap(:,:,l,frame)).^2;
        %intensity_var(:,l,frame) = variance(unwrapped_intensity);
        %variance(unwrapped_intensity)

        %%%%%%%simulated%%%%%%%
        %corr_func = corr2_ft(unwrap(:,:,l,frame), unwrap(:,:,l,frame), delta1);
        %t_st_func = phz_st_func(corr_func);

%{
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
%}
        %clear U unwrap 

    end
end
%an_corr =  mean(an_corr,3);
%st_func =  mean(st_func,3);

%avg_st_func = mean(all_st_func,3);
%avg_corr_func = mean(all_corr_func,3);
%avg_var = mean(all_var(:));

p = linspace(0,Dz,n);
rytov = zeros(1,n_lam);
iso_theta = zeros(1,n_lam);
for l = 1:n_lam
    rytov(l) =  0.563 * k(l)^(7/6) * Dz* sum( Cn2 .* (1-p./Dz).^(5/6) .* (p(2)-p(1)));
    iso_theta(l) =  (2.91 * k(l)^(2) * sum( Cn2 .* (p).^(5/3) .* (p(2)-p(1))))^(-3/5);
end
iso_theta = rad2deg(iso_theta).*3600; %arc-seconds
r0 = (0.423 *(cos(z)^-1) * (k.^2) * sum(Cn2(1:n)*del_z)).^(-3.0/5);
%r0_n = ( 0.423 * (k^2) .*Cn2 .*del_z ).^(-3.0/5.0);
%r0pw = (0.423 * k^2 * Cn2 * Dz)^(-3/5);

%1.03*(D2/r0)^(5.0/3.0);

%{
a = all_corr_func;
an = an_corr;
for l = 1:n_lam
    a(:,l) = a(:,l)./mean(a(:,l).^2);
    an(:,l) = an(:,l)./mean(an(:,l).^2);
    %st_func(:,l) = st_func(:,l)./15;
end
norm= max(a(:,11));
an_norm = max(an(:,11));
for l = 1:n_lam
    a(:,l) = a(:,l)./norm;
    an_corr(:,l) = an_corr(:,l)./an_norm;
end
%}

%{
x=1:(N_final/2);
figure
subplot(2,2,1)
plot(x,all_st_func(:,11),'b')
title("Simulated Str Func (400 nm)")
subplot(2,2,2)
plot(x,st_func(:,11),'b--')
title("Theory Str Func (400 nm)")
subplot(2,2,3)
plot(x,all_st_func(:,1),'r')
title("Simulated Str Func (1000 nm)")
subplot(2,2,4)
plot(x,st_func(:,1),'r--')
title("Theory Str Func (1000 nm)")
figure
subplot(2,2,1)
plot(x,all_corr_func(:,11),'b')
title("Simulated Corr Func (400 nm)")
subplot(2,2,2)
plot(x,an_corr(:,11),'b--')
title("Theory Corr Func (400 nm)")
subplot(2,2,3)
plot(x,all_corr_func(:,1),'r')
title("Simulated Corr Func (460 nm)")
subplot(2,2,4)
plot(x,an_corr(:,1),'r--',x,an_corr(:,11),'b--')
title("Theory Corr Func (460 nm)")
%}

%n = s^2 / lam*z;

%figure
%plot(fliplr(dz),Cn2,'b')
%imagesc(unwrap), colormap(gray), colorbar, axis square, axis off, title('GS Unwrapped phase');

fitswrite(unwrap,'phases_DR0_32_nl_40_np_11_2.fits');fitswrite(amp_U,'amplitude_DR0_32_nl_40_np_11_2.fits')

