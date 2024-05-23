hs = mk_hs();
d_lc = 4e-06; % model lc thickness
thetai = 2/180*pi*0; % light
theta0 = grating_angle(2, -1, thetai):0.005:grating_angle(2, +1, thetai); % outgoing angles
theta0d = [grating_angle(2, -1, thetai), grating_angle(4, -1, thetai), grating_angle(6, -1, thetai), ...
    grating_angle(6, 1, thetai), grating_angle(4, 1, thetai), grating_angle(2, 1, thetai)]; % outgoing angles
alpha = (0:1:90)/180*pi; % molecule rotation wrt slm plane
alpha0 = 2/180*pi; %pre-rotation 

dphase = pi*2*d_lc/hs.wvl_rgb(1)*(nref_lc(theta0, alpha)-nref_lc(0, pi/2));
dphased = pi*2*d_lc/hs.wvl_rgb(1)*(nref_lc(theta0d, alpha)-nref_lc(0, pi/2));
figure; surf(theta0, alpha, dphase);
xlabel('light ray angle, rad'); ylabel('LC angle, rad'); zlabel('phase modulation, rad');
dphasedi = zeros(size(dphase));
for i = 1:size(dphased,1)
    dphasedi(i, :) = interp1(theta0d, dphased(i,:), theta0, 'linear');
end
figure; surf(theta0, alpha, dphasedi);
xlabel('light ray angle, rad'); ylabel('LC angle, rad'); zlabel('phase modulation, rad');
disp(['Max interpolation error: ' num2str(max(abs(dphase-dphasedi),[], 'all'))]);
for i = 1:size(dphase,1)
 dphase(i,:) = dphase(i,:) - min(dphase(i,:));
end
for i = 1:size(dphased,1)
 dphased(i,:) = dphased(i,:) - min(dphased(i,:));
end
figure; surf(theta0, alpha, dphase);
xlabel('light ray angle, rad'); ylabel('LC angle, rad'); zlabel('phase modulation difference, rad');
hold on;
%surf(theta0d, alpha, dphased);
hold off;

figure;
hold on;
for i = 1:1:91
    dphase = pi*2*d_lc/hs.wvl_rgb(1)*(nref_lc(theta0, alpha(i))-nref_lc(0, pi/2));
    dphase = dphase - min(dphase(:));
    plot(theta0, dphase);
    dphased = pi*2*d_lc/hs.wvl_rgb(1)*(nref_lc(theta0d, alpha(i))-nref_lc(0, pi/2));
    %dphasedi = interp1(theta0d, dphased - min(dphased(:)), theta0, 'spline');
    dphasedi = interp1(theta0d, dphased, theta0, 'spline');
    dphasedi = dphasedi-min(dphasedi(:));
    plot(theta0,dphasedi);
end
hold off;

