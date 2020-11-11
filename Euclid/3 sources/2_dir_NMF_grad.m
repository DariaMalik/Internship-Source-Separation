close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im0 = fitsread('myCMC_config_4sources_scn_vrai_ss_bruit_NISP_GRED_0.0_10_IMG.fits','image');
im90 = fitsread('myCMC_config_4sources_scn_vrai_ss_bruit_NISP_GRED_90.0_10_IMG.fits','image');

load simulation_4sources.mat

figure
subplot(121)
imagesc(im0)
title('0°')
subplot(122)
imagesc(im90)
title('90°')

figure
imagesc(im90')
title('Image originale, 90°')

im90 = im90';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X0, Y0] = find(im0==max(max(im0)));
[X90, Y90] = find(im90==max(max(im90(:,1:1170))));

x10 = im0(X0-1,600:1650);
% x20 = im0(X0,600:1650);

x190 = im90(600:1650,Y90+1)';
x290 = im90(600:1650,Y90)';

% figure
% imagesc([x10; x190; x290])

figure
plot (x10)
hold on
plot (x190)
plot (x290)
legend('x10', 'x190', 'x290')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = [x10; x190; x290];

iter = 1;

W = rand(3,3);
H = rand(size(V));

W(7) = 0;
W(5) = 0;
W(6) = 0;

mu = 0.0001;
err = 1;

while err > 1e-6
    v = W*H;
    
    gradW = -(V - v) * H';
    W = W - (mu*gradW);
    W = max(W,eps);
    
    W(7) = 0;
    W(5) = 0;
    W(6) = 0;

    gradH = -W' * (V-v);
    H = H - (mu*gradH);
    H = max(H,eps);
    
    e = ((V-W*H).^2)/2;
    err = sum(sum(e))
    
    iter = iter+1;
    
    if iter>5000
        break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Comparer les résultats

figure
plot(H(1,:))
hold on
plot(H(2,:))
plot(H(3,:))
axis tight

