clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vraies sources pour comparaison
obj2=cell2mat(fitsread('obj2_amplified_spectrum.spc.fits','BinTable'));
obj50=cell2mat(fitsread('obj50_amplified_spectrum.spc.fits','BinTable'));

figure
plot(obj2(:,1),obj2(:,2))
hold on
plot(obj50(:,1),obj50(:,2))
legend('Objet 2', 'Objet 50')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% lire l'image source
im = fitsread('myCMC_config_2sources_fastica_NISP_GRED_0.0_10_IMG.fits','image');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% afficher l'image source
figure
imagesc(im)
title('2 sources melangées dans 1 direction')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% isoler l'objet et les lignes à étudier
[x,y]=find(im==max(max(im)));

figure
imagesc(im(x-5:x+5, y-100+17:y+550-8))
title('2 sources melangées dans 1 direction')

x1 = im(x,y-100:y+550);
x2 = im(x+1,y-100:y+550);

figure
subplot(211)
imagesc(x1)
title('Spectre 2D - ligne/observation 1')
subplot(212)
imagesc(x2)
title('Spectre 2D - ligne/observation 2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% passer en spectre 1D
figure
plot(x1)
title('Spectres 1D melangés')
hold on
plot(x2)
legend('Observation 1', 'Observation 2')
axis tight


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% preliminaires pour ICA
x1 = x1 - mean(x1);
x2 = x2 - mean(x2);

X = [x1' x2'];

X = whitening(X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% estimations

w = min_kurt(X);

E = w'*X;

if abs(min(E(1,:)))>abs(max(E(1,:)))                        % les rendre positives
    E(1,:) = -E(1,:);
end
if abs(min(E(2,:)))>abs(max(E(2,:)))
    E(2,:) = E(2,:)*(-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% caler les raies obtenues sur les vrais spectres

s2 = obj2(:,2)';
s50 = obj50(:,2)';

l2 = obj2(:,1);
l50 = obj50(:,1);

e1 = E(1,:);
e2 = E(2,:);

e1 = e1 + abs(min(e1));
e2 = e2 + abs(min(e2));


x2 = find(s2==max(s2));
x50 = find(s50==max(s50));
xe1 = find(e1==max(e1));
xe2 = find(e2==max(e2));

if max(e1)>max(e2)
    pix2 = l2(x2-xe1+1 : x2+(length(e1)-xe1));
    pix50 = l50(x50-xe2+1 : x50+(length(e2)-xe2));
    e22 = e1;
    e50 = e2;
else
    pix2 = l2(x2-xe2+1 : x2+(length(e2)-xe2));
    pix50 = l50(x50-xe1+1 : x50+(length(e1)-xe1));
    e22 = e2;
    e50 = e1;
end

figure
subplot(211)%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis left
plot(l2, s2)
ylabel ('Source 1')
hold on
yyaxis right
plot(pix2, e22)
ylabel('Estimation 1')
subplot (212)%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis left
plot(l50, s50)
ylabel ('Source2')
hold on
yyaxis right
plot(pix50, e50)
ylabel('Estimation 2')

s2 = s2/max(s2);
s50 = s50/max(s50);

e22 = e22/max(e22);
e50 = e50/max(e50);

pix21 = linspace(14000, l2(x2), find(e22==max(e22)));
pix22 = linspace(l2(x2), 20080, length(e22)-find(e22==max(e22))+1);
pix2 = [pix21 pix22(2:end)];

pix501 = linspace(14850, l50(x50), find(e50==max(e50)));
pix502 = linspace(l50(x50), 20700, length(e50)-find(e50==max(e50))+1);
pix50 = [pix501 pix502(2:end)];

figure
subplot(211)%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis left
plot(l2, s2, '*r', 'MarkerSize', 4)
ylabel ('Source 1')
hold on
yyaxis right
plot(pix2, e22, 'k', 'linewidth', 3)
ylabel('Estimation 1')
legend('Source 1', 'Estimation de la source 1')
axis ([14000 20000 0 1])
subplot (212)%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis left
plot(l50, s50, '*r', 'MarkerSize', 4)
ylabel ('Source2')
hold on
yyaxis right
plot(pix50, e50, 'k', 'linewidth', 3)
ylabel('Estimation 2')
legend('Source 2', 'Estimation de la source 2')
axis ([14000 20000 0 1])
