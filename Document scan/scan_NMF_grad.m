close all;
clc;

fileID = fopen('Scan recto verso - NMF gradient.txt', 'a');
fprintf(fileID, '*********************************************************\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X1, map1] = imread('image_mix1.bmp');
[X2, map2] = imread('image_mix2.bmp');

x1c = double(X1(:));
x2c = double(X2(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = [x1c'; x2c'];
n = min(size(V));

mu = 1e-7;

W = rand(n,n) + 1;
H = rand(size(V)) + 1;

iter = 0;
err = 1;

figure
tic
while err > 1e-6
    v = W*H;
    
    gradW = -(V - v) * H';
    W = W - (mu*gradW);
    W = max(W,eps);
    
    gradH = -W' * (V-v);
    H = H - (mu*gradH);
    H = max(H,eps);

    e = ((V-W*H).^2)./2;
    err = sum(sum(e));
    
    if iter>300
        break
    end
    
    iter = iter+1
    
    fprintf (fileID, 'iteration %d\terreur = %e\n', iter, err);
    
    plot(iter, err, '*', 'MarkerSize', 10)
    hold on
    
end

fprintf (fileID, 'temps de calcul : %f\n', toc);

e1 = H(1,:);
e1 = (e1-min(e1))./(max(e1)-min(e1));       % rescale estimation n°1
e1 = reshape(e1,hX,[]);
imwrite(e1, 'estimate_1.png');

e2 = H(2,:);
e2 = (e2-min(e2))./(max(e2)-min(e2));       % rescale
e2 = reshape(e2,hX,[]);                     % reshape in matric
imwrite(e2, 'estimate_2.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% les données sont en format double, l'image est en grayscale
% imwrite prend les données de 0 à 1 et les étale sur 0-255
% les valeurs négatives sont mises à 0
% le valeurs > 1 sont mises à 1

sizeR = size(X1);
hX = sizeR(1); hY = sizeR(2);

l = [1:hX]';    l = repmat(l,1,hY);
c = [1:hY];     c = repmat(c,hX,1);

e1 = H(1,:);                                % estimation
e1 = e1 + min(e1);                          % make all values positive
e1 = (e1-min(e1))./(max(e1)-min(e1));       % rescale
e1 = reshape(e1,hX,[]);                     % reshape in matric
%e1 = fliplr(e1);
%e1 = imcomplement(e1);
imwrite(e1, 'estimate_1.png');

e2 = H(2,:);                                % estimation
e2 = e2 + min(e2);                          % make all values positive
e2 = (e2-min(e2))./(max(e2)-min(e2));       % rescale
e2 = reshape(e2,hX,[]);                     % reshape in matric
imwrite(e2, 'estimate_2.png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%