close all;
clc;

fileID = fopen('Scan recto verso - NMF ALS.txt', 'a');
fprintf(fileID, '*********************************************************\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X1, map1] = imread('image_mix1.bmp');
[X2, map2] = imread('image_mix2.bmp');

x1c = double(X1(:));
x2c = double(X2(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = [x1c'; x2c'];
n = min(size(V));

W = rand(n,n);

iter = 0;
err = 1;

figure
tic
while err > 1e-6

    H = inv(W'*W)*W'*V;
    H = max(H,eps);
    
    W = V*H'*inv(H*H');
    W = max(W,eps);
    
    e = ((V-W*H).^2)./2;
    err = sum(sum(e))
    
    if iter>400
        break
    end
    
    iter = iter+1;
    
end
fprintf (fileID, 'temps de calcul : %f\n', toc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sizeR = size(X1);
hX = sizeR(1); hY = sizeR(2);

e1 = H(1,:);
e1 = (e1-min(e1))./(max(e1)-min(e1));       % rescale estimation n°1
e1 = reshape(e1,hX,[]);
imwrite(e1, 'estimate_1.png');

e2 = H(2,:);
e2 = (e2-min(e2))./(max(e2)-min(e2));       % rescale
e2 = reshape(e2,hX,[]);                     % reshape in matric
imwrite(e2, 'estimate_2.png');
%%
% les données sont en format double, l'image est en grayscale
% imwrite prend les données de 0 à 1 et les étale sur 0-255
% les valeurs négatives sont mises à 0
% le valeurs > 1 sont mises à 1



l = [1:hX]';    l = repmat(l,1,hY);
c = [1:hY];     c = repmat(c,hX,1);

e1 = H(1,:);                                % estimation
% e1 = e1 + min(e1);                          % make all values positive
e1 = (e1-min(e1))./(max(e1)-min(e1));       % rescale
e1 = reshape(e1,hX,[]);                     % reshape in matric
%e1 = fliplr(e1);
%e1 = imcomplement(e1);
imwrite(e1, 'estimate_1.png');

e2 = H(2,:);                                % estimation
% e2 = e2 + min(e2);                          % make all values positive
e2 = (e2-min(e2))./(max(e2)-min(e2));       % rescale
e2 = reshape(e2,hX,[]);                     % reshape in matric
imwrite(e2, 'estimate_2.png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%