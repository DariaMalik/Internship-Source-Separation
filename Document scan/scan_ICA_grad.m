clear all;
close all;
clc;

fileID = fopen('Traitement du scan recto-verso.txt', 'w');
fprintf(fileID, '*********************************************************\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X1, map1] = imread('image_mix1.bmp');
[X2, map2] = imread('image_mix2.bmp');

x1c = double(X1(:));
x2c = double(X2(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1 = x1c - mean(x1c);
x2 = x2c - mean(x2c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID, '\nOriginal data \t\t\t\t Centered data');
fprintf(fileID, '\n************************ \t\t ************************\n');
fprintf(fileID, '\nx1 \t\t x2 \t\t\t x1 \t\t x2\n');
fprintf(fileID, '%d \t\t %d \t\t\t %.2f \t\t %.2f\n', [x1c x2c x1 x2]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = [x1 x2];
% size_M = size(M)

C = cov(M);                         % matrice de covariance des observ. centrées
[E, D] = eig(C);                    % val et vect propres de matrice de cov

V = D^(-1/2) * E';                  % matrice de whitening

M = V*M';                           % whitening
C = cov(M');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M1 = M;
M = M';

mu = 0.01;                          % "pas" de gradient

w1 = rand(2,1);                     % vecteur initial

for iter=1:100
    
    y = w1'*M1;
     
    y3 = y.^3;
        
    grad = 4 * sign(kurtosis(y)) * ( mean([ M(:,1).*y3', M(:,2).*y3']) - (3*w1') );
    
    w1 = w1 + (mu*grad');
    w1 = w1 / norm(w1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w2 = [ w1(2); -w1(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% les données sont en format double, l'image est en grayscale
% imwrite prend les données de 0 à 1 et les étale sur 0-255
% les valeurs négatives sont mises à 0
% le valeurs > 1 sont mises à 1

sizeR = size(X1);   hX = sizeR(1);  hY = sizeR(2);

l = [1:hX]';    l = repmat(l,1,hY);
c = [1:hY];     c = repmat(c,hX,1);

e1 = w1'*M1;                                % estimation
e1 = e1 + min(e1);                          % make all values positive
e1 = (e1-min(e1))./(max(e1)-min(e1));       % rescale
e1 = reshape(e1,hX,[]);                     % reshape in matric
e1 = fliplr(e1);
e1 = imcomplement(e1);
imwrite(e1, 'estimate_1.png');

e2 = w2'*M1;                                % estimation
e2 = e2 + min(e2);                          % make all values positive
e2 = (e2-min(e2))./(max(e2)-min(e2));       % rescale
e2 = reshape(e2,hX,[]);                     % reshape in matric
imwrite(e2, 'estimate_2.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID, '*********************************************************\n');
fprintf(fileID, '\nEstimated data');
fprintf(fileID, '\n****************************\n');
fprintf(fileID, '\ne1 \t\t e2\n');
fprintf(fileID, '%8.2f\t%8.2f\n', [e1 e2]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%