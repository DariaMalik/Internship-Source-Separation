%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Daria MALIK, juillet-août 2019, stage de fin de L3 EEA CMI

% le programme detecte lui-même les fichiers .wav qui se trouve dans le
% dossier avec ce main.m et leur quantité, il compose lui-même la matrice
% X avec ces fichiers .wav consideré comme les observations

% dans le dossier avec ce fichier et les audios de mélange il faut créer un
% dossier "estimations" où le programme enregistrera automatiquelment tous
% les fichiers d'éstimations

% le fichier kurtosis.m est une fonction, il faudra donc l'avoir dans le
% dossier avec ce fichier main.m

close all;
clear all;
clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lecture autoimatique de tous les fichiers .wav qui se trouve dans le
% dossier et le stockage des échantillons dans la matrice X

folder_info = dir('*.wav');
n = length(folder_info);

F = zeros(1,n);

for j = 1 : n
  fileID = folder_info(j).name;
  [x, F(j)] = audioread(fileID);
  if j==1
      X = zeros(size(x));
      X = repmat(X,1,n);
  end
  X(:,j) = x;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% centering and whitening

% centering
for j=1:n
    X(:,j) = X(:,j) - mean(X(:,j));
end

C = cov(X);                         % matrice de covariance des observ. centrées
[E, D] = eig(C);                    % val et vect propres de matrice de cov

V = D^(-1/2) * E';                  % matrice de whitening

X = V*X';                           % whitening
C = cov(X');                        % nouvelle matrice de covariance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% algorithme du gradient

mu = 0.01;                          % "pas" de gradient

w = ones(n,n);

for p=1:n
    wp = rand (n,1);
    
    if p==1
        for iter=1:200
            y = wp' * X;
            y3 = repmat((y.^3)',1,n);
            
            X1 = X'.*y3;
            
            grad = 4 * sign(kurtosis(y)) * ( mean(X1) - (3*wp') );
            
            wp = wp + (mu*grad');
            wp = wp / norm(wp);
        end
    end

    if p>1
        for iter = 1:200
            y = wp' * X;
            y3 = repmat((y.^3)',1,n);
            
            X1 = X'.*y3;
            
            grad = 4 * sign(kurtosis(y)) * ( mean(X1) - (3*wp') );
            
            wp = wp + (mu*grad');
            
            sigma = 0;
            for j=1:p-1
                sigma = sigma + ((wp'*w(:,j))*w(:,j));
            end
            wp = wp - sigma;
            wp = wp /  norm(wp);
        end
    end
    
    w(:,p) = wp;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for j=1:n
    estim = w(:,j)'*X;
    estim = estim/max(estim);
    fileID = sprintf('estimations/estimation%d.wav',j);
    audiowrite(fileID, estim, F(1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%