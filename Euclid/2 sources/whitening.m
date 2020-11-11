function W = whitening(M)

C = cov(M);                         % matrice de covariance des observ. centrées
[E, D] = eig(C);                    % val et vect propres de matrice de cov

V = D^(-1/2) * E';                  % matrice de whitening

W = V*M';                           % whitening