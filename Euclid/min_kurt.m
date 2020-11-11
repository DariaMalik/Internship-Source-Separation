function w = min_kurt(X)

    mu = 0.01;                          % "pas" de gradient
    
    w = ones(min(size(X)));
    n = length(w);
    
    if n==2
        wp = rand (n,1);
        for iter=1:200
            y = wp' * X;
            y3 = repmat((y.^3)',1,n);

            X1 = X'.*y3;

            grad = 4 * sign(kurtosis(y)) * ( mean(X1) - (3*wp') );

            wp = wp + (mu*grad');
            wp = wp / norm(wp);
        end
        w(:,1) = wp;
        w(:,2) = [ wp(2); -wp(1)];
    end
    
    if n>2
        for p=1:n

            wp = rand (n,1);

            if p==1
                for iter=1:150
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
    end

    