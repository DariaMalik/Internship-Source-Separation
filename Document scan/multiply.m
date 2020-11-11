function C = multiply(A,B)

a = size(A);
b = size(B);

c = [a(1) b(2)];

for j=1:c(1)
    for k=1:c(2)
        C(j,k) = A(j,:)*B(:,k);
    end
end