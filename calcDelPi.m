function Dp = calcDelPi(A,B)

n = size(A,1);
Dp = zeros(n);

for i = 1:n
    for j = 1:n
        partSum = 0;
        for p = 1:n
            partSum = partSum + (A(p,i)-A(p,j))*(B(p,j)-B(p,i)) ...
                + (A(i,p)-A(j,p))*(B(j,p)-B(i,p));
        end
        Dp(i,j) = partSum + (A(i,i)-A(j,j))*(B(j,j)-B(i,i)) ...
            + (A(i,j)-A(j,i))*(B(j,i)-B(i,j));
    end
end