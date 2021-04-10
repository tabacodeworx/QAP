function B = calcDist(P)

n = size(P,1);

for i=1:n
    for j=1:n
        B(i,j)=sqrt((P(i,1)-P(j,1))^2 + (P(i,2)-P(j,2))^2);
    end
end