% Tabu search for the Qaudratic Assignment Problem

close all; clc; clear all;

% Set paramters
raw = graph1();              % Get initial locations.
P = raw(:,2:3);
n = size(P,1);          % Graph size
permOrder = (1:n)';
T = n;              % Tabu length
As = 5*n^2;        % Aspiration length
maxit = 400;    % Maximum iterations

%Make arbitrary flow matrix
A = zeros(n); max_aij = 20;
for i = 1:n-1
    A(i,i+1:n) = randi(max_aij, 1, n-i); 
    A(i+1:n,i) = A(i,i+1:n)';
end

% Calculate distance matrix and position matrix Z
load case14
B = calcDist(P);
Z = zeros(n);

% Initial objective function
OF = 0.5*sum(sum(A.*B));
DpVec = [];

Dp = calcDelPi(A,B);

[val ind] = sort(Dp(:));
jOpt = rem(ind(1),n); if jOpt==0; jOpt = n; end
iOpt = 1 + (ind(1)-jOpt)/n;

Z(iOpt,jOpt) = Z(iOpt,jOpt) + 1;
Z(jOpt,iOpt) = Z(jOpt,iOpt) + 1;

DpOpt = val(1);
DpVec = [DpVec; DpOpt];
DpBest = DpOpt;

itr = 0;
tabu = Z~=0 & Z+T >= itr & Z'+T >= itr;
aspiratn = Dp <= DpBest | max(max([Z+A Z'+A]))<=itr;

PUp = P; PUp(iOpt,:) = P(jOpt,:); PUp(jOpt,:) = P(iOpt,:);

temp = permOrder; temp(iOpt) = permOrder(jOpt); temp(jOpt) = permOrder(iOpt); 
permOrder = temp; clear temp

BUp = calcDist(PUp);
P1 = PUp;


% Begin Search
for itr = 2:maxit
   i = iOpt; j = jOpt;
   DpUpd = zeros(n);
   for u = 1:n-1 
       for v = u+1:n
           if i==u || i==v || j==u || j==v 
               continue
           else
               DpUpd(u,v) = (A(u,i)-A(u,j)+A(v,j)-A(v,i))*(BUp(v,i)-BUp(v,j)+BUp(u,j)-BUp(u,i)) ...
                  +(A(i,u)-A(j,u)+A(j,v)-A(i,v))*(BUp(i,v)-BUp(j,v)+BUp(j,u)-BUp(i,u));
            end
        end
    end
%     Dp = DpOpt + (DpUpd+DpUpd');
    Dp = (DpUpd+DpUpd');
    Dp_fsble = 10*ones(size(Dp));
    Dp_fsble(~tabu) = Dp(~tabu);
    Dp_fsble(aspiratn) = Dp(aspiratn);
    [val ind] = sort(Dp_fsble(:));
    
    rejectMove = 1; count = 0;
    
    while rejectMove
        count = count+1;      
        jtemp = rem(ind(count),n); if jtemp==0; jtemp = n; end
        itemp = 1 + (ind(count)-jtemp)/n;       
        rejectMove = tabu(itemp,jtemp)==1 && aspiratn(itemp,jtemp)==0;       
    end
    
    iOpt = itemp; jOpt = jtemp;
    Z(iOpt,jOpt) = itr;
    Z(jOpt,iOpt) = itr;  
    tabu = Z~=0 & Z+T >= itr & Z'+T >= itr;
    aspiratn = Dp <= DpBest | max(max([Z+A Z'+A]))<=itr;
    
    DpOpt = val(count);
    DpVec = [DpVec; DpOpt];
    if DpOpt <= DpBest; DpBest = DpOpt; end

    PUp1 = P1; PUp(iOpt,:) = P1(jOpt,:); PUp(jOpt,:) = P1(iOpt,:);

    temp = permOrder; temp(iOpt) = permOrder(jOpt); temp(jOpt) = permOrder(iOpt); 
    permOrder = temp; clear temp

    BUp = calcDist(PUp); 
    P1 = PUp;
   
end  
plot(DpVec)