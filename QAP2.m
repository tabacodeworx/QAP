% Tabu search for the Qaudratic Assignment Problem
%% AUTHOR: FEYI OLALOTITI-LAWAL
close all; clc; clear all;

% load case14
load case48
n = size(P,1);          % Graph size
permOrder = (1:n)';
T = n;              % Tabu length
As = 5*n^2;        % Aspiration length
maxit = 200;    % Maximum iterations

% Calculate distance matrix and position matrix Z
B = calcDist(P);
Z = zeros(n);

% Initial objective function
OF = sum(sum(A.*B));
DpVec = [];
OF_vec = [];
Dp = zeros(size(A));
for i = 1:n-1
    for j = i+1:n
        P2 = P; P2(i,:) = P(j,:); P2(j,:) = P(i,:);
        B2 = calcDist(P2);
        Dp(i,j) = sum(sum(A.*B2))-sum(sum(A.*B));
    end
end
Dp = Dp+Dp';      
[val ind] = sort(Dp(:)); jOpt = rem(ind(1),n); if jOpt==0; jOpt = n; end
iOpt = 1 + (ind(1)-jOpt)/n;
Z(iOpt,jOpt) = Z(iOpt,jOpt) + 1; Z(jOpt,iOpt) = Z(jOpt,iOpt) + 1;
DpOpt = val(1); DpVec = [DpVec; DpOpt];
OFvec = OF + DpOpt; DpBest = DpOpt;
itr = 0; 
tabu = Z~=0 & Z+T >= itr & Z'+T >= itr;
aspiratn = Dp <= DpBest | max(max([Z+A Z'+A]))<=itr;
PUp = P; PUp(iOpt,:) = P(jOpt,:); PUp(jOpt,:) = P(iOpt,:);
temp = permOrder; temp(iOpt) = permOrder(jOpt); temp(jOpt) = permOrder(iOpt); 
permOrder = temp; clear temp
BUp = calcDist(PUp);
P1 = PUp;

%Visualize Initial condition
figure(1); hold on;
for ii = 1:n-1; for jj = ii+1:n; coord = [P(ii,:);P(jj,:)];plot(coord(:,1),coord(:,2));end; end
for ii = 1:n 
    plot(P(ii,1),P(ii,2),'o','MarkerSize',18,'MarkerFaceColor','c','MarkerEdgeColor','k')
    text(P(ii,1)-40,P(ii,2), num2str(permOrder(ii)));
end
axis tight;axis off

% Begin Search
for itr = 2:maxit
    Dp = zeros(size(A));
    for i = 1:n-1
        for j = i+1:n
            P2 = P1; P2(i,:) = P1(j,:); P2(j,:) = P1(i,:);
            B2 = calcDist(P2);
            Dp(i,j) = sum(sum(A.*B2))-sum(sum(A.*BUp));
        end
    end
    Dp = Dp+Dp';
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
    check = ceil(maxit/10); if itr>check+1 && length(find(DpVec(itr-check:itr)~=0)) <=3; break; end
    OFvec = [OFvec; OFvec(itr-1) + DpOpt];
    if DpOpt <= DpBest; DpBest = DpOpt; end
    
    PUp = P1; PUp(iOpt,:) = P1(jOpt,:); PUp(jOpt,:) = P1(iOpt,:);
    
    temp = permOrder; temp(iOpt) = permOrder(jOpt); temp(jOpt) = permOrder(iOpt);
    permOrder = temp; clear temp
    
    BUp = calcDist(PUp);
    P1 = PUp;
    
end

%Visualize Final condition
figure(2); hold on;
for ii = 1:n-1 
    for jj = ii+1:n 
        coord = [P(ii,:);P(jj,:)];
        plot(coord(:,1),coord(:,2));
    end
end
for ii = 1:n 
    plot(P(ii,1),P(ii,2),'o','MarkerSize',18,'MarkerFaceColor','c','MarkerEdgeColor','k')
    text(P(ii,1)-40,P(ii,2), num2str(permOrder(ii)));
end
axis tight;axis off

%Objective Function
figure(3); plot(OFvec,'o-'); xlabel('Iterations'); ylabel('Objective Function');