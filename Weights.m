% Roxana Hess, September 2017

% Calculate the weights of the dirac measure

% may cause problems and give stupid output (without giving warnings) as
% linear system might be overdifined or underdefined

function [w] = Weights(expl,d,pts,momv)

% Define dimension n
if expl == 1, n = 1;
elseif expl == 2 || expl == 3 || expl == 4 || expl == 5
    n = 2;
elseif expl == 6, n = 3;
end

% generate a matrix where each column is the monomial vector evaluated at
% a point of the support of the dirac
pow = genpow(n+1,2*d);
B = zeros(length(pow),size(pts,2));
for i = 1 : size(pts,2)
    l = ones(length(pow),1);
    for j = 1 : size(pts,1)
        l = l.*(pts(j,i)*ones(length(pow),1)).^pow(:,j+1);
    end
    B(:,i) = l;
end

% solve the system of linear equations B*w = momv
% B*w is the moment vector of the dirac, so the solution w is the vector of
% the weights of the dirac measure
w = linsolve(B,momv(1:nchoosek(2*d+n,n)))'; % potentially overdefined
%w = mldivide(B(1:size(pts,2),1:size(pts,2)),momv(1:size(pts,2)))'; % potentially underdefined
%w = (B\momv)';
end