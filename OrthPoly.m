% Roxana Hess, September 2017

% Compute orthonormal polynomials associated with a measure given by its
% truncated moment matrix M
% Function called in RecoverChristoffel.m to compute the Christoffel
% polynomials without the inverse

% output: 
% P       ... coefficient vector of the orthogonal polynomial, meaning
%             P*v is the orthogonal polynomial for v the monomial vector  
%             of correct length
% csquare ... square of the normalising coefficient, meaning
%             sqrt(csquare)*P*v is the orthonormal polynomial

function [P,csquare] = OrthPoly(Mr)

% coefficient vector
nM = length(Mr);
P = zeros(1,nM);
for i = 1 : nM
    Mi = [Mr(1:(nM-1),1:(i-1)), Mr(1:(nM-1),(i+1):nM)];
    P(i) = (-1)^(i+nM).*det(Mi);
end

% normalising coefficient
csquare = 1/sum(sum((P'*P).*Mr));
end