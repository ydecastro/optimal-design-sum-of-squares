% Roxana Hess, September 2017

% Recover the atomic representing measure of the moments found in Step 1
% using the method by Nie
% Instead of generating a random polynomial strictly positive on the
% support, we minimize the trace of the moment matrix

% output
% pts ... points on which the representing dirac measure is supported

function [pts] = RecoverNie(expl,d,momv)

% Define dimension n
if expl == 1, n = 1;
elseif expl == 2 || expl == 3 || expl == 4 || expl == 5
    n = 2;
elseif expl == 6, n = 3;
end

mpol('x',n)

% Support of the measure (design space)

% Univariate unit interval
if expl == 1
    n = 1;
    spt = 1-x^2>=0;

% Wynn's polygon
elseif expl == 2
    n = 2;
    % polygon scaled to fit the unit circle
    spt = [x(1) >= -1/2/sqrt(2), x(2) >= -1/2/sqrt(2), ...
        x(1)^2+x(2)^2 <= 1, ... % Need it, otherwise problem maybe unbounded
        x(1) <= 1/3*(x(2)+4/2/sqrt(2)), x(2) <= 1/3*(x(1)+4/2/sqrt(2))];
    
% Ellipse with hole
elseif expl == 3
    n = 2;
    % Bigger ellipse with a smaller ellipse as a hole
    spt = [9*x(1)^2 + 13*x(2)^2 <= 7.3,...
        5*x(1)^2 + 13*x(2)^2 >= 2];
    
% Moon
elseif expl == 4
    n = 2;
    % Bigger disk without an intersecting smaller disk
    spt = [(x(1)+.2)^2 + x(2)^2 <= .36,...
        (x(1)-.6)^2 + x(2)^2 >= .16];

% Folium
elseif expl == 5
    n = 2;
    spt = [-x(1)*(x(1)^2-2*x(2)^2)-(x(1)^2+x(2)^2)^2 >= 0,...
        x(1)^2+x(2)^2 <= 1]; % Need it, otherwise problem maybe unbounded

% The 3-dimensional unit sphere
elseif expl == 6
    n = 3;
    spt = x'*x==1;
    
end

% relaxation order
% k = half degree for relaxation
if n == 1, k = d+1;
elseif n == 3 && d == 3, k = d+4;
else k = d+3;
end

% Solve problem 
v = mmon(x,k);
f = v'*v; % objective function to minimize
%R = msdp(mom(mmon(x,0,2*d))-momv(1:nchoosek(2*d+n,n))==0, spt, k);
R = msdp(min(-f),0==mom(mmon(x,0,2*d))-momv(1:nchoosek(2*d+n,n)), spt, k); 
% model problem with Gloptipoly
% Attention: By default GloptiPoly minimizes the trace, when no objective
% function is found
mset('yalmip',false);
msol(R); % solve problem

% output: support of the dirac
pts = double(x);
pts = reshape(pts,n,size(pts,3));

end
