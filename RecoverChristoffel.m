% % Roxana Hess, September 2017

% Recover the atomic representing measure of the moments found in Step 1
% using Christoffel polynomials
% Among the two options, meaning minimizing the Riesz functional of the
% Christoffel polynomial or the trace of the moment matrix, we choose
% whichever works better for the given example

% output
% pts ... points on which the representing dirac measure is supported
% ChP ... matrix for plotting the Christoffel polynomial

function [pts,ChP] = RecoverChristoffel(expl,d,q,Mr)

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

% Compute polynomial p* = s(d) - Christoffel
nM = nchoosek(n+d,n);
pow = genpow(n+1,d); e = ones(length(pow),1); 
v = e;
for i = 1 : n
    v = v.*(x(i)*e).^pow(:,i+1); % vector of monomials
end
% Ch = nM - v'*inv(Mr)*v; % Christoffel polynomial
% Since the computation of the inverse may cause numerical problems, we 
% calculate the Christoffel polynomial via the orthonormal polynomials 
% associated with our measure.
if q == 1
    Ch = trace(Mr)-v'*v;
elseif q == 0
    Ch = nM;
    for i = 1 : nM
        [P,c] = OrthPoly(Mr(1:i,1:i));
        Ch = Ch - c*(P*v(1:i))^2;
    end
end

% relaxation order
% k = half degree for relaxation
if n == 1, k = d+1;
elseif n == 2 && d == 3, k = d+5;
elseif n == 3 && d == 1, k = d+2;
else k = d+3;
end

% Solve problem

if n == 3
    R = msdp(min(Ch), spt, k); % minimizing the Riesz functional of
    % the Christoffel polynomial
else
    R = msdp(mom(Ch)==0, spt, k); % minimizing the trace
end
msol(R); % solve problem

% output: support of the dirac
pts = [];
try
    pts = double(x);
    pts = reshape(pts,n,size(pts,3));
catch 
    warning('Probably no solution extracted')
end

% For plotting the Christoffel polynomial
if n == 1
    X = linspace(-1,1,200); %#ok
    ChP = eval(vectorize(Ch,'X'));
elseif n == 2
    [X1,X2] = meshgrid(linspace(-1,1,100)); %#ok
    ChP = eval(vectorize(Ch,'X1','X2'));
elseif n == 3
    ChP = [];
end

end