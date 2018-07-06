% Roxana Hess, September 2017

% Find the truncated moment sequence of the approximate optimal
% design by solving the SDP problem resulting from relaxing the
% ideal problem on measures with the Lasserre hierarchy

% output
% Mr   ... moment matrix
% momv ... moment vector

function [Mr, momvd] = SDPApproxOptDesign(expl,d,q)

% Define dimension n
if expl == 1, n = 1;
elseif expl == 2 || expl == 3 || expl == 4 || expl == 5
    n = 2;
elseif expl == 6, n = 3;
end

mpol('x',n)
nM = nchoosek(n+d,n);

% Support of the measure (design space)

% Univariate unit interval
if expl == 1
    spt = 1-x^2>=0;

% Wynn's polygon
elseif expl == 2
    % polygon scaled to fit the unit circle
    spt = [x(1) >= -1/2/sqrt(2), x(2) >= -1/2/sqrt(2), ...
        x(1)^2+x(2)^2 <= 1, ... % Need it, otherwise problem maybe unbounded
        x(1) <= 1/3*(x(2)+4/2/sqrt(2)), x(2) <= 1/3*(x(1)+4/2/sqrt(2))];
    
% Ellipse with hole
elseif expl == 3
    % Bigger ellipse with a smaller ellipse as a hole
    spt = [9*x(1)^2 + 13*x(2)^2 <= 7.3,...
        5*x(1)^2 + 13*x(2)^2 >= 2];
    
% Moon
elseif expl == 4
    % Bigger disk without an intersecting smaller disk
    spt = [(x(1)+.2)^2 + x(2)^2 <= .36,...
        (x(1)-.6)^2 + x(2)^2 >= .16];

% Folium
elseif expl == 5
    spt = [-x(1)*(x(1)^2-2*x(2)^2)-(x(1)^2+x(2)^2)^2 >= 0,...
        x(1)^2+x(2)^2 <= 1]; % Need it, otherwise problem maybe unbounded

% The 3-dimensional unit sphere
elseif expl == 6
    spt = x'*x==1;
    
end

% relaxation order
% k = half degree for relaxation
if n == 2, k = d+3;
else k = d; end

% solve SDP problem
if q == 0 % D-optimal design
    
    P = msdp(spt,k); % model problem with Gloptipoly 3
    [F,~,momv] = myalmip(P); % give it to YALMIP
    
    if (n == 1 && d == 1) || n == 3
        momM = sdpvar(F(2)); % moment matrix for relaxation
    % Attention: YALMIP orders matrices in F according to type of constraint
    else
        momM = sdpvar(F(1)); % moment matrix for relaxation
    end
    regM = momM(1:nM,1:nM); % moment matrix for regression
    f = -geomean(regM); % objective function to minimize (D-optimal)
    % Attention: geomean overloaded function in YALMIP (geomean of eig for
    % sdpvar)
    %f = -trace(regM); % objective function to minimize (T-optimal, q = 1)
    optimize(F,f,sdpsettings('solver','sedumi')); %solve problem
    
%     % For fixing moments
% 
%     Sigma = [1 .01 .95; .01 2 0; .95 0 1];
%     u = (2*pi)^(3/2)*sqrt(det(Sigma));
%     S = Sigma./u;
%     optimize([F;momv(7)==S(2,2);momv(5)==S(1,2);momv(6)==S(1,3);momv(9)==S(3,3)],f,sdpsettings('solver','sedumi')); %solve problem
% 
%     Sigma = [1 .1 .95; .1 1 0; .95 0 1];
%     u = (2*pi)^(3/2)*sqrt(det(Sigma));
%     S = Sigma./u;
%     optimize([F;momv(5)==S(1,2);momv(6)==S(1,3)],f,sdpsettings('solver','sedumi')); %solve problem
    
    momvd = [1;double(momv)]; % vector of moments
    M = double(momM);

%     % Doing it by hand (without gloptipoly) for experimental reasons
%
%     sdpvar x1 x2;
%     
%     % Psatz multipliers
%     [q0,q0c] = polynomial([x1 x2],k);
%     [q1,q1c] = polynomial([x1 x2],k-2);
%     [q2,q2c] = polynomial([x1 x2],k-2);
%     [q3,q3c] = polynomial([x1 x2],k-2);
%     [q4,q4c] = polynomial([x1 x2],k-2);
%     [q5,q5c] = polynomial([x1 x2],k-2);
%     
%     sol = solvesos([sos(trace(M)-1-x1^2-x2^2-q0-q1*(x1+1/2/sqrt(2))-q2*(x2+1/2/sqrt(2))-q3*(1-x1^2-x2^2)-q4*(1/3*(x2+4/2/sqrt(2))-x1)-q5*(1/3*(x1+4/2/sqrt(2))-x2));...
%             sos(q0);sos(q1);sos(q2);sos(q3);sos(q4);sos(q5)],...
%             0 , [], [momv;q0c;q1c;q2c;q3c;q4c;q5c]);
%     p = sol.problem;
        
elseif q == 1 % T-optimal design
    
    v = mmon(x,d);
    f = -v'*v; % objective function to minimize
    
    P = msdp(min(f),spt,k); % model problem with Gloptipoly 3
    msol(P); % solve problem
    
    momvd = double(mvec(meas)); % moment vector
    M = double(mmat(meas)); % moment matrix
    
end
Mr = M(1:nM,1:nM);

end