% Roxana Hess, September 2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPROXIMATE OPTIMUM DESIGNS ON SEMI-ALGEBRAIC DESIGN SPACES %
% Examples                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Requires GloptiPoly3, YALMIP and SeDuMi 

clear all; close all; clc
mset clear

% set parameters

% Six examples: One univariate, four in two dimensions, and one three-dimensional
% example
% expl = 1: Univariate unit interval
% expl = 2: Wynn's polygon
% expl = 3: Ring of ellipses
% expl = 4: Moon
% expl = 5: Folium
% expl = 6: The 3-dimensional unit sphere
expl = 2;

% Regression order
% d = half degree of regression order
% Consider d = 1 up to 7 for the univariate problem and d = 1 up to 3 for 
% the two- and three-dimensional examples
d = 3;

% Kiefer's \phi_q-criteria
% q = 0: D-optimal design ( \phi_q(M) = log det M )
% q = 1: T-optimal design ( \phi_q(M) = trace(M) )
q = 0;

% Choose method for recovering the representing measure in Step 2
% recover = 0: Use method by Nie
% recover = 1: Use Christoffel polynomials - makes only sense for q = 0
recover = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: Find the truncated moment sequence of the approximate optimal
% design

[M, momv] = SDPApproxOptDesign(expl,d,q);

% Step 2: Find representing atomic measures of the moments found in Step 1

% Find the support of the dirac measure
if recover == 0
    pts = RecoverNie(expl,d,momv); % Use method by Nie
    Ch = [];
elseif recover == 1
    [pts,Ch] = RecoverChristoffel(expl,d,q,M); % Use Christoffel polynomials
end

% Calculate the weights of the dirac
w = Weights(expl,d,pts,momv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot solution
 
disp('Support of the dirac measure and corresponding weights:')
pts %#ok
w %#ok
PlotSolution(expl,d,q,pts,w,Ch,M);
% Design space in black
% Support of the dirac measure: red points
% Size of the points represents their respective weights
% In case Christoffel polynomaials are used: 
% - one-dimensional: Christoffel polynomial in blue
% - two-dimensional: zero level set of Christoffel polynomial in blue
%                    mesh plot in figure 2