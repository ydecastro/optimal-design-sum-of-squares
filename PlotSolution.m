% Roxana Hess, September 2017

% Plot solution

function PlotSolution(expl,d,q,pts,w,Ch,M)
axes('FontSize',14)
hold on

if expl == 1
    
    % plot interval
    plot([-1;1], [0;0], 'k', 'linewidth',3);
    % plot support of the dirac and corresponding weights
    if ~isempty(pts)
        plot(pts(1,:),0,'ro','MarkerFaceColor','r', 'MarkerSize',10);
        plot([pts(1,:);pts(1,:)],[zeros(length(w),1)';w], 'g', 'linewidth',2);
        axis([-1.2 1.2 -.1 max(w)+.1])
    end
    % plot polynomial p* = s(d) - Christoffel
    if ~isempty(Ch)
        X = linspace(-1,1,200);
        plot(X,Ch,'b','linewidth',2);
        axis([-1.2 1.2 -.1 max(Ch)+.1])
    end
    
%     % Christoffel-like polynomial for T-optimal design
%     % It is independent of the design space and can be computed without
%     % any optimization as done below.
%     if q == 1 
%         X = linspace(-2,2,100);
%         e = 2*[0:d]; Y = zeros(size(X));
%         for i = 1 : 100
%             Y(i) = sum((X(i)*ones(1,length(e))).^e);
%         end
%         plot(X,trace(M)-Y,'b','linewidth',2)
%         axis([-2 2 -1 1])
%     end

elseif expl == 2 || expl == 3 || expl == 4 || expl == 5
    
    if expl == 2
    % plot polygon
    plot([-1,1]/2/sqrt(2),[-1,-1]/2/sqrt(2),'k','linewidth',3);
    plot([-1,-1]/2/sqrt(2),[-1,1]/2/sqrt(2),'k','linewidth',3);
    plot([-1,2]/2/sqrt(2),[1,2]/2/sqrt(2),'k','linewidth',3);
    plot([1,2]/2/sqrt(2),[-1,2]/2/sqrt(2),'k','linewidth',3);
    
    elseif expl == 3
    % plot ellipses
    [X1,X2] = meshgrid(linspace(-1,1,100));
    contour(X1,X2, 9*X1.^2 + 13*X2.^2 - 7.3, [0 0], 'k','linewidth',3);
    contour(X1,X2, 5*X1.^2 + 13*X2.^2 - 2, [0 0], 'k','linewidth',3);
    
    elseif expl == 4
    % plot moon
    % plot section of bigger circle
    ang=0.5:0.01:2*pi-.5; xp=.6*cos(ang); yp=.6*sin(ang);
    plot(-.2+xp,yp,'k','linewidth',3);
    % plot section of smaller circle
    ang=pi-.8:0.01:pi+.8; xp=.4*cos(ang); yp=.4*sin(ang);
    plot(.6+xp,yp,'k','linewidth',3);
    
    elseif expl == 5
    % plot folium
    [X1,X2] = meshgrid(linspace(-1,1,100));
    contour(X1,X2, -X1.*(X1.^2-2*X2.^2)-(X1.^2+X2.^2).^2, [0 0], 'k','linewidth',3);
    end
    
    % plot support and weights
    if ~isempty(pts)
        % plot support of dirac with weights visualized by marker size
        scatter(pts(1,:),pts(2,:),'o', 'MarkerFaceColor','r', 'MArkerEdgeColor', 'none', 'SizeData',w.*1500);
%         % plot support of dirac with weights visualized by marker size and
%         % intensity of marker color
%         t = 0 : pi/10 : 2*pi;
%         for i = 1 : length(w)
%             point = patch((sin(t).*w(i)/2+ pts(1,i)),(cos(t).*w(i)/2+pts(2,i)),'r','edgecolor','none');
% %             point = patch((sin(t).*.1/2+ pts(1,i)),(cos(t).*.1/2+pts(2,i)),'r','edgecolor','none'); % same marker size
%             alpha(point,w(i)/max(w)); % Matlab2014a does not yet have the option MarkerFaceAlpha, so I do it with patch
%         end
%         % plot support of dirac with weights visualized by green bars
%         plot(pts(1,:),pts(2,:),'ro','MarkerFaceColor','r', 'MarkerSize',10);
%         plot3([pts(1,:);pts(1,:)],[pts(2,:);pts(2,:)],[zeros(length(w),1)';w], 'g', 'linewidth',2);
    end
    
    % plot polynomial p* = s(d) - Christoffel
    if ~isempty(Ch)
        [X1,X2] = meshgrid(linspace(-1,1,100));
        contour(X1,X2,Ch,[0 0],'b','linewidth',2);
    end
    
    xlabel('x_1','FontSize',14);
    ylabel('x_2','FontSize',14);
    if expl == 2
        axis([-.5 1 -.5 1])
    elseif expl == 3
        axis([-1 1 -1 1])
    elseif expl == 4
        axis([-.8 .8 -.8 .8])
    elseif expl == 5
        axis([-1.2 .8 -.8 .8])
    end
    axis equal
    box on

    % For a mesh plot of p* = s(d) - Christoffel polynomial
    if ~isempty(Ch)
        % set all negative values of the matrix representing p*to 0 in
        % order to see what happens for non-negative values
        for i = 1 : length(X1)
            for j = 1 : length(X1)
                if Ch(i,j) < 0
                    Ch(i,j) = 0;
                end
            end
        end
        
      figure(2)
      axes('FontSize',14)
      hold on
      
      mesh(X1,X2,Ch);
      White(:,:,1) = ones(size(Ch)); % spezifying color
      White(:,:,2) = ones(size(Ch));
      White(:,:,3) = ones(size(Ch));
      mesh(X1,X2,zeros(size(Ch)),White);
        
      xlabel('x_1','FontSize',14);
      ylabel('x_2','FontSize',14);
      if expl == 2
          axis([-.5 1 -.5 1])
      elseif expl == 3
          axis([-1 1 -1 1])
      elseif expl == 4
          axis([-.8 .8 -.8 .8])
      elseif expl == 5
          axis([-1.2 .8 -.8 .8])
      end
      view(50,50)
      grid on
    end
    
%     % Christoffel-like polynomial for T-optimal design
%     % It is independent of the design space and can be computed without
%     % any optimization as done below
%     if q == 1 && isempty(Ch)
%         [X1,X2] = meshgrid(linspace(-1,1,100));
%         g = 2*genpow(2+1,d); Y = zeros(size(X1)); e = ones(size(g,1),1);
%         for i = 1 : 100
%             for j = 1 : 100
%                 Y(i,j) = sum((X1(i,j)*e).^g(:,2).*(X2(i,j)*e).^g(:,3));
%             end
%         end
%         %mesh(X1,X2,trace(M)-Y);
%         contour(X1,X2,trace(M)-Y,[0 0],'b','linewidth',2);
%     end
    
elseif expl == 6
    
    % plot sphere
    [X1,X2,X3]=sphere(50);
    surf(X1,X2,X3);
    
    % plot support of dirac
    if ~isempty(pts)
        plot3(pts(1,:),pts(2,:),pts(3,:),'.r','markersize',30);
    end
    
    %axis vis3d
    axis equal
    view(115,40)
    %camlight(0,0)
    %camlight(45,45)
    xlabel('x_1','FontSize',14)
    ylabel('x_2','FontSize',14)
    zlabel('x_3','FontSize',14)
    colormap gray
    grid on
 
%     % Christoffel-like polynomial for T-optimal design
%     % It is independent of the design space and can be computed without
%     % any optimization as done below
%     if q == 1 && isempty(Ch)
%         [X1,X2,X3] = meshgrid(linspace(-1,1,100));
%         Y = zeros(size(X1));
%         g = 2*genpow(n+1,d); e = ones(size(g,1),1);
%         for i = 1 : 100
%             for j = 1 : 100
%                 for k = 1 : 100
%                     Y(i,j,k) = sum((X1(i,j,k)*e).^g(:,2).*(X2(i,j,k)*e).^g(:,3).*(X3(i,j,k)*e).^g(:,4));
%                 end
%             end
%         end
%     figure(2)
%     isosurface(X1,X2,X3,Y,2);
%     end
    
end

end