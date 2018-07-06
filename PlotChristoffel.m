load wynn_christoffel3
close;
mesh(X1,X2,-Ch);
axis([-0.5 0.9 -0.5 0.9 -max(max(Ch)) 0])
hold on;
    % plot polygon
    plot([-1,1]/2/sqrt(2),[-1,-1]/2/sqrt(2),'k','linewidth',3);
    plot([-1,-1]/2/sqrt(2),[-1,1]/2/sqrt(2),'k','linewidth',3);
    plot([-1,2]/2/sqrt(2),[1,2]/2/sqrt(2),'k','linewidth',3);
    plot([1,2]/2/sqrt(2),[-1,2]/2/sqrt(2),'k','linewidth',3);
contour(X1,X2,Ch,linspace(-2e-2,2e-2,10),'b','linewidth',2)
    % Plot support points
    scatter(pts(1,:),pts(2,:),'o', 'MarkerFaceColor','r', 'MArkerEdgeColor', 'none', 'SizeData',w.*1500);
xlabel('x_1');ylabel('x_2');
view(-30,45);
