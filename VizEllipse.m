%% Ellipse - plots

h0 = 1;
alpha = 30*pi/180;
g = [0;2*pi/3;-2*pi/3];
% g = [0;0;0];
J = computeJac(g,h0,alpha);

% SVD4444
[S,Sigma,V] = svd(J);

xc = 0;
yc = 0;
zc = 0;

xr = Sigma(1,1);
yr = Sigma(2,2);
zr = Sigma(3,3);
[X,Y,Z] = ellipsoid(xc,yc,zc,xr,yr,zr);
Xc = reshape(X,[],1);
Yc = reshape(Y,[],1);
Zc = reshape(Z,[],1);
coord = [Xc, Yc, Zc];
Rot = S*coord';

Xr = Rot(1,:);
Yr = Rot(2,:);
Zr = Rot(3,:);

Xs = reshape(Xr, 21,21);
Ys = reshape(Yr, 21,21);
Zs = reshape(Zr, 21,21);
C = sqrt(Xs.^2+Ys.^2+Zs.^2);

mesh(Xs,Ys,Zs,C);

xlabel('$$x$$ axis','Interpreter','latex','FontSize',15)
ylabel('$$y$$ axis','Interpreter','latex','FontSize',15)
zlabel('$$z$$ axis','Interpreter','latex','FontSize',15)
title('Principal Components Visualization','FontSize',15)
subtitle('Triangular array','FontSize',12)
grid off
box off
axis square
xlim([-1.5 1.5])
ylim([-1.5 1.5])
c = colorbar;
c.Label.String = 'Magnitude';
c.Label.FontSize = 12;

view(0, 90)

%% Ellipsoid - movie

% Initialize video
myVideo = VideoWriter('ellipsoid_movie','Uncompressed AVI'); % open video file
myVideo.FrameRate = 20;  
open(myVideo)

xc = 0;
yc = 0;
zc = 0;

alpha = 30*pi/180;
h0 = 1;

t = 0:0.1:90;
for i = 1:2:length(g)
    J = computeJac(g(:,i),h0,alpha);
    % SVD
    [S,Sigma,V] = svd(J);
    
    xr = Sigma(1,1);
    yr = Sigma(2,2);
    zr = Sigma(3,3);
    [X,Y,Z] = ellipsoid(xc,yc,zc,xr,yr,zr);
    Xc = reshape(X,[],1);
    Yc = reshape(Y,[],1);
    Zc = reshape(Z,[],1);
    coord = [Xc, Yc, Zc];
    Rot = S*coord';

    Xr = Rot(1,:);
    Yr = Rot(2,:);
    Zr = Rot(3,:);

    Xs = reshape(Xr, 21,21);
    Ys = reshape(Yr, 21,21);
    Zs = reshape(Zr, 21,21);
    C = sqrt(Xs.^2+Ys.^2+Zs.^2);
    
    mesh(Xs,Ys,Zs,C);

    xlabel('$$x$$ axis','Interpreter','latex','FontSize',15)
    ylabel('$$y$$ axis','Interpreter','latex','FontSize',15)
    zlabel('$$z$$ axis','Interpreter','latex','FontSize',15)
    title('Principal Components Visualization','FontSize',15)
    subtitle('Triangular array','FontSize',12)
    grid off
    box off
    axis square
    xlim([-2 2])
    ylim([-2 2])
    c = colorbar;
    c.Label.String = 'Magnitude';
    c.Label.FontSize = 12;
    
    view(0, 90)

    pause(0.025)
    
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end

close(myVideo)

%% Determinant - movie

% Initialize video
myVideo = VideoWriter('determinant_movie','Uncompressed AVI'); % open video file
myVideo.FrameRate = 20; 
open(myVideo)

t = 0:0.1:90;
for i = 1:2:length(manipulability)
    plot(t(1:i),manipulability(1:i),'r','LineWidth',1)
    xlabel('time [s]','Interpreter','latex','FontSize',15);
    ylabel('$$\sqrt{det(JJ^T)}$$','Interpreter','latex','FontSize',15);
    title('Manipulability index','Interpreter','latex','FontSize',15);
    ylim([0 30])
    xlim([0 90])
    % grid on
    box off
     
    pause(0.025)
    
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)

%% Jacobian 

function J = computeJac(g,h0,alpha)
J = h0*[-sin(alpha + g(1))   cos(g(2))  -sin(alpha - g(3))  0;
         cos(alpha + g(1))   sin(g(2))  -cos(alpha - g(3))  0;
                        0           0                  0  0];
end
