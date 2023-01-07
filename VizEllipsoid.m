%% Ellipsoid - plots

h0 = 1;
% g = [pi/2;0;0;-pi/2];
g = [pi/2;0;0;0];
% g = [0;0;0;0];
J = computeJac(g,h0);

% SVD
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
subtitle('Roof array','FontSize',12)
grid off
box off
axis square
xlim([-2 2])
ylim([-2 2])
zlim([-2 2])
c = colorbar;
c.Label.String = 'Magnitude';
c.Label.FontSize = 12;

% view(0, 90)

%% Ellipsoid - movie

% Initialize video
myVideo = VideoWriter('ellipsoid_movie','Uncompressed AVI'); %open video file
myVideo.FrameRate = 20;  
open(myVideo)

xc = 0;
yc = 0;
zc = 0;

for i = 1:2:length(out.gimbal.Data)
    g = out.gimbal.Data(i,:)';
    J = [-cos(beta)*cos(g(1)) sin(g(2)) cos(beta)*cos(g(3)) -sin(g(4));
        -sin(g(1)) -cos(beta)*cos(g(2)) sin(g(3)) cos(beta)*cos(g(4));
        sin(beta)*cos(g(1)) sin(beta)*cos(g(2)) sin(beta)*cos(g(3)) sin(beta)*cos(g(4))];
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
    colorbar 
    caxis([0 2])
    colormap default
    axis([-2 2 -2 2 -2 2])
    xlabel('x direction')
    ylabel('y direction')
    zlabel('z direction')
    title('Principal Components Visualization')
    grid on

    pause(0.025)
    
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end

close(myVideo)

%% Determinant - movie

% Initialize video
myVideo = VideoWriter('determinant_movie','Uncompressed AVI'); %open video file
myVideo.FrameRate = 20; 
open(myVideo)

for i = 1:2:length(out.det.Data)
    plot(out.det.Time(1:i),out.det.Data(1:i),'r','LineWidth',1)
    axis([0 45 0 2])
    xlabel('Time [s]')
    ylabel('Jacobian Determinant')
    title('Singularity measure')
    grid on
     
    pause(0.025)
    
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end

close(myVideo)

%% Jacobian 

function J = computeJac(g,h0)
J = h0*[ cos(g(1))   -sin(g(2))             0            0;
                 0            0     sin(g(3))    cos(g(4));
         sin(g(1))    cos(g(2))     cos(g(3))   -sin(g(4))];
end