%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-CMG momentum visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hugo Pereira
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("..\funcs\")

%% Parameters
h0 = 1;                  % Angular momentum of each CMG
lambda = 30*pi/180;      % Triangle inner angle

%% Angular momentum - 2D

samples = 30;
g1 = linspace(-pi,pi,samples);
g2 = linspace(-pi,pi,samples);
g3 = linspace(-pi,pi,samples);

j = 1;
h = zeros(2,samples^3);
for i1 = 1:length(g1)
    for i2 = 1:length(g2)
        for i3 = 1:length(g3)
            g = [g1(i1);g2(i2);g3(i3)];
            h(:,j) = TriangleMomentum(g,h0,lambda);
            j = j + 1;
        end
    end
end
C = sqrt(h(1,:).^2+h(2,:).^2);
scatter(h(1,:),h(2,:),[],C,'filled')
xlabel('$$h_x/h_0$$','Interpreter','latex','FontSize',15)
ylabel('$$h_y/h_0$$','Interpreter','latex','FontSize',15)
title('Angular momentum envelope','FontSize',15)
subtitle('Triangular array','FontSize',12)
grid off
box off
axis square
xlim([-3.5 3.5])
ylim([-3.5 3.5])
c = colorbar;
c.Label.String = 'Magnitude';
c.Label.FontSize = 12;