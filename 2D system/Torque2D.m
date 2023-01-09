%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-CMG torque visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hugo Pereira
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("..\funcs\")

%% Parameters
h0 = 1;                  % Angular momentum of each CMG
lambda = 30*pi/180;      % Triangle inner angle

%% Torque envelope - 2D
g = [0;2*pi/3;-2*pi/3];     % Non-singular configuration
% g = [0;0;0];                % Singularity
J = TriangleJacobian(g,h0,lambda);

samples = 20;
gr1 = linspace(-1.5,1.5,samples);
gr2 = linspace(-1.5,1.5,samples);
gr3 = linspace(-1.5,1.5,samples);

hDot = zeros(2,samples^4);
j = 1;
for i1 = 1:length(gr1)
    for i2 = 1:length(gr2)
        for i3 = 1:length(gr3)
            gr = [gr1(i1);gr2(i2);gr3(i3)];
            hDot(:,j) = J*gr;
            j = j + 1;
        end
    end
end

C = sqrt(hDot(1,:).^2+hDot(2,:).^2);
scatter(hDot(1,:),hDot(2,:),10,C,'filled')
xlabel('$$\dot{h}_x/h_0$$ [$$s^{-1}$$]','Interpreter','latex','FontSize',15)
ylabel('$$\dot{h}_y/h_0$$ [$$s^{-1}$$]','Interpreter','latex','FontSize',15)
title('Torque envelope','FontSize',15)
subtitle('Triangular array','FontSize',12)
grid off
box off
axis square
xlim([-3.5 3.5])
ylim([-3.5 3.5])
c = colorbar;
c.Label.String = 'Magnitude';
c.Label.FontSize = 12;