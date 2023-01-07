%% Singular momentum surface

g1 = sym('g1', 'real');
g2 = sym('g2', 'real');

beta = 54.73*pi/180;
h0 = 0.2288;
g3 = linspace(-pi,pi,10);
g4 = linspace(-pi,pi,10);
s = [];
for i3 = 1:length(g3)
    for i4 = 1:length(g4)
        D1 = sin(beta)*(sin(g2)*sin(g3(i3))*cos(g4(i4)) + cos(g2)*sin(g3(i3))*sin(g4(i4)) + ...
            cos(beta)*(cos(g2)*cos(g3(i3))*sin(g4(i4)) - sin(g2)*cos(g3(i3))*cos(g4(i4))) + ...
            2*cos(beta)^2*cos(g2)*cos(g3(i3))*cos(g4(i4)));
        
        D2 = sin(beta)*(sin(g3(i3))*sin(g4(i4))*cos(g1) + cos(g3(i3))*sin(g4(i4))*sin(g1) + ...
            cos(beta)*(cos(g3(i3))*cos(g4(i4))*sin(g1) - sin(g3(i3))*cos(g4(i4))*cos(g1)) + ...
            2*cos(beta)^2*cos(g3(i3))*cos(g4(i4))*cos(g1));
        eq1 = D1 == 0;
        eq2 = D2 == 0;
        [solg1,solg2] = solve(eq1,eq2,g1,g2);
        s1 = [double(solg1(1)) double(solg1(2)) double(solg1(3)) double(solg1(4))];
        s2 = [double(solg2(1)) double(solg2(2)) double(solg1(3)) double(solg1(4))];
        s = [s [s1; s2]];    
    end
end


%% Plot - singular momentum surface

h_t = [];
g32 = [];
for i = 1:10
    for j = 1:10
        g32 = [g32 [g3(i); g4(j)]];
    end
end
for i = 1:100
    g = [s(1,4*i-3),s(2,4*i-3),g32(1,i),g32(2,i)];
    h = compute_h(h0, beta, g);
    h_t = [h_t h];
end
for i = 1:100
    g = [s(1,4*i-2),s(2,4*i-2),g32(1,i),g32(2,i)];
    h = compute_h(h0, beta, g);
    h_t = [h_t h];
end
C = sqrt(h_t(1,:).^2+h_t(2,:).^2+h_t(3,:).^2);
figure
scatter3(h_t(1,:),h_t(2,:),h_t(3,:), 10, C, 'fill')
c = colorbar;
colormap default;
c.Label.String = 'Magnitude |h|';
xlabel('x direction [Nms]')
ylabel('y direction [Nms]')
zlabel('z direction [Nms]')
title('Angular momentum singularities')
grid on

%% Gamma

g1 = sym('g1', 'real');
g2 = sym('g2', 'real');


beta = 54.73*pi/180;   % 54.73
% h0 = 4;
n = 100;
g3 = linspace(-pi,pi,n);
% g4 = pi/3;
g4 = 100*pi/180;
P = SetupParameters;
s = [];
for i3 = 1:length(g3)
    D1 = sin(beta)*(sin(g2)*sin(g3(i3))*cos(g4) + cos(g2)*sin(g3(i3))*sin(g4) + ...
        cos(beta)*(cos(g2)*cos(g3(i3))*sin(g4) - sin(g2)*cos(g3(i3))*cos(g4)) + ...
        2*cos(beta)^2*cos(g2)*cos(g3(i3))*cos(g4));
    
    D2 = sin(beta)*(sin(g3(i3))*sin(g4)*cos(g1) + cos(g3(i3))*sin(g4)*sin(g1) + ...
        cos(beta)*(cos(g3(i3))*cos(g4)*sin(g1) - sin(g3(i3))*cos(g4)*cos(g1)) + ...
        2*cos(beta)^2*cos(g3(i3))*cos(g4)*cos(g1));
    eq1 = D1 == 0;
    eq2 = D2 == 0;
    [solg1,solg2] = solve(eq1,eq2,g1,g2);
    s1 = [double(solg1(1)) double(solg1(2)) double(solg1(3)) double(solg1(4))];
    s2 = [double(solg2(1)) double(solg2(2)) double(solg1(3)) double(solg1(4))];
    s = real([s [s1; s2]]);
%     sol = [double(solg1), double(solg2)];
%     m = Manipulability([sol(1), sol(2), g3, g4],P);
%     disp(sol)
%     disp(m),
    disp(i3)
end

%% Plot - Gamma plot

figure
g1 = s(1,1:4:end);
g2 = s(2,1:4:end);
scatter3(g1,g2,g3, 10,'fill','b')
hold on
g1 = s(1,2:4:end);
g2 = s(2,2:4:end);
scatter3(g1,g2,g3, 10,'fill','b')
hold on
g1 = s(1,3:4:end);
g2 = s(2,3:4:end);
scatter3(g1,g2,g3, 10,'fill','b')
hold on
g1 = s(1,4:4:end);
g2 = s(2,4:4:end);
scatter3(g1,g2,g3, 10,'fill','b')
hold on
plotcube([2*pi 2*pi 2*pi],[-pi -pi -pi],0,[0 0 1]);
xlabel('$\gamma_1$ [rad]','Interpreter','latex','FontSize',15);
ylabel('$\gamma_2$ [rad]','Interpreter','latex','FontSize',15);
zlabel('$\gamma_3$ [rad]','Interpreter','latex','FontSize',15);
title('Singularity envelope ($\gamma_4 = \pi/2$)','Interpreter','latex','FontSize',15);
subtitle('Pyramid array','FontSize',12);
xlim([-pi pi])
ylim([-pi pi])
zlim([-pi pi])
grid off
axis equal
box on
% view(335,20)
% view(70,30)
% view(0,0)
% view(90,0)
% view(0,90)

%%
% hold on
% g1 = s10(1,1:4:end);
% g2 = s10(2,1:4:end);
% scatter3(g1,g2,g3, 10,'fill','r')
% hold on
% g1 = s10(1,2:4:end);
% g2 = s10(2,2:4:end);
% scatter3(g1,g2,g3, 10,'fill','r')
% hold on
% g1 = s10(1,3:4:end);
% g2 = s10(2,3:4:end);
% scatter3(g1,g2,g3, 10,'fill','r')
% hold on
% g1 = s10(1,4:4:end);
% g2 = s10(2,4:4:end);
% scatter3(g1,g2,g3, 10,'fill','r')

% hold on
% g1 = s_5(1,1:4:end);
% g2 = s_5(2,1:4:end);
% scatter3(g1,g2,g3, 10,'fill','g')
% hold on
% g1 = s_5(1,2:4:end);
% g2 = s_5(2,2:4:end);
% scatter3(g1,g2,g3, 10,'fill','g')
% hold on
% g1 = s_5(1,3:4:end);
% g2 = s_5(2,3:4:end);
% scatter3(g1,g2,g3, 10,'fill','g')
% hold on
% g1 = s_5(1,4:4:end);
% g2 = s_5(2,4:4:end);
% scatter3(g1,g2,g3, 10,'fill','g')

% hold on
% g1 = s_10(1,1:4:end);
% g2 = s_10(2,1:4:end);
% scatter3(g1,g2,g3, 10,'fill','k')
% hold on
% g1 = s_10(1,2:4:end);
% g2 = s_10(2,2:4:end);
% scatter3(g1,g2,g3, 10,'fill','k')
% hold on
% g1 = s_10(1,3:4:end);
% g2 = s_10(2,3:4:end);
% scatter3(g1,g2,g3, 10,'fill','k')
% hold on
% g1 = s_10(1,4:4:end);
% g2 = s_10(2,4:4:end);
% scatter3(g1,g2,g3, 10,'fill','k')

% hold on
% g1 = s5(1,1:4:end);
% g2 = s5(2,1:4:end);
% scatter3(g1,g2,g3, 10,'fill','m')
% hold on
% g1 = s5(1,2:4:end);
% g2 = s5(2,2:4:end);
% scatter3(g1,g2,g3, 10,'fill','m')
% hold on
% g1 = s5(1,3:4:end);
% g2 = s5(2,3:4:end);
% scatter3(g1,g2,g3, 10,'fill','m')
% hold on
% g1 = s5(1,4:4:end);
% g2 = s5(2,4:4:end);
% scatter3(g1,g2,g3, 10,'fill','m')
% 
% hold on
% g1 = s0(1,1:4:end);
% g2 = s0(2,1:4:end);
% scatter3(g1,g2,g3, 10,'fill','y')
% hold on
% g1 = s0(1,2:4:end);
% g2 = s0(2,2:4:end);
% scatter3(g1,g2,g3, 10,'fill','y')
% hold on
% g1 = s0(1,3:4:end);
% g2 = s0(2,3:4:end);
% scatter3(g1,g2,g3, 10,'fill','y')
% hold on
% g1 = s0(1,4:4:end);
% g2 = s0(2,4:4:end);
% scatter3(g1,g2,g3, 10,'fill','y')

% hold on
% g1 = s180(1,1:4:end);
% g2 = s180(2,1:4:end);
% scatter3(g1,g2,g3, 10,'fill','g')
% hold on
% g1 = s180(1,2:4:end);
% g2 = s180(2,2:4:end);
% scatter3(g1,g2,g3, 10,'fill','g')
% hold on
% g1 = s180(1,3:4:end);
% g2 = s180(2,3:4:end);
% scatter3(g1,g2,g3, 10,'fill','g')
% hold on
% g1 = s180(1,4:4:end);
% g2 = s180(2,4:4:end);
% scatter3(g1,g2,g3, 10,'fill','g')
% 
% hold on
% g1 = s170(1,1:4:end);
% g2 = s170(2,1:4:end);
% scatter3(g1,g2,g3, 10,'fill','y')
% hold on
% g1 = s170(1,2:4:end);
% g2 = s170(2,2:4:end);
% scatter3(g1,g2,g3, 10,'fill','y')
% hold on
% g1 = s170(1,3:4:end);
% g2 = s170(2,3:4:end);
% scatter3(g1,g2,g3, 10,'fill','y')
% hold on
% g1 = s170(1,4:4:end);
% g2 = s170(2,4:4:end);
% scatter3(g1,g2,g3, 10,'fill','y')

hold on
g1 = s90(1,1:4:end);
g2 = s90(2,1:4:end);
scatter3(g1,g2,g3, 10,'fill','m')
hold on
g1 = s90(1,2:4:end);
g2 = s90(2,2:4:end);
scatter3(g1,g2,g3, 10,'fill','m')
hold on
g1 = s90(1,3:4:end);
g2 = s90(2,3:4:end);
scatter3(g1,g2,g3, 10,'fill','m')
hold on
g1 = s90(1,4:4:end);
g2 = s90(2,4:4:end);
scatter3(g1,g2,g3, 10,'fill','m')

hold on
g1 = s120(1,1:4:end);
g2 = s120(2,1:4:end);
scatter3(g1,g2,g3, 10,'fill','k')
hold on
g1 = s120(1,2:4:end);
g2 = s120(2,2:4:end);
scatter3(g1,g2,g3, 10,'fill','k')
hold on
g1 = s120(1,3:4:end);
g2 = s120(2,3:4:end);
scatter3(g1,g2,g3, 10,'fill','k')
hold on
g1 = s120(1,4:4:end);
g2 = s120(2,4:4:end);
scatter3(g1,g2,g3, 10,'fill','k')

hold on
g1 = s110(1,1:4:end);
g2 = s110(2,1:4:end);
scatter3(g1,g2,g3, 10,'fill','r')
hold on
g1 = s110(1,2:4:end);
g2 = s110(2,2:4:end);
scatter3(g1,g2,g3, 10,'fill','r')
hold on
g1 = s110(1,3:4:end);
g2 = s110(2,3:4:end);
scatter3(g1,g2,g3, 10,'fill','r')
hold on
g1 = s110(1,4:4:end);
g2 = s110(2,4:4:end);
scatter3(g1,g2,g3, 10,'fill','r')

%%
g = sym('g',[4 1]);
assume(g,'positive')
assumeAlso(g <= pi)
h0 = 1;
beta = 54.73*pi/180;
% beta = sym('s');
% Compute Jacobian
J = h0*[-cos(beta)*cos(g(1)) sin(g(2)) cos(beta)*cos(g(3)) -sin(g(4));
        -sin(g(1)) -cos(beta)*cos(g(2)) sin(g(3)) cos(beta)*cos(g(4));
         sin(beta)*cos(g(1)) sin(beta)*cos(g(2)) sin(beta)*cos(g(3)) sin(beta)*cos(g(4))];

simplify(det(J*J'))
% eqn = det(J*J') == 0;
% sol = solve(eqn,g);

%% Functions
function momentum = compute_h(h0,beta,g)
    momentum = h0*[-cos(beta)*sin(g(1))-cos(g(2))+cos(beta)*sin(g(3))+cos(g(4));
                    cos(g(1))-cos(beta)*sin(g(2))-cos(g(3))+cos(beta)*sin(g(4));
                    sin(beta)*sin(g(1))+sin(beta)*sin(g(2))+sin(beta)*sin(g(3))+sin(beta)*sin(g(4))];
end




