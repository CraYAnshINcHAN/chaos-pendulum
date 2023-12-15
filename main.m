L1 = 0.15; L2 = 0.1;
m1 = 0.6; m2 = 0.5;
J1 = 1/12 * m1 * L1^2; J2 = 1/12 * m2 * L2^2;
c = 0.05;
g = 9.8;

%%
syms l1 l2 c1 c2 'positive'
syms theta1 theta2 omega1 omega2 'real'
C = [0 0.5];
CS1 = 1/2 * L1 - l1;
S1 = C  + CS1 * [cos(theta1) sin(theta1)];
D = C - l1 * [cos(theta1) sin(theta1)];
DS2 = 1/2 * L2 - l2;
S2 = D + DS2 * [cos(theta2) sin(theta2)];
V1 = diff(S1,theta1,1) * omega1 + diff(S1,theta2,1) * omega2;
V2 = diff(S2,theta1,1)*omega1 + diff(S2,theta2,1)*omega2;
EK = 1/2 * (V1 * m1 * V1' + V2 * m2 * V2' + ...
    omega1 * J1 * omega1 + omega2 * J2 * omega2);
EP = m1 * S1(2)*g + m2 * S2(2)*g;

EK_1 = diff(EK,theta1,1);
EK_2 = diff(EK,theta2,1);
EK_11 = diff(EK,omega1,1);
EK_22 = diff(EK,omega2,1);
EK_1_11 = diff(EK_1,omega1,1);
EK_1_22 = diff(EK_1,omega2,1);
EK_2_11 = diff(EK_2,omega1,1);
EK_2_22 = diff(EK_2,omega2,1);
EK_11_11 = diff(EK_11,omega1,1);
EK_11_22 = diff(EK_11,omega2,1);
EK_22_22 = diff(EK_22,omega2,1);



EP_1 = diff(EP,theta1,1);
EP_2 = diff(EP,theta2,1);


S1y = S1(2);
S1y_1 = diff(S1y,theta1,1);
S1y_2 = diff(S1y,theta2,1);
S2y = S2(2);
S2y_1 = diff(S2y,theta1,1);
S2y_2 = diff(S2y,theta2,1);


F1 = -c1*omega1 - c2*(omega1 - omega2) - m1*g*S1y_1 - m2*g*S2y_1;
F2 = -c2*(omega2 - omega1) - m1*g*S1y_2 - m2*g*S2y_2;

%%
%赋值
A = [EK_11_11,EK_11_22;EK_11_22,EK_22_22];
B = [F1 + EK_1 -EP_1 -EK_1_11*omega1 - EK_2_11*omega2;F2 + EK_2 -EP_2 -EK_1_22*omega1 - EK_2_22*omega2];
alpha = A\B;
alpha1 = simplify(alpha(1));
alpha2 = simplify(alpha(2));


E_final = C(2)*(m1 + m2)*g - abs(((0.5*L1 - l1)*m1-m2*l1))*g - abs(m2*(0.5*L2 - l2))*g;
%%
%构造matlab函数
matlabFunction(eval([omega1;omega2;alpha1;alpha2]),'Vars', {[theta1 theta2 omega1 omega2]', [l1 l2 c1 c2]}, 'File', 'hodefcn');
matlabFunction(EK + EP - 1.001 * E_final, 'Vars', {[theta1 theta2 omega1 omega2]', [l1 l2 c1 c2]}, 'File', 'engyfcn');
rehash

%%
%题1

% l1 = 0.05;
% l2 = 0.025;
% [t, y, te] = cal_ode([l1 l2 0.05 0.05]);

%%
%画视频
% G = D - l2*[cos(theta2) sin(theta2)];
% F = D + (L2-l2)*[cos(theta2) sin(theta2)];
% E = C + (L1 - l1)*[cos(theta1) sin(theta1)];
% VE = diff(E,theta1,1) * omega1 + diff(E,theta2,1) * omega2;
% VF = diff(F,theta1,1) * omega1 + diff(F,theta2,1) * omega2;
% figure(1)
% axis equal
% grid on
% hold on
% axis([-0.2 0.2, 0.25 0.65])
% set(gca, 'xtick', -0.2:0.05:0.2, 'ytick', 0.25:0.05:65)
% plot([-0.05;0.05], [0.25;0.25], 'k-')
% plot([-0.05;C(1)], [0.25;C(2)], 'k-')
% plot([0.05;C(1)], [0.25;C(2)], 'k-')
% plot(C(1), C(2), 'ok', 'MarkerSize', 4);
% Dr = eval(subs(D,[theta1 theta2 omega1 omega2],y(1,:)));
% Er = eval(subs(E,[theta1 theta2 omega1 omega2],y(1,:)));
% Fr = eval(subs(F,[theta1 theta2 omega1 omega2],y(1,:)));
% S1r = eval(subs(S1,[theta1 theta2 omega1 omega2],y(1,:)));
% S2r = eval(subs(S2,[theta1 theta2 omega1 omega2],y(1,:)));
% Gr = eval(subs(G,[theta1 theta2 omega1 omega2],y(1,:)));
% 
% DEr = [Dr;Er];
% GFr = [Gr;Fr];
% DEplot = plot(DEr(:,1),DEr(:,2), '-', 'Color', '#61E161', 'LineWidth', 2);
% GFplot = plot(GFr(:,1),GFr(:,2), '-', 'Color', '#61E161', 'LineWidth', 2);
% Dplot = plot(Dr(:,1),Dr(:,2), 'ok', 'MarkerSize', 4);
% S1plot = plot(S1r(:,1),S1r(:,2), '.r', 'MarkerSize', 8);
% S2plot = plot(S2r(:,1),S2r(:,2), '.r', 'MarkerSize', 8);
% video = VideoWriter(['out/', mfilename]);
% video.FrameRate = 60;
% open(video);
% FF = griddedInterpolant(t, y);
% tq = 0:(1 / 60):te;
% for i = 1:length(tq)
%     yq = FF(tq(i));
%     Dr = eval(subs(D, [theta1, theta2, omega1, omega2], yq));
%     Er = eval(subs(E, [theta1, theta2, omega1, omega2], yq));
%     Fr = eval(subs(F, [theta1, theta2, omega1, omega2], yq));
%     Gr = eval(subs(G, [theta1, theta2, omega1, omega2], yq));
%     S1r = eval(subs(S1, [theta1, theta2, omega1, omega2], yq));
%     S2r = eval(subs(S2, [theta1, theta2, omega1, omega2], yq));
%     DEr = [Dr; Er];
%     GFr = [Gr; Fr];
%     set(DEplot, 'XData', DEr(:,1), 'YData', DEr(:,2));
%     set(GFplot, 'XData', GFr(:,1), 'YData', GFr(:,2));
%     set(Dplot, 'XData', Dr(:, 1), 'YData', Dr(:, 2));
%     set(S1plot, 'XData', S1r(:, 1), 'YData', S1r(:, 2));
%     set(S2plot, 'XData', S2r(:, 1), 'YData', S2r(:, 2));
%     drawnow
%     curframe = getframe;
%     writeVideo(video, curframe);
% end
% close(video);
%%
%画图


% n = length(t);
% Erx = zeros(n,1);
% Ery = zeros(n,1);
% VErx = zeros(n,1);
% VEry = zeros(n,1);
% for i =1:n
%     Er = eval(subs(E,[theta1 theta2 omega1 omega2],y(i,:)));
%     Ery(i) = Er(2);
%     Erx(i) = Er(1);
%     VEr = eval(subs(VE,[theta1 theta2 omega1 omega2],y(i,:)));
%     VErx(i) = VEr(:,1);
%     VEry(i) = VEr(:,2);
% end
% t_ = t(1 : n-1);
% aErx = diff(VErx)./diff(t);
% aEry = diff(VEry)./diff(t);
% 
% subplot(1,3,1);
% plot3(Erx,Ery,t,"LineWidth",2,"Color","red");
% xlabel('xE(m)');
% ylabel('yE(m)');
% zlabel('t(s)');
% grid on
% 
% subplot(1,3,2);
% plot3(VErx,VEry,t,"LineWidth",2,"Color","red");
% xlabel('vxE(m/s)');
% ylabel('vyE(m/s)');
% zlabel('t(s)');
% grid on
% 
% subplot(1,3,3);
% plot3(aErx,aEry,t_,"LineWidth",2,"Color","red");
% xlabel('axE(m/s^2)');
% ylabel('ayE(m/s^2)');
% zlabel('t(s)');
% grid on

% Frx = zeros(n,1);
% Fry = zeros(n,1);
% VFrx = zeros(n,1);
% VFry = zeros(n,1);
% for i =1:n
%     Fr = eval(subs(F,[theta1 theta2 omega1 omega2],y(i,:)));
%     Fry(i) = Fr(2);
%     Frx(i) = Fr(1);
%     VFr = eval(subs(VF,[theta1 theta2 omega1 omega2],y(i,:)));
%     VFrx(i) = VFr(:,1);
%     VFry(i) = VFr(:,2);
% end
% aFrx = diff(VFrx)./diff(t);
% aEry = diff(VFry)./diff(t);
% 
% subplot(1,3,1);
% plot3(Frx,Fry,t,"LineWidth",2,"Color","red");
% xlabel('xF(m)');
% ylabel('yF(m)');
% zlabel('t(s)');
% grid on
% 
% subplot(1,3,2);
% plot3(VErx,VEry,t,"LineWidth",2,"Color","red");
% xlabel('vxF(m/s)');
% ylabel('vyF(m/s)');
% zlabel('t(s)');
% grid on
% 
% subplot(1,3,3);
% plot3(aErx,aEry,t_,"LineWidth",2,"Color","red");
% xlabel('axF(m/s^2)');
% ylabel('ayF(m/s^2)');
% zlabel('t(s)');
% grid on
%%
%题2-1
% nx = 101;ny = 201;
% x1 = linspace(0,L1,nx);
% y1 = linspace(0,L2,ny);
% T = zeros(nx,ny);
% for i =1: nx
%     for j = 1:ny
%         [t11,y11,te1] = cal_ode([x1(i) y1(j) 0.05 0.05]);
%         T(i,j) = te1;
%     end
% end
% mesh(y1,x1,T);

%%
%题2-2
% ny = 10000;
% y1 = linspace(0.045,0.052,ny);
% TT = zeros(ny,1);
% for i =1:ny
%     [t11,y11,te1] = cal_ode([0.15 y1(i)  0.05 0.05]);
%     TT(i)=te1;
% % end
% plot(y1,TT);
%l1 = 0.15 ; l2 = 0.0495309;


%改阻尼
nc1 = 41;nc2 = 61;
xc1 = linspace(0.005,1,nc1);
yc2 = linspace(0.005,1,nc2);
T = zeros(nc1,nc2);
for i =1: nc1
    for j = 1:nc2
        [t11,y11,te1] = cal_ode([0.15 0.0495309 xc1(i) yc2(j)]);
        T(i,j) = te1;
    end
end
mesh(yc2,xc1,T);

%%
%
function [t ,y, te] = cal_ode(params)
    options = odeset("Events", @odestopfunction);
    if engyfcn([0;0;0;0], params) < 0
        t = 0; y = [0 0 0 0]; te = NaN;
    else
        [t , y, te, ~, ~] = ode45(@(~, y) hodefcn(y, params), [0 10000], [0 0 0 0], options);
    end

    function [value, isterminal, direction] = odestopfunction(~, y)
        value = engyfcn(y, params);
        isterminal = true;
        direction = 0;
    end
end
