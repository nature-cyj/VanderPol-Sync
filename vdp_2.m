% James Panteleone
% Nondemensionalized equation
% for solving
% y1''+(1+d)sin(y1)+u((y1/y0)^2 - 1)y1'-b*cos(y1)*(sin(y1)+sin(y2))'' == 0 
% y2''+(1-d)sin(y1)+u((y2/y0)^2 - 1)y1'-b*cos(y2)*(sin(y1)+sin(y2))'' == 0 

g = 9.81;
theta0 = 0.45; Van_e = 0.425;      % theta_0는 단위 rad
m_pen = 16.3e-3; m_met = 87.66e-3; % 위에 다는 추 제외
m_plt = 37.8e-3; % 2.7mm 작은 아크릴 판 기준
r_cm = 22.3e-3; % 질량중심까지의 거리

w = 208*pi/60; % 기본 각진동수
I = m_pen*g*r_cm / w^2;

% dimensionless parameters
beta = (m_pen * r_cm)^2 / (2*m_met + m_plt + 2* m_pen) / I;
delta = 0;   % (w1 - w2)/w
mu = 0.01;    
y_0 = 0.39;

syms y1(tau) y2(tau)
eqn1 = diff(y1,2)+(1+delta)*sin(y1)+mu*((y1/y_0)^2-1)*diff(y1)-beta*cos(y1)*diff(sin(y1)+sin(y2),2) == 0;
eqn2 = diff(y2,2)+(1-delta)*sin(y2)+mu*((y2/y_0)^2-1)*diff(y2)-beta*cos(y2)*diff(sin(y1)+sin(y2),2) == 0;
eqns = [eqn1 eqn2];

[V,S] = odeToVectorField(eqns)
K = matlabFunction(V,'vars',{'tau','Y'});

interval = [0 1500]; % 푸는 범위
% 초기조건
% [y1(0) y1'(0) y2(0) y2'(0)]
conds = [1 0 -1.0001 0];

[t,y] = ode45(K,interval,conds);
y_1 = y(:,1);
y_2 = y(:,3);

%% Plotting the angle
figure;
plot(t/w,y_1,t/w,y_2);
title('Angle of 2 Metronomes','FontSize',18);
xlabel('time (s)','FontSize',13);
ylabel('Angle (rad)','FontSize',13);
legend('\theta_1','\theta_2')


%% hilbert transformation of result
y_1_h = hilbert(y_1); y_2_h = hilbert(y_2);
phase_1 = unwrap(angle(y_1_h)); phase_2 = unwrap(angle(y_2_h));

% figure;subplot(2,1,1);plot(t,phase_1);subplot(2,1,2);plot(t,phase_2);
phase_diff = phase_1 - phase_2;

% figure; subplot(2,1,1);plot(t/w,y_1_r,t/w,y_1_i);subplot(2,1,2);plot(t/w,y_2_r,t/w,y_2_i);
figure; plot(t/w, phase_diff); title('Phase difference of two metronomes','FontSize',17);
xlabel('time (s)','FontSize',13); ylabel('Phase difference (rad)','FontSize',13);

clear eqn eqn1 eqn2 y1 y2 cond1 cond2