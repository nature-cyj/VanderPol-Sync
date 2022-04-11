% Remove spikes of hilbert
% Multiple plots
phase_list = [];
tic
for i = 0:0.005:0.015
w = 208*pi/60;
beta = 0.018;
delta = 0;   % (w1 - w2)/w
mu = 0.002+i;    
y_0 = 0.39;

syms y1(tau) y2(tau)
eqn1 = diff(y1,2)+(1+delta)*sin(y1)+mu*((y1/y_0)^2-1)*diff(y1)-beta*cos(y1)*diff(sin(y1)+sin(y2),2) == 0;
eqn2 = diff(y2,2)+(1-delta)*sin(y2)+mu*((y2/y_0)^2-1)*diff(y2)-beta*cos(y2)*diff(sin(y1)+sin(y2),2) == 0;
eqns = [eqn1 eqn2]; [V,S] = odeToVectorField(eqns);
K = matlabFunction(V,'vars',{'tau','Y'});

interval = [0 700]; % 푸는 범위
% 초기조건 [y1(0) y1'(0) y2(0) y2'(0)]
conds = [1 0 -1.04 0];

[t,y] = ode45(K,interval,conds);
y_1 = y(:,1); y_2 = y(:,3);
y_1_h = hilbert(y_1); y_2_h = hilbert(y_2);
phase_1 = unwrap(angle(y_1_h)); phase_2 = unwrap(angle(y_2_h));

% figure;subplot(2,1,1);plot(t,phase_1);subplot(2,1,2);plot(t,phase_2);
phase_diff = phase_1 - phase_2;
phase_list = [phase_list phase_diff(1:3300)];
%phase_m = medfilt1(phase_diff,100); % 튀는 값 날리기
% phase_list = [phase_list phase_m(1:3300)];
% medfilt를 쓰면 인덱스 수가 달라져서 하나로 만들려고

clear eqn eqn1 eqn2 y1 y2 cond1 cond2
end
toc

% 3d plotting

surf(phase_list,'Edgecolor','None')