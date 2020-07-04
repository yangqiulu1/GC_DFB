clear all;

%Na = 3.22; Nb = 2.888839004576706; %由κL=3.32558+0j算得一个组合 【例1】
%Na = 3.22; Nb = 2.596922074398683; %由κL=6.451633+0j算得一个组合 【例2】
%Na = 3.22; Nb = 2.891556951572583 - 0.031044851809367i;%由κL=3.225859+0.322079j算得一个组合 【例三】
Na = 3.22; Nb = 3.098876636825481 - 0.310924771189149i;%由κL=1+3i算得一个组合【另一文献1989Nakano】

r = (Na-Nb)/(Na+Nb);
ta_COE = 2*Na/(Na+Nb); % t本该包含 场的重叠比例【式10a】，这里用t_COE表示t/重叠比例【比例是个z向无关的常数，不影响后面阈值的位置，仅影响发射谱幅值】
tb_COE = 2*Nb/(Na+Nb);
kappaLAM = log(Na/Nb); % 文章中这个形式比较精确
% 解锁以下程序，获得ln(Na/Nb)与2r 近似差异的对比
% kappaLAM1 = 2*r
% kappaLAM2 = log(Na/Nb)
lambda0 = 1.55; % 以 1.55为中心波长/布拉格波长 设计
LAM_a = lambda0/(4*real(Na)); LAM_b = lambda0/(4*real(Nb));
LAM = LAM_a + LAM_b;
m = 30; % 先设计成30对光栅
Lg = LAM*m;
kappaLg = kappaLAM*m;
deltaLg = (-15:0.05:15)';
gL_scan = 0:0.02:1; 

% beta_s = pi/(2*LAM_s) + delta_s + j*g

for a=1:1:length(gL_scan)
    for b=1:1:length(deltaLg)
        beta_a = pi/(2*LAM_a) + deltaLg(b)/Lg + j*gL_scan(a)/Lg;
        beta_b = pi/(2*LAM_b) + deltaLg(b)/Lg + j*gL_scan(a)/Lg;
%         T11 = exp(-j*beta_b*LAM_b)*(exp(-j*beta_a*LAM_a) - r^2*exp(j*beta_a*LAM_a));
%         T12 = r*(exp(-j*beta_a*LAM_a) - exp(j*beta_a*LAM_a));
%         T21 = -r*(exp(-j*beta_a*LAM_a) - exp(j*beta_a*LAM_a));
%         T22 = exp(j*beta_b*LAM_b)*(exp(j*beta_a*LAM_a) - r^2*exp(-j*beta_a*LAM_a));
        T11 = exp(j*(beta_b*LAM_b+beta_a*LAM_a)) - r^2*exp(-j*(beta_a*LAM_a-beta_b*LAM_b));
        T12 = r*(exp(-j*(beta_b*LAM_b+beta_a*LAM_a) )- exp(j*(beta_a*LAM_a-beta_b*LAM_b)));
        T21 = r*(exp(j*(beta_b*LAM_b+beta_a*LAM_a)) - exp(-j*(beta_a*LAM_a-beta_b*LAM_b)));
        T22 = exp(-j*(beta_b*LAM_b+beta_a*LAM_a)) - r^2*exp(j*(beta_a*LAM_a-beta_b*LAM_b));
        T = (1/(ta_COE*tb_COE))*[T11 T12;T21 T22];
        Tg{a,b} = T^m;
        S21(a,b) = (Tg{a,b}(1))^(-1);
    end
end

%% 依次画线(画出每个kappaLg(即m)下 扫描增益的S21结果)，为了观察寻峰 
figure(1);
for a=1:1:length(gL_scan)     
    plot(deltaLg,abs(S21(a,:)));
    ylim([0 100]);% 对应文献中三组数据，依此纵坐标范围100 20 100
    text(-13.5,95,['κL_g=' num2str(kappaLg)],'FontSize',12);
    text(-13.5,90,['(Γg-α_i)L_g=' num2str(gL_scan(a)) ' ' 'a=' num2str(a)],'FontSize',12);
    pause;
end

%% 返回与阈值点坐标（gL,δL）
M = max(abs(S21(42,:))); % 29是观察最大发射（对应特征方程成立，分母为0，极点）的结果【根据观察结果更改数值】
% 三个例子分别29(简并)； 11(简并)；[16(0阶) 42(1阶)]
[row col] = find(abs(S21)==M);
disp(['κL=' num2str(kappaLg) '的DFB激光器，阈值增益与对应归一化频率位置为：']);
disp(['(gL,δL)=(' num2str(gL_scan(row)) ',' num2str(deltaLg(col)) ')']);

%% 归一化复耦和发射谱，将多条线画在一起
A=abs(S21(27,:));
B=abs(S21(22,:));
C=abs(S21(17,:));
max=max(A);A=A./max;clear max;
max=max(B);B=B./max;clear max;
max=max(C);C=C./max;clear max;
result=zeros(length(B),3);result(:,1)=A;result(:,2)=B;result(:,3)=C;
plot(deltaLg,result);legend('gL','gL-0.1','gL-0.2');
