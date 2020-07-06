% 本程序用于由GC_DFB_TMM.m复制而来，以供其他仿真

% 复数参数参考：Toshihiko Makino--Threshold Condition of DFB Semiconductor Lasers by the Local-Normal-Mode Transfer-Matrix Method: Correspondence to the Coupled-Wave Method
% 画图取值参考：Nakano, Luo, Tada--Facet reflection independent, single longitudinal mode oscillation in a GaAlAs/GaAs distributed feedback laser equipped with a gain~coupling mechanism
clear all;

%由κL=[2 2+0.2i 2+0.4i 2+0.6i 2+0.8i 2+i]算得6组[Na Nb]
Na_group = 3.6*ones(1,6); 
Nb_group = [3.3678 3.3678-0.011226i 3.3678-0.022452i 3.3677-0.033678i 3.3675-0.044903i 3.3674-0.056128i];
N_mean_group = 2./(1./real(Na_group)+1./real(Nb_group));

r_group = (Na_group - Nb_group)./(Na_group + Nb_group);
ta_COE_group = 2*Na_group./(Na_group + Nb_group); % t本该包含 场的重叠比例【式10a】，这里用t_COE表示t/重叠比例【比例是个z向无关的常数，不影响后面阈值的位置，仅影响发射谱幅值】
tb_COE_group = 2*Nb_group./(Na_group + Nb_group);
r_cleavage1 = (Na_group-1)./(Na_group+1); % 满足材料1 Na=3.6,r_cleavage = 0.565的设定(写成Na,Nb的形式是为了引用虚数)

kappaLAM_group = log(Na_group./Nb_group); % 文章中这个形式比较精确
% 解锁以下程序，获得ln(Na/Nb)与2r 近似差异的对比
% kappaLAM1 = 2*r
% kappaLAM2 = log(Na/Nb)
lambda0 = 1.55; % 以 1.55为中心波长/布拉格波长 设计
LAM_a_group = lambda0./(4*real(Na_group)); LAM_b_group = lambda0./(4*real(Nb_group));
LAM_group = LAM_a_group + LAM_b_group;
m = 30; % 先设计成30对光栅
Lg_group = LAM_group*m;
kappaLg_group = kappaLAM_group*m;

deltaLg_group = (-10:0.1:10)';
gL_scan_group = 0:0.01:1.2; 

% 传常
% beta_s = pi/(2*LAM_s) + delta_s + j*(g-αi)
% 内损（表达式中含波长[即deltaLg相关]，所以得写到循环里）
%由δ=2pi*n/lambda-2pi*n/lambda0获得λ表达式
% lambda=(deltaLg/(Lg*2*pi*N_mean) + 1/lambda0).^(-1);
% alpha_i_a = 2*pi*imag(Na)/lambda; 
% alpha_i_b = 2*pi*imag(Nb)/lambda; 
%%
for a=1:1:length(gL_scan_group)
    for b=1:1:length(deltaLg_group)
        for c=1:1:6 % 6组κLg取值
            gL_scan = gL_scan_group(a);
            deltaLg = deltaLg_group(b);
            Lg = Lg_group(c); LAM_a = LAM_a_group(c); LAM_b = LAM_b_group(c);
            r = r_group(c);Na = Na_group(c);Nb = Nb_group(c);N_mean = N_mean_group(c);
            ta_COE = ta_COE_group(c); tb_COE = tb_COE_group(c);
            r_cleavage = r_cleavage1(c);
            
            lambda=(deltaLg/(Lg*2*pi*N_mean) + 1/lambda0).^(-1);
            alpha_i_a = 2*pi*imag(Na)/lambda;
            alpha_i_b = 2*pi*imag(Nb)/lambda;
            beta_a = pi/(2*LAM_a) + deltaLg/Lg + j*(gL_scan/Lg-alpha_i_a);
            beta_b = pi/(2*LAM_b) + deltaLg/Lg + j*(gL_scan/Lg-alpha_i_b);
%         T11 = exp(-j*beta_b*LAM_b)*(exp(-j*beta_a*LAM_a) - r^2*exp(j*beta_a*LAM_a));
%         T12 = r*(exp(-j*beta_a*LAM_a) - exp(j*beta_a*LAM_a));
%         T21 = -r*(exp(-j*beta_a*LAM_a) - exp(j*beta_a*LAM_a));
%         T22 = exp(j*beta_b*LAM_b)*(exp(j*beta_a*LAM_a) - r^2*exp(-j*beta_a*LAM_a));
            T11 = exp(j*(beta_b*LAM_b+beta_a*LAM_a)) - r^2*exp(-1i*(beta_a*LAM_a-beta_b*LAM_b));
            T12 = r*(exp(-j*(beta_b*LAM_b+beta_a*LAM_a))- exp(j*(beta_a*LAM_a-beta_b*LAM_b)));
            T21 = r*(exp(j*(beta_b*LAM_b+beta_a*LAM_a)) - exp(-j*(beta_a*LAM_a-beta_b*LAM_b)));
            T22 = exp(-j*(beta_b*LAM_b+beta_a*LAM_a)) - r^2*exp(j*(beta_a*LAM_a-beta_b*LAM_b));
            T = (1/(ta_COE*tb_COE))*[T11 T12;T21 T22];
            
            % 解理面相位是随机的【应该】，这里先取左右解理面端面差pi/4【这是根据coldren书上AR/HR镀膜的图中，其中一个面λ/8相移的图取的值】
            T_left = Tmatrix_interface(0.05);%r_cleavage);
            T_right = Tmatrix_interface(-0.05);%-r_cleavage);
            
            Tg{a,b,c} = T^m;
            TT{a,b,c} = T_left*Tg{a,b,c}*Tmatrix_line((beta_b*LAM_b)/2)*T_right;
            S21(a,b,c) = (TT{a,b,c}(1))^(-1);
        end
    end
end

%% 依次画线(画出每个kappaLg(即m)下 扫描增益的S21结果)，为了观察寻峰 
figure(1);
for c=6:1:6
    for a=1:1:length(gL_scan_group)
        plot(deltaLg_group,abs(S21(a,:,c)));
        ylim([0 100]);% 对应文献中三组数据，依此纵坐标范围100 20 100
        text(-8.5,95,['κL_g=' num2str(kappaLg_group(c)) ' ' 'c=' num2str(c)],'FontSize',12);
        text(-8,90,['ΓgL_g=' num2str(gL_scan_group(a)) ' ' 'a=' num2str(a)],'FontSize',12);
        pause;
    end
end

%% 返回+1/-1阶模阈值点坐标（gL,δL）
% 观察S21得到阈值坐标[row col]
row1 = [100 75 53 29 5]; %+1模 阈值增益gL_scan(a)【短波记做+1阶】
col1 = [68 69 69 70 70]; %+1模 失谐量deltaLg(b)
row2 = [89 80 69 62 52]; %-1模 阈值增益gL_scan(a)
col2 = [135 135 136 136 136]; %-1模 失谐量deltaLg(b)
plot(deltaLg_group(col1),gL_scan_group(row1),'*b',deltaLg_group(col2),gL_scan_group(row2),'*r');
xlim([-10 10]);ylim([0 1.2]);
xlabel('\delta L_g');
ylabel('gL');
title('+1/-1阶模阈值');

%% 绘制 +1/-1阶模阈值差与κL的关系
x=imag(kappaLg_group(1:5));
y=abs(row2-row1)*0.01;%0.01是gL_scan扫描精度
plot(x,y);
xlim([0 0.5]);ylim([0 0.6]);
xlabel('κ_g_a_i_n L_g');
ylabel('ΔgL');
title('+1/-1阶模阈值差与增益耦合κ_g_a_i_n 的关系（κ_i_n_d_e_sx = 2）');


