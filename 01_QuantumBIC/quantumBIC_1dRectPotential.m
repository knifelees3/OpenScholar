clear; clc; close all;

%% ===== 参数设置 =====
L = 20;              % 空间范围 [-L, L]
N = 1200;            % 网格数（越大越接近连续谱）
x = linspace(-L, L, N)';
dx = x(2) - x(1);

V0 = 5;              % 势阱深度
a = 2;               % 势阱半宽

%% ===== 构造势能 =====
V = zeros(N,1);
V(abs(x) < a) = -V0;

%% ===== 构造哈密顿量 =====
e = ones(N,1);

% 二阶导数（有限差分）
T = (-2*diag(e) + diag(e(1:end-1),1) + diag(e(1:end-1),-1)) / dx^2;

% 哈密顿量 H = -d^2/dx^2 + V
H = -T + diag(V);

%% ===== 求本征值和本征函数 =====
[psi, E] = eig(H);
E = diag(E);

% 排序
[E, idx] = sort(E);
psi = psi(:, idx);

%% ===== 图1：能谱 =====
figure;
hold on;

% ===== 能谱绘制 =====
plot(E(1:50), 's', 'MarkerSize', 15, 'MarkerFaceColor','k');

% ===== 设置纵轴范围 =====
ylim([min(E(1:50))-1, max(E(1:50))+1]);
y_limits = ylim;  % 现在就正确了

% ===== 背景颜色 =====
% E < 0 背景灰色
patch([1 50 50 1], [y_limits(1) y_limits(1) 0 0], ...
      [0.9 0 0], 'EdgeColor','none', 'FaceAlpha',0.5);

% E > 0 背景淡蓝
patch([1 50 50 1], [0 0 y_limits(2) y_limits(2)], ...
      [0 0 0.8], 'EdgeColor','none', 'FaceAlpha',0.2);

% ===== E=0 虚线 =====
yline(0, 'r--', 'LineWidth', 1.5);

% ===== 坐标设置 =====
xlabel('State index');
ylabel('Energy E');
title('Spectrum: Discrete (E<0) vs Quasi-continuum (E>0)');
grid on;

set(gca,'FontName','times new roman','Fontsize',25,'XColor','k','YColor','k','LineWidth',1.3);

% 确保 patch 在最底层
uistack(findobj(gca,'Type','patch'),'bottom');
%% ===== 找到几个典型态 =====
% 选几个 E<0（束缚态）
bound_idx = find(E < 0);

% 选几个 E>0（连续谱中的态）
cont_idx = find(E > 0);

%% ===== 图2：束缚态波函数 =====
figure;
hold on;

% ===== 波函数 =====
for i = 1:min(3, length(bound_idx))
    n = bound_idx(i);
    plot(x, psi(:,n)/max(abs(psi(:,n))) + i, 'LineWidth', 1.5);
end

% ===== 势阱边界 =====
xline(-a, '--k', 'LineWidth', 1.5);
xline(a, '--k', 'LineWidth', 1.5);

% ===== 坐标设置 =====
xlabel('x');
ylabel('Wavefunction (offset)');
title('Bound States (E < 0): Localized');
grid on;

set(gca,'FontName','times new roman','Fontsize',25,...
    'XColor','k','YColor','k','LineWidth',1.3);
%% ===== 图3：连续谱态波函数 =====
figure;
hold on;
for i = 1:3
    n = cont_idx(i*10);  % 隔开选几个
    E(n)
    plot(x, psi(:,n)/max(abs(psi(:,n)))/2 + i, 'LineWidth', 1.5);
end
xlabel('x');
ylabel('Wavefunction (offset)');
title('Scattering States (E > 0): Extended');
grid on;
% ===== 势阱边界 =====
xline(-a, '--k', 'LineWidth', 1.5);
xline(a, '--k', 'LineWidth', 1.5);

set(gca,'FontName','times new roman','Fontsize',25,...
    'XColor','k','YColor','k','LineWidth',1.3);

%% 态密度
figure;
hold on;

% bound states
histogram(E(E<0), 20, 'FaceColor','r', 'EdgeColor','none');

% continuum
histogram(E(E>0), 100, 'FaceColor','k', 'EdgeColor','none');

xline(0,'--r','LineWidth',1.5);

legend('Bound states','Continuum');

xlabel('Energy E');
ylabel('DOS');
xlim([-5,10])

%%
figure;
hold on;

% 只画前50个态
num_states = 50;

% ===== 背景颜色 =====
E_min = min(E(1:num_states)) - 1;
E_max = max(E(1:num_states)) + 1;

% E < 0 背景灰色
patch([0 num_states num_states 0], [E_min E_min 0 0], [0.9 0.9 0.9], ...
      'EdgeColor','none', 'FaceAlpha',0.3);

% E > 0 背景淡蓝
patch([0 num_states num_states 0], [0 0 E_max E_max], [0.8 0.9 1], ...
      'EdgeColor','none', 'FaceAlpha',0.2);

% ===== 横线绘制 =====
for i = 1:num_states
    y = E(i);
    plot([0 1], [y y], 'r', 'LineWidth', 3);  % 每个态一条横线
end

% ===== E=0 虚线 =====
yline(0, 'r--', 'LineWidth', 2);

% ===== 坐标设置 =====
ylabel('Energy E');
%xlabel('State index (arbitrary)');
%title('Spectrum: Discrete (E<0) vs Quasi-continuum (E>0)');
grid on;

xlim([0 1]);  % 横坐标不重要
ylim([E_min E_max]);

set(gca,'FontName','times new roman','Fontsize',25,'XColor','k','YColor','k','LineWidth',1.3);

% 背景在最底层
uistack(findobj(gca,'Type','patch'),'bottom');