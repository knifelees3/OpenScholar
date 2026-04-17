clear; clc; close all;
% 势函数加微弱势阱

%% ===== 参数 =====
L = 50;
N = 2000;
x = linspace(1e-3, L, N)';   % 避开 r=0
dx = x(2) - x(1);

k = 1.0;
c = 1.0;

%% ===== von Neumann–Wigner 构造 =====
int_term = cumtrapz(x, sin(k*x).^2);
f = 1 ./ (1 + c * int_term);
u = sin(k*x) .* f;

u_xx = gradient(gradient(u, dx), dx);
V = k^2 + u_xx ./ u;

% 平滑
V = smoothdata(V, 'gaussian', 15);

%% ===== 加一个弱势阱（用于产生 E<0）=====
V0 = 10;
V = V - V0 * exp(-x.^2);

%% ===== 哈密顿量 =====
e = ones(N,1);
T = (-2*diag(e) + diag(e(1:end-1),1) + diag(e(1:end-1),-1)) / dx^2;
H = -T + diag(V);

%% ===== 本征求解 =====
[psi, E] = eig(H);
E = diag(E);
[E, idx] = sort(E);
psi = psi(:, idx);

%% ===== 选三类态 =====
idx_bound = find(E < 0, 1);       % bound
idx_cont  = find(E > 1, 1);       % continuum
[~, idx_bic] = min(abs(E - k^2)); % BIC

%% =========================================================
%% 🔹 图1：散点谱
figure;
plot(E(1:80), 'bs', 'MarkerSize', 12, 'MarkerFaceColor','k'); hold on;
yline(0,'r--','LineWidth',1.5);

xlabel('State index');
ylabel('Energy E');
title('Spectrum (scatter)');
grid on;

set(gca,'FontName','times new roman','Fontsize',20,'LineWidth',1.3);

%% =========================================================
%% 🔹 图2：横线谱
figure; hold on;

num_states = 2000;
E_min = min(E(1:num_states)) - 1;
E_max = max(E(1:num_states)) + 1;

% 背景
patch([0 1 1 0], [E_min E_min 0 0], [0.9 0.9 0.9], ...
    'EdgeColor','none','FaceAlpha',0.3);
patch([0 1 1 0], [0 0 E_max E_max], [0.8 0.9 1], ...
    'EdgeColor','none','FaceAlpha',0.2);

% 横线
for i = 1:num_states
    plot([0 1], [E(i) E(i)], 'r', 'LineWidth', 2);
end

yline(0,'r--','LineWidth',1.5);
xlim([0 1]); ylim([E_min E_max]);

%title('Spectrum (line representation)');
set(gca,'FontSize',20);
uistack(findobj(gca,'Type','patch'),'bottom');

%% =========================================================
%% 🔹 图3：三类波函数（对称延拓）

figure;

titles = {'Continuum state (E>0)', ...
          'BIC (E>0, localized)'};

indices = [idx_bound, idx_cont, idx_bic];

for i = 1:2
    psi_r = psi(:, 10);

    % 对称延拓 (-r, r)
    x_full = [-flipud(x); x];
    psi_full = [flipud(psi_r); psi_r];

    subplot(2,1,i);
    plot(x_full, psi_full / max(abs(psi_full)), ...
        'k', 'LineWidth', 1.5);
    grid on;

    xlabel('x');
    ylabel('\psi(x)');
    title([titles{i}, ',  E = ', num2str(E(indices(i)))]);
    
    set(gca,'FontName','times new roman','Fontsize',18,'LineWidth',1.2);
end

%% 势能函数
%% ===== 对称扩展到 (-r, r) =====
x_full = [-flipud(x); x];
V_full = [flipud(V); V];

% ===== 绘制 =====
figure;
plot(x_full, V_full, 'k', 'LineWidth', 1.8); hold on;

yline(0, 'r--', 'LineWidth', 1.5);

xlabel('x');
ylabel('V(x)');
title('von Neumann–Wigner Potential');

grid on;

set(gca,'FontName','times new roman','Fontsize',25,...
    'XColor','k','YColor','k','LineWidth',1.3);

%% ===== 计算 IPR =====
IPR = zeros(length(E),1);

for n = 1:length(E)
    psi_n = psi(:,n);
    IPR(n) = trapz(x, abs(psi_n).^4);
end

% ===== 在连续谱中找最局域的态 =====
cont_idx = find(E > 0);   % 只看连续谱

[~, max_idx] = max(IPR(cont_idx));

idx_bic = cont_idx(max_idx);

% ===== 对称延拓 =====
psi_r = psi(:, idx_bic);

x_full = [-flipud(x(2:end)); x];
psi_full = [-flipud(psi_r(2:end)); psi_r];

% ===== 绘图 =====
figure;
plot(x_full, psi_full / max(abs(psi_full)), ...
    'k', 'LineWidth', 2);

grid on;
xlabel('x');
ylabel('\psi(x)');
title(['BIC (localized), E = ', num2str(E(idx_bic))]);

set(gca,'FontName','times new roman','Fontsize',25,'LineWidth',1.3);

%%
%% ===== 选择 BIC + 4 个普通态 =====
% 计算 IPR（用于挑 BIC）
IPR = zeros(length(E),1);
for n = 1:length(E)
    psi_n = psi(:,n);
    IPR(n) = trapz(x, abs(psi_n).^4);
end

% BIC：连续谱中 IPR 最大
cont_idx = find(E > 0);
[~, max_idx] = max(IPR(cont_idx));
idx_bic = cont_idx(max_idx);

% 普通连续态：挑几个不同能量的态
% 可以均匀选连续谱区间
num_cont = 4;
cont_samples = round(linspace(1, length(cont_idx), num_cont));
idx_cont = [1,10,20,30];

% ===== 绘图：BIC + 4 个普通态 =====
figure;

% 对称延拓函数
mirror_extend = @(xvec, psi_r) deal([-flipud(psi_r(2:end)); psi_r], [-flipud(xvec(2:end)); xvec]);

% % % BIC
% [psi_full, x_full] = mirror_extend(x, psi(:, idx_bic));
% subplot(4,1,1);
% plot(x_full, psi_full / max(abs(psi_full)), 'k', 'LineWidth', 2);
% grid on;
% title(['BIC (localized), E = ', num2str(E(idx_bic))]);
% xlabel('x'); ylabel('\psi(x)');

% 普通连续态
for i = 1:4
    [psi_full, x_full] = mirror_extend(x, psi(:, idx_cont(i)));
    subplot(4,1,i);
    plot(x_full, psi_full / max(abs(psi_full)), 'k', 'LineWidth', 1.5);
    grid on;
    if i>1
    title(['Continuum state ', num2str(i), ', E = ', num2str(E(idx_cont(i)))]);
    else
    title(['Bound state ', num2str(i), ', E = ', num2str(E(idx_cont(i)))]); 
    end
    xlabel('x'); ylabel('\psi(x)');
end

% 格式
set(gcf,'Color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');