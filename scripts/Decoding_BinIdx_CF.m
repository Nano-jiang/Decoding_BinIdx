%% Decoding_BinIdx_CF

% 脚本描述：计算小鼠在 CircleField 中每个 bin 上的解码正确率

% 作者：WXW

% 创建时间：2026/02/17

% 更新时间：2026/02/19

%%
% close all
% clear
% clc

%%  Initial Settings
MiniscopeID = 0;
behavCamID = 1;

% Circle field size
width_cm = 30;
height_cm = 30;

% frame rate
fps = 30;

%% Import Data
currentFile = pwd;
targetFolder = fullfile(currentFile, 'FrameCount2');

load(fullfile(targetFolder, 'BehaviorResults.mat'));
load(fullfile(targetFolder, 'Region_OF.mat'));
load(fullfile(targetFolder, 'mouseTrajectory.mat'));

cnmfeFile = dir('*_CNMFE_Cleaned.mat');
matFileName = cnmfeFile(1).name;
load(matFileName);

%% Create the results keeping folder
tokens = regexp(currentFile, '[^\\/]+$', 'match');
if ~isempty(tokens)
    currentFileName = tokens{1};
    disp(currentFileName);

    resultsFolderName = 'DecAcc_BinIdx_CF';
    resultsFolderPath = fullfile(currentFile, resultsFolderName);

    if ~exist(resultsFolderPath, 'dir')
        mkdir(resultsFolderPath);
    end
end

%% 检验 trials 数是否满足条件
min_trialNum = 10;

% 判断 trials 是否不小于 min_trialNum
if size(BehaviorTable,1) >= min_trialNum
    fprintf('满足条件，开始计算...\n');  
else
    % 如果不满足条件，抛出错误并停止运行
    fprintf('Trials 不满足条件, 无法进行下一步计算');
    return
end

%% 小鼠的足迹和速率

% 高斯滤波器
sigma = 3; 
positions_smoothed(:, 1) = imgaussfilt(positions(:, 1), sigma);
positions_smoothed(:, 2) = imgaussfilt(positions(:, 2), sigma);

width_pixel = Region_OF(3);
height_pixel = Region_OF(4);

X_cm = abs(positions_smoothed(:,1) - Region_OF(1))* width_cm / width_pixel;
Y_cm = abs(positions_smoothed(:,2) - Region_OF(2))* height_cm / height_pixel;
positions_cm = [X_cm, Y_cm];

behavFramesNum = length(X_cm);
speed = zeros(behavFramesNum, 1);

for i = 1:behavFramesNum-1
    dx = X_cm(i+1) - X_cm(i);
    dy = Y_cm(i+1) - Y_cm(i);
    speed(i) = sqrt(dx^2 + dy^2) * fps;
end
speed(end) = speed(end-1);

% mouse's speed should no less than 3 cm/s
speed_ok_fms = find(speed >= 3);

%% 预处理 NeuronP1
[n,m] = size(NeuronP);

tsdata = importdata('timestamp.dat');

% behavCam frames → CaCam frames 的映射
FMS = correct_frames2(tsdata,MiniscopeID,behavCamID); 
% 构建速率阈值掩码，且掩码帧不能为nan，且限定在钙成像范围内
mapped = FMS(speed_ok_fms);
mapped = mapped(~isnan(mapped) & mapped>=1 & mapped<=m);
valid_fms = unique(mapped);

NeuronP1 = nan(size(NeuronP));
valid_data = NeuronP(:, valid_fms);
% origin
NeuronP1(:, valid_fms) = valid_data;

%% 整理 trials
trialNum = size(BehaviorTable,1);
trialFrame = cell(trialNum,1);
trialStamp = cell(trialNum,1);
trialStamp_cm = cell(trialNum,1);

for ti = 1:trialNum
    % frames
    trialFrame{ti} = BehaviorTable(ti,1):BehaviorTable(ti,2);
    trialFrame{ti} = trialFrame{ti}(:);
    % (x,y) in pixels
    trialStamp{ti} = positions_smoothed(BehaviorTable(ti,1):BehaviorTable(ti,2),:);
    % (x,y) in cm
    trialStamp_cm{ti} = positions_cm(BehaviorTable(ti,1):BehaviorTable(ti,2),:);
end

%% 将 Circle Field 等分成 bins
bin_angle = 36;
binNum = 360/bin_angle;
CF_center = [Region_OF(1) + Region_OF(3)/2, Region_OF(2) + Region_OF(4)/2];

trialStamp_bin = cell(trialNum,1); 

for ti = 1:trialNum
    coords = trialStamp{ti};              % n x 2，每一行为 [x, y]
    rel_coords = coords - CF_center;      % 相对于中心点的坐标

    % 1. 计算顺时针角度（范围 [0, 360)）
    angles = mod(360 - atan2d(rel_coords(:,2), rel_coords(:,1)), 360);  
    
    % 2. 将312°设为bin起点（偏移角度坐标）
    angles = mod(angles - 312, 360);      % 现在 bin 1 对应原来的 312° (奖励点位于右下角)

    % 3. 初步按逆时针编号 bin（1 到 30）
    bin_idx = floor(angles / bin_angle) + 1;

    % 4. 将 bin 编号反转为顺时针编号
    bin_idx = mod(31 - bin_idx, binNum);      % 把逆时针编号翻转成顺时针（0–29）
    bin_idx(bin_idx == 0) = binNum;           % 修正 0 为 30（模运算后）

    % 5. 存储结果
    trialStamp_bin{ti} = bin_idx;
end

% 对 trialStamp_bin 进行清洗 (防止首尾越界的情况)
for i = 1:trialNum
    currentData = trialStamp_bin{i};
    
    % 重要：确保数据是 double 类型，否则 NaN 无法正常存储
    if ~isfloat(currentData)
        currentData = double(currentData);
    end
    
    dataLen = length(currentData);
    
    % 1. 处理前 20 帧：出现 30 赋值为 NaN
    checkRange_front = 1:min(20, dataLen);
    % 在前 20 帧的范围内，找到 bin id == binNum 的逻辑索引
    idx30 = currentData(checkRange_front) == binNum;
    % 注意：修改的是 currentData 在 checkRange_front 里的对应位置
    currentData(checkRange_front(idx30)) = NaN;
    
    % 2. 处理最后 20 帧：出现 1 赋值为 NaN
    % 确定最后 20 帧的起始位置
    startIdx_back = max(1, dataLen - 19); 
    checkRange_back = startIdx_back:dataLen;
    
    % 在最后 20 帧的范围内，找到 bin id == 1 的逻辑索引
    idx1 = currentData(checkRange_back) == 1;
    currentData(checkRange_back(idx1)) = NaN;
    
    % 将修改后的数据放回 cell
    trialStamp_bin{i} = currentData;
end

%% 绘制分 bin 后的 Circle Field 以及运动轨迹

% %% 运动轨迹
figure('Color', 'w');
hold on;

x_position = positions_smoothed(:,1);
y_position = positions_smoothed(:,2);

% 1. 绘制小鼠路径：蓝色半透明 (Alpha = 0.2)
% 使用 RGBA 矢量 [R, G, B, Alpha]
for ti = 1:trialNum
    plot(x_position(BehaviorTable(ti,1):BehaviorTable(ti,2)), ...
         y_position(BehaviorTable(ti,1):BehaviorTable(ti,2)), ...
         'Color', [0, 0, 1, 0.2], 'LineWidth', 1);
end

% 2. 基础参数准备
xmin = Region_OF(1); xmax = xmin + Region_OF(3);
ymin = Region_OF(2); ymax = ymin + Region_OF(4);

% 3. 计算并绘制辐射线 (与矩形边界的交点)
theta_bins = 0:bin_angle:359; % 0 到 360 度
for i = 1:length(theta_bins)
    theta = deg2rad(theta_bins(i));
    dx = cos(theta);
    dy = sin(theta);
    
    % 与矩形四条边求交点的参数 t
    t = [];
    if dx > 0, t(end+1) = (xmax - CF_center(1)) / dx; end
    if dx < 0, t(end+1) = (xmin - CF_center(1)) / dx; end
    if dy > 0, t(end+1) = (ymax - CF_center(2)) / dy; end
    if dy < 0, t(end+1) = (ymin - CF_center(2)) / dy; end
    
    % 最小正 t 即为射线与边界的最近交点
    t_final = min(t(t > 0));
    x_end = CF_center(1) + t_final * dx;
    y_end = CF_center(2) + t_final * dy;
    
    % 绘制灰色辐射线
    plot([CF_center(1), x_end], [CF_center(2), y_end], 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
end

% 4. 绘制中心点、外框及界面美化
scatter(CF_center(1), CF_center(2), 80, 'ro', 'Filled', 'MarkerEdgeColor', 'k');
rectangle('Position', Region_OF, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1.5);

axis equal tight;
set(gca, 'YDir', 'reverse'); % 匹配图像坐标系
axis off;
title('Binned Circle Field', 'FontSize', 12);
hold off;

% 定义保存名称
gridSaveName = 'Circle_Field_Grid_Trajectory.png';
gridSavePath = fullfile(resultsFolderPath, gridSaveName);

% 使用 exportgraphics 保存，Resolution 设置为 300 DPI 保证清晰度
% 'BackgroundColor', 'current' 保持白色背景
exportgraphics(gcf, gridSavePath, 'Resolution', 300);

fprintf('Circle Field 轨迹网格图已保存至: %s\n', gridSavePath);

%% 分配每一个 bin 上的 frames
% fms_bin_behav_CF: 每一个 bin 上的 behav frames
% fms_bin_scope_CF: 每一个 bin 上的 scope frames

% 计算每个 bin 上存储的 frames
fms_bin_behav_CF = cell(binNum,1);
fms_bin_scope_CF = cell(binNum,1);

for ti = 1:trialNum
    trialfms = trialFrame{ti};
    trialfms_b2c = FMS(trialFrame{ti});

    for bi = 1:binNum
        fms_bin_behav = unique(trialfms(trialStamp_bin{ti} == bi));
        fms_bin_scope = unique(trialfms_b2c(trialStamp_bin{ti} == bi));

        fms_bin_behav_CF{bi} = [fms_bin_behav_CF{bi}; fms_bin_behav(:)];
        fms_bin_scope_CF{bi} = [fms_bin_scope_CF{bi}; fms_bin_scope(:)];
    end
end

% %% 用彩色散点核验所有 bin 上的 frames 是否正确
figure('Color', 'w', 'Name', 'Bin Verification');
hold on;

% 1. 定义颜色映射
colors = jet(binNum); 

% 2. 绘制背景轮廓 (Region_OF)
rectangle('Position', Region_OF, 'EdgeColor', [0.8 0.8 0.8], 'LineStyle', '--', 'LineWidth', 1);

% 3. 绘制足迹散点
for bi = 1:length(fms_bin_scope_CF)
    frames = fms_bin_behav_CF{bi};
    
    if ~isempty(frames)
        pts_x = positions_smoothed(frames, 1);
        pts_y = positions_smoothed(frames, 2);
        
        scatter(pts_x, pts_y, 15, colors(bi, :), 'filled', 'MarkerFaceAlpha', 0.6);
    end
end

% 4. 绘制中心点
% 标记 Circle Field 的圆心，便于验证旋转逻辑
plot(CF_center(1), CF_center(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'LineWidth', 1.5);

% 5. 界面美化与坐标设置
cb = colorbar;
colormap(jet(binNum)); 
ylabel(cb, 'Bin Index (1-30)');
clim([1 binNum]); % 锁定颜色轴范围为 1 到 30

axis equal tight;
set(gca, 'YDir', 'reverse'); % 匹配图像像素坐标系（Y轴向下）
axis off; 

title('Spatial Bin Mapping', 'FontSize', 12);

hold off;

% 定义保存文件路径
binVerifySaveName = 'Bin_Verification_Scatter.png';
binVerifySavePath = fullfile(resultsFolderPath, binVerifySaveName);

% 使用 exportgraphics 进行保存
% 'Resolution', 300 确保散点不会因为分辨率低而变成模糊的点
% 'BackgroundColor', 'none' 如果你想保持背景透明可以选 none，通常 'current' (白色) 即可
exportgraphics(gcf, binVerifySavePath, 'Resolution', 300, 'BackgroundColor', 'current');

fprintf('Bin 核验散点图已成功保存至: %s\n', binVerifySavePath);

%% 计算 NeuronDataBinAll (n * bin * trial)
NeuronDataBinAll = zeros(n,binNum,trialNum);

for ti = 1:trialNum
    trialfms_b2c = FMS(trialFrame{ti});

    for bi = 1:binNum
        NeuronData = mean(NeuronP1(:, unique(trialfms_b2c(trialStamp_bin{ti} == bi))), 2, 'omitnan');
        NeuronDataBinAll(:, bi, ti) = NeuronData;
    end
end

%% 对 NeuronDataBinAll 进行全局 z-score
[nz, bin_num_z, trial_num_z] = size(NeuronDataBinAll);

% 1. 重塑为 [n, binNum * trialNum]
% 这样每一行包含了一个神经元在整个实验过程中的所有放电观测值
data_reshaped = reshape(NeuronDataBinAll, nz, []);

% 2. 沿着第二维度（行方向）计算 Z-score
% zscore(X, flag, dim) -> dim=2 表示对每一行进行独立标准化
% 结果 data_z_reshaped 的大小依然是 [n, binNum * trialNum]
data_z_reshaped = zscore(data_reshaped, 0, 2);

% 3. 处理异常值（Silent Neurons）
% 如果某颗细胞全程不放电，其标准差为 0，Z-score 会产生 NaN
data_z_reshaped(isnan(data_z_reshaped)) = 0;

% 4. 还原回三维形状 [n, binNum, trialNum]
NeuronDataBinAll_Z = reshape(data_z_reshaped, n, bin_num_z, trial_num_z);

NeuronDataBinAll = NeuronDataBinAll_Z;

%% 对神经群体进行解码准确率的计算
% --- 第一部分：实际解码准确率 (Actual) ---
ACCbin_all_trials = zeros(trialNum, binNum);
ErrBin_all_trials = zeros(trialNum, binNum);

for ti = 1:trialNum
    % 1. 划分索引
    idxTest = ti;
    idxTrain = setdiff(1:trialNum, ti);
    
    % 2. 准备训练数据：
    % 结果维度：binNum x n (特征)
    % XTrain = mean(NeuronDataBinAll(:,:,idxTrain), 3)'; 
    % YTrain = (1:binNum)'; 
    
    % 原始维度: [nNeuron, binNum, trainTrialNum]
    % 转置并重排为: [binNum, trainTrialNum, nNeuron]
    tempTrain = permute(NeuronDataBinAll(:,:,idxTrain), [2, 3, 1]);
    % 展平为: [(binNum * trainTrialNum), nNeuron]
    XTrain = reshape(tempTrain, [binNum * length(idxTrain), size(NeuronDataBinAll, 1)]);
    YTrain = repmat((1:binNum)', length(idxTrain), 1);

    % 3. 准备测试数据：即当前的这个 trial
    % 结果维度：binNum x n
    XTest = NeuronDataBinAll(:,:,idxTest)'; 
    YTest = (1:binNum)';
    
    % 4. 训练与预测
    mdl = fitcecoc(XTrain, YTrain);
    YHat = predict(mdl, XTest);
    
    % 5. 记录当前 trial 的准确率 (每个 bin 预测对了吗？)
    ACCbin_all_trials(ti, :) = (YHat == YTest)'; 

    % 6. 记录解码误差
    bin_dist = abs(YHat - YTest); 
    angle_dist = bin_dist * bin_angle;     
    ErrBin_all_trials(ti, :) = angle_dist';
end

% 得到每个 Bin 的平均准确率 (1 x binNum)
ACCbin_Avg = mean(ACCbin_all_trials, 1);
ACC_Avg = mean(ACCbin_Avg);

ErrBin_Avg = mean(ErrBin_all_trials, 1);
Err_Avg = mean(ErrBin_Avg);              

% --- 第二部分：计算 Shuffle 解码准确率 ---
sft = 100;
rng("shuffle"); 
ACCbinSf = zeros(sft, binNum);
ErrBinSf = zeros(sft, binNum);

for sfi = 1:sft
    fprintf('sfi = %d\n', sfi);
    ACCbinSf_all_trials = zeros(trialNum, binNum);
    ErrBinSf_all_trials = zeros(trialNum, binNum);
    for ti = 1:trialNum
        idxTest = ti;
        idxTrain = setdiff(1:trialNum, ti);
       
        % XTrain = mean(NeuronDataBinAll(:,:,idxTrain), 3)';
        % YTrain_shf = randperm(binNum)'; % 打乱训练标签
        tempTrain = permute(NeuronDataBinAll(:,:,idxTrain), [2, 3, 1]);
        XTrain = reshape(tempTrain, [binNum * length(idxTrain), size(NeuronDataBinAll, 1)]);
        YTrain_shf = repmat((1:binNum)', length(idxTrain), 1);
        YTrain_shf = YTrain_shf(randperm(length(YTrain_shf)));


        XTest  = NeuronDataBinAll(:,:,idxTest)';
        YTest = (1:binNum)';
        
        mdl = fitcecoc(XTrain, YTrain_shf);
        YHat = predict(mdl, XTest);
        
        ACCbinSf_all_trials(ti, :) = (YHat == YTest)';
        ErrBinSf_all_trials(ti, :) = (abs(YHat - YTest) * bin_angle)';
    end
    % 记录这一次 shuffle 实验的平均结果
    ACCbinSf(sfi, :) = mean(ACCbinSf_all_trials, 1);
    ErrBinSf(sfi, :) = mean(ErrBinSf_all_trials, 1);
end

% 得到 Shuffle 的平均基准线 (1 x binNum)
ACCbinSf_Avg = mean(ACCbinSf, 1);
ACCSf_Avg = mean(ACCbinSf_Avg);
ErrBinSf_Avg = mean(ErrBinSf, 1);
ErrSf_Avg = mean(ErrBinSf_Avg);

%% 指定要保存的变量名
varsToSave = {
    'ACC_Avg', 'ACCSf_Avg', ...
    'ACCbin_Avg', 'ACCbinSf_Avg',...
    'Err_Avg', 'ErrSf_Avg',...
    'ErrBin_Avg', 'ErrBinSf_Avg'
    };

saveFileName = 'DecAcc_BinIdx_CF.mat';
finalSavePath = fullfile(resultsFolderPath, saveFileName);

save(finalSavePath, varsToSave{:});
fprintf('数据已成功保存至: %s\n', finalSavePath);

%% 绘制解码正确率散点图
acc = ACCbin_Avg(:);
accSf = ACCbinSf_Avg(:);
binNum = length(acc);

mu = [mean(acc), mean(accSf)]; 
sem = [std(acc)/sqrt(binNum), std(accSf)/sqrt(binNum)]; 

p_val = signrank(acc, accSf);

% 显著性标签判断
if p_val < 0.001, sig_label = '***';
elseif p_val < 0.01, sig_label = '**';
elseif p_val < 0.05, sig_label = '*';
else, sig_label = 'n.s.';
end

figure('Color', 'w', 'Units', 'pixels', 'Position', [200, 200, 400, 500]);
hold on;

% 颜色定义
color_blue = [0, 0.447, 0.741];
color_grey = [0.6, 0.6, 0.6];

% --- 绘制透明散点 (带 Jitter) ---
x_jit1 = 1 + (rand(binNum, 1) - 0.5) * 0.15; % 均匀分布抖动
x_jit2 = 2 + (rand(binNum, 1) - 0.5) * 0.15;

scatter(x_jit1, acc, 30, color_blue, 'filled', ...
    'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', 'none');
scatter(x_jit2, accSf, 30, color_grey, 'filled', ...
    'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', 'none');

% --- 绘制 SEM 误差棒 (无 Marker) ---
errorbar(1, mu(1), sem(1), 'Color', 'b', 'LineWidth', 1.5, 'CapSize', 8);
errorbar(2, mu(2), sem(2), 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'CapSize', 8);

% --- 绘制均值短线 ---
% 在 mu 处画一条短横线，宽度为 0.2 个单位
line_width_half = 0.1;
line([1-line_width_half, 1+line_width_half], [mu(1), mu(1)], 'Color', 'b', 'LineWidth', 3);
line([2-line_width_half, 2+line_width_half], [mu(2), mu(2)], 'Color', [0.3 0.3 0.3], 'LineWidth', 3);

% --- 显著性标注 (位于 SEM 之上) ---
% 找到两组中 Mean+SEM 的最高点
y_top = max(mu + sem);
y_sig_line = y_top + max(sem)*0.5; % 在误差棒上方留出间隙
y_text = y_sig_line + max(sem)*0.2;

line([1, 2], [y_sig_line, y_sig_line], 'Color', 'k', 'LineWidth', 1);
text(1.5, y_text, sig_label, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');

set(gca, ...
    'XTick', [1, 2], ...
    'XTickLabel', {'Actual', 'Shuffle'}, ...
    'XLim', [0.5, 2.5], ...
    'YLim', [0, 1], ...              % 设置 Y 轴范围
    'YTick', 0:0.25:1, ...           % 设置 Y 轴刻度间隔为 0.25
    'TickDir', 'in', ...             % 刻度线朝内
    'LineWidth', 1, ...              % 坐标轴线宽
    'Box', 'off', ...                % 关闭上边和右边框
    'TickLength', [0.02, 0.02], ...  % 刻度线长度
    'FontSize', 10);                 % 字体大小
  
ylabel('Decoding accuracy');
grid off;
hold off;

saveName = 'DecAcc_BinIdx_CF_Plot.png';
savePath = fullfile(resultsFolderPath, saveName);

if ~exist(resultsFolderPath, 'dir')
    mkdir(resultsFolderPath);
end

exportgraphics(gcf, savePath, 'Resolution', 300);
fprintf('解码准确率曲线图已保存至: %s\n', savePath);

%% 绘制解码误差 (Decoding Error) 散点图
% 提取数据：将每个 Bin 的平均误差作为散点
err = ErrBin_Avg(:);     % 实际每个 bin 的平均误差
errSf = ErrBinSf_Avg(:); % Shuffle 每个 bin 的平均误差
binNum_plot = length(err);

% 计算均值和 SEM (用于误差棒显示)
mu_err = [mean(err), mean(errSf)]; 
sem_err = [std(err)/sqrt(binNum_plot), std(errSf)/sqrt(binNum_plot)]; 

% 统计检验 (符号秩检验)
p_val_err = signrank(err, errSf);

% 显著性标签判断
if p_val_err < 0.001, sig_label = '***';
elseif p_val_err < 0.01, sig_label = '**';
elseif p_val_err < 0.05, sig_label = '*';
else, sig_label = 'n.s.';
end

% --- 开始绘图 ---
figure('Color', 'w', 'Units', 'pixels', 'Position', [250, 250, 400, 500]);
hold on;

% 颜色定义 (误差图建议使用红色或暖色调区分于准确率的蓝色)
color_actual = [0.850, 0.325, 0.098]; % 橙红色
color_shuffle = [0.6, 0.6, 0.6];     % 灰色

% 1. 绘制透明散点 (带 Jitter 抖动)
x_jit1 = 1 + (rand(binNum_plot, 1) - 0.5) * 0.15; 
x_jit2 = 2 + (rand(binNum_plot, 1) - 0.5) * 0.15;
scatter(x_jit1, err, 35, color_actual, 'filled', ...
    'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');
scatter(x_jit2, errSf, 35, color_shuffle, 'filled', ...
    'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');

% 2. 绘制 SEM 误差棒
errorbar(1, mu_err(1), sem_err(1), 'Color', color_actual, 'LineWidth', 2, 'CapSize', 8);
errorbar(2, mu_err(2), sem_err(2), 'Color', [0.3 0.3 0.3], 'LineWidth', 2, 'CapSize', 8);

% 3. 绘制均值短线
line_w = 0.12;
line([1-line_w, 1+line_w], [mu_err(1), mu_err(1)], 'Color', color_actual, 'LineWidth', 4);
line([2-line_w, 2+line_w], [mu_err(2), mu_err(2)], 'Color', [0.3 0.3 0.3], 'LineWidth', 4);

% 4. 显著性标注
y_max_data = max([err; errSf]);
y_sig_line = y_max_data + max(sem_err)*1.2; 
y_text = y_sig_line + max(sem_err)*0.3;

line([1, 2], [y_sig_line, y_sig_line], 'Color', 'k', 'LineWidth', 1);
text(1.5, y_text, sig_label, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize', 13, 'FontWeight', 'bold');

% 5. 坐标轴美化
set(gca, ...
    'XTick', [1, 2], ...
    'XTickLabel', {'Actual', 'Shuffle'}, ...
    'XLim', [0.5, 2.5], ...
    'YLim', [0, y_text * 1.1], ...     % 动态设置 Y 轴高度
    'TickDir', 'in', ...
    'LineWidth', 1.2, ...
    'Box', 'off', ...
    'FontSize', 11);

ylabel('Decoding Error (degrees)'); % 纵坐标改为误差单位
title('Decoding Distance Error', 'FontSize', 12);
grid off;
hold off;

% 保存图像
errPlotSaveName = 'DecError_BinIdx_CF_Plot.png';
errPlotSavePath = fullfile(resultsFolderPath, errPlotSaveName);
exportgraphics(gcf, errPlotSavePath, 'Resolution', 300);
fprintf('解码误差散点图已保存至: %s\n', errPlotSavePath);
