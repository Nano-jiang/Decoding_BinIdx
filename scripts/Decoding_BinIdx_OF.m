%% Decoding_BinIdx_OF

% 脚本描述：计算小鼠在 OpenField 中每个 bin 上的解码正确率

% 作者：WXW

% 创建时间：2026/02/09

% 更新时间：2026/02/19

%%
% close all
% clear
% clc

%%  Initial Settings
MiniscopeID = 0;
behavCamID = 1;

% Open field size
width_cm = 30;
height_cm = 30;

% frame rate
fps = 30;

%% Import Data
currentFile = pwd;
targetFolder = fullfile(currentFile, 'FrameCount2');

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

    resultsFolderName = 'DecAcc_BinIdx_OF';
    resultsFolderPath = fullfile(currentFile, resultsFolderName);

    if ~exist(resultsFolderPath, 'dir')
        mkdir(resultsFolderPath);
    end
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

for ki = 1:behavFramesNum-1
    dx = X_cm(ki+1) - X_cm(ki);
    dy = Y_cm(ki+1) - Y_cm(ki);
    speed(ki) = sqrt(dx^2 + dy^2) * fps;
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

%% 将 Open Field 等分成 4*4 的bins

bin_num = [4 4];
binNum = bin_num(1)*bin_num(2);

bin_width = Region_OF(:,3)/bin_num(1);
bin_height = Region_OF(:,4)/bin_num(2);

% 假设轨迹数据范围
x_OF = [Region_OF(:,1), Region_OF(:,1)+Region_OF(:,3)];
y_OF = [Region_OF(:,2),Region_OF(:,2)+Region_OF(:,4)];

% 生成网格中心点
x_centers_OF = x_OF(1)+bin_width/2 : bin_width : x_OF(2)-bin_width/2;
y_centers_OF = y_OF(1)+bin_height/2 : bin_height : y_OF(2)-bin_height/2;

% 生成网格矩阵（每一对 (X,Y) 表示一个 bin 的中心）
[Xgrid_OF, Ygrid_OF] = meshgrid(x_centers_OF, y_centers_OF);

%% 绘制分 bin 后的 Open Field 以及运动轨迹

figure('Color', 'w');
hold on;

% 1. 绘制运动轨迹：蓝色半透明 (Alpha = 0.2)
% 放在最底层，避免遮挡网格线
plot(positions_smoothed(:,1), positions_smoothed(:,2), ...
     'Color', [0, 0, 1, 0.2], 'LineWidth', 1);

% 2. 绘制 Open Field 网格和外框
% 使用黑色 ('k') 绘制每一个 Bin 的边界
for ti = 1:numel(Xgrid_OF)
    x_rect = Xgrid_OF(ti) - bin_width/2;
    y_rect = Ygrid_OF(ti) - bin_height/2;
    rectangle('Position', [x_rect, y_rect, bin_width, bin_height], ...
              'EdgeColor', 'k', 'LineWidth', 0.8);
end

% 3. 绘制特定点标记
% 起始点 (第1帧) 和 标记点 (第20帧)
scatter(positions_smoothed(1,1), positions_smoothed(1,2), 40, 'ro', 'filled');
scatter(positions_smoothed(20,1), positions_smoothed(20,2), 60, 'r+', 'LineWidth', 1.5);

% 4. 界面美化与坐标设置
axis equal tight;
set(gca, 'YDir', 'reverse'); % 匹配图像坐标系
axis off;
title('Open Field Binned Regions', 'FontSize', 12);

hold off;

% 定义保存名称
gridSaveName = 'Open_Field_Grid_Trajectory.png';
gridSavePath = fullfile(resultsFolderPath, gridSaveName);

% 使用 exportgraphics 保存，Resolution 设置为 300 DPI 保证清晰度
% 'BackgroundColor', 'current' 保持白色背景
exportgraphics(gcf, gridSavePath, 'Resolution', 300);

fprintf('Open Field 轨迹网格图已保存至: %s\n', gridSavePath);

%% 分配每一个 bin 上的 frames
% fms_bin_behav_OF: 每一个 bin 上的 behav frames
% fms_bin_scope_OF: 每一个 bin 上的 scope frames

% 计算每个 bin 上存储的 frames
fms_bin_behav_OF = cell(binNum,1);
fms_bin_scope_OF = cell(binNum,1);

centers_OF = [Xgrid_OF(:), Ygrid_OF(:)];

coords = positions_smoothed;  % n x 2, 每一行是 [x, y]

centers = centers_OF;
% 计算每个点与所有 bin 中心的距离（用 pdist2）
dists = pdist2(coords, centers);  % [nFrame x nBin]

% 距离小鼠位置最近的bin中心编号被当作小鼠此时所在的bin，min函数总是返回靠前的idx
[~, minIdx] = min(dists, [], 2); 

% 小鼠每一帧所在的 BinID
position_bin_idx = minIdx;

behav_frames = 1:size(positions,1);
for bi = 1:binNum
    fms_bin_behav = unique(behav_frames(position_bin_idx == bi));
    fms_bin_scope = unique(FMS(behav_frames(position_bin_idx == bi)));

    fms_bin_behav_OF{bi} = [fms_bin_behav_OF{bi}; fms_bin_behav(:)];
    fms_bin_scope_OF{bi} = [fms_bin_scope_OF{bi}; fms_bin_scope(:)];
end

% %% 用彩色散点核验所有 bin 上的 frames 是否正确
figure('Color', 'w', 'Name', 'Bin Verification');
hold on;

% 1. 定义颜色映射
colors = jet(binNum); 

% 2. 绘制背景轮廓 (Region_OF)
rectangle('Position', Region_OF, 'EdgeColor', [0.8 0.8 0.8], 'LineStyle', '--', 'LineWidth', 1);

% 3. 绘制足迹散点
for bi = 1:length(fms_bin_scope_OF)
    frames = fms_bin_behav_OF{bi};
    
    if ~isempty(frames)
        pts_x = positions_smoothed(frames, 1);
        pts_y = positions_smoothed(frames, 2);
        
        scatter(pts_x, pts_y, 15, colors(bi, :), 'filled', 'MarkerFaceAlpha', 0.6);
    end
end

% 4. 界面美化与坐标设置
cb = colorbar;
colormap(jet(binNum)); 
ylabel(cb, 'Bin Index (1-225)');
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

%% 对神经群体进行解码准确率的计算

kfolds = 5; % divide the total session into 5 folds
fold_length = round(m/kfolds);
folds = cell(kfolds,1);
for ki = 1:kfolds
    if ki < kfolds
        folds{ki} = (ki-1)*fold_length +1 : ki*fold_length;
    else
        folds{ki} = (ki-1)*fold_length +1 : m;
    end
end

% --- 准备工作：计算每个 bin 的中心点坐标 ---
bin_width_cm = 30/bin_num(1);
bin_height_cm = 30/bin_num(2);

bin_centers = zeros(binNum, 2); % [x_cm, y_cm]
for bi = 1:binNum
    [row, col] = ind2sub(bin_num, bi); 
    bin_centers(bi, 1) = (col - 0.5) * bin_width_cm;
    bin_centers(bi, 2) = (row - 0.5) * bin_height_cm;
end

% --- 第一部分：计算实际解码准确率 (Actual) ---
ACCbin_kfolds = zeros(kfolds, binNum);
ErrBin_kfolds = zeros(kfolds, binNum);

for ki = 1:kfolds
    fmsTest = folds{ki}; 
    fmsTrain = setdiff(1:m, fmsTest); 
    
    % 构建数据 (X: binNum x n)
    [tempXTrain, tempXTest] = deal(zeros(binNum, n));
    for bi = 1:binNum
        v_train = intersect(fms_bin_scope_OF{bi}, fmsTrain, 'stable');
        v_test  = intersect(fms_bin_scope_OF{bi}, fmsTest, 'stable');
        if ~isempty(v_train)
            tempXTrain(bi,:) = mean(NeuronP1(:, v_train), 2, 'omitnan')'; 
        end
        if ~isempty(v_test)
            tempXTest(bi,:)  = mean(NeuronP1(:, v_test), 2, 'omitnan')'; 
        end
    end
    
    YTrain = (1:binNum)'; 
    mdl = fitcecoc(tempXTrain, YTrain);
    YHat = predict(mdl, tempXTest);
    ACCbin_kfolds(ki, :) = (YHat' == (1:binNum)); 

    % 计算解码误差 (物理距离)
    pred_coords = bin_centers(YHat, :);        % 预测的 [x, y]
    true_coords = bin_centers(1:binNum, :);    % 真实的 [x, y]
    
    % 计算欧几里得距离: sqrt((x1-x2)^2 + (y1-y2)^2)
    dist_error = sqrt(sum((pred_coords - true_coords).^2, 2));
    ErrBin_kfolds(ki, :) = dist_error';
end
% 最终得到 1 x binNum 的实际平均准确率 (转为百分比)
ACCbin_Avg = mean(ACCbin_kfolds, 1); 
ACC_Avg = mean(ACCbin_Avg);
ErrBin_Avg = mean(ErrBin_kfolds, 1); 
Err_Avg = mean(ErrBin_Avg);

% --- 第二部分：计算 Shuffle 解码准确率 ---
sft = 100;
rng("shuffle"); 
ACCbinSf = zeros(sft, binNum);
ErrBinSf = zeros(sft, binNum);

for sfi = 1:sft
    fprintf('sfi = %d\n', sfi);
    ACCbinSf_kfolds = zeros(kfolds, binNum);
    ErrBinSf_kfolds = zeros(kfolds, binNum);

    for ki = 1:kfolds
        fmsTest = folds{ki}; 
        fmsTrain = setdiff(1:m, fmsTest); 
    
        % 构建数据 (X: binNum x n)
        [tempXTrain, tempXTest] = deal(zeros(binNum, n));
        for bi = 1:binNum
            v_train = intersect(fms_bin_scope_OF{bi}, fmsTrain, 'stable');
            v_test  = intersect(fms_bin_scope_OF{bi}, fmsTest, 'stable');
            if ~isempty(v_train), tempXTrain(bi,:) = mean(NeuronP1(:, v_train), 2, 'omitnan')'; end
            if ~isempty(v_test),  tempXTest(bi,:)  = mean(NeuronP1(:, v_test), 2, 'omitnan')'; end
        end

        % 打乱标签
        YTrain_shuffled = randperm(binNum)'; % 随机打乱 1:binNum 的顺序
        
        mdl = fitcecoc(tempXTrain, YTrain_shuffled);
        YHat = predict(mdl, tempXTest);
        
        % 注意：测试集的标签也要对应打乱，或者保持 YHat 与原标签比对，效果是一样的
        % 这里常用逻辑是：看模型在乱序训练后，是否还能“瞎猜”对原始标签
        ACCbinSf_kfolds(ki, :) = (YHat' == (1:binNum)); 

        % 计算 Shuffle 物理误差 (计算预测 Bin 中心与真实 Bin 中心的距离)
        pred_coords_sf = bin_centers(YHat, :);
        true_coords_sf = bin_centers(1:binNum, :);
        dist_error_sf = sqrt(sum((pred_coords_sf - true_coords_sf).^2, 2));
        ErrBinSf_kfolds(ki, :) = dist_error_sf';
    end
    ACCbinSf(sfi, :) = mean(ACCbinSf_kfolds, 1);
    ErrBinSf(sfi, :) = mean(ErrBinSf_kfolds, 1);
end

% 得到用于绘图的 Shuffle 平均值
ACCbinSf_Avg = mean(ACCbinSf, 1); 
ACCSf_Avg = mean(ACCbinSf_Avg);
ErrBinSf_Avg = mean(ErrBinSf, 1); 
ErrSf_Avg = mean(ErrBinSf_Avg);   

%% 指定要保存的变量名
varsToSave = {
    'ACC_Avg', 'ACCSf_Avg',...
    'ACCbin_Avg', 'ACCbinSf_Avg',...
    'Err_Avg', 'ErrSf_Avg',...
    'ErrBin_Avg', 'ErrBinSf_Avg'
    };

saveFileName = 'DecAcc_BinIdx_OF.mat';
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

saveName = 'DecAcc_BinIdx_OF_Plot.png';
savePath = fullfile(resultsFolderPath, saveName);

if ~exist(resultsFolderPath, 'dir')
    mkdir(resultsFolderPath);
end

exportgraphics(gcf, savePath, 'Resolution', 300);
fprintf('解码准确率曲线图已保存至: %s\n', savePath);

%% 绘制解码误差 (Decoding Error) 散点图
err = ErrBin_Avg(:);      % 实际每个 bin 的平均误差 (cm)
errSf = ErrBinSf_Avg(:);  % Shuffle 后每个 bin 的平均误差 (cm)
binNum = length(err);

% 计算均值与标准误 (SEM)
mu_err = [mean(err), mean(errSf)]; 
sem_err = [std(err)/sqrt(binNum), std(errSf)/sqrt(binNum)]; 

% 统计检验 (Wilcoxon 符号秩检验)
p_val_err = signrank(err, errSf);

% 显著性标签判断
if p_val_err < 0.001, sig_label = '***';
elseif p_val_err < 0.01, sig_label = '**';
elseif p_val_err < 0.05, sig_label = '*';
else, sig_label = 'n.s.';
end

figure('Color', 'w', 'Units', 'pixels', 'Position', [650, 200, 400, 500]); % 偏移位置，不与准确率图重叠
hold on;

% 颜色定义 (误差图建议换一种颜色以示区分，比如橙色)
color_actual = [0.85, 0.325, 0.098]; % 砖红色/橙色
color_grey = [0.6, 0.6, 0.6];

% --- 绘制透明散点 (带 Jitter) ---
x_jit1 = 1 + (rand(binNum, 1) - 0.5) * 0.15;
x_jit2 = 2 + (rand(binNum, 1) - 0.5) * 0.15;

scatter(x_jit1, err, 30, color_actual, 'filled', ...
    'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', 'none');
scatter(x_jit2, errSf, 30, color_grey, 'filled', ...
    'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', 'none');

% --- 绘制 SEM 误差棒 ---
errorbar(1, mu_err(1), sem_err(1), 'Color', color_actual, 'LineWidth', 1.5, 'CapSize', 8);
errorbar(2, mu_err(2), sem_err(2), 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'CapSize', 8);

% --- 绘制均值短横线 ---
lw_half = 0.1;
line([1-lw_half, 1+lw_half], [mu_err(1), mu_err(1)], 'Color', color_actual, 'LineWidth', 3);
line([2-lw_half, 2+lw_half], [mu_err(2), mu_err(2)], 'Color', [0.3 0.3 0.3], 'LineWidth', 3);

% --- 显著性标注 (位于最高点上方) ---
y_max_val = max([err; errSf]); % 找到所有点中的最高值
y_sig_line = y_max_val * 1.05;  % 线上方留 5% 空间
y_text = y_sig_line * 1.02;

line([1, 2], [y_sig_line, y_sig_line], 'Color', 'k', 'LineWidth', 1);
text(1.5, y_text, sig_label, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');

% --- 坐标轴美化 ---
set(gca, ...
    'XTick', [1, 2], ...
    'XTickLabel', {'Actual', 'Shuffle'}, ...
    'XLim', [0.5, 2.5], ...
    'YLim', [0, ceil(y_text/5)*5], ... % 动态调整 Y 轴上限，按 5 取整
    'TickDir', 'in', ...
    'LineWidth', 1, ...
    'Box', 'off', ...
    'TickLength', [0.02, 0.02], ...
    'FontSize', 10);

ylabel('Decoding Error (cm)');
grid off;
hold off;

% 保存图片
saveNameErr = 'DecError_BinIdx_OF_Plot.png';
savePathErr = fullfile(resultsFolderPath, saveNameErr);
exportgraphics(gcf, savePathErr, 'Resolution', 300);
fprintf('解码误差散点图已保存至: %s\n', savePathErr);
