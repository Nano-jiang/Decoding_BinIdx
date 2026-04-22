%% Decoding_BinIdx_Tmaze

% 脚本描述：计算小鼠在 T-maze 中每个 bin 上的解码正确率

% 作者：WXW

% 创建时间：2026/02/17

% 更新时间：2026/02/20
% Left & Right 分别解码

%%
% close all
% clear
% clc

%% Initial Settings
MiniscopeID = 0;
behavCamID = 1;

% T-maze size
width_cm = 47;
height_cm = 15+21+16;

% Reward bowl diameter (cm)
Reward_diameter = 3;

% frame rate
fps = 30;

% binNum = [center stem, left/right arm]
binNum = [5 5];

%% Import Data
currentFile = pwd;
targetFolder = fullfile(currentFile, 'FrameCount2');

load(fullfile(targetFolder, 'Boundary.mat'));
load(fullfile(targetFolder, 'TmazeRegion.mat'));
load(fullfile(targetFolder, 'mouseTrajectory.mat'));
load(fullfile(targetFolder, 'BehaviorResults.mat'));

% For overexp mouse data
if exist('overexpFramesConfirm.mat', 'file')
    load('overexpFramesConfirm.mat');   % 变量直接进工作区
    BehaviorTable = BehaviorTable_noOverexpTrials;
end

cnmfeFile = dir('*_CNMFE_Cleaned.mat');
matFileName = cnmfeFile(1).name;
load(matFileName);

%% Create the results keeping folder
tokens = regexp(currentFile, '[^\\/]+$', 'match');
if ~isempty(tokens)
    currentFileName = tokens{1};
    disp(currentFileName);

    resultsFolderName = 'DecAcc_BinIdx_Tmaze';
    resultsFolderPath = fullfile(currentFile, resultsFolderName);

    if ~exist(resultsFolderPath, 'dir')
        mkdir(resultsFolderPath);
    end
end

%% 检验 trials 数是否满足条件
min_trialNum = 8;

binNum_all = binNum(1) + binNum(2)*2;

left_trials_id = find(BehaviorTable(:,8) == 1);
right_trials_id = find(BehaviorTable(:,8) == 2);

left_trials_num = length(left_trials_id);
right_trials_num = length(right_trials_id);

% 判断 trials 是否不小于 min_trialNum
if left_trials_num >= min_trialNum && right_trials_num >= min_trialNum
    fprintf('满足条件，开始计算...\n');  
else
    % 如果不满足条件，抛出错误并停止运行
    fprintf('至少有一侧的 trials 不满足条件, 无法进行下一步计算');
    return
end

%% 小鼠的足迹和速率

% 高斯滤波器
sigma = 3; 
positions_smoothed(:, 1) = imgaussfilt(positions(:, 1), sigma);
positions_smoothed(:, 2) = imgaussfilt(positions(:, 2), sigma);

X1 = positions_smoothed(:,1);
Y1 = positions_smoothed(:,2);

X_origin = Region_S3(1);
Y_origin = Region_S1(2) + Region_S1(4);

width_pixel = Region_S3(3);
height_pixel = Region_S1(4) + Region_S2(4) + Region_S3(4);

X_cm = abs(X1 - X_origin) ./ width_pixel .* width_cm;
Y_cm = abs(Y1 - Y_origin) ./ height_pixel .* height_cm;
positions_cm = [X_cm, Y_cm];

behavFramesNum = length(X_cm);
speed = zeros(behavFramesNum, 1);

for ni = 1:behavFramesNum-1
    dx = X_cm(ni+1) - X_cm(ni);
    dy = Y_cm(ni+1) - Y_cm(ni);
    speed(ni) = sqrt(dx^2 + dy^2) * fps;
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
trialNum = length(BehaviorTable);
trialFrame = cell(trialNum,3);
trialStamp = cell(trialNum,3);
trialStamp_cm = cell(trialNum,3);

for si = 1:3
    for ti = 1:trialNum
        % frames
        trialfms = BehaviorTable(ti,si*2-1):BehaviorTable(ti,si*2);
        trialFrame{ti,si} = trialfms(:);
        % (x,y) in pixels
        trialStamp{ti,si} = positions_smoothed(trialfms,:);
        % (x,y) in cm
        trialStamp_cm{ti,si} = positions_cm(trialfms,:);
    end
end

%% 将 T-maze 划分成不同的 bins
% Start box
start_box_w = linspace(Region_S1(1), Region_S1(1) + Region_S1(3), 4);
start_box_h = linspace(boundary_S1S2, boundary_S1S2 + Region_S1(4), 4);

% Center stem
center_stem = linspace(boundary_S2S3, boundary_S1S2, binNum(1) + 1);

% Tmaze_axis (计算中心轴)
Tmaze_axis = Region_S2(1) + Region_S2(3)/2;

% Left L
x_left_end = Region_S3(1) + Reward_diameter ./ width_cm .* width_pixel;

% Right L
x_right_end = Region_S3(1) + Region_S3(3) - Reward_diameter ./ width_cm .* width_pixel;

left_arm = linspace(x_left_end, Tmaze_axis, binNum(2) + 1);
right_arm = linspace(Tmaze_axis, x_right_end, binNum(2) + 1);

trialStamp_bin = cell(trialNum, 3);

for ti = 1:trialNum
    % 均保留靠近物体侧边界上的点
    % center stem
    si = 2;
    y_coords = trialStamp{ti,si}(:,2);
    bin_idx = discretize(y_coords, center_stem, 'IncludedEdge', 'left'); 
    trialStamp_bin{ti,si} = binNum(1) - bin_idx + 1;

    % left-arm and right-arm
    si = 3;
    x_coords = trialStamp{ti,si}(:,1);
    direction = BehaviorTable(ti,8);

    if direction == 1
        bin_idx = discretize(x_coords, left_arm, 'IncludedEdge', 'right');
        bin_idx = binNum(2) - bin_idx + 1;
    elseif direction == 2 
        bin_idx = discretize(x_coords, right_arm, 'IncludedEdge', 'left');
        bin_idx = bin_idx + binNum(2); % 计算 Decoding_BinIdx_Tmaze 使用
    end
    bin_idx = bin_idx + binNum(1);
    trialStamp_bin{ti,si} = bin_idx;
end

%% 绘制分 bin 后的 T-maze 以及 left/right trials 轨迹 

% 绘制 T-maze 轮廓
figure;
rectangle('Position', Region_S1, 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1.2); hold on;
rectangle('Position', Region_S2, 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1.2);
rectangle('Position', Region_S3, 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1.2);

% --- 绘制 Center Stem (Region_S2 的水平等分线) ---
% center_stem 定义了 y 轴的刻度
x_s2_range = [Region_S2(1), Region_S2(1) + Region_S2(3)];
for ni = 1:length(center_stem)
    % 绘制水平线，限制在 Region_S2 的宽度内
    line(x_s2_range, [center_stem(ni), center_stem(ni)], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
end

% --- 绘制 Left Arm (Region_S3 左侧部分的垂直等分线) ---
% y 轴限制在 Region_S3 的高度内
y_s3_range = [Region_S3(2), Region_S3(2) + Region_S3(4)];
for ni = 1:length(left_arm)
    % 绘制蓝色垂直线
    line([left_arm(ni), left_arm(ni)], y_s3_range, 'Color', 'b', 'LineWidth', 0.8);
end

% --- 绘制 Right Arm (Region_S3 右侧部分的垂直等分线) ---
for ni = 1:length(right_arm)
    % 绘制红色垂直线
    line([right_arm(ni), right_arm(ni)], y_s3_range, 'Color', 'r', 'LineWidth', 0.8);
end

% --- 绘制 Start Box 网格 (可选) ---
for ni = 1:length(start_box_w)
    line([start_box_w(ni), start_box_w(ni)], [Region_S1(2), Region_S1(2)+Region_S1(4)], 'Color', [0.8 0.8 0.8]);
end
for ni = 1:length(start_box_h)
    line([Region_S1(1), Region_S1(1)+Region_S1(3)], [start_box_h(ni), start_box_h(ni)], 'Color', [0.8 0.8 0.8]);
end

% 界面美化
axis equal;
grid off; % 关闭默认网格，只看我们画的自定义网格
title('Binned T-maze');
xlabel('Pixels (X)'); ylabel('Pixels (Y)');

% 标注区域名称
% text(Region_S1(1)+5, Region_S1(2)+5, 'Start Box', 'Color', [0.5 0.5 0.5]);
% text(Region_S2(1)+5, Region_S2(2)+5, 'Center Stem', 'Color', 'k');
% text(Region_S3(1)+5, Region_S3(2)+5, 'Arms', 'Color', 'k');

for ti = 1:size(BehaviorTable,1)
    direction = BehaviorTable(ti,8);
    x_trial = positions_smoothed(BehaviorTable(ti,3):BehaviorTable(ti,6),1);
    y_trial = positions_smoothed(BehaviorTable(ti,3):BehaviorTable(ti,6),2);
    
    if direction == 1
        lineColor = [0, 0, 1, 0.15]; % 透明蓝色 (Left)
    elseif direction == 2
        lineColor = [1, 0, 0, 0.15]; % 透明红色 (Right)
    else
        lineColor = [0.5, 0.5, 0.5, 0.1]; % 其他情况用灰色
    end
    
    plot(x_trial, y_trial, 'Color', lineColor, 'LineWidth', 1);
    hold on
end

axis equal;
grid off;
set(gca, 'Color', 'w'); % 背景设为白色，使透明度更清晰
set(gca, 'YDir', 'reverse');

% 定义保存名称
gridSaveName = 'Tmaze_Grid_Trajectory.png';
gridSavePath = fullfile(resultsFolderPath, gridSaveName);

% 使用 exportgraphics 保存，Resolution 设置为 300 DPI 保证清晰度
% 'BackgroundColor', 'current' 保持白色背景
exportgraphics(gcf, gridSavePath, 'Resolution', 300);

fprintf('Tmaze 轨迹网格图已保存至: %s\n', gridSavePath);

%% 分配 所有 trials 每一个 bin 上的 frames
% fms_bin_behav: 用于后续 pii 的计算
% fms_bin_scope: 用于后续 mean firing rate 的计算

% 计算 所有 trials 每个 bin 上存储的 frames
fms_bin_behav = cell(binNum_all,1);
fms_bin_scope = cell(binNum_all,1);

for ti = 1:trialNum
    for si = 2:3
        if si == 2
            bin_num = binNum(1);
            bin_idxs = 1:binNum(1);
        elseif si == 3
            bin_num = binNum(2)*2;
            bin_idxs = (binNum(1) + 1) : binNum_all;
        end

        trialfms = trialFrame{ti,si};
        trialfms_b2c = FMS(trialFrame{ti,si});

        for bi = 1:bin_num
            bii = bin_idxs(bi);
            temp_bin_behav = unique(trialfms(trialStamp_bin{ti,si} == bii));
            temp_bin_scope = unique(trialfms_b2c(trialStamp_bin{ti,si} == bii));

            fms_bin_behav{bii} = [fms_bin_behav{bii}; temp_bin_behav(:)];
            fms_bin_scope{bii} = [fms_bin_scope{bii}; temp_bin_scope(:)];
        end
    end
end

% %% 用彩色散点核验所有 bin 上的 frames 是否正确
figure; hold on;

% 定义颜色映射
colors = jet(binNum_all); 

% --- 左侧子图：Left Trials ---
for bi = 1:length(fms_bin_behav)
    frames = fms_bin_behav{bi};
    
    if ~isempty(frames)
        pts_x = positions_cm(frames, 1);
        pts_y = positions_cm(frames, 2);
        
        % 绘制散点：同一Bin颜色相同
        scatter(pts_x, pts_y, 15, colors(bi, :), 'filled', 'MarkerFaceAlpha', 0.6);
    end
end
axis equal; grid on; 
xlabel('X (cm)'); ylabel('Y (cm)');
title('Left Trials: Binned Coordinates');
cb1 = colorbar; 
colormap(jet); 
ylabel(cb1, 'Bin Index');

% 统一坐标轴范围（可选，方便对比）
linkaxes(findall(gcf, 'Type', 'axes'), 'xy');

% 定义保存文件路径
binVerifySaveName = 'Bin_Verification_Scatter.png';
binVerifySavePath = fullfile(resultsFolderPath, binVerifySaveName);

% 使用 exportgraphics 进行保存
% 'Resolution', 300 确保散点不会因为分辨率低而变成模糊的点
% 'BackgroundColor', 'none' 如果你想保持背景透明可以选 none，通常 'current' (白色) 即可
exportgraphics(gcf, binVerifySavePath, 'Resolution', 300, 'BackgroundColor', 'current');

fprintf('Bin 核验散点图已成功保存至: %s\n', binVerifySavePath);

%% 计算 NeuronDataBinAll (n * bin * trial)
NeuronDataBinS2 = zeros(n,binNum(1),trialNum);
NeuronDataBinS3 = zeros(n,binNum(2)*2,trialNum);

for ti = 1:trialNum
    for si = 2:3
        if si == 2
            bin_num = binNum(1);
            bin_idxs = 1:binNum(1);
        elseif si == 3
            bin_num = binNum(2) *2;
            bin_idxs = (binNum(1) + 1) : (binNum(1)+binNum(2)*2);
        end

        trialfms_b2c = FMS(trialFrame{ti,si});

        for bi = 1:bin_num
            bii = bin_idxs(bi);
            NeuronData = mean(NeuronP1(:, unique(trialfms_b2c(trialStamp_bin{ti,si} == bii))), 2, 'omitnan');

            if si == 2
                NeuronDataBinS2(:, bi, ti) = NeuronData;
            elseif si == 3
                NeuronDataBinS3(:, bi, ti) = NeuronData;
            end
        end
    end
end

% NeuronDataBinS2: center stem 5 bins
% NeuronDataBinS3: Left arm 5 bins + Right arm 5 bins = 10 bins
% NeuronDataBinAll: 15 bins
NeuronDataBinAll = cat(2, NeuronDataBinS2, NeuronDataBinS3);

%% 对 NeuronDataBinAll 进行全局 z-score (手动处理 NaN 掩码)
[nz, bin_num_z, trial_num_z] = size(NeuronDataBinAll);
data_reshaped = reshape(NeuronDataBinAll, nz, []);

% 初始化一个全为 NaN 的结果矩阵，保证原始缺失位置依然是 NaN
data_z_reshaped = nan(size(data_reshaped));

% 遍历每一个神经元（行）
for i = 1:nz
    row_data = data_reshaped(i, :);
    
    % 创建掩码：找到当前行中不是 NaN 的位置
    mask = ~isnan(row_data);
    
    % 如果该行有有效数据
    if any(mask)
        valid_data = row_data(mask);
        
        % 计算该行有效数据的标准差
        sd = std(valid_data);
        
        if sd > 0
            % 只有标准差不为 0 时才计算 Z-score
            % zscore(valid_data) 默认对列操作，如果 row_data 是行向量，它会自动处理
            data_z_reshaped(i, mask) = zscore(valid_data);
        else
            % 如果标准差为 0（即非 NaN 数据全部相等），设为 0
            data_z_reshaped(i, mask) = 0;
        end
    end
end

% 4. 还原回三维形状 [n, binNum, trialNum]
NeuronDataBinAll_Z = reshape(data_z_reshaped, nz, bin_num_z, trial_num_z);
NeuronDataBinAll = NeuronDataBinAll_Z;

%% 对神经群体进行解码准确率的计算
binNum_frm = binNum(1)+binNum(2);
bin_idx_left = 1:binNum_frm;
bin_idx_right = [1:binNum(1), binNum_frm + (1:binNum(2))];

%% --- Left：实际解码准确率 (Actual) ---

% T-maze size
% width_cm = 47;
% height_cm = 15+21+16;

% Reward bowl diameter (cm)
% Reward_diameter = 3;

CenterStem_cm = 21;
Arm_cm = (width_cm - Reward_diameter*2)/2;

% 计算每一个 Bin 的大小用于解码误差的计算
bin_size_S2 = CenterStem_cm / binNum(1);
bin_size_S3 = Arm_cm / binNum(2);
bin_size = (bin_size_S2 + bin_size_S3) /2; % 4.15 cm

bin_num = length(bin_idx_left);
ACCbin_left_trials = zeros(left_trials_num, bin_num);
ErrBin_left_trials = zeros(left_trials_num, bin_num);

for ti = 1:left_trials_num
    % 1. 划分索引
    tii = left_trials_id(ti);
    idxTest = tii;
    idxTrain = setdiff(left_trials_id, idxTest);

    % 2. 准备训练数据：
    % 结果维度：binNum x n (特征)
    % XTrain = mean(NeuronDataBinAll(:,bin_idx_left,idxTrain), 3)';
    % YTrain = bin_idx_left';

    % 原始维度: [nNeuron, binNum, trainTrialNum]
    % 转置并重排为: [binNum, trainTrialNum, nNeuron]
    tempTrain = permute(NeuronDataBinAll(:,bin_idx_left,idxTrain), [2, 3, 1]);
    % 展平为: [(binNum * trainTrialNum), nNeuron]
    XTrain = reshape(tempTrain, [length(bin_idx_left) * length(idxTrain), size(NeuronDataBinAll, 1)]);
    YTrain = repmat((bin_idx_left)', length(idxTrain), 1);
    
    % 3. 准备测试数据：即当前的这个 trial
    % 结果维度：binNum x n 
    XTest = NeuronDataBinAll(:,bin_idx_left,idxTest)'; 
    YTest = bin_idx_left';
    
    % 4. 训练与预测
    mdl = fitcecoc(XTrain, YTrain);
    YHat = predict(mdl, XTest);
    
    % 5. 记录当前 trial 的准确率 (每个 bin 预测对了吗？)
    ACCbin_left_trials(ti, :) = (YHat == YTest)'; 

    % 6. 记录解码误差
    ErrBin_left_trials(ti, :) = (abs(YHat - YTest) * bin_size)';
end

% 得到每个 Bin 的平均准确率 (1 x binNum)
ACCbin_left_Avg = mean(ACCbin_left_trials, 1);
ACC_left_Avg = mean(ACCbin_left_Avg);

ErrBin_left_Avg = mean(ErrBin_left_trials, 1);
Err_left_Avg = mean(ErrBin_left_Avg);         

% --- Left：计算 Shuffle 解码准确率 ---
sft = 50;
rng("shuffle"); 
ACCbinSf_left = zeros(sft, bin_num);
ErrBinSf_left = zeros(sft, bin_num);

for sfi = 1:sft
    fprintf('sfi = %d\n', sfi);
    ACCbinSf_left_trials = zeros(left_trials_num, bin_num);
    ErrBinSf_left_trials = zeros(left_trials_num, bin_num);
    for ti = 1:left_trials_num
        tii = left_trials_id(ti);
        idxTest = tii;
        idxTrain = setdiff(left_trials_id, idxTest);
       
        % XTrain = mean(NeuronDataBinAll(:,bin_idx_left,idxTrain), 3)';
        % YTrain_shf = bin_idx_left(randperm(length(bin_idx_left)))'; % 打乱训练标签

        tempTrain = permute(NeuronDataBinAll(:,bin_idx_left,idxTrain), [2, 3, 1]);
        XTrain = reshape(tempTrain, [length(bin_idx_left) * length(idxTrain), size(NeuronDataBinAll, 1)]);
        YTrain_shf = repmat((bin_idx_left)', length(idxTrain), 1);
        YTrain_shf = YTrain_shf(randperm(length(YTrain_shf)));

        XTest  = NeuronDataBinAll(:,bin_idx_left,idxTest)';
        YTest = bin_idx_left';
        
        mdl = fitcecoc(XTrain, YTrain_shf);
        YHat = predict(mdl, XTest);
        
        ACCbinSf_left_trials(ti, :) = (YHat == YTest)';
        ErrBinSf_left_trials(ti, :) = (abs(YHat - YTest) * bin_size)';
    end
    % 记录这一次 shuffle 实验的平均结果
    ACCbinSf_left(sfi, :) = mean(ACCbinSf_left_trials, 1);
    ErrBinSf_left(sfi, :) = mean(ErrBinSf_left_trials, 1);
end

% 得到 Shuffle 的平均基准线 (1 x binNum)
ACCbinSf_left_Avg = mean(ACCbinSf_left, 1);
ACCSf_left_Avg = mean(ACCbinSf_left_Avg);
ErrBinSf_left_Avg = mean(ErrBinSf_left, 1);
ErrSf_left_Avg = mean(ErrBinSf_left_Avg);

%% --- Right：实际解码准确率 (Actual) ---
bin_num = length(bin_idx_right);
ACCbin_right_trials = zeros(right_trials_num, bin_num);
ErrBin_right_trials = zeros(right_trials_num, bin_num);

for ti = 1:right_trials_num
    % 1. 划分索引
    tii = right_trials_id(ti);
    idxTest = tii;
    idxTrain = setdiff(right_trials_id, idxTest);

    % 2. 准备训练数据：
    % 结果维度：binNum x n (特征)
    % XTrain = mean(NeuronDataBinAll(:,bin_idx_right,idxTrain), 3)';
    % YTrain = bin_idx_right';

    % 原始维度: [nNeuron, binNum, trainTrialNum]
    % 转置并重排为: [binNum, trainTrialNum, nNeuron]
    tempTrain = permute(NeuronDataBinAll(:,bin_idx_right,idxTrain), [2, 3, 1]);
    % 展平为: [(binNum * trainTrialNum), nNeuron]
    XTrain = reshape(tempTrain, [length(bin_idx_right) * length(idxTrain), size(NeuronDataBinAll, 1)]);
    YTrain = repmat((bin_idx_right)', length(idxTrain), 1);

    
    % 3. 准备测试数据：即当前的这个 trial
    % 结果维度：binNum x n 
    XTest = NeuronDataBinAll(:,bin_idx_right,idxTest)'; 
    YTest = bin_idx_right';
    
    % 4. 训练与预测
    mdl = fitcecoc(XTrain, YTrain);
    YHat = predict(mdl, XTest);
    
    % 5. 记录当前 trial 的准确率 (每个 bin 预测对了吗？)
    ACCbin_right_trials(ti, :) = (YHat == YTest)'; 
    
    % 6. 记录解码误差 (注意right arm编号继于left arm之后，在计算时候需要考虑减去)
    mask = YHat > binNum(1);
    adjustedYHat = YHat;
    adjustedYHat(mask) = adjustedYHat(mask) - binNum(2);
    bin_dist = abs(adjustedYHat - YTest);
    actual_dist = bin_dist * bin_size;
    ErrBin_right_trials(ti, :) = actual_dist';
end

% 得到每个 Bin 的平均准确率 (1 x binNum)
ACCbin_right_Avg = mean(ACCbin_right_trials, 1);
ACC_right_Avg = mean(ACCbin_right_Avg);

ErrBin_right_Avg = mean(ErrBin_right_trials, 1);
Err_right_Avg = mean(ErrBin_right_Avg);     

% --- Right：计算 Shuffle 解码准确率 ---
% rng("shuffle"); 
ACCbinSf_right = zeros(sft, bin_num);
ErrBinSf_right = zeros(sft, bin_num);

for sfi = 1:sft
    fprintf('sfi = %d\n', sfi);
    ACCbinSf_right_trials = zeros(right_trials_num, bin_num);
    ErrBinSf_right_trials = zeros(right_trials_num, bin_num);
    for ti = 1:right_trials_num
        tii = right_trials_id(ti);
        idxTest = tii;
        idxTrain = setdiff(right_trials_id, idxTest);
       
        % XTrain = mean(NeuronDataBinAll(:,bin_idx_right,idxTrain), 3)';
        % YTrain_shf = bin_idx_right(randperm(length(bin_idx_right)))'; % 打乱训练标签

        tempTrain = permute(NeuronDataBinAll(:,bin_idx_right,idxTrain), [2, 3, 1]);
        XTrain = reshape(tempTrain, [length(bin_idx_right) * length(idxTrain), size(NeuronDataBinAll, 1)]);
        YTrain_shf = repmat((bin_idx_right)', length(idxTrain), 1);
        YTrain_shf = YTrain_shf(randperm(length(YTrain_shf)));

        XTest  = NeuronDataBinAll(:,bin_idx_right,idxTest)';
        YTest = bin_idx_right';
        
        mdl = fitcecoc(XTrain, YTrain_shf);
        YHat = predict(mdl, XTest);
        
        ACCbinSf_right_trials(ti, :) = (YHat == YTest)';
        % (注意right arm编号继于left arm之后，在计算时候需要考虑减去)
        tempYHat = YHat;
        mask = YHat > binNum(1);
        tempYHat(mask) = tempYHat(mask) - binNum(2);
        ErrBinSf_right_trials(ti, :) = (abs(tempYHat - YTest) * bin_size)';
    end
    % 记录这一次 shuffle 实验的平均结果
    ACCbinSf_right(sfi, :) = mean(ACCbinSf_right_trials, 1);
    ErrBinSf_right(sfi, :) = mean(ErrBinSf_right_trials, 1);
end

% 得到 Shuffle 的平均基准线 (1 x binNum)
ACCbinSf_right_Avg = mean(ACCbinSf_right, 1);
ACCSf_right_Avg = mean(ACCbinSf_right_Avg);
ErrBinSf_right_Avg = mean(ErrBinSf_right, 1);
ErrSf_right_Avg = mean(ErrBinSf_right_Avg);

% 计算平均值
ACC_Avg = (ACC_left_Avg + ACC_right_Avg)/2;
ACCSf_Avg = (ACCSf_left_Avg + ACCSf_right_Avg)/2;
ACCbin_Avg = (ACCbin_left_Avg + ACCbin_right_Avg) ./2;
ACCbinSf_Avg = (ACCbinSf_left_Avg + ACCbinSf_right_Avg) ./2;

Err_Avg = (Err_left_Avg + Err_right_Avg)/2;
ErrSf_Avg = (ErrSf_left_Avg + ErrSf_right_Avg)/2;
ErrBin_Avg = (ErrBin_left_Avg + ErrBin_right_Avg) ./2;
ErrBinSf_Avg = (ErrBinSf_left_Avg + ErrBinSf_right_Avg) ./2;

%% 指定要保存的变量名
varsToSave = {
    'ACC_Avg', 'ACCSf_Avg', ...
    'ACCbin_Avg', 'ACCbinSf_Avg',...
    'Err_Avg', 'ErrSf_Avg',...
    'ErrBin_Avg', 'ErrBinSf_Avg'
    };

saveFileName = 'DecAcc_BinIdx_Tmaze.mat';
finalSavePath = fullfile(resultsFolderPath, saveFileName);

save(finalSavePath, varsToSave{:});
fprintf('数据已成功保存至: %s\n', finalSavePath);

%% 绘制解码正确率曲线
[~, bin_num] = size(ACCbin_Avg);

x = 1:bin_num;

figure('Color', 'w');
hold on;

patch([4.75 5.25 4.25 5.75], [0 0 100 100], [1 0 0], ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% --- 绘制均值曲线 ---
plot(x, ACCbin_Avg, 'b', 'LineWidth', 2, 'DisplayName', 'Actual');
plot(x, ACCbinSf_Avg, 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'DisplayName', 'Shuffle');

% --- 图形修饰 ---
ylim([0 1]);
xlim([1 11]);
set(gca, 'XTick', 1:10);
set(gca, 'YTick', 0:0.25:1);

xlabel('Spatial bin');
ylabel('Decoding accuracy');

set(gca, 'TickDir', 'in');
grid off;                  
box off;                  

legend('Location', 'southeast');

hold off;

saveName = 'DecAcc_BinIdx_Curve.png';
savePath = fullfile(resultsFolderPath, saveName);

if ~exist(resultsFolderPath, 'dir')
    mkdir(resultsFolderPath);
end

exportgraphics(gcf, savePath, 'Resolution', 300);
fprintf('解码准确率曲线图已保存至: %s\n', savePath);

%% 绘制解码正确率散点图
acc = ACCbin_Avg(:);
accSf = ACCbinSf_Avg(:);
binNum_plot = length(acc);

mu = [mean(acc), mean(accSf)]; 
sem = [std(acc)/sqrt(binNum_plot), std(accSf)/sqrt(binNum_plot)]; 

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
x_jit1 = 1 + (rand(binNum_plot, 1) - 0.5) * 0.15; % 均匀分布抖动
x_jit2 = 2 + (rand(binNum_plot, 1) - 0.5) * 0.15;

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

saveName = 'DecAcc_BinIdx_Tmaze_Plot.png';
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

ylabel('Decoding Error (cm)'); % 纵坐标改为误差单位
title('Decoding Distance Error', 'FontSize', 12);
grid off;
hold off;

% 保存图像
errPlotSaveName = 'DecError_BinIdx_CF_Plot.png';
errPlotSavePath = fullfile(resultsFolderPath, errPlotSaveName);
exportgraphics(gcf, errPlotSavePath, 'Resolution', 300);
fprintf('解码误差散点图已保存至: %s\n', errPlotSavePath);
