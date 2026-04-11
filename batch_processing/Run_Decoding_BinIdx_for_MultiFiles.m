%% Run_Decoding_BinIdx_for_MultiFiles

% 脚本描述：批量运行 Decoding_BinIdx 脚本，并保存DecAcc_BinIdx*文件至有序结构

% 作者：WXW

% 创建时间：2026/02/17

% 更新时间：2026/02/19

%%
close all
clear
clc

%% 初始化路径
sourceBase = 'D:\LAB\Miniscope\Basic_Simplified'; % 请确保这是你的实际根目录
targetBase = 'D:\LAB\Miniscope\Decoding_BinIdx_Results_260219';

% 定义配置：将名字作为 key，方便精准匹配
% config.CircleField = struct('script', 'Decoding_BinIdx_CF', 'genDir', 'DecAcc_BinIdx_CF');
% config.OpenField   = struct('script', 'Decoding_BinIdx_OF', 'genDir', 'DecAcc_BinIdx_OF');
config.Tmaze       = struct('script', 'Decoding_BinIdx_Tmaze', 'genDir', 'DecAcc_BinIdx_Tmaze');

expNames = fieldnames(config);

%% 开始遍历
for i = 1:numel(expNames)
    thisExp = expNames{i};
    expPath = fullfile(sourceBase, thisExp);
    
    if ~exist(expPath, 'dir'), continue; end
    
    % 获取当前实验类型的配置
    currentCfg = config.(thisExp);
    
    % 遍历基因型 (KO, WT)
    genotypes = {'KO', 'WT'};
    for g = 1:length(genotypes)
        genoName = genotypes{g};
        genoPath = fullfile(expPath, genoName);
        if ~exist(genoPath, 'dir'), continue; end
        
        % 获取小鼠文件夹
        mouseDirs = dir(genoPath);
        for m = 1:length(mouseDirs)
            if strcmp(mouseDirs(m).name, '.') || strcmp(mouseDirs(m).name, '..') || ~mouseDirs(m).isdir
                continue; 
            end
            
            mouseName = mouseDirs(m).name;
            mousePath = fullfile(genoPath, mouseName);
            
            % 获取每天数据文件夹
            dayDirs = dir(mousePath);
            for d = 1:length(dayDirs)
                if strcmp(dayDirs(d).name, '.') || strcmp(dayDirs(d).name, '..') || ~dayDirs(d).isdir
                    continue;
                end
                
                dayName = dayDirs(d).name;
                currentDataPath = fullfile(mousePath, dayName);
                
                % --- 执行处理 ---
                % 记录初始目录，确保报错后能切回来
                originalDir = pwd; 
                
                try
                    cd(currentDataPath);
                    fprintf('正在处理 [%s]: %s\n', thisExp, dayName);
                    
                    % 运行处理脚本
                    feval(currentCfg.script); 
                    
                    % --- 移动生成的文件夹 ---
                    % 1. 确定源文件夹位置
                    generatedFolderPath = fullfile(currentDataPath, currentCfg.genDir);
                    
                    % 2. 确定目标父目录：D:\...\Pick_up_Place_Cells\OpenField\KO\Mouse\Day
                    finalDestParent = fullfile(targetBase, thisExp, genoName, mouseName, dayName);
                    
                    if exist(generatedFolderPath, 'dir')
                        % 3. 强制递归创建目标父路径 (确保 mkdir 不会因为层级深而失败)
                        if ~exist(finalDestParent, 'dir')
                            [status, msg] = mkdir(finalDestParent);
                            if ~status, error('创建目录失败: %s', msg); end
                        end
                        
                        % 4. 执行移动
                        movefile(generatedFolderPath, finalDestParent, 'f'); 
                        fprintf('  [成功] 已移动至目标目录\n');
                    else
                        warning('  [警告] 脚本已运行但未发现文件夹: %s', currentCfg.genDir);
                    end
                    
                catch ME
                    fprintf('  [错误] 处理 %s 时出错: %s\n', dayName, ME.message);
                end
                
                % 回到初始目录，防止路径错乱
                cd(originalDir);
                close all
            end
        end
    end
end
fprintf('### 批量任务结束 ###\n');
