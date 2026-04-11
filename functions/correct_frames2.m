function [fmsout] = correct_frames2(tsdata,CaCamID,BehavCamID)
% CORRECT_FRAMES
% Have the same function as correct_frames1 in Tmaze (Previous)
% Rename this function as correct_frames2 to prevente citation error cause
% from the duplication from CellX_analysis

data1 = tsdata.data;
ts0 = data1(data1(:,1)==CaCamID,:);
ts1 = data1(data1(:,1)==BehavCamID,:);

ts0(1,3) = 1;
ts1(1,3) = 1;

% Find the repeat frames in CaCam and make them unique
% But actually this condition is infrequent
% It is more like an insurance
t = ts0(:,3);
while length(unique(t)) ~= length(t)
    [~,ic,~] = unique(t);
    t(~ismember((1:length(t)),ic),:) = t(~ismember((1:length(t)),ic),:)+1;
end

fms = ts0(:,2);
timeq = ts1(:,3);
fmsout = interp1(t,fms,timeq,'nearest', 'extrap');
fmsout = round(fmsout);

end

