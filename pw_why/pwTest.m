clear,clc,close all
phase = -pi + 2*pi*rand(16,16);
% load 'IM.mat'; IM = IM(1:10,1:10); phase = angle(IM);

N = size(phase);
figure, imshow(phase, [-pi, pi]) ,title('Original phase');
[ resMap ] = phaseResidue(phase);
% resMap3 = PhaseResidues(phase, ones(size(phase))); % 如有很小的差别, 应该是机器精度误差, 由于一定是2pi整数倍,所以不影响
figure, imshow(resMap) ,title('Residue flag');


%% a test for the branch cut generation
close all
resMap = zeros(10);
resMap(5,5) = -1;
resMap(6,7) = 1;
resMap(4,3) = -1;
resMap(4,5) = 1;

IM_mask = ones(size(resMap)); maxBoxRadius = 2;
brCut=BranchCuts(resMap, maxBoxRadius, IM_mask);            %Place branch cuts
brCut2 = branchCut( resMap ,maxBoxRadius );

% figure, imshow([brCut, ones(10,10), brCut2]);
% imtool(brCut);
% imtool(brCut2);
