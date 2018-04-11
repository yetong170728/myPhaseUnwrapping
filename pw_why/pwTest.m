clear,clc,close all
phase = -pi + 2*pi*rand(128,128);
% load 'IM.mat'; IM = IM(1:10,1:10); phase = angle(IM);

N = size(phase);
figure, imshow(phase, [-pi, pi]) ,title('Original phase');
[ resMap ] = phaseResidue(phase);
% resMap3 = PhaseResidues(phase, ones(size(phase))); % 如有很小的差别, 应该是机器精度误差, 由于一定是2pi整数倍,所以不影响
figure, imshow(resMap) ,title('Residue flag');

y=zeros(N,'logical');