function [ resMap ] = phaseResidue( phaseData, flag)
% function [ resMap, resMap2 ] = phaseResidue( phaseData, flag)
%calculate phase residues in phase unwrapping
% Input:
%    phaseData: 
%    flag: 0 for 2D or four-point residue, 1 for 3D residue, 2 for 4D residues
% Ref. to Book, "Two dimensional phase unwrapping: Theory, algorithm and software"
% Last modified by Hanyu@cbir(c), 4/11/2018 

if nargin < 2
    flag = 0;
end

switch flag
    case 0
        [nx,ny] = size(phaseData);

        
%         An awkard serial procedure
%         resMap = zeros(size(phaseData));
%         for j =1:ny-1
%             for i = 1:nx-1
% %                 four-point residue
%                 resMap(i,j) = wrap( phaseData(i+1,j) - phaseData(i,j) ) + ...
%                                      wrap( phaseData(i+1,j+1) - phaseData(i+1,j) ) + ...
%                                      wrap( phaseData(i,j+1) - phaseData(i+1, j+1) ) + ...
%                                      wrap( phaseData(i,j) - phaseData(i,j+1) );    
%                 
%             end
%         end

%         Another matrix procedure, more cost on memory but fast
        phaseUp = phaseData(1:end-1,1:end-1);
        phaseBelow = phaseData(2:end,1:end-1);
        phaseRight = phaseData(1:end-1,2:end);
        phaseBelowRight = phaseData(2:end,2:end);
        
        d1 = wrap(phaseBelow - phaseUp);
        d2 = wrap(phaseBelowRight - phaseBelow);
        d3 = wrap(phaseRight - phaseBelowRight );
        d4 = wrap(phaseUp - phaseRight);
        res = d1+d2+d3+d4;
%         resMap(resMap > 6) = 1; % residue for 2pi
%         resMap(resMap < -6) = -1; % residue for  -2pi
%         If finish by now, some tiny error number(machine accuracy) will
%         disturb the resMap, so sth is needed for a pure bit-based resMap
%         thus,
        resMap = (res > 6) - (res < -6);
        
%         Now append zeros to make resMap2  a nx*ny size matrix
        resMap = [resMap, zeros(nx-1,1)];
        resMap = [resMap; zeros(1,ny)];
        
    case 1
        return;
    case 2
        return;
    otherwise
        error('Input arg ''flag'' must be 0(2D), 1(3D), or 2(4D)!');
end       
end