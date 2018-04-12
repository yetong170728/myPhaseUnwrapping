function [ output_args ] = branchCut( resMap ,MaxBoxRadius )
%Form branch cuts from residue map
%   Last modified by Hanyu@cbir ,4/11/2018
%   Input: resMap 
%               
S = size(resMap); 
% assume 2 dimensional case

% flag the active residue 
resActMap = zeros(S,'logical');

% if the residue is balanced
isBalance = zeros(S,'logical');

% charge to balance the residue
charge = zeros(S);

% is border pixel
isBorder = zeros(S,'logical');
isBorder(:,1) = true;
isBorder(1,:) = true;
isBorder(:,end) = true;
isBorder(end,:) = true;

for i = 1:S(1)
    for j = 1:S(2)
        if (resMap(i,j) ~= 0) && (~isBalance(i,j))
%             isResidue + unbalanced

            resActMap(i,j) = true;
%             mark the residue active


            charge(i,j) = 1*(resMap(i,j)>0) - 1*(resMap(i,j)<0);
%             let charge = 1/-1 for positive/negative residue
            
            for r = 1:MaxBoxRadius
%                 r for boxradius, box size n = 1+2r
                
               
            end
                
        end
    end
end



end

