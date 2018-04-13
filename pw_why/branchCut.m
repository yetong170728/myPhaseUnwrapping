function [ brCut ] = branchCut( resMap ,MaxBoxRadius )
%Form branch cuts from residue map
%   Last modified by Hanyu@cbir ,4/11/2018
%   Input: 
%            resMap, a map for polarities of residues(+1, -1 and 0)
%            MaxBoxRadius, maximum of searching radius
%   Output:
%            brCut, branch-cut map
% modified by Hanyu@cbir(c), 4/13/2018
%% Below are some pre-processing work
% if nargin < 3
% %     if  the pixel is a border pixel
% %     can be assigned outside the function
%     isBorder = zeros(S,'logical');
%     isBorder(:,1) = true;
%     isBorder(1,:) = true;
%     isBorder(:,end) = true;
%     isBorder(end,:) = true;
% end

% update the resMap bordered by isBorder

% get the indices of residue in resMap, resRow has same length with resCol
[resRow, resCol]  = find(resMap);

% assume 2 dimensional case
S = size(resMap); 

% intialize the branch-cut matrix
brCut = zeros(S, 'logical');

% if the residue is active
% isResAct = zeros(S,'logical');

% if the residue is balanced
isBalance = zeros(S,'logical');

% if the residue is isolated
isIso = zeros(size(resRow), 'logical');

% charge to balance the residue
charge = zeros( size(resRow) );

%% Belows are some persudo code for setting branch cuts
%             what to do? -- for each UNBALANCED residue
%                 1. get the outermost pixels (named as P) around the active residue
%                 2. get the UNBALANCED residues in P, stored as Pres
%                 3. get the border pixel in P, stored as Pborder
%                 4. while( charge(i) )
%                     4.1 if charge(i) is non-zero (this is not necessary)
%                         5. first for each UNBALANCED residue in Pres (except the active residue itself)
%                             pick some of the residues to balance the charge(i)
%                             mark THESE residues as BALANCED
%                             set charge(i) = 0
%                             mark the residue BALANCED
%                             set a minimum-length branch cut connecting all the residues
%                             break;
%                         6. if all residues in Pres don't work( charge(i) can not reach zero regardless of nearby residues ):
%                             the question now is how to choose: Expand radius or set branch cut to border ?
%                             if the Pborder is not empty, try border
%                                 pick the nearest border pixel(it should be the first one in Pborder) to the active residue
%                                 set a branch cut
%                                 set charge(i) = 0               
%                             else expand radius
%                                 update the Pres and Pborder
%                     4.2 if charge of the residue( refers to charge(i) ) is zero, the balance is finished, mark the residue BALANCED
%   


%% Below are the real part of this function

% ----------------------------------------------------------------------------------------------------------------
% This is a manner to realize the working flow, very easy-to-think
% when we find border first, we don't continuously find residues to balance
% the active residue, and that's a necessity in Goldstrein's Branch Cut
% Method
for i =1:length(resRow)
%     for all residues in resMap
    resX = resRow(i); resY=resCol(i);
    if ~isBalance(resX,resY)
%     for each unbalanced residue
%         isResAct(resX,resY) = true;
        charge(i) = resMap(resX, resY);
%         mark the residue active and give the intial charge
        while( charge(i) )
%            Loop while the charge is NOT zero
            for radius = 1:MaxBoxRadius
%               for n = 3 to MaxBoxSize} step 2
                for pX = resX-radius : resX+radius
                    for pY = resY-radius: resY+radius
                        tt_res = zeros(size(resMap)); tt_res(resX,resY) = resMap(resX,resY); % 1 for central residue
                        tt_res(pX,pY) = 2; % 2 for box pixel
                        if charge(i) && ( (abs(pX-resX) == radius) || (abs(pY-resY) == radius) ) 
%                         for the outermost pixels of search box centered at
%                         (resX, resY), and  charge(i) now is non-zero
                            [pX1,pY1,borderFlag] = getBorder(S, [pX,pY]);
                            if borderFlag
%                                 if the box pixel reaches border of resMap
                                charge(i) = 0;
                                isBalance(resX,resY) = true;
                                brCut = setBranchCut(brCut, [resX, resY], [pX1,pY1]);
                            elseif resMap(pX1,pY1)
%                             elseif resMap(pX1,pY1)  && ~isResAct(pX1,pY1)
%                                 if the box pixel is an inactive residue
                                brCut = setBranchCut(brCut, [resX, resY], [pX1,pY1]);
%                                 connect all residues to the center residue under whatever conditions
                                if ~isBalance(pX1,pY1)
                                    charge(i) = charge(i) + resMap(pX1,pY1);
%                                     isResAct(pX1,pY1) = true;
                                    isBalance(pX1,pY1) = true;             
                                end
                                if ~charge(i)
%                             if the charge of residue is 0
                                    isBalance(resX,resY) = true;
                                end 
                            end % end of the borderflag
                        end % end of outermost shell
                    end % end of pY loop
                end  % end of pX loop            
                if ~charge(i)
                    break;
                end
            end % end of radius loop
            if charge(i)
%             if not break the loop, so the charge(i) can not be zero even MaxBoxRadius is reached
                isIso(i) = true; % set this residue as an isolated residue
%             get the nearest border
                bX = abs(resX - S(1)/2); bY = abs(resY-S(2)/2);
                if bX <= bY
                    bX = resX;
                    bY = 1*(resY < S(2)/2) + S(2)*(resY >= S(2)/2);
                else
                    bY = resY;
                    bX = 1*(resX < S(1)/2) + S(1)*(resX >= S(1)/2);
                end
                charge(i) = 0;
                brCut = setBranchCut(brCut, [resX, resY], [bX,bY]);
                isBalance(resX,resY) = true;
            end % end of ~charge(i) loop   
         end   
    end % end of judge if the residue is BALANCE
end % end of all residues
% ----------------------------------------------------------------------------------------------------------------

end % End of the function


function [brIn] = setBranchCut(brIn, pS, pE)
% this function is to set the branch cut in brCut from start point to end point
% brOut = brIn;
% if sum(pS ~= pE)<1 % if same points are put in
%     return;
% end
% dX = pE(1) - pS(1);
% dY = pE(2) - pS(2);
% r = sqrt(dX^2+dY^2);
% orient = [dX, dY] ./ r;
% steps = floor(r); % we need #steps pixels to draw a branch cut
% 
% for r = 1:steps
%     brOut( pS(1)+round(r*orient(1)), pS(2)+round(r*orient(2)) ) = 1;
% end

r1 = pS(1);c1= pS(2);
r2 = pE(1);c2 = pE(2);

brIn(r1,c1)=1;                       %Fill the initial points
brIn(r2,c2)=1;
radius=sqrt((r2-r1)^2 + (c2-c1)^2);         %Distance between points
warning off MATLAB:divideByZero;            
m=(c2-c1)/(r2-r1);                          %Line gradient
theta=atan(m);                              %Line angle

if c1==c2                                   %Cater for the infinite gradient
    if r2>r1
        for i=1:radius
            r_fill=r1 + i;  
            brIn(r_fill, c1)=1;
        end
        return;
    else
        for i=1:radius
            r_fill=r2 + i;   
            brIn(r_fill, c1)=1;
        end
        return;
    end
end

%Use different plot functions for the different quadrants (This is very clumsy, 
%I'm sure there is a better way of doing it...)
if theta<0 && c2<c1                           
    for i=1:radius                          %Number of points to fill in
        r_fill=r1 + round(i*cos(theta));    
        c_fill=c1 + round(i*sin(theta));
        brIn(r_fill, c_fill)=1;
    end
    return;
elseif theta<0 && c2>c1 
    for i=1:radius                         %Number of points to fill in
        r_fill=r2 + round(i*cos(theta));
        c_fill=c2 + round(i*sin(theta));
        brIn(r_fill, c_fill)=1;
     end
    return;
end

if theta>0 && c2>c1                          
    for i=1:radius                          %Number of points to fill in
        r_fill=r1 + round(i*cos(theta));    
        c_fill=c1 + round(i*sin(theta));
        brIn(r_fill, c_fill)=1;
    end
    return;
elseif theta>0 && c2<c1 
     for i=1:radius                         %Number of points to fill in
        r_fill=r2 + round(i*cos(theta));
        c_fill=c2 + round(i*sin(theta));
        brIn(r_fill, c_fill)=1;
    end
    return;   
end
end

function [pX,pY,flag] = getBorder(SizeOfMap, index)
% set the correct index (pX,pY) if meets border
% flag =1 for the point defined by index meets border
% default for 2-D situation
flagX = true;
flagY = true;

if index(1) <= 1
    pX = 1;  
elseif index(1) >= SizeOfMap(1)
    pX = SizeOfMap(1);    
else
    pX = index(1);
    flagX = false;
end

if index(2) <= 1
    pY = 1;
elseif index(2) >= SizeOfMap(2)
    pY = SizeOfMap(2);
else
    pY = index(2);
    flagY = false;
end

flag = flagX || flagY;
end