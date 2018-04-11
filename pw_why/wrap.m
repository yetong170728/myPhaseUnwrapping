function [ phaseOut ] = wrap( phaseIn )
%wrap-operator in phase unwrapping
% phaseOut = mod(phaseIn + pi, 2pi) - pi
% Ref. to Book, "Two dimensional phase unwrapping: Theory, algorithm and software"
% Last modified by Hanyu@cbir(c), 4/11/2018 
phaseOut = mod(phaseIn+pi, 2*pi) - pi;
end

