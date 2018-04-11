function [curl ] = curl4(data)
% test for cal of curl
[nx,ny] = size(data);
if nx*ny > 4
    error('Wrong input!');
end
d1 = wrap(data(2,1) - data(1,1) );
d2 = wrap(data(2,2) - data(2,1) );
d3 = wrap(data(1,2) - data(2,2) );
d4 = wrap(data(1,1) - data(1,2) );

curl = d1+d2+d3+d4;

end

