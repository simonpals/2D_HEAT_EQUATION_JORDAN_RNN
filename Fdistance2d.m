function [ Max_Deviation, MLS ] = Fdistance2d( U1, U2 )

assert(ndims(U1) == ndims(U2));
assert(ndims(U2) == 2);

dim_size = size(U1);
max_dev = 0;
% max_dev_t = 0;
max_dev_x = 0;
max_dev_y = 0;
MLS = 0;

% for t = 1:1:dim_size(3)
    for y = 1:1:dim_size(1)
        for x = 1:1:dim_size(2)
            
            MLS = MLS + (U1(y,x) - U2(y,x))^2;
            
            if max_dev < abs(U1(y,x) - U2(y,x))
                max_dev = abs(U1(y,x) - U2(y,x));
%                 max_dev_t = t;
                max_dev_x = x;
                max_dev_y = y;
            end
        end
    end
% end

Max_Deviation = [max_dev_x max_dev_y max_dev];

end

