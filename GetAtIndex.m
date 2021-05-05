function [ res ] = GetAtIndex( array, index1, index2, index3, index4 )

size_arr = size(array);
res = 0;
assert(nargin-1 == length(size(array)));

switch nargin
    case 2
        ndim1 = size_arr(1);
        if index1 > 0 && index1 <= ndim1
            res = array(index1);
        end
    case 3
        ndim1 = size_arr(1);
        ndim2 = size_arr(2);
        if index1 > 0 && index1 <= ndim1...
            && index2 > 0 && index2 <= ndim2
            res = array(index1, index2);
        end
    case 4
        ndim1 = size_arr(1);
        ndim2 = size_arr(2);
        ndim3 = size_arr(3);
        if index1 > 0 && index1 <= ndim1...
            && index2 > 0 && index2 <= ndim2...
            && index3 > 0 && index3 <= ndim3
            res = array(index1, index2, index3);
        end
    case 5
        ndim1 = size_arr(1);
        ndim2 = size_arr(2);
        ndim3 = size_arr(3);
        ndim4 = size_arr(4);
        if index1 > 0 && index1 <= ndim1...
            && index2 > 0 && index2 <= ndim2...
            && index3 > 0 && index3 <= ndim3...
            && index4 > 0 && index4 <= ndim4
            res = array(index1, index2, index3, index4);
        end        
    otherwise
        assert(0);
end

end

