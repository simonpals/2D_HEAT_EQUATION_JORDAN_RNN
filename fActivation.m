function [ out_fx ] = fActivation(in_x,type)

a = 1.7159;
b = 2/3;

if strcmp(type,'none')
    out_fx = in_x;
else if strcmp(type,'tanh') 
    out_fx = a*tanh(b*in_x);
elseif strcmp(type, 'sigmoid')
    out_fx = 1 / (1 + exp(-in_x));
end

end

