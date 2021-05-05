function [ out_fx_derr ] = fActivation_derivative(in_fx,type)

a = 1.7159;
b = 2/3;

if strcmp(type,'none')
    out_fx_derr = 1.0;
elseif strcmp(type,'tanh')
    out_fx_derr = (b/a)*(a - in_fx)*(a + in_fx);
elseif strcmp(type, 'sigmoid')
    out_fx_derr = in_fx*(1-in_fx);
end

end

