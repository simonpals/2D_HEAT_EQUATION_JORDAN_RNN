function [ Ures, W, Wrec ] = Test_RNN2( W, Wrec, ExactSol, task_number, optimizing_method )
 
% look to CostFunction2
if task_number > -1
    test_number = task_number;
else
    test_number = 6;
end    
use_new_weights = length(W)<=1 || length(Wrec)<=1; %true;
if optimizing_method > -1
    optimization_method = optimizing_method;
else
    optimization_method = 5;
end
maximum_passes = 2000;


domain_len = 1;
% if test_number <= 3
%     domain_len = 2;
% else
%     domain_len = 1;
% end
learning_rate = 0.001; % 5^0.5; % 0.01125; % 
T = 0.05;
deltaT = 0.005;
NTnodes = T / deltaT + 1;
Mx = 1 * domain_len;
deltaX = 0.05 * domain_len;
NXnodes = Mx / deltaX + 1;
My = 1 * domain_len;
deltaY = 0.05 * domain_len;
NYnodes = My / deltaY + 1;
Nhiden = NXnodes; %(NXnodes-1) * 1.5 + 1;
activation_type = 'none';

Uleft = repmat(0.0, [1 NYnodes]); % zeros(1, NYnodes);
Uright = repmat(0.0, [1 NYnodes]); 
Utop = repmat(0.0, [1 NXnodes]);
Ubottom = repmat(0.0, [1 NXnodes]);
Uinit = repmat(0.0, [NYnodes NXnodes]); %zeros(NXnodes,NYnodes);

biases = repmat(0.00, [NYnodes NXnodes]);

% thermal_diffusivity_factor = 0.2; %sqrt( 2*deltaX^2 / (3*deltaT) ) + 0.1;

switch(test_number) 
    case 0 
        %==============Test0================================================
        Uinit(ceil(NYnodes*0.2):ceil(NYnodes*0.8),ceil(NXnodes*0.2):ceil(NXnodes*0.8)) = 0.95;
    case 1
        %==============Test1================================================
        for i = 1:1:NYnodes
            for j = 1:1:NXnodes
                Uinit(i,j) = -0.1*(1+2*pi)*sin(2*pi*(deltaY*i+deltaX*j));
            end
        end
    case 2
        %==============Test2================================================
        for i = 1:1:NYnodes
            for j = 1:1:NXnodes
                Uinit(i,j) = exp(-((deltaX*j - 0.5*Mx)^2 / (2*0.3*Mx^2) + (deltaY*i-0.5*My)^2 / (2*0.1*My^2) ));
            end
        end
    case 3
        %==============Test3================================================
        thermal_diffusivity_factor=(1.0/8.0);
        Uinit(1:ceil(NYnodes*0.5),1:NXnodes) = 0.5;
    case 4
        %==============Test4================================================
        sigma = 1.0;
        thermal_diffusivity_factor = 1.0;
        Tmax = 1.0;
        x=0:deltaX:2;
        y=0:deltaY:2;
        for i = 1:1:NYnodes
            for j = 1:1:NXnodes
                Uinit(i,j) = Tmax * exp((-(x(j)^2 + y(i)^2))/(sigma^2));
            end
        end
        Uleft(:) = Uinit(:,1);
        Uright(:) = Uinit(:,NXnodes);
        Utop(:) = Uinit(NYnodes,:);
        Ubottom(:) = Uinit(1,:);
    case 5
        %==============Test5================================================
        thermal_diffusivity_factor=(1.0/5.0); %1.0;
        for i = 1:1:NYnodes
            for j = 1:1:NXnodes
                Uinit(i,j) = (deltaY*(i-1)*deltaX*(j-1)); % x*y
            end
        end
        Uleft(:) = 0;
        Uright(:) = 0;
        Utop(:) = Uinit(NYnodes,:);
        Ubottom(:) = Uinit(1,:);
    case 6
        %==============Test6================================================
        thermal_diffusivity_factor= 1.0/2.0; %(1.0/3.0); %1.0;
        for i = 1:1:NYnodes
            for j = 1:1:NXnodes
                Uinit(i,j) = (deltaY*(i-1)*(1-deltaY*(i-1))*deltaX*(j-1)*(1-deltaX*(j-1))); % x*(1-x)*y*(1-y)
            end
        end
        %===================================================================
    case 7
        %==============Test7================================================
        thermal_diffusivity_factor=(1.0/3.0); %1.0;
        Uinit(:,:) = 0.1;
        Utop(:) = Uinit(NYnodes,:);
        %===================================================================
end    

Uinit(NYnodes,:) = Utop(:);
Uinit(1,:) = Ubottom(:);
Uinit(:,1) = Uleft(:);
Uinit(:,NXnodes) = Uright(:);

[W,Wrec,biases] = HEAT_RNN2_LEARNING(optimization_method,T,maximum_passes,learning_rate,NTnodes,deltaT,NXnodes,deltaX,NYnodes,deltaY,...
 Nhiden,Uinit,Uleft,Uright,Utop,Ubottom,thermal_diffusivity_factor,activation_type, biases, test_number, use_new_weights, W, Wrec, ExactSol);
[ Yhat_T, XT ] = HEAT_RNN2( NTnodes, NXnodes, NYnodes, Nhiden, Uinit, W, Wrec, activation_type, biases, -1, test_number );                            
Ures = Yhat_T;

end

