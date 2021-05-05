function [ error ] = CostFunction2( Yhat_T, deltaX, deltaY, deltaT, thermal_diffusivity_factor, boundary_time )

[NYnodes, NXnodes, NTnodes] = size(Yhat_T);
error = 0;
rfx = 0.5*(thermal_diffusivity_factor*deltaT)/(deltaX^2);
rfy = 0.5*(thermal_diffusivity_factor*deltaT)/(deltaY^2);

if -1 == boundary_time
    boundary_time = NTnodes;
end

for y = 2:1:NYnodes-1
    for x = 2:1:NXnodes-1
        for t = 2:1:boundary_time % !!!
%         for t = boundary_time:1:boundary_time
            error = error + 0.5*( (Yhat_T(y,x,t)-Yhat_T(y,x,t-1)) ...
            - rfx*(Yhat_T(y,x+1,t)-2*Yhat_T(y,x,t)+Yhat_T(y,x-1,t)+Yhat_T(y,x-1,t-1)+Yhat_T(y,x+1,t-1)-2*Yhat_T(y,x,t-1)) ...
            - rfy*(Yhat_T(y+1,x,t)-2*Yhat_T(y,x,t)+Yhat_T(y-1,x,t)+Yhat_T(y-1,x,t-1)+Yhat_T(y+1,x,t-1)-2*Yhat_T(y,x,t-1)) )^2;            
        end
    end
end

end

