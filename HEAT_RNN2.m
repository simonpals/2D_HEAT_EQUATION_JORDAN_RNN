function [ Yhat_T, XT ] = HEAT_RNN2( NTnodes, NXnodes, NYnodes, Nhiden, Uinit, W, Wrec, activation_type, biases, boundary_time, test_number )
%                                 ( NTnodes, N, m, Uinit, W, Wrec, activation_type, isLearningMode )

    Yhat_T = repmat(0, [NYnodes NXnodes NTnodes]);
    XT = repmat(0, [NYnodes Nhiden+NXnodes NTnodes]);
    x_t = zeros(NYnodes, Nhiden+NXnodes);
    Yhat = zeros(NYnodes,NXnodes);
    
    Yhat(:,:) = Uinit(:,:);
    x_t(:,Nhiden+1:Nhiden+NXnodes) = Uinit(:,:);
    xt_1 = x_t;
    
    if -1 == boundary_time
        boundary_time = NTnodes;
    end

    XT(:,:,1) = x_t(:,:);
    Yhat_T(:,:,1) = Yhat(:,:);
    
    for t = 2:1:boundary_time    
        [x_t,Yhat] = NET2d(NXnodes,NYnodes,Nhiden,W,Wrec,xt_1,activation_type, biases); 
        if test_number == 5
            Yhat(1,:) = Yhat(2,:);
            Yhat(NYnodes,:) = Yhat(NYnodes-1,:);
        elseif test_number == 7    
            Yhat(1,:) = 0;
            Yhat(NYnodes,:) = Yhat(NYnodes-1,:);
        else
            Yhat(1,:) = 0;
            Yhat(NYnodes,:) = 0;
        end
        Yhat(:,1) = 0;
        Yhat(:,NXnodes) = 0;
        x_t(:,Nhiden+1:Nhiden+NXnodes) = Yhat(:,:);
        xt_1 = x_t;
        
        XT(:,:,t) = x_t(:,:);
        Yhat_T(:,:,t) = Yhat(:,:);
    end %t

end

