function [W, PW] = Optimization_method2d( optimization_method, learning_rate, W, F_W, PW, NXnodes, NYnodes, Nhiden, pass_number, F_W_1, Wdirect, deltaX, deltaY, deltaT, thermal_diffusivity_factor, Uinit, activation_type, NTnodes, biases, F_biases, test_number )
%                                       ( optimization_method, pass_number, learning_rate, W, F_W, F_W_1, PW, Wrow, Wcol )

neuron_learning_rate = learning_rate;

switch(optimization_method)
    case 1
        %Steepest descent
        for ye = 1:1:NYnodes            
            for xe = 1:1:NXnodes
                for yb = 1:1:NYnodes
                    for xb = 1:1:Nhiden
                        W(ye,xe,yb,xb) = W(ye,xe,yb,xb) - neuron_learning_rate * F_W(ye,xe,yb,xb);
                    end
                end
            end
        end
               
    case 2
        %Conjugate gradient
        k = pass_number;
        if k == 1
            PW = -F_W;
        else                                                            
            %beta = norm(F_W) / norm(F_W_1); %transpose(F_W) * F_W / (transpose(F_W_1) * F_W_1);
            %PW = -F_W + beta * PW;
            
            last_beta = inf;
            for ye = 1:1:NYnodes
                for xe = 1:1:NXnodes
                                        
                    fw = repmat(0, [NYnodes Nhiden]);
                    fw_1 = repmat(0, [NYnodes Nhiden]);
                    
                    for yb = 1:1:NYnodes
                        for xb = 1:1:Nhiden
                            fw(yb,xb) = F_W(ye,xe,yb,xb);
                            fw_1(yb,xb) = F_W_1(ye,xe,yb,xb);
                        end
                    end
                    
                    beta = (norm(fw)) / (norm(fw_1)); %(fw*fw') / (fw_1*fw_1'); %
                    
                    if beta == inf || isnan(beta)
                        beta = last_beta;
                    end
                    
                    if beta ~= inf && ~isnan(beta)
                        for yb = 1:1:NYnodes
                            for xb = 1:1:Nhiden
                                PW(ye,xe,yb,xb) = -F_W(ye,xe) + beta * PW(ye,xe,yb,xb);
                            end
                        end
                        last_beta = beta;
                    end
                                        
                end
            end
            
        end
        
%------------------------------------------------------                
            [ exact ] = Exact_solution(  );
            alphas = zeros(NYnodes,NXnodes);
%             Wrec = Wrec + alpha2*PW;
            Wrec = W;

            for ye = 1:1:NYnodes
                for xe = 1:1:NXnodes                                        
                              
                    alpha1 = 0.0005/(1.0*pass_number);  
                    local_1d_iteration_count = 21;
                    error1 = 1000000000;
                    alpha2 = alpha1; %(alpha3+alpha1)*0.5;                    
                    err_curve = 1:1:local_1d_iteration_count;                    
                    start_step = 0.005;                                      
                                        
                    xb = xe;  yb = ye;
                    for iter1d = 1:1:local_1d_iteration_count 
%                         Wrec = W;
                        
                        Wrec(ye,xe,yb,xb) = W(ye,xe,yb,xb);
                        if yb+1 <= NYnodes
                            Wrec(ye,xe,yb+1,xb) = W(ye,xe,yb+1,xb);
                        end
                        if yb-1 > 0
                            Wrec(ye,xe,yb-1,xb) = W(ye,xe,yb-1,xb);
                        end
                        if xb+1 <= NXnodes
                            Wrec(ye,xe,yb,xb+1) = W(ye,xe,yb,xb+1);
                        end
                        if xb-1 > 0
                            Wrec(ye,xe,yb,xb-1) = W(ye,xe,yb,xb-1);
                        end    


                        %                     for yb = 1:1:NYnodes
                        %                         for xb = 1:1:Nhiden
                        Wrec(ye,xe,yb,xb) = Wrec(ye,xe,yb,xb) + alpha2 * PW(ye,xe,yb,xb);
                        if yb+1 <= NYnodes
                            Wrec(ye,xe,yb+1,xb) = Wrec(ye,xe,yb+1,xb) + alpha2 * PW(ye,xe,yb+1,xb);
                        end
                        if yb-1 > 0
                            Wrec(ye,xe,yb-1,xb) = Wrec(ye,xe,yb-1,xb) + alpha2 * PW(ye,xe,yb-1,xb);
                        end
                        if xb+1 <= NXnodes
                            Wrec(ye,xe,yb,xb+1) = Wrec(ye,xe,yb,xb+1) + alpha2 * PW(ye,xe,yb,xb+1);
                        end
                        if xb-1 > 0
                            Wrec(ye,xe,yb,xb-1) = Wrec(ye,xe,yb,xb-1) + alpha2 * PW(ye,xe,yb,xb-1);
                        end                        
                        %                         end
                        %                     end

                        [Yhat_T, XT] = HEAT_RNN2( NTnodes, NXnodes, NYnodes, Nhiden, Uinit, Wdirect, Wrec, activation_type, biases, 3, test_number );
                        error = CostFunction2(Yhat_T, deltaX, deltaY, deltaT, thermal_diffusivity_factor, 3);
                        
%                         [ Max_Deviation, MLS ] = Fdistance2d( exact(:,:,11), rot90(Yhat_T(:,:,11)) );
%                         err_curve(iter1d) = MLS;
                        
                        if error < error1
                            error1 = error;
                            alpha1 = alpha2;
                        end
                        alpha2 = alpha2 + start_step;
                                                
                    end
                    
                    alphas(ye,xe) = alpha1;
                    
                     Wrec(ye,xe,yb,xb) = Wrec(ye,xe,yb,xb) + alpha1 * PW(ye,xe,yb,xb);
                        if yb+1 <= NYnodes
                            Wrec(ye,xe,yb+1,xb) = Wrec(ye,xe,yb+1,xb) + alpha1 * PW(ye,xe,yb+1,xb);
                        end
                        if yb-1 > 0
                            Wrec(ye,xe,yb-1,xb) = Wrec(ye,xe,yb-1,xb) + alpha1 * PW(ye,xe,yb-1,xb);
                        end
                        if xb+1 <= NXnodes
                            Wrec(ye,xe,yb,xb+1) = Wrec(ye,xe,yb,xb+1) + alpha1 * PW(ye,xe,yb,xb+1);
                        end
                        if xb-1 > 0
                            Wrec(ye,xe,yb,xb-1) = Wrec(ye,xe,yb,xb-1) + alpha1 * PW(ye,xe,yb,xb-1);
                        end                                            

                end
            end
            
                     
     
%------------------------------------------------------        
%         alpha2 = alpha1;        
        for ye = 1:1:NYnodes            
            for xe = 1:1:NXnodes
                for yb = 1:1:NYnodes
                    for xb = 1:1:Nhiden
                        W(ye,xe,yb,xb) = W(ye,xe,yb,xb) + alphas(ye,xe) * PW(ye,xe,yb,xb);
                    end
                end
            end
        end
        alphas
    case 3       
        %Annealing
%         for i = 1:1:Wrow   
%             neuron_learning_rate = neuron_learning_rate * 0.995;
%             for j = 1:1:Wcol
%                 W(i,j) = W(i,j) - neuron_learning_rate * F_W(i,j);
%             end
%         end    
        
        for ye = 1:1:NYnodes    
            neuron_learning_rate = neuron_learning_rate * 0.995;
            for xe = 1:1:NXnodes
                for yb = 1:1:NYnodes
                    for xb = 1:1:Nhiden
                        W(ye,xe,yb,xb) = W(ye,xe,yb,xb) - neuron_learning_rate * F_W(ye,xe,yb,xb);
                    end
                end
            end
        end
       
   case 4
      %for biases optimization
            
      [Yhat_T, XT] = HEAT_RNN2( NTnodes, NXnodes, NYnodes, Nhiden, Uinit, Wdirect, W, activation_type, biases, -1, test_number );
      error1 = CostFunction2(Yhat_T, deltaX, deltaY, deltaT, thermal_diffusivity_factor, -1);
            
      min_lr=0;
      for lr = 0.001 : 0.001 : 0.2
          Wb=biases;
                    
          for y = 1:1:NYnodes
              for x = 1:1:NXnodes
                  Wb(y,x) = Wb(y,x) - lr * F_biases(y,x);
              end
          end    
          
          [Yhat_T, XT] = HEAT_RNN2( NTnodes, NXnodes, NYnodes, Nhiden, Uinit, Wdirect, W, activation_type, Wb, -1, test_number );
          error = CostFunction2(Yhat_T, deltaX, deltaY, deltaT, thermal_diffusivity_factor, -1);
          
          if(error<error1)
              error1 = error;
              min_lr = lr;
          end
          
      end
                       
      min_lr
      for y = 1:1:NYnodes
          for x = 1:1:NXnodes
              biases(y,x) = biases(y,x) - min_lr * F_biases(y,x);
          end
      end
      
      
   case 5
        %Conjugate gradient
        k = pass_number;
        if k == 1
            PW = -F_W;
        else                                                            
            %beta = norm(F_W) / norm(F_W_1); %transpose(F_W) * F_W / (transpose(F_W_1) * F_W_1);
            %PW = -F_W + beta * PW;
            
            last_beta = inf;
            for ye = 1:1:NYnodes
                for xe = 1:1:NXnodes
                                        
                    fw = repmat(0, [NYnodes Nhiden]);
                    fw_1 = repmat(0, [NYnodes Nhiden]);
                    
                    for yb = 1:1:NYnodes
                        for xb = 1:1:Nhiden
                            fw(yb,xb) = F_W(ye,xe,yb,xb);
                            fw_1(yb,xb) = F_W_1(ye,xe,yb,xb);
                        end
                    end
                    
                    beta = (norm(fw)) / (norm(fw_1)); %(fw*fw') / (fw_1*fw_1'); %
                    
                    if beta == inf || isnan(beta)
                        beta = last_beta;
                    end
                    
                    if beta ~= inf && ~isnan(beta)
                        for yb = 1:1:NYnodes
                            for xb = 1:1:Nhiden
                                PW(ye,xe,yb,xb) = -F_W(ye,xe) + beta * PW(ye,xe,yb,xb);
                            end
                        end
                        last_beta = beta;
                    end
                                        
                end
            end
            
        end
%------------------------------------------------------        
%         local_1d_iteration_count = 21;
%         start_step = 0.005;
%         alpha1 = 0.001/(1.0*pass_number);
                
%         local_1d_iteration_count = 65;
%         start_step = 0.035;
%         alpha1 = 0.001/(1000);

        local_1d_iteration_count = 5501;
        start_step = 0.005;
        alpha1 = 0.001/(10000);
%------------------------------------------------------
        Wrec = W;
        for ye = 1:1:NYnodes
            for xe = 1:1:NXnodes
%                 for yb = 1:1:NYnodes
%                     for xb = 1:1:Nhiden
                yb = ye; xb = xe;
                        Wrec(ye,xe,yb,xb) = Wrec(ye,xe,yb,xb) + alpha1 * PW(ye,xe,yb,xb);
                        if yb+1 <= NYnodes
                            Wrec(ye,xe,yb+1,xb) = Wrec(ye,xe,yb+1,xb) + alpha1 * PW(ye,xe,yb+1,xb);
                        end
                        if yb-1 > 0
                            Wrec(ye,xe,yb-1,xb) = Wrec(ye,xe,yb-1,xb) + alpha1 * PW(ye,xe,yb-1,xb);
                        end    
                        if xb+1 <= NXnodes
                            Wrec(ye,xe,yb,xb+1) = Wrec(ye,xe,yb,xb+1) + alpha1 * PW(ye,xe,yb,xb+1);
                        end    
                        if xb-1 > 0
                            Wrec(ye,xe,yb,xb-1) = Wrec(ye,xe,yb,xb-1) + alpha1 * PW(ye,xe,yb,xb-1);
                        end    
%                     end
%                 end
            end
        end
        
        [Yhat_T, XT] = HEAT_RNN2( NTnodes, NXnodes, NYnodes, Nhiden, Uinit, Wdirect, Wrec, activation_type, biases, -1, test_number );
        error1 = CostFunction2(Yhat_T, deltaX, deltaY, deltaT, thermal_diffusivity_factor, -1);
                                     
        alpha2 = alpha1; %(alpha3+alpha1)*0.5;
        
        for iter1d = 1:1:local_1d_iteration_count  
            Wrec = W;
            for ye = 1:1:NYnodes
                for xe = 1:1:NXnodes
                    yb = ye; xb = xe;
%                     for yb = 1:1:NYnodes
%                         for xb = 1:1:Nhiden
                            Wrec(ye,xe,yb,xb) = Wrec(ye,xe,yb,xb) + alpha2 * PW(ye,xe,yb,xb);
                            if yb+1 <= NYnodes
                                Wrec(ye,xe,yb+1,xb) = Wrec(ye,xe,yb+1,xb) + alpha2 * PW(ye,xe,yb+1,xb);
                            end
                            if yb-1 > 0
                                Wrec(ye,xe,yb-1,xb) = Wrec(ye,xe,yb-1,xb) + alpha2 * PW(ye,xe,yb-1,xb);
                            end
                            if xb+1 <= NXnodes
                                Wrec(ye,xe,yb,xb+1) = Wrec(ye,xe,yb,xb+1) + alpha2 * PW(ye,xe,yb,xb+1);
                            end    
                            if xb-1 > 0
                                Wrec(ye,xe,yb,xb-1) = Wrec(ye,xe,yb,xb-1) + alpha2 * PW(ye,xe,yb,xb-1);
                            end    
%                         end
%                     end
                end
            end
            
            [Yhat_T, XT] = HEAT_RNN2( NTnodes, NXnodes, NYnodes, Nhiden, Uinit, Wdirect, Wrec, activation_type, biases, -1, test_number);
            error = CostFunction2(Yhat_T, deltaX, deltaY, deltaT, thermal_diffusivity_factor, -1);    
            
            if error==0
                break;
            elseif error < error1
                error1 = error;
                alpha1 = alpha2;
            end
            alpha2 = alpha2 + start_step;
                       
        end
        alpha2 = alpha1;
        
        for ye = 1:1:NYnodes            
            for xe = 1:1:NXnodes
%                 for yb = 1:1:NYnodes
%                     for xb = 1:1:Nhiden
                        yb = ye; xb = xe;                        
                        W(ye,xe,yb,xb) = W(ye,xe,yb,xb) + alpha2 * PW(ye,xe,yb,xb);                        
                        if yb+1 <= NYnodes
                            W(ye,xe,yb+1,xb) = W(ye,xe,yb+1,xb) + alpha2 * PW(ye,xe,yb+1,xb);
                        end
                        if yb-1 > 0
                            W(ye,xe,yb-1,xb) = W(ye,xe,yb-1,xb) + alpha2 * PW(ye,xe,yb-1,xb);
                        end
                        if xb+1 <= NXnodes
                            W(ye,xe,yb,xb+1) = W(ye,xe,yb,xb+1) + alpha2 * PW(ye,xe,yb,xb+1);
                        end    
                        if xb-1 > 0
                            W(ye,xe,yb,xb-1) = W(ye,xe,yb,xb-1) + alpha2 * PW(ye,xe,yb,xb-1);
                        end
%                     end
%                 end
            end
        end
        alpha2      
       
    otherwise
        assert(0);
end


end

