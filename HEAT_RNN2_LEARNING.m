function [W,Wrec,biases] = HEAT_RNN2_LEARNING(optimization_method,T,maximum_passes,learning_rate,NTnodes,deltaT,NXnodes,deltaX,NYnodes,deltaY,Nhiden,Uinit,Uleft,Uright,Utop,Ubottom,thermal_diffusivity_factor,activation_type,biases, test_number, new_init_W, W, Wrec, ExactSol)
              

if new_init_W == true  
    
    noiseValue = 0; % rand() * 0.04; % 0.05; % 0.1942; %   0.0689; %  0.0554;
    noiseValue
    WfactorMLS = 0;
    Wfactor = zeros(NYnodes,NXnodes); % Wrec
    W2factor = zeros(NYnodes,NXnodes); % W
    ye = 2;
    xe = 2;

    Wfactor(ye,xe) = 0.01;
    W2factor(ye,xe) = 0.1;

    if ye+1 <= NYnodes
        Wfactor(ye+1,xe) = 0.1; %0.125;
        W2factor(ye+1,xe) = 0.2; % 0.125;
    end
    if ye-1 > 0
        Wfactor(ye-1,xe) = 0.1; % 0.125;
        W2factor(ye-1,xe) = 0.2; % 0.125;
    end
    if xe+1 <= NXnodes
        Wfactor(ye,xe+1) = 0.1; % 0.125;
        W2factor(ye,xe+1) = 0.2; % 0.125;
    end 
    if xe-1 > 0
        Wfactor(ye,xe-1) = 0.1; % 0.125;
        W2factor(ye,xe-1) = 0.2; % 0.125;
    end


    [W] = Init4d_Distribution(NYnodes,NXnodes,NYnodes,Nhiden, false, W2factor);
    [Wrec] = Init4d_Distribution(NYnodes,Nhiden,NYnodes,NXnodes, true, Wfactor);
end

PW = repmat(0, [NYnodes NXnodes NYnodes Nhiden]);
PWrec = repmat(0, [NYnodes Nhiden NYnodes NXnodes]);

F_W_1 = repmat(0, [NYnodes Nhiden NYnodes NXnodes]);
F_W_DIR_1 = repmat(0, [NYnodes Nhiden NYnodes NXnodes]);

if length(ExactSol) > 1 % && ExactSol > -1 
    exact = ExactSol;
else    
    [exact] = Exact_solution( test_number );
end 
 
 Yhat_T_1 = 0;
 MLS3rec = 10000000;
 
for pass_number = 1:1:maximum_passes
    
    disp('Pass number:');
    disp(pass_number);
    
    [ Yhat_T, XT ] = HEAT_RNN2( NTnodes, NXnodes, NYnodes, Nhiden, Uinit, W, Wrec, activation_type, biases, -1, test_number );    
    if test_number == 3
        [ Max_Deviation3, MLS3 ] = Fdistance2d( rot90(Yhat_T(:,:,11)), exact(:,:,11) );
        if MLS3 < 0.058
            break;
        end
    elseif test_number == 4
        [ Max_Deviation3, MLS3 ] = Fdistance2d( (Yhat_T(:,:,11)), exact(:,:,11) );
    elseif test_number == 5
        [ Max_Deviation3, MLS3 ] = Fdistance2d( (Yhat_T(:,:,11))', exact(:,:,11) );
        if MLS3 < 0.025
            break;
        end        
    elseif test_number == 6
        [ Max_Deviation3, MLS3 ] = Fdistance2d( (Yhat_T(:,:,11))', exact(:,:,11) );
        if MLS3 < 0.00005
            break;
        end        
    elseif test_number == 7
        [ Max_Deviation3, MLS3 ] = Fdistance2d( (Yhat_T(:,:,11))', exact(:,:,11) );
        if MLS3 < 0.001
            break;
        end
    end
    MLS3
    if MLS3rec < MLS3        
        break;
    end
    MLS3rec = MLS3;
    
    error = CostFunction2(Yhat_T, deltaX, deltaY, deltaT, thermal_diffusivity_factor, -1);
    error
    
    if pass_number > 1
        [ Max_Deviation_Yhat,Yhat_T_MLS ] = Fdistance2d( Yhat_T_1(:,:,11), Yhat_T(:,:,11));
        Yhat_T_MLS
        Max_Deviation_Yhat
    end
    Yhat_T_1 = Yhat_T;
       
    
    F_W = repmat(0, [NYnodes NXnodes NYnodes Nhiden]);
    F_NET = repmat(0, [NYnodes Nhiden+NXnodes]);
    F_NETrec = repmat(0, [NYnodes Nhiden+NXnodes]);
    F_Wrec = repmat(0, [NYnodes Nhiden NYnodes NXnodes]);
    Yhat_t_1 = zeros(NYnodes,NXnodes);    
    W3 = W; Wrec3 = Wrec;
    F_biases = repmat(0, [NYnodes NXnodes]);
        
    init_exit_cond = 0;
    factorPWrec = 1;
    for t = NTnodes:-1:1
        xt = XT(:,:,t);
        Yhat_t = Yhat_T(:,:,t);
               
        
        if t == 1
            [F_Yhat] = zeros(NYnodes,NXnodes);
        else
            Yhat_t_1 = Yhat_T(:,:,t-1);                        
            [F_Yhat] = CostFunctionDerivarive2d(NTnodes,NXnodes,NYnodes,t,Yhat_t,Yhat_t_1,deltaX,deltaY,deltaT,thermal_diffusivity_factor,Uleft,Uright,Utop,Ubottom,false);
            
            if t < NTnodes
                [F_Yhat_1] = CostFunctionDerivarive2d(NTnodes,NXnodes,NYnodes,t,Yhat_T(:,:,t+1),Yhat_T(:,:,t),deltaX,deltaY,deltaT,thermal_diffusivity_factor,Uleft,Uright,Utop,Ubottom,true);
                F_Yhat = F_Yhat + F_Yhat_1;
            end
        end                        
                
        exit_cond = norm(F_Yhat); % F_Yhat' * F_Yhat;
        if exit_cond > init_exit_cond
            init_exit_cond = exit_cond;
        end

        
        [ F_NET,F_W,F_Wrec,F_biases ] = F_NET2d(Nhiden,NXnodes,NYnodes,F_Yhat,W,Wrec,xt,F_NET,F_NETrec,F_W,F_Wrec,activation_type,factorPWrec,biases,F_biases);
        
        F_NETrec = F_NET;
    end
    
    init_exit_cond

Wdirect = W;
if max(max(biases)) ~= 0
    [biases, pb] = Optimization_method2d( 4, learning_rate, Wrec, F_Wrec, PWrec, NXnodes, NYnodes, Nhiden,  pass_number, F_W_1, Wdirect, deltaX, deltaY, deltaT, thermal_diffusivity_factor, Uinit, activation_type, NTnodes, biases, F_biases, test_number );
end
% uncomment for W optimization
%[W,PW] = Optimization_method2d( 1, learning_rate, W, F_W, PW, NXnodes, NYnodes, Nhiden );
[W,PW] = Optimization_method2d( optimization_method, learning_rate, W, F_W, PW, NXnodes, NYnodes, Nhiden,  pass_number, F_W_DIR_1, Wdirect, deltaX, deltaY, deltaT, thermal_diffusivity_factor, Uinit, activation_type, NTnodes, biases, F_biases, test_number );
[Wrec,PWrec] = Optimization_method2d( optimization_method, learning_rate, Wrec, F_Wrec, PWrec, NXnodes, NYnodes, Nhiden,  pass_number, F_W_1, Wdirect, deltaX, deltaY, deltaT, thermal_diffusivity_factor, Uinit, activation_type, NTnodes, biases, F_biases, test_number );
F_W_1 = F_Wrec;
F_W_DIR_1 = F_W;

end %pass_number
end

