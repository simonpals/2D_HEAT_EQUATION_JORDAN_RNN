%Calculate derivatives
function [ F_NET,F_W,F_Wrec,F_biases ] = F_NET2d(Nhiden,NXnodes,NYnodes,F_Yhat,W,Wrec,xt,F_NET,F_NETrec,F_W,F_Wrec,activation_type, factorPW,biases,F_biases)

F_x = zeros(NYnodes, Nhiden+NXnodes);

for y = 1:1:NYnodes
    for x = 1:1:NXnodes
        F_x(y,x+Nhiden) = F_Yhat(y, x);
    end
end

for ye = NYnodes:-1:1
    for xe = Nhiden+NXnodes:-1:Nhiden+1                
        F_x(ye,xe) = F_x(ye,xe) + Wrec(ye,xe-NXnodes,ye,xe-Nhiden) * F_NETrec(ye,xe-NXnodes);
        if ye-1 > 0 
            F_x(ye,xe) = F_x(ye,xe) + Wrec(ye-1,xe-NXnodes,ye,xe-Nhiden) * F_NETrec(ye-1,xe-NXnodes);
        end        
        if ye+1 <= NYnodes                        
            F_x(ye,xe) = F_x(ye,xe) + Wrec(ye+1,xe-NXnodes,ye,xe-Nhiden) * F_NETrec(ye+1,xe-NXnodes);
        end        
        if xe-1-NXnodes > 0 %if xe-1 > 0   
            F_x(ye,xe) = F_x(ye,xe) + Wrec(ye,xe-1-NXnodes,ye,xe-Nhiden) * F_NETrec(ye,xe-1-NXnodes);
        end        
        if xe+1-NXnodes <= Nhiden %if xe+1 <= NXnodes           
            F_x(ye,xe) = F_x(ye,xe) + Wrec(ye,xe+1-NXnodes,ye,xe-Nhiden) * F_NETrec(ye,xe+1-NXnodes);
        end                  
        
        F_NET(ye,xe) = F_x(ye,xe) * fActivation_derivative(xt(ye,xe),'none');        
        F_biases(ye,xe-Nhiden) = F_biases(ye,xe-Nhiden) + F_x(ye,xe);
                       
        
        F_W(ye,xe-Nhiden,ye,xe-Nhiden) = F_W(ye,xe-Nhiden,ye,xe-Nhiden) + F_NET(ye,xe) * xt(ye,xe-Nhiden);
        
        if ye-1 > 0
            F_W(ye,xe-Nhiden,ye-1,xe-Nhiden) = F_W(ye,xe-Nhiden,ye-1,xe-Nhiden) + F_NET(ye,xe) * xt(ye-1,xe-Nhiden);            
        end
        
        if ye+1 <= NYnodes            
            F_W(ye,xe-Nhiden,ye+1,xe-Nhiden) = F_W(ye,xe-Nhiden,ye+1,xe-Nhiden) + F_NET(ye,xe) * xt(ye+1,xe-Nhiden);            
        end
        
        if xe-1-Nhiden > 0 %if xe-1 > 0
            F_W(ye,xe-Nhiden,ye,xe-1-Nhiden) = F_W(ye,xe-Nhiden,ye,xe-1-Nhiden) + F_NET(ye,xe) * xt(ye,xe-1-Nhiden);            
        end
        
        if xe+1-Nhiden <= NXnodes %if xe+1 <= NXnodes
            F_W(ye,xe-Nhiden,ye,xe+1-Nhiden) = F_W(ye,xe-Nhiden,ye,xe+1-Nhiden) + F_NET(ye,xe) * xt(ye,xe+1-Nhiden);                        
        end                
        
    end
end

Yhat = zeros(NYnodes,NXnodes);
for y = 1:1:NYnodes
    for x = 1:1:NXnodes
        Yhat(y,x) = xt(y,Nhiden+x);
    end
end

for ye = NYnodes:-1:1
    for xe = Nhiden:-1:1
       
%Commented because only one neuron in hidden correspond 
%only one neuron in output layer
%         for yb = 1:1:NYnodes
%             for xb = Nhiden+1:1:Nhiden+NXnodes
%                 F_x_term_2 = W(yb,xb-Nhiden,ye,xe) * F_NET(yb,xb);
%                 F_x(ye,xe) = F_x(ye,xe) + F_x_term_2;
%             end
%         end
%         F_x_term_2 = W(ye,xe,ye,xe) * F_NET(ye,xe+Nhiden);
%         F_x(ye,xe) = F_x(ye,xe) + F_x_term_2;

        F_x(ye,xe) = F_x(ye,xe) + W(ye,xe,ye,xe) * F_NET(ye,xe+Nhiden);
        if ye-1 > 0
            F_x(ye,xe) = F_x(ye,xe) + W(ye-1,xe,ye,xe) * F_NET(ye-1,xe+Nhiden);
        end        
        if ye+1 <= NYnodes
            F_x(ye,xe) = F_x(ye,xe) + W(ye+1,xe,ye,xe) * F_NET(ye+1,xe+Nhiden);
        end
        if xe-1 > 0
            F_x(ye,xe) = F_x(ye,xe) + W(ye,xe-1,ye,xe) * F_NET(ye,xe-1+Nhiden);
        end
        if xe+1 <= NXnodes
            F_x(ye,xe) = F_x(ye,xe) + W(ye,xe+1,ye,xe) * F_NET(ye,xe+1+Nhiden);
        end        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%

        F_NET(ye,xe) = F_x(ye,xe) * fActivation_derivative(xt(ye,xe),activation_type);
        
%Commented because we train only finite difference patern weigths         
%         for yb = 1:1:NYnodes
%             for xb = Nhiden+1:1:Nhiden+NXnodes
%                 F_Wrec_term1 = F_NETrec(ye,xe) * xt(yb,xb);
%                 F_Wrec_term1 = F_Wrec_term1 + F_Wrec(ye,xe,yb,xb-Nhiden);
%                 F_Wrec(ye,xe,yb,xb-Nhiden) = F_Wrec_term1;
%             end
%         end        

        F_Wrec(ye,xe,ye,xe) = F_Wrec(ye,xe,ye,xe) + F_NETrec(ye,xe) * Yhat(ye,xe);
        
        if ye-1 > 0
            F_Wrec(ye,xe,ye-1,xe) = F_Wrec(ye,xe,ye-1,xe) + F_NETrec(ye,xe) * Yhat(ye-1,xe);
        end
        
        if ye+1 <= NYnodes
            F_Wrec(ye,xe,ye+1,xe) = F_Wrec(ye,xe,ye+1,xe) + F_NETrec(ye,xe) * Yhat(ye+1,xe);
        end
        
        if xe-1 > 0
            F_Wrec(ye,xe,ye,xe-1) = F_Wrec(ye,xe,ye,xe-1) + F_NETrec(ye,xe) * Yhat(ye,xe-1);
        end

        if xe+1 <= NXnodes
            F_Wrec(ye,xe,ye,xe+1) = F_Wrec(ye,xe,ye,xe+1) + F_NETrec(ye,xe) * Yhat(ye,xe+1);
        end
        
    end
end



end

