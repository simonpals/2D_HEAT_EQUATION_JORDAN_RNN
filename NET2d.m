function [x_t,Yhat] = NET2d(NXnodes,NYnodes,Nhiden,W,Wrec,xt_1,activation_type,biases)

x_t = zeros(NYnodes, Nhiden+NXnodes);
Yhat = zeros(NYnodes,NXnodes);

for ye = 1:1:NYnodes
    for xe = 1:1:Nhiden
        net = 0;
        
        yb = ye; xb = xe;        
        net = net + Wrec(ye, xe, yb, xb) * xt_1(yb, xb+Nhiden);
        if yb+1 <= NYnodes
            net = net + Wrec(ye, xe, yb+1, xb) * xt_1(yb+1, xb+Nhiden);
        end
        if yb-1 > 0
            net = net + Wrec(ye, xe, yb-1, xb) * xt_1(yb-1, xb+Nhiden);
        end
        if xb+1 <= Nhiden
            net = net + Wrec(ye, xe, yb, xb+1) * xt_1(yb, xb+Nhiden+1);
        end
        if xb-1 > 0
            net = net + Wrec(ye, xe, yb, xb-1) * xt_1(yb, xb+Nhiden-1);
        end        
                
        x_t(ye,xe) = fActivation(net, activation_type);
    end
end

for ye = 1:1:NYnodes
    for xe = Nhiden+1:1:Nhiden+NXnodes
        net = 0;
        
        yb = ye; xb = xe;
        net = net + W(ye,xe-Nhiden,yb, xb-NXnodes) * x_t(yb, xb-NXnodes);
                
        if yb+1 <= NYnodes
            net = net + W(ye,xe-Nhiden,yb+1, xb-NXnodes) * x_t(yb+1, xb-NXnodes);
        end
        if yb-1 > 0
            net = net + W(ye,xe-Nhiden,yb-1, xb-NXnodes) * x_t(yb-1, xb-NXnodes);
        end
        if xb+1 <= Nhiden+NXnodes
           net = net + W(ye,xe-Nhiden,yb, xb-NXnodes+1) * x_t(yb, xb-NXnodes+1); 
        end
        if xb-1 > Nhiden
            net = net + W(ye,xe-Nhiden,yb, xb-NXnodes-1) * x_t(yb, xb-NXnodes-1);
        end                        
        
        x_t(ye,xe) = fActivation(net, 'none') + biases(ye,xe-Nhiden);
    end
end

for y = 1:1:NYnodes
    for x = 1:1:NXnodes
        Yhat(y,x) = x_t(y,Nhiden+x);
    end
end

end