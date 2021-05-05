%Initialize weights matrix with normal distributed random values
function [ W ] = Init4d_Distribution( in_y_count, in_x_count, out_y_count, out_x_count, is_fdm_arch, Wfactor )

W = repmat(0, [in_y_count in_x_count out_y_count out_x_count]);
small_value = 0.00;

if true == is_fdm_arch
    assert(in_x_count == out_x_count);
    
    for ye = 1:1:in_y_count
        for xe = 1:1:in_x_count
            for yb = 1:1:out_y_count
                for xb = 1:1:out_x_count
                    
                    
                    if ye == yb && xe == xb
                        
                        W(ye,xe,yb,xb) = Wfactor(2,2);                         
                        
                        if yb+1 <= out_y_count
                            W(ye,xe,yb+1,xb) = Wfactor(3,2);
                        end
                        if yb-1 > 0
                            W(ye,xe,yb-1,xb) = Wfactor(1,2);
                        end
                        
                        if xb+1 <= in_x_count                           
                            W(ye,xe,yb,xb+1) =  Wfactor(2,3);
                        end
                        
                        if xb-1 > 0
                            W(ye,xe,yb,xb-1) = Wfactor(2,1);
                        end

                    end
                    
                end
            end
            
        end
    end
    
else
    for ye = 1:1:in_y_count
        for xe = 1:1:in_x_count
          for yb = 1:1:out_y_count
                for xb = 1:1:out_x_count                                        
                    if ye == yb && xe == xb                                           
                        W(ye,xe,yb,xb) = Wfactor(2,2);                                           
                        if yb+1 <= out_y_count                            
                            W(ye,xe,yb+1,xb) = Wfactor(3,2);
                        end
                        if yb-1 > 0                           
                            W(ye,xe,yb-1,xb) = Wfactor(1,2);
                        end                        
                        if xb+1 <= in_x_count                           
                            W(ye,xe,yb,xb+1) = Wfactor(2,3);
                        end                        
                        if xb-1 > 0
                            W(ye,xe,yb,xb-1) = Wfactor(2,1);
                        end
                    end                    
                end
          end

        end
    end
end

end

