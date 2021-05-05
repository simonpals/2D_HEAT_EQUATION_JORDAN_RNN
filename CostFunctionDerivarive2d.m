function [ F_Yhat ] = CostFunctionDerivarive2d(NTnodes,NXnodes,NYnodes,t,Yhat_t,Yhat_t_1,deltaX,deltaY,deltaT,thermal_diffusivity_factor,Uleft,Uright,Utop,Ubottom, FY_t_1)

F_Yhat = zeros(NYnodes,NXnodes);

assert(NXnodes>=6);
assert(NYnodes>=6);

rfx = 0.5*(thermal_diffusivity_factor*deltaT)/(deltaX^2);
rfy = 0.5*(thermal_diffusivity_factor*deltaT)/(deltaY^2);

y_sign = 1;

if true == FY_t_1
    y_sign = -1;
end


%(2,2)
[left_scheme_side, right_scheme_side1, right_scheme_side2] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, 2, 2 );
F_Yhat(2, 2) = (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (y_sign*1+2*rfx+2*rfy);
[left_scheme_side, right_scheme_side1, right_scheme_side2] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, 2, 3 );
F_Yhat(2, 2) = F_Yhat(2, 2) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfx);
[left_scheme_side, right_scheme_side1, right_scheme_side2] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, 3, 2 );
F_Yhat(2, 2) = F_Yhat(2, 2) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfy);

%(2,NYnodes-1)
[left_scheme_side, right_scheme_side1, right_scheme_side2] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, NYnodes-1, 2 );
F_Yhat(NYnodes-1, 2) = (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (y_sign*1+2*rfx+2*rfy);
[left_scheme_side, right_scheme_side1, right_scheme_side2] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, NYnodes-2, 2 );
F_Yhat(NYnodes-1, 2) = F_Yhat(NYnodes-1, 2) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfy);
[left_scheme_side, right_scheme_side1, right_scheme_side2] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, NYnodes-1, 3 );
F_Yhat(NYnodes-1, 2) = F_Yhat(NYnodes-1, 2) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfx);

%(NXnodes-1,NYnodes-1)
[left_scheme_side, right_scheme_side1, right_scheme_side2] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, NYnodes-1, NXnodes-1 );
F_Yhat(NYnodes-1, NXnodes-1) = (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (y_sign*1+2*rfx+2*rfy);
[left_scheme_side, right_scheme_side1, right_scheme_side2] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, NYnodes-2, NXnodes-1 );
F_Yhat(NYnodes-1, NXnodes-1) = F_Yhat(NYnodes-1, NXnodes-1) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfy);
[left_scheme_side, right_scheme_side1, right_scheme_side2] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, NYnodes-1, NXnodes-2 );
F_Yhat(NYnodes-1, NXnodes-1) = F_Yhat(NYnodes-1, NXnodes-1) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfx);

%(NXnodes-1,2)
[left_scheme_side, right_scheme_side1, right_scheme_side2] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, 2, NXnodes-1 );
F_Yhat(2, NXnodes-1) = (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (y_sign*1+2*rfx+2*rfy);
[left_scheme_side, right_scheme_side1, right_scheme_side2] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, 3, NXnodes-1 );
F_Yhat(2, NXnodes-1) = F_Yhat(2, NXnodes-1) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfy);
[left_scheme_side, right_scheme_side1, right_scheme_side2] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, 2, NXnodes-2 );
F_Yhat(2, NXnodes-1) = F_Yhat(2, NXnodes-1) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfx);


%(2,3:NYnodes-2)
for yn = 3:1:NYnodes-2
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, yn, 2 );
    F_Yhat(yn, 2) = (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (y_sign*1+2*rfx+2*rfy);    
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, yn, 3);
    F_Yhat(yn, 2) = F_Yhat(yn, 2) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfx);    
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, yn-1, 2);
    F_Yhat(yn, 2) = F_Yhat(yn, 2) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfy);
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, yn+1, 2);
    F_Yhat(yn, 2) = F_Yhat(yn, 2) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfy);
end

%(NXnodes-1,3:NYnodes-2)
for yn = 3:1:NYnodes-2
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, yn, NXnodes-1 );
    F_Yhat(yn, NXnodes-1) = (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (y_sign*1+2*rfx+2*rfy);    
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, yn, NXnodes-2);
    F_Yhat(yn, NXnodes-1) = F_Yhat(yn, NXnodes-1) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfx);    
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, yn-1, NXnodes-1);
    F_Yhat(yn, NXnodes-1) = F_Yhat(yn, NXnodes-1) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfy);
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, yn+1, NXnodes-1);
    F_Yhat(yn, NXnodes-1) = F_Yhat(yn, NXnodes-1) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfy);
end

%(3:NXnodes-2,2)
for xn = 3:1:NXnodes-2
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, 2, xn );
    F_Yhat(2, xn) = (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (y_sign*1+2*rfx+2*rfy);        
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, 3, xn );
    F_Yhat(2, xn) = F_Yhat(2, xn) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfy);    
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, 2, xn-1 );
    F_Yhat(2, xn) = F_Yhat(2, xn) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfx);
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, 2, xn+1 );
    F_Yhat(2, xn) = F_Yhat(2, xn) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfx);
end

%(3:NXnodes-2,NYnodes-1)
for xn = 3:1:NXnodes-2
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, NYnodes-1, xn );
    F_Yhat(NYnodes-1, xn) = (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (y_sign*1+2*rfx+2*rfy);        
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, NYnodes-2, xn );
    F_Yhat(NYnodes-1, xn) = F_Yhat(NYnodes-1, xn) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfy);    
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, NYnodes-1, xn-1 );
    F_Yhat(NYnodes-1, xn) = F_Yhat(NYnodes-1, xn) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfx);
    [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, NYnodes-1, xn+1 );
    F_Yhat(NYnodes-1, xn) = F_Yhat(NYnodes-1, xn) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfx);
end


%(3:NXnodes-2,3:NYnodes-2)
for yn = 3:1:NYnodes-2
    for xn = 3:1:NXnodes-2
        [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, yn, xn );
        F_Yhat(yn, xn) = (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (y_sign*1+2*rfx+2*rfy);  
        
        [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, yn-1, xn );
        F_Yhat(yn, xn) = F_Yhat(yn, xn) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfy);
        [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, yn+1, xn );
        F_Yhat(yn, xn) = F_Yhat(yn, xn) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfy);
        
        [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, yn, xn-1 );
        F_Yhat(yn, xn) = F_Yhat(yn, xn) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfx);
        [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, yn, xn+1 );
        F_Yhat(yn, xn) = F_Yhat(yn, xn) + (left_scheme_side - rfx*right_scheme_side1 - rfy*right_scheme_side2) * (-rfx);
    end
end

end

