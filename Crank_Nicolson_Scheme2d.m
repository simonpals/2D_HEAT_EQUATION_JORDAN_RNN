function [ left_scheme_side, right_scheme_side1, right_scheme_side2 ] = Crank_Nicolson_Scheme2d( Yhat_t, Yhat_t_1, yn, xn )

left_scheme_side = (Yhat_t(yn,xn) - Yhat_t_1(yn,xn));
right_scheme_side1 = ( Yhat_t(yn,xn+1) - 2.0*Yhat_t(yn,xn) + Yhat_t(yn,xn-1)+ Yhat_t_1(yn,xn+1) - 2.0*Yhat_t_1(yn,xn) + Yhat_t_1(yn,xn-1) );
right_scheme_side2 = ( Yhat_t(yn+1,xn) - 2.0*Yhat_t(yn,xn) + Yhat_t(yn-1,xn)+ Yhat_t_1(yn+1,xn) - 2.0*Yhat_t_1(yn,xn) + Yhat_t_1(yn-1,xn) );

end

