function [ exact ] = Exact_solution( task_number )

nx=21;
ny=21; 
nt=11; 

a = 1;
b = 1;
T = 0.05;

dt=T/(nt-1);   %0.005;   
dx=a/(nx-1);
dy=b/(ny-1);

t=0:dt:T;
x=0:dx:a;                      
y=0:dy:b;

Yhat_T = repmat(0.0, [nx ny nt]);
series_part_cnt = 100;

switch(task_number)    
    case 3
        for it=1:11
            for yi=1:ny
                for xi=1:nx
                    Yhat_T(xi,yi,it) = 0.0;
                    k = ((1.0/8.0));
                    assert(a==b);
                    coef = k / a^2;
                    
                    for m=1:series_part_cnt
                        for n=1:series_part_cnt
                            Yhat_T(xi,yi,it) = Yhat_T(xi,yi,it) + 4 * 0.5 * ( ( -1/(pi*n)*cos(n*pi/2) + 1/(pi*n) )*( -1/(pi*m)*cos(m*pi/1) + 1/(pi*m) ) * sin(m*pi/a*x(xi)) * sin(n*pi/b*y(yi))*exp((-pi^2)*(m^2+n^2)*t(it)*coef) );
                        end
                    end
                    
                end
            end
        end
   
    case 4 % invalid
        sigma = 1.0;
        tdf = 0.1;
        Tmax = 1.0;
        for it=1:11
            for yi=1:ny
                for xi=1:nx        
                    Yhat_T(xi,yi,it) = Tmax/((1+4.0*t(it)*tdf/(sigma^2))^0.5) * exp((-(x(xi)^2 + y(yi)^2))/(sigma^2+4*t(it)*tdf));
                end
            end
        end
        
    case 5
        for it=1:11
            for yi=1:ny
                for xi=1:nx
                    Yhat_T(xi,yi,it) = 0.0; %x(xi)*y(yi); %  x(xi)*(1-x(xi))*y(yi)*(1-y(yi));% 
                    k = ((1.0/5.0));
                    
                    for m=1:series_part_cnt
                        Yhat_T(xi,yi,it) = Yhat_T(xi,yi,it) + (a*b*((-1)^(m+1)))/(m*pi) *sin((m*pi*x(xi))/a)*exp(-((m*pi/a)^2)*k*t(it));
                        for n=1:series_part_cnt                            
                            Yhat_T(xi,yi,it) = Yhat_T(xi,yi,it) + (4*a*b*((-1)^m)*(1-(-1)^n))/(m*n^2*pi^3) * sin((m*pi*x(xi))/a) * cos((n*pi*y(yi))/b)* exp(-((m*pi/a)^2+(n*pi/b)^2)*k*t(it));
                        end
                    end
                    
                end
            end
        end   
        
    case 6
        for it=1:11
            for yi=1:ny
                for xi=1:nx
                    Yhat_T(xi,yi,it) = 0.0; 
                    k = 1/2; % ((1.0/3.0));
                    
                    for m=1:series_part_cnt                        
                        for n=1:series_part_cnt   
                            lambda = -pi^2*(n^2+m^2);
                            Yhat_T(xi,yi,it) = Yhat_T(xi,yi,it) +  (16.0 / (pi^6)) * ( ((-1)^n-1)*((-1)^m-1)*exp(k*lambda*t(it)) / (n^3*m^3) * sin(n*pi*x(xi)) * sin(m*pi*y(yi)) );
                        end
                    end
                    
                end
            end
        end
        
    case 7
        L = 1;
        H = 1;
        k = 1/3;
        for it=1:11
            for yi=1:ny
                for xi=1:nx
                    Yhat_T(xi,yi,it) = 0.0;
                    
                    for m=1:series_part_cnt
                        for n=1:series_part_cnt
                            Yhat_T(xi,yi,it) = Yhat_T(xi,yi,it) + 1.6/((2*n-1)*(2*m-1)*pi^2)*sin(((2*n-1)*pi*x(xi))/L)*sin(((2*m-1)*pi*y(yi))/(2*H)) * exp( -((((2*m-1)^2)/(4*H^2)) + (((2*n-1)^2)/(L^2)))*k*pi^2*t(it) );
                        end
                    end
                    
                end
            end
        end
    otherwise
        assert(0);        
end

exact = Yhat_T;

end

