function[] = animateU2d(U, thermal_diffusivity_factor, X0, deltaX, Y0, deltaY, T0, deltaT)

    dim = size(U);
    dim_size = size(dim');
    assert(dim_size(1)==3); % third dim it's time
    
    for i=1:1:dim(1)
        x(i) = X0+(i-1)*deltaX;
    end
    
    for i=1:1:dim(2)
        y(i) = Y0+(i-1)*deltaY;
    end
    
    for nTime = 1:1:dim(3)
        Ut = U(:,:,nTime);
        h=surf(x,y,Ut,'DisplayName','Ut');
        shading interp
        axis ([0 2 0 2 0 2])
          pause(0.2)
        title({['2-D Diffusion with {\nu} = ',num2str(thermal_diffusivity_factor)];['time (\itt) = ',num2str(T0+nTime*deltaT)]})
        xlabel('Spatial co-ordinate (x) \rightarrow')
        ylabel('{\leftarrow} Spatial co-ordinate (y)')
        zlabel('Transport property profile (u) \rightarrow')
        drawnow;
        refreshdata(h)
    end

end