function ParasitesFD(parasite)
%clear all
%close all
showplot = 1; %boolean for graphics

parasite = 'me'; % 'me' mesh 
                % 'mo' model 

scheme = 'FTCS'; % 'FTBS' order 1 time, order 1 space upwinding
                 % 'FT2S' order 1 time, order 2 space upwinding
                 % 'FT3S' order 1 time, order 3 space upwinding
                 % 'FTCS' order 1 time, order 2 space central
                 
T = 1.; % final time
maxit = 4; % max refinement iterations
nq = 2; %number of quad points in 1d per cell

ex = @(x,t) exp(-100*(x-.25 - t).^2); %smooth exact solution 
%ex = @(x,t) (x-t > 0.25).*(x-t <= 0.75); %discontinuous exact solution 
%ex = @(x,t) (x-t > 0.25).*(x-t <= 0.75) .* (1-2*abs(2*(x-t)-1)); %C0 exact solution 


%fprintf(1,'Using the %s scheme\n',scheme)


for m = 4:4

    h = 2^(-(m+2)); % mesh size
    if parasite == 'me'
        x = [-.5:h:0.5, 0.5+h/8:h/8:1.5]; % spatial partition
    else
        x = [-.5:h:1.5]; % spatial partition
        ep = h*(x>0.5);
    end
    dt = 0.1*h^2; % timestep size
    uold = ex(x,0); % set initial condition
    if showplot
        fig = figure();
        subplot(1,2,1)
        lnh = plot(x,uold,'LineWidth',2); %approximate in blue
        hold on 
        if parasite == 'me'
        plot(x,-0.1*ones(length(x)),'. k'); %approximate in blue
        else
        plot(x,ep,'. k'); %approximate in blue
        end
        ylim([-0.2,1]);
        subplot(1,2,2)
        lnh2 = plot(x,uold,'LineWidth',2); %approximate in blue
        hold on 
        xlim([-0.4,0.7]);
        ylim([-0.05,0.05]);
        hold on
        pause
    
    end
    t = dt; %initialise time
    i = 0; 
    L2L2R = 0; %L2L2 accumulation of ||R||^2
    while t < T+eps
        i = i + 1;
        if scheme == 'FTCS'
            if parasite == 'me'
                u = uold + dt*(circshift(uold,1) - circshift(uold,-1))./(circshift(x,-1) - circshift(x,1));
            else
                u = uold + dt*(circshift(uold,1) - circshift(uold,-1))./(circshift(x,-1) - circshift(x,1)) ...
                    + dt*ep.*(circshift(uold,1) - 2*circshift(uold,0) + circshift(uold,-1))./(circshift(x,1) - circshift(x,0)).^2;
            end
        elseif scheme == 'FTBS'
            u = uold + dt/h*(circshift(uold,1) - circshift(uold,0));
        else
            fprintf(1,'Invalid scheme choice');
            break
        end
        if showplot && mod(i,100) == 0
            set(lnh,'XData',x,'YData',u)
            set(lnh2,'XData',x,'YData',u)
            %set(exh,'XData',x,'YData',ex(x,t+dt))
            drawnow
        end
        
        t = t + dt; %move in time
        uold = u; 

    end
end

end