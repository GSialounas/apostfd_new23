function ParasitesFD(parasite)
%clear all
close all
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

ex = @(x,t) 0.5*(sech(0.5*((x+5)-t))).^2;

for m = 4:4

    h = 0.1; % mesh size
    if parasite == 'me'
        x = [-20:h:20];%, 5+h/2:h/2:7, 7+h:h:20]; % spatial partition
    else
        x = [-.5:h:1.5]; % spatial partition
        ep = h*(x>0.5);
    end
    dt = h^2; % timestep size
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
        %ylim([-0.2,1]);
        subplot(1,2,2)
        lnh2 = plot(x,uold,'LineWidth',2); %approximate in blue
        hold on 
        xlim([-0.4,0.7]);
        ylim([-0.05,0.05]);
        hold on6*circshift(uold,0).*
        pause
    
    end
    t = dt; %initialise time
    i = 0; 
    L2L2R = 0; %L2L2 accumulation of ||R||^2
    while t < T+eps
        i = i + 1;
        
        u = uold - dt*6*circshift(uold,0).*(circshift(uold,1) - circshift(uold,-1))./(circshift(x,1) - circshift(x,-1)) ...
            - dt * 0.5*(circshift(uold,2) - 2*circshift(uold,1)+2*circshift(uold,-1)-circshift(uold,-2))./(circshift(x,0) - circshift(x,-1)).^3;
  
        if showplot && mod(i,1) == 0
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