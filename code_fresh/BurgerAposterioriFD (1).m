clear all
close all
showplot = 0; %boolean for graphics

scheme = 'FTBS'; % 'FTBS' order 1 time, order 1 space upwinding
                 % 'FT2S' order 1 time, order 2 space upwinding
                 % 'FT3S' order 1 time, order 3 space upwinding
                 % 'FTCS' order 1 time, order 2 space central
                 
interpolant = 'pwl'; % 'pwl' the standard piecewise bilinear lagrange interpolant
                     % 'pwq' the piecewise linear(time)quadratic(space) hermite interpolant
                     % 'spl' cubic spline
T = 1.; % final time
maxit = 5; % max refinement iterations
global nq xq wq % gauss quadrature
nq = 2; %number of quad points in 1d per cell
gauss();

%ex = @(x,t) exp(-100*(x-.5 - t).^2); %smooth exact solution 
ex = @(x,t) (x-t > 0.25).*(x-t <= 0.75); %discontinuous exact solution 
%ex = @(x,t) (x-t > 0.25).*(x-t <= 0.75) .* (1-2*abs(2*(x-t)-1)); %C0 exact solution 

fprintf(1,'Using the %s scheme\n',scheme)

for m = 1:maxit

    h = 2^(-(m+2)); % mesh size
    x = 0:h:1; % spatial partition    
    dt = 0.1*h; % timestep size
    uold = ex(x,0); % set initial condition
    if showplot
        fig = figure();
        lnh = plot(x,uold); %approximate in blue
        hold on
        exh = plot(x,ex(x,0),'r'); %exact in red
    end
    t = dt; %initialise time
    i = 0; 
    etat = 0; etac = 0; % accumulation of est
    L2L2R = 0; %L2L2 accumulation of ||R||^2
    while t < T+eps
        i = i + 1;
        if scheme == 'FTCS'
            %u = uold + dt/(2*h)*(circshift(uold,1) - circshift(uold,-1));
        elseif scheme == 'FTBS'
            u = uold + dt/h*(0.5*circshift(uold,1).^2 - 0.5*circshift(uold,0).^2);
        elseif scheme == 'FT2S'
            %u = uold + dt/(2*h)*(u8;
        elseif scheme == 'FT3S'
        else
            fprintf(1,'Invalid scheme choice');
            break
        end
        if showplot
            set(lnh,'XData',x,'YData',u)
            set(exh,'XData',x,'YData',ex(x,t+dt))
            drawnow
        end
        
        %Loop over quad points in time 
        for iq = 1 : nq
            tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
            
            if interpolant == 'pwl'
                [RL1iq,RL2iq] = compute_Rl(x,uold,u,tq(iq),t,dt); %compute R at gauss point
            else
                fprintf(1,'Invalid interpolant')
                break
            end
            
            L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
        end
        etat = etat + sum(abs(u - uold))*dt*h;
        etac = etac + 2*(h+dt)*dt*sum(abs(0.5*circshift(uold,1).^2 - 0.5*circshift(uold,0).^2))...
            + (h+dt)^2*dt*h*(sum(abs(circshift(uold,0))+ abs(circshift(uold,1))));
        t = t + dt; %move in time
        uold = u; 

    end
    [L1, L2] = space_int(x,ex(x,0),0,ex);
    initL1err(m) = abs(L1); %compute final time L2 error
    initL2err(m) = sqrt(L2); %compute final time L2 error
    RL2(m) = sqrt(L2) + sqrt((L2L2R)); %bound
    RL1(m) = T*(space_int(x,ex(x,0),0,ex) + etat + etac) + ...
        2*T * sqrt(etat + etac); %bound

    EOCeL1(1) = 0;
    EOCeL2(1) = 0;
    EOCRL1(1) = 0;
    EOCRL2(1) = 0;
    if m > 1
        EOCeL1(m) = log(initL1err(m-1)/initL1err(m))/log(2);
        EOCeL2(m) = log(initL2err(m-1)/initL2err(m))/log(2);
        EOCRL2(m) = log(RL2(m-1)/RL2(m))/log(2);
        EOCRL1(m) = log(RL1(m-1)/RL1(m))/log(2);
    end
    EI(m) = RL1(m)/initL2err(m);
    %fprintf(1,'||(u - IU)(0)||_L1 = %.5f  EOC = %1.2f\n',initL1err(m),EOCeL1(m))
    %fprintf(1,'||(u - IU)(0)||_L2 = %.5f  EOC = %1.2f\n',initL2err(m),EOCeL2(m))
    fprintf(1,'           est(L1) = %.5f  EOC = %1.2f\n',RL1(m),EOCRL1(m))
    fprintf(1,'      ||R||_L2(L2) = %.5f  EOC = %1.2f\n',RL2(m),EOCRL2(m))
    %fprintf(1,'                EI = %.5f  \n',EI(m))
end


function [L1Rt, L2Rt] = compute_Rl(x,uold,u,evalt,tj,dt)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% for any t it can be represented as a pw linear function in space on each
% spatial element
L2Rt = 0;
L1Rt = 0;
h = x(2) - x(1);

IU_t = (u - uold)/dt; %defined as a constant for any t
IU = (evalt - tj)*uold/dt + (dt+tj-evalt)*u/dt; %Lagrange interpolant of IU over the time interval evaluated at evalt

c1 = IU;
c2 = (circshift(IU,-1) - IU)/h;
c1t = IU_t;
c2t = (circshift(IU_t,-1) - IU_t)/h;

for j = 1 : length(x)-1
    for iq = 1:nq
        xiq = 0.5*h*xq(iq) + x(j) + h/2;
        IUTS = c1(j);
        IUTS_x = c2(j);
        IUTS_t = c1t(j)+c2t(j)*(xiq-x(j));
        L2Rt = L2Rt + wq(iq)*h*(IUTS_t + IUTS.*IUTS_x)^2;
        L1Rt = L1Rt + wq(iq)*h*abs(IUTS_t + IUTS.*IUTS_x);

    end
end

end

%At a specific time level compute L2 spatial error use 2 point gauss quadrature elementwise
function [intL1, intL2] = space_int(x,v,evalt,ex)
global nq xq wq

h = x(2)-x(1);
intL1 = 0;
intL2 = 0;
for j = 1 : length(x)-1
    for iq = 1:nq
        xiq = 0.5*h*xq(iq) + x(j) + h/2;
        v(iq) = (xiq - x(j))*v(j+1)/h + (x(j+1)-xiq)*v(j)/h;
        
        intL1 = intL1 + wq(iq)*h*abs(v(iq) - ex(xiq,evalt));
        intL2 = intL2 + wq(iq)*h*(v(iq) - ex(xiq,evalt))^2.;
    end
end

end

function [] = gauss
% For n point gauss quadrature, return evaluation points and weights for
% gauss quadrature over [-1,1]
global nq xq wq
if nq == 1
    xq = 0;
    wq = 1;
elseif nq == 2
    xq = [-0.57735026918962576451, 0.57735026918962576451];
    wq = [0.5, 0.5];
elseif nq == 3
    xq = [-0.7745966692414834, 0, 0.7745966692414834];
    wq = 0.5*[0.5555555555555556, 0.8888888888888888, 0.5555555555555556];
elseif nq > 3
    fprintf(1,'No Gauss quadrature implemented of this degree\n')
    return
end
end