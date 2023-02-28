clear all
close all
showplot = 0; %boolean for graphics

scheme = 'FTBS'; % 'FTBS' order 1 time, order 1 space upwinding
                 % 'FT2S' order 1 time, order 2 space upwinding
                 % 'FT3S' order 1 time, order 3 space upwinding
                 % 'FTCS' order 1 time, order 2 space central
                 
interpolant = 'pwl'; % 'pwl' the standard piecewise bilinear lagrange interpolant
                     % 'spl' cubic spline
T = .1; % final time
maxit = 7; % max refinement iterations
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
    dt = 0.1*h^2; % timestep size
    uold = ex(x,0); % set initial condition
    if showplot
        fig = figure();
        lnh = plot(x,uold); %approximate in blue
        hold on
        exh = plot(x,ex(x,0),'r'); %exact in red
    end
    t = dt; %initialise time
    i = 0; 
    L2L2R = 0; %L2L2 accumulation of ||R||^2
    while t < T+eps
        i = i + 1;
        if scheme == 'FTCS'
            u = uold + dt/(2*h)*(circshift(uold,1) - circshift(uold,-1));
        elseif scheme == 'FTBS'
            u = uold + dt/h*(circshift(uold,1) - circshift(uold,0));
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
                RL2iq = compute_R(x,uold,u,tq(iq),t,dt); %compute R at gauss point
            elseif interpolant == 'spl'
                RL2iq = compute_Rs(x,uold,u,tq(iq),t,dt); %compute R at gauss point
            else
                fprintf(1,'Invalid interpolant')
                break
            end
            
            L2L2R = L2L2R + wq(iq)*dt*(RL2iq*exp(-tq(iq))); %quadrature formula
        end
        t = t + dt; %move in time
        uold = u; 

    end
    finalL2err(m) = sqrt(space_int(x,u,T,ex)); %compute final time L2 error
    R(m) = sqrt(exp(T)*(space_int(x,ex(x,0),0,ex) + L2L2R)); %bound
    EOCe(1) = 0;
    EOCR(1) = 0;
    if m > 1
        EOCe(m) = log(finalL2err(m-1)/finalL2err(m))/log(2);
        EOCR(m) = log(R(m-1)/R(m))/log(2);
    end
    EI(m) = R(m)/finalL2err(m);
    fprintf(1,'||(u - IU)(T)||_L2 = %.5f  EOC = %1.2f\n',finalL2err(m),EOCe(m))
    fprintf(1,'      ||R||_L2(L2) = %.5f  EOC = %1.2f\n',R(m),EOCR(m))
    fprintf(1,'                EI = %.5f  \n',EI(m))
end

function L2Rt = compute_Rs(x,uold,u,evalt,tj,dt)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% Take IU as a cubic spline
% for any t it can be represented as a pw linear function in space on each
% spatial element
L2Rt = 0;
h = x(2) - x(1);
u = spline(x,u);
[breaks,coefs,l,k,d] = unmkpp(u);
du = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);

uold = spline(x,uold);
[breaks,coefs,l,k,d] = unmkpp(uold);
duold = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
%xx=0:h/100:1;
%yy = ppval(u, xx);
%plot(xx,yy)
for j = 1 : length(x)-1
    for iq = 1:nq
        xl(iq) = 0.5*h*xq(iq) + x(j) + h/2;
        
        IUt = (ppval(u,xl(iq)) - ppval(uold,xl(iq)))/dt; %defined as a constant for any t
        IUx = (evalt - tj)*(ppval(duold,xl(iq)))/dt + (dt+tj-evalt)*(ppval(du,xl(iq)))/dt; %Lagrange interpolant of IU over the time interval evaluated at evalt
        
        %IUt(iq) = (xl(iq) - x(j))*IUt(j+1)/h + (x(j+1)-xl(iq))*IUt(j)/h;
        
        L2Rt = L2Rt + wq(iq)*h*(IUt + IUx)^2;
    end
end

end

function L2Rt = compute_R(x,uold,u,evalt,tj,dt)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% for any t it can be represented as a pw linear function in space on each
% spatial element
L2Rt = 0;
h = x(2) - x(1);

IUt = (u - uold)/dt; %defined as a constant for any t
IU = (evalt - tj)*uold/dt + (dt+tj-evalt)*u/dt; %Lagrange interpolant of IU over the time interval evaluated at evalt
for j = 1 : length(x)-1
    for iq = 1:nq
        xiq = 0.5*h*xq(iq) + x(j) + h/2;
        vl = IU(j);
        vr = IU(j+1);
        IUx = (vr - vl)/h;
        IUtx(iq) = (xiq - x(j))*IUt(j+1)/h + (x(j+1)-xiq)*IUt(j)/h;
        
        L2Rt = L2Rt + wq(iq)*h*(IUtx(iq) + IUx)^2;
    end
end

end

%At a specific time level compute L2 spatial error use 2 point gauss quadrature elementwise
function int = space_int(x,v,evalt,ex)
global nq xq wq

h = x(2)-x(1);
int = 0;
for j = 1 : length(x)-1
    for iq = 1:nq
        xiq = 0.5*h*xq(iq) + x(j) + h/2;
        v(iq) = (xiq - x(j))*v(j+1)/h + (x(j+1)-xiq)*v(j)/h;
        
        int = int + wq(iq)*h*(v(iq) - ex(xiq,evalt))^2.;
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