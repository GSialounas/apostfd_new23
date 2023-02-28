function A = createA_advection_TRISTAN_nonu_central(U, dt, x,Lx)

% This function will create the BTCS matrix

% _____INPUT_______
% U  : velocity (take it as constant here)
% dt : time-step
% x  : x-grid
% Lx : final point of the x-grid(which coincides with first point in a
% periodic domain
% Length of the solution vector - same as number of grid-points
N = length(x);

% Comments
% Bear in mind that the grid (x) which is being fed to this function is
% actually x(1:end-1) of the original grid.  The reason is because the
% domain is periodic so the point x(end) and the point x(1) are actually
% the same point.  In order to work around this when constructing the
% operator matrix we use Lx as a stand in for x(end) of THE ENTIRE GRID 
% where in x(end) is NOT x(1) e.g. x = 0:0.1:10, x(1) = 0, x(end)= 10! 
% For example, x(end) used below is actually x(end-1).

% Transpose x if necessary:
% The functionality in this function works if x is a column vector, so if
% it's not just transpose it:
if (size(x,2) > size(x,1))
    x = x';
end

% Creation of vectors that will hold diagonal entries
% Note that for f we don't care about the entries where they don't appear
% so we just put them in the same vector instead of creating separate
% vectors for the penultimate super and sub diagonals
x_dplus1  =     U * dt * [0 ; ((x(2)-x(1)).^(-1)-(x(2)-x(1) + (Lx-x(end))).^(-1)) ; ((x(3:end)-x(2:(end-1))).^(-1)- (x(3:end)-x(1:(end-2))).^(-1))]; % First superdiagonal

x_d       =    U*dt* [ ((Lx-x(end))^(-1) - (x(2)-x(1))^(-1) ) ; (((x(2:end-1)-x(1:(end-2))).^(-1)- (x(3:end)-x(2:(end-1))).^(-1))); ( (x(end)-x(end-1))^(-1) - (Lx-x(end))^(-1) )] ; % Main Diagonal Entries
% x_dplus1  =     U * dt * [0 ; (x(2)-x(1) + (x(end)-x(end-1))).^(-1) ; (x(3:end)-x(1:(end-2))).^(-1)]; % First superdiagonal
x_dminus1 =      U * dt * [((x(3:end)-x(1:(end-2))).^(-1) - (x(2:end-1)-x(1:(end-2))).^(-1)); ((x(end)-x(end-1) + (Lx-x(end)))^(-1) - (x(end)-x(end-1))^(-1) ); 0]; % First subdiagonal
% x_dminus1 =     - U * dt * [(x(3:end)-x(1:(end-2))).^(-1) ; (x(2)-x(1) + (x(end)-x(end-1)))^(-1) ; 0]; % First subdiagonal
f         =     zeros(size(x_d)); % penultimate sup and sub diagonals
f(1)      =     U * dt * ((Lx-x(end))^(-1) - ((Lx-x(end)+x(end)-x(end-1))^-1  )); % penultimate subdiagonal
% f(1)      =     U * dt * (x(2)-x(1) + (x(end)-x(end-1)))^(-1); % penultimate subdiagonal
f(end)  =     U * dt * (-(Lx-x(end))^(-1) + ((Lx-x(end)+x(end)-x(end-1))^-1  )); % penultimate superdiagonal
% f(end)  =     -U * dt *  (x(2)-x(1) + (x(end)-x(end-1)))^(-1); % penultimate superdiagonal

% Creation of A
A = spdiags([f x_dminus1 x_d x_dplus1 f], [-(N-1) -1:1 N-1], N, N);


