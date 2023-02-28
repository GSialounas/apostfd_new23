function A = createA_diffusion_TRISTAN(dt,x,K,Lx)

N = length(x);


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

F = (x(3:end) - x(1:(end-2))) .* (x(3:end) - x(2:(end-1))) .* (x(2:(end-1)) - x(1:(end-2)));
F  = dt./(0.5*F);


% This is the factor for the first row, which happens to be the same as the
% factor to the last row
F1 = (x(2) - x(1) + x(end) - x(end-1)) * (x(2) - x(1)) * (x(end) - x(end-1));

% Entries of the first row
F1_1    =  (dt/(0.5*F1)) * (-1) * (x(2) - x(1) +x(end) - x(end-1)); % diagonal entry
F1_2    =  (dt/(0.5*F1)) * (x(end) - x(end-1)); % second entry
F1_end  =  (dt/(0.5*F1)) * (x(2) - x(1)); % penultimate entry



% Entries of the last row (the factor F1 is equal to Fend so no renaming, yolo.)
Fend_end       =   (dt/(0.5*F1)) * (-1) * (x(2) - x(1) +x(end) - x(end-1)); % diagonal entry
Fend_end_min1  =   (dt/(0.5*F1)) * (x(2) - x(1)); % penultimate entry
Fend_1         =   (dt/(0.5*F1)) * (x(end)-x(end-1)); % second entry

x_d        =  F* (-1) .* (x(3:end) - x(1:(end-2))); 
x_d        =  [F1_1 ; x_d; Fend_end];
x_dplus1   =  [0; F1_2; F.* (x(2:(end-1)) - x(1:(end-2)))]; % First superdiagonal
x_dminus1  =  [F .* (x(3:(end)) - x(2:(end-1))); Fend_end_min1; 0]; % First subdiagonal

f         =   zeros(size(x_d)); % penultimate sup and sub diagonals
f(1)      =   Fend_1; % penultimate subdiagonal
f(end)  =   F1_end; % penultimate superdiagonal

% Creation of A
A = K * spdiags([f x_dminus1 x_d x_dplus1 f], [-(N-1) -1:1, N-1], N, N);