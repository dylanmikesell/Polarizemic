%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% convsm1d.m
%
% PROGRAMMER:
% Matt Haney
%
% Last revision date:
% 22 May 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program convsm1d is a Matlab function to convolutionally smooth a 
% vector or time series
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
% A             vector or time series
% nx            number of gridpoints on one side of the smoothing arm
%
% Output:
% K33           a smoothed version of the original vector A
%
% Notes:
% The arm of the convolutional smoother is symmetric; for example, 
% for nx=2, with an input vector A and an output vector K33, the 
% convolutional smoother gives:       
%                                                           
% K33_(i) = A_(i-2) + A_(i-1) + A_(i) + A_(i+1) + A_(i+2)     
%                                                           
% so that nx=2 means that the sum ranges from (i-2) to (i+2).        
%                                                           
% Since the smoothing is done using Fourier transforms, edge effects are 
% avoided by padding at the edges before smoothing and then removing the 
% padding. Reflecting boundary conditions are used at the edges of the
% vector A to define the padding.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K33=convsm1d(A,nx)

% see if the input vector is a row or a column
% if it is a row, change it to a column, do everything, 
% then change it back to a row at the end for output
ln = size(A);
isflip = 0; % a flag which says if vector was flipped
if (ln(1) == 1)
    A = A';
    isflip = 1;
else
end

% start
ln = size(A);
Nx = ln(1);

% the padded input %
Ap = zeros(1,Nx+2*nx);

% put A in the middle of it %
Ap((nx+1):(nx+Nx)) = A;

% now pad the input with reflective boundaries the length of an arm %
Ap(1:nx) = Ap((nx+2):(2*nx+1));
Ap((nx+Nx+1):(Nx+2*nx)) = Ap(Nx:(nx+Nx-1));

% the padded filter %
Nxp = Nx+2*nx;

% expand input do even/odd cases %
kx = (1:Nxp);

% two cases - one for even number of samples, one for odd number
if (ceil(Nx/2) == floor(Nx/2))
    K3x = ((1/(2*nx+1))*((((2*sin(.5*(nx+1)*...
                ((kx-((Nxp+2)/2)+eps)*((2*pi)/Nxp))))./...
                sin(.5*(kx-((Nxp+2)/2)+eps)*((2*pi)/Nxp))).*...
                cos(.5*nx*(kx-((Nxp+2)/2))*((2*pi)/Nxp)))-1));
else
    K3x = ((1/(2*nx+1))*((((2*sin(.5*(nx+1)*...
                ((kx-((Nxp+1)/2)+eps)*((2*pi)/Nxp))))./...
                sin(.5*(kx-((Nxp+1)/2)+eps)*((2*pi)/Nxp))).*...
                cos(.5*nx*(kx-((Nxp+1)/2))*((2*pi)/Nxp)))-1));
end

%filter %
Ap = real(ifftn(ifftshift(fftshift(fftn(Ap)).*K3x)));

% pull out the middle part and stick it in a vector 
K33 = Ap((nx+1):(nx+Nx)); 

% if we were originally handed a column vector, flip to give column back
if (isflip == 1)
else
    K33 = K33';
end



