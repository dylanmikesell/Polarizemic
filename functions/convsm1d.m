function Ktt=convsm1d(A,nx)

% for a nice demo, do the following - smooth a small vector
% using a smoothing window which is 2*2+1 = 5 samples long:
%
% avs = convsm1d(1:6,2)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
% convsm1d                                              %
%                                                       %
% convolutionally smooth a vector                       %
% mmhaney 11/15/2007                                    %
%                                                       %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%   EXPLANATION OF INPUT VARIABLES                      %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A:    a data vector
% nx:   number of data points on one side of the smoothing window %

% NOTES
%
% The window of the convolutional smoother is symmetric;    %
% for example, for nx=2, with an input array f and          % 
% an output array F, the convolutional smoother gives       %
%                                                           %
% F_(i) = (1/5)*[f_(i-2)+f_(i-1)+f_(i)+f_(i+1)+f_(i+2)]     %
%                                                           %
% nx=2 means that the sum ranges from (i-2) to (i+2).       %
%                                                           %
% Since the smoothing is done using Fourier                 %
% transforms, edge effects are avoided by padding at the    %
% edges before smoothing and then removing the padding.     %
% The padding is done by mirroring the data closest to      %
% the edges.                                                %
%                                                           %
% For a nice demo, do the following - smooth a small vector %
% using a smoothing window which is 2*2+1 = 5 samples long: %
%                                                           %
% avs = convsm1d(1:6,2)                                     %
%                                                           %
%
%



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
Ktx = zeros(1,Nxp); 

% expand input do even/odd cases %
kx = [1:Nxp];

if (ceil(Nx/2) == floor(Nx/2))
    Ktx = ((1/(2*nx+1))*((((2*sin(.5*(nx+1)*...
                ((kx-((Nxp+2)/2)+eps)*((2*pi)/Nxp))))./...
                sin(.5*(kx-((Nxp+2)/2)+eps)*((2*pi)/Nxp))).*...
                cos(.5*nx*(kx-((Nxp+2)/2))*((2*pi)/Nxp)))-1));
else
    Ktx = ((1/(2*nx+1))*((((2*sin(.5*(nx+1)*...
                ((kx-((Nxp+1)/2)+eps)*((2*pi)/Nxp))))./...
                sin(.5*(kx-((Nxp+1)/2)+eps)*((2*pi)/Nxp))).*...
                cos(.5*nx*(kx-((Nxp+1)/2))*((2*pi)/Nxp)))-1));
end

%filter %
Ap = real(ifftn(ifftshift(fftshift(fftn(Ap)).*Ktx)));

% pull out the middle part and stick it in vector %
Ktt = zeros(1,Nx);
Ktt = Ap((nx+1):(nx+Nx)); 

% if we were originally handed a column vector, flip to give column back
if (isflip == 1)
else
    Ktt = Ktt';
end

