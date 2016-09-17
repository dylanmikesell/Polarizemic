%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% polar_coherency.m
%
% PROGRAMMER:
% Matt Haney
%
% Last revision date:
% 22 May 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program polar_coherency is a Matlab function to perform complex-valued 
% polarization analysis using a coherency matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
% dtac          a matrix of size (3 x Nt), where Nt is the number of time 
%               samples; the three rows of the matrix are the vertical, 
%               east, and north components in that order
% wndo          the length of the window to calculate polarization 
%               parameters in samples
%
% Output:
% azim          time series of azimuth angles in clockwise degrees from 
%               north
% incd          time series of incidence angles in degrees from vertical
% ellip         time series of ellipticity - intermediate axis of ellipsoid 
%               divided by major axis
%
% Dependency:
% Uses program convsm1d.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [azim, incd, ellip] = polar_coherency(dtac,wndo)

% unpack the matrix of three-component motion
pv6z = dtac(1,:);
pv6e = dtac(2,:);
pv6n = dtac(3,:);

% length of the time series
tln = length(pv6z);

% remove local mean over the window length
pv6z = pv6z-convsm1d(pv6z,wndo);
pv6n = pv6n-convsm1d(pv6n,wndo);
pv6e = pv6e-convsm1d(pv6e,wndo);

% take hilbert transform
pv6zh = hilbert(pv6z);
pv6nh = hilbert(pv6n);
pv6eh = hilbert(pv6e);

% form the diagonal components of the coherency matrix and 
% smooth over a window length - these are all real-valued
mvec11 = convsm1d(real(pv6zh.*conj(pv6zh)),wndo);
mvec22 = convsm1d(real(pv6nh.*conj(pv6nh)),wndo);
mvec33 = convsm1d(real(pv6eh.*conj(pv6eh)),wndo);

% form the off-diagonal components of the coherency matrix and 
% smooth over a window length - these are all complex-valued
mvec12 = convsm1d(real(pv6zh.*conj(pv6nh)),wndo) +...
       1i*convsm1d(imag(pv6zh.*conj(pv6nh)),wndo);
mvec13 = convsm1d(real(pv6zh.*conj(pv6eh)),wndo) +...
       1i*convsm1d(imag(pv6zh.*conj(pv6eh)),wndo);
mvec23 = convsm1d(real(pv6nh.*conj(pv6eh)),wndo) +...
       1i*convsm1d(imag(pv6nh.*conj(pv6eh)),wndo);

% initialize azim, ellip, and incd
azim = zeros(1,tln);
ellip = azim; 
incd = azim;

% loop over all time samples
for ii=1:tln

    % solve for the eigenvalues and eigenvectors of the 
    % coherency matrix
    [v2,d] = eig([ mvec11(ii) mvec12(ii) mvec13(ii);
                   conj(mvec12(ii)) mvec22(ii) mvec23(ii);
                   conj(mvec13(ii)) conj(mvec23(ii)) mvec33(ii) ]);         

    % sort the eigenvalues into ascending order
    [ds, indx] = sort(diag(d),1,'ascend');
    %d = d(:,indx);
    v2 = v2(:,indx);
    
    % phase angle associated with maximum real part of largest eigenvector
    % Vidale (1986, BSSA) recommends a line search to do this, but it 
    % can be done analytically, without need for a line search
    nmr = -2*((real(v2(1,3))*imag(v2(1,3)))+...
                  (real(v2(2,3))*imag(v2(2,3)))+...
                  (real(v2(3,3))*imag(v2(3,3))));
    dmr = -(((imag(v2(1,3))^2)-(real(v2(1,3))^2)) + ...
                ((imag(v2(2,3))^2)-(real(v2(2,3))^2)) + ...
                ((imag(v2(3,3))^2)-(real(v2(3,3))^2)));
    alph = atan2(nmr,dmr);
    
    % rotate the largest eigenvector to maximize the real part
    v2rv = v2(1:3,3)*exp(1i*alph);
    % the maximum real part of largest eigenvector
    v2r = sqrt((real(v2rv(1))^2) + (real(v2rv(2))^2) + (real(v2rv(3))^2));

    % calculate the azimuth and incidence angle
    v2 = real(v2rv); %v2 = real(v2(1:3,3));
%    azim(ii) = atan2(v2(3),v2(2))*(180/pi);
    incd(ii) = atan2(sqrt(v2(2)^2 + v2(3)^2),v2(1))*(180/pi);
    azim(ii) = atan(v2(3)/v2(2))*(180/pi);
%    incd(ii) = atan(sqrt(v2(2)^2 + v2(3)^2)/v2(1))*(180/pi);    
    
    % calculate the ellipticity as in Vidale (1986, BSSA)
    ellip(ii) = sqrt(1-(real(v2r)^2))/real(v2r);

end
