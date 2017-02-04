%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% polar_covariance.m
%
% PROGRAMMER:
% Matt Haney
%
% Last revision date:
% 22 May 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program polar_covariance is a Matlab function to perform real-valued 
% polarization analysis using a covariance matrix
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
function [azim incd ellip] = polar_covariance(dtac,wndo)

% unpack the matrix of three-component motion
pv6z = dtac(1,:);
pv6e = dtac(2,:);
pv6n = dtac(3,:);

% length of the time series
tln = length(pv6z);

% remove local mean over the window length
pv6zh = pv6z-convsm1d(pv6z,wndo);
pv6nh = pv6n-convsm1d(pv6n,wndo);
pv6eh = pv6e-convsm1d(pv6e,wndo);

% form the diagonal components of the covariance matrix and 
% smooth over a window length
mvec11 = convsm1d((pv6zh.*pv6zh),wndo);
mvec22 = convsm1d((pv6nh.*pv6nh),wndo);
mvec33 = convsm1d((pv6eh.*pv6eh),wndo);

% form the off-diagonal components of the covariance matrix and 
% smooth over a window length
mvec12 = convsm1d((pv6zh.*pv6nh),wndo);
mvec13 = convsm1d((pv6zh.*pv6eh),wndo);
mvec23 = convsm1d((pv6nh.*pv6eh),wndo);

% initialize azim, ellip, and incd
azim = zeros(1,tln);
ellip = azim; 
incd = azim;

% loop over all time samples
for ii=1:tln

    % solve for the eigenvalues and eigenvectors of the 
    % covariance matrix
    [v2,d] = eig([ mvec11(ii) mvec12(ii) mvec13(ii);
                   mvec12(ii) mvec22(ii) mvec23(ii);
                   mvec13(ii) mvec23(ii) mvec33(ii) ]);         

    % sort the eigenvalues into ascending order
    [ds indx] = sort(diag(d),1,'ascend');
    d = d(:,indx);
    v2 = v2(:,indx);
    
    % calculate the azimuth and incidence angle
    azim(ii) = atan2((v2(3,3)*1),(v2(2,3)*1))*(180/pi);
    incd(ii) = atan2(sqrt(v2(2,3)^2 + v2(3,3)^2),v2(1,3))*(180/pi);
    
    % calculate the ellipticity
    ellip(ii) = real(sqrt(d(2,2)/d(3,3)));
    
end
