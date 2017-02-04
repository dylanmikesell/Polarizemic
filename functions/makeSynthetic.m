function [data,tt] = makeSynthetic(tln,dt,omga)


% make a synthetic time series of a rectilinear signal followed by an
% elliptical signal


%--------------------------------------------------------------------------
% a vector of time values for the time series
tt = (0:tln-1) .* dt; % [s]

%--------------------------------------------------------------------------
% two tapers - one for the first half of the time series,
% the other for the second half
tpr = [hanning(tln/2)' zeros(1,tln/2)];
tpr2 = [zeros(1,tln/2) hanning(tln/2)'];

%--------------------------------------------------------------------------
% the rectilinear signal - all 3 components in phase
%
% Ampplitude of each component
ampZ = 1;
ampE = 2;
ampN = 1;
% Phase of each component [rad]
phsZ = 0 * (pi/180);
phsE = 0 * (pi/180);
phsN = 0 * (pi/180);
% Build time series and taper the edges with first taper
w1Z = ampZ * cos( omga*tt + phsZ ) .* tpr;
w1E = ampE * cos( omga*tt + phsE ) .* tpr;
w1N = ampN * cos( omga*tt + phsN ) .* tpr;

%--------------------------------------------------------------------------
% the elliptical signal - vertical component 90 degrees out of phase
%
% Ampplitude of each component
ampZ = 1;
ampE = 3;
ampN = 1;
% Phase of each component [rad]
phsZ = -90 * (pi/180);
phsE =   0 * (pi/180);
phsN =   0 * (pi/180);
% Build time series and taper the edges with second taper
w2Z = ampZ * cos( omga*tt + phsZ ) .* tpr2;
w2E = ampE * cos( omga*tt + phsE ) .* tpr2;
w2N = ampN * cos( omga*tt + phsN ) .* tpr2;

%--------------------------------------------------------------------------
% Add rectilinear and elliptical signals
Z = w1Z + w2Z;
E = w1E + w2E;
N = w1N + w2N;

%--------------------------------------------------------------------------
% pack the complete three-component signal into a matrix with the ordering
% 1,2,3 = Z,E,N
data(1,:) = Z;
data(2,:) = E;
data(3,:) = N;

end