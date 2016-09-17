%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% polariz.m
%
% PROGRAMMER:
% Matt Haney
%
% Last revision date:
% 22 May 2009
% Last modified by Dylan Mikesell (13 September 2016)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program polariz is a Matlab script that compares polarization 
% estimates based on covariance and coherence matrices 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is a script
clear all
close all
clc

addpath('./functions');

%% Step 1: Make synthetic 3-component times series

% make a synthetic time series of a rectilinear signal followed by an
% elliptical signal 
%--------------------------------------------------------------------------
% length of synthetic time series
tln = 10000;
%--------------------------------------------------------------------------
% dominant frequency of synthetic time series in Hz
f = 2; % [Hz]
omga = 2*pi*f; % [radians/s]
%--------------------------------------------------------------------------
% time sample interval
dt = 0.001; % [s]
% window to calculate polarization over, in samples
cycs = 2; % number of cycles
wndo = floor( (1/f) * (1/dt) ) * cycs; % samples per cycle times # of cycles
%--------------------------------------------------------------------------
% a vector of time values for the time series
tt = (1:tln) .* dt; % [s]
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
ampZ = 1;
% Phase of each component [rad]
phsZ = 0 * (pi/180);
phsE = 0 * (pi/180);
phsN = 0 * (pi/180);
% Build time series and taper the edges with first taper
w1Z = ampZ * cos( omga*tt + phsZ ) .* tpr;
w1E = ampE * cos( omga*tt + phsE ) .* tpr;
w1N = ampZ * cos( omga*tt + phsN ) .* tpr;
%--------------------------------------------------------------------------
% the elliptical signal - vertical component 90 degrees out of phase
%
% Ampplitude of each component
ampZ = 1;
ampE = 3;
ampZ = 1;
% Phase of each component [rad]
phsZ = -90 * (pi/180);
phsE =   0 * (pi/180);
phsN =   0 * (pi/180);
% Build time series and taper the edges with second taper
w2Z = ampZ * cos( omga*tt + phsZ ) .* tpr2;
w2E = ampE * cos( omga*tt + phsE ) .* tpr2;
w2N = ampZ * cos( omga*tt + phsN ) .* tpr2;
%--------------------------------------------------------------------------
% Add rectilinear and elliptical signals 
Z = w1Z + w2Z;
E = w1E + w2E;
N = w1N + w2N;
%--------------------------------------------------------------------------
% pack the complete three-component signal into a matrix with the ordering
% 1,2,3 = Z,E,N 
dtac(1,:) = Z;
dtac(2,:) = E;
dtac(3,:) = N;
%--------------------------------------------------------------------------
% plot the three-component synthetic data
lSize = 2;
h = figure('Color','white');
subplot(3,1,1)
plot(tt,dtac(1,:),'LineWidth',lSize); axis([tt(1) tt(end) -3.5 3.5]); grid on;
title('Three-component synthetic data');
ylabel('Vertical [Z]');
subplot(3,1,2)
plot(tt,dtac(2,:),'LineWidth',lSize); axis([tt(1) tt(end) -3.5 3.5]); grid on;
ylabel('East [E]');
subplot(3,1,3)
plot(tt,dtac(3,:),'LineWidth',lSize); axis([tt(1) tt(end) -3.5 3.5]); grid on;
ylabel('North [N]'); xlabel('Time [s]');
%--------------------------------------------------------------------------
% Fix figure properties
fsize = 16;
% set( findall( h, '-property', 'FontWeight' ), 'FontWeight', 'Bold' );
set( findall( h, '-property', 'Fontsize' ), 'Fontsize', fsize );
%
%% Apply polar coherency method
%
% calculate the polarization parameters using the complex-valued coherency
% method 
[azim, incd, ellip] = polar_coherency( dtac, wndo );
%--------------------------------------------------------------------------
% plot the coherency-based polarization parameters
lSize = 2;
h = figure('Color','white');
subplot(3,1,1);
plot(tt,azim,'LineWidth',lSize); axis([tt(1) tt(end) 60 75]); grid on;
title('Coherency matrix method (Vidale, BSSA, 1986)');
ylabel('Azimuth [deg]');
subplot(3,1,2)
plot(tt,incd,'LineWidth',lSize); axis([tt(1) tt(end) 50 100]); grid on;
ylabel('Incididence [deg]');
subplot(3,1,3)
plot(tt,ellip,'LineWidth',lSize); axis([tt(1) tt(end) 0 1]); grid on;
xlabel('Time [s]'); ylabel('Ellipticity');
%--------------------------------------------------------------------------
% Fix figure properties
fsize = 16;
% set( findall( h, '-property', 'FontWeight' ), 'FontWeight', 'Bold' );
set( findall( h, '-property', 'Fontsize' ), 'Fontsize', fsize );
%
%% Apply polar covariance method
%
% calculate the polarization parameters using the real-valued covariance
% method 
[azim, incd, ellip] = polar_covariance( dtac, wndo );
%--------------------------------------------------------------------------
% plot the covariance-based polarization parameters
lSize = 2;
h = figure('Color','white');
subplot(3,1,1)
plot(tt,azim,'LineWidth',lSize); axis([tt(1) tt(end) 60 75]); grid on;
title('Covariance matrix method');
ylabel('Azimuth [deg]');
subplot(3,1,2)
plot(tt,incd,'LineWidth',lSize); axis([tt(1) tt(end) 50 100]); grid on;
ylabel('Incididence [deg]');
subplot(3,1,3)
plot(tt,ellip,'LineWidth',lSize); axis([tt(1) tt(end) 0 1]); grid on;
xlabel('Time [s]'); ylabel('Ellipticity');
%--------------------------------------------------------------------------
% Fix figure properties
fsize = 16;
% set( findall( h, '-property', 'FontWeight' ), 'FontWeight', 'Bold' );
set( findall( h, '-property', 'Fontsize' ), 'Fontsize', fsize );



