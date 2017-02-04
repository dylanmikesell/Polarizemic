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
% Last modified by Dylan Mikesell (19 December 2016)
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

%--------------------------------------------------------------------------
% Determine window to calculate polarization over, in samples
cycs = 2; % number of cycles (2 to 3 is usually sufficient)
wndo = floor( (1/f) * (1/dt) ) * cycs; % samples per cycle times # of cycles

%--------------------------------------------------------------------------
[dtac,tt] = makeSynthetic(tln,dt,omga); % compute the synthetics

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
plot(tt,azim,'LineWidth',lSize); %axis([tt(1) tt(end) 60 75]); 
grid on;
axis([tt(1) tt(end) 50 100]); 
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
plot(tt,azim,'LineWidth',lSize); %axis([tt(1) tt(end) 60 75]); 
grid on;
axis([tt(1) tt(end) 50 100]); 
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



