%=============================================================*
%                                                             *
% run_1D.m                                                    *
% Version: 2.0 (created on 6 July 2021)                       *
% AUTHORS: M.J. Colebank, M.U. Qureshi,  M.S. Olufsen         *
%          Department of Mathematics                          *
%          North Carolina State University, Raleigh, USA      *
% DATE UPDATED: 31 Dec. 2019.                                 *
%                                                             *
% DESCRIPTION: This script creates the inteface with the C++  *
% code for the 1D fluid dynamics model. The script creates    *
% and runs the executable by passing selected model           *
% parameters and plots the hemodynamic waveforms              *
%                                                             *
%=============================================================*



% Model usines a simple linear elastic wall model and a three element Windkessl outflow
% boundary condition. The system is driven by imposing an inflow profile.

% clc; clear all; close all;
% !chmod +x sor06
% !make clean
% !make
format shortg;
%% Define the parameters
% The stiffness is assumed to follow an exponential curve, i.e.
% Eh/r0 = k1*exp(k2*r0) + k3
% where E, h, and r0 are the Youngs Modulus, wall thickness, and reference
% radius, respectively. To set a single stiffness value for all the
% vessels, set k1=0;
k1 = 0;%1e+5;
k2 = -25;
k3 =  8e4;%1e4%2.3e3;%5e+4;
%% Create an Inflow waveform for pulse wave simlation in a single vessel
% Comment this section out if you have a flow waveform you wish to provide
% as an input to the model
Qmin = 0.0;    %d minumum flow rate during diastole (ml/s)
Qmax = 500;  % maximum flow rate during systole (ml/s)


t0 = 0.0; % initial time (s)
T  = 0.9; % Length of the cardiac cycle (s)

Tm   = 0.15; % Peak systole time (s)
Td   = 0.4;  % End systole time (s)
Tneg = 0; % Period of negative flow

tt = Td-Tm; % time difference between peak and end of systole
ttneg = Tneg-Td;
[Qin,tQ] = InFlow(Qmin,Qmax,t0,Tm,tt,ttneg,T); % Generates a desired inflow profile and saves it to a file Qin_*.dat
% OR load in a flow waveform
% Qin = load('MYFLOWPROFILE.dat');

% or use exponential flow
% x = linspace(0,T,length(Qin));
% f = @(par) par(1).*exp(-(par(2)-x).^2./0.005);
% Qin = f([300,0.3]);
fname = strcat('Qin.dat');
dlmwrite(fname,Qin');
    

%% Load in a pressure profile
% This code was developed for parameter estimation purposes. If no pressure
% profile is available, simply provide a systolic, diastolic, and mean
% pressure value below
% Pdat = load('MYPRESSUREPROFILE.dat');
% Psys = max(Pdat); Pdia = min(Pdat); Pmean = mean(Pdat);

% Or hard code the systolic, diastolic, and mean value
Psys = 25; Pdia = 5; Pmean = (25+2.*5)/3;
P_vein = 0;
Pdat = [Psys, Pmean, Pdia,P_vein];
%% Specify Windkessel boundary
% Mihaela: these are still passed in but not used
Rp_1 = 0.7;
Rd_1 = 1.5;
CT_1  = 1.8;

% These are the new parameters. P0 is the reference pressure in the
% pressure-area relationship (I would just fix it to a value between 0 and
% 8) and p_out is the pressure boundary condition, which gets converted to
% area
p0 =8.0;
p_out = 16.0;
%% Define scaling factors for resistance/compliance if needed
pars = [k3,Rp_1,Rd_1,CT_1,p0,p_out,1];
pars_str = mat2str(pars);
%% call the model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%unix(sprintf('sor06.exe  %s',pars_str(2:end-1)));
unix(sprintf('./sor06  %s',pars_str(2:end-1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting Simulated data
% Plot inlet/outlet waveforms in the whole tree
% data = load('art_ALL.2d');
data = dlmread('output_1.2d');
[time,x,p,q,A,C] = gnuplot(data);
t = time(:,1)-time(1,1); % Time starts from 0
% Proximal predictions are 1->N
% Midpoint predictions are N+1->2N
% Distal   prdictions are 2N+1->3N (or end)
figure(1); clf; hold on;
% plot(linspace(0,T,length(Pdat)),Pdat,'k','LineWidth',3);
plot(t,p(:,1),'LineWidth',2);
set(gca, 'FontSize',30);grid on;
xlabel('Time (s)');
ylabel('Pressure (mmHg)');
plot(t,t.*0+Psys,'--k','LineWidth',2);
plot(t,t.*0+Pdia,'--k','LineWidth',2)
hold off;

figure(2); hold on;
plot(p(:,[1 end]));
legend('Inlet','Outlet')


