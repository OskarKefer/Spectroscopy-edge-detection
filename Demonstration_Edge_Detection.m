%% Demonstration of the edge-detection algorithm based on simulated Pump-Probe Data
% Basic introduction on using the edgeDet-function for edge-detection,
% demonstrated on simulated Pump-Probe (PP-) Data. The PP data simulated
% contains a dispserion curve, typical for white-light continuum probes 
% used in Pump-probe experiments. The dispersion is characerized after convolution
% of the PP data with a wavelet and mapping of the maximum value of the
% resulting gaussian curves.
clear
clc
close all

% Initialize time- and wavelength parameters for PP-Data simulation
PP.time   = (-20:.05:50)';
PP.lambda =  400:1:800;

% Simulation of the dispersion curve via polynomial functions
dispersion = polyval([-4e-5 2e-2 0],PP.lambda-580);

% Characteristics of the gaussians used to simulate the spectra of the
% PP-Data. Arranged as follows: A0, Lambda0, width
gausses = [ 1 520 70;...
           .9 540 50;...
            1 580 80;...
           -1 650 50;...
          -.5 720 50];

% Exponential decay characteristics to simulate some form of dynamics
decays  = 100;

% Functions to simulate the gaussian and exponential curves. The
% edge-light(step-like) onset of PP-data is simulated by an error-function,
% which also results from convolution of gaussian-function with exponential
% decay functions.
sigma    = @(x) 2*sqrt(2*x)/log(2);% modifier for width
gaus     = @(x) x(:,1).*exp(-((PP.lambda-x(:,2))./sigma(x(:,3))).^2);% gauss-function
expdecay = @(x) erfc(-(PP.time-x)/.55)/2.*exp(-(PP.time-x)/decays);% exponential-decay function + onset

% Generate PP-data and normalize it
PP.Data = sum(gaus(gausses),1).*expdecay(dispersion); 
PP.Data = PP.Data./max(PP.Data,[],"all");

% Input-parameter for wavelet-construction
width  = .65;% Width of gaussian envelope (Gamma_w)
period = width/1.2;% Period of oscillation (T_w)

% Construct wavelet and perform convolution with PP-data. Read documentation of edgeDet for more information
[convo, time_interp, data_interp, wavelet] = edgeDet(PP.time,PP.Data,[-10 10],length(PP.time(abs(PP.time)<=15)),period,width);


figure(1) % plots "raw"-PP-data and shows relation to wavelet
contourf(PP.lambda,time_interp,data_interp,100,"linecolor","none")
hold on
contour(PP.lambda,time_interp-6,real(wavelet)*ones(1,length(PP.lambda))*.8,100)
hold off
colormap(redblue2)
axis([450 780 -8 8])
caxis([-1 1]*.95)
yline(0,"linewidth",2)
title("Simulated PP Difference Spectrum")
xlabel("Wavelength, nm")
ylabel("Delay, ps")
setfunc(gca,11,"off",false)



figure(2) % plots result from convolution of wavelet with PP-data
contourf(PP.lambda,time_interp,abs(convo)./max(abs(convo)),100,"linecolor","none")
colormap(redblue2)
axis([450 780 -8 8])
caxis([-1 1]*.95)
yline(0,"linewidth",2)
title("Convolution of PP and Wavelet")
xlabel("Wavelength, nm")
ylabel("Delay, ps")
setfunc(gca,11,"off",false)


% Here, we exemplarily characterize the dispersion curve. In this case, it
% is enough to just look at the maximum value of the resulting gaussians
% (from the convolution), which are obtained by taking the absolute value
% of the convolution. (Dispersion can also be characterized by actually
% fitting the gaussian function, which should be more precise)
[~,max_position] = max(abs(convo));

figure(3) % plots simulated dispersion curve and retrieved one. 
% Retrieved by taking the temporal position of the maximum of the 
% absolute-valued result from the convolution
plot(PP.lambda,dispersion,"b","LineWidth",4)
hold on
plot(PP.lambda,time_interp(max_position),"r:","LineWidth",4)
hold off
title("Dispersion Simulated vs. Retrieved")
yline(0,"linewidth",2)
legend(["simulated","retrieved"],"Location","best")
xlabel("Wavelength, nm")
ylabel("Delay, ps")
setfunc(gca,11,"off",false)
axis([450 780 -8 8])

% Take the retrieved dispersion curve to correct the time-zero position at
% each wavelength, resulting in a dispersion-free PP-dataset, ready for
% further analysis
corrected_PP = PP.Data;
for ij = 1:length(PP.lambda)
corrected_PP(:,ij) = interp1(PP.time-time_interp(max_position(ij)),PP.Data(:,ij),PP.time);
end
reduced_time = PP.time;
reduced_time(any(isnan(corrected_PP),2)) = [];  % After time-zero correction, some data-points become NaNs, remove them
corrected_PP(any(isnan(corrected_PP),2),:) = [];% After time-zero correction, some data-points become NaNs, remove them

figure(4) % Plots the corrected PP-dataset
contourf(PP.lambda,reduced_time(abs(PP.time)<=8),corrected_PP(abs(PP.time)<=8,:),100,"linecolor","none")
colormap(redblue2)
axis([450 780 -3 8])
caxis([-1 1]*.95)
yline(0,"linewidth",2)
title("Simulated - and Corrected- PP Difference Spectrum")
xlabel("Wavelength, nm")
ylabel("Delay, ps")
setfunc(gca,11,"off",false)



function setfunc(x,sizefont,minor,major)
set(x,"FontSize",sizefont,"BoxStyle","full","TickDir","out","XMinorTick","on","YMinorTick","on","XMinorGrid",minor,"YMinorGrid",minor,"FontName","Tahoma","Xgrid",major,"Ygrid",major);
end

function c = redblue2(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.
%   Adam Auton, 9th October 2009
if nargin < 1
    m = size(get(gcf,'colormap'),1); 
end
red   = [ linspace(0,1,11)' [.0196 .1294 .2627 .5725 .8196 .9680 .9921 .9569 .8392 .698  .4039]'];
green = [ linspace(0,1,11)' [.1882    .4 .5765 .7725 .8980 .9686 .8588 .6470 .3764 .0941     0]'];
blue  = [ linspace(0,1,11)' [.3803 .6745 .7647 .8705 .9412 .9686 .7803 .5098 .3019 .1686 .1216]'];

interpolator = linspace(0,1,m)';

r = interp1(red(:,1),red(:,2),interpolator);
g = interp1(green(:,1),green(:,2),interpolator);
b = interp1(blue(:,1),blue(:,2),interpolator);
c = [r g b]; 
end
