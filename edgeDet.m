function [output, time_interp, data_interp, wavelets] = edgeDet(time, data, timerange, num_points, periods, widths)
% This function can be used to detect the onset of data (reimagined as an
% edge-detection problem) in any kind of time-resolved data (1-|2-D) by
% convolution of wavelets (periodically oscillating function with a
% gaussian envelope (in time domain) with the original dataset.
% If you are using this method, please cite the corresponding publication (DOI:XXXXXX.XXX)

% Author: Oskar Kefer, for any inquiry, please contact me!
% E-Mail: oskar.kefer@uni-heidelberg.de

% Inputs:
%        - time:        Vector of time-axis corresponding to the dataset
%        - data:        Dataset of interest, {1-D | N-D} representation
%        - timerange:   Part of the time-axis to be considered (2-elements)
%        - num_points:  How many points to interpolate within timerange (1-element)
%        - periods:     Periods of wavelets, {1-|N-} elements
%        - widths:      Widths of gaussian envelope, {1-|N-} elements 

time = time(:);% time runs along first column
data = convert3D(data,time);% permute data to max with dimension of time-axis

% Generate time axis for interpolation (constant-step time necessary for FFT)
time_interp = linspace(timerange(1),timerange(2),num_points); 
time_interp = time_interp(:);

data_interp = interp1(time,data,time_interp);% Interpolate data onto newly generated time-axis
data_interp(isnan(data_interp)) = 0;% Remove NaNs

% Data is concatenated to ensure artifact-free Fourier-Transform later
data_interp2 = cat(1,ones(size(data_interp)).*data_interp(1,:,:,:),data_interp,ones(size(data_interp)).*data_interp(end,:,:,:));

%optional: Pass data over low-pass filter
% data_interp2 = lowpass(data_interp.',1e-15).'; 

% optional: Smooth data before convolution (gets rid of some high-frequency noise)
% for ij=1:size(data_interp,2)
%     data_interp2(:,ij) = smooth(data_interp(:,ij),15);
% end

% Initialize periods and widths of the wavelet and permute to correct
% dimensionality. Periods in 3-rd position and Widths in 4-th position
periods = permute(periods(:),[3 2 1]);
widths  = permute(widths(:),[4 3 2 1]);

% Generate wavelets and zero-pad to match size of dataset
wavelets = exp(-1i*(pi/2+(time_interp-mean(time_interp))./periods)-((time_interp-mean(time_interp))./widths).^2);
wavelets = cat(1,zeros(size(wavelets)),wavelets,zeros(size(wavelets)));


% Perform covolution of data and wavelets. Output can be 4-D array if
% multiple periods & widths are applied. Afterwards remove concatenated
% parts to have size match with newly generated time-axis
output = ifft(fft(data_interp2,[],1).*fft(fftshift(wavelets,1),[],1));
output = output(length(time_interp)+1:end-length(time_interp),:,:,:);

% Reduce wavelet array to match length of time-axis accordingly
wavelets = wavelets(length(time_interp)+1:end-length(time_interp),:,:,:);
end

function Matrix3D = convert3D(Input, varargin)
% Converts the dataset to be in line with the axes (Axis1,Axis2) specified. 
% Can be used with 1,2 and 3...x dimensional datasets.
% Input dataset will always be sorted in descending mannor when nothing is
% specified.

prelimsize       = size(Input);
dimensionality   = length(prelimsize);
decision         = zeros(dimensionality,1);
permutationorder = zeros(dimensionality,1);

if isempty(varargin)
    decision = sort(prelimsize,"descend");
else
    for ij = 1:length(varargin)
        if all(length(varargin{ij}) ~= 1) && ~isempty(varargin{ij})
            decision(ij) = length(varargin{ij});
        elseif all(size(varargin{ij}) == 1)
            decision(ij) = varargin{ij};
        end
    end
    decision(decision==0) = sort(prelimsize(~ismember(prelimsize,decision)),"descend");
end

for ij = 1:dimensionality
    permutationorder(ij) = find(floor(decision(ij)) == prelimsize,1);
end

Matrix3D = permute(Input,permutationorder);
end