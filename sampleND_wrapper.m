%% sampleND_wrapper.m
%
% DISTRIBUTION STATEMENT A. Approved for public release. Distribution is 
% unlimited.
% 
% This material is based upon work supported by the Federal Aviation 
% Administration under Air Force Contract No. FA8702-15-D-0001. Any 
% opinions, findings, conclusions or recommendations expressed in this 
% material are those of the author(s) and do not necessarily reflect the 
% views of the Federal Aviation Administration.
% 
% © 2024 Massachusetts Institute of Technology.
% 
% Subject to FAR52.227-11 Patent Rights - Ownership by the contractor 
% (May 2014)
% 
% The software/firmware is provided to you on an As-Is basis
% 
% Delivered to the U.S. Government with Unlimited Rights, as defined in 
% DFARS Part 252.227-7013 or 7014 (Feb 2014). Notwithstanding any copyright
% notice, U.S. Government rights in this work are defined by DFARS 
% 252.227-7013 or DFARS 252.227-7014 as detailed above. Use of this work 
% other than as specifically authorized by the U.S. Government may violate 
% any copyrights that exist in this work.
%
% sample_ND_wrapper.m
%
% Wrapper script illustating the usage of inverseTransformSampleND.m
% function using 3-dimension synthetic data, including
% - data preparation
% - bin edge assignment
% - probablility matrix construction
% - diagnostic plots comparing the original data to generated sample

clear
close all

%% generate synthetic data to illustrate functionality; replace with desired dataset
rng(215);
% Define the length of each vector
n = 1000;
% Generate a base signal
base_signal = linspace(0, 10, n)'; % A simple linear ramp
% Generate noise
noise = 0.25 * randn(n,1);
% Create three correlated vectors
% Data vectors must be of identical length, but can contain NaNs
data1 = 10*(sin(base_signal) + noise);
data2 = -5*(-cos(base_signal) + noise);
data3  = 100*(exp(0.1 * base_signal) + noise);
% Replace a random 10% of each vector with NaNs
data1(randperm(n, round(0.1 * n))) = NaN;
data2(randperm(n, round(0.1 * n))) = NaN;
data3(randperm(n, round(0.1 * n))) = NaN;

%% parameters to define the sampling
% - nBins: number of bins over which to divide the range of input (one per
%          data element)
% - nSamples: How many output samples to produce
nBins = [30 30 30];
nSamples = 100000;

%% data pre-processing and setup
% place desired data vectors in 'var' structure, one element per dimension.
% Each vector should have the same length.
var.data1 = data1;
var.data2 = data2;
var.data3 = data3;

names = fieldnames(var);
nVar = numel(names);

% initializations
sampT = table();
varList = '';
edgeList = '';
trimList = '';

% determine bins, edges, and associated values for each input variable.
% This is governed by nBins.
for i = 1:nVar
    maxVal = max(var.(names{i}));
    minVal = min(var.(names{i}));
    range = maxVal - minVal;
    bin.(names{i}) = range / nBins(i);
    edges.(names{i}) = [0:nBins(i)]*bin.(names{i}) + minVal;
    vals.(names{i}) = edges.(names{i})(1:end-1) + bin.(names{i})/2;
    % slightly increase outer edge of last bin so all will fall in the n-1
    % bin of histcnd function
    edges.(names{i})(end) = edges.(names{i})(end) * 1.001;
    % accumulate variable and edge lists to feed to counting function
    varList = sprintf('%s,var.%s', varList, names{i});
    edgeList = sprintf('%s,edges.%s', edgeList, names{i});
    trimList = [trimList,',1:nBins'];
end

% trim leading ',' from these lists
varList(1) = [];
edgeList(1) = [];
trimList(1) = [];

%% compute random sample distributions
% generate histogram counts matrix
eval(['counts = createNDHistogram(' varList ',' edgeList ');'])
% for counts, get rid of each dimension's extraneous last bin
% (quirk of histcnd function)
eval(['counts = counts(' trimList ');'])

% generate probabililty matrix from counts
prob   = counts ./ numel(var.(names{1}));
% generate random samples via inverse transform method
samp = inverseTransformSampleND(prob, vals, nSamples);

% place output samples from structure to table
for kk = 1:nVar
    sampT.(names{kk}) = samp(:,kk);
end

%% generate diagnostic plots
% Plot 1D histograms of marginal distributions for each variable, orginal vs sample
figure( 'pos', [75 50 1500 750]);
for k = 1:nVar
    subplot(1, nVar, k);
    histogram(var.(names{k}),[-Inf edges.(names{k}) Inf],'Normalization','probability');
    histogram(var.(names{k}),[-Inf edges.(names{k}) Inf],'Normalization','probability');
    hold on;
    histogram(samp(:,k),[-Inf edges.(names{k}) Inf],'Normalization','probability');
    grid on;
    xlabel(names{k});
    ylabel('probability');
    legend('Source','Sample');
    hold off;
end
sgtitle(sprintf('1D Source (N=%d) vs Sample (N=%d) Comparison',...
    numel(var.(names{1})), nSamples),'FontWeight','bold');

% Plot 2D histograms for all 2-way combinations of input variables
% - find the unique 2-way combinations
comb = combvec([1:nVar], [1:nVar])';
% - remove combinations of idential elements
comb(comb(:,1)==comb(:,2),:) = [];
% - remove duplicate combinations (after sorting)
comb = sort(comb,2);
comb =  unique(comb, 'rows');

% plotting: top row is originals, bottom is samples, for each
% combination
[nPlotCols,~] = size(comb);
nPlot = 0;
figure( 'pos', [75 50 1500 750]);

for j = 1:nPlotCols
    subplot(2,nPlotCols,j);
    h1 = histogram2(var.(names{comb(j,1)}), var.(names{comb(j,2)}),...
        edges.(names{comb(j,1)}), edges.(names{comb(j,2)}),...
        'EdgeColor', [0 0 0], ...
        'Normalization', 'probability', 'DisplayStyle', 'tile');
    xlabel(names{comb(j,1)},'FontWeight','bold');
    ylabel(names{comb(j,2)},'FontWeight','bold');
    title(sprintf('Source N=%d', numel(var.(names{comb(j,1)}))));
    hc1 = colorbar;
    hc1.Label.String = 'Probability';
    ax1 = gca();

    subplot(2,nPlotCols,j+nPlotCols);
    h2 = histogram2(samp(:,comb(j,1)), samp(:,comb(j,2)),...
        edges.(names{comb(j,1)}), edges.(names{comb(j,2)}),...
        'EdgeColor', [0 0 0], ...
        'Normalization', 'probability', 'DisplayStyle', 'tile');
    xlabel(names{comb(j,1)},'FontWeight','bold');
    ylabel(names{comb(j,2)},'FontWeight','bold');
    title(sprintf('Sample N=%d', nSamples));
    hc2 = colorbar;
    hc2.Label.String = 'Probability';
    ax2 = gca();
    ax2.CLim = ax1.CLim;
end

sgtitle('2D Source vs Sample Comparisons', 'FontWeight', 'bold');

%% function to compute an N-dimensional histogram count of the input data
function histND = createNDHistogram(varargin)
    % Validate input: must have an even number of arguments
    if mod(length(varargin), 2) ~= 0
        error('Inputs must consist of pairs of data vectors and their associated edges.');
    end

    % Number of dimensions
    numDims = length(varargin) / 2;

    % Extract data and edges
    data = varargin(1:numDims);
    edges = varargin(numDims+1:end);

    % Validate consistency
    dataLengths = cellfun(@length, data);
    if any(dataLengths ~= dataLengths(1))
        error('All data vectors must have the same length.');
    end

    % Initialize histogram with zeros
    binSizes = cellfun(@(e) length(e)-1, edges);
    histND = zeros(binSizes);

    % Process each data point
    numDataPoints = dataLengths(1);
    for i = 1:numDataPoints
        % Find the bin indices for each dimension
        binIndices = zeros(1, numDims);
        valid = true;
        for dim = 1:numDims
            currentValue = data{dim}(i);
            binIndex = find(currentValue >= edges{dim}(1:end-1) & currentValue < edges{dim}(2:end), 1);

            % Special handling for the last edge
            if currentValue == edges{dim}(end)
                binIndex = binSizes(dim);
            end

            % If no valid bin is found, invalidate this point
            if isempty(binIndex)
                valid = false;
                break;
            else
                binIndices(dim) = binIndex;
            end
        end

        % Increment the histogram if the bin indices are valid
        if valid
            histND = incrementBin(histND, binIndices);
        end
    end
end

function histND = incrementBin(histND, binIndices)
    % Use dynamic indexing to increment the appropriate bin
    subs = num2cell(binIndices);
    histND(subs{:}) = histND(subs{:}) + 1;
end
