%%
clear all
close all

% Specify dataset
filepath = './supermag_data/csv/2010-April-06.csv';

% Load data
[jd2000, IAGA, MLT, MLAT, IGRF_DECL, SZA, N, E, Z] = supermagdataall2matrix(filepath);

% IAGA variable contains relative indexed information for all available
% stations (alphabetized). Each station index corresponds to a column in any of the data
% vectors.

