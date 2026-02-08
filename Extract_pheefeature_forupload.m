%{
This script was written by Junfeng Huang.
Originally created in 2021-09.
Last modified on 2025-03-04.

Purpose:
This script extracts acoustic features of Phee calls from audio recordings
based on manually annotated labels. It supports batch processing across
multiple groups and subjects, and computes both single-call and multi-pip
Phee statistics.

The script is intended for research use and public release.
%}

clear;
close all;
clc;

currentpath = cd;
cd(currentpath);

%% Load data and configuration

% Read batch configuration table (starting from cell B2)
[~,~,Sheet_table] = xlsread('批量提取pheefeature参数excel表_20250717.xlsx','B2:D7');

global call_count Plot_or_not Save_figure_name
call_count = 0;
Plot_or_not = 0;   % If 1, plot and save spectrograms and F0 curves

jmarmoset = 0;

%% Extract Phee features for all groups
group_number = size(Sheet_table,1);
mkdir Phee_Feature_All_Groups
cd('Phee_Feature_All_Groups');

for group_index = 1:group_number

    label_folder = Sheet_table{group_index,1};   % Path to label files
    audio_folder = Sheet_table{group_index,2};   % Path to audio files
    group_file   = Sheet_table{group_index,3};   % Group folder name

    label_file = dir(fullfile(label_folder,'*.txt'));
    audio_file = dir(fullfile(audio_folder,'*.wav'));

    begin_file = 1;
    marmoset_number = length(audio_file);

    mkdir(group_file);
    cd(group_file);

    %% Process each subject/session
    for session_index = begin_file:marmoset_number

        jmarmoset = jmarmoset + 1;
        Name{jmarmoset} = audio_file(session_index).name;

        disp(['Group ',group_file, ...
              ': processing ',num2str(session_index), ...
              ' / ',num2str(marmoset_number)]);
        disp(audio_file(session_index).name);

        audio_file_path = fullfile(audio_folder,audio_file(session_index).name);
        marmoset_file = audio_file(session_index).name(1:end-4);

        mkdir(marmoset_file);
        cd(marmoset_file);

        %% Match label file
        label_flag = 0;
        for ifile = 1:length(label_file)
            if strcmp(marmoset_file,label_file(ifile).name(1:end-4))
                label_flag = 1;
                break;
            end
        end

        label_file_path = fullfile(label_folder,[marmoset_file,'.txt']);

        if label_flag == 0
            warning(['No label file found for ',audio_file(session_index).name]);
            data_name = [marmoset_file,'_pheedata'];
            having_call = 0;
            save(data_name,'having_call');
            cd ..
            continue;
        end

        %% Load label file and keep only Phee calls
        call_data = Inset_loading_label_data_HJF(label_file_path);
        non_phee = call_data.type ~= 1;

        call_data.number(non_phee)     = [];
        call_data.type(non_phee)       = [];
        call_data.start_time(non_phee) = [];
        call_data.end_time(non_phee)   = [];

        %% Load audio
        [audio_data_origin,fs_origin] = audioread(audio_file_path);
        fs = 48000;

        if fs_origin ~= fs
            audio_data = resample(audio_data_origin,fs,fs_origin);
        else
            audio_data = audio_data_origin;
        end
        clear audio_data_origin

        %% Extract features for all Phee calls
        Save_figure_name = marmoset_file;
        [feature,call_data] = Inset_AllPheeCall_feature_HJF( ...
            call_data,audio_data,fs,audio_file(session_index).name,Human_modified_call);

        %% Calculate n-pip Phee statistics
        interval_threshold = 1;   % seconds
        nPhee_call = Inset_calculate_npip_Phee_HJF(call_data,interval_threshold);
        nPhee_feature = Inset_Calculate_npip_Phee_feature(nPhee_call,feature);

        audio_duration = length(audio_data)/fs;
        clear audio_data call_data

        %% Save results
        data_name = [marmoset_file,'_pheedata'];
        having_call = 1;
        save(data_name);

        cd ..
    end
    cd ..
end
cd ..

%% Collect all MAT files into one folder
current_file = cd;
mkdir('All_mat');

for group_index = 1:group_number
    Group_folder = fullfile(current_file,'Phee_Feature_All_Groups',Sheet_table{group_index,3});
    New_group_folder = fullfile(current_file,'All_mat',Sheet_table{group_index,3});
    mkdir(New_group_folder);

    marmoset_folder = dir(Group_folder);
    for i = 3:length(marmoset_folder)
        if marmoset_folder(i).isdir
            cd(fullfile(Group_folder,marmoset_folder(i).name));
            copyfile('*pheedata.mat',New_group_folder);
        end
    end
end
cd(current_file);

%% ========================= Function 1 =========================
function call_data = Inset_loading_label_data_HJF(filename)
% Load call label information from a Praat-exported text file.
% Labeling convention:
%   1: Phee call

warning off  % Disable warning messages

% Read label file exported from Praat (TextGrid converted to text)
data = readtable(filename,'Delimiter','=');

% Note:
% The first 12 rows contain header information and are ignored.
data_call_start_index = 13;

%% Find call onset times (xmin)
call_num = 0;
for i = data_call_start_index:length(data.FileType)
    if isequal(data.FileType{i},'xmin')
        call_num = call_num + 1;
        call_data.start_index(call_num) = i;
        call_data.start_time(call_num)  = data.ooTextFile(i);
    end
end

%% Find call offset times (xmax)
call_num = 0;
for i = data_call_start_index:length(data.FileType)
    if isequal(data.FileType{i},'xmax')
        call_num = call_num + 1;
        call_data.end_index(call_num) = i;
        call_data.end_time(call_num)  = data.ooTextFile(i);
    end
end

%% Find call labels
call_num = 0;
for i = data_call_start_index:length(data.FileType)
    if isequal(data.FileType{i},'text')
        call_num = call_num + 1;
        call_data.label_index(call_num) = i;
        call_data.label_value(call_num) = data.ooTextFile(i);
    end
end

%% Consistency check
if ~isequal(length(call_data.start_index), ...
            length(call_data.end_index), ...
            length(call_data.label_index))
    error('HJF: Inconsistent number of starts, ends, and labels.')
end

%% Simplify call_data structure
for i = 1:length(call_data.start_index)
    call_data_temp.number(i)     = i;
    call_data_temp.type(i)       = call_data.label_value(i);
    call_data_temp.start_time(i) = call_data.start_time(i);
    call_data_temp.end_time(i)   = call_data.end_time(i);
end

call_data = call_data_temp;
disp('Nice! Loading call data completed.');

clearvars -except call_data
end


%% ========================= Function 2 =========================
function nPhee_call = Inset_calculate_npip_Phee_HJF(call_data, interval_threshold)
% Identify n-pip Phee calls based on inter-call interval threshold (seconds)

Phee_call_data = call_data;
Phee_pip_index = find(Phee_call_data.type == 1);

if isempty(Phee_pip_index)
    disp('No Phee calls detected.');
    return
end

Phee_start_index(1) = Phee_pip_index(1);
Phee_num = 1;

for i = 1:length(Phee_pip_index)-1
    if Phee_call_data.start_time(Phee_pip_index(i+1)) - ...
       Phee_call_data.end_time(Phee_pip_index(i)) > interval_threshold
        Phee_end_index(Phee_num)   = Phee_pip_index(i);
        Phee_start_index(Phee_num+1) = Phee_pip_index(i+1);
        Phee_num = Phee_num + 1;
    end
end

Phee_end_index(Phee_num) = Phee_pip_index(end);

for i = 1:Phee_num
    nPhee(i) = Phee_end_index(i) - Phee_start_index(i) + 1;
end

%% Count different n-pip Phee calls
nPhee_call.num_Phee_1 = sum(nPhee == 1);
nPhee_call.num_Phee_2 = sum(nPhee == 2);
nPhee_call.num_Phee_3 = sum(nPhee == 3);
nPhee_call.num_Phee_4 = sum(nPhee == 4);
nPhee_call.num_Phee_Bigger_than_4 = sum(nPhee > 4);

%% Index of first pip of each n-pip Phee
nPhee_call.index_Phee_1 = Phee_start_index(nPhee == 1);
nPhee_call.index_Phee_2 = Phee_start_index(nPhee == 2);
nPhee_call.index_Phee_3 = Phee_start_index(nPhee == 3);
nPhee_call.index_Phee_4 = Phee_start_index(nPhee == 4);
nPhee_call.index_Phee_Bigger_than_4 = Phee_start_index(nPhee > 4);
end


%% ========================= Function 3 =========================
function [feature, call_data] = Inset_AllPheeCall_feature_HJF( ...
    call_data, audio_data, fs, audio_name, Human_modified_call)
% Extract acoustic features for all Phee calls in one audio file

global Phee_index
j = 0;
h_waitbar = waitbar(0,'Calculating call features...');

tmp = audio_name(1:end-4);
tmp_index = [];

for i = 1:length(Human_modified_call)
    if strcmp(tmp, Human_modified_call(i).fileName)
        tmp_index = [tmp_index,i];
    end
end

for i = 1:length(call_data.number)
    if floor(fs * call_data.end_time(i)) > length(audio_data)
        break
    end

    j = j + 1;
    Phee_index = j;
    waitbar(j / length(call_data.number));

    single_call_signal = audio_data( ...
        max(ceil(fs * call_data.start_time(i)),1) : ...
        floor(fs * call_data.end_time(i)));

    F0_corrected_index = [];
    for k = 1:length(tmp_index)
        if Phee_index == Human_modified_call(tmp_index(k)).call_index
            F0_corrected_index = tmp_index(k);
        end
    end

    [single_call_feature, single_F0_signal, new_time, mfcc_feature] = ...
        Inset_extract_single_Phee_feature_HJF( ...
        single_call_signal, fs, F0_corrected_index, Human_modified_call);

    feature.mfcc{Phee_index}        = mfcc_feature;
    feature.F0{Phee_index}          = single_F0_signal;
    feature.call_feature(Phee_index,:) = single_call_feature;
    feature.call_type(Phee_index)   = call_data.type(Phee_index);

    call_data.start_time(Phee_index) = ...
        call_data.start_time(Phee_index) + (new_time(1)-1)/fs;
    call_data.end_time(Phee_index) = ...
        call_data.start_time(Phee_index) + new_time(2)/fs;

    feature.start_time(Phee_index) = call_data.start_time(Phee_index);
    feature.end_time(Phee_index)   = call_data.end_time(Phee_index);
end

close(h_waitbar);
end

%% ========================= Function 4 (Optimized) =========================
function [call_feature,F0_signal,new_time] = ...
    Inset_extract_single_Phee_feature_HJF(call_signal,fs,F0_corrected_index,Human_modified_call)
% Optimized version:
% Only compute call_duration, call_peak_frequency, and F0-related features.

%% Ensure column vector
if size(call_signal,2) > 1
    call_signal = call_signal(:);
end

%% Band-pass filtering for F0 & peak frequency
filter_lf = 4000;
filter_hf = 12000;

signal_filt = Inset_Butterbp_hjf(call_signal,filter_lf,filter_hf,fs);
signal_filt = Inset_Butterbp_hjf(flipud(signal_filt),filter_lf,filter_hf,fs);
signal_filt = flipud(signal_filt);

%% Trim silence (same logic, but applied once)
wlen = round(fs*0.03);
inc  = round(fs*0.01);
win  = hanning(wlen);

frames = enframe(signal_filt,win,inc)';
frame_amp = sum(abs(frames),1);
threshold = mean(frame_amp)/20;

idx = find(frame_amp > threshold);
if isempty(idx)
    new_time = [1,length(signal_filt)];
else
    start_frame = idx(1);
    end_frame   = idx(end);
    frameTime = Inset_frame2time(size(frames,2),wlen,inc,fs);
    new_start = floor(frameTime(start_frame)*fs);
    new_end   = min(floor(frameTime(end_frame)*fs),length(signal_filt));
    signal_filt = signal_filt(new_start:new_end);
    new_time = [new_start,new_end];
end

%% Call duration
call_duration = length(signal_filt)/fs;

%% Peak frequency (FFT-based)
N = length(signal_filt);
Y = abs(fft(signal_filt));
Y = Y(1:floor(N/2)+1);
freq = (0:length(Y)-1)*fs/N;
[~,idx_peak] = max(Y);
call_peak_frequency = freq(idx_peak);

%% F0 extraction (spectrogram ridge)
[s,f,t] = spectrogram(signal_filt,hanning(wlen),wlen-inc,wlen,fs);
P = abs(s);
[~,idx] = max(P,[],1);
DF = f(idx);
DF = smooth(DF,0.05,'loess');

%% Manual F0 correction if available
if ~isempty(F0_corrected_index)
    DF = Human_modified_call(F0_corrected_index).F0_data;
end

F0_signal = DF;

%% F0 features
F0_nanmean = nanmean(DF);
F0_Min  = min(DF);
F0_Max  = max(DF);
F0_BW   = F0_Max - F0_Min;
F0_Std  = std(DF);
F0_Start = DF(1);
F0_End   = DF(end);
F0_Slope_All = (F0_End - F0_Start)/(t(end)-t(1));

%% Collect selected features ONLY
call_feature = [ ...
    call_duration, ...
    call_peak_frequency, ...
    F0_nanmean, ...
    F0_Min, ...
    F0_Max, ...
    F0_BW, ...
    F0_Std, ...
    F0_Start, ...
    F0_End, ...
    F0_Slope_All ...
    ];
end




%% ========================= Function 5 =========================
function [fmaxb,fminb,fmax,peak_freq_magratio] = ...
    Inset_Phee_dominant_frequency(call_signal,fs)
% Extract dominant frequency and bandwidth from FFT spectrum.

boundary_threshold = 0.05;

wlen = length(call_signal);
W2   = floor(wlen/2) + 1;
Y    = fft(call_signal);
Y_mag= abs(Y(1:W2));
freq = (0:W2-1)*fs/wlen;

Y_power_filtered = smooth(Y_mag,ceil(wlen/200));
[~,fmax_index] = max(Y_power_filtered);
fmax = freq(fmax_index);
peak_freq_magratio = Y_power_filtered(fmax_index)/sum(Y_power_filtered);

bw_idx = find(Y_power_filtered/max(Y_power_filtered) < boundary_threshold);
bw_idx = sort(bw_idx);

if length(bw_idx) > 2
    fmaxb = freq(bw_idx(find(bw_idx>fmax_index,1)));
    fminb = freq(bw_idx(find(bw_idx<fmax_index,1,'last')));
else
    error('Spectrum boundary not detected correctly.');
end
end


%% ========================= Function 6 =========================
function nPhee_feature = Inset_Calculate_npip_Phee_feature(nPhee_call,feature)
% Compute statistics for n-pip Phee calls (n = 1–4).

nPhee_feature.Phee_1   = nPhee_call.num_Phee_1;
nPhee_feature.Phee_1_1 = nanmean(feature.call_feature(nPhee_call.index_Phee_1,:),1);

nPhee_feature.Phee_2(1) = nPhee_call.num_Phee_2;
nPhee_feature.Phee_2(2) = nanmean(feature.end_time(nPhee_call.index_Phee_2+1) - ...
                                  feature.start_time(nPhee_call.index_Phee_2));
nPhee_feature.Phee_2(3) = nanmean(feature.start_time(nPhee_call.index_Phee_2+1) - ...
                                  feature.end_time(nPhee_call.index_Phee_2));
nPhee_feature.Phee_2_1  = nanmean(feature.call_feature(nPhee_call.index_Phee_2,:),1);
nPhee_feature.Phee_2_2  = nanmean(feature.call_feature(nPhee_call.index_Phee_2+1,:),1);
end


%% ========================= Function 7 =========================
function frameTime = Inset_frame2time(frameNum, framelen, inc, fs)
% Convert frame index to time (center of each frame)
frameTime = (((1:frameNum)-1)*inc + framelen/2) / fs;
end


%% ========================= Function 8 =========================
function signal_filtered = Inset_Butterbp_hjf(signal, Flow, Fhigh, Fs)
% Apply a 3rd-order Butterworth band-pass filter
[b,a] = butter(3,[Flow/(Fs/2), Fhigh/(Fs/2)],'bandpass');
signal_filtered = filter(b,a,signal);
end
