%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function out = extract_raw_data (in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (c) 2014-2019, Infineon Technologies AG
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification,are permitted provided that the
% following conditions are met:
%
% Redistributions of source code must retain the above copyright notice, this list of conditions and the following
% disclaimer.
%
% Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
% disclaimer in the documentation and/or other materials provided with the distribution.
%
% Neither the name of the copyright holders nor the names of its contributors may be used to endorse or promote
% products derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE  FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
% WHETHER IN CONTRACT, STRICT LIABILITY,OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% This simple example demos the acquisition of data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cleanup and init
% Before starting any kind of device the workspace must be cleared and the
% MATLAB Interface must be included into the code. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
disp('******************************************************************');
addpath('..\..\RadarSystemImplementation'); % add Matlab API
clear all %#ok<CLSCR>
close all
resetRS; % close and delete ports

% 1. Create radar system object
szPort = findRSPort; % scan all available ports
oRS = RadarSystem(szPort); % setup object and connect to board

disp('Connected RadarSystem:');
oRS %#ok<*NOPTS>

% 2. Enable automatic trigger with frame time 1s
% oRS.oEPRadarBase.set_automatic_frame_trigger(1000000);
% 
oRS.oEPRadarBase.stop_automatic_frame_trigger; % stop it to change values 
 
  oRS.oEPRadarFMCW.lower_frequency_kHz = 24025000; % lower FMCW frequency   
  oRS.oEPRadarFMCW.upper_frequency_kHz = 24225000; % upper FMCW frequency   
  oRS.oEPRadarFMCW.tx_power = oRS.oEPRadarBase.max_tx_power;   
  oRS.oEPRadarBase.num_chirps_per_frame = 16;   
  oRS.oEPRadarBase.num_samples_per_chirp = 128; % [32, 64, 128, 256]   
  oRS.oEPRadarBase.rx_mask = bin2dec('0011'); % enable two RX antennas   
  oRS.oEPRadarFMCW.direction = 'Up Only';

  
% Tc = cast(oRS.oEPRadarBase.chirp_duration_ns * (10^-9),'double'); % chirp time
oRS.oEPTargetDetection.min_range_cm = 90; % set max distance
oRS.oEPTargetDetection.max_range_cm = 300;
oRS.oEPTargetDetection.max_speed_kmh = 10; % set max speed
oRS.oEPTargetDetection.min_speed_kmh = 0;
i = 0;

chirp_time = oRS.oEPRadarBase.chirp_duration_ns
frame_interval = oRS.oEPRadarBase.min_frame_interval_us
min_rf_freq = oRS.oEPRadarBase.min_rf_frequency_kHz
max_rf_freq = oRS.oEPRadarBase.max_rf_frequency_kHz
% disp(chirp_time)
% disp(frame_interval)
% disp(min_rf_freq)
% disp(max_rf_freq)

% % Variables needed to write ydata to text file
% Change variables accordingly
chirps = 16; % chirps per frame
samples = 128; % samples per chirp
receivers = 1; % text write and read doesn't work for receivers = 2
frames = 600; % how many frames do you want recorded?
current_line = 0;
j = 1;

% % % % Variables needed for audio capability
            max_range = 3; % in meters
            range_interval = (max_range-0.9)/4;
            s1_range = max_range - (range_interval);
            s2_range = max_range - (2*range_interval);
            s4_range = max_range - (3*range_interval);
            s3_range = 0.9;


fileID1 = fopen('T1Single_Target_notesRx1.txt','w'); % Receiver 1 ydata
fileID2 = fopen('T1Single_Target_notesRx2.txt','w'); % Receiver 2 ydata

while true
    % 3. Trigger radar chirp and get the raw data
    [mxRawData, sInfo] = oRS.oEPRadarBase.get_frame_data;
    ydata = mxRawData; % get raw data

    disp(j)


 % % % % Write ydata matrix to text file (can only upload data from 1 receiver at a time)
    % % Write data from Rx1
    if j <= frames
        for i = 1:chirps
            fprintf(fileID1,'%f %f\n',[real(ydata(:,1,i)),imag(ydata(:,1,i))].'); 
        end
    else
        break
    end
    % % Write data from Rx2
    if j <= frames
        for i = 1:chirps
            fprintf(fileID2,'%f %f\n',[real(ydata(:,2,i)),imag(ydata(:,2,i))].'); 
        end
    else
        break
    end
    j = j+1;
    
    disp("target position")
    t_radius = (oRS.oEPTargetDetection.get_targets.radius)*(0.01)
    t_angle_az = oRS.oEPTargetDetection.get_targets.azimuth
    t_angle_elevation = oRS.oEPTargetDetection.get_targets.azimuth
    r_speed = oRS.oEPTargetDetection.get_targets.radial_speed
    a_speed = oRS.oEPTargetDetection.get_targets.azimuth_speed
    e_speed = oRS.oEPTargetDetection.get_targets.elevation_speed
    

    if (t_radius <= max_range) & (t_radius >= 0.9) % limits audio output when target is out of bounds
                if t_angle_az >= 0

                        if (t_radius <= max_range) & (t_radius > s1_range)
                            Fo = 349; % F
                        elseif t_radius <= s1_range & t_radius > s2_range
                            Fo = 392; % G
                        elseif t_radius <= s2_range & t_radius >= s4_range
                            Fo = 440; % A
                        elseif t_radius <= s4_range & t_radius >= s3_range
                            Fo = 494; % B
                        end
                                               
                elseif t_angle_az < 0 

                        if t_radius <= max_range & t_radius > s1_range
                            Fo = 349; % F
                        elseif t_radius <= s1_range & t_radius > s2_range
                            Fo = 330; % E
                        elseif t_radius <= s2_range & t_radius >= s4_range
                            Fo = 294; % D
                        elseif t_radius <= s4_range & t_radius >= s3_range
                            Fo = 262; % C
                        end

                end
%%%% Piano   
                  vector = 1./logspace(1,3,(4*7000)+1);
                  t = 0:1/7000:4;
                  signal = vector.*(5*(sin(2*pi*Fo*t))+(3*sin((2*pi*Fo*2*t)))+(.5*sin((2*pi*Fo*3*t)))+(0.05*sin(2*pi*Fo*4*t)));
                  soundsc(signal)
                  pause(5)
%%%% Guitar             
%                 Fs       = 44100;
%                 x = zeros(Fs*4, 1);
%                 delay = round(Fs/Fo);
%                 b  = firls(42, [0 1/delay 2/delay 1], [0 0 1 1]);
%                 a  = [1 zeros(1, delay) -0.5 -0.5];
% 
%                 zi = rand(max(length(b),length(a))-1,1);
%                 note = filter(b, a, x, zi);
%                 note = note-mean(note);
%                 note = note/max(abs(note));
%                 hplayer = audioplayer(note, Fs); play(hplayer)
%                 pause(3)

    end


end

% % % Conclude write to text file for receiver 1 and 2
fclose(fileID1); 
fclose(fileID2);

% disp(j)
% PhDiff_avg_array(1,j) = sum(PhDiff_array)/(chirps/2);
% vr_avg_array(1,j) = sum(vr_array)/(chirps/2);
% PhDiff_avg_array(1,j)
% vr_avg_array(1,j)