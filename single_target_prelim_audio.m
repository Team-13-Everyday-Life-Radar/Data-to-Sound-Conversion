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
oRS.oEPRadarADCXMC.samplerate_Hz = 640000; % 640000 affects the chirp time
  
% Tc = cast(oRS.oEPRadarBase.chirp_duration_ns * (10^-9),'double'); % chirp time
oRS.oEPTargetDetection.min_range_cm = 90; % set max distance
oRS.oEPTargetDetection.max_range_cm = 300;
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
frames = 2000; % how many frames do you want recorded?
current_line = 0;
j = 1;


% 
c = 3e8;
snum = double(oRS.oEPRadarBase.num_samples_per_chirp); % sample per chirp
cnum = double(oRS.oEPRadarBase.num_chirps_per_frame); % chirp per frame
fs = double(oRS.oEPRadarADCXMC.samplerate_Hz);  % samppling freq
bw = double(oRS.oEPRadarFMCW.bandwidth_per_second); % chcirpslope in microsec
ts =  double(oRS.oEPRadarBase.chirp_duration_ns); % chirp duration(ns) 
t =  0 : 1/fs : (1e-9*ts)-(1/fs); 
f1= double(oRS.oEPRadarFMCW.lower_frequency_kHz)*1e3; % min freq 
f2= double(oRS.oEPRadarFMCW.upper_frequency_kHz)*1e3; % max freq
f= f1: (f2-f1)/snum: f2-((f2-f1)/snum); 
chirpslope = (f2-f1)/(t(end));

threshhold = 15;

minrange = 0.9; maxrange= 3;
minbeatfreq = ((chirpslope)*2*minrange)/c ;
maxbeatfreq = ((chirpslope)*2*maxrange)/c ;
% % % % Variables needed for audio capability
            maxrange = 3; % in meters
            range_interval = (maxrange-minrange)/4;
            s1_range = maxrange - (range_interval);
            s2_range = maxrange - (2*range_interval);
            s4_range = maxrange - (3*range_interval);
            s3_range = minrange;
            instrument = "piano";
            angle_buffer = 5;
            range_buffer = 0.2;
            prev_range = 0;

fC = (f2+f1)/2;
lambda = c / fC;
antenna_spacing = 6.22e-3; % in meters

fileID1 = fopen('T2Single_Target_notesRx1.txt','w'); % Receiver 1 ydata
fileID2 = fopen('T2Single_Target_notesRx2.txt','w'); % Receiver 2 ydata
 t_radius = 0;
while true
    % 3. Trigger radar chirp and get the raw data
    [mxRawData, sInfo] = oRS.oEPRadarBase.get_frame_data;
    ydata = mxRawData; % get raw data
    
    disp(j)

% % % % % % % % % % Angle Detection Algorithm
  %     Amp1(:,1,(1:cnum))= sqrt(real(ydata(:,1,(1:cnum))).^2 + imag(ydata(:,1,(1:cnum))).^2); 
%     Amp2(:,1,(1:cnum))= sqrt(real(ydata(:,2,(1:cnum))).^2 + imag(ydata(:,2,(1:cnum))).^2); 
    %angle= atan(imag(ydata(:,1))./ real(ydata(:,1)));
    
    Amp1(:,1,(1:cnum))= ydata(:,1,(1:cnum)); 
    Amp2(:,1,(1:cnum))= ydata(:,2,(1:cnum)); 
    nfft = 2048;
    
    Af1(:,1,(1:cnum)) = abs(fftshift(fft(Amp1(:,1,(1:cnum)),nfft)));
    Af2(:,1,(1:cnum)) = abs(fftshift(fft(Amp2(:,1,(1:cnum)),nfft)));
    
    Pf1(:,1,(1:cnum)) = (fftshift(fft(Amp1(:,1,(1:cnum)),nfft)));
    Pf2(:,1,(1:cnum)) = (fftshift(fft(Amp2(:,1,(1:cnum)),nfft)));
      
    freq = (fs/nfft)*(-nfft/2:nfft/2 -1);
    fr = freq((nfft/2 +1):(end));
    
   y1 = mean(Af1,3);
   y2 = mean(Af2,3);
   
   x1 = mean(Pf1,3);
   x2 = mean(Pf2,3);

   y1_half = y1((nfft/2 +1):end);
   y2_half = y2((nfft/2 +1):end);
   
   x1_half = x1((nfft/2 +1):end);
   x2_half = x2((nfft/2 +1):end);
   
       stepfreq = fr(2) - fr(1); % affected by fs and nfft
       minindx = round((minbeatfreq - fr(1))/stepfreq);
       maxindx = round((maxbeatfreq - fr(1))/stepfreq); % rounded up or down ?? 
       
    frequ = fr((minindx+1):(maxindx+1));
    y1_final = y1_half((minindx+1):(maxindx+1));
    y2_final = y2_half((minindx+1):(maxindx+1));
    
    x1_final = x1_half((minindx+1):(maxindx+1));
    x2_final = x2_half((minindx+1):(maxindx+1));
     
%  [p1 , ind1] = findpeaks(y1_final); %pass the averaged spectrum of all chirps
%     
%     for i= 1 : length(p1)
%         if (p1(i) > threshhold)
%             pks1(i) = p1(i);
%             df1(i)= frequ(ind1(i)); %df(i) = (ind(i)-1).*(freq(2)-freq(1));
%             target1(i)= round((c*df1(i))/(2*(chirpslope)),2,'significant');
%         end
%       
%     end
%      
%      [p2 , ind2] = findpeaks(y2_final); %pass the averaged spectrum of all chirps
%     
%     for i= 1 : length(p2)
%         if (p2(i) > threshhold)
%             pks2(i) = p2(i);
%             df2(i)= frequ(ind2(i)); %df(i) = (ind(i)-1).*(freq(2)-freq(1));
%             target2(i)= round((c*df2(i))/(2*(chirpslope)),2,'significant');
%         end
%       
%      end
%     
%    if (~exist('pks1') || ~exist('df1') || ~exist('target1'))
%            pks1 = 0;
%            df1 = 0 ; 
%            target1 = 0;
%         disp(' No Object Detected ')
%    else
           
           [val_max_1, ind_max_1] = max(y1_final);
           [val_max_2, ind_max_2] = max(y2_final);

            ang_rx_1 = angle(x1_final(ind_max_1));
            ang_rx_2 = angle(x2_final(ind_max_2));

            d_phi = ang_rx_1 - ang_rx_2;

            
            if (d_phi <= 0)
                d_phi = d_phi + 2*pi;
            end
             d_phi = d_phi - pi;

            target_angle = (asin((d_phi * lambda) / (antenna_spacing * (2*pi)))); % AOA in radians

            target_angle_deg = ((target_angle) * 180 / pi); % AOA in degrees


    t_angle_az = target_angle_deg
% % % % % % % % % % END Angle Detection Algorithm
    
% % % % % % % % % % Range Detection Algorithm
    Amp_range(:,1,(1:cnum))= sqrt(real(ydata(:,1,(1:cnum))).^2 + imag(ydata(:,1,(1:cnum))).^2); 
    Amp2_range(:,1,(1:cnum))= sqrt(real(ydata(:,2,(1:cnum))).^2 + imag(ydata(:,2,(1:cnum))).^2); 

    
    nfft = 2048;
    Af_range(:,1,(1:cnum)) = abs(fftshift(fft(Amp_range(:,1,(1:cnum)),nfft)));
    Af2_range(:,1,(1:cnum)) = abs(fftshift(fft(Amp2_range(:,1,(1:cnum)),nfft)));
    
%     freq = (fs/nfft)*(-nfft/2:nfft/2 -1);
%     fr = freq((nfft/2 +1):(end));
    
   y1_range = mean(Af_range,3);
   y2_range = mean(Af2_range,3);
   y_range = 0.5*(y1_range+y2_range);
  

   yhalf_range = y_range((nfft/2 +1):end);
   

%        stepfreq = fr(2) - fr(1); % affected by fs and nfft
%        minindx = round((minbeatfreq - fr(1))/stepfreq);
%        maxindx = round((maxbeatfreq - fr(1))/stepfreq); % rounded up or down ?? 
%        
%     frequ = fr((minindx+1):(maxindx+1));
    yfinal_range = yhalf_range((minindx+1):(maxindx+1));
%    

 [p , ind] = findpeaks(yfinal_range); %pass the averaged spectrum of all chirps
    
% p = p.*(p >= threshhold);
 
    for i= 1 : length(p)
        if p(i) >= threshhold 
            pks(i) = p(i);
            df(i)= frequ(ind(i)); %df(i) = (ind(i)-1).*(freq(2)-freq(1));
            target(i)= round((c*df(i))/(2*(chirpslope)),4,'significant');
           % target(i)= ceil(tgt(i)*4)/4;
           
        end
      
    end
    
   if (~exist('pks') || ~exist('df') || ~exist('target'))
           pks = 0;
          df = 0 ; 
          % target = 0;
        disp(' No Object Detected ')
        t_radius = 0;
   else
      
    rang = (c/(2*chirpslope)).*frequ;
    disp(max(target))
    tgt_range = target;
           t_radius = max(tgt_range) % current range
           
      clear pks; 
     clear target ;
    clear df ;
   end
% % % % % % % % END Range Detection Algorithm
    
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
    
    

    if (t_radius <= maxrange) & (t_radius > 1.07) % limits audio output when target is out of bounds
          if abs(t_radius-prev_range) <= range_buffer      
        if t_angle_az >= angle_buffer
angle_exist = 1;
                        if (t_radius <= maxrange) & (t_radius > s1_range)
                            Fo = 349; % F
                        elseif t_radius <= s1_range & t_radius > s2_range
                            Fo = 392; % G
                        elseif t_radius <= s2_range & t_radius >= s4_range
                            Fo = 440; % A
                        elseif t_radius <= s4_range & t_radius >= s3_range
                            Fo = 494; % B
                        end
                                               
                elseif t_angle_az <= -angle_buffer
angle_exist = 1;
                        if t_radius <= maxrange & t_radius > s1_range
                            Fo = 349; % F
                        elseif t_radius <= s1_range & t_radius > s2_range
                            Fo = 330; % E
                        elseif t_radius <= s2_range & t_radius >= s4_range
                            Fo = 294; % D
                        elseif t_radius <= s4_range & t_radius >= s3_range
                            Fo = 262; % C
                        end
        else
            angle_exist = 0
                end
                if instrument == "piano" && angle_exist == 1
%%%% Piano   
                  vector = 1./logspace(1,3,(4*7000)+1);
                  t = 0:1/7000:4;
                  signal = vector.*(5*(sin(2*pi*Fo*t))+(3*sin((2*pi*Fo*2*t)))+(.5*sin((2*pi*Fo*3*t)))+(0.05*sin(2*pi*Fo*4*t)));
                  soundsc(signal)
                  pause(5)
                elseif instrument == "guitar" && angle_exist == 1
%%%% Guitar             
                Fs       = 44100;
                x = zeros(Fs*4, 1);
                delay = round(Fs/Fo);
                b  = firls(42, [0 1/delay 2/delay 1], [0 0 1 1]);
                a  = [1 zeros(1, delay) -0.5 -0.5];

                zi = rand(max(length(b),length(a))-1,1);
                note = filter(b, a, x, zi);
                note = note-mean(note);
                note = note/max(abs(note));
                hplayer = audioplayer(note, Fs); play(hplayer)
                pause(3)

                end
          else
          end
    end

prev_range = t_radius;
end

% % % Conclude write to text file for receiver 1 and 2
fclose(fileID1); 
fclose(fileID2);

% disp(j)
% PhDiff_avg_array(1,j) = sum(PhDiff_array)/(chirps/2);
% vr_avg_array(1,j) = sum(vr_array)/(chirps/2);
% PhDiff_avg_array(1,j)
% vr_avg_array(1,j)
