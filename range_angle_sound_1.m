
    % Trigger radar chirp and get the raw data
    [mxRawData, sInfo] = oRS.oEPRadarBase.get_frame_data;
    ydata = mxRawData; % get raw data
    
% % % % % % % % % % Angle Detection Algorithm
    
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
     
    [p1 , ind1] = findpeaks(y1_final); %pass the spectrum of all chirps in RX1
    
    for i= 1 : length(p1)
        if (p1(i) > threshhold)
            pks1(i) = p1(i);
        end
        if (~exist('pks1'))
 
        disp(' No Object Detected ')
        else
            ang_rx_1 = angle(max(x1_final));
            ang_rx_2 = angle(max(x2_final));
            d_phi = ang_rx_1 - ang_rx_2;
            
            if (d_phi <= 0)
                d_phi = d_phi + 2*pi;
            end
             d_phi = d_phi - pi;

            target_angle = (asin((d_phi * lambda) / (antenna_spacing * (2*pi)))); % AOA in radians

            target_angle_deg = ((target_angle) * 180 / pi); % AOA in degrees
        end
    end
     
    t_angle_az = target_angle_deg
% % % % % % % % % % END Angle Detection Algorithm
    
% % % % % % % % % % Range Detection Algorithm
    Amp_range(:,1,(1:cnum))= sqrt(real(ydata(:,1,(1:cnum))).^2 + imag(ydata(:,1,(1:cnum))).^2); 
    Amp2_range(:,1,(1:cnum))= sqrt(real(ydata(:,2,(1:cnum))).^2 + imag(ydata(:,2,(1:cnum))).^2); 

    
    nfft = 2048;
    Af_range(:,1,(1:cnum)) = abs(fftshift(fft(Amp_range(:,1,(1:cnum)),nfft)));
    Af2_range(:,1,(1:cnum)) = abs(fftshift(fft(Amp2_range(:,1,(1:cnum)),nfft)));
    
    freq = (fs/nfft)*(-nfft/2:nfft/2 -1);
    fr = freq((nfft/2 +1):(end));
    
   y1_range = mean(Af_range,3);
   y2_range = mean(Af2_range,3);
   y_range = 0.5*(y1_range+y2_range);
  

   yhalf_range = y_range((nfft/2 +1):end);
   

       stepfreq = fr(2) - fr(1); % affected by fs and nfft
       minindx = round((minbeatfreq - fr(1))/stepfreq);
       maxindx = round((maxbeatfreq - fr(1))/stepfreq); % rounded up or down ?? 
%        
    frequ = fr((minindx+1):(maxindx+1));
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
    max_target = max(target)
    tgt_range = target;
    t_radius = max(tgt_range); % current range
           
    clear pks; 
    clear target ;
    clear df ;
   end
% % % % % % % % END Range Detection Algorithm
    
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
