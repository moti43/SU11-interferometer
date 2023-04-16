% Simulation of SU11

function EoutS=Simulation_SU11()
    w_EOM=10e9;   %frequency of the EOM
    w_laser=100.001e6; %periodicity of the laser
    
    Pw=193.4e12;  %cental frequency of the pump
    Ptau=0.5e-12; %temporal pump width
    
    Sw=190.95e12; %central frequency of the signal
    Stau=1e-12;   %temporal signal width
    
    n=0; %just a boring index
    EoutS=zeros(100001,1000); %all the signal output
    EoutI=EoutS;              %all the idler output
    for t0=0:1/w_laser:1/w_laser*1000
        n=n+1;
        t=(-50e-12:1e-15:50e-12);
        % setting the input signal and pump pulses
        P=exp(-(t/Ptau).^2).*exp(1i*Pw*t);
        S=exp(-(t/Stau).^2).*exp(1i*Sw*t);
        
        % pump dispersion
        P=dispersion(P,18e-3*0.4);
        P2=P;
        
        % pump EOM
        P1=P.*exp(1i*cos(2*pi*w_EOM*3*(t+t0)));
                
        % signal dispersion
        S0=dispersion(S,18e-3*0.2);
        
        % FWM interaction
        I1=FWM(P1,P1,conj(S0));
        S1=S0;
        
        
        % some phase which I can add to one arm 
        %I1=I1.*exp(1i*cos(2*pi*w_EOM*10*(t+t0)));
        
        
        % FWM interaction
        I2=FWM(P2,P2,conj(S1));
        S2=FWM(P2,P2,conj(I1));
        
        % Normalizing all the fields (asuming high nonlinear efficiency)
        S1=S1/max(S1);
        I1=I1/max(I1);
        S2=S2/max(S2);
        I2=I2/max(I2);
        
        % Combining the arms of the interferometer
        S3=S2+S1;
        I3=I2+I1;
        
        %More dispersion for after the time-lens
        S4=dispersion(S3,-18e-3*0.2);
        I4=dispersion(I3,-18e-3*0.2);
        
        % Time stretch - We dont need to simulate the time-stretch since we
        % we are not limited by the detector resolution.
%        S4=dispersion(S3,140e-3);
%        I4=dispersion(I3,140e-3);
        
        %For showing each pulse
        if 1==0
            figure
            hold on
            plot(abs(S1));
            plot(abs(S2));
            plot(abs(S3));
            plot(abs(S4));
        end
        EoutI(:,n)=I4;
        EoutS(:,n)=S4;
    end
    figure
    imagesc(abs(EoutS));
    figure
    imagesc(abs(EoutI));
    EoutSS=zeros(size(EoutS));
    for n=1:1:1001
        [val,ind]=max(EoutS(:,n));
        inds=floor(length(EoutSS)/2);
        EoutSS(inds-3000:inds+3000,n)=EoutS(ind-3000:ind+3000,n);
    end
end


function d=dispersion(E,D)
    w=(1:1:length(E))-1;
    
    Ef=fftshift(fft(E));
    [w_max,w_ind]=max(abs(Ef));
    
    w=w-w_ind;
    w=w.*1e10;
    
    Ef=Ef.*exp(1i.*D.*w.^2/4e19);
    d=ifft(Ef);
end

function s=addPhase(E)
    w=(1:1:length(E))-1;
    
    Ef=fftshift(fft(E));
    [w_max,w_ind]=max(abs(Ef));
    
    w=w-w_ind;
    w=w.*1e10;
    
    Ef=Ef.*exp(1i.*w/4e10);
    s=ifft(Ef);
end


function E4=FWM(E1,E2,E3)
    Chi3=0.1;
    E4=Chi3.*E1.*E2.*E3;
end




