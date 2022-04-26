function [P,f]=spec1(x,window,overlap,fs);
% Estimate the PSD for the signal x, using FFTs of length(window), and the overlap (0-1)
% fs is the sampling rate. 
Nfft=length(window);
N=length(x);
x=reshape(x,N,1);
window=reshape(window,Nfft,1);
P=nan;
f=nan;
shift=round((1-overlap)*Nfft);
if (shift>0) & (Nfft<N)
    y=[];
    weight=[];
    weight0=sum(window.^2);
    weight1=sum(window.^2);
    i_start=1;
    i_end=Nfft;
    no_windows=0;
    while i_end<=N
        x_temp=x(i_start:i_end).*window;
        weight=[weight,weight0];
        y=[y,x_temp];
        no_windows=no_windows+1;
        i_start=i_start+shift;
        i_end=i_end+shift;
    end          
    Y=fft(y)';
    Y2=abs(Y).^2;
    if no_windows>1
        P=mean(Y2);
    else
        P=Y2;
    end
    P=P/weight0/fs;
    if mod(Nfft,2)==0
        P(2:Nfft/2+1)=P(2:Nfft/2+1)+P(Nfft:-1:Nfft/2+1);
        P=P(1:Nfft/2+1);
    else
        P(2:(Nfft-1)/2+1)=P(2:(Nfft+1)/2)+P(Nfft:-1:(Nfft+1)/2+1);
        P=P(1:(Nfft+1)/2);
    end
    P(1)=2*P(1);
    f=[0:length(P)-1]/Nfft*fs;
end

% keyboard

