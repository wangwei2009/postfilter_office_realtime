function [ DS, x1] = DelaySumURA( x,fs,N,frameLength,inc,r,angle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%frequency-domain delay-sum beamformer using circular array
%   
%      input :
%          x : input signal ,samples * channel
%          fs: sample rate
%          N : fft length,frequency bin number
%frameLength : frame length,usually same as N
%        inc : step increment
%          r : array element radius
%      angle : incident angle
%
%     output :
%         DS : delay-sum output
%         x1 : presteered signal,same size as x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 340;
Nele = size(x,2);
omega = zeros(frameLength,1);
H = ones(N/2+1,Nele);

theta = angle(2);%0*pi/180; %固定一个俯仰角
gamma = [0 90 180 270]*pi/180;%麦克风位置
%gamma = [30 90 150 210 270 330]*pi/180;%麦克风位置
tao = r*cos(theta)*cos(angle(1)-gamma)/c;     %方位角 0 < angle <360
yds = zeros(length(x(:,1)),1);
x1 = zeros(size(x));

% frequency bin weights
for k = 2:1:N/2+1
    omega(k) = 2*pi*(k-1)*fs/N;   
    % steering vector
    H(k,:) = exp(-1j*omega(k)*tao);
end
    
for i = 1:inc:length(x(:,1))-frameLength/2-1

    d = fft(bsxfun(@times, x(i:i+frameLength-1,:),hamming(frameLength)));
%     d = fft(x(i:i+frameLength-1,:).*hamming(frameLength)');
%     x_fft = d(1:129,:).*H;%./abs(d(1:129,:));
    x_fft=bsxfun(@times, d(1:N/2+1,:),H);
    
    % phase transformed
    %x_fft = bsxfun(@rdivide, x_fft,abs(d(1:N/2+1,:)));
    yf = sum(x_fft,2);
    Cf = [yf;conj(flipud(yf(2:N/2)))];
    
    % 恢复延时累加的信号
    yds(i:i+frameLength-1) = yds(i:i+frameLength-1)+(ifft(Cf));
    
    
    % 恢复各路对齐后的信号
    xf  = [x_fft;conj(flipud(x_fft(2:N/2,:)))];
    x1(i:i+frameLength-1,:) = x1(i:i+frameLength-1,:)+(ifft(xf));
end
DS = real(yds)/Nele;  
x1 = real(x1);


end

