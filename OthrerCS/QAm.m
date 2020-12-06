clc

M=16;
k=log2(M);
n=100000;                          %比特序列长度
samp=1;                            %过采样率
x=randi([0 1], n, 1);                    %生成随机二进制比特流
stem(x(1:50),'filled');            %画出相应的二进制比特流信号
title('二进制随机比特流');
xlabel('比特序列');ylabel('信号幅度');
x4=reshape(x,k,length(x)/k);       %将原始的二进制比特序列每四个一组分组，并排列成k行length(x)/k列的矩阵
xsym=bi2de(x4.','left-msb');       %将矩阵转化为相应的16进制信号序列
figure;
stem(xsym(1:50));                  %画出相应的16进制信号序列
title('16进制随机信号');
xlabel('信号序列');ylabel('信号幅度');
y= qammod(xsym, M);  %用16QAM调制器对信号进行调制
scatterplot(y);                    %画出16QAM信号的星座图
text(real(y)+0.1,imag(y),dec2bin(xsym));  %text(x,y,'string')：在二维图形中指定的位置(x,y)上显示字符串string
axis([-5 5 -5 5]);
EbNo=15;
snr=EbNo+10*log10(k)-10*log10(samp); %信噪比
yn=awgn(y,snr,'measured');         %加入高斯白噪声
h=scatterplot(yn,samp,0,'b.');     %经过信道后接收到的含白噪声的信号星座图
hold on;
scatterplot(y,1,0,'k+',h);         %加入不含白噪声的信号星座图
title('接收信号星座图');
legend('含噪声接收信号','不含噪声信号');
axis([-5 5 -5 5]);
hold on;
eyediagram(yn,2);                  %眼图
yd = qamdemod(yn, M); %此时解调出来的是16进制信号
z=de2bi(yd,'left-msb');            %转化为对应的二进制比特流
z=reshape(z.',numel(z),1');
[number_of_errors,bit_error_rate]=biterr(x,z)
% berfit 用平滑的曲线尽量把离散点连接起来
