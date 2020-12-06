clc

M=16;
k=log2(M);
n=100000;                          %�������г���
samp=1;                            %��������
x=randi([0 1], n, 1);                    %������������Ʊ�����
stem(x(1:50),'filled');            %������Ӧ�Ķ����Ʊ������ź�
title('���������������');
xlabel('��������');ylabel('�źŷ���');
x4=reshape(x,k,length(x)/k);       %��ԭʼ�Ķ����Ʊ�������ÿ�ĸ�һ����飬�����г�k��length(x)/k�еľ���
xsym=bi2de(x4.','left-msb');       %������ת��Ϊ��Ӧ��16�����ź�����
figure;
stem(xsym(1:50));                  %������Ӧ��16�����ź�����
title('16��������ź�');
xlabel('�ź�����');ylabel('�źŷ���');
y= qammod(xsym, M);  %��16QAM���������źŽ��е���
scatterplot(y);                    %����16QAM�źŵ�����ͼ
text(real(y)+0.1,imag(y),dec2bin(xsym));  %text(x,y,'string')���ڶ�άͼ����ָ����λ��(x,y)����ʾ�ַ���string
axis([-5 5 -5 5]);
EbNo=15;
snr=EbNo+10*log10(k)-10*log10(samp); %�����
yn=awgn(y,snr,'measured');         %�����˹������
h=scatterplot(yn,samp,0,'b.');     %�����ŵ�����յ��ĺ����������ź�����ͼ
hold on;
scatterplot(y,1,0,'k+',h);         %���벻�����������ź�����ͼ
title('�����ź�����ͼ');
legend('�����������ź�','���������ź�');
axis([-5 5 -5 5]);
hold on;
eyediagram(yn,2);                  %��ͼ
yd = qamdemod(yn, M); %��ʱ�����������16�����ź�
z=de2bi(yd,'left-msb');            %ת��Ϊ��Ӧ�Ķ����Ʊ�����
z=reshape(z.',numel(z),1');
[number_of_errors,bit_error_rate]=biterr(x,z)
% berfit ��ƽ�������߾�������ɢ����������
