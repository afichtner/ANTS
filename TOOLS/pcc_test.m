t=load('test_pcc_t.txt');
u=load('test_pcc_u.txt');
pcc=load('test_pcc.txt');


figure(1)
subplot(2,1,1)
hold on
plot(t(1:450),u(1:450),'g','LineWidth',1.5);
xlabel('Time (s)')
ylabel('Amplitude')
legend('Test cosine')

subplot(2,1,2)
hold on
plot(linspace(-1000,1000,length(pcc)),pcc,'g','LineWidth',1.5);
xlabel('Lag (s)')
ylabel('Amplitude')
legend('Phase cross correlation')

