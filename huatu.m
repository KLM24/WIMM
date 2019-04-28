load('cvct_u.mat','-mat')
load('cvct_o.mat','-mat')
load('cvct_n.mat','-mat')
RMSE_IMM_u = sqrt(mean(RMSE_IMM1_u,1));
RMSE_WIMM_u = sqrt(mean(RMSE_WIMM1_u,1));
RMSE_KF_u = sqrt(mean(RMSE_KF1_u,1));
RMSE_WKF_u = sqrt(mean(RMSE_WKF1_u,1));

RMSE_IMM_o = sqrt(mean(RMSE_IMM1,1));
RMSE_WIMM_o = sqrt(mean(RMSE_WIMM1,1));
RMSE_KF_o = sqrt(mean(RMSE_KF1,1));
RMSE_WKF_o = sqrt(mean(RMSE_WKF1,1));

RMSE_IMM_n = sqrt(mean(RMSE_IMM1_n,1));
RMSE_WIMM_n = sqrt(mean(RMSE_WIMM1_n,1));
RMSE_KF_n = sqrt(mean(RMSE_KF1_n,1));
RMSE_WKF_n = sqrt(mean(RMSE_WKF1_n,1));

varIMM = (RMSE_IMM_n-RMSE_IMM_u)./RMSE_IMM_u;
varWIMM = (RMSE_WIMM_n-RMSE_WIMM_u)./RMSE_WIMM_u;
varKF = (RMSE_KF_n-RMSE_KF_u)./RMSE_KF_u;
varWKF = (RMSE_WKF_n-RMSE_WKF_u)./RMSE_WKF_u;

figure(1)
plot(1:K+1,RMSE_IMM_u,'r');
hold on
plot(1:K+1,RMSE_WIMM_u,'g');
hold on
plot(1:K+1,RMSE_KF_u,'b');
hold on
plot(1:K+1,RMSE_WKF_u,'y');
xlabel('State number')
ylabel('RMSE')
title('RMSE for different algorithms under good situation')
legend('IMM','WIMM','KF','WKF')
hold off

figure(2)
plot(xk1_u{rN}(1,:),xk1_u{rN}(3,:),'k');
hold on
plot(xk_IMM_L1_u{rN}(1,:),xk_IMM_L1_u{rN}(3,:),'r');
hold on
plot(xk_WIMM1_u{rN}(1,:),xk_WIMM1_u{rN}(3,:),'g');
hold on
plot(xk_KF1_u{rN}(1,:),xk_KF1_u{rN}(3,:),'b');
hold on
plot(xk_WKF1_u{rN}(1,:),xk_WKF1_u{rN}(3,:),'y');
xlabel('x-axis')
ylabel('y-axis')
title('The last time paths under good situation')
legend('True', 'IMM','WIMM','KF','WKF')
hold off

figure(3)
plot(1:K+1,RMSE_IMM_n,'r');
hold on
plot(1:K+1,RMSE_WIMM_n,'g');
hold on
plot(1:K+1,RMSE_KF_n,'b');
hold on
plot(1:K+1,RMSE_WKF_n,'y');
xlabel('State number')
ylabel('RMSE')
title('RMSE for different algorithms under bad situation')
legend('IMM','WIMM','KF','WKF')
hold off

figure(4)
plot(xk1_n{rN}(1,:),xk1_n{rN}(3,:),'k');
hold on
plot(xk_IMM_L1_n{rN}(1,:),xk_IMM_L1_n{rN}(3,:),'r');
hold on
plot(xk_WIMM1_n{rN}(1,:),xk_WIMM1_n{rN}(3,:),'g');
hold on
plot(xk_KF1_n{rN}(1,:),xk_KF1_n{rN}(3,:),'b');
hold on
plot(xk_WKF1_n{rN}(1,:),xk_WKF1_n{rN}(3,:),'y');
xlabel('x-axis')
ylabel('y-axis')
title('The last time paths under bad situation')
legend('True', 'IMM','WIMM','KF','WKF')
hold off

figure(5)
plot(1:K+1,RMSE_IMM_u,'r')
hold on
plot(1:K+1,RMSE_IMM_n,'g')
xlabel('State number')
ylabel('RMSE')
title('IMM RMSE good VS bad')
legend('good','bad')
hold off

figure(6)
plot(1:K+1,RMSE_WIMM_u,'r')
hold on
plot(1:K+1,RMSE_WIMM_n,'g')
xlabel('State number')
ylabel('RMSE')
title('WIMM RMSE good VS bad')
legend('good','bad')
hold off

figure(7)
plot(1:K+1,RMSE_KF_u,'r')
hold on
plot(1:K+1,RMSE_KF_n,'g')
xlabel('State number')
ylabel('RMSE')
title('KF RMSE good VS bad')
legend('good','bad')
hold off

figure(8)
plot(1:K+1,RMSE_WKF_u,'r')
hold on
plot(1:K+1,RMSE_WKF_n,'g')
xlabel('State number')
ylabel('RMSE')
title('WKF RMSE good VS bad')
legend('good','bad')
hold off

figure(9)
plot(1:K+1,varIMM,'r')
hold on
plot(1:K+1,varWIMM,'g')
hold on
plot(1:K+1,varKF,'b')
hold on
plot(1:K+1,varWKF,'y')
xlabel('State number')
ylabel('RMSE varying rate')
title('RMSE varying rate for all 4 different algorithms')
legend('IMM','WIMM','KF','WKF')
hold off
figure(10)
plot(1:K+1,varIMM,'r')
hold on
plot(1:K+1,varWIMM,'g')
xlabel('State number')
ylabel('RMSE varying rate')
title('RMSE varying rate for IMM and WIMM')
legend('IMM','WIMM')
hold off
figure(11)
plot(1:K+1,varKF,'r')
hold on
plot(1:K+1,varWKF,'g')
xlabel('State number')
ylabel('RMSE varying rate')
title('RMSE varying rate for KF and WKF')
legend('KF','WKF')
hold off

figure(12)
plot(1:K+1,RMSE_IMM_o,'r');
hold on
plot(1:K+1,RMSE_WIMM_o,'g');
hold on
plot(1:K+1,RMSE_KF_o,'b');
hold on
plot(1:K+1,RMSE_WKF_o,'y');
xlabel('State number')
ylabel('RMSE')
title('RMSE for different algorithms under medium situation')
legend('IMM','WIMM','KF','WKF')
hold off

figure(13)
plot(xk1{rN}(1,:),xk1{rN}(3,:),'k');
hold on
plot(xk_IMM_L1{rN}(1,:),xk_IMM_L1{rN}(3,:),'r');
hold on
plot(xk_WIMM1{rN}(1,:),xk_WIMM1{rN}(3,:),'g');
hold on
plot(xk_KF1{rN}(1,:),xk_KF1{rN}(3,:),'b');
hold on
plot(xk_WKF1{rN}(1,:),xk_WKF1{rN}(3,:),'y');
xlabel('x-axis')
ylabel('y-axis')
title('The last time paths under medium situation')
legend('True', 'IMM','WIMM','KF','WKF')
hold off