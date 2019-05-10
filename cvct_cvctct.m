clc
clear

%parpool(4)
%%
%%Uniform motion:
T=5; sigma_v=0.1;sigma_w=100;
Fu=[1  T  0  0; 0  1  0  0; 0  0  1  T; 0  0  0  1];
Au=[1/2*T^2   0; T   0; 0   1/2*T^2; 0   T];
Hu=[1  0  0  0; 0  0  1  0];
Quv=sigma_v^2*eye(2,2);
Qu=Au*Quv*Au';
Ru=sigma_w^2*eye(2,2);

%%Maneuver:
% omega=pi/180;
% Fm=Fm_CT(omega,1);
Am=[1/2*T^2  0  0; T  0   0; 0  1/2*T^2   0; 0  T   0; 0   0   T];
Hm=[1  0  0  0  0; 0  0  1  0  0];
Qmv=sigma_v^2*eye(3,3);
%Qm0=Am*Qmv*Am';
Fa=[1 T 0 0 0.5*T^2 0;0 1 0 0 T 0;0 0 1 T 0 0.5*T^2;0 0 0 1 0 T;0 0 0 0 1 0;0 0 0 0 0 1];
Ha=[1 0 0 0 0 0;0 0 1 0 0 0];
Aa=[T^3/6 0;T^2/2 0;0 T^3/6;0 T^2/2;T 0;0 T];
Ra=sigma_w^2*eye(2,2);
Qa=Aa*Quv*Aa';

Rm=sigma_w^2*eye(2,2);
rN=20;unormal_vector_IMM_CT=ones(1,rN);unormal_vector_IMM_CT_PF=ones(1,rN);
K=125/T+90/T+125/T+30/T+125/T+1;

RMSE_KF1_n=zeros(rN,K+1);RMSE_IMM1_n=zeros(rN,K+1);RMSE_WIMM1_n=zeros(rN,K+1);RMSE_WKF1_n=zeros(rN,K+1);
rho = 0.1;
xk1_n=cell(rN,1);
xk_IMM_L1_n=cell(rN,1);
xk_WIMM1_n=cell(rN,1);
xk_KF1_n=cell(rN,1);
xk_WKF1_n=cell(rN,1);

parfor    rn=1:rN
    RMSE_KF=zeros(1,K+1);RMSE_IMM=zeros(1,K+1);RMSE_WIMM=zeros(1,K+1);RMSE_WKF=zeros(1,K+1);
    %%The ATC Scenario
    %states of target
    randn('state',rn);% 100 300 400 2000
    xk=zeros(7,K+1);
    
    xk(:,1)=[25000  -120   10000  0   0   0   0]';% the initial value of state
    
    for   k=1:(125/T+1)  %the aircraft flies westward for 125 s at 120 m/s
        xk(1:4,k+1)=Fu*xk(1:4,k)+Au*mvnrnd([0 0],Quv)';
        
    end
    
    omega=pi/180;%executing a 1 degree/s coordinated turn for 90s
    Fm=Fm_CT(omega,T);
    for   k=(125/T+2):(125/T+90/T+1)
        
        xk(1:5,k+1)=Fm*xk(1:5,k)+Am*mvnrnd([0 0 0],Qmv)';
        
    end
    
    for   k=(125/T+90/T+2):(125/T+90/T+125/T+1)%then it flies southward for another 125s
        xk(1:4,k+1)=Fu*xk(1:4,k)+Au*mvnrnd([0 0],Quv)';
        
    end
    
    omega=-pi/60;%executing a 3 degree/s coordinated turn for 30s
    Fm=Fm_CT(omega,T);
    for   k=(125/T+90/T+125/T+2):(125/T+90/T+125/T+30/T+1)
        xk(1:5,k+1)=Fm*xk(1:5,k)+Am*mvnrnd([0 0 0],Qmv)';
        
    end
    
    for   k=(125/T+90/T+125/T+30/T+2):(125/T+90/T+125/T+30/T+125/T+1)%then it flies southward(i think it should be westward) for another 125s
        xk(1:4,k+1)=Fu*xk(1:4,k)+Au*mvnrnd([0 0],Quv)';
        
    end
    
    %          figure(1)
    %          plot(xk(1,1:K),xk(3,1:K))
    %          hold on
    %          axis([-2*10^4 3*10^4  -2.5*10^4  2.5*10^4])
    %          plot(0,0,'x')
    %          text(-1000,1000,'R')
    %          text(xk(1,1)+500,xk(3,1)+1500,'I')
    %          text(xk(1,K)+500,xk(3,K)+1500,'F')
    
    %measurements of target
    zk=zeros(2,K+1);
    for    k=1:(K+1)
        zk(1,k)=xk(1,k)+ sigma_w*randn;
        zk(2,k)=xk(3,k)+ sigma_w*randn;
    end
    
    %%  Kalman filter
         Quv_KF=1^2*eye(2,2);%compare 0.1^2
         Qu_KF=Au*Quv_KF*Au';
         Fu_KF=Fu;
         xk_KF=zeros(4,K+1);

         xk_plus=zeros(4,1);xk_plus(1,1)=zk(1,1); xk_plus(3,1)=zk(2,1);xk_KF(:,1)=xk_plus;%initial value of KF
         Pk_plus=100*eye(4,4);
         for   k=1:K
             [xk_plus,Pk_plus]=Kalman_filter(Fu_KF,Qu_KF,Hu,Ru,xk_plus,Pk_plus, zk(:,k+1));%Kalman filter
             xk_KF(:,k+1)=xk_plus;
             RMSE_KF(1,k+1)=(xk_plus(1)-xk(1,k+1))^2+(xk_plus(3)-xk(3,k+1))^2;%RMSE position estimation errors
             fprintf('run unpurt %d simulation %d times KF\n',rn,k)
         end
    RMSE_KF1_n(rn,:)=RMSE_KF;
    %%  WKF
    Quv_WKF=1^2*eye(2,2);%compare 0.1^2
    Qu_WKF=Au*Quv_WKF*Au';
    Fu_WKF=Fu;
    xk_WKF=zeros(4,K+1);
    
    xk_plus=zeros(4,1);xk_plus(1,1)=zk(1,1); xk_plus(3,1)=zk(2,1);xk_WKF(:,1)=xk_plus;%initial value of KF
    Pk_plus=100*eye(4,4);
    for k=1:K
        [xk_plus,Pk_plus]=My_WKF(rho,Fu_WKF,Qu_WKF,Hu,Ru,xk_plus,Pk_plus,zk(:,k+1));
        xk_WKF(:,k+1)=xk_plus;
        RMSE_WKF(1,k+1)=(xk_plus(1)-xk(1,k+1))^2+(xk_plus(3)-xk(3,k+1))^2;%RMSE position estimation errors
        fprintf('run unpurt %d simulation %d times WKF\n',rn,k)
    end
    RMSE_WKF1_n(rn,:)=RMSE_WKF;
    %% IMM-L_u
    %IMM-L: An IMM estimator with two second-order linear kinematic
    %models(WNA) with two noise levels
    Quv_L1=0.1^2*eye(2,2); Quv_L2=0.1^2*eye(3,3);Quv_L3=0.1^2*eye(3,3); %The one with lower noise level with standard deviation 0.1 m/s^2 and the other one with 2 m/s^2
    Qu_L1=Au*Quv_L1*Au'; Qu_L2=Am(1:4,:)*Quv_L2*Am(1:4,:)';Qu_L3=Am(1:4,:)*Quv_L3*Am(1:4,:)';
    Pi_L=[0.9  0.05 0.05; 0.05  0.9 0.05;0.05 0.05 0.9];%The mode transition probability matrix
    Fu_1=Fu;Fu_22=Fm_CT(pi/90,T);Fu_2=Fu_22(1:4,1:4);Fu_33=Fm_CT(-pi/45,T);Fu_3=Fu_33(1:4,1:4);
    
    xk_IMM_L=zeros(4,K+1);
    
    xk_plus=zeros(4,1);xk_plus(1,1)=zk(1,1); xk_plus(3,1)=zk(2,1);xk_IMM_L(:,1)=xk_plus;%initial value of IMM_L
    Pk_plus=100*eye(4,4);muk_plus=[1/2 1/4 1/4]';
    xk_plus_1=xk_plus;xk_plus_2=xk_plus;xk_plus_3=xk_plus;
    Pk_plus_1=Pk_plus;Pk_plus_2=Pk_plus;Pk_plus_3=Pk_plus;
    for   k=1:K
        [xk_plus,Pk_plus,xk_plus_1,Pk_plus_1,xk_plus_2,Pk_plus_2,xk_plus_3,Pk_plus_3, muk_plus]...
            =IMM_L_3(Fu_1,Qu_L1,Fu_2,Qu_L2,Fu_3,Qu_L3,Hu,Ru, Pi_L,xk_plus_1,Pk_plus_1,xk_plus_2,Pk_plus_2, xk_plus_3,Pk_plus_3,muk_plus, zk(:,k+1));%IMM_L with two models
        xk_IMM_L(:,k+1)=xk_plus;
        RMSE_IMM(1,k+1)=(xk_plus(1)-xk(1,k+1))^2+(xk_plus(3)-xk(3,k+1))^2;%RMSE position estimation errors
        fprintf('run unpurt %d simulation %d times IMM\n',rn,k)
    end
    RMSE_IMM1_n(rn,:)=RMSE_IMM;
    %% WIMM_u
    Quv_L1=0.1^2*eye(2,2); Quv_L2=0.1^2*eye(3,3);Quv_L3=0.1^2*eye(3,3); %The one with lower noise level with standard deviation 0.1 m/s^2 and the other one with 2 m/s^2
    Qu_L1=Au*Quv_L1*Au'; Qu_L2=Am(1:4,:)*Quv_L2*Am(1:4,:)';Qu_L3=Am(1:4,:)*Quv_L3*Am(1:4,:)';
    Pi_L=[0.9  0.05 0.05; 0.05  0.9 0.05;0.05 0.05 0.9];%The mode transition probability matrix
    Fu_1=Fu;Fu_22=Fm_CT(pi/90,T);Fu_2=Fu_22(1:4,1:4);Fu_33=Fm_CT(-pi/45,T);Fu_3=Fu_33(1:4,1:4);
    
    xk_WIMM=zeros(4,K+1);
    
    xk_plus=zeros(4,1);xk_plus(1,1)=zk(1,1); xk_plus(3,1)=zk(2,1);xk_WIMM(:,1)=xk_plus;%initial value of IMM_L
    Pk_plus=100*eye(4,4);muk_plus=[1/2 1/4 1/4]';
    xk_plus_1=xk_plus;xk_plus_2=xk_plus;xk_plus_3=xk_plus;
    Pk_plus_1=Pk_plus;Pk_plus_2=Pk_plus;Pk_plus_3=Pk_plus;
    for   k=1:K
        [xk_plus,Pk_plus,muk_plus]...
            =WIMM_3(rho,Fu_1,Qu_L1,Fu_2,Qu_L2,Fu_3,Qu_L3,Hu,Ru, Pi_L,xk_plus,Pk_plus, muk_plus, zk(:,k+1));%IMM_L with two models
        xk_WIMM(:,k+1)=xk_plus;
        RMSE_WIMM(1,k+1)=(xk_plus(1)-xk(1,k+1))^2+(xk_plus(3)-xk(3,k+1))^2;%RMSE position estimation errors
        fprintf('run unpurt %d simulation %d times WIMM\n',rn,k)
    end
    RMSE_WIMM1_n(rn,:)=RMSE_WIMM;
    
    xk1_n{rn}=xk;
    xk_IMM_L1_n{rn}=xk_IMM_L;
    xk_WIMM1_n{rn}=xk_WIMM;
    xk_KF1_n{rn}=xk_KF;
    xk_WKF1_n{rn}=xk_WKF;
end

