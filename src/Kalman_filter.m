function  [xk_plus,Pk_plus]=Kalman_filter(F,Q,H,R,xk_plus,Pk_plus, zk)
                
 % prediction 
 xk_minus=F*xk_plus;  %һ��Ԥ��
 Pk_minus=F*Pk_plus*F'+Q; %һ��Ԥ��������
 
 %update
 S=H*Pk_minus*H'+R;%inovation covariance
 W=Pk_minus*H'*inv(S);%Filter gain
 
 xk_plus=xk_minus+W*(zk-H*xk_minus);%updated state estimate
 Pk_plus=Pk_minus-W*S*W';%updated state covariance

end