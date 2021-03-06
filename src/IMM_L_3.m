function  [xk_plus,Pk_plus,xk_plus_1,Pk_plus_1,xk_plus_2,Pk_plus_2,xk_plus_3,Pk_plus_3, muk_plus] = IMM_L_3(Fu_1,Qu_L1,Fu_2,Qu_L2,Fu_3,Qu_L3,Hu,Ru, Pi_L,xk_plus_1,Pk_plus_1,xk_plus_2,Pk_plus_2,xk_plus_3,Pk_plus_3, muk_plus, zk)
    
                 r=3;
                %1. Calculation of the mixing probabilities
                muk_mix=zeros(r,r);
                
                c_bar=Pi_L'*muk_plus;
                for   i=1:r
                       for    j=1:r
                               muk_mix(i,j)=Pi_L(i,j)*muk_plus(i)/c_bar(j);
                       end
                end
                %2. Mixing
                xk_01=xk_plus_1*muk_mix(1,1)+xk_plus_2*muk_mix(2,1)+xk_plus_3*muk_mix(3,1);
                
                xk_02=xk_plus_1*muk_mix(1,2)+xk_plus_2*muk_mix(2,2)+xk_plus_3*muk_mix(3,2);
                
                xk_03=xk_plus_1*muk_mix(1,3)+xk_plus_2*muk_mix(2,3)+xk_plus_3*muk_mix(3,3)
                
                Pk_01=muk_mix(1,1)*(Pk_plus_1+(xk_plus_1-xk_01)*(xk_plus_1-xk_01)')...
                           +muk_mix(2,1)*(Pk_plus_2+(xk_plus_2-xk_01)*(xk_plus_2-xk_01)')+muk_mix(3,1)*(Pk_plus_3+(xk_plus_3-xk_01)*(xk_plus_3-xk_01)')
                
                Pk_02=muk_mix(1,2)*(Pk_plus_1+(xk_plus_1-xk_02)*(xk_plus_1-xk_02)')...
                           +muk_mix(2,2)*(Pk_plus_2+(xk_plus_2-xk_02)*(xk_plus_2-xk_02)')+muk_mix(3,2)*(Pk_plus_3+(xk_plus_3-xk_02)*(xk_plus_3-xk_02)');
                
                Pk_03=muk_mix(1,3)*(Pk_plus_1+(xk_plus_1-xk_03)*(xk_plus_1-xk_03)')...
                           +muk_mix(2,3)*(Pk_plus_2+(xk_plus_2-xk_03)*(xk_plus_2-xk_03)')+muk_mix(3,3)*(Pk_plus_3+(xk_plus_3-xk_03)*(xk_plus_3-xk_03)');
                       
                       
                %3. Mode-matched filtering
                xk_minus_1=Fu_1*xk_01;
                xk_minus_2=Fu_2*xk_02;
                xk_minus_3=Fu_3*xk_03;
                
                Pk_minus_1=Fu_1*Pk_01*Fu_1'+Qu_L1
                Pk_minus_2=Fu_2*Pk_02*Fu_2'+Qu_L2;
                Pk_minus_3=Fu_3*Pk_03*Fu_3'+Qu_L3;
                
                S_1=Hu*Pk_minus_1*Hu'+Ru;
                S_2=Hu*Pk_minus_2*Hu'+Ru;
                S_3=Hu*Pk_minus_3*Hu'+Ru;
                
                W_1=Pk_minus_1*Hu'*inv(S_1);
                W_2=Pk_minus_2*Hu'*inv(S_2);
                W_3=Pk_minus_3*Hu'*inv(S_3);
                
                nu_1=zk-Hu*xk_minus_1;
                nu_2=zk-Hu*xk_minus_2;
                nu_3=zk-Hu*xk_minus_3; 
                
                xk_plus_1=xk_minus_1+W_1*nu_1;
                xk_plus_2=xk_minus_2+W_2*nu_2;
                xk_plus_3=xk_minus_3+W_3*nu_3;
                
                Pk_plus_1=Pk_minus_1-W_1*S_1*W_1';
                Pk_plus_2=Pk_minus_2-W_2*S_2*W_2';
                Pk_plus_3=Pk_minus_3-W_3*S_3*W_3';
                

                Lambda_1=mvnpdf(nu_1, [0 0]', (S_1+S_1')/2);
                %Lambda_1=mvnpdf(zk,Hu*xk_minus_1,S_1);
                Lambda_2=mvnpdf(nu_2, [0 0]', (S_2+S_2')/2);
                %Lambda_2=mvnpdf(zk,Hu*xk_minus_2,S_2);
                Lambda_3=mvnpdf(nu_3, [0 0]', (S_3+S_3')/2);
                %4. Mode probability update
                c=Lambda_1*c_bar(1)+Lambda_2*c_bar(2)+Lambda_3*c_bar(3);
                muk_plus(1)=Lambda_1*c_bar(1)/c;
                muk_plus(2)=Lambda_2*c_bar(2)/c;
                muk_plus(3)=Lambda_3*c_bar(3)/c;
                
                %5. Estimate and covariance combination

                
                xk_plus=xk_plus_1*muk_plus(1)+xk_plus_2*muk_plus(2)+xk_plus_3*muk_plus(3);
                Pk_plus=muk_plus(1)*(Pk_plus_1+(xk_plus_1-xk_plus)*(xk_plus_1-xk_plus)')...
                             +muk_plus(2)*(Pk_plus_2+(xk_plus_2-xk_plus)*(xk_plus_2-xk_plus)')+muk_plus(3)*(Pk_plus_3+(xk_plus_3-xk_plus)*(xk_plus_3-xk_plus)');

end