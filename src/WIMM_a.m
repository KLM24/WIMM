function [xk_plus,Pk_plus,xk_plus_1,Pk_plus_1,xk_plus_2,Pk_plus_2, muk_plus] = WIMM_a(rho,Fu_1,Qu_L1,Fu_2,Qu_L2,Hu,Ru, Pi_L,xk_plus,Pk_plus,muk_plus,zk)
r=2;
muk_mix=zeros(r,r);

c_bar=Pi_L'*muk_plus;
for   i=1:r
    for    j=1:r
        muk_mix(i,j)=Pi_L(i,j)*muk_plus(i)/c_bar(j);
    end
end

A_aug_1 = [Fu_1; Hu * Fu_1];
B_aug_1 = [eye(6) zeros(6,2);Hu eye(2)];
mu_t_1 = A_aug_1 * xk_plus;
sigma_t_1 = A_aug_1 * Pk_plus * A_aug_1' + B_aug_1 * [Qu_L1 zeros(6,2);zeros(2,6) Ru] * B_aug_1';
A_aug_2 = [Fu_2; Hu * Fu_2];
B_aug_2 = [eye(6) zeros(6,2);Hu eye(2)];
mu_t_2 = A_aug_2 * xk_plus;
sigma_t_2 = A_aug_2 * Pk_plus * A_aug_2' + B_aug_2 * [Qu_L2 zeros(6,2);zeros(2,6) Ru] * B_aug_2';


[phi_star_1, Q_star_1] = F_W(mu_t_1, sigma_t_1, rho, 6);
G_t_1 = phi_star_1.G;
S_t_1 = Q_star_1.Sigma;
Pk_plus_1 = S_t_1(1:6, 1:6) - G_t_1 * S_t_1(7:end, 1:6);
xk_plus_1 = G_t_1*(zk - mu_t_1(7:end,1)) + mu_t_1(1:6,1);

[phi_star_2, Q_star_2] = F_W(mu_t_2, sigma_t_2, rho, 6);
G_t_2 = phi_star_2.G;
S_t_2 = Q_star_2.Sigma;
Pk_plus_2 = S_t_2(1:6, 1:6) - G_t_2 * S_t_2(7:end, 1:6);
xk_plus_2 = G_t_2*(zk - mu_t_2(7:end,1)) + mu_t_2(1:6,1);

Lambda_1 = mvnpdf(zk-mu_t_1(7:end,1),[0 0]',S_t_1(7:end,7:end));
Lambda_2 = mvnpdf(zk-mu_t_2(7:end,1),[0 0]',S_t_2(7:end,7:end));
c=Lambda_1*c_bar(1)+Lambda_2*c_bar(2);
muk_plus(1)=Lambda_1*c_bar(1)/c;
muk_plus(2)=Lambda_2*c_bar(2)/c;
%5. Estimate and covariance combination


xk_plus=xk_plus_1*muk_plus(1)+xk_plus_2*muk_plus(2);
Pk_plus=muk_plus(1)*(Pk_plus_1+(xk_plus_1-xk_plus)*(xk_plus_1-xk_plus)')...
    +muk_plus(2)*(Pk_plus_2+(xk_plus_2-xk_plus)*(xk_plus_2-xk_plus)');
end