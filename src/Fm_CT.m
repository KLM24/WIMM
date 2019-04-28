function   Fm=Fm_CT(omega,T)

if  omega==0
        Fm=[1   T  0  0   0
                 0   1  0  0   0
                 0   0  1  T   0
                 0   0  0  1   0
                 0   0  0  0  1];
else
Fm=[1  sin(omega*T)/omega          0           -(1-cos(omega*T))/omega   0
         0  cos(omega*T)                     0           -sin(omega*T)                    0 
         0  (1-cos(omega*T))/omega   1           sin(omega*T)/omega          0 
         0  sin(omega*T)                      0           cos(omega*T )                   0
         0        0                                    0                    0                                1]; 
end
end