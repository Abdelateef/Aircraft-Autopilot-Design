function [F_vec,M_vec]=F_M_Cal(SD, s0,V,dc,s,wd)
  global m g I
    % Force and moment matrix, based on stability derivatives


    phi=s(7);theta=s(8);psi=s(9); 
    u=s(1);v=s(2);w=s(3);
    p=s(4);q=s(5);r=s(6);


    phi_0=s0(7);theta_0=s0(8);psi_0=s0(9);
    u0=s0(1);v0=s0(2);w0=s0(3); wd0=0;
    p0=s0(4);q0=s0(5);r0=s0(6);


    Delta=[u-u0;v-v0;w-w0;wd-wd0;p-p0;q-q0;r-r0;dc(1);dc(2);dc(3);dc(4)];
    F_M_Matrix =  [SD(1) 0 SD(4) 0 0 0 0 0 0 SD(11) SD(14);
                   0 SD(17) 0 0 0 0 0 0 SD(26) 0 0;
                   SD(2) 0 SD(5) SD(7) 0 SD(8) 0 0 0 SD(12) SD(15);
                   0 SD(19)/V 0 0 SD(21) 0 SD(23) SD(27) SD(29) 0 0;
                   SD(3) 0 SD(6) SD(9) 0 SD(10) 0 0 0 SD(13) SD(16);
                   0 SD(20)/V 0 0 SD(22) 0 SD(24) SD(28) SD(30) 0 0];
    
    Delta_F_M_0=F_M_Matrix*Delta;
    F_vec=Delta_F_M_0(1:3)*m +[m * g * sin(theta_0); -m * g * cos(theta_0) * sin(phi_0); -m * g * cos(theta_0) * cos(phi_0)]+[-m * g * sin(theta); m * g * cos(theta) * sin(phi); m * g * cos(theta) * cos(phi)];
    M_vec=[Delta_F_M_0(4)*I(1,1);Delta_F_M_0(5)*I(2,2);Delta_F_M_0(6)*I(3,3)];

end