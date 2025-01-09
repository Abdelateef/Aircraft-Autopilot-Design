function SD_Lat=LateralFunc(SD_Lat_dash,Ixz,Izz,Ixx,V)
    
    G=1/(1-Ixz^2/(Ixx*Izz));
    A=Ixz/Ixx;
    B=Ixz/Izz;
    dashed=[SD_Lat_dash(3);SD_Lat_dash(5);SD_Lat_dash(7);SD_Lat_dash(13);SD_Lat_dash(11)
            SD_Lat_dash(4);SD_Lat_dash(6);SD_Lat_dash(8);SD_Lat_dash(14);SD_Lat_dash(12)];
    Matrix=[G 0 0 0 0 G*A 0 0 0 0;
            0 G 0 0 0 0 G*A 0 0 0;
            0 0 G 0 0 0 0 G*A 0 0;
            0 0 0 G 0 0 0 0 G*A 0;
            0 0 0 0 G 0 0 0 0 G*A;
            G*B 0 0 0 0 G 0 0 0 0;
            0 G*B 0 0 0 0 G 0 0 0;
            0 0 G*B 0 0 0 0 G 0 0;
            0 0 0 G*B 0 0 0 0 G 0;
            0 0 0 0 G*B 0 0 0 0 G];


    star=[SD_Lat_dash(9);SD_Lat_dash(10)];
    non_dashed=inv(Matrix)*dashed;
    non_stared=star*V;
    SD_Lat=[SD_Lat_dash(1);SD_Lat_dash(2);non_dashed(1);non_dashed(6);
            non_dashed(2);non_dashed(7);non_dashed(3);non_dashed(8);
            non_stared(1);non_stared(2);non_dashed(5);non_dashed(10);
            non_dashed(4);non_dashed(9)];
end