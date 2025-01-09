function diff=getstates(t, y)
    global F
    global m
    global I
    global M
    %% Allocating the States Vector
    velocities=y(1:3); % Allocating velocity positions
    angular_velocities=y(4:6); % Allocating angular velocity positions
    angles=y(7:9); % Allocating angles positions
    position=y(10:12); % Allocating xyz positions
    %% Differential Equations Setup
    velocities_dot=(1/m)*F - cross(angular_velocities ,velocities); % Velocity derivative equations
    angular_velocities_dot=I\(M-cross(angular_velocities,I*angular_velocities)); % Angular Velocity derivative equations 
    angles_dot=[1, sin(angles(1)) * tan(angles(2)), cos(angles(1)) * tan(angles(2));
         0, cos(angles(1)), -sin(angles(1));
        0, sin(angles(1)) / cos(angles(2)), cos(angles(1)) / cos(angles(2))] * angular_velocities; % Angles derivative equations
    positions_dot=eul2rotm([angles(3) angles(2) angles(1)])*velocities; % Position derivatives equations
    diff=[velocities_dot; angular_velocities_dot; angles_dot; positions_dot]; % Creating the differential equation matrix
end