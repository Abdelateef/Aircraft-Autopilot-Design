function [t, y] = RK4_solver(ode, y0, t0, tf, n,SD,Vto, dc)
    global M 
    global F 
    h = (tf - t0) / n;         % Step size
    t = linspace(t0, tf, n+1); % Time vector
    y = zeros(length(y0), n+1);% Solution matrix
    y(:,1) = y0;               % Set initial condition

    % RK4 Loop
    for i = 1:n
        k1 =  h*(ode(t(i), y(:,i)));
        k2 =  h*(ode(t(i) + 0.5*h, y(:,i) + 0.5*k1));
        k3 =  h*(ode(t(i) + 0.5*h, y(:,i) + 0.5*k2));
        k4 =  h*(ode(t(i) + h, y(:,i) + k3));
        y(:,i+1) = y(:,i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
        if i == 1
            wd = (y(3,i+1)-y0(3))/h;
        else
            wd = (y(3,i+1)-y(3, i))/h;
        end
        [F,M]=F_M_Cal(SD,y0,Vto,dc,y(:,i+1),wd);
     
    end
    y=y(:,1:n); % to remove the last row 
    t=t(:,1:n);
end
