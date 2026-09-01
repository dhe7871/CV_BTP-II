function plotFigures(params, state)
    dt = params.dt;
    simTime = params.simTime;
    
    figure;
    subplot(2, 2, 1);
    plot(0:dt:(simTime-dt), state(7, 1:end-1)*180/pi,"--", linewidth=1.5);
    hold on;
    grid on; grid minor;
    plot(0:dt:(simTime-dt), state(6, 1:end-1)*180/pi, lineWidth=1.5);
    xlabel("Time [s]");
    ylabel("\chi, [deg]");
    legend("Desired","Actual");
    title("Course angle variation with time");
    hold off;
    
    subplot(2, 2, 2);
    plot(0:dt:(simTime-dt), state(10, 1:end-1), lineWidth=1.5);
    hold on;
    grid on; grid minor;
    xlabel("Time [s]");
    ylabel("Cross Track Error");
    title("Cross Track Error variation with time");
    hold off;
    
    subplot(2, 2, 3);
    plot(0:dt:(simTime-dt), state(5, 1:end-1), linewidth=1.5);
    hold on;
    grid on; grid minor;
    xlabel("Time [s]");
    ylabel("V_g, [m/s]");
    title("Ground speed variation with time");
    hold off;
    
    subplot(2, 2, 4);
    plot(0:dt:(simTime-dt), state(8, 1:end-1),"--", linewidth=1.5);
    hold on;
    grid on; grid minor;
    plot(0:dt:(simTime-dt), state(9, 1:end-1), linewidth=1.5);
    xlabel("Time [s]");
    ylabel("d\chi/dt, [rad/s]");
    legend("Desired", "commanded");
    title("Course rate variation with time");
    hold off;
end