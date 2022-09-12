function [] = Plotposition(X_hat_InEKF,X_true,title_str)
N = length(X_hat_InEKF);
figure
for i=1:3
    subplot(3,1,i)
    temp_a = reshape(X_hat_InEKF(1:3,5,:),[3,N]);
    plot(1:N,temp_a(i,:));
    hold on
    temp_b = reshape(X_true(1:3,5,:),[3,N]);
    plot(1:N,temp_b(i,:));
    hold off
    legend('Estiamtion','GroundTruth')
    xlabel('time step')
    switch i
        case 1
            ylabel('x(m)')
        case 2
            ylabel('y(m)')
        case 3
            ylabel('z(m)')
    end
end
sgtitle(title_str+": Estimation of linear position")
disp(title_str+" --- ATE of position: " +num2str(norm(temp_b-temp_a)))
end

