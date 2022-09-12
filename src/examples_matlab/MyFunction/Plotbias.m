function [] = Plotbias(bias_hat,bg_true,ba_true)
figure
for i=1:2
    for j=1:3
        subplot(3,2,i+2*(j-1))     
        if i==1
            plot(bias_hat(j,:))
            hold on
            plot(bg_true(j,:))
            hold off 
            legend('estimation','groudtruth')
            title('bg')
        else
            plot(bias_hat(3+j,:))
            hold on
            plot(ba_true(j,:))
            hold off 
            legend('estimation','groudtruth')
            title('ba')
        end
    end
end
sgtitle("Estimation of Bias")
end

