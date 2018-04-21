function loss_perc = info_loss(baseline_tc,amp_tc)
%clear all

nu_neuron = 101;
theta_pref = [-(nu_neuron-1)/2:(nu_neuron-1)/2];

sigma_tc = 30;

if nargin < 2
    baseline_tc=10;
    amp_tc = 80;
end

time_max=500;
time = [1:time_max];
sigma_t = 100;
speed_t = exp(-(time-time_max/2).^2/(2*sigma_t^2));
speed_t = speed_t-min(speed_t);
speed_t = speed_t/max(speed_t);

for j=1:time_max
    resp(j,:) = speed_t(j)*amp_tc*exp(-(theta_pref).^2/(2*sigma_tc^2))+ (1- 0.99*(speed_t(j)/max(speed_t)))* baseline_tc;
    resp_der(j,:) = -speed_t(j)*amp_tc*(theta_pref)/sigma_tc^2.*exp(-(theta_pref).^2/(2*sigma_tc^2));
    %Optimal weights
    weights(j,:) = abs(resp_der(j,:))./resp(j,:);
end




for j=1:nu_neuron
    info_opt_n(j) = sum(resp_der(:,j).^2./resp(:,j));
    if sum(weights(:,j).^2.*resp(:,j))==0
        info_weighted_n(j) = 0;
    else
        info_weighted_n(j) = sum(weights(:,j).*resp_der(:,j))^2 / sum(weights(:,j).^2.*resp(:,j));
    end
    if sum(weights(:,j).^2.*resp(:,j))==0
        info_sum_n(j) = 0;
    else
        info_sum_n(j) = sum(resp_der(:,j))^2 / sum(resp(:,j));
    end
end


info_opt = sum(info_opt_n);
%weights = resp_der_t./resp_t;
info_weighted = sum(info_weighted_n);
info_sum = sum(info_sum_n);
loss = info_opt-info_sum;
loss_perc = loss/info_opt*100;

if nargin < 2
    
    fprintf('info_otp %f  ',info_opt)
    fprintf('info_sum %f  ',info_sum)
    fprintf('info_weighted %f  ',info_weighted)
    fprintf('Percentage loss %f \n',loss_perc)
    
    subplot(221)
    surf(resp);
    shading INTERP
    
    subplot(222)
    surfl(resp_der)
    shading INTERP
    
    subplot(223)
    plot(info_sum_n,'b');
    hold on
    plot(info_opt_n,'.r');
    hold off
    
    subplot(224)
    plot(weights(:,10))
end
