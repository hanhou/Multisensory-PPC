
true_info= [0.726];
info = [0.059 0.320 0.373 0.418];
info_shuffled=[0.059 0.592 0.837 1.136];
output_rate = [10 100 150 200];
mean_true_info = mean(true_info);
nu_point = length(output_rate);

true_info_sharp = [0.726];
info_sharp = [0.330 0.320  0.267   0.239];
info_shuffled_sharp = [0.601  0.592   0.444    0.365];
output_rate_sharp = [20.016 21.260   23.974 25.121];
nu_point2 = length(output_rate_sharp);

subplot(221);
plot(output_rate,info,'-*r');
hold on 
plot(output_rate,info_shuffled,'-ob');
plot([output_rate(1) output_rate(nu_point)],[mean_true_info mean_true_info]);
%plot(output_rate_sharp,info_sharp,'-oy');
hold off

subplot(222);
plot(output_rate_sharp,info_sharp,'-*r');
hold on 
plot(output_rate_sharp,info_shuffled_sharp,'-ob');
plot([output_rate_sharp(1) output_rate_sharp(nu_point2)],[mean_true_info mean_true_info]);
hold off

subplot(111);
plot(output_rate,info,'-*r');
hold on 
plot(output_rate,info_shuffled,'-.*r');
plot([19 26],[mean_true_info mean_true_info]);
plot(output_rate_sharp,info_sharp,'-og');
plot(output_rate_sharp,info_shuffled_sharp,'-.og');
hold off
%axis([8 35 0 4]);
axis([19 26 0 0.8]);
