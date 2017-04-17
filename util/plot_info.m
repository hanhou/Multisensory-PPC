

clear

time_scale = [];
load info64-2
MLinfo(1,:)=ML_info;
info_outx(2:5)=info_outx;
info_outx(1)=0;
info_train_outx(2:5)=info_train_outx;
info_train_outx(1)=0;
info(1,:)=(info_outx+info_train_outx)/2;
info_test(1,:)= info_outx;
info_train(1,:)= info_train_outx;

[aux nu_time_step] = size(ML_info);


time_scale = [0:50:50*(nu_time_step-1)];


load info128
MLinfo(2,:)=ML_info;
info_outx(2:5)=info_outx;
info_outx(1)=0;
info_train_outx(2:5)=info_train_outx;
info_train_outx(1)=0;
info(2,:)=(info_outx+info_train_outx)/2;
info_test(2,:)= info_outx;
info_train(2,:)= info_train_outx;

load info256-2
MLinfo(3,:)=ML_info;
info_outx(2:5)=info_outx;
info_outx(1)=0;
info_train_outx(2:5)=info_train_outx;
info_train_outx(1)=0;
info(3,:)=(info_outx+info_train_outx)/2;
info_test(3,:)= info_outx;
info_train(3,:)= info_train_outx;

load info512
MLinfo(4,:)=ML_info;
info_outx(2:5)=info_outx;
info_outx(1)=0;
info_train_outx(2:5)=info_train_outx;
info_train_outx(1)=0;
info(4,:)=(info_outx+info_train_outx)/2;
info_test(4,:)= info_outx;
info_train(4,:)= info_train_outx;

load info256512
MLinfo(5,:)=ML_info;
info_outx(2:5)=info_outx;
info_outx(1)=0;
info_train_outx(2:5)=info_train_outx;
info_train_outx(1)=0;
info(5,:)=(info_outx+info_train_outx)/2;
info_test(5,:)= info_outx;
info_train(5,:)= info_train_outx;


% load info0512
% MLinfo(6,:)=ML_info;
% info_outx(2:5)=info_outx;
% info_outx(1)=0;
% info_train_outx(2:5)=info_train_outx;
% info_train_outx(1)=0;
% info(6,:)=(info_outx+info_train_outx)/2;
% info_test(6,:)= info_outx;

load info_in64
info_inx(2:5)=info_inx;
info_inx(1)=0;
info_train_inx(2:5)=info_train_inx;
info_train_inx(1)=0;
info_in(2,:)=(info_inx+info_train_inx)/2;
info_in_test(2,:)= info_inx;

load info_in128
info_inx(2:5)=info_inx;
info_inx(1)=0;
info_train_inx(2:5)=info_train_inx;
info_train_inx(1)=0;
info_in(3,:)=(info_inx+info_train_inx)/2;
info_in_test(3,:)= info_inx;

load info_in256
info_inx(2:5)=info_inx;
info_inx(1)=0;
info_train_inx(2:5)=info_train_inx;
info_train_inx(1)=0;
info_in(4,:)=(info_inx+info_train_inx)/2;
info_in_test(4,:)= info_inx;

load info_in512
info_inx(2:5)=info_inx;
info_inx(1)=0;
info_train_inx(2:5)=info_train_inx;
info_train_inx(1)=0;
info_in(5,:)=(info_inx+info_train_inx)/2;
info_in_test(5,:)= info_inx;


subplot(221)
plot(time_scale,info(5,:),'s:')
hold on
plot(time_scale,info(1:4,:)','s-')
hold off
axis([0 (nu_time_step-1)*50 0 1.1*max(max(info))]);
plot(time_scale,[info(5,:)' info(5,:)' info(1,:)' info(2,:)' info(3,:)' info(4,:)'],'s-')


subplot(222)
plot(time_scale,[MLinfo(5,:)' MLinfo(5,:)' MLinfo(1,:)' MLinfo(2,:)' MLinfo(3,:)' MLinfo(4,:)'],'s-')

hold on

% for j=2:5
%    errorbar(time_scale,info(j,:)',info_test(j,:)',info_train(j,:)');
% end
plot(time_scale,[info(5,:)' info(5,:)' info(1,:)' info(2,:)' info(3,:)' info(4,:)'],'o-')
plot(time_scale,[info_test(5,:)' info_test(5,:)' info_test(1,:)'...
    info_test(2,:)' info_test(3,:)' info_test(4,:)'],':')
plot(time_scale,[info_train(5,:)' info_train(5,:)' info_train(1,:)'...
    info_train(2,:)' info_train(3,:)' info_train(4,:)'],':')


hold off
axis([0 (nu_time_step-1)*50 0 1.3*max(max(MLinfo))]);



subplot(223)
plot(time_scale,info_test','s-')

subplot(224)
plot(time_scale,info_in','s-')