function cost = CostFunc(model_y,real_ys,weight,alpha,increase_step_trick)
% global increase_step_trick alpha;
weight = weight-increase_step_trick;
proj_PSTH = reshape((weight) * reshape(real_ys,size(real_ys,1),[]),[],3);
cost = nanmean((proj_PSTH(:) - model_y(:)).^2) + alpha * norm(weight); % Mean Squared error
end
