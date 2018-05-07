function circ_aver = circAverage(A)
    step = round(size(A,1)/2);
    B = nan(size(A));
    for i = 1:size(A,1)
        B(i,:) = circshift(A(i,:),[0 step-i]);
    end
    circ_aver = mean(B,1);
end
