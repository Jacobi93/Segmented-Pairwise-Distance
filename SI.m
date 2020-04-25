% derivative dynamic time warping of two signals
function mean_sil = SI(data)
% Silhouette index for the CM dataset
% there are 21 time series in the CM dataset
% there are 7 surgeons with each performing 3 surgeries
% on the same patient (bone)
% this function should be modified if the dataset is different

sil= zeros(21,1);
for i = 1:21
    x = data(i,:);
    if mod(i,3) == 0
        a = (x(i-1)+x(i-2))/2;
        x((i-2):i)=[];
    elseif mod(i,3) == 1
        a = (x(i+1)+x(i+2))/2;
        x(i:(i+2))=[];
    else
        a = (x(i-1)+x(i+1))/2;
        x((i-1):(i+1))=[];
    end   
    % 18 numbers are left
    diff_clus = zeros(6,1);
    for j = 1:6
        diff_clus(j) = (x(3*j-2)+x(3*j-1)+x(3*j))/3;
    end
    c = min(diff_clus);
    sil(i) = (c-a)/max(a,c);
end
mean_sil = sum(sil)/21;

end
