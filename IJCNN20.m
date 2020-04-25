%% z-normalization
for m = 1:21
    matFileName = sprintf('ps%d', m);
    data = eval(matFileName);
    data_mean = mean(data);
    std_data = std(data);
    z = sprintf('psz%d', m);
    eval([z '= (data-data_mean)./std_data;']);
end

%% dataset
dataset = {pz1, pz2, pz3, pz4, pz5, pz6, pz7, pz8, pz9, pz10, pz11, pz12, ...
    pz13, pz14, pz15, pz16, pz17, pz18, pz19, pz20, pz21};
dataset2 = {psz1, psz2, psz3, psz4, psz5, psz6, psz7, psz8, psz9, psz10, ...
    psz11, psz12, psz13, psz14, psz15, psz16, psz17, psz18, psz19, psz20, psz21};

%% DTW
tic
result1=zeros(21,21);
for m = 1:21
    for n = 1:21
        if m ~= n
            matFileName1 = sprintf('psz%d', m);
            matFileName2 = sprintf('psz%d', n);
            matData = eval(matFileName1);
            matData2 = eval(matFileName2);
            result1(m,n) = dtw(matData', matData2')/max(length(matData),length(matData2));
        end
    end
end
toc

b = result1;
% Silhouette
sil_dtw = zeros(21,1);
for i = 1:21
    x = b(i,:);
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
    sil_dtw(i) = (c-a)/max(a,c);
end
mean_sil_dtw = sum(sil_dtw)/21;
mean_sil_dtw
% 0.0540

%% SDTW
tic
quan = 0.99;
result_ol99min = zeros(21,21);
result_ol99max = zeros(21,21);
for m = 1:21
    for n = 1:21
        if m < n
            m1 = dataset2{1, m};
            m2 = dataset2{1, n};
            % segment when distance > quantile
            diff1=sqrt(sum((diff(m1)).^2, 2));
            diff2=sqrt(sum((diff(m2)).^2, 2));
            % set the quantile
            q_1 = quantile(diff1, quan);
            q_2 = quantile(diff2, quan);
            loc1 = 1;
            loc2 = 1;
            for i = 1:length(diff1)
                if diff1(i) > q_1
                    loc1 = [loc1, i];
                end
            end
            loc1 = [loc1, length(m1)];
            for i = 1:length(diff2)
                if diff2(i) > q_2
                    loc2 = [loc2, i];
                end
            end
            loc2 = [loc2, length(m2)];
            % complete segmentation
            m1_seg = cell(1,(length(loc1)-1));
            m2_seg = cell(1,(length(loc2)-1));
            for i = 1:(length(loc1)-1)
                m1_seg{1,i} = m1((loc1(i)+1):loc1(i+1),:);
            end
            for i = 1:(length(loc2)-1)
                m2_seg{1,i} = m2((loc2(i)+1):loc2(i+1),:);
            end
            
            % delete very short segments (noise)
            del1 = [];
            for i = 1:(length(loc1)-1)
                a = size(m1_seg{1,i});
                if a(1) < 3
                    del1 = [del1, i];
                end
            end
            m1_seg(:,del1) = [];
            
            del2 = [];
            for i = 1:(length(loc2)-1)
                a = size(m2_seg{1,i});
                if a(1) < 3
                    del2 = [del2, i];
                end
            end
            m2_seg(:,del2) = [];
            
            % calculate dtwd matrix for two segmented trajectories
            inter_result = zeros(length(m1_seg),length(m2_seg));
            for i = 1:(length(m1_seg))
                for j = 1:(length(m2_seg))
                        inter_result(i,j) = dtw(m1_seg{1,i}',m2_seg{1,j}');
                end
            end

            % delete overlap
            % policy 1
            % calculate dtwd matrix for two segmented trajectories
            inter_result1 = inter_result;

            % find minimal numbers for p1
            loc3 = [];
            sum_seg1 = 0;
            for i = 1:(length(m1_seg))
                [M, I] = min(inter_result1(i,:));
                sum_seg1 = sum_seg1 + M;
                loc3 = [loc3, I];
            end
            loc4 = unique(loc3);
            inter_result1(:,loc4) = [];
            % find minimal numbers for rest p2
            new = size(inter_result1);
            for i = 1:new(2)
                sum_seg1 = sum_seg1 + min(inter_result1(:,i));
            end

            % policy 2
            % calculate dtwd matrix for two segmented trajectories
            inter_result2 = inter_result;

            % find minimal numbers for p2
            loc3 = [];
            sum_seg2 = 0;
            for i = 1:(length(m2_seg))
                [M, I] = min(inter_result2(:,i));
                sum_seg2 = sum_seg2 + M;
                loc3 = [loc3, I];
            end
            loc4 = unique(loc3);
            inter_result2(loc4,:) = [];
            % find minimal numbers for rest p2
            new = size(inter_result2);
            for i = 1:new(1)
                sum_seg2 = sum_seg2 + min(inter_result2(i,:));
            end

            % dtw
            result_ol99min(m,n) = min(sum_seg1/(length(m1)+length(m2)), ...
                sum_seg2/(length(m1)+length(m2)));
            result_ol99max(m,n) = max(sum_seg1/(length(m1)+length(m2)), ...
                sum_seg2/(length(m1)+length(m2)));
        end
    end
end
toc

for m = 1:21
    for n = 1:21
        if m > n
            result_ol99min(m,n) = result_ol99min(n,m);
            result_ol99max(m,n) = result_ol99max(n,m);
        end
    end
end

b = result_ol99min;
% Silhouette
sil_ol99min= zeros(21,1);
for i = 1:21
    x = b(i,:);
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
    sil_ol99min(i) = (c-a)/max(a,c);
end
mean_sil_ol99min = sum(sil_ol99min)/21;
mean_sil_ol99min
% 0.1394

%% CIDTW
tic
result1=zeros(21,21);
for m = 1:21
    for n = 1:21
        if m ~= n
            matFileName1 = sprintf('psz%d', m);
            matFileName2 = sprintf('psz%d', n);
            matData = eval(matFileName1);
            matData2 = eval(matFileName2);
            % complexity estimate
            diff1=sqrt(sum((diff(matData)).^2));
            diff2=sqrt(sum((diff(matData2)).^2));
            result1(m,n) = dtw(matData', matData2')/max(length(matData),length(matData2))*...
                (max(diff1,diff2)/min(diff1,diff2));
        end
    end
end
toc

b = result1;
% Silhouette
sil_dtw = zeros(21,1);
for i = 1:21
    x = b(i,:);
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
    sil_dtw(i) = (c-a)/max(a,c);
end
mean_sil_dtw = sum(sil_dtw)/21;
mean_sil_dtw
% 0.1562

%% SCIDTW
tic
quan = 0.99;
result_ol99min = zeros(21,21);
result_ol99max = zeros(21,21);
for m = 1:21
    for n = 1:21
        if m < n
            m1 = dataset2{1, m};
            m2 = dataset2{1, n};
            % segment when distance > quantile
            diff1=sqrt(sum((diff(m1)).^2, 2));
            diff2=sqrt(sum((diff(m2)).^2, 2));
            % set the quantile
            q_1 = quantile(diff1, quan);
            q_2 = quantile(diff2, quan);
            loc1 = 1;
            loc2 = 1;
            for i = 1:length(diff1)
                if diff1(i) > q_1
                    loc1 = [loc1, i];
                end
            end
            loc1 = [loc1, length(m1)];
            for i = 1:length(diff2)
                if diff2(i) > q_2
                    loc2 = [loc2, i];
                end
            end
            loc2 = [loc2, length(m2)];
            % complete segmentation
            m1_seg = cell(1,(length(loc1)-1));
            m2_seg = cell(1,(length(loc2)-1));
            for i = 1:(length(loc1)-1)
                m1_seg{1,i} = m1((loc1(i)+1):loc1(i+1),:);
            end
            for i = 1:(length(loc2)-1)
                m2_seg{1,i} = m2((loc2(i)+1):loc2(i+1),:);
            end
            
            % delete very short segments (noise)
            del1 = [];
            for i = 1:(length(loc1)-1)
                a = size(m1_seg{1,i});
                if a(1) < 3
                    del1 = [del1, i];
                end
            end
            m1_seg(:,del1) = [];
            
            del2 = [];
            for i = 1:(length(loc2)-1)
                a = size(m2_seg{1,i});
                if a(1) < 3
                    del2 = [del2, i];
                end
            end
            m2_seg(:,del2) = [];
            
            % calculate dtwd matrix for two segmented trajectories
            inter_result = zeros(length(m1_seg),length(m2_seg));
            for i = 1:(length(m1_seg))
                for j = 1:(length(m2_seg))
                    % complexity estimate
                    diff1=sqrt(sum((diff(m1_seg{1,i})).^2));
                    diff2=sqrt(sum((diff(m2_seg{1,j})).^2));
                    inter_result(i,j) = dtw(m1_seg{1,i}',m2_seg{1,j}')*...
                        (max(diff1,diff2)/min(diff1,diff2));
                end
            end

            % delete overlap
            % policy 1
            % calculate dtwd matrix for two segmented trajectories
            inter_result1 = inter_result;

            % find minimal numbers for p1
            loc3 = [];
            sum_seg1 = 0;
            for i = 1:(length(m1_seg))
                [M, I] = min(inter_result1(i,:));
                sum_seg1 = sum_seg1 + M;
                loc3 = [loc3, I];
            end
            loc4 = unique(loc3);
            inter_result1(:,loc4) = [];
            % find minimal numbers for rest p2
            new = size(inter_result1);
            for i = 1:new(2)
                sum_seg1 = sum_seg1 + min(inter_result1(:,i));
            end

            % policy 2
            % calculate dtwd matrix for two segmented trajectories
            inter_result2 = inter_result;

            % find minimal numbers for p2
            loc3 = [];
            sum_seg2 = 0;
            for i = 1:(length(m2_seg))
                [M, I] = min(inter_result2(:,i));
                sum_seg2 = sum_seg2 + M;
                loc3 = [loc3, I];
            end
            loc4 = unique(loc3);
            inter_result2(loc4,:) = [];
            % find minimal numbers for rest p2
            new = size(inter_result2);
            for i = 1:new(1)
                sum_seg2 = sum_seg2 + min(inter_result2(i,:));
            end

            % dtw
            result_ol99min(m,n) = min(sum_seg1/(length(m1)+length(m2)), ...
                sum_seg2/(length(m1)+length(m2)));
            result_ol99max(m,n) = max(sum_seg1/(length(m1)+length(m2)), ...
                sum_seg2/(length(m1)+length(m2)));
        end
    end
end
toc

for m = 1:21
    for n = 1:21
        if m > n
            result_ol99min(m,n) = result_ol99min(n,m);
            result_ol99max(m,n) = result_ol99max(n,m);
        end
    end
end

b = result_ol99min;
% Silhouette
sil_ol99min= zeros(21,1);
for i = 1:21
    x = b(i,:);
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
    sil_ol99min(i) = (c-a)/max(a,c);
end
mean_sil_ol99min = sum(sil_ol99min)/21;
mean_sil_ol99min
% 0.2260

%% DDTW
tic
result1=zeros(21,21);
for m = 1:21
    for n = 1:21
        if m ~= n
            matFileName1 = sprintf('psz%d', m);
            matFileName2 = sprintf('psz%d', n);
            matData = eval(matFileName1);
            matData2 = eval(matFileName2);
            ddtw1 = zeros(length(matData),3);
            ddtw2 = zeros(length(matData2),3);
            for i = 2:length(matData)-1
                ddtw1(i,:) = (matData(i,:)-matData(i-1,:) + ...
                    (matData(i+1,:)-matData(i-1,:))./2)./2;               
            end
            ddtw1(1,:)=ddtw1(2,:);
            ddtw1(length(matData),:)=ddtw1(length(matData)-1,:);
            for i = 2:length(matData2)-1
                ddtw2(i,:) = (matData2(i,:)-matData2(i-1,:) + ...
                    (matData2(i+1,:)-matData2(i-1,:))./2)./2;               
            end
            ddtw2(1,:)=ddtw2(2,:);
            ddtw2(length(matData2),:)=ddtw2(length(matData2)-1,:);
            result1(m,n) = dtw(ddtw1', ddtw2')/max(length(ddtw1),length(ddtw2));
        end
    end
end
toc

b = result1;
% Silhouette
sil_dtw = zeros(21,1);
for i = 1:21
    x = b(i,:);
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
    sil_dtw(i) = (c-a)/max(a,c);
end
mean_sil_dtw = sum(sil_dtw)/21;
mean_sil_dtw
% -0.2341

%% SDDTW
tic
quan = 0.99;
result_ol99min = zeros(21,21);
result_ol99max = zeros(21,21);
for m = 1:21
    for n = 1:21
        if m < n
            m1 = dataset2{1, m};
            m2 = dataset2{1, n};
            % segment when distance > quantile
            diff1=sqrt(sum((diff(m1)).^2, 2));
            diff2=sqrt(sum((diff(m2)).^2, 2));
            % set the quantile
            q_1 = quantile(diff1, quan);
            q_2 = quantile(diff2, quan);
            loc1 = 1;
            loc2 = 1;
            for i = 1:length(diff1)
                if diff1(i) > q_1
                    loc1 = [loc1, i];
                end
            end
            loc1 = [loc1, length(m1)];
            for i = 1:length(diff2)
                if diff2(i) > q_2
                    loc2 = [loc2, i];
                end
            end
            loc2 = [loc2, length(m2)];
            % complete segmentation
            m1_seg = cell(1,(length(loc1)-1));
            m2_seg = cell(1,(length(loc2)-1));
            for i = 1:(length(loc1)-1)
                m1_seg{1,i} = m1((loc1(i)+1):loc1(i+1),:);
            end
            for i = 1:(length(loc2)-1)
                m2_seg{1,i} = m2((loc2(i)+1):loc2(i+1),:);
            end
            
            % delete very short segments (noise)
            del1 = [];
            for i = 1:(length(loc1)-1)
                a = size(m1_seg{1,i});
                if a(1) < 3
                    del1 = [del1, i];
                end
            end
            m1_seg(:,del1) = [];
            
            del2 = [];
            for i = 1:(length(loc2)-1)
                a = size(m2_seg{1,i});
                if a(1) < 3
                    del2 = [del2, i];
                end
            end
            m2_seg(:,del2) = [];
            
            % calculate dtwd matrix for two segmented trajectories
            inter_result = zeros(length(m1_seg),length(m2_seg));
            for i = 1:(length(m1_seg))
                for j = 1:(length(m2_seg))
                    md1 = m1_seg{1,i};
                    md2 = m2_seg{1,j};
                    ddtw1 = zeros(length(md1),3);
                    ddtw2 = zeros(length(md2),3);
                    for k = 2:length(md1)-1
                        ddtw1(k,:) = (md1(k,:)-md1(k-1,:) + ...
                            (md1(k+1,:)-md1(k-1,:))./2)./2;               
                    end
                    ddtw1(1,:)=ddtw1(2,:);
                    ddtw1(length(md1),:)=ddtw1(length(md1)-1,:);
                    for k = 2:length(md2)-1
                        ddtw2(k,:) = (md2(k,:)-md2(k-1,:) + ...
                            (md2(k+1,:)-md2(k-1,:))./2)./2;               
                    end
                    ddtw2(1,:)=ddtw2(2,:);
                    ddtw2(length(md2),:)=ddtw2(length(md2)-1,:);                   
                    inter_result(i,j) = dtw(ddtw1', ddtw2');
                end
            end

            % delete overlap
            % policy 1
            % calculate dtwd matrix for two segmented trajectories
            inter_result1 = inter_result;

            % find minimal numbers for p1
            loc3 = [];
            sum_seg1 = 0;
            for i = 1:(length(m1_seg))
                [M, I] = min(inter_result1(i,:));
                sum_seg1 = sum_seg1 + M;
                loc3 = [loc3, I];
            end
            loc4 = unique(loc3);
            inter_result1(:,loc4) = [];
            % find minimal numbers for rest p2
            new = size(inter_result1);
            for i = 1:new(2)
                sum_seg1 = sum_seg1 + min(inter_result1(:,i));
            end

            % policy 2
            % calculate dtwd matrix for two segmented trajectories
            inter_result2 = inter_result;

            % find minimal numbers for p2
            loc3 = [];
            sum_seg2 = 0;
            for i = 1:(length(m2_seg))
                [M, I] = min(inter_result2(:,i));
                sum_seg2 = sum_seg2 + M;
                loc3 = [loc3, I];
            end
            loc4 = unique(loc3);
            inter_result2(loc4,:) = [];
            % find minimal numbers for rest p2
            new = size(inter_result2);
            for i = 1:new(1)
                sum_seg2 = sum_seg2 + min(inter_result2(i,:));
            end

            % dtw
            result_ol99min(m,n) = min(sum_seg1/(length(m1)+length(m2)), ...
                sum_seg2/(length(m1)+length(m2)));
            result_ol99max(m,n) = max(sum_seg1/(length(m1)+length(m2)), ...
                sum_seg2/(length(m1)+length(m2)));
        end
    end
end
toc

for m = 1:21
    for n = 1:21
        if m > n
            result_ol99min(m,n) = result_ol99min(n,m);
            result_ol99max(m,n) = result_ol99max(n,m);
        end
    end
end

b = result_ol99min;
% Silhouette
sil_ol99min= zeros(21,1);
for i = 1:21
    x = b(i,:);
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
    sil_ol99min(i) = (c-a)/max(a,c);
end
mean_sil_ol99min = sum(sil_ol99min)/21;
mean_sil_ol99min
% -0.0444

%% WDTW
tic
g=0.01;
result1=zeros(21,21);
for m = 1:21
    for n = 1:21
        if m ~= n
            matFileName1 = sprintf('psz%d', m);
            matFileName2 = sprintf('psz%d', n);
            matData = eval(matFileName1);
            matData2 = eval(matFileName2);
            result1(m,n) = wdtw(matData', matData2',g)/max(length(matData),length(matData2));
        end
    end
end
toc

b = result1;
% Silhouette
sil_dtw = zeros(21,1);
for i = 1:21
    x = b(i,:);
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
    sil_dtw(i) = (c-a)/max(a,c);
end
mean_sil_dtw = sum(sil_dtw)/21;
mean_sil_dtw
% -0.6268

%% SWDTW
tic
g=0.01;
quan = 0.99;
result_ol99min = zeros(21,21);
result_ol99max = zeros(21,21);
for m = 1:21
    for n = 1:21
        if m < n
            m1 = dataset2{1, m};
            m2 = dataset2{1, n};
            % segment when distance > quantile
            diff1=sqrt(sum((diff(m1)).^2, 2));
            diff2=sqrt(sum((diff(m2)).^2, 2));
            % set the quantile
            q_1 = quantile(diff1, quan);
            q_2 = quantile(diff2, quan);
            loc1 = 1;
            loc2 = 1;
            for i = 1:length(diff1)
                if diff1(i) > q_1
                    loc1 = [loc1, i];
                end
            end
            loc1 = [loc1, length(m1)];
            for i = 1:length(diff2)
                if diff2(i) > q_2
                    loc2 = [loc2, i];
                end
            end
            loc2 = [loc2, length(m2)];
            % complete segmentation
            m1_seg = cell(1,(length(loc1)-1));
            m2_seg = cell(1,(length(loc2)-1));
            for i = 1:(length(loc1)-1)
                m1_seg{1,i} = m1((loc1(i)+1):loc1(i+1),:);
            end
            for i = 1:(length(loc2)-1)
                m2_seg{1,i} = m2((loc2(i)+1):loc2(i+1),:);
            end
            
            % delete very short segments (noise)
            del1 = [];
            for i = 1:(length(loc1)-1)
                a = size(m1_seg{1,i});
                if a(1) < 3
                    del1 = [del1, i];
                end
            end
            m1_seg(:,del1) = [];
            
            del2 = [];
            for i = 1:(length(loc2)-1)
                a = size(m2_seg{1,i});
                if a(1) < 3
                    del2 = [del2, i];
                end
            end
            m2_seg(:,del2) = [];
            
            % calculate dtwd matrix for two segmented trajectories
            inter_result = zeros(length(m1_seg),length(m2_seg));
            for i = 1:(length(m1_seg))
                for j = 1:(length(m2_seg))
                    inter_result(i,j) = wdtw(m1_seg{1,i}',m2_seg{1,j}',g);
                end
            end

            % delete overlap
            % policy 1
            % calculate dtwd matrix for two segmented trajectories
            inter_result1 = inter_result;

            % find minimal numbers for p1
            loc3 = [];
            sum_seg1 = 0;
            for i = 1:(length(m1_seg))
                [M, I] = min(inter_result1(i,:));
                sum_seg1 = sum_seg1 + M;
                loc3 = [loc3, I];
            end
            loc4 = unique(loc3);
            inter_result1(:,loc4) = [];
            % find minimal numbers for rest p2
            new = size(inter_result1);
            for i = 1:new(2)
                sum_seg1 = sum_seg1 + min(inter_result1(:,i));
            end

            % policy 2
            % calculate dtwd matrix for two segmented trajectories
            inter_result2 = inter_result;

            % find minimal numbers for p2
            loc3 = [];
            sum_seg2 = 0;
            for i = 1:(length(m2_seg))
                [M, I] = min(inter_result2(:,i));
                sum_seg2 = sum_seg2 + M;
                loc3 = [loc3, I];
            end
            loc4 = unique(loc3);
            inter_result2(loc4,:) = [];
            % find minimal numbers for rest p2
            new = size(inter_result2);
            for i = 1:new(1)
                sum_seg2 = sum_seg2 + min(inter_result2(i,:));
            end

            % dtw
            result_ol99min(m,n) = min(sum_seg1/(length(m1)+length(m2)), ...
                sum_seg2/(length(m1)+length(m2)));
            result_ol99max(m,n) = max(sum_seg1/(length(m1)+length(m2)), ...
                sum_seg2/(length(m1)+length(m2)));
        end
    end
end
toc

for m = 1:21
    for n = 1:21
        if m > n
            result_ol99min(m,n) = result_ol99min(n,m);
            result_ol99max(m,n) = result_ol99max(n,m);
        end
    end
end

b = result_ol99min;
% Silhouette
sil_ol99min= zeros(21,1);
for i = 1:21
    x = b(i,:);
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
    sil_ol99min(i) = (c-a)/max(a,c);
end
mean_sil_ol99min = sum(sil_ol99min)/21;
mean_sil_ol99min
% 0.0085

%% WDDTW
tic
g=0.01;
result1=zeros(21,21);
for m = 1:21
    for n = 1:21
        if m ~= n
            matFileName1 = sprintf('psz%d', m);
            matFileName2 = sprintf('psz%d', n);
            matData = eval(matFileName1);
            matData2 = eval(matFileName2);
            ddtw1 = zeros(length(matData),3);
            ddtw2 = zeros(length(matData2),3);
            for i = 2:length(matData)-1
                ddtw1(i,:) = (matData(i,:)-matData(i-1,:) + ...
                    (matData(i+1,:)-matData(i-1,:))./2)./2;               
            end
            ddtw1(1,:)=ddtw1(2,:);
            ddtw1(length(matData),:)=ddtw1(length(matData)-1,:);
            for i = 2:length(matData2)-1
                ddtw2(i,:) = (matData2(i,:)-matData2(i-1,:) + ...
                    (matData2(i+1,:)-matData2(i-1,:))./2)./2;               
            end
            ddtw2(1,:)=ddtw2(2,:);
            ddtw2(length(matData2),:)=ddtw2(length(matData2)-1,:);
            result1(m,n) = wdtw(ddtw1', ddtw2',g)/max(length(ddtw1),length(ddtw2));
        end
    end
end
toc

b = result1;
% Silhouette
sil_dtw = zeros(21,1);
for i = 1:21
    x = b(i,:);
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
    sil_dtw(i) = (c-a)/max(a,c);
end
mean_sil_dtw = sum(sil_dtw)/21;
mean_sil_dtw
% -0.6727

%% SWDDTW
tic
g=0.01;
quan = 0.99;
result_ol99min = zeros(21,21);
result_ol99max = zeros(21,21);
for m = 1:21
    for n = 1:21
        if m < n
            m1 = dataset2{1, m};
            m2 = dataset2{1, n};
            % segment when distance > quantile
            diff1=sqrt(sum((diff(m1)).^2, 2));
            diff2=sqrt(sum((diff(m2)).^2, 2));
            % set the quantile
            q_1 = quantile(diff1, quan);
            q_2 = quantile(diff2, quan);
            loc1 = 1;
            loc2 = 1;
            for i = 1:length(diff1)
                if diff1(i) > q_1
                    loc1 = [loc1, i];
                end
            end
            loc1 = [loc1, length(m1)];
            for i = 1:length(diff2)
                if diff2(i) > q_2
                    loc2 = [loc2, i];
                end
            end
            loc2 = [loc2, length(m2)];
            % complete segmentation
            m1_seg = cell(1,(length(loc1)-1));
            m2_seg = cell(1,(length(loc2)-1));
            for i = 1:(length(loc1)-1)
                m1_seg{1,i} = m1((loc1(i)+1):loc1(i+1),:);
            end
            for i = 1:(length(loc2)-1)
                m2_seg{1,i} = m2((loc2(i)+1):loc2(i+1),:);
            end
            
            % delete very short segments (noise)
            del1 = [];
            for i = 1:(length(loc1)-1)
                a = size(m1_seg{1,i});
                if a(1) < 3
                    del1 = [del1, i];
                end
            end
            m1_seg(:,del1) = [];
            
            del2 = [];
            for i = 1:(length(loc2)-1)
                a = size(m2_seg{1,i});
                if a(1) < 3
                    del2 = [del2, i];
                end
            end
            m2_seg(:,del2) = [];
            
            % calculate dtwd matrix for two segmented trajectories
            inter_result = zeros(length(m1_seg),length(m2_seg));
            for i = 1:(length(m1_seg))
                for j = 1:(length(m2_seg))
                    md1 = m1_seg{1,i};
                    md2 = m2_seg{1,j};
                    ddtw1 = zeros(length(md1),3);
                    ddtw2 = zeros(length(md2),3);
                    for k = 2:length(md1)-1
                        ddtw1(k,:) = (md1(k,:)-md1(k-1,:) + ...
                            (md1(k+1,:)-md1(k-1,:))./2)./2;               
                    end
                    ddtw1(1,:)=ddtw1(2,:);
                    ddtw1(length(md1),:)=ddtw1(length(md1)-1,:);
                    for k = 2:length(md2)-1
                        ddtw2(k,:) = (md2(k,:)-md2(k-1,:) + ...
                            (md2(k+1,:)-md2(k-1,:))./2)./2;               
                    end
                    ddtw2(1,:)=ddtw2(2,:);
                    ddtw2(length(md2),:)=ddtw2(length(md2)-1,:);                   
                    inter_result(i,j) = wdtw(ddtw1', ddtw2',g);
                end
            end

            % delete overlap
            % policy 1
            % calculate dtwd matrix for two segmented trajectories
            inter_result1 = inter_result;

            % find minimal numbers for p1
            loc3 = [];
            sum_seg1 = 0;
            for i = 1:(length(m1_seg))
                [M, I] = min(inter_result1(i,:));
                sum_seg1 = sum_seg1 + M;
                loc3 = [loc3, I];
            end
            loc4 = unique(loc3);
            inter_result1(:,loc4) = [];
            % find minimal numbers for rest p2
            new = size(inter_result1);
            for i = 1:new(2)
                sum_seg1 = sum_seg1 + min(inter_result1(:,i));
            end

            % policy 2
            % calculate dtwd matrix for two segmented trajectories
            inter_result2 = inter_result;

            % find minimal numbers for p2
            loc3 = [];
            sum_seg2 = 0;
            for i = 1:(length(m2_seg))
                [M, I] = min(inter_result2(:,i));
                sum_seg2 = sum_seg2 + M;
                loc3 = [loc3, I];
            end
            loc4 = unique(loc3);
            inter_result2(loc4,:) = [];
            % find minimal numbers for rest p2
            new = size(inter_result2);
            for i = 1:new(1)
                sum_seg2 = sum_seg2 + min(inter_result2(i,:));
            end

            % dtw
            result_ol99min(m,n) = min(sum_seg1/(length(m1)+length(m2)), ...
                sum_seg2/(length(m1)+length(m2)));
            result_ol99max(m,n) = max(sum_seg1/(length(m1)+length(m2)), ...
                sum_seg2/(length(m1)+length(m2)));
        end
    end
end
toc

for m = 1:21
    for n = 1:21
        if m > n
            result_ol99min(m,n) = result_ol99min(n,m);
            result_ol99max(m,n) = result_ol99max(n,m);
        end
    end
end

b = result_ol99min;
% Silhouette
sil_ol99min= zeros(21,1);
for i = 1:21
    x = b(i,:);
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
    sil_ol99min(i) = (c-a)/max(a,c);
end
mean_sil_ol99min = sum(sil_ol99min)/21;
mean_sil_ol99min
% -0.0937
