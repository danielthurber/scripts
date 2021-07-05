clear; close all;
colors = [
         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
opt.normalize = 0;    
opt.watershed = 1;
switch opt.watershed
    % 1: Logan River, 2: cub, 3: blacksmith, 4: little cottonwood, 
                    % 5: big cottonwood,  6: Red Butte Creek, 7: Parleys
                    % 8: east canyon? 9: Milcreek
    case 1
        load obs_LR
        id_incomplete = [38];
        opt.snotel = 5;
        opt.flow = 4; % 1:low flow, 2: high flow, 3: mean flow, 4: summer (Jul-Sep) flow average, 5: Jan flow
    case 2
        load obs_cub
        id_incomplete = [1,2,9:27,33];
        opt.snotel = 1;
        opt.flow = 4;
    case 3
        load obs_blacksmith
        id_incomplete = [40];
        opt.snotel = 3;
        opt.flow = 4; % 1:low flow, 2: high flow, 3: mean flow, 4: summer (Jul-Sep) flow average, 5: Jan flow
    case 4
        load obs_LittleCotton
        id_incomplete = [26];
        opt.snotel = 1;
        opt.flow = 4;
    case 5
        load obs_BigCotton
        id_incomplete = [29];
        opt.snotel = 1;
        opt.flow = 4;
    case 6
        load obs_RBC 
        id_incomplete = [1:21,40];
        opt.snotel = 3;
        opt.flow = 4;
    case 7
        load obs_parleys
        id_incomplete = [37];
        opt.snotel = 3;
        opt.flow = 4;
    case 8
        load obs_east_canyon
        id_incomplete = [15 16];
        opt.snotel = 4;
        opt.flow = 4;
    case 9
        load obs_millcreek
        id_incomplete = 29;
        opt.snotel = 2;
        opt.flow = 4;
end

swe_year = swe_year(:,:,opt.snotel);
data_snow = data_snow(:,opt.snotel);
years(id_incomplete) = [];
swe_year(:,id_incomplete) = [];
Q_year(:,id_incomplete) = [];
p_year(:,id_incomplete) = [];
T_year(:,id_incomplete,:) = [];
data_Q(id_incomplete,:) = [];
data_snow(id_incomplete,:) = [];
data_p(id_incomplete,:) = [];
data_T(id_incomplete,:) = [];

% display(['cov: lowflow=',num2str(nanstd(data_Q(:,1))/nanmean(data_Q(:,1)))]);
% display(['cov: summer=',num2str(nanstd(data_Q(:,4))/nanmean(data_Q(:,4)))]);
display(['std/mu: lowflow=',num2str(nanstd(data_Q(:,1))/nanmean(data_Q(:,1)))]);
display(['std/mu: summer=',num2str(nanstd(data_Q(:,4))/nanmean(data_Q(:,4)))]);
display(['std/mu: annual=',num2str(nanstd(data_Q(:,3))/nanmean(data_Q(:,3)))]);
display(['std/mu: annual precip=',num2str(nanstd(data_p)/nanmean(data_p))]);

%% Q_t vs Q_t-1
Q_id = 1;  % 1:low flow, 2: high flow, 3: mean flow, 4: summer (Jul-Sep) flow average, 5: Jan flow
figure; hold on
cmap = colormap;
for i = 1:size(data_Q,1)-1
    plot(data_Q(i,Q_id),data_Q(i+1,Q_id),'o','MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:));
    text(data_Q(i,Q_id)+2,data_Q(i+1,Q_id)-2, num2str(i));
end

%% get monthly statistics
p_month = zeros(12,length(years));
Q_month = p_month;
swe_month = p_month;
for j = 1:length(years)
    for i = 1:12
        mm = mod(i+9,12);
        if mm == 0
            mm = 12;
        end
        if mm > 9
            year = years(j);
        else
            year = years(j) + 1;
        end
        id1 = datenum([num2str(mm),'/01/',num2str(year)]) - datenum(['10/01/',...
            num2str(years(j))]) + 1;
        id2 = id1 - 1 + eomday(year,mm);
        p_month(i,j) = p_year(id2,j) - p_year(id1,j);
        Q_month(i,j) = mean(Q_year(id1:id2,j));
        swe_month(i,j) = mean(swe_year(id1:id2,j));
    end
end

%% plot multi-year average monthly flow (signature)
figure;
plot(1:12,nanmean(Q_month,2));
xticks(1:12);
xticklabels({'10','11','12','1','2','3','4','5','6','7','8','9'});

%% remove trend and average seasonal signal
var_name = {'Q_month','p_month','swe_month'};
for i = 1:length(var_name)
    eval(['y = mean(',var_name{i},');']);
    id = ~isnan(y);
    X = [ones(length(years(id)),1), years(id)'];
    b = X\y(id)';
    X = [ones(length(years),1), years'];
    yhat = (X * b)';
    yhat = ones(12,1) * yhat;
    eval([var_name{i}, ' = ', var_name{i}, ' - yhat - mean(', var_name{i}, ',2,''omitnan'') * ones(1,length(years));']);
end

%% ACF and cross ACF
lags = 0:48;
rho = zeros(3,length(lags)); % ACF of Q, cross Q & P, cross Q & SWE
Q_month = Q_month(:);
p_month = p_month(:);
swe_month = swe_month(:);
for i = 1:length(lags)
    tmp = corrcoef(Q_month(1:end-lags(i)),Q_month(lags(i)+1:end),'Rows','complete');
    rhos(1,i) = tmp(1,2);
    tmp = corrcoef(p_month(1:end-lags(i)),Q_month(lags(i)+1:end),'Rows','complete');
    rhos(2,i) = tmp(1,2);
    tmp = corrcoef(swe_month(1:end-lags(i)),Q_month(lags(i)+1:end),'Rows','complete');
    rhos(3,i) = tmp(1,2);
end
figure('Position',[0 0 400 700]);
for i = 1:3
    subplot(3,1,i);
    bar(lags,rhos(i,:));
    axis([-1 max(lags)+1 -0.3 1]);
end
clear Q_month p_month swe_month

%% dry vs wet plot
% pick dry and wet years
swe_annual = nanmax(swe_year);
[~,id] = sort(swe_annual);
id_dry = id(1:5);
id_wet = id(end-4:end);

swe_dry = mean(swe_year(:,id_dry),2);
Q_dry = mean(Q_year(:,id_dry),2);

swe_wet = mean(swe_year(:,id_wet),2);
Q_wet = mean(Q_year(:,id_wet),2);

% plot
figure('Position',[0 0 800 300]);
subplot(1,2,1);
plot((1:366)/30.4+10,swe_wet);
hold on
plot((1:366)/30.4+10,swe_dry);
ylabel('SWE [mm]');
xlim([10 22]);
xticks(10:22);
xticklabels({'10','11','12','1','2','3','4','5','6','7','8','9'});
xlabel('Month');
legend('Wet','Dry','Location','NorthWest');
legend boxoff

subplot(1,2,2);
plot((1:366)/30.4+10,Q_wet);
hold on
plot((1:366)/30.4+10,Q_dry);
ylabel('Discharge [mm/d]');
xlim([10 22]);
xticks(10:22);
xticklabels({'10','11','12','1','2','3','4','5','6','7','8','9'});
xlabel('Month');
legend('Wet','Dry','Location','NorthWest');
legend boxoff

%% correlation SWE, P & Q
opt.lag = [0]; 

figure('Position',[0 0 800 350*length(opt.lag)]);
for i = 1:length(opt.lag)
    for j = 1:2
        lag = opt.lag(i);
        y = data_Q(1+lag:end,4);
        y2 = data_Q(2+lag:end,1);
        switch j
            case 1
                x = data_snow(1:end-lag);
                x2 = data_snow(1:end-lag-1);
            case 2
                x = data_p(1:end-lag);
                x2 = data_p(1:end-lag-1);
        end
        % remove nans for linear regression
        id_keep = ~isnan(x) & ~isnan(y);
        x = x(id_keep);
        y = y(id_keep);

        id_keep = ~isnan(x2) & ~isnan(y2);
        x2 = x2(id_keep);
        y2 = y2(id_keep);   

        if opt.normalize == 1
            x = x / mean(x) * 100;
            y = y / mean(y) * 100;
            x2 = x2 / mean(x2) * 100;
            y2 = y2 / mean(y2) * 100;
        end
        
        x_max = max(x) + range(x)/10;

        tmp_corr = corrcoef(x,y);
        rho1 = tmp_corr(1,2);
        tmp_corr = corrcoef(x2,y2);
        rho2 = tmp_corr(1,2);

        subplot(length(opt.lag),2,2*i-2+j);
        scatter(x,y,[],colors(1,:));
        hold on
        scatter(x2,y2,[],colors(2,:));

        % fit a straight line
        X = [ones(length(x),1), x];
        b = X\y;
        beta = b(2);
        X = [ones(2,1), [0;x_max]];
        yhat = X*b;
        plot([0;x_max], yhat,'-.', 'color',colors(1,:));

        % fit a straight line
        X = [ones(length(x2),1), x2];
        b = X\y2;
        beta2 = b(2);
        X = [ones(2,1), [0;x_max]];
        yhat2 = X*b;
        plot([0;x_max], yhat2,'-.', 'color',colors(2,:));

    %     if b(1) > 0
    %         eq_str = sprintf('Q=%.2f*SWE+%.2f',b(2),b(1));
    %     else
    %         eq_str = sprintf('Q=%.2f*SWE%.2f',b(2),b(1));
    %     end
    %     text(min(x),max(y),eq_str);

        text(max(x)-range(x)/5,yhat(end)+range(y)/10,sprintf('\\rho=%.2f,\\beta=%.2f',rho1,beta));
        text(max(x)-range(x)/5,yhat2(end)-range(y)/10,sprintf('\\rho=%.2f,\\beta=%.2f',rho2,beta2));
        box on;
        xlim([0 x_max]);
        ylabel('Discharge [mm/year]');
        switch j
            case 1
                xlabel('Max SWE [mm]');
            case 2
                xlabel('Precipitation [mm]');
        end
        legend('Summer flow','Low flow','Location','NorthWest');
        legend boxoff
    end
end

%% DOY melt out, DOY max SWE, Melt rate, mean spring T
lag = 0;
% melt out day
tmelt = zeros(length(years),1);
for i = 1:length(years)
    tmp1 = swe_year(1:end-3,i);
    tmp2 = swe_year(2:end-2,i);
    tmp3 = swe_year(3:end-1,i);
    tmp4 = swe_year(4:end,i);
    id = find( tmp1 > 0 & tmp2 > 0 & tmp3 > 0 & tmp4 == 0);
    if length(id) > 1
        id2 = find(id >= 200 & id < 300);
        tmelt(i) = id(id2(1)) + 1 - 92;
    else
        tmelt(i) = id + 1 - 92;
    end
end

% max snowpack day
tmax = zeros(length(years),1);
for i = 1:length(years)
    [~,tmax(i)] = nanmax(swe_year(:,i));
    tmax(i) = tmax(i) - 92;
end

% melt rate
rmelt = zeros(length(years),1);
for i = 1:length(years)
    dswe = swe_year(2:end,i) - swe_year(1:end-1,i);
    id = find(dswe < 0);
    rmelt(i) = - sum(dswe(id))/length(id);
end

figure('Position',[0 0 800 700]);
for j = 1:4
    y = data_Q(1+lag:end,4);
    y2 = data_Q(2+lag:end,1);
    switch j
        case 1
            x = tmelt(1:end-lag);
            x2 = tmelt(1:end-lag-1);
        case 2
            x = tmax(1:end-lag);
            x2 = tmax(1:end-lag-1);
        case 3
            x = rmelt(1:end-lag);
            x2 = rmelt(1:end-lag-1);
        case 4
            x = data_T(1:end-lag,2);
            x2 = data_T(1:end-lag-1,2);
    end
    
    % remove nans for linear regression
    id_keep = ~isnan(x) & ~isnan(y);
    x = x(id_keep);
    y = y(id_keep);
    
    id_keep = ~isnan(x2) & ~isnan(y2);
    x2 = x2(id_keep);
    y2 = y2(id_keep);  
    
    if opt.normalize == 1
        x = x / mean(x) * 100;
        y = y / mean(y) * 100;
        x2 = x2 / mean(x2) * 100;
        y2 = y2 / mean(y2) * 100;
    end
    x_min = min(x) - range(x)/10;
    x_max = max(x) + range(x)/10;
    tmp_corr = corrcoef(x,y);
    rho = tmp_corr(1,2);
    
    if length(x) > 5
        subplot(2,2,j);
        scatter(x,y,[],colors(1,:));
        hold on
        scatter(x2,y2,[],colors(2,:));

        % fit a straight line
        X = [ones(length(x),1), x];
        b = X\y;
        if opt.normalize == 1
            beta = b(2);
        else
            beta = b(2) / mean(y) * mean(x);
        end
        X = [ones(2,1), [x_min;x_max]];
        yhat = X*b;
        plot([x_min;x_max], yhat,'-.', 'color',colors(1,:));

        % fit a straight line
        X = [ones(length(x2),1), x2];
        b = X\y2;
        if opt.normalize == 1
            beta = b(2);
        else
            beta = b(2) / mean(y) * mean(x);
        end
        X = [ones(2,1), [x_min;x_max]];
        yhat2 = X*b;
        plot([x_min;x_max], yhat2,'-.', 'color',colors(2,:));

        tmp_corr = corrcoef(x,y);
        rho1 = tmp_corr(1,2);
        tmp_corr = corrcoef(x2,y2);
        rho2 = tmp_corr(1,2);    
        text(max(x)-range(x)/5,yhat(end)+range(y)/10,sprintf('\\rho=%.2f,\\beta=%.2f',rho1,beta));
        text(max(x)-range(x)/5,yhat2(end)-range(y)/10,sprintf('\\rho=%.2f,\\beta=%.2f',rho2,beta2));

        box on;
        xlim([x_min x_max]);
        ylabel('Discharge [mm/year]');
        switch j
            case 1 
                xlabel('DOY melt out');
            case 2
                xlabel('DOY max SWE');
            case 3
                xlabel('Melt rate [mm/d]');
            case 4
                xlabel('Spring average daily maximum temperature [C]');
        end
        legend('Summer flow','Low flow','Location','NorthWest');
        legend boxoff
    end
end

%% trend
figure('Position',[0 0 1200 700]);
subplot(2,3,1);
plot(years,data_snow,'s');
xlabel('Year');
ylabel('Max SWE [mm]');
% fit a straight line
fit1 = fitlm(years,data_snow);
p = coefTest(fit1);
ypred = predict(fit1,years');
hold on
plot(years, ypred, '-r');
text(years(end-4),ypred(end)+range(ypred)/5,sprintf('p=%.2f',p));

subplot(2,3,2);
plot(years,data_p,'s');
xlabel('Year');
ylabel('Precipitation [mm]');
% fit a straight line
fit1 = fitlm(years,data_p);
p = coefTest(fit1);
ypred = predict(fit1,years');
hold on
plot(years, ypred, '-r');
text(years(end-4),ypred(end)+range(ypred)/5,sprintf('p=%.2f',p));

tmp0 = [1 2 4];
for i = 1:length(tmp0)
    opt.flow = tmp0(i);
    subplot(2,3,3+i);
    tmp = data_Q(:,opt.flow);
    plot(years,tmp,'s');
    xlabel('Year');
    switch opt.flow 
        case 1
            ylabel('Low flow [mm]');
        case 2
            ylabel('High flow [mm]');
        case 3
            ylabel('Mean flow [mm]');
        case 4
            ylabel('Summer average flow [mm]');      
        case 5
            ylabel('Winter average flow [mm]');
    end
    % fit a straight line
    fit1 = fitlm(years,tmp);
    p = coefTest(fit1);
    ypred = predict(fit1,years');
    hold on
    plot(years, ypred, '-r');
    text(years(end-4),ypred(end)+range(ypred)/5,sprintf('p=%.2f',p));
end

%% wet periods
% figure('Position',[0 0 1000 300]);
% subplot(1,3,1); 
% hist(data_snow);
% subplot(1,3,2);
% hist(data_p);
% subplot(1,3,3);
% hist(data_Q(:,opt.flow));
% 
% flag_wet = nan(size(data_snow));
% flag_wet(data_snow > quantile(data_snow,0.75)) = 1;
% flag_wet(data_snow < 800) = 0;
% 
% id_dry = [3 9];
% id_wet = [4:8];
% 
% swe_wet = swe_year(:,id_wet);
% swe_dry = swe_year(:,id_dry);
% Q_wet = Q_year(:,id_wet);
% Q_dry = Q_year(:,id_dry);
% 
% % smoothing, default window = 5
% for i = 1:size(Q_wet,2)
%     Q_wet(:,i) = smooth(Q_wet(:,i));
% end
% for i = 1:size(Q_dry,2)
%     Q_dry(:,i) = smooth(Q_dry(:,i));
% end
% 
% figure('Position',[0 0 800 300]);
% subplot(1,2,1);
% plot((1:366)/30.4+10,swe_wet);
% hold on
% plot((1:366)/30.4+10,swe_dry);
% ylabel('SWE [mm]');
% xlim([10 22]);
% xticks(10:22);
% xticklabels({'10','11','12','1','2','3','4','5','6','7','8','9'});
% xlabel('Month');
% 
% subplot(1,2,2);
% plot((1:366)/30.4+10,Q_wet);
% hold on
% plot((1:366)/30.4+10,Q_dry);
% ylabel('Discharge [mm/d]');
% xlim([10 22]);
% xticks(10:22);
% xticklabels({'10','11','12','1','2','3','4','5','6','7','8','9'});
% xlabel('Month');

%% droughts
% id_dry = [12,13,14];
% id_wet = [11 15];
% 
% swe_wet = swe_year(:,id_wet);
% swe_dry = swe_year(:,id_dry);
% Q_wet = Q_year(:,id_wet);
% Q_dry = Q_year(:,id_dry);
% 
% % smoothing, default window = 5
% for i = 1:size(Q_wet,2)
%     Q_wet(:,i) = smooth(Q_wet(:,i));
% end
% for i = 1:size(Q_dry,2)
%     Q_dry(:,i) = smooth(Q_dry(:,i));
% end
% 
% figure('Position',[0 0 800 300]);
% subplot(1,2,1);
% plot((1:366)/30.4+10,swe_wet);
% hold on
% plot((1:366)/30.4+10,swe_dry);
% ylabel('SWE [mm]');
% xlim([10 22]);
% xticks(10:22);
% xticklabels({'10','11','12','1','2','3','4','5','6','7','8','9'});
% xlabel('Month');
% 
% subplot(1,2,2);
% plot((1:366)/30.4+10,Q_wet);
% hold on
% plot((1:366)/30.4+10,Q_dry);
% ylabel('Discharge [mm/d]');
% xlim([10 22]);
% xticks(10:22);
% xticklabels({'10','11','12','1','2','3','4','5','6','7','8','9'});
% xlabel('Month');

%% memory effect
figure('Position',[0 0 800 700]);
for j = 4:5 % 1:low flow, 2: high flow, 3: mean flow, 4: summer (Jul-Sep) flow average, 5: Jan flow
    opt.flow = j;
    % first check out autocorrelation
    lags = 0:1:5;
    rhos = zeros(size(lags));
    for i = 1:length(lags)
        tmp = corrcoef(data_Q(1:end-lags(i),opt.flow),data_Q(lags(i)+1:end,opt.flow),'Rows','complete');
        rhos(i) = tmp(1,2);
    end
    subplot(2,2,(j-4)*2+1);
    bar(lags,rhos);
    subplot(2,2,(j-4)*2+2);
    for i = 1:length(opt.lag)
        switch opt.flow 
            case 1
                lag = 1;
            case 4
                lag = 0;
        end
        tmpy = data_Q(1+lag:end,opt.flow);
        tmpx = data_snow(1:end-lag);

        id_dry = find(data_snow < mean(data_snow));
        id_dry_prev = id_dry + 1;
        id_dry_prev(id_dry_prev > length(tmpx)) = [];
        id_wet_prev = setdiff(1:length(tmpx),id_dry_prev);

        x2 = tmpx(id_dry_prev);
        y2 = tmpy(id_dry_prev);
        x = tmpx(id_wet_prev);
        y = tmpy(id_wet_prev);

        % remove nans for linear regression
        id_keep = ~isnan(x) & ~isnan(y);
        x = x(id_keep);
        y = y(id_keep);

        id_keep = ~isnan(x2) & ~isnan(y2);
        x2 = x2(id_keep);
        y2 = y2(id_keep);   

        x_max = max(x) + range(x)/10;
        tmp_corr = corrcoef(x,y);
        rho = tmp_corr(1,2);

        scatter(x,y,[],colors(1,:));
        hold on
        scatter(x2,y2,[],colors(2,:));

        % fit a straight line
        X = [ones(length(x),1), x];
        b = X\y;
        X = [ones(2,1), [0;x_max]];
        yhat = X*b;
        plot([0;x_max], yhat,'-.', 'color',colors(1,:));

        % fit a straight line
        X = [ones(length(x2),1), x2];
        b = X\y2;
        X = [ones(2,1), [0;x_max]];
        yhat2 = X*b;
        plot([0;x_max], yhat2,'-.', 'color',colors(2,:));

        if b(1) > 0
            eq_str = sprintf('Q=%.2f*SWE+%.2f',b(2),b(1));
        else
            eq_str = sprintf('Q=%.2f*SWE%.2f',b(2),b(1));
        end
    %     text(min(x),max(y),eq_str);
        box on;
        xlim([0 x_max]);
        ylabel('Discharge [mm/year]');
        xlabel('Max SWE [mm]');
        legend('Previous year above average','Previous year below average','Location','NorthWest');
        legend boxoff
    end
end

