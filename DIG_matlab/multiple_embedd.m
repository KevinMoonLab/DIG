%% Multiple 
load('n10_8-12hz.mat');

series_distancev = ["EIG"];
ndimv = [2];
tv = [10];
iegv = [0];
pot_methodv = ['dm'];
alphav = [20];
maxfeature = 1:1;
gammav = [0];
% 'L', 3840

i = 1;


for gamma = gammav
j = 1;    
for diff_dis = iegv
for series_distance = series_distancev 
for t = tv
for alpha = alphav
[Y, z_hist] = DIG(data(4:8,:), 'L',3840, 't', t, 'ieg', diff_dis, 'pot_method', 'gamma', 'series_distance', ... 
    series_distance, 'ndim', ndimv, 'a', alpha, 'gamma', gamma, 'weighted', 1);
% LabelGraph{j,i} = "SeriesDistance: " + series_distance + newline ...
%       + "DiffDistance: " + "GAMMA" + newline + "alpha: " + alpha + " and t: " + t + newline + ...
%       "Channels: " + "1" + "to " + 10 +  "Gamma: " + gamma; 

end
j=j+1;
end
end 
end 
i=i+1;
end 

%% 2D plots
PlotPFig = 6;
NumberOfGraphs = ceil(length(Y)/PlotPFig);
% time progress
% color = linspace(1,10,length(hyp));
color = 1:2000;

for Nf = 1:NumberOfGraphs
figure(Nf)
for sub = 1:PlotPFig
subplot(2,3, sub);
color = 1:length(Y{sub ,i});
scatter(Y{sub ,i}(:,1), Y{sub ,i}(:,2), 10, color, 'filled');
%scatter(Y{sub+PlotPFig*(Nf-1)}(:,1), Y{sub+PlotPFig*(Nf-1)}(:,2), 10, hyp, 'filled');
%title(LabelGraph{sub+PlotPFig*(Nf-1)})
colorbar;
colormap(jet); 
end 
end 

%% Filter data to just 1 sleep stage


newdata{3} = data(1:10,:);
L=3840;
hyp2{1}=hyp;
index=[];
k=1;
for i=1:860
if (hyp(i) ~= 5)
index = [index i];
newdata{3}(:,1+(k-1)*L:k*L) = []; 
k = k-1;
end 
i
k = k+1;
end 


%% Video 
l =  length(Y);

% color = zeros(1,860);
% color(610) = 1;
i=1;
for j=1:6
    figure(j)
    color = 1:length(Y{j,i});
scatter(Y{j,i}(:,1), Y{j,i}(:,2), 10, color, 'filled');
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'fontname','Times New Roman','FontWeight','bold')
%title(LabelGraph{i})
colorbar;
colormap(parula); 
drawnow;
pause(1.5)

end 
%% Video gscatter
l =  length(Y);
color = 1:length(Y{1});
% color = zeros(1,860);
% color(610) = 1;
j=1;
for i=1:13
    figure(i)
gscatter(Y{j,i}(:,1), Y{j,i}(:,2), hyp5,parula(4),".ox*",[15 5 5 5], 'off');
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'fontname','Times New Roman','FontWeight','bold')
%title(LabelGraph{i})
%colorbar;
% colormap(jet); 
% colormap([0 0 .5; 0 1 0; 1 1 0; 1 0 0])
drawnow;
pause(1.5)

end 
%%
for ii=1:N
    text(Y{i}(ii),Y{i}(ii),textCell(ii))
end

%% Video TrustW

l =  length(T);
for i=1:1
figure(i)  

imagesc(T{i})

% c = gray;
% c = flipud(c);
% colormap(c);

%colormap(jet);
%caxis([0.8 1]);
xticks([1:13])
yticks([1:13])
xticklabels({'DM','-1','-0.8','-0.6','-0.4','-0.2','0','0.2','0.4','0.6','0.8','PHATE','Geo. Distance'})   
yticklabels({'DM','-1','-0.8','-0.6','-0.4','-0.2','0','0.2','0.4','0.6','0.8','PHATE','Geo. Distance'})
xtickangle(45)
%h = colorbar('XTick', linspace(min(T{i},[],'all'), max(T{i},[],'all'),5));
daspect([1 1 1])
% h = colorbar('XTick', linspace(0.85, 1,5));
caxis([0.9 1])
% tix=get(h,'XTick')';
% set(h,'xticklabel',num2str(tix,'%.2f'))
if i ==3
h = colorbar('XTick', linspace(0.9, 1,5));
tix=get(h,'XTick')';
set(h,'xticklabel',num2str(tix,'%.2f'))   
end 

% Create ylabel
if i ==1
ylabel({'Gamma'});
else
ylabel({''});
% Create xlabel
end 
xlabel({'Gamma'});
set(gca,'fontname','Times New Roman','FontWeight','bold')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times New Roman','fontsize',9)
drawnow;
pause(2.5)
  
end 
%%


index = (hyp ~= 3 & hyp ~= 4);
Y2 = Y;
for j = 1:size(Y2,1)
for i = 1:size(Y2,2)

Y2{j,i}(index,:) = [];

end 
end 
%% subset the data


index = (hyp ~= 2);
Y5 = Y2;

for i=1:length(Y)

Y5{i}(index,:) = [];

end 

%% 3D

PlotPFig = 1;
NumberOfGraphs = ceil(length(Y)/PlotPFig);
% time progress
color = linspace(1,10,length(hyp));


for Nf = 1:NumberOfGraphs
figure(Nf)
for sub = 1:6
subplot(2,3, sub);
scatter3(Y{sub+6*(Nf-1)}(:,1), Y{sub+6*(Nf-1)}(:,2), Y{sub+6*(Nf-1)}(:,3), 10, hyp, 'filled');
colormap(jet); 
title(LabelGraph{sub+6*(Nf-1)})
end 
end 

%% Multiple EIG bins and cov windows



series_distancev = ["IEG"];
ndimv = [2];
tv = [10];
iegv = [0];
pot_methodv = ['dm'];
alphav = [10];
maxfeature = 1:1;
bins = [5 10 50];
covw = [10 15]; 


i = 1;
for diff_dis = iegv
for series_distance = series_distancev 
for t = tv
for alpha = alphav
    for maxF = 2:2
    for nbin = bins
    for ncov = covw
        data = data2(1:maxF,:);
Y{i} = main(data, 'L', 3840, 't', t, 'ieg', diff_dis, 'pot_method', 'log', 'series_distance', series_distance, 'ndim', ndimv, 'a', alpha, ...
       'nbin', nbin, 'ncov', ncov);
LabelGraph{i} = "SeriesDistance: " + series_distance + newline ...
      + "DiffDistance: " + "PHATE" + newline + "alpha: " + alpha + " and t: " + t + newline + ...
      "Channels: " + "1" + "to " + maxF; 
i = i+1;
    end 
    end 
    end 
end
end
end 
end 

%% Video multiple emb vs individual
l =  length(Y);
color = 1:length(Y{1});
color2=(hyp ~= 1 & hyp ~= 2);
colorMap = zeros(length(color2),3);
for k = 1 : length(color2)
  if color2(k) == 1
    colorMap(k,:) = [0 0 1]; % Red
  else
    colorMap(k,:) = [1 0 0]; % Blue
  end
end
% color = zeros(1,860);
% color(610) = 1;
j=3;
for i=1:l
figure(i)
subplot(1,2,1)
scatter(Y{j,i}(:,1), Y{j,i}(:,2), 10, colorMap, 'filled');
subplot(1,2,2)
scatter(Y2{j,i}(:,1), Y2{j,i}(:,2), 10, "red", 'filled');
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'fontname','Times New Roman','FontWeight','bold')
%title(LabelGraph{i})
drawnow;
pause(1.5)

end 


%% trim the data

data_trim = data_ori(2:5, (data_ori(5,:) < 3*var(data_ori(5,:)) & data_ori(5,:) > -3*var(data_ori(5,:)))...
    & (data_ori(2,:) < 3*var(data_ori(2,:)) & data_ori(2,:) > -3*var(data_ori(2,:))) ... 
    & (data_ori(3,:) < 3*var(data_ori(3,:)) & data_ori(3,:) > -3*var(data_ori(3,:)))...
    & (data_ori(4,:) < 3*var(data_ori(4,:)) & data_ori(4,:) > -3*var(data_ori(4,:))));

data_trimN =  data_trim + normrnd(0,3,[4,length(data_trim)]);
%data_trimN =  data_trim + 6*rand(4,length(data_trim));
