clear all;
% close all;

% ------------ Image directory  -------------------------%
%%
im_dir = 'D:\Dropbox\Dropbox\Acads (UP)\5th Yr  (1st Sem 2019)\AP 186\Act 12\test_zip\apple\'; roi = imread('roi_a.jpg'); leg = 'x'; %apple
% im_dir = 'D:\Dropbox\Dropbox\Acads (UP)\5th Yr  (1st Sem 2019)\AP 186\Act 12\test_zip\banana\'; roi = imread('roi_b.jpg'); leg = 'd'; %banana
% im_dir = 'D:\Dropbox\Dropbox\Acads (UP)\5th Yr  (1st Sem 2019)\AP 186\Act 12\test_zip\orange\';roi = imread('roi_o.jpg'); leg = '*'; %orange

im_list = 'All';        % Images or frames to read. 1st image sets reference plane

% im_dir = 'D:\Kim\Pictures\FreeVideoToJPGConverter\P1 (10-15-2019 6-24-24 PM)\'; %Projectile
% im_list = [22:29]; %Projectile
% tthres = 0.095;
% case_ = 2; 
% % filename = 'Projectile.gif';

% ------------ Read images and conversion ------------------------------%
%%
files = dir(im_dir);
if ischar(im_list)       
    im_list = 1:(length(files)-2);	% Create a list of all images in directory
end
num_images = length(im_list);
for i = 1:num_images
%     im(i).image = rgb2gray(imread(strcat(im_dir, files(im_list(i)+2).name))); 
    im(i).image = (imread(strcat(im_dir, files(im_list(i)+2).name))); 
    im(i).gray = rgb2gray(im(i).image);
    im(i).double = double(im(i).image);
    
    im(i).bw = im2bw(im(i).image);
    im(i).bw = ~im(i).bw;
    im(i).bw = 1- im(i).bw;
    im(i).bw = (im(i).bw == 0);

%     figure(1);
%     imshow(im(i).image);

%     im(i).name = files(im_list(i)+2).name;
end


% ------------ Variables --------------------------------%
%%

% 
norm = 1000000; %To avoid division by zero for NCC
% size_x = size(im(1).image,1);
% size_y = size(im(1).image,2);



roi = double(imresize(roi, [size(im(1).image,1) size(im(1).image,2)]));
% imshow(roi);

% ------------ L*a*b Extraction --------------------------------%
%%

for i = 1:num_images
    im(i).lab = rgb2lab(im(i).image);
    im(i).L = mean(mean(im(i).lab(:,:,1)));
    im(i).a = mean(mean(im(i).lab(:,:,2)));
    im(i).b = mean(mean(im(i).lab(:,:,3)));
    im(i).Hue = (abs(atan(im(i).b/im(i).a)));

end


% ------------ NCC convertion ------------------------------%
%%

for i = 1:num_images
%     im(i).image = rgb2gray(imread(strcat(im_dir, files(im_list(i)+2).name))); 
    
    
    im(i).NCC = im(i).double(:,:,1) + im(i).double(:,:,2) + im(i).double(:,:,3);
    im(i).NCC(find(im(i).NCC==0))=norm;
    
    im(i).r_img = im(i).double(:,:,1)./im(i).NCC;
    im(i).g_img = im(i).double(:,:,2)./im(i).NCC;
%     figure(1);
%     imshow(im(i).image);

%     im(i).name = files(im_list(i)+2).name;
end

% imshow(im(1).image);




% ------------ Segmentation ------------------------------%
%%

% % //Convert to normalized chromaticity coordinates
I=roi(:,:,1)+roi(:,:,2)+roi(:,:,3);
% I(find(I==0))=1000000;
I(find(I==0))=norm;
r=roi(:,:,1)./I;
g=roi(:,:,2)./I;
% 
% //2D Histogram generation
BINS=32;
rint=round(r*(BINS-1)+1);
gint=round(g*(BINS-1)+1);

% FOR NCC Graph
colors=gint(:)+(rint(:)-1)*BINS;

ncc = zeros(BINS, BINS, 3);
for row = 1:1:BINS
    for col = 1:1:(BINS - row + 1)
        ncc(row,col,:) = [round(row*255/BINS) round(col*255/BINS) (255 - round(row*255/BINS) - round(col*255/BINS))];
%         ncc(row,col)=length(find(colors==(((col_+(row_-1)*BINS)))));
    end
end


% ------------ Parametric Segmentation ------------------------------%
%%
% Parametric
% for i = 1:num_images
%     im(i).p_r =(1 ./(std(r)*sqrt(2*pi))).*exp(-((im(i).r_img -mean(r)) ./ (2*std(r))).^2);
%     im(i).p_g =(1 ./(std(g)*sqrt(2*pi))).*exp(-((im(i).g_img -mean(g)) ./ (2*std(g))).^2);
%     im(i).para = im(i).p_r.*im(i).p_g;
% %     imshow(im(i).para);
% end

% ------------ Non-Parametric Segmentation ------------------------------%
%%
% 
hist=zeros(BINS,BINS);
for row = 1:BINS
    for col = 1:(BINS-row+1) 
        hist(row,col)=length(find(colors==(((col+(row-1)*BINS)))));
    end
end
hist=mat2gray(hist);

% Non-Parametric
for i = 1:num_images
    im(i).r_int_img = round(im(i).r_img*(BINS-1)+1);
    im(i).g_int_img = round(im(i).g_img*(BINS-1)+1);
    NP(i).n = im(i).r_int_img + (im(i).g_int_img -1).*size(hist,1);
    NP(i).hist_list = hist(NP(i).n);
    NP(i).nonpara = reshape(NP(i).hist_list,size(im(i).r_img,1),size(im(i).r_img,2));
    NP(i).nonpara(find(0<NP(i).nonpara<0.5))=1;
    
    %make values into 1 and 0's
    im(i).bw_ = im2bw(NP(i).nonpara);
    im(i).bw_ = ~im(i).bw_;
    im(i).bw_ = 1- im(i).bw_;
    im(i).bw_ = (im(i).bw_ == 0);
%     imshow(NP(i).nonpara);
end

% ------------ Threshold Segmentation ------------------------------%
%%
% threshold = 20;
% num_max = 0;
% tthres = 0.5;
% 
% for i = 1:num_images
%   [im(i).count, im(i).cells] = imhist(im(i).image);
%   im(i).bw = im2bw(im(i).image,tthres);
%   [bw(i).L bw(i).num] = bwlabel(im(i).bw);
%   
% %   num_max = bw(i).num;
%   if bw(i).num > num_max;
%       num_max = bw(i).num;
%   end
% %   bw(i).counts = sum(bsxfun(@eq,bw(i).L(:),1:max(bw(i).num)));
% 
%   bw(i).counts = sum(bsxfun(@eq,bw(i).L(:),1:num_max));
%   bw(i).B1 = bsxfun(@eq,bw(i).L,permute(find(bw(i).counts>threshold),[1 3 2]));
%   stats(i).stats = regionprops(bw(i).B1,'Centroid','MajorAxisLength','Area','MinorAxisLength','BoundingBox','Area');
%   
% %   conv_fy = stats(1).stats.BoundingBox(4)/(2*r); %Pixel/mm , Radius of pinpong is 40 mm
% %   conv_fx = stats(1).stats.BoundingBox(3)/(2*r); %Pixel/mm , Radius of pinpong is 40 mm
%   
%   %Adjust the centroid from the smudge
% %   if i > 1 && case_ == 1
% %       stats(i).stats.Centroid(2) = stats(i).stats.Centroid(2) + (stats(i).stats.BoundingBox(4)/2) - (stats(1).stats.BoundingBox(4)/2);
% %   elseif i > 1 && case_ == 2
% %       stats(i).stats.Centroid(1) = stats(i).stats.Centroid(1) + (stats(i).stats.BoundingBox(4)/2) - (stats(1).stats.BoundingBox(4)/2);
% %       stats(i).stats.Centroid(2) = stats(i).stats.Centroid(2) + (stats(i).stats.BoundingBox(4)/2) - (stats(1).stats.BoundingBox(4)/2);
% %       
% %   end
% %   stats(i).xcoords = stats(i).stats.Centroid(1)*(1/conv_fx);
% %   stats(i).ycoords = -(stats(i).stats.Centroid(2)-size(im(1).image,1))*(1/conv_fy);
% %   imshow(bw(i).B1);
% end


% ------------ Blob analysis ------------------------------%
%%

for i = 1:num_images
% for i = 1:2
    temp = [];
    im(i).intstats = regionprops(im(i).bw_,'Centroid','MajorAxisLength','Area','MinorAxisLength','BoundingBox','Area','Eccentricity');
    
    o = size(im(i).intstats,1); %gets the length
    temp = [temp,im(i).intstats.Area];
    temp = temp';
    max_ = max(temp);
    
    for j = 1:o %finds the max largest area
        if max_ == im(i).intstats(j).Area;
            im(i).stats = im(i).intstats(j,:);
%             im(i).max_ind = j;
%             disp(j);
        end
    end    
    temp = [];

end


% max_ = max(temp(:));

% ------------ Preview ------------------------------%
%%

% i = 1;
% j = 1;
% 
% for i = 1:num_images
% % %     hold on;
%     figure;
% %     imshow(NP(i).nonpara);
% %     imshow(im(i).image);
% %     rectangle('Position', [stats(i).stats.BoundingBox],...
% %           'EdgeColor', 'r',...
% %           'LineWidth', 3,...
% %           'LineStyle','-');
% %     hold on;
% %     scatter(stats(i).stats.Centroid(1),stats(i).stats.Centroid(2),'+','g','LineWidth',2.5);
% 
% 
%   imshow(im(i).bw_);
% %   imshow(NP(i).nonpara);
%   
% % % % PLACES BOUNDING BOX
%   rectangle('Position', [im(i).stats(1).BoundingBox(1,:)],...
%           'EdgeColor', 'r',...
%           'LineWidth', 3,...
%           'LineStyle','-');
%       
%       
% %   imshow(im(i).lab(:,:,3));
% %   imshow(bw(i).B1);
% end

%SAMPLE OF POST PROCESSING

% figure(1);
% row =1;
% col = 5;
% subplot(row,col,1)
% imshow(im(1).image);
% title('Original Image');
% 
% subplot(row,col,2)
% imshow(im(1).gray);
% title('Grayscale of Image');
% 
% subplot(row,col,3)
% imshow(im(1).lab);
% title('L*a*b of Image');
% 
% subplot(row,col,4)
% imshow(im(1).bw_);
% title('Segmentation of Image');
% 
% subplot(row,col,5)
% imshow(im(1).bw_);
% rectangle('Position', [im(1).stats(1).BoundingBox(1,:)],...
%           'EdgeColor', 'r',...
%           'LineWidth', 3,...
%           'LineStyle','-');
% title('Blob Analysis of Image');
% ------------ Graphs ------------------------------%
%%
Hue = [];
Area = [];
Ecc = [];
% xcoords = [];
% ycoords = [];
% 
for i = 1:num_images
    Hue = [Hue;im(i).Hue];
    Area = [Area;im(i).stats.Area];
    Ecc = [Ecc;im(i).stats.Eccentricity];
%     xcoords = [xcoords, stats(i).xcoords];
%     ycoords = [ycoords, stats(i).ycoords];
end
% 
%Normalizes data points
% Hue = normalize(Hue,'range');
Area = normalize(Area,'range');

figure(1);
scatter(Hue,Ecc, leg);
% legend('apple','banana','orange');
legend('apple','orange');
title('Feature Extraction');
xlabel('Hue');
ylabel('Eccentricity');
hold on;

% banana = [Hue,Ecc];

% % cent_x = (centx');
% % cent_y = (centy');
% 
% % 
% % 
% fps=29.6;
% t_int=1./fps;
% t=0:t_int:(length(centx)-1)*t_int;
% t_ = linspace(0,(length(centx)-1)*t_int);
% 
% %Polyfit
% p_x = polyfit(t,xcoords,2);
% f_x = polyval(p_x,t_);
% 
% p_y = polyfit(t,ycoords,2);
% f_y = polyval(p_y,t_);
% 
% 
% 


% row = 2 ;
% col = 2 ;
% figure(2);
% subplot(row,col,1);
% scatter(t, -((centy)-size(im(1).image,1))); %720 (image height) corrects orientation of Y
% % scatter(t, ((centy)));
% title('y direction');
% ylabel('y (pixels)');
% xlabel('t (secs)');
% 
% subplot(row,col,2);
% scatter(t, centx);
% title('x direction');
% xlabel('t (secs)');
% ylabel('x (pixels)');
% 
% subplot(row,col,3);
% scatter(t, ycoords); %720 (image height) corrects orientation of Y
% plot(t,ycoords,'ko',t_,f_y,'r--');
% legend('data','linear fit');
% title('y direction');
% ylabel('y (m)');
% xlabel('t (secs)');
% 
% subplot(row,col,4);
% % scatter(t, xcoords);
% plot(t,xcoords,'ko',t_,f_x,'r--');
% legend('data','linear fit');
% % hold on;
% title('x direction');
% xlabel('t (secs)');
% ylabel('x (m)');

% ------------ Data Export ------------------------------%
%%

% %Capture a series of plots for increasing values
fig = figure;
for i = 1:num_images
% %     hold on;
%   imshow(im(i).bw);
%   imshow(bw(i).B1);

%     figure(1);
%     imshow(NP(i).nonpara);
    imshow(im(i).image);
    title(['Image ' num2str(i) ' of ' num2str((num_images)) ]);
%     rectangle('Position', [stats(i).stats.BoundingBox],...
%           'EdgeColor', 'r',...
%           'LineWidth', 3,...
%           'LineStyle','-');
%     hold on;
%     scatter(stats(i).stats.Centroid(1),stats(i).stats.Centroid(2),'+','g','LineWidth',2.5);
    drawnow
    frame = getframe(fig);
    im_{i} = frame2im(frame);
% 
% 
end
close;
% 
% 
% %Display the series of images in one figure.
% figure;
% for idx = 1:num_images
%     subplot(4,3,idx);
%     imshow(im_{idx});
% end
% 
% % Save the nine images into a GIF file. Because three-dimensional data is not supported for GIF files, call rgb2ind to 
% % convert the RGB data in the image to an indexed image A with a colormap map. To append multiple images to the first image, 
% % call imwrite with the name-value pair argument 'WriteMode','append'.
% 
% % filename = 'testAnimated.png'; % Specify the output file name
% % t_delay = 0.4;
% % 
% % for idx = 1:num_images
% %     [A,map] = rgb2ind(im_{idx},256);
% %     if idx == 1
% %         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',t_delay);
% %     else
% %         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',t_delay);
% %     end
% % end
