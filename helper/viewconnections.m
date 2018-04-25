function img = viewconnections(img,spotEvents,objs_link)

% Read video into 3D array for quick access
if ischar(img)
    filename = img;
    iminf = imfinfo(img);
    N = length(iminf);
    m = iminf(1).Width; n = iminf(1).Height;
    img = zeros(m,n,N,'uint16');
    TifLink = Tiff(filename,'r');
    for k = 1:N
        TifLink.setDirectory(k);
        img(:,:,k) = TifLink.read();
    end
    TifLink.close();
else
    N = size(img,3);
end

% Gather coordinates of linked trajectories
L = length(spotEvents);
trajx = nan(L,N);
trajy = nan(L,N);
for k = 1:L
    TOI = spotEvents(k).trajectory;
    t = length(TOI);
    inds = TOI(~isnan(TOI));
    if t >= 1
        trajframe = objs_link(5,inds);
        trajx(k,trajframe) = objs_link(1,inds);
        trajy(k,trajframe) = objs_link(2,inds);
    end
end

% Display points
figure
hAxesCM = axes;
axis(hAxesCM,'off')
hAxesP = axes;
hAxesP.YLim = [1 512]; hAxesP.XLim = [1 512];
for k = 1:N
%     imshow(imadjust(img(:,:,k))), hold on
%     c = linspace(1,65535,L);
%     scatter(trajx(:,k),trajy(:,k),20,c), hold off
%     pause(1)
% colormap(hAxesCM,gray);
% im = imshow(imadjust(img(:,:,k)));
% imgmin = min(min(img(:,:,k)));
% imgmax = max(max(img(:,:,k)));
% caxis(hAxesCM,[imgmin imgmax]);
colormap(hAxesP,hsv);
c = uint16(linspace(1,65535,L));
scatter(trajx(:,k),trajy(:,k),20,c); hold on
pause(1)
clear im
end

% %create data
% % [X,Y,Z] = peaks(30);
% % s = scatter(trajx(:,k),trajy(:,k),20);
% % Zprime = del2(Z);
% 
% % pcolormin = min(Zprime(:));
% % pcolormax = max(Zprime(:));
% %create figure and store handle
% hF = figure;
% %create axes for the countourm axes
% hAxesCM = axes;
% %set visibility for axes to 'off' so it appears transparent
% axis(hAxesCM,'off')
% %set colormap for contourm axes
% colormap(hAxesCM,gray);
% %plot contourm
% im = imshow(imadjust(img(:,:,k)));
% % contourmPlot = contourm(X,Y,Z,20);
% %create color bar and set range for color
% % cbCM = colorbar(hAxesCM,'Location','east');
% imgmin = min(img(:,:,k));
% imgmax = max(img(:,:,k));
% caxis(hAxesCM,[imgmin imgmax]);
% %create axes for pcolor and store handle
% hAxesP = axes;
% %set colormap for pcolor axes
% colormap(hAxesP,hsv);
% %plot pcolor for gradient
% c = uint16(linspace(1,65535,L));
% scatterPlot = scatter(trajx(:,k),trajy(:,k),20,c);
% %pcolorPlot = pcolor(X,Y,s);
% %set(pcolorPlot,'FaceColor','interp','EdgeColor','interp');
% %create color bar and set range for color
% %cbP = colorbar(hAxesP,'Location','west');
% caxis(hAxesP,[1 65535]);
% 
% %link the two overlaying axes so they match at all times to remain accurate
% linkaxes([hAxesP,hAxesCM]);

% figure,imshow(img_1)
% hold on, plot(pts(:,2),pts(:,1),'co'),hold off
% for i = 1:length(pts)
%     if currentTraj(i) ~= 0
%         text((double(pts(i,2))+4),double(pts(i,1)),sprintf('%d',currentTraj(i)),'Color','c');
%     end
% end
    