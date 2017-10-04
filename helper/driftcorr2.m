% Drift Correction 2 - center of mass shift

% h = waitbar(0,'Analyzing...');
% % Convert cell indexing to matrices of coordinates
% trkID = objs_link(6,:); numTrajs = max(trkID); 
% trks = cell(numTrajs,1); frms = trks;
% X = NaN(length(trks),Nframes); Y = X; b = X; d = X; s = X;
% for k = 1:numTrajs
%     trks{k} = objs_link(1:2,trkID==k);
%     frms{k} = objs_link(5,trkID==k);
%     X(k,frms{k}) = objs_link(1,trkID==k);
%     Y(k,frms{k}) = objs_link(2,trkID==k);
%     b(k,frms{k}) = objs_link(3,trkID==k);
%     d(k,frms{k}) = objs_link(8,trkID==k);
%     s(k,frms{k}) = objs_link(7,trkID==k);
%     waitbar(k/numTrajs)
% end
% close(h)
% Get center of mass displacements

% dx = zeros(numTrajs,Nframes-1); dy = dx;
% dCOM = zeros(Nframes-1,2); % COM = zeros(Nframes,2); COM2 = zeros(Nframes,2);
dCOM = cell(Nframes-1,1);
for k = 2:Nframes
    ind = sum(X(:,k-1:k+T-2),2)==T & sum(X(:,k:k+T-1),2)==T;
    for l = 1:T
        COM1 = [mean(X(ind,k+l-2)),mean(Y(ind,k+l-2))];
        COM2 = [mean(X(ind,k+l-1)),mean(Y(ind,k+l-1))];
        dCOM{k-2+l} = [dCOM{k-2+l},COM2 - COM1];
    end
    
%     ind = ~isnan(X(:,k-1)) & ~isnan(X(:,k));
    
%     dx(k-1,:) = diff([XC(ind,k-1),XC(ind,k)],1,2);
%     dy(k-1,:) = diff([YC(ind,k-1),YC(ind,k)],1,2);
%     if sum(ind)>0
%         w1 = B(ind,k-1)/sum(ind);%./d(ind,k-1)/sum(ind); 
%         w2 = B(ind,k)/sum(ind);%./d(ind,k)/sum(ind);
%         W1 = sum(w1); W2 = sum(w2);
%         COM1 = [sum(w1.*X(ind,k-1))/W1,sum(w1.*Y(ind,k-1))/W1];
%         COM2 = [sum(w2.*X(ind,k))/W2,sum(w2.*Y(ind,k))/W2];
%     else
%         COM1(k-1,:) = [0,0]; COM2(k,:) = [0,0];
%     end
%     if sum(ind)>0
%         w1 = 1./s(ind,k-1).^2/sum(ind); 
%         w2 = 1./s(ind,k).^2/sum(ind);
% %         w1 = B(ind,k-1)/sum(ind);%./d(ind,k-1)/sum(ind); 
% %         w2 = B(ind,k)/sum(ind);%./d(ind,k)/sum(ind);
%         W1 = sum(w1); W2 = sum(w2);
%         if k == 2
%             % COM1(k-1,:) = [sum(w1.*X(ind,k-1))/W1,sum(w1.*Y(ind,k-1))/W1];
%             COM(k-1,:) = [sum(w1.*X(ind,k-1))/W1,sum(w1.*Y(ind,k-1))/W1];
%             COM2(k-1,:) = [mean(X(ind,k-1)),mean(Y(ind,k-1))];
%         else
%             COM(k,:) = [sum(w2.*X(ind,k))/W2,sum(w2.*Y(ind,k))/W2];
%             COM2(k,:) = [mean(X(ind,k)),mean(Y(ind,k))];
%         end
%     else
%         if k == 2
%             COM(k-1,:) = [NaN,NaN]; COM2(k-1,:) = [NaN,NaN];
%         else
%             COM(k,:) = [NaN,NaN]; COM2(k,:) = [NaN,NaN];
%         end
%     end
% %     if sum(ind)>0
%         COM1 = [sum(1/sum(ind)*XC(ind,k-1)),sum(1/sum(ind)*YC(ind,k-1))];
%         COM2 = [sum(1/sum(ind)*XC(ind,k)),sum(1/sum(ind)*YC(ind,k))];
%     else
%         COM1 = [0,0]; COM2 = [0,0];
%     end
%     COM1 = [mean(XC(ind,k-1)),mean(YC(ind,k-1))];
%     COM2 = [mean(XC(ind,k)),mean(YC(ind,k))];
    dCOM{k-1,:} = [COM2(1)-COM1(1),COM2(2)-COM1(2)];
    % waitbar(k/Nframes)
end
% cumDiff = [[0;0],cumsum(dCOM,1)'];
% close(h)

% Plot COM
figure,plot(COM1(:,1),COM1(:,2),'r.'),set(gca,'Xlim',[0 512],'YLim',[0 512])
hold on, plot(COM2(:,1),COM2(:,2),'k.'),plot(C{57}(1,:),C{57}(2,:),'b.'),hold off

% % Apply correction
% Xcorr = X - cumDiff(1,:); Ycorr = Y - cumDiff(2,:);
% 
% % Convert into individual trajectory cell form
% Ccorr = cell(size(trks));
% for k = 1:length(Ccorr)
%     xs = Xcorr(k,:); ys = Ycorr(k,:);
%     xs(isnan(xs)) = []; ys(isnan(ys)) = [];
%     if ~isempty(xs)
%         frames = frms{k};
%         Ccorr{k} = [xs',ys',frames'];
%     end
% end
% Ccorrmat = cell2mat(Ccorr);
% trksmat = cell2mat(cellfun(@(x) x',trks,'UniformOutput',0));
% frmsmat = cell2mat(cellfun(@(x) x',frms,'UniformOutput',0));
% orig = [trksmat,frmsmat];
% 
% % Display corrected points
% figure, % h = gca; h.XLim = [1 512]; h.YLim = [1 512]; 
% hold on
% for k = 1:Nframes
%     centers = Ccorrmat(Ccorrmat(:,3) == k,1:2);
%     plot(centers(:,1),512-centers(:,2),'b.')
%     origcenters = orig(orig(:,3) == k,1:2);
%     plot(origcenters(:,1),512-origcenters(:,2),'r.')
%     pause(0.025);
% end
