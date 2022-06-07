% load('/hdd2/RecBCD2/codedev/Analysis/EXP-22-BY4440/therun/prod/Pos11/trackedCells.mat');
radconst = (pi/180);
% ok, thats fast enought - 8ms per cell on average
tic
% figure
% hold on

for ci = 1:200%:length(mCells)
        
    for ti = 1%:mCells(ci).lifeTime
        x = mCells(ci).boundaries{ti}(:,1);
        y = mCells(ci).boundaries{ti}(:,2);
        
        refpoint = [mean(x),mean(y)];
        
        x = x - mean(x);
        y = y - mean(y);
        
        tmp = zeros(180,1);
        for theta = 1:180
            radtheta = theta * radconst;
            x1 = x*cos(radtheta) - y*sin(radtheta);
            y1 = x*sin(radtheta) + y*cos(radtheta);
            tmp(theta) = max(x1)-min(x1);
        end
        
        [~,minthetaloc] = min(tmp);
        mintheta = min(minthetaloc) * radconst;
        x1 = x*cos(mintheta) - y*sin(mintheta);
        y1 = x*sin(mintheta) + y*cos(mintheta);
%         nexttile
%         plot(x1,y1)
%         axis equal
%         drawnow
    end
end
toc