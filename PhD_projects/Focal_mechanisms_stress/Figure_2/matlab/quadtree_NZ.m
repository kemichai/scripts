clear
data=load('input/matlab_quadtree_input.dat');

plottype='patch';

% Grid data (do 'help grid_Iceland' for explanation of parameters);
% Minimum number of events per bin
min_number=15;
% Maximum number of events per bin
max_number=50;
% Limits for region to be divided
limits=[-3.55*10^6 -3.3*10^6 4.05*10^6 4.3*10^6 ];
% Minimum size for the bins
min_size=0.005*10^6;% NOT SURE IN WHAT UNITS THIS SIZE IS...


[centres, contents, boxes, ind]=grid_Iceland([data 9*ones(length(data),4)],max_number,min_number,0,limits,min_size);


% Compute and plot mean (or median) depths (Helps with readjusting the size
% and distributions of the bins
hold on
switch lower(plottype)
    case 'scatter'
        for i=1:max(ind)
            % Check if box contains enough eqs
            if sum(ind==i)>=min_number
                % Compute average epicentre of eqs in box
                avX(i)=mean(data(ind==i,1));
                avY(i)=mean(data(ind==i,2));
                % Compute average depth of eqs in box (using mean or median)
                % avdepth(i)=mean(data(ind==i,3));
                avdepth(i)=prctile(data(ind==i,3),90);
                scatter(avX,avY,20,avdepth,'filled')
            end
        end
    case 'patch'
        for i=1:max(ind)
            % create lists even if there are not enough events in the box
            seismic_moment(i)=0;
            median_dep(i)=0;
            % create boxes too
            X=boxes(i,[1,2]);
            Y=boxes(i,[3,4]);
            X1_(i)=X(1);
            X2_(i)=X(2);
            Y1_(i)=Y(1);
            Y2_(i)=Y(2);
            % Check if box contains enough EQs
            if sum(ind==i)>=min_number
                % avdep=prctile(data(ind==i,3),95);
                % seismic moment rate
                seis_mom=sum(data(ind==i,4))
                X=boxes(i,[1 2]);
                Y=boxes(i,[3 4]);
                X1_(i)=X(1);
                X2_(i)=X(2);
                Y1_(i)=Y(1);
                Y2_(i)=Y(2);
                patch([X(1) X(2) X(2) X(1) X(1)],[Y(1) Y(1) Y(2) Y(2) Y(1)],seis_mom)
                %avdepth2(i)=prctile(data(ind==i,3),95);
                % seismic moment rate...
                % Create a list of the seismic moment sums
                % and a list of the median depths for each box
                seismic_moment(i)=sum(data(ind==i,4))
                median_dep(i)=median(data(ind==i,3));
                
            end
        end
    otherwise
        disp('Unknown plot type')
end

axis tight
colorbar

% Plot grid boxes along with coast line and faults to re-adjust the grids
hold on
fault=load('faults.xy');
map = load('nz_coast.xy');
% scatter(map(:,1),map(:,2),10,'k','o')
% hold on
plot(map(:,1),map(:,2),'k-','MarkerSize',10)
hold on
samba = load('SAMBA.dat')
scatter(samba(:,1),samba(:,2),50,'r','v','filled','MarkerEdgeColor','k')
geo = load('GEONET.dat')
scatter(geo(:,1),geo(:,2),50,[0.9290, 0.6940, 0.1250],'v','filled','MarkerEdgeColor','k')
hold off


av=seismic_moment';
dep=median_dep';
c = [data,ind]


% focal mechanism info and index of box that each event is in.
csvwrite('output/cat_box_num.csv',c)
% Contains details on the eqz within each bin (lat, lon, median depth and
% number of events within the bin).
csvwrite('output/cluster_details.csv',[centres,contents,dep])
% ------------------------------------------------------------------------
% The two files below are created in order to plot the seismic moment
% release within each box.
% Sum of seismic moment for each bin
csvwrite('output/final_values.csv',av) 
% Coordinates of the four edges of each bin
csvwrite('output/final_boxes.csv',boxes)


