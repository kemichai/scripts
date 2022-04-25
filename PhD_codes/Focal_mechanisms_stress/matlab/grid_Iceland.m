function [centers,contents,boxes,ind]=socal(xydata, maxnumber, minnumber, minmagnitude, lim, minsize, saveflag, rootname);
%function [centers,contents,boxes,ind]=socal(xydata, maxnumber, minnumber, minmagnitude, lim, minsize, saveflag, rootname);
%
% Bins [X Y] data using the QUADTREE algorithm in the Saga toolbox
% Plotting data, 'centers' and 'b' are returned to Matlab, and the actual datapoints in each box
% are saved to disk (if selected).

% Customised for use with Icelandic data set in July 2013
% John Townend, UW-Madison, 15 July 2013

% Extract a subset 'd' of 'xydata' meeting the minimum magnitude criterion

% Handle input

if nargin < 1
  error(' Not enough input arguments.')
end

if nargin<8, rootname='test'; end %default output name
if nargin<7, saveflag=0; end %default is to not save
if nargin<6, minsize=0; end	%no limit to minimum box size
if nargin<5, lim=[min(xydata(:,1)) max(xydata(:,1)) min(xydata(:,2)) max(xydata(:,2))]; end
if nargin<4, minmagnitude=2; end
if nargin<3, minnumber=10; end
if nargin<2, maxnumber=100; end



% Decide whether to begin algorithm

currentbox=min(abs(lim(1)-lim(2)),abs(lim(3)-lim(4)));
s=xydata(:,6)>minmagnitude;

if sum(s)<=maxnumber
   error(['Fewer than ' int2str(maxnumber) ' events larger than ' int2str(minmagnitude)])
end
if currentbox<=minsize
   error('Quadrature will produce boxes smaller than minimum allowable size')
end

% Extract working data and report number of events

d=xydata(s,:);
disp([int2str(sum(s)) ' events to analyze'])

% NOTE:  FROM HERE ON, 'd' IS THE DATA MATRIX FOR EVERYTHING
% ('xydata' IS NO LONGER USED)

% Call QUADTREE

[ind,bx,by,nb,lx,ly] = quadtree(d(:,1),d(:,2),ones(length(d),1),maxnumber,lim,minsize);

% Extract points in each box and save for all non-empty boxes

b=zeros(size(nb,2),1);	%initialize a vector to be used in finding empty boxes
box=zeros(size(nb,2),2);

for jj = 1:size(nb,2)
   nn=find(ind==jj);
   subset=d(nn,[3 4 5]);  %data in box:  extract strike, dip, rake -> cols. 1-3
   b(jj)=size(subset,1);
   centers(jj,:)=0.5*[lx(jj,1)+lx(jj,2) ly(jj,1)+ly(jj,2)];
   depths(jj)=nanmean(d(nn,end));
   if (b(jj)>0)
       mags(jj,:)=[min(d(nn,6)) max(d(nn,6))];
   else
       mags(jj,:)=[NaN NaN];
   end
   % Save data if 'saveflag' is set to 1
   if saveflag==1
      
      % Set up output format 'of' to write input files for Michael's inversion 
      
      of=char('%7.2f %7.2f %7.2f\n');	
           
      % Only save data in boxes with more than 'minnumber' of events
      
      if b(jj)>=minnumber 
         filename=[rootname '_' sprintf('%0.4d',jj)];
         header=['%' int2str(minnumber) '<n<' int2str(maxnumber) ' M' int2str(minmagnitude) ...
               '+ events per box, ' int2str(minsize) ' km min. size, box ' int2str(jj) ' of ' ...
               int2str(size(nb,2)) ];         
         fid=fopen(filename,'w');
         fprintf(fid,'%s\n',header);
         
         %convert strike to dip direction

         subset(:,1)=mod(subset(:,1)+90,360);
         
         %write columns 1, 2, 3 (dip direction, dip, rake)
         
         fprintf(fid,of,subset');
         %the transpose makes it save each column of the data 
         %in a different column of the file (otherwise it works row-wise)
         
         fclose(fid);
                 
      end
   end
end

% Plot data and boxes
xb = lx(:,[1 2 2 1 1])';
yb = ly(:,[1 1 2 2 1])';
plot(d(:,1),d(:,2),'.','MarkerSize',1,'MarkerEdgeColor',[0.59 0.59 0.59])		%data

axis equal
hold on
boxes=[lx ly];
plot(xb,yb,':','color','k')		%all boxes

for i=1:length(xb)
   if b(i)>=minnumber
      plot(xb(:,i),yb(:,i),'-','color','k')		%filled boxes
   end
end
title(['Criteria:  ' int2str(minnumber) '-' int2str(maxnumber) ' (inc.) M' int2str(minmagnitude) ...
      '+ events/box, ' int2str(minsize) ' km min. size'],'fontsize',10);
xlabel([int2str(sum(b>=minnumber)) '/' int2str(size(nb,2)) ' boxes found for ' int2str(sum(s)) ' events'])
hold off

depths=depths';
contents=b;