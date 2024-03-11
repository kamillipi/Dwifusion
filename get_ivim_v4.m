function [results] = get_ivim_v4(bvals,data,varargin)
%out = get_ivim(bvals,data,method,bsplit)
%provide this function with:
%bvals - b values vector
%data - 2D (x, b values) 3D (x,y,b values) or 4D (x,y,z,b values) with IVIM data
%method - "1step" Fits all parameters at once
%         "segmented" fit for D and S0 and then for Dstar and f
%         "grids" 2 step grid search for D and S0 and Dstar

p = inputParser; %branch test

allowedmethods = {'1step','segmented','grid'};
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p,'bvals');
addRequired(p,'data');
addParameter(p,'method','1step',@(x) any(validatestring(x,allowedmethods)));
addParameter(p,'bsplit',250,validScalarPosNum);
addParameter(p,'mask',0);
parse(p,bvals,data,varargin{:});
dimensions=size(data);
ndimensions=numel(dimensions);
number_of_points1=int16(360);
number_of_points2=int16(300);
Dstar_min=0.002;
Dstar_max=0.085;
D_min=5e-4;
D_max=0.002;
f_min=0.001;
f_max=0.5;
switch ndimensions
    case 2
        table_cols=dimensions(2);
        table_rows=dimensions(1);
        D=zeros(dimensions(1),1);
    case 3
        table_cols=dimensions(3);
        table_rows=dimensions(1)*dimensions(2);
        D=zeros(dimensions(1),dimensions(2));
    case 4
        table_cols=dimensions(4);
        table_rows=dimensions(1)*dimensions(2)*dimensions(3);
        D=zeros(dimensions(1),dimensions(2),dimensions(3));
    otherwise
        disp("Please provide proper data: 2,3,4 dimensional where last dimension corresponds with bvals" + newline)
        results=-1;
        return
end

sorted_data=reshape(data,table_rows,table_cols);

if isequal(p.Results.mask,0)
    mask=max(data,[],ndimensions)>0.01*max(data,[],"all");
    sorted_mask=reshape(mask,table_rows,1);
    [values_location,~]=find(sorted_mask>0);
    disp("Using thresholded mask, calculating " + sum(mask,'all')/numel(mask)*100 + "% of data, which is " + sum(mask,'all') + " voxels");
elseif isequal(numel(p.Results.mask),1)
    mask=max(data,[],ndimensions)>p.Results.mask;
    sorted_mask=reshape(mask,table_rows,1);
    [values_location,~]=find(sorted_mask>0);
    disp("Using thresholded mask, calculating " + sum(mask,'all')/numel(mask)*100 + "% of data, which is " + sum(mask,'all') + " voxels");

else
    sorted_mask=reshape(p.Results.mask,numel(p.Results.mask),1);
    sorted_mask=repmat(sorted_mask,table_rows/numel(p.Results.mask),1);
    [values_location,~]=find(sorted_mask>0);
    disp("Using provided mask, calculating " + sum(p.Results.mask,'all')/numel(p.Results.mask)*100 + "% of data, which is " + sum(p.Results.mask,'all') + " voxels");

end
disp("Started at " + datestr(datetime));
to_calculation=double(sorted_data(values_location,:)');
calculated_values=zeros(numel(values_location),4);
if ~iscolumn(bvals)
    bvals=bvals';
end
top_signal=double(max(data,[],'all'));
switch p.Results.method
    %     case "lsqf"
    %         sorted_data=sorted_data';
    %         xes=repmat(bvals,1,length(to_calculation));
    %         fun=@(x,xdata)x(1).*x(2).*exp(-xdata.*x(3))+(1-x(2)).*x(1).*exp(-xdata.*x(4));
    %         x0=[1000,0.05,0.001,0.01];
    %         res=lsqcurvefit(fun,x0,xes,to_calculation);

    case "1step"
        %         disp("Performing all parameters fitting");
        %         Fit all parameters
        fo3 = fitoptions('Method','NonlinearLeastSquares',... % D Dstar S0 fblood
            'StartPoint',[1e-3 1e-4 top_signal 0.1 ],...
            'Lower', [D_min Dstar_min 0 f_min], ...
            'Upper', [D_max Dstar_max top_signal f_max]);
        %                'TolFun', 1e-18, ...
        %                'MaxIter', 10000, ...
        %                'MaxFunEvals', 10000);
        ft3 = fittype('S0*f*exp(-x*Dstar)+(1-f)*S0*exp(-x*D)','options',fo3);

    progbar= progressBar(size(to_calculation,2),'pname','Calculating 1step');

        parfor i=1:numel(values_location)
            [fit3,~,~]=fit(bvals,squeeze(to_calculation(:,i)),ft3);
            calculated_values(i,:)=[fit3.S0 fit3.f fit3.Dstar fit3.D];
            progbar.progress
        end
    case "segmented"
        progbar= progressBar(size(to_calculation,2),'pname','Calculating segmented');
        %         disp("Performing segmented fitting");
        ivimx=bvals(bvals<p.Results.bsplit);
        nonivimx=bvals(bvals>p.Results.bsplit);

        if size(ivimx,2)>1
            ivimx=ivimx(:);
        end
        if size(nonivimx,2)>1
            nonivimx=nonivimx(:);
        end
        
        fo1 = fitoptions('Method','NonlinearLeastSquares',...
            'StartPoint',[1e-3 0.6*top_signal],...
            'Lower', [D_min 0], ...
            'Upper', [D_max top_signal]);

        ft1 = fittype('S0*exp(-x*D)','options',fo1);

        %fo2 = fitoptions('Method','NonlinearLeastSquares',...
            % 'StartPoint',[0.04 0.01],...
            % 'Upper', [0.33 0.1],...
            % 'Lower', [0 0]);
        fo2 = fitoptions('Method','NonlinearLeastSquares',...
            'StartPoint',[1e-2 0.2*top_signal],...
            'Lower', [Dstar_min 0], ...
            'Upper', [Dstar_max 0.5*top_signal]);
        %ft2 = fittype(@(f,Dstar,S0,x)(f/(1-f)*S0*exp(-x*Dstar)),'problem','S0','options',fo2);
        ft2 = fittype('S0*exp(-x*Dstar)','options',fo2);
        parfor i = 1:numel(values_location)
            progbar.progress
            nonivimy=to_calculation(bvals>p.Results.bsplit,i);
            tissueandivimy=to_calculation(bvals<p.Results.bsplit,i);
            [fit1, ~, ~]=fit(nonivimx,nonivimy,ft1);
            ivimy=tissueandivimy-fit1.S0*exp(-fit1.D*ivimx);
            %[fit2, ~, ~]=fit(ivimx,ivimy,ft2,'problem',fit1.S0);
            [fit2, ~, ~]=fit(ivimx,ivimy,ft2);
            calculated_values(i,:)=[(fit1.S0+fit2.S0) (fit2.S0/(fit1.S0+fit2.S0)) fit2.Dstar fit1.D];
            %calculated_values(i,:)=[(fit1.S0/(1-fit2.f)) fit2.f fit2.Dstar fit1.D];
        end
    case "grid"
        %         disp("Performing bayesian fitting");
        bsplit=p.Results.bsplit;


        nonivimx=bvals(bvals>bsplit);
        ivimx=bvals(bvals<bsplit);
        n1 = numel(nonivimx);
        n2 = numel(ivimx);

        deviation=std(to_calculation(find(bvals==0),:),[],"all");
        if deviation == 0
        deviation=0.01*max(data,[],"all");
        end
  
        progbar= progressBar(size(to_calculation,2),'pname','Calculating grid search');
        parfor i = 1:numel(values_location)
            progbar.progress
            nonivimy=to_calculation(bvals>bsplit,i);
            % ivimy=to_calculation(bvals<bsplit,i);
            
            %S0pred=1/exp(-min(nonivimx)*0.001)*max(nonivimy);

            %w1 = linspace(0.5*S0pred,2*S0pred,round(number_of_points*1.2)); %S0
            w1 = linspace(0.2*top_signal,top_signal,number_of_points1); %S0
            %w2 = linspace(0.0001,0.01,number_of_points); %D
            %w1 = linspace(0.9/exp(-min(nonivimx)*5e-4)*max(nonivimy),...
            % 1.1/exp(-min(nonivimx)*2e-3)*max(nonivimy),number_of_points); %S0
            w2 = linspace(D_min,D_max,number_of_points1); %D
            [vw1,vw2] = meshgrid(w1,w2);

            N = length(vw1(:));
            Y = repmat(nonivimy,1,N);
            S = repmat(vw1(:)',n1,1).*exp(-nonivimx*(vw2(:)'));


            mu = sum((Y-S).^2,1)'/2/deviation^2;
            li = exp(-mu);
            li = li/sum(li(:));
            li = reshape(li,size(vw1));
            ind = find(li==max(li(:)));
            ind = round(median(ind));
            [xindex,yindex]=ind2sub(size(vw1),ind);

            S0=w1(yindex);
            Dp=w2(xindex);


            onlyivimy=to_calculation(bvals<bsplit,i)-S0*exp(-Dp*ivimx);
            % w1 = linspace(0.02*max(to_calculation(bvals<bsplit,i)),0.08*max(to_calculation(bvals<bsplit,i)),100); %S0 %S0 ivim part
           % w1 = linspace(0,1.5*(abs(max(onlyivimy))),number_of_points); %S0 %S0 ivim part
            %w2 = linspace(5e-3,5e-2,number_of_points); %Dstar
            %w1 = linspace(1,0.5*max(to_calculation(bvals<bsplit,i)),number_of_points); %S0 %S0 ivim part
            w1 = linspace(0.001*top_signal,0.4*top_signal,number_of_points2); %S0 %S0 ivim part
            w2 = linspace(Dstar_min,Dstar_max,number_of_points2); %Dstar

            [vw1,vw2] = meshgrid(w1,w2);

            N = length(vw1(:));
            Y = repmat(onlyivimy,1,N);
            S = repmat(vw1(:)',n2,1).*exp(-ivimx*vw2(:)');
            mu = sum((Y-S).^2,1)'/2/deviation^2;
            li = exp(-mu);
            li = li/sum(li(:));
            li = reshape(li,size(vw1));
            ind = find(li==max(li(:)));
            ind = round(median(ind));
            [xindex,yindex]=ind2sub(size(vw1),ind);

            f=w1(yindex)/(S0+w1(yindex));
            Dstar=w2(xindex);
            calculated_values(i,:)=[(S0+w1(yindex)) f Dstar Dp];
        end
end
imags=zeros(size(sorted_data,1),4);
imags(values_location,:)=calculated_values(1:numel(values_location),:);

S0=reshape(imags(:,1),size(D));
f=reshape(imags(:,2),size(D));
Dstar=reshape(imags(:,3),size(D));
Dstar(f<0.01)=0;
D=reshape(imags(:,4),size(D));
results=cat(ndimensions,S0,f,Dstar,D);
disp("Ended at " + datestr(datetime));
end
