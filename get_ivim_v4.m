function [results] = get_ivim_v4(bvals,data,varargin)
%out = get_ivim(bvals,data,method,bsplit)
%provide this function with:
%bvals - b values vector
%data - 2D (x, b values) 3D (x,y,b values) or 4D (x,y,z,b values) with
%method - "1step" Fits all parameters at once
%         "segmented" fit for D and S0 and then for Dstar and f
%         "bayes" 2 step grid search for D and S0 and Dstar

p = inputParser; %branch test

allowedmethods = {'1step','segmented','bayesian'};
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p,'bvals');
addRequired(p,'data');
addParameter(p,'method','1step',@(x) any(validatestring(x,allowedmethods)));
addParameter(p,'bsplit',250,validScalarPosNum);
addParameter(p,'mask',0);
parse(p,bvals,data,varargin{:});

dimensions=size(data);
ndimensions=numel(dimensions);

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

switch p.Results.mask
    case 0
        [values_location,~]=find(sum(sorted_data,2)>1000);
         disp("Calculating data with sum >1000. " + numel(values_location)/length(sorted_data)*100 + "% of proivded data, which is " + numel(values_location) + " voxels");
    otherwise
        sorted_mask=reshape(p.Results.mask,table_rows,1);
        [values_location,~]=find(sorted_mask>0);
        disp("Using mask, calculating " + numel(values_location)/length(sorted_data)*100 + "% of data, which is " + numel(values_location) + " voxels");
end

to_calculation=double(sorted_data(values_location,:))';
calculated_values=zeros(numel(values_location),4);


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
            'StartPoint',[0.0008 0.009 max(data,[],'all') 0.05 ],...
            'Lower', [0 0 0 0.01], ...
            'Upper', [0.02 0.2 Inf 0.33]);
        %                'TolFun', 1e-18, ...
        %                'MaxIter', 10000, ...
        %                'MaxFunEvals', 10000);
        ft3 = fittype('S0*f*exp(-x*Dstar)+(1-f)*S0*exp(-x*D)','options',fo3);
        
        parfor i=1:numel(values_location)
            [fit3,~,~]=fit(bvals,squeeze(to_calculation(:,i)),ft3);
            calculated_values(i,:)=[fit3.S0 fit3.f fit3.Dstar fit3.D];
        end
    case "segmented"
        %         disp("Performing segmented fitting");
        ivimx=bvals(bvals<p.Results.bsplit)';
        nonivimx=bvals(bvals>p.Results.bsplit)';
        fo1 = fitoptions('Method','NonlinearLeastSquares',...
            'StartPoint',[0.0008 max(data,[],'all')],...
            'Lower', [0.0001 0.5*max(data,[],'all')], ...
            'Upper', [0.01 max(data,[],'all')]);
        
        ft1 = fittype('S0*exp(-x*D)','options',fo1);
        
        fo2 = fitoptions('Method','NonlinearLeastSquares',...
            'StartPoint',[0.04 0.01],...
            'Upper', [0.33 0.1],...
            'Lower', [0 0]);
        ft2 = fittype(@(f,Dstar,S0,x)(f/(1-f)*S0*exp(-x*Dstar)),'problem','S0','options',fo2);
        nonivimy=to_calculation(bvals>p.Results.bsplit,:);
        tissueandivimy=to_calculation(bvals<p.Results.bsplit,:);
        parfor i = 1:length(to_calculation)
            [fit1, ~, ~]=fit(nonivimx,nonivimy(:,i),ft1);
            ivimy=tissueandivimy(:,i)-fit1.S0*exp(-fit1.D*ivimx);
            [fit2, ~, ~]=fit(ivimx,ivimy,ft2,'problem',fit1.S0);
            calculated_values(i,:)=[(fit1.S0/(1-fit2.f)) fit2.f fit2.Dstar fit1.D];
        end
    case "bayesian"
        %         disp("Performing bayesian fitting");
        bsplit=p.Results.bsplit;
        nonivimx=bvals(bvals>bsplit)';
        ivimx=bvals(bvals<bsplit)';
        n1 = numel(nonivimx);
        n2 = numel(ivimx);
        
        deviation=std(to_calculation(bvals==0,:),[],'all');
        
        if deviation==0
            deviation=1000;
        end
        nonivimy=to_calculation(bvals>p.Results.bsplit,:);
        parfor i = 1:numel(values_location)
%           ivimy=to_calculation(bvals<bsplit,i);
            S0pred=1/exp(-min(nonivimx)*0.001)*max(nonivimy(:,i));
            
            w1 = linspace(0.9*S0pred,1.1*S0pred,100); %S0
            w2 = linspace(0.0001,0.01,200); %D
            [vw1,vw2] = meshgrid(w1,w2);
            
            N = length(vw1(:));
            Y = repmat(nonivimy(:,i),1,N);
            S = repmat(vw1(:)',n1,1).*exp(-nonivimx*vw2(:)');
            
            
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
            w1 = linspace(1,0.25*max(to_calculation(bvals<bsplit,i)),200); %S0 %S0 ivim part
            w2 = linspace(0.001,0.02,400); %Dstar
            
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
imags=zeros(length(sorted_data),4);
imags(values_location,:)=calculated_values(1:numel(values_location),:);
S0=reshape(imags(:,1),size(D));
f=reshape(imags(:,2),size(D));
Dstar=reshape(imags(:,3),size(D));
D=reshape(imags(:,4),size(D));
results=cat(4,S0,f,Dstar,D);
end