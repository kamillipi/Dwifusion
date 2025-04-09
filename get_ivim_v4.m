function [results] = get_ivim_v4(bvals,signal_data,varargin)
%out = get_ivim(bvals,data,method,bsplit)
%provide this function with:
%bvals - b values vector
%data - 2D (x, b values) 3D (x,y,b values) or 4D (x,y,z,b values) with IVIM data
%method - "1step" Fits all parameters at once
%         "segmented" fit for D and S0 and then for Dstar and f
%         "grid" 2 step grid search for D and S0 and Dstar

p = inputParser; %branch test

allowedmethods = {'1step','seg','grid','v_grid','v_seg'};
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p,'bvals');
addRequired(p,'signal_data');
addParameter(p,'method','1step',@(x) any(validatestring(x,allowedmethods)));
addParameter(p,'bsplit',250,validScalarPosNum);
addParameter(p,'mask',0);
addParameter(p,'mask_action',"calculate");
parse(p,bvals,signal_data,varargin{:});
dimensions=size(signal_data);
ndimensions=numel(dimensions);
number_of_points1=int16(300);
number_of_points2=int16(360);
Dstar_min=0;
Dstar_max=0.4;
D_min=0;
D_max=0.0025;
f_min=0.001;
f_max=1;
v_min=0;
v_max=3;

if p.Results.mask_action=="average"
    avmask=double(p.Results.mask);
    avmask(p.Results.mask==0)=NaN;
    switch ndimensions
        case 2
            signal_data=squeeze(mean(signal_data.*single(avmask),[1],"omitnan"));
        case 3
            signal_data=squeeze(mean(signal_data.*single(avmask),[1 2],"omitnan"));
        case 4
            signal_data=squeeze(mean(signal_data.*single(avmask),[1 2 3],"omitnan"));
        otherwise
            disp("Please provide proper data: 2,3,4 dimensional where last dimension corresponds with bvals" + newline)
            results=-1;
            return
    end
    if ~isrow(signal_data)
        signal_data=signal_data';
    end
    dimensions=size(signal_data);
    ndimensions=numel(dimensions);
end

switch ndimensions
    case 2
        
        table_cols=dimensions(2);
        table_rows=dimensions(1);
        D=zeros(1,1);
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


sorted_data=reshape(signal_data,table_rows,table_cols);
if p.Results.mask_action=="calculate"
    if isequal(p.Results.mask,0)
        mask=max(signal_data,[],ndimensions)>0.001*max(signal_data,[],"all");
        sorted_mask=reshape(mask,table_rows,1);
        [values_location,~]=find(sorted_mask>0);
        disp("Using thresholded mask, calculating " + sum(mask,'all')/numel(mask)*100 + "% of data, which is " + sum(mask,'all') + " voxels");
    elseif isequal(numel(p.Results.mask),1) && p.Results.mask>1
        mask=max(signal_data,[],ndimensions)>p.Results.mask;
        sorted_mask=reshape(mask,table_rows,1);
        [values_location,~]=find(sorted_mask>0);
        disp("Using thresholded mask, calculating " + sum(mask,'all')/numel(mask)*100 + "% of data, which is " + sum(mask,'all') + " voxels");
    else
        sorted_mask=reshape(p.Results.mask,numel(p.Results.mask),1);
        if(table_rows==1||table_cols==1)
            sorted_data=sorted_data';
        else
            sorted_mask=repmat(sorted_mask,table_rows/numel(p.Results.mask),1);
        end
        [values_location,~]=find(sorted_mask>0);
        disp("Using provided mask, calculating " + sum(p.Results.mask,'all')/numel(p.Results.mask)*100 + "% of data, which is " + sum(p.Results.mask,'all') + " voxels");
    end
else
    values_location=1:dimensions();
end
disp("Started at " + datestr(datetime));
to_calculation=double(sorted_data(values_location,:));
if iscolumn(to_calculation)
else
    to_calculation=to_calculation';
end
calculated_values=zeros(numel(values_location),4);
if ~iscolumn(bvals)
    bvals=bvals';
end
top_signal=double(max(signal_data,[],'all'));
switch p.Results.method

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
    case "seg"
        progbar= progressBar(size(to_calculation,2),'pname','Calculating segmented');
        %         disp("Performing segmented fitting");
        ivimb=bvals(bvals<p.Results.bsplit);
        nonivimb=bvals(bvals>p.Results.bsplit);

        if size(ivimb,2)>1
            ivimb=ivimb(:);
        end
        if size(nonivimb,2)>1
            nonivimb=nonivimb(:);
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
            [fit1, ~, ~]=fit(nonivimb,nonivimy,ft1);
            ivimy=tissueandivimy-fit1.S0*exp(-fit1.D*ivimb);
            [fit2, ~, ~]=fit(ivimb,ivimy,ft2);
            calculated_values(i,:)=[(fit1.S0+fit2.S0) (fit2.S0/(fit1.S0+fit2.S0)) fit2.Dstar fit1.D];
            %calculated_values(i,:)=[(fit1.S0/(1-fit2.f)) fit2.f fit2.Dstar fit1.D];
        end
    case "v_seg"
        Delta=0.043192;
        delta=0.023848;
        progbar= progressBar(size(to_calculation,2),'pname','Calculating segmented');
        %         disp("Performing segmented fitting");
        ivimb=bvals(bvals<p.Results.bsplit);
        nonivimb=bvals(bvals>p.Results.bsplit);

        if size(ivimb,2)>1
            ivimb=ivimb(:);
        end
        if size(nonivimb,2)>1
            nonivimb=nonivimb(:);
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
            'StartPoint',[0.2*top_signal 1e-2 ],...
            'Lower', [0 v_min], ...
            'Upper', [0.5*top_signal v_max]);
        %ft2 = fittype(@(f,Dstar,S0,x)(f/(1-f)*S0*exp(-x*Dstar)),'problem','S0','options',fo2);
        ft2 = fittype('S0*exp(-sqrt(x/(Delta - delta/3))*V)', ...
              'independent', 'x', ...
              'dependent', 'y', ...
              'problem', {'Delta', 'delta'}, ...
              'coefficients', {'S0', 'V'}, ...
              'options', fo2);
        for i = 1:numel(values_location)
            progbar.progress
            nonivimy=to_calculation(bvals>p.Results.bsplit,i);
            tissueandivimy=to_calculation(bvals<p.Results.bsplit,i);
            [fit1, ~, ~]=fit(nonivimb,nonivimy,ft1);
            ivimy=tissueandivimy-fit1.S0*exp(-fit1.D*ivimb);
            [fit2, ~, ~]=fit(ivimb,ivimy,ft2,'problem', {Delta, delta});
            calculated_values(i,:)=[(fit1.S0+fit2.S0) (fit2.S0/(fit1.S0+fit2.S0)) fit2.V fit1.D];
            %calculated_values(i,:)=[(fit1.S0/(1-fit2.f)) fit2.f fit2.Dstar fit1.D];
        end
    case "grid"
        %         disp("Performing bayesian fitting");
        bsplit=p.Results.bsplit;
        nonivimb=bvals(bvals>bsplit);
        ivimb=bvals(bvals<bsplit);
        n1 = numel(nonivimb);
        n2 = numel(ivimb);
        first_ivim_bval=min(nonivimb);
        %max_attenuation=exp(-first_ivim_bval*1.4e-3); % diffusion of healthy grey matter
        deviation=std(to_calculation(find(bvals==0),:),[],"all");
        if deviation == 0
            deviation=0.01*max(signal_data,[],"all");
        end


        w2Dstar = linspace(Dstar_min,Dstar_max,number_of_points2); %Dstar
        w2D = linspace(D_min,D_max,number_of_points1); %D
        progbar= progressBar(size(to_calculation,2),'pname','Calculating grid search');
        parfor i = 1:numel(values_location)
            progbar.progress
            nonivimy=to_calculation(bvals>bsplit,i);
            ivimy=to_calculation(bvals<bsplit,i);
            max_attenuation=max(nonivimy)/max(ivimy);
            top_signal=max(nonivimy)/max_attenuation;

            w1 = linspace(0.75*top_signal,top_signal,number_of_points1); %S0

            [vw1,vw2] = meshgrid(w1,w2D);

            N = length(vw1(:));
            Y = repmat(nonivimy,1,N);
            S = repmat(vw1(:)',n1,1).*exp(-nonivimb*(vw2(:)'));


            mu = sum((Y-S).^2,1)'/2/deviation^2;
            li = exp(-mu);
            %li = exp(-mu)/sum(li(:));
            li = reshape(li,size(vw1));
            ind = find(li==max(li(:)));
            ind = round(median(ind));
            [xindex,yindex]=ind2sub(size(vw1),ind);

            S0=w1(yindex);
            Dp=w2D(xindex);


            onlyivimy=ivimy-S0*exp(-Dp*ivimb);
            max_ivim_signal=max(onlyivimy);
            if max_ivim_signal>0.02*S0
          
                w1 = linspace(0.5*max_ivim_signal,2*max_ivim_signal,number_of_points2); %S0 %S0 ivim part

                [vw1,vw2] = meshgrid(w1,w2Dstar);

                N = length(vw1(:));
                Y = repmat(onlyivimy,1,N);
                S = repmat(vw1(:)',n2,1).*exp(-ivimb*vw2(:)');
                mu = sum((Y-S).^2,1)'/2/deviation^2;
                li = exp(-mu);
                %li = exp(-mu)/sum(li(:));
                li = reshape(li,size(vw1));
                ind = find(li==max(li(:)));
                ind = round(median(ind));
                [xindex,yindex]=ind2sub(size(vw1),ind);

                f=w1(yindex)/(S0+w1(yindex));
                Dstar=w2Dstar(xindex);
                calculated_values(i,:)=[(S0+w1(yindex)) f Dstar Dp];
            else
                Dstar=0;
                f=0;
                calculated_values(i,:)=[S0 f Dstar Dp];
            end

        end
    case "v_grid"
        %         disp("Performing bayesian fitting");
        bsplit=p.Results.bsplit;
        nonivimb=bvals(bvals>bsplit);
        ivimb=bvals(bvals<bsplit);
        n1 = numel(nonivimb);
        n2 = numel(ivimb);
        first_ivim_bval=min(nonivimb);
        %max_attenuation=exp(-first_ivim_bval*1.4e-3); % diffusion of healthy grey matter
        deviation=std(to_calculation(find(bvals==0),:),[],"all");
        if deviation == 0
            deviation=0.01*max(signal_data,[],"all");
        end


        w2V = linspace(v_min,v_max,number_of_points2); %V
        w2D = linspace(D_min,D_max,number_of_points1); %D
        progbar= progressBar(size(to_calculation,2),'pname','Calculating grid search');
        parfor i = 1:numel(values_location)
            progbar.progress
            nonivimy=to_calculation(bvals>bsplit,i);
            ivimy=to_calculation(bvals<bsplit,i);
            max_attenuation=max(nonivimy)/max(ivimy);
            top_signal=max(nonivimy)/max_attenuation;

            w1 = linspace(0.75*top_signal,top_signal,number_of_points1); %S0

            [vw1,vw2] = meshgrid(w1,w2D);

            N = length(vw1(:));
            Y = repmat(nonivimy,1,N);
            S = repmat(vw1(:)',n1,1).*exp(-nonivimb*(vw2(:)'));


            mu = sum((Y-S).^2,1)'/2/deviation^2;
            li = exp(-mu);
            %li = exp(-mu)/sum(li(:));
            li = reshape(li,size(vw1));
            ind = find(li==max(li(:)));
            ind = round(median(ind));
            [xindex,yindex]=ind2sub(size(vw1),ind);

            S0=w1(yindex);
            Dp=w2D(xindex);


            onlyivimy=ivimy-S0*exp(-Dp*ivimb);
            max_ivim_signal=max(onlyivimy);
            if max_ivim_signal>0.02*S0
          
                w1 = linspace(0.5*max_ivim_signal,2*max_ivim_signal,number_of_points2); %S0 %S0 ivim part

                [vw1,vw2] = meshgrid(w1,w2V);
                Delta=0.043192;
                delta=0.023848;
                N = length(vw1(:));
                Y = repmat(onlyivimy,1,N);
                S = repmat(vw1(:)',n2,1).*exp(-sqrt(ivimb*Delta)*vw2(:)');
                %S = repmat(vw1(:)',n2,1).*exp(-Delta*sqrt(ivimb/(Delta-delta/3))*vw2(:)');
                mu = sum((Y-S).^2,1)'/2/deviation^2;
                li = exp(-mu);
                %li = exp(-mu)/sum(li(:));
                li = reshape(li,size(vw1));
                ind = find(li==max(li(:)));
                ind = round(median(ind));
                [xindex,yindex]=ind2sub(size(vw1),ind);

                f=w1(yindex)/(S0+w1(yindex));
                V=w2V(xindex);
                calculated_values(i,:)=[(S0+w1(yindex)) f V Dp];
            else
                V=0;
                f=0;
                calculated_values(i,:)=[S0 f V Dp];
            end

        end
end
imags=zeros(size(sorted_data,1),4);
imags(values_location,:)=calculated_values(1:numel(values_location),:);
S0=reshape(imags(:,1),size(D));
f=reshape(imags(:,2),size(D));
Dstar=reshape(imags(:,3),size(D));
Dstar(Dstar==0.2)=0;
D=reshape(imags(:,4),size(D));
results=cat(ndimensions,S0,f,Dstar,D);
disp("Ended at " + datestr(datetime));

end
