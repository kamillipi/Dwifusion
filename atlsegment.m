function signal = atlsegment(data,nbvals)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
atlas=niftiread('ROI_MNI_V7.nii');
signal=zeros(170,nbvals);
for j = 1:nbvals
    vol=squeeze(data(:,:,:,j));
    for i = 1:170
        signal(i,j)=sum(vol(atlas==i),"all");
    end
end
end

