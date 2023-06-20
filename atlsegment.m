function signal = atlsegment(data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
atlas=niftiread('ROI_MNI_V7.nii');
data_size=size(data,1:4);
if data_size(1:3)~=size(atlas)
    atlas=imresize3(atlas,data_size(1:3),method="nearest");
end
signal=zeros(170,data_size(4));
for j = 1:data_size(4)
    vol=squeeze(data(:,:,:,j));
    for i = 1:170
        signal(i,j)=sum(vol(atlas==i),"all");
    end
end
end

