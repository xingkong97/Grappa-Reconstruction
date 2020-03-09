%example for Grappa reconstruction
%Bo Li 1.29.2020
close all
clear

load ('kdata.mat');
load ('refscan.mat');
res=DoGrappa(KSpaceDATA,refscan);
subplot(2,2,1);
imshow(sqrt(sum(abs(fft2c(KSpaceDATA)),3)'),[]);
title('image before Grappa recon');
subplot(2,2,2);
imshow(sqrt(sum(abs(fft2c(res)),3)'),[]);
title('image after Grappa recon');
subplot(2,2,3);
imshow(abs(KSpaceDATA(:,:,1))',[]);
title('k space before Grappa recon');
subplot(2,2,4);
imshow(abs(res(:,:,1))',[]);
title('k space after Grappa recon');

function res=DoGrappa(kdata, refdata)%kdata:nx,ny,ncoil,refdata:nx,ny,ncoil
    kdata=permute(kdata,[2,1,3]);
    R=[1,2];%the undersampling rate of readout(nx) and phase(ny) directions
    kernel=[3,2];%one is odd and one is even
    kcabliData=permute(refdata,[2,1,3]);
    res=GrappaRecon(kdata, kcabliData, kernel, R);
    res=permute(res,[2,1,3]);
end

function res=GrappaRecon(kdata,kcabliData,kernel,R)
%Fast Grappa reconstruction for undersampling k-space
%Bo Li 1.29.2020
%kdata: undersampling k-space data
%kcabliData: calibration Data (full sampling data)
%kernel: grappa reconstruction kernel like [3,2], the second one is even
%R: the undersampling rate of readout and phase directions

pad=grappa_get_pad_size(kernel, R);
kdata_pad=padarray(kdata(:,:,1:end), [pad(2) pad(1)]);%pad(2) has to be even

kcabliData_pad=padarray(kcabliData(:,:,1:end), [pad(1) pad(1)]);
kcabliData_pad_1=permute(kcabliData_pad,[2,1,3]);

[sx,sy,sz] = size(kcabliData_pad_1);
winSize=[kernel(1),kernel(1)];

calib0 = zeros((sx-winSize(1)+1)*(sy-winSize(2)+1),prod(winSize),sz);
count=0;
for y=1:winSize(2)
    for x=1:winSize(1)
        count = count+1;
        calib0(:,count,:) = conj(reshape(kcabliData_pad_1(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:),...
            (sx-winSize(1)+1)*(sy-winSize(2)+1),1,sz));
    end
end

calib0_S=calib0(:,[1:kernel(1),end-kernel(1)+1:end],:);
S=reshape(calib0_S,size(calib0_S,1),size(calib0_S,2)*size(calib0_S,3));
calib0_T=calib0(:, floor(end/2+0.5),:);
T=reshape(calib0_T,size(calib0_T,1),size(calib0_T,2)*size(calib0_T,3));
W = pinv(S) * T;

kdata_pad_1=permute(kdata_pad,[2,1,3]);
[sx,sy,sz] = size(kdata_pad_1);
kdata0 = zeros((sx-winSize(1)+1)*(sy-winSize(2)+1),prod(winSize),sz);
count=0;
for y=1:winSize(2)
    for x=1:winSize(1)
        count = count+1;
        kdata0(:,count,:) = conj(reshape(kdata_pad_1(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:),...
            (sx-winSize(1)+1)*(sy-winSize(2)+1),1,sz));
    end
end

kdata0=reshape(kdata0,sx-winSize(1)+1,sy-winSize(2)+1,prod(winSize),sz);
kdata0=kdata0(:,2:2:end,:,:);
kdata0=reshape(kdata0,size(kdata0,1)*size(kdata0,2),prod(winSize),sz);
S_ch_new=kdata0(:,[1:3,7:9],:);
S=reshape(S_ch_new,size(S_ch_new,1),size(S_ch_new,2)*size(S_ch_new,3));
T_new=S * W;

T_ch_new = reshape(T_new,size(kdata, 2),size(T_new, 1)/size(kdata, 2),size(T_new, 2));
res = kdata;
for icoil=1:size(kdata, 3)
    res(1:2:end,:,icoil)=T_ch_new(:,1:end-1,icoil)';%2:end for even line undersampling, 1:end-1 for odd line undersampling
end

end

function pad = grappa_get_pad_size(kernel, R)
%   Compute size of padding needed in each direction
    pad =   floor(R.*kernel/2);
end


