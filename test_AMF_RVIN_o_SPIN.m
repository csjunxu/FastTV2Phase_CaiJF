clear;
Original_image_dir  =    'C:\Users\csjunxu\Desktop\TWSCGIN\cleanimages\';
Sdir = regexp(Original_image_dir, '\', 'split');
fpath                     =   fullfile(Original_image_dir, '*.png');
im_dir                    =   dir(fpath);
im_num                    =   length(im_dir);

method = 'TVL1';
write_MAT_dir = ['C:/Users/csjunxu/Desktop/TWSCGIN/'];
write_sRGB_dir = [write_MAT_dir method];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end


r=0; % radius of Gaussian Kernel is 3
% beta_m=[0.00001 0.00002 0.00002 0.00002 0.00002];
beta_m = [0.1, 0.3, 0.3, 0.3, 0.3]; % for full variation

tt=[];
pp=[];
tol=1e-4;
eta=1;

beta2=beta_m(2);

disp([num2str(beta2)]);
Type = 0;
for nSig = [10 20 30]
    for sp = [0.1 0.3 0.5]
        imPSNR = [];
        imSSIM = [];
        for i = 1:im_num
            Par.I                =   double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
            IMname = regexp(im_dir(i).name, '\.', 'split');
            IMname = IMname{1};
            % beta_m=[0.005 0.01 0.02 0.1 0.2];
            % %beta_m=[0.01 0.02 0.02 0.1 0.2]; %For full variation
            %     blu=blur_gauss_impulse(lena,r,0,0);
            
            t = cputime;
            %% add Gaussian noise
            %                 randn('seed', 0);
            %                 par.nim = par.I + nSig*randn(size(par.I ));
            %% add "salt and pepper" noise
            %                 rand('seed', 0)
            %                 par.nim = 255*imnoise(par.nim/255, 'salt & pepper', SpikyRatio);
            %% add "salt and pepper" noise : 0 or RVIN noise : 1
            %                 rand('seed',Sample-1)
            %                 [par.nim,Narr]          =   impulsenoise(par.nim,SpikyRatio,1);
            if Type == 0
                imname = sprintf([write_MAT_dir 'noisyimages/G' num2str(nSig) '_SPIN' num2str(sp) '_' im_dir(i).name]);
                %                                     imwrite(Par.nim/255,imname);
                Par.nim = double( imread(imname));
            elseif Type == 1
                imname = sprintf([write_MAT_dir 'noisyimages/G' num2str(nSig) '_RVIN' num2str(sp) '_' im_dir(i).name]);
                %                                     imwrite(Par.nim/255,imname);
                Par.nim = double( imread(imname));
            else
                break;
            end
            [rec,ind] = adpmedft(Par.nim,19);
            %             ind=(rec~=par.nim)&((par.nim==255)|(par.nim==0));
            %             rec(~ind)=par.nim(~ind);
            
            %ind=false(size(lena)); % for full variation.
            
            [rec1,v]=deblur_TV_L1_inc(rec,r,~ind,tol,beta2,eta);
            t=cputime-t;
            rec1(rec1>255)=255;rec1(rec1<0)=0;
            ps2 = csnr( rec1,Par.I , 0, 0 );
            ss2 = cal_ssim( rec1,Par.I , 0, 0 );
            disp(['The two phase method: PSNR = ', num2str(ps2),' time = ',num2str(t)]);
            tt=[tt t];
            %% output
            imPSNR = [imPSNR ps2];
            imSSIM  = [imSSIM ss2];
            imwrite(rec1/255, [write_sRGB_dir '/Cai2010_AMF2_GSPIN_' IMname '_' num2str(nSig) '_' num2str(sp) '.png']);
            fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name,ps2,ss2);
        end
        %% save output
        mPSNR = mean(imPSNR);
        mSSIM = mean(imSSIM);
        result = sprintf([write_MAT_dir '/' method '_AMF2_GSPIN_p_' Sdir{end-1} '_nSig' num2str(nSig) '_sp' num2str(sp) '.mat']);
        save(result,'nSig','sp','imPSNR','imSSIM','mPSNR','mSSIM');
    end
end
