clear;
rand('seed',0);
randn('seed',0);
Original_image_dir  =    'C:\Users\csjunxu\Desktop\ECCV2016\grayimages\';
fpath                     =   fullfile(Original_image_dir, '*.png');
im_dir                    =   dir(fpath);
im_num                    =   length(im_dir);

r=0; % radius of Gaussian Kernel is 3
beta_m=[0.00001 0.00002 0.00002 0.00002 0.00002];
%beta_m=[0.01 0.02 0.05 0.2 0.5]; % for full variation

tt=[];
pp=[];
tol=1e-4;
eta=1;

beta2=beta_m(2);

disp([num2str(beta2)]);

for nSig = [10 20]
    for SpikyRatio = [0.15 1.3]
        for Sample = 1:1
            imPSNR{Sample} = [];
            imSSIM{Sample} = [];
            for i = 1:im_num
                par.I                =   double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
                IMname = regexp(im_dir(i).name, '\.', 'split');
                IMname = IMname{1};
                % beta_m=[0.005 0.01 0.02 0.1 0.2];
                % %beta_m=[0.01 0.02 0.02 0.1 0.2]; %For full variation
                %     blu=blur_gauss_impulse(lena,r,0,0);
                
                t = cputime;
                %% add Gaussian noise
                randn('seed',0);
                par.nim = par.I + nSig*randn(size(par.I ));
                %% add spiky noise or "salt and pepper" noise 1
                rand('seed',Sample-1)
                par.nim = 255*imnoise(par.nim/255, 'salt & pepper', SpikyRatio);
                %% add spiky noise or "salt and pepper" noise 1
                %                 rand('seed',0)
                %                 par.nim  =   impulsenoise(par.nim,SpikyRatio,0);
                
                % tic;
                % rec=blu;
                % for i=1:4
                %     rec=acwmf(rec,lena,i);
                % end
                % toc;
                % ind=(rec~=blu);
                
                %                 tic;
                %                 [rec,ind] = adpmedft(par.nim,19);
                %                 toc;
                %                 ind=(rec~=par.nim)&((par.nim==255)|(par.nim==0));
                %                 rec(~ind)=par.nim(~ind);
                
                ind=false(size(par.I)); % for full variation.
                
                [rec1,v]=deblur_TV_L1_inc(par.nim,r,~ind,tol,beta2,eta);
                t=cputime-t;
                rec1(rec1>255)=255;rec1(rec1<0)=0;
                ps2 = csnr( rec1,par.I , 0, 0 );
                ss2 = cal_ssim( rec1,par.I , 0, 0 );
                disp(['The two phase method: PSNR = ', num2str(ps2),' time = ',num2str(t)]);
                tt=[tt t];
                %% output
                imPSNR{Sample} = [imPSNR{Sample} ps2];
                imSSIM{Sample}  = [imSSIM{Sample} ss2];
                imwrite(rec1/255, ['./results/Cai2010_GauSpi_' IMname '_' num2str(nSig) '_' num2str(SpikyRatio) '_' num2str(Sample) '.png']);
                fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name,ps2,ss2);
            end
            %% save output
            SmPSNR(Sample) = mean(imPSNR{Sample});
            SmSSIM(Sample) = mean(imSSIM{Sample});
            result = sprintf('Cai2010_AMF_GauSpi_%d_%2.2f.mat',nSig,SpikyRatio);
            save(result,'nSig','SpikyRatio','imPSNR','imSSIM','SmPSNR','SmSSIM');
        end
        %% save output
        mPSNR = mean(SmPSNR);
        mSSIM = mean(SmSSIM);
        result = sprintf('Cai2010_GauSpi_%d_%2.2f.mat',nSig,SpikyRatio);
        save(result,'nSig','SpikyRatio','mPSNR','mSSIM','imPSNR','imSSIM','SmPSNR','SmSSIM');
    end
end
