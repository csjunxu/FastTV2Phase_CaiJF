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

for nSigG = 10
    for nSigL = [15 30 45]
        for Sample = 1:1
            imPSNR{Sample} = [];
            imSSIM{Sample} = [];
            for i = 1:im_num
                par.I = double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
                IMname = regexp(im_dir(i).name, '\.', 'split');
                IMname = IMname{1};
                % beta_m=[0.005 0.01 0.02 0.1 0.2];
                % %beta_m=[0.01 0.02 0.02 0.1 0.2]; %For full variation
                %     blu=blur_gauss_impulse(lena,r,0,0);
                t = cputime;
                %% add Gaussian noise
                randn('seed',Sample-1);
                par.nim = par.I + nSigG*randn(size(par.I ));
                %% add Laplacian noise
                rand('seed',Sample-1)
                par.nim = par.nim + nSigL*randl(size(par.I ));
                %% add "salt and pepper" noise : 0 or RVIN noise : 1
                %                 rand('seed',Sample-1)
                %                 [par.nim,Narr]          =   impulsenoise(par.nim,SpikyRatio,1);
                [rec,ind] = adpmedft(par.nim,19);
                ind=(rec~=par.nim)&((par.nim==255)|(par.nim==0));
                rec(~ind)=par.nim(~ind);
                
                %ind=false(size(lena)); % for full variation.
                
                [rec1,v]=deblur_TV_L1_inc(rec,r,~ind,tol,beta2,eta);
                t=cputime-t;
                rec1(rec1>255)=255;rec1(rec1<0)=0;
                ps2 = csnr( rec1,par.I , 0, 0 );
                ss2 = cal_ssim( rec1,par.I , 0, 0 );
                disp(['The two phase method: PSNR = ', num2str(ps2),' time = ',num2str(t)]);
                tt=[tt t];
                %% output
                imPSNR{Sample} = [imPSNR{Sample} ps2];
                imSSIM{Sample}  = [imSSIM{Sample} ss2];
                imwrite(rec1/255, ['C:/Users/csjunxu/Desktop/ECCV2016/1_Results/Caietal/GauLap/Cai2010_AMF_GauLap_' IMname '_' num2str(nSigG) '_' num2str(nSigL) '.png']);
                fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name,ps2,ss2);
            end
            %% save output
            SmPSNR(Sample) = mean(imPSNR{Sample});
            SmSSIM(Sample) = mean(imSSIM{Sample});
            result = sprintf('C:/Users/csjunxu/Desktop/ECCV2016/1_Results/Caietal/Cai2010_AMF_GauLap_%d_%d.mat',nSigG,nSigL);
            save(result,'imPSNR','imSSIM','SmPSNR','SmSSIM');
        end
        %% save output
        mPSNR = mean(SmPSNR);
        mSSIM = mean(SmSSIM);
        result = sprintf('C:/Users/csjunxu/Desktop/ECCV2016/1_Results/Caietal/Cai2010_AMF_GauLap_%d_%d.mat',nSigG,nSigL);
        save(result,'mPSNR','mSSIM','imPSNR','imSSIM','SmPSNR','SmSSIM');
    end
end
