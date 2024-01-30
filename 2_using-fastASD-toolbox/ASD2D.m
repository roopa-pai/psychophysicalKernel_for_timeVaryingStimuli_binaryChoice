% ASD 2D modified
clearvars
diary(sprintf('run_%s_%d.txt',datestr(now,'yyyy_mm_dd_HH_MM_SS'),randi([1,10000],1)));

% get cued and uncued data from experiment datafile and set up analysis structure("mat").
% "matsize" set to 2 because 2 conditions(cued and uncued) were analyzed.
load ki09012015.mat;
matsize = 2;
mat = struct('matrix', cell(1, matsize),'s', cell(1, matsize),'l', cell(1, matsize), 'l2', cell(1, matsize),'strength', cell(1, matsize),'ls', cell(1, matsize), 'y', cell(1, matsize),'st', cell(1, matsize),'ASD2DW', cell(1, matsize),'PK_orig', cell(1, matsize), 'PK_model', cell(1, matsize));
mat(1).matrix = res.m_cued;
mat(2).matrix = res.m_uncued;

folder = '4 - ASD2D';
figfolder = strcat('C:\Users\roopa\Documents\caesardata\18.1.24\',folder,'\');
mkdir (figfolder);

for mat_i = 1:matsize
    
    %  caption for naming figures etc.
    if mat_i == 1
        caption = 'cued';
    elseif mat_i == 2
        caption = 'uncued';
    end
    
    m = mat(mat_i).matrix;
    s = m(:,4:93);  % stimulus sequence (disparities) on each trial
    mat(mat_i).s = s;
    l = size(mat(mat_i).s,1); % number of trials
    mat(mat_i).l = l;
    l2 = size(mat(mat_i).s,2); % number of timepoints per trial
    mat(mat_i).l2 = l2;
    
    % get unique signal strengths
    if mat_i == 1
        strength = unique(m(:,2));
    else
        strength = unique(m(:,3));
    end
    mat(mat_i).strength = strength;
    ls = length(mat(mat_i).strength);
    mat(mat_i).ls = ls;
    
    % get decision vector: (-1 is near, 1 is far)
    y = m(:,1);
    mat(mat_i).y = y;

    % separate stimulus and decision arrays by signal strength (1:ls), get unique stimulus values and psychophysical kernel (PK) for all
    clear st;
    st = struct('strength', cell(1, ls+1),'stim', cell(1, ls+1), 'dec', cell(1, ls+1),'u', cell(1, ls+1),'pk', cell(1, ls+1),'stim_c_temp', cell(1, ls+1),'stim_c', cell(1, ls+1), 'dec_c', cell(1, ls+1), 'k_asd', cell(1, ls+1),'sigma_ASD', cell(1, ls+1),'tau_ASD', cell(1, ls+1), 'asdstats', cell(1, ls+1), 'modelpk', cell(1, ls+1), 'model_mcr', cell(1, ls+1));
    for i = 1:ls
        if mat_i == 1
            ss = find(abs(m(:,2) - strength(i))<1e-5);
        elseif mat_i == 2
            ss = find(abs(m(:,3) - strength(i))<1e-5);
        end
        s_matrix = zeros(length(ss),l2);
        st(i).strength = strength(i);
        for j = 1:length(ss)
            s_matrix(j,:) = s(ss(j),:);
        end
        st(i).stim = s_matrix;
        if i == 1
            st(i).stim_c_temp = st(i).stim;
            st(i).dec = y(ss);
            st(i).dec_c = st(i).dec;
            st(i).u = unique(st(i).stim(:,:));
            d = hist(st(i).stim',st(i).u);
            if all(st(i).dec == -1)
                st(i).pk = mean(d(:,st(i).dec == -1),2);
            elseif all(st(i).dec == 1)
                st(i).pk = - mean(d(:,st(i).dec == 1),2);
            else
                st(i).pk = mean(d(:,st(i).dec == -1),2) - mean(d(:,st(i).dec == 1),2);
            end
            st(i).stim_c = st(i).stim;
            two_d = st(i).stim;
            two_d_c = st(i).stim;
            for td_i = 1:ls %add a page of strength(i) on the page number corresponding to "i".
                if td_i == i
                    two_d(:,:,td_i) = st(i).stim;
                else
                    two_d(:,:,td_i) = zeros(size(st(i).stim));
                end
            end
            st(i).stim = two_d; %add second dimension to stimulus analysis. (first is time, second is stim strength) [separate]
            st(i).stim_c = two_d;%add second dimension to stimulus analysis. (first is time, second is stim strength) [combined]

        elseif 1<i<ls+1
            st(i).stim_c_temp = vertcat(st(i-1).stim_c_temp, st(i).stim);
            st(i).dec = y(ss);
            st(i).dec_c = vertcat(st(i-1).dec_c, st(i).dec);
            st(i).u = unique(st(i).stim(:,:));
            d = hist(st(i).stim',st(i).u);
            if all(st(i).dec == -1)
                st(i).pk = mean(d(:,st(i).dec == -1),2);
            elseif all(st(i).dec == 1)
                st(i).pk = - mean(d(:,st(i).dec == 1),2);
            else
                st(i).pk = mean(d(:,st(i).dec == -1),2) - mean(d(:,st(i).dec == 1),2);
            end
            two_d = st(i).stim;
            for td_i = 1:ls % add a page of strength(i) on the page number corresponding to "i".
                if td_i == i
                    two_d(:,:,td_i) = st(i).stim;
                else
                    two_d(:,:,td_i) = zeros(size(st(i).stim));
                end
            end
            st(i).stim = two_d; % add second dimension to stimulus analysis. (first is time, second is stim strength) [separate]
            st(i).stim_c = vertcat(st(i-1).stim_c, st(i).stim); % add second dimension to stimulus analysis. (first is time, second is stim strength) [combined]
        end
    end
    i = i + 1;
    st(i).stim = st(i-1).stim_c;
    st(i).dec = st(i-1).dec_c;
    st(i).u = unique(st(i).stim(:,:));

    % Generate true filter vector k
    nks = [l2 ls];  % number of filter pixels along [cols, rows]
    nk = prod(nks); % total number of filter coeffs
    len = [1 1];  % length scale along each dimension
    rho = 1;  % marginal prior variance

    % Generate factored ASD prior covariance matrix in Fourier domain ->
    % creates my sigma*B
    [cdiag1,U1,wvec1] = mkcov_ASDfactored([len(1),1],nks(1)); % columns
    [cdiag2,U2,wvec2] = mkcov_ASDfactored([len(2),1],nks(2)); % rows
    nf1 = length(cdiag1); % number of frequencies needed
    nf2 = length(cdiag2);

    % Draw true regression coeffs 'k' by sampling from ASD prior ->
    %creates my k_init
    kh = sqrt(rho)*randn(nf1,nf2).*(sqrt(cdiag1)*sqrt(cdiag2)'); % Fourier-domain kernel
    fprintf('Filter has: %d pixels, %d significant Fourier coeffs\n',nk,nf1*nf2);

    % Inverse Fourier transform
    kim = U1*(U2*kh')'; % convert to space domain (as 2D image )
    k = kim(:);  % as vector

    % Make full covariance matrix (for inspection purposes only; will cause
    % out-of-memory error if filter dimensions too big!)
    C1 = U1*diag(cdiag1)*U1';
    C2 = U2*diag(cdiag2)*U2';
    Cprior = rho*kron(C2,C1);

    nsamps = size(st(i).stim,1);
    
    % st(i).stim needs to be reshaped so that all stimuli are in a 2D matrix.
    %so, weights from 1 to 91, then weights from back leaf/page, i.e. for ls.
    % => st(i).stim = reshape(st(i).stim, nsamps, nk)

    ASD2D_x = reshape(st(i).stim, nsamps, nk); % combined matrix of stimuli in 3D array
        
    figure ('Visible','off');
    subplot(221);
    imagesc(Cprior)
    set(gca,'YDir','normal')
    colorbar;
    title('prior covariance');
    mat(mat_i).Cprior = Cprior;
    
    % Compute ridge regression estimate 
    fprintf('\n...Running ridge regression with fixed-point updates...\n');
    % Sufficient statistics (old way of doing it, not used for ASD)
    dd.xx = ASD2D_x'*ASD2D_x;   % stimulus auto-covariance
    dd.xy = (ASD2D_x'*st(i).dec); % stimulus-response cross-covariance
    dd.yy = st(i).dec'*st(i).dec;   % marginal response variance
    dd.nx = nk;     % number of dimensions in stimulus
    dd.ny = nsamps;  % total number of samples

    % Run ridge regression using fixed-point update of hyperparameters
    maxiter = 100;
    kridge = autoRidgeRegress_fp(dd,maxiter);

    % Compute ASD estimate
    fprintf('\n\n...Running ASD_2D...\n');

    minlens = [1;1];  % minimum length scale along each dimension
    [kasd,asdstats] = fastASD(ASD2D_x,st(i).dec,nks,minlens);
    
    st(i).k_asd = kasd;
    st(i).asdstats = asdstats;
    st(i).sigma_ASD = st(i).asdstats.rho;
    st(i).tau_ASD = st(i).asdstats.len;

    %  ---- Make Plots ----

    subplot(222);
    sc_x = [mat(mat_i).strength];
    sc_y = [1 l2*16];
    imagesc(sc_x, sc_y, reshape(kridge,nks))
    set(gca,'YDir','normal')
    colorbar;
    title('ridge');
    xlabel ('signal strength')
    ylabel ('time (ms)')


    subplot(224);
    sc_x = [mat(mat_i).strength];
    sc_y = [1 l2*16];
    imagesc(sc_x, sc_y, reshape(kasd,nks))
    set(gca,'YDir','normal')
    colorbar;
    title(['ASD 2D Weights: ', caption]);
    xlabel ('signal strength')
    ylabel ('time (ms)')
    mat(mat_i).ASD2DW = reshape(kasd,nks);
    
    filenameExtensionpk = sprintf(['ASD2D_Weights, ridge, covariance; ',caption,'.bmp']);
    saveas(gcf,strcat(figfolder, filenameExtensionpk), 'bmp');
    
    
    % plot all final weights as subplots
    figure ('Visible','off');
    fig = gcf;
    ASD2DW = mat(mat_i).ASD2DW;
    for i = 1:ls
        ncol = round(ls/2);
        subplot(2,ncol,i)
        plot(1:length(ASD2DW(:,i)), ASD2DW(:,i));
        title(['ASD 2D Weights: ',num2str(strength(i))])
        xlabel ('timepoints')
        ylabel('weights')
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 17 4.5];
        filenameExtensionwsub = sprintf(['ASD 2D Weights_subplots(', caption,')']);
        saveas(gcf,strcat(figfolder, filenameExtensionwsub), 'jpeg');
    end
    
    % simulate data and plot PK1, PK2 matrices
    reps = 100;
    figure ('Visible','off');
    fig = gcf;
    PK1 = zeros(length(st(i).u), mat(mat_i).ls);
    PK2 = zeros(length(st(i).u), mat(mat_i).ls);
    for i = 1:mat(mat_i).ls
        PK1(:,i) = st(i).pk;
        stimpk = st(i).stim;
        [st(i).modelpk, st(i).model_mcr] = modelPKfutn (st(i).u, stimpk(:,:,i), st(i).dec, ASD2DW(:,i), reps);
        PK2(:,i) = st(i).modelpk;
    end 
    mat(mat_i).PK_orig = PK1;
    mat(mat_i).PK_model = PK2;
    subplot(1,2,1)
    sc_x = [mat(mat_i).strength];
    sc_y = [st(i).u];
    imagesc(sc_x, sc_y, PK1)
    set(gca,'YDir','normal')
    colorbar;
    title(['ASD 2D PK: Data (', caption,')'])
    xlabel ('signal strength')
    ylabel ('disparity ( ^{\circ} )')
    subplot(1,2,2)
    sc_x = [mat(mat_i).strength];
    sc_y = [st(i).u];
    imagesc(sc_x, sc_y, PK2)
    set(gca,'YDir','normal')
    colorbar;
    title(['ASD 2D PK:  Model (', caption,')'])
    xlabel ('signal strength')
    ylabel ('disparity ( ^{\circ} )')
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 11 4];
    % note: amplitude = #frames/trial (near-far)
    filenameExtensionPK12 = sprintf(['ASD 2D, PK Data vs Model (',caption,')']);
    saveas(gcf,strcat(figfolder, filenameExtensionPK12), 'jpeg');

    % Display facts about estimate
    ci = asdstats.ci;
    fprintf('\nHyerparam estimates (+/-1SD)\n-----------------------\n');
    fprintf('     l: %5.1f  %5.1f (+/-%.1f)\n',len(1),asdstats.len,ci(1));
    fprintf('   rho: %5.1f  %5.1f (+/-%.1f)\n',rho(1),asdstats.rho,ci(2));

   mat(mat_i).st = st; 
   
end 

% plot cued and uncued ASDD2D weights in one figure
figure ('Visible','off');
subplot(1,2,1)
sc_x = [mat(1).strength]; % iteration of mat_i has already ended. check prev files to see if this is why xaxis labeling in PK12 plot was messed up
sc_y = [1 l2*16];
imagesc(sc_x, sc_y, mat(1).ASD2DW)
set(gca,'YDir','normal')
colorbar;
title('ASD 2D Weights (cued)');
xlabel ('signal strength')
ylabel ('time (ms)')
subplot(1,2,2)
sc_x = [mat(2).strength];
sc_y = [1 l2*16];
imagesc(sc_x, sc_y, mat(2).ASD2DW)
set(gca,'YDir','normal')
colorbar;
title('ASD 2D Weights (uncued)');
xlabel ('signal strength')
ylabel ('time (ms)')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 4];
filenameExtensionpk = sprintf('ASD2D_Weights');
saveas(gcf,strcat(figfolder, filenameExtensionpk), 'bmp');


% plot Cpriors
figure ('Visible','off');
subplot(1,2,1)
imagesc(mat(1).Cprior)
set(gca,'YDir','normal')
colorbar;
title('ASD 2D Prior Covariance (cued)');
xlabel ('index')
ylabel ('index')
subplot(1,2,2)
imagesc(mat(2).Cprior)
set(gca,'YDir','normal')
colorbar;
title('ASD 2D Prior Covariance (uncued)');
xlabel ('index')
ylabel ('index')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 4];
filenameExtensionCp = sprintf('ASD2D _Cpriors');
saveas(gcf,strcat(figfolder, filenameExtensionCp), 'bmp');

% save workspace as ".mat" file
stname = sprintf(['analysis - ',folder,'.mat']);
save([figfolder, stname],'mat')

%-----------------functions------------------------------------------------
 
function[modelpk, model_mcr] = modelPKfutn (unique, stim, dec, finalweights, reps)
    stim_mcr = repmat(stim, reps, 1);
    dec_mcr = repmat(dec,reps,1);
    comb = repmat(stim*finalweights, reps, 1);
    p_mcr = sigmoid(comb);
    y_mcr = rand(size(p_mcr));
    cond_mcr1 = (y_mcr<=p_mcr);
    y_mcr(cond_mcr1) = -1;
    y_mcr(~cond_mcr1) = 1;
    yp_mcr = zeros(size(p_mcr)); 
    cond_mcr2 = (y_mcr == dec_mcr);
    yp_mcr(cond_mcr2) = 1; % nnz if model decision == real decision
    d_mcr = hist(stim_mcr',unique);
    if all(y_mcr == -1)
        modelpk = mean(d_mcr(:,y_mcr == -1),2);
    elseif all(y_mcr == 1)
        modelpk = - mean(d_mcr(:,y_mcr == 1),2);
    else
        modelpk = mean(d_mcr(:,y_mcr == -1),2) - mean(d_mcr(:,y_mcr == 1),2);
    end
    model_mcr = (length(yp_mcr) - nnz(yp_mcr))/length(yp_mcr); %misclassification rate of model = number of mismatches/number of fake trials
end

function g = sigmoid(z)
    g = 1./(1+exp(-z));        
end