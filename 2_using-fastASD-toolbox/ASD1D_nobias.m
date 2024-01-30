% ASD 1D modified
clearvars
diary(sprintf('run_%s_%d.txt',datestr(now,'yyyy_mm_dd_HH_MM_SS'),randi([1,10000],1)));

% get cued and uncued data from experiment datafile and set up analysis structure("mat").
% "matsize" set to 2 because 2 conditions(cued and uncued) were analyzed.
load ki09012015.mat;
matsize = 2;
mat = struct('matrix', cell(1, matsize),'s', cell(1, matsize),'l', cell(1, matsize), 'l2', cell(1, matsize),'strength', cell(1, matsize),'ls', cell(1, matsize), 'y', cell(1, matsize),'st', cell(1, matsize),'ASD1DW', cell(1, matsize),'PK_orig', cell(1, matsize), 'PK_model', cell(1, matsize));
mat(1).matrix = res.m_cued;
mat(2).matrix = res.m_uncued;

folder = '3 - ASD1D';
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
    st = struct('strength', cell(1, ls),'stim', cell(1, ls), 'dec', cell(1, ls),'u', cell(1, ls),'pk', cell(1, ls), 't', cell(1, ls),'k_ASD1D', cell(1, ls),'ASD1D_sigma', cell(1, ls),'ASD1D_tau', cell(1, ls),'modelpk', cell(1, ls),'model_mcr', cell(1, ls));
    for i = 1:ls
        % ss: indexes of all trials of signal strength "i"
        if mat_i == 1
            ss = find(abs(m(:,2) - strength(i))<1e-5);
        elseif mat_i == 2
            ss = find(abs(m(:,3) - strength(i))<1e-5);
        end
        s_matrix = zeros(length(ss),l2);
        st(i).strength = strength(i);
        % s_matrix = collection of all trials of signal strength "i"
        for j = 1:length(ss)
            s_matrix(j,:) = s(ss(j),:);
        end
        
        st(i).stim = s_matrix;
        st(i).dec = y(ss);
        st(i).u = unique(st(i).stim(:,:));
        d = hist(st(i).stim',st(i).u);
        
        % PK ( = average number of times a particular stimulus appeared when
        % subject chose "near" - av. number of times it appeared when the
        % subject chose "far")
        if all(st(i).dec == -1)
            st(i).pk = mean(d(:,st(i).dec == -1),2);
        elseif all(st(i).dec == 1)
            st(i).pk = - mean(d(:,st(i).dec == 1),2);
        else
            st(i).pk = mean(d(:,st(i).dec == -1),2) - mean(d(:,st(i).dec == 1),2);
        end
        
        % run fast ASD 1D to get optimal sigma and tau
        sigma = 1;
        tau = 10;
        nk = size(st(i).stim,2);
        B = calcB(nk, tau);
        Binv = pinv(B);
        k_init = mvnrnd(zeros(1,nk),sigma*B)';
        st(i).t = 1:nk;
        minlen = 1;
        [st(i).k_ASD1D,ASD1Dstats] = fastASD(st(i).stim, st(i).dec, nk, minlen);
        st(i).ASD1D_sigma = ASD1Dstats.rho;
        st(i).ASD1D_tau = ASD1Dstats.len;
        
        % generate simulated data from model, get modelPK and compare to original
        % NOTE: reps is how many times you use each of the actual trials to
        % simulate a response. i.e. reps = 2 means cycling through all
        % trials twice to give <2*trials> simulated trials and responses.
        reps = 100;
        [st(i).modelpk, st(i).model_mcr] = modelPKfutn (st(i).u, st(i).stim, st(i).dec, st(i).k_ASD1D, reps);
    
        mat(mat_i).st = st;
    end
    
    % plot all final weights as subplots
    figure ('Visible','off');
    fig = gcf;
    for i = 1:ls
        ncol = round(ls/2);
        subplot(2,ncol,i)
        plot(1:length(st(i).k_ASD1D), st(i).k_ASD1D);
        title(['ASD 1D Weights: ',num2str(strength(i))])
        xlabel ('timepoints')
        ylabel('weights')
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 17 4.5];
        filenameExtensionwsub = sprintf(['ASD 1D Weights_subplots(', caption,')']);
        saveas(gcf,strcat(figfolder, filenameExtensionwsub), 'jpeg');
    end

    % plot ASD1D W
    figure('Visible','off');
    W = zeros(nk,ls);
    for i = 1:ls
        W(:,i) = st(i).k_ASD1D;
    end
    mat(mat_i).ASD1DW = W;
    sc_x = [mat(mat_i).strength];
    sc_y = [1 nk];
    imagesc(sc_x, sc_y*16, W)
    set(gca,'YDir','normal')
    colorbar;
    title(['ASD 1D Weights (', caption,')'])
    xlabel ('signal strength')
    ylabel ('time (ms)')
    filenameExtensionW = sprintf(['ASD 1D Weights(', caption,')']);
    saveas(gcf,strcat(figfolder, filenameExtensionW), 'jpeg');

    % plot PKs on one graph
    figure('Visible','off');
    fig = gcf;
    for i = 1:ls
       ncol = round(ls/2);
       subplot(2,ncol,i)
       plot(st(i).u, st(i).pk); hold on;
       plot(st(i).u, st(i).modelpk); hold off;   
       title(['ASD 1D PK: ',num2str(strength(i))])
       xlabel ('disparity ( ^{\circ} )')
       ylabel('PK amp. [frames/trial]')
       fig.PaperUnits = 'inches';
       fig.PaperPosition = [0 0 17 4.5];
       filenameExtensionPK = sprintf(['ASD 1D - PKs_subplots (',caption,')']);
       saveas(gcf,strcat(figfolder, filenameExtensionPK), 'jpeg');
    end

    % plot PK1, PK2 matrices
    figure ('Visible','off');
    fig = gcf;
    PK1 = zeros(length(st(i).u), ls);
    PK2 = zeros(length(st(i).u), ls);
    for i = 1:ls
        PK1(:,i) = st(i).pk;
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
    title(['ASD 1D PK: Data (', caption,')'])
    xlabel ('signal strength')
    ylabel ('disparity ( ^{\circ} )')
    subplot(1,2,2)
    sc_x = [mat(mat_i).strength];
    sc_y = [st(i).u];
    imagesc(sc_x, sc_y, PK2)
    set(gca,'YDir','normal')
    colorbar;
    title(['ASD 1D PK:  Model (', caption,')'])
    xlabel ('signal strength')
    ylabel ('disparity ( ^{\circ} )')
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 11 4];
    %note: amplitude = #frames/trial (near-far)
    filenameExtensionPK12 = sprintf(['PK Data vs Model (',caption,')']);
    saveas(gcf,strcat(figfolder, filenameExtensionPK12), 'jpeg');

end

% save analysis results as ".mat" file
stname = sprintf(['analysis - ',folder,'.mat']);
save([figfolder, stname],'mat')


%-----------------functions------------------------------------------------

function[B] = calcB(nk, tau)
    i_B = (1:nk)';
    num = bsxfun(@minus,i_B, i_B').^2;
    B = exp(-0.5*num/tau.^2);
end

function g = sigmoid(z)
    g = zeros(size(z));
    cond_sig = (z >= 0);
    g(cond_sig) = 1./(1+exp(-z(cond_sig)));        
    g(~cond_sig) = exp(z(~cond_sig))./(1+exp(z(~cond_sig)));
end

% (1)misclassification rate of model (2) comparison of model & original PK
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
