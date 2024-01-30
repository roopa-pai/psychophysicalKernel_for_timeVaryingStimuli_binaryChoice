% CV2D
% (1) "fminunc" and "crossval" used to maximize log likelihood
% (2) smoothing along two dimensions, time and signal strength

clearvars
diary(sprintf('run_%s_%d.txt',datestr(now,'yyyy_mm_dd_HH_MM_SS'),randi([1,10000],1)));

% get cued and uncued data from experiment datafile and set up analysis structure("mat").
% "matsize" set to 2 because 2 conditions(cued and uncued) were analyzed.
load ki09012015.mat;
matsize = 2;
mat = struct('s', cell(1, matsize),'l', cell(1, matsize), 'l2', cell(1, matsize),'strength', cell(1, matsize),'ls', cell(1, matsize), 'y', cell(1, matsize),'st', cell(1, matsize),'CVW', cell(1, matsize),'PK_orig', cell(1, matsize), 'PK_model', cell(1, matsize));
mat(1).matrix = res.m_cued;
mat(2).matrix = res.m_uncued;

folder = '1 - CV2D section';
figfolder = strcat('C:\Users\roopa\Documents\caesardata\18.1.16\',folder,'\');
mkdir (figfolder);

for mat_i = 1
    
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
    
% convert decision vector: (-1,1) -> (1,0)
    raw_y = m(:,1);
    y = raw_y;
    cond = (y<0);
    y(cond) = 1;
    y(~cond) = 0;
    mat(mat_i).y = y;

    % separate stimulus and decision arrays by signal strength (1:ls), get unique stimulus values and psychophysical kernel (PK) for all
    clear st;
    st = struct('stim', cell(1, ls+1), 'dec', cell(1, ls+1),'u', cell(1, ls+1),'pk', cell(1, ls+1),'stim_c_temp', cell(1, ls+1),'stim_c', cell(1, ls+1), 'dec_c', cell(1, ls+1),'u_c', cell(1, ls+1),'pk_c', cell(1, ls+1), 't', cell(1, ls+1),'opt_L', cell(1, ls+1),'opt_ws', cell(1, ls+1),'opt_sigma', cell(1, ls+1),'opt_tau', cell(1, ls+1),'STM',cell(1, ls+1), 'modelpk', cell(1, ls+1),'model_mcr', cell(1, ls+1), 'CVW', cell(1, ls+1));
    for i = 4:8
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
        if i == 1
            st(i).stim_c_temp = st(i).stim;
            st(i).dec = y(ss);
            st(i).dec_c = st(i).dec;
            st(i).u = unique(st(i).stim(:,:));
            d = hist(st(i).stim',st(i).u);
            
            % PK ( = average number of times a particular stimulus appeared when
            % subject chose "near" - av. number of times it appeared when the
            % subject chose "far")
            if all(st(i).dec == 1)
                st(i).pk = mean(d(:,st(i).dec == 1),2);
            elseif all(st(i).dec == 0)
                st(i).pk = - mean(d(:,st(i).dec == 0),2);
            else
                st(i).pk = mean(d(:,st(i).dec == 1),2) - mean(d(:,st(i).dec == 0),2);
            end
            two_d = st(i).stim;
            for td_i = 1:ls %add a page of strength(i) on the page number corresponding to "i".
                if td_i == i
                    two_d(:,:,td_i) = st(i).stim;
                else
                    two_d(:,:,td_i) = zeros(size(st(i).stim));
                end
            end
            st(i).stim = two_d; %add second dimension to stimulus analysis. (first is time, second is stim strength) [separate]
            st(i).stim_c = two_d; %add second dimension to stimulus analysis. (first is time, second is stim strength) [combined]

        elseif 1<i<ls+1
            st(i).stim_c_temp = vertcat(st(i-1).stim_c_temp, st(i).stim);
            st(i).dec = y(ss);
            st(i).dec_c = vertcat(st(i-1).dec_c, st(i).dec);
            st(i).u = unique(st(i).stim(:,:));
            d = hist(st(i).stim',st(i).u);
            if all(st(i).dec == 1)
                st(i).pk = mean(d(:,st(i).dec == 1),2);
            elseif all(st(i).dec == 0)
                st(i).pk = - mean(d(:,st(i).dec == 0),2);
            else
                st(i).pk = mean(d(:,st(i).dec == 1),2) - mean(d(:,st(i).dec == 0),2);
            end
            two_d = st(i).stim;
            for td_i = 1:ls %add a page of strength(i) on the page number corresponding to "i".
                if td_i == i
                    two_d(:,:,td_i) = st(i).stim;
                else
                    two_d(:,:,td_i) = zeros(size(st(i).stim));
                end
            end
            st(i).stim = two_d; %add second dimension to stimulus analysis. (first is time, second is stim strength) [separate]
            st(i).stim_c = vertcat(st(i-1).stim_c, st(i).stim); %add second dimension to stimulus analysis. (first is time, second is stim strength) [combined]
        end
    end
    i = i + 1;
    st(i).stim = st(i-1).stim_c;
    st(i).dec = st(i-1).dec_c;
    st(i).u = unique(st(i).stim(:,:));
%     stim_section = st(i).stim;
%     stim_section = stim_section(:,:,:);
%     stim_dec = st(i).dec;
%     stim_dec = stim_dec(:,:,:);

    % find combination of hyperparameters (sigma and tau) and weights that give highest likelihood
    sigmas = [0.1;1;10];
    taus1 = [1;5];
    taus2 = 2;
    STM = zeros(length(sigmas), length(taus1), length(taus2));
    L = -Inf;
    for index_s = 1:length(sigmas)
        sigma = sigmas(index_s);
        for index_t1 = 1:length(taus1)
            tau1 = taus1(index_t1);
            for index_t2 = 1:length(taus2)
                tau2 = taus2(index_t2);
                [L_STM, ws_STM] = CVfutn (st(i).stim, st(i).dec, sigma, tau1, tau2);
                STM(index_s, index_t1, index_t2) = L_STM;
                disp('did at least one STM')
                if L_STM>L
                    st(i).opt_L = L_STM;
                    st(i).opt_ws = ws_STM';
                    st(i).opt_sigma = sigma;
                    st(i).opt_tau1 = tau1;
                    st(i).opt_tau2 = tau2;                        
                    L = L_STM;
                end
            end
        end
    end
    st(i).STM = STM;
    W = reshape(st(i).opt_ws,[l2 ls]);
    mat(mat_i).CVW = W;
    
    % simulate data
    
    %reps is how many times you use each of the actual trials to
    %simulate a response. i.e. reps = 2 means cycling through all
    %trials twice to give <2*trials> simulated trials and responses.
    reps = 1;
    
    for i = 1:ls
        [st(i).modelpk, st(i).model_mcr] = modelPKfutn (st(i).u, st(i).stim, st(i).dec, W(:,i), reps); 
    end
    
    mat(mat_i).st = st;

    % PLOTS
    % (1) plot W matrix
    figure('Visible','off');
    sc_x = [strength(1) strength(ls)];
    sc_y = [1 l2];
    imagesc(sc_x, sc_y, W)
    set(gca,'YDir','normal')
    colorbar;
    title(['CV 2D Weights (', caption,')'])
    xlabel ('signal strength')
    ylabel ('time (ms)')
    filenameExtensionW = strcat('CV 2D Weights(', caption, ')');
    saveas(gcf,strcat(figfolder, filenameExtensionW), 'bmp');
    
    % (2) plot W subplots
    figure ('Visible','off');
    fig = gcf;
    for i = 1:ls
        ncol = round(ls/2);
        subplot(2,ncol,i)
        plot(1:l2*16, W(:,i),'r');
        title(['CV 2D Weights: ',num2str(strength(i))])
        xlabel ('time (ms)')
        ylabel('weights')
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 17 4.5];
        filenameExtensionwsub = strcat('CV 2D Weights_subplots(', caption,')');
        saveas(gcf,strcat(figfolder, filenameExtensionwsub), 'bmp');
    end

    % (3) plot PK subplots
    figure('Visible','off');
    fig = gcf;
    for i = 1:ls
       ncol = round(ls/2);
       subplot(2,ncol,i)
       plot(st(i).u, st(i).pk,'k:'); hold on; % originalPK
       plot(st(i).u, st(i).modelpk, 'r--', 'LineWidth', 2); hold off; % modelPK
       title(['PK comparison(',num2str(strength(i)),')'])
       xlabel ('disparity ( ^{\circ} )')
       ylabel('PK amp. [frames/trial]')
       legend('k_o', 'k_mod')
    end
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 17 4.5];
    filenameExtensionpk = sprintf(['PK comparison (%s)', caption]);
    saveas(gcf,strcat(figfolder, filenameExtensionpk), 'bmp');
    
    % (4) plot PK (data vs model) comparison
    figure ('Visible','off');
    fig = gcf;
    PK1 = zeros(ls, length(st(i).u));
    PK2 = zeros(ls, length(st(i).u));
    for i = 1:ls
        PK1(i,:) = st(i).pk;
        PK2(i,:) = st(i).modelpk;
    end
    subplot(1,2,1)
    sc_x = [strength(1) strength(length(strength))];
    sc_y = [1 length(st(i).u)];
    imagesc(sc_x, sc_y, PK1)          % note: amplitude = #frames/trial (near-far)
    set(gca,'YDir','normal')
    colorbar;
    title(['CV PK: Data (', caption,')'])
    xlabel ('signal strength')
    ylabel ('disparity ( ^{\circ} )')
    subplot(1,2,2)
    sc_x = [strength(1) strength(length(strength))];
    sc_y = [1 length(st(i).u)];
    imagesc(sc_x, sc_y, PK2)
    set(gca,'YDir','normal')
    colorbar;
    title(['CV PK:  Model (', caption,')'])
    xlabel ('signal strength')
    ylabel ('disparity ( ^{\circ} )')
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 11 4];
    filenameExtensionPK12 = strcat('PK Data vs Model (',caption,')');
    saveas(gcf,strcat(figfolder, filenameExtensionPK12), 'bmp');

    % save analysis results as ".mat" file
    stname = sprintf(['analysis - ',folder,'.mat']);
    save([figfolder, stname],'mat')
end


%-----------------functions------------------------------------------------

function g = sigmoid(z)
% sigmoid function
    g = zeros(size(z));
    cond_sig = (z >= 0);
    g(cond_sig) = 1./(1+exp(-z(cond_sig)));        
    g(~cond_sig) = exp(z(~cond_sig))./(1+exp(z(~cond_sig)));
end

function[modelpk, model_mcr] = modelPKfutn (unique, stim, dec, finalweights, reps)
% gives: (1)misclassification rate of model (2) comparison of model & original PK
    stim_mcr = repmat(stim, reps, 1);
    dec_mcr = repmat(dec,reps,1);
    comb = repmat(stim*finalweights, reps, 1);
    p_mcr = sigmoid(comb);
    y_mcr = rand(size(p_mcr));
    cond_mcr1 = (y_mcr<=p_mcr);
    y_mcr(cond_mcr1) = 1;
    y_mcr(~cond_mcr1) = 0;
    yp_mcr = zeros(size(p_mcr)); 
    cond_mcr2 = (y_mcr == dec_mcr);
    yp_mcr(~cond_mcr2) = 1; % nnz if model decision is wrong
    d_mcr = hist(stim_mcr',unique);
    % model PK
    if all(y_mcr == 1)
        modelpk = mean(d_mcr(:,y_mcr == 1),2);
    elseif all(y_mcr == 0)
        modelpk = - mean(d_mcr(:,y_mcr == 0),2);
    else
        modelpk = mean(d_mcr(:,y_mcr == 1),2) - mean(d_mcr(:,y_mcr == 0),2);
    end
    % misclassification rate of model = number of mismatches/number of fake trials
    model_mcr = nnz(yp_mcr)/length(yp_mcr);
end

function [L_modelCV, ws_model_CV] = CVfutn (stim, dec, sigma, tau1, tau2)
% performs tenfold cross-validation of the likelihood for the chosen pair of sigma and tau values
% over all the trials.
    [~,sz2,sz3] = size(stim);
    ws_modelCV_vec = [];
    classf = @(XTRAIN, YTRAIN, XTEST, YTEST)(model(XTRAIN, YTRAIN, XTEST, YTEST, sigma, tau1, tau2));
    L_modelCV_vec = crossval(classf, stim, dec);
    disp('finished crossval')
    L_modelCV = sum(L_modelCV_vec)/10;
    ws_model_CV = sum(ws_modelCV_vec)/10; %averaged weights for all 10 subsets for 1 STT condition set

    function [L] = model(XTRAIN, YTRAIN, XTEST, YTEST, sigma, tau1, tau2)
    % model takes weights generated by using 90% of the trials ("XTRAIN")and tests it
    % on the remaining 10% ("XTEST")
        % XTRAIN already (automatically) arranged such that all stimuli for one trial is on the same row - pages aligned L to R
        U = calcB(sz2, tau1);
        disp('U done')
        Uinv = pinv(U);
        disp('Uinv done')
        V = calcB (sz3, tau2);
        disp('V done')
        Vinv = pinv(V);
        disp('Vinv done')
        k_X = randn(sz2,sz3); % <- initialize random 90x11 weights
        k_A = chol(U, 'upper');
        disp('k_A done')
        k_B = chol(V, 'lower');
        disp('k_B done')
        k_init = sigma*k_A*k_X*k_B; % <- include covariance info
        k_init = reshape(k_init, [sz2*sz3 1]); % <- make k_init a 990x1 vector. order doesn't matter because just initialization.
        % ws_model gives weights for one of the ten CV subsets for one STT combination.
        options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true, 'CheckGradients', false, 'FiniteDifferenceType','central');
        [ws_model, ~] = fminunc(@(w) (maxLLw(XTRAIN, w, YTRAIN, sigma, Uinv, Vinv)), k_init, options);
        disp('exited maxLLw')
        ws_model = reshape(ws_model, [sz2*sz3 1]);
        comb = XTEST*ws_model;
        p = sigmoid(comb);
        L = sum(YTEST.* log(p) + (1 - YTEST).* log(1 - p));
        ws_modelCV_vec = vertcat(ws_modelCV_vec, ws_model.');
     
        function[B] = calcB(nk, tau)
        % squared exponential covariance function
            i_B = (1:nk)';
            num = bsxfun(@minus,i_B, i_B').^2;
            B = exp(-0.5*num/tau.^2);
        end

        function g = sigmoid(z)
        % sigmoid function
        g = zeros(size(z));
        cond_sig = (z >= 0);
        g(cond_sig) = 1./(1+exp(-z(cond_sig)));        
        g(~cond_sig) = exp(z(~cond_sig))./(1+exp(z(~cond_sig)));
        end
  
        function[f,g] = maxLLw(stim, weights, y, sigma, Uinv, Vinv)
        % function ("f") and gradients ("g") used to minimize the
        % negative of the log likelihood
            comb_w = stim*weights;
            weights = reshape(weights, [sz2 sz3]);
            pw = (0.5*trace(Uinv*weights*Vinv*weights')/(sigma^2));
            pwg = (Uinv*weights*Vinv)/(sigma^2);
            pwg = reshape(pwg, [1 sz2* sz3]);
            f = -sum(y.*comb_w - comb_w - log(1 + exp(-comb_w)) - pw); % 1 value overall
            g = -sum((y -  1./(1+exp(-comb_w))).*stim - pwg); % each row/trial has 990 gradient values
        end
    end
end