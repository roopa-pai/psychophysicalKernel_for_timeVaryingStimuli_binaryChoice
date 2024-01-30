% CV 1D
% tau (length scale) in log scale, round 2
clearvars
diary(sprintf('run_%s_%d.txt',datestr(now,'yyyy_mm_dd_HH_MM_SS'),randi([1,10000],1)));

% get cued and uncued data from experiment datafile and set up analysis structure("mat").
% "matsize" set to 2 because 2 conditions(cued and uncued) were analyzed.
load ki09012015.mat;
matsize = 2;
mat = struct('s', cell(1, matsize),'l', cell(1, matsize), 'l2', cell(1, matsize),'strength', cell(1, matsize),'ls', cell(1, matsize), 'y', cell(1, matsize),'st', cell(1, matsize),'CVW', cell(1, matsize),'PK_orig', cell(1, matsize), 'PK_model', cell(1, matsize));
mat(1).matrix = res.m_cued;
mat(2).matrix = res.m_uncued;

folder = '1 - CV 1D, second round, 100 reps';
figfolder = strcat('C:\Users\roopa\Documents\caesardata\18.1.27\',folder,'\');
mkdir (figfolder);

for mat_i = 1:matsize
    
    %  caption for naming figures etc.
    if mat_i == 1
        caption = 'cued';
    elseif mat_i == 2
        caption = 'uncued';
    end
    
    m = mat(mat_i).matrix;
    s = m(:,4:93);  % stimulus sequence (disparities) for all trials
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
    st = struct('strength', cell(1, ls),'stim', cell(1, ls), 'dec', cell(1, ls),'u', cell(1, ls),'pk', cell(1, ls), 'opt_L', cell(1, ls),'opt_ws', cell(1, ls),'opt_sigma', cell(1, ls),'opt_tau1', cell(1, ls), 'modelpk', cell(1, ls),'model_mcr', cell(1, ls));
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
        if all(st(i).dec == 1)
            st(i).pk = mean(d(:,st(i).dec == 1),2);
        elseif all(st(i).dec == 0)
            st(i).pk = - mean(d(:,st(i).dec == 0),2);
        else
            st(i).pk = mean(d(:,st(i).dec == 1),2) - mean(d(:,st(i).dec == 0),2);
        end
        
        % set range of sigmas and taus to test. loop through all and print
        %map of log likelihood to find best st value pair. range of st in
        %this script determined by wider (log-2 to log2) previous run.
        
        if mat_i == 1
            if (i == 1)|| (i ==11)
                sigmas = 100:100:1000;
                taus1 = 1:1:10;
            elseif i == 2
                sigmas = 100:100:1000;
                taus1 = 0.01:0.01:0.1;
            elseif (i == 3)|| (i == 4)||(i == 5)|| (i == 6)||(i == 9)|| (i == 10)
                sigmas = 10:10:100;
                taus1 = 1:1:10;
            elseif i == 7
                sigmas = 1:1:10;
                taus1 = 1:1:10;
            elseif i == 8
                sigmas = 0.1:0.1:1;
                taus1 = 100:100:1000;
            end
        elseif mat_i == 2
            if i == 1
                sigmas = 0.01:0.01:0.1;
                taus1 = 0.1:0.1:1;       
            elseif i == 2
                sigmas = 0.01:0.01:0.1;
                taus1 = 0.01:0.01:0.1;      
            elseif i == 3
                sigmas = 0.1:0.1:1;
                taus1 = 0.1:0.1:1; 
            end
        end          
        STM = zeros(length(sigmas), length(taus1));
        L = -Inf;
        for index_s = 1:length(sigmas)
            sigma = sigmas(index_s);
            for index_t1 = 1:length(taus1)
                tau1 = taus1(index_t1);
                [L_STM, ws_STM] = CVfutn (st(i).stim, st(i).dec, sigma, tau1);
                STM(index_s, index_t1) = L_STM; 
                if L_STM>L
                    st(i).opt_L = L_STM;
                    st(i).opt_ws = ws_STM';
                    st(i).opt_sigma = sigma;
                    st(i).opt_tau1 = tau1;
                    L = L_STM;
                end
            end
        end
        
        figure;
        sc_x = [taus1(1) taus1(length(taus1))];
        sc_y = [sigmas(1) sigmas(length(sigmas))];
        imagesc(sc_x, sc_y, STM)
        set(gca,'YDir','normal')
        title(['Log likelihood (',num2str(strength(i)),') = ',num2str(st(i).opt_L),'; sigma = ',num2str(st(i).opt_sigma),', tau = ', num2str(st(i).opt_tau1)]);
        xlabel ('tau')
        ylabel ('sigma')
        colorbar;
        filenameExtensionL = strcat('Log likelihood map (', num2str(strength(i)),') (',num2str(caption),').bmp');
        saveas(gcf,strcat(figfolder, filenameExtensionL), 'bmp');
        
        %reps is how many times you use each of the actual trials to
        %simulate a response. i.e. reps = 2 means cycling through all
        %trials twice to give <2*trials> simulated trials and responses.
        reps = 100;
        [st(i).modelpk, st(i).model_mcr] = modelPKfutn (st(i).u, st(i).stim, st(i).dec, st(i).opt_ws, reps); 
      
        mat(mat_i).st = st;
    end
    
    % plot final weights as subplots
    figure ('Visible','off');
    fig = gcf;
    for i = 1:ls
        ncol = round(ls/2);
        subplot(2,ncol,i)
        plot(1:length(st(i).opt_ws), st(i).opt_ws);
        title(['Weights: ',num2str(strength(i))])
        xlabel ('timepoints')
        ylabel('weights')
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 17 4.5];
        filenameExtensionwsub = strcat('CV1D Weights subplots(', caption, ')');
        saveas(gcf,strcat(figfolder, filenameExtensionwsub), 'bmp');
    end
    
    % plot and display weights in imagesc matrix format
    figure('Visible','off');
    W = zeros(l2,ls);
    for i = 1:ls
        W(:,i) = st(i).opt_ws;
    end
    mat(mat_i).CVW = W;
    sc_x = [mat(mat_i).strength];
    sc_y = [1 l2];
    imagesc(sc_x, sc_y*16, W)
    set(gca,'YDir','normal')
    colorbar;
    title(['Weights smoothed over time (', caption,')'])
    xlabel ('signal strength')
    ylabel ('time (ms)')
    filenameExtensionW = strcat('CV1D Weights(', caption, ')');
    saveas(gcf,strcat(figfolder, filenameExtensionW), 'bmp');

    % plot PKs as subplots
    figure('Visible','off');
    fig = gcf;
    for i = 1:ls
       ncol = round(ls/2);
       subplot(2,ncol,i)
       plot(st(i).u, st(i).pk); hold on;
       plot(st(i).u, st(i).modelpk); hold off;   
       title(['PK: ',num2str(strength(i))])
       xlabel ('disparity ( ^{\circ} )')
       ylabel('PK amp. [frames/trial]')
       fig.PaperUnits = 'inches';
       fig.PaperPosition = [0 0 17 4.5];
       filenameExtensionPK = strcat('CV1D - PK subplots (',caption,')');
       saveas(gcf,strcat(figfolder, filenameExtensionPK), 'bmp');
    end

    % subplot original vs model PK
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
    title(['PK: Data (', caption,')'])
    xlabel ('signal strength')
    ylabel ('disparity ( ^{\circ} )')
    subplot(1,2,2)
    sc_x = [mat(mat_i).strength];
    sc_y = [st(i).u];
    imagesc(sc_x, sc_y, PK2)
    set(gca,'YDir','normal')
    colorbar;
    title(['PK: Model (', caption,')'])
    xlabel ('signal strength')
    ylabel ('disparity ( ^{\circ} )')
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 11 4];
    %note: amplitude = #frames/trial (near-far)
    filenameExtensionPK12 = sprintf([' CV1D - PK Data vs Model (',caption,')']);
    saveas(gcf,strcat(figfolder, filenameExtensionPK12), 'bmp');
end

% save analysis results as ".mat" file
stname = strcat('analysis - ',folder,'.mat');
save([figfolder, stname],'mat')

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

function [L_STM, ws_STM] = CVfutn (stim, dec, sigma, tau)
% does tenfold cross-validation of the likelihood for the chosen pair of sigma and tau values
%over all the trials.
    ws_STM_vec=[];
    classf = @(XTRAIN, YTRAIN, XTEST, YTEST)(model(XTRAIN, YTRAIN, XTEST, YTEST, sigma, tau));
    L_STM_vec = crossval(classf, stim, dec);
    L_STM = sum(L_STM_vec)/10;
    ws_STM = sum(ws_STM_vec)/10; %averaged weights for all 10 subsets for 1 ST pair

    function [L] = model(XTRAIN, YTRAIN, XTEST, YTEST, sigma, tau)
    % model takes weights generated by using 90% of the trials ("XTRAIN")and tests it
    % on the remaining 10% ("XTEST")
        nk = size(XTRAIN,2);
        B = calcB(nk, tau);
        Binv = pinv(B);
        k_init = mvnrnd(zeros(1,nk),sigma*B)';
        %ws_stmodel gives weights for one of the ten CV subsets for one ST pair.
        options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true, 'CheckGradients', false, 'FiniteDifferenceType','central');
        [ws_stmodel, ~] = fminunc(@(w) (maxLL(XTRAIN, w, YTRAIN, sigma, Binv)), k_init, options);
        comb = XTEST*ws_stmodel;
        p = sigmoid(comb);
        L = sum(YTEST.* log(p) + (1 - YTEST).* log(1 - p));
        ws_STM_vec = vertcat(ws_STM_vec, ws_stmodel.');
     
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
  
        function[f,g] = maxLL(stim, weights, y, sigma, Binv)
        % function ("f") and gradients ("g") used to minimise the
        % negative of the log likelihood
            comb_w = stim*weights;
            pw = (0.5*weights.'*Binv*weights)/(sigma^2);
            pwg = (weights.'*Binv)/(sigma^2);
            f = -sum(y.*comb_w - comb_w - log(1 + exp(-comb_w)) - pw);
            g = -sum((y -  1./(1+exp(-comb_w))).*stim - pwg);
        end
    end
end