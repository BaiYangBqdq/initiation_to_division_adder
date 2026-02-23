close all ;
clear all ;


alist(1:20) = 1:20;
% lambdalist = 0.3*ones(size(alist));

DTlist = 30*ones(size(alist));
Clist = 40*ones(size(alist));
% DTlist = round(log(2)./lambdalist * 60);
% Clist = round((log(2)./lambdalist*60).*(0.8801*lambdalist+0.2375));
% alpha_m = 0.28;       beta_m=0.99;      m0 = 1;
% mi=m0/log(2)*(alpha_m+beta_m*lambdalist).*exp(-alpha_m-beta_m*lambdalist);
% md=m0/log(2)*(0.28+0.99*lambdalist);
dVidlist = 4*ones(size(alist));

exp_lambda = [1.60 1.27 1.49 1.14 0.8 0.66 0.68 0.59 0.25 0.44];
exp_kdeltaid = [0.49 0.47 0.36 0.37 0.34 0.30 0.29 0.20 0.13 0.06];
auto_corr_list = zeros(size(DTlist));

% h0 = figure;
% plot(lambdalist,auto_corr_list,'-','linewidth',2);  hold all
% scatter(exp_lambda,exp_kdeltaid ,50,'filled'); hold all
% legend('exp data','fitted');
% set(gca,'fontsize',14);
% xlim([0 2]);
% ylim([0 1]);
% xlabel('growth rate (hr^-^1)');
% ylabel('auto correlation of \delta_i_d');
% drawnow;
%%

ii_sizer_slope = nan(size(DTlist));
dd_sizer_slope = nan(size(DTlist));
id_sizer_slope = nan(size(DTlist));
ii_adder_slope = nan(size(DTlist));
dd_adder_slope = nan(size(DTlist));
id_adder_slope = nan(size(DTlist));
ii_timer_slope = nan(size(DTlist));
dd_timer_slope = nan(size(DTlist));
id_timer_slope = nan(size(DTlist));
rlist = nan(size(DTlist));

MI = nan(size(DTlist));
for ir = 1:length(DTlist)

    tic;

    rep_Vi1 = [];
    rep_Vi2 = [];
    rep_dVi = [];
    rep_dTi = [];

    rep_Vd1 = [];
    rep_Vd2 = [];
    rep_dVd = [];
    rep_dTd = [];

    rep_dVid = [];
    rep_dTid = [];

    for rep = 1:100
    %% global parameters

    DT = DTlist(ir) ; % doubling rate (min)
    TT = DT*100 ; % total calculation time (min)
    dt = 0.2 ; % calculation time step (min)

    gamma0 = log(2)/DT ; % cell mass growth rate (min^-1)
    epsilon_A = 0.1;
    epsilon_V = 0;
    epsilon_id = 0;
    auto_corr = auto_corr_list(ir);

    gamma = gamma0*(epsilon_V * randn(1)+1);
    Cperiod = round((log(2)./gamma).*(0.8801*gamma*60+0.2375));

    % Cperiod = Clist(ir);
    alpha0 = alist(ir)*gamma*60;
    % alpha0 = alist(ir)*lambdalist(ir);

    datAn = 250 ; % dnaA binding site on datA
    Chns = 150 ; % number of hns binding to each chromosome

    DNAlength = ceil(Cperiod/dt) ; % DNA length is set to fix V fork = 1 ;
    dVidc0 = dVidlist(ir);


    %% initial condition
    DNAcopy = ones(1,DNAlength); % starting DNAcopy = 1 ;
    DNAcopy(1,1) = 1 ;
    DNA = 1;
    alpha = alpha0;

    Nd = 1 ; % copy number of datA
    nrep = 1 ; % number of DNA replication. starting from 1, not 0
    V = 5 ; % initial cell mass
    A = 1 ; % initial dnaA protein
    nd = 1;

    Abox = Nd*datAn + sum(DNAcopy(1:round(0.8*Cperiod/dt)))/numel(DNAcopy(1:round(0.8*Cperiod/dt)))*Chns;
    % Abox = Nd*datAn + DNA*Chns;

    rpt = NaN(1,ceil(TT/DT)+10) ; % DNA replication time
    tert = NaN(1,ceil(TT/DT)+10)  ; % DNA termination time
    divt =  NaN(1,ceil(TT/DT)+10)  ; % division time
    Prep = NaN(1,ceil(TT/DT)+10)  ; % position of DNA replication  
    Vi =  NaN(1,ceil(TT/DT)+10)  ; % initiation mass
    Vt =  NaN(1,ceil(TT/DT)+10)  ; % termination mass
    Vd =  NaN(1,ceil(TT/DT)+10)  ; % division mass

    % gammalist =  NaN(1,ceil(TT/DT)+10)  ; % growth rate
    Aboxlist =  NaN(1,ceil(TT/DT)+10) ; % threshold at initiation
    dVid_set =  NaN(1,ceil(TT/DT)+10) ; %
    dVid_set(1) =  dVidc0 ; %
    OriC =  NaN(1,ceil(TT/DT)+10) ; % 
    Num =  NaN(1,ceil(TT/DT)+10) ; % 
    % gamma = gamma0;


    OriC(1) = 1;
    Num(1) = 1; % cell number

    div_mass = [inf];

    %% main loop simulation

    for t = 1 : TT/dt
        DNA = sum(DNAcopy)/DNAlength ;
        Nd = DNAcopy(1); % copy number of datA site
        A = A + alpha*V*dt;
        V = V + gamma*V*dt ;
    %     Abox = Nd*datAn + DNA*Chns;
        Abox = Nd*datAn + sum(DNAcopy(1:round(0.8*Cperiod/dt)))/numel(DNAcopy(1:round(0.8*Cperiod/dt)))*Chns;

        % check new replication initiation condition
        if A>Abox  
            nrep = nrep +1; % replication number +1
            Prep(nrep) = 1 ; % generate new replication position
            rpt(nrep) = t*dt; % replication time +dt
            Vi(nrep) = V; % note initiation mass
    %         gammalist(nrep) = gamma; % 
    %         gamma = gamma0*(epsilon_V * randn(1)+1);
    %         gamma = gamma0;
            alpha = alpha0*(epsilon_A * randn(1)+1);
            dVid_set(nrep) = (1-auto_corr)*dVidc0 + auto_corr*dVid_set(nrep-1) + dVidc0*(epsilon_id * randn(1));
            OriC(nrep) = DNAcopy(1);
            div_mass(nrep) = Vi(nrep) + dVid_set(nrep)*OriC(nrep);

        end
        % DNA replication
        for f = 1 : length(Prep) 
            if Prep(f) <= Cperiod/dt
                DNAcopy(Prep(f)) = DNAcopy(Prep(f))*2 ; % all DNA replicates
                Prep(f) = Prep(f) +1 ; % replication position moves
            else if Prep(f) > Cperiod/dt % DNA replication terminates
                Prep(f) = NaN ;   
                tert(nrep) = t*dt ; % DNA termination time
                Vt(nrep) = V/DNAcopy(end) ; % termination mass
                end
            end
        end

        % Division
        if any(V >= div_mass) 
            nd = nd +1;
            Vd(nd) = V;
            Num(nd) = Num(nd-1)*2;
            divt(nd) = t*dt;
            div_mass(V>div_mass) = inf;
        end

    end


    %% arrange data
    OpC = OriC(nrep)/Num(nd)*2;

    Vi1 = Vi(10:nd-1)./OriC(10:nd-1);
    Vi2 = Vi(11:nd)./OriC(11:nd);
    dVi = 2*Vi2-Vi1;
    dTi = rpt(11:nd)-rpt(10:nd-1);

    Vd1 = Vd(10:nd-1)./Num(10:nd-1)*2;
    Vd2 = Vd(11:nd)./Num(11:nd)*2;
    dVd = Vd2*2-Vd1;
    dTd = divt(11:nd)-divt(10:nd-1);

    OriC = OriC(10:nd-1);
    Num = Num(10:nd-1);


    dVid = Vd1 - Vi1;
    dTid = divt(10:nd-1) - rpt(10:nd-1);


    rep_Vi1 = [rep_Vi1 Vi1];
    rep_Vi2 = [rep_Vi2 Vi2];
    rep_dVi = [rep_dVi dVi];
    rep_dTi = [rep_dTi dTi];

    rep_Vd1 = [rep_Vd1 Vd1];
    rep_Vd2 = [rep_Vd2 Vd2];
    rep_dVd = [rep_dVd dVd];
    rep_dTd = [rep_dTd dTd];

    rep_dVid = [rep_dVid dVid];
    rep_dTid = [rep_dTid dTid];

    end
    %%
    scaleVi1 = rep_Vi1/nanmean(rep_Vi1);
    scaleVi2 = rep_Vi2/nanmean(rep_Vi2);
    scaledVi = rep_dVi/nanmean(rep_dVi);
    scaledTi = rep_dTi/nanmean(rep_dTi);
    scaleVd1 = rep_Vd1/nanmean(rep_Vd1);
    scaleVd2 = rep_Vd2/nanmean(rep_Vd2);
    scaledVd = rep_dVd/nanmean(rep_dVd);
    scaledTd = rep_dTd/nanmean(rep_dTd);

    scaledVid = rep_dVid/mean(rep_dVid);
    scaledTid = rep_dTid/mean(rep_dTid);

    %% plot initiation correlation
    bin = 0.6:0.1:1.4;

    h1 = figure('position',[100 100 900 300]);
    subplot(1,3,1)
    disc = discretize(scaleVi1,bin);
    meanVi2 = nan(size(bin));    stdVi2 = nan(size(bin));
    for i = 1:length(bin)
        meanVi2(i) = nanmean(scaleVi2(disc==i));
        stdVi2(i) = nanstd(scaleVi2(disc==i));
    end
    scatter(scaleVi1,scaleVi2,10,'filled'); hold all
    plot(bin,meanVi2,'ro-','linewidth',2);
    xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);
    pa=polyfit(scaleVi1,scaleVi2,1);
    isla = pa(1);
    xlabel('Rescaled l_i^n');
    ylabel('Rescaled l_i^n^+^1');
    title(sprintf('slope = %.2f',isla));
    set(gca,'fontsize',14);

    subplot(1,3,2)
    meandV = nan(size(bin));    stddV = nan(size(bin));
    for i = 1:length(bin)
        meandV(i) = nanmean(scaledVi(disc==i));
        stddV(i) = nanstd(scaledVi(disc==i));
    end
    scatter(scaleVi1,scaledVi,10,'filled'); hold all
    plot(bin,meandV,'ro-','linewidth',2);
    xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);
    pb=polyfit(scaleVi1,scaledVi,1);
    islb = pb(1);
    xlabel('Rescaled l_i');
    ylabel('Rescaled \deltal_i');
    title(sprintf('slope = %.2f',islb));
    set(gca,'fontsize',14);

    subplot(1,3,3)
    meandTi = nan(size(bin));    stddTi = nan(size(bin));
    for i = 1:length(bin)
        meandTi(i) = nanmean(scaledTi(disc==i));
        stddTi(i) = nanstd(scaledTi(disc==i));
    end
    scatter(scaleVi1,scaledTi,10,'filled'); hold all
    plot(bin,meandTi,'ro-','linewidth',2);
    xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);
    p=polyfit(scaleVi1,scaledTi,1);
    islc = p(1);
    xlabel('Rescaled l_i');
    ylabel('Rescaled \tau_i');
    title(sprintf('slope = %.2f',islc));
    set(gca,'fontsize',14);
    drawnow;

    %% ii adder plot for pub
    figure;
    finbin = 0.7:0.05:1.3;
    findisc = discretize(scaleVi1,finbin);
    finbin_cent = (finbin(1:end-1)+finbin(2:end))/2;
    finbin_meandV = nan(size(finbin_cent));    % stddV = nan(size(finbin_cent));
    for i = 1:length(finbin_cent)
        finbin_meandV(i) = nanmean(scaledVi(findisc==i));
    %     stddV(i) = nanstd(scaledVi(disc==i));
    end
    scatter(scaleVi1,scaledVi,10,'filled'); hold all
    plot(finbin_cent,finbin_meandV,'ro-','linewidth',2);
    xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);
    pb=polyfit(scaleVi1,scaledVi,1);
    islb = pb(1);
    xlabel('Rescaled l_i');
    ylabel('Rescaled \deltal_i');
    title(sprintf('slope = %.2f',islb));

    set(gca,'fontsize',14);
    axis('square');

    %% division correlation

    bin = 0.6:0.1:1.4;
    h2 = figure('position',[100 400 900 300]); 

    subplot(1,3,1)
    disc = discretize(scaleVd1,bin);
    for i = 1:length(bin)
        meanscaleVd2(i) = nanmean(scaleVd2(disc==i));
        stdscaleVd2(i) = nanstd(scaleVd2(disc==i));
    end
    scatter(scaleVd1,scaleVd2,10,'filled'); hold all
    plot(bin,meanscaleVd2,'ro-','linewidth',2);
    xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);
    pa=polyfit(scaleVd1,scaleVd2,1);
    dsla = pa(1);
    xlabel('Rescaled l_d^n');
    ylabel('Rescaled l_d^n^+^1');
    title(sprintf('slope = %.2f',dsla));
    set(gca,'fontsize',14);

    subplot(1,3,2)
    disc = discretize(scaleVd1,bin);
    for i = 1:length(bin)
        meanscaledVd(i) = nanmean(scaledVd(disc==i));
        stddscaledVd(i) = nanstd(scaledVd(disc==i));
    end
    scatter(scaleVd1,scaledVd,10,'filled'); hold all
    plot(bin,meanscaledVd,'ro-','linewidth',2);
    xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);
    pb=polyfit(scaleVd1,scaledVd,1);
    dslb = pb(1);
    xlabel('Rescaled l_d');
    ylabel('Rescaled \deltal_d');
    title(sprintf('slope = %.2f',dslb));
    set(gca,'fontsize',14);

    subplot(1,3,3)
    disc = discretize(scaleVd1,bin);
    for i = 1:length(bin)
        meanscaledTd(i) = nanmean(scaledTd(disc==i));
        stdscaledTd(i) = nanstd(scaledTd(disc==i));
    end
    scatter(scaleVd1,scaledTd,10,'filled'); hold all
    plot(bin,meanscaledTd,'ro-','linewidth',2);
    xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);
    p=polyfit(scaleVd1,scaledTd,1);
    dslc = p(1);
    xlabel('Rescaled l_d');
    ylabel('Rescaled \tau_d');
    title(sprintf('slope = %.2f',dslc));
    set(gca,'fontsize',14);

    drawnow;

    %% initiation to division correlation
    bin = 0.6:0.1:1.4;
    h3 = figure('position',[900 100 900 300]); 
    subplot(1,3,1)
    disc = discretize(scaleVi1,bin);
    for i = 1:length(bin)
        meanscaleVd1(i) = nanmean(scaleVd1(disc==i));
    end
    scatter(scaleVi1,scaleVd1,10,'filled'); hold all
    plot(bin,meanscaleVd1,'ro-','linewidth',2);
    p=polyfit(scaleVi1,scaleVd1,1);
    idsla = p(1);
    xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);
    xlabel('Rescaled l_i');
    ylabel('Rescaled l_d');
    title(sprintf('slope = %.2f',idsla));
    set(gca,'fontsize',14);

    subplot(1,3,2)
    disc = discretize(scaleVi1,bin);
    for i = 1:length(bin)
        meanscaledVid(i) = nanmean(scaledVid(disc==i));
    end
    scatter(scaleVi1,scaledVid,10,'filled'); hold all
    plot(bin,meanscaledVid,'ro-','linewidth',2);
    p=polyfit(scaleVi1,scaledVid,1);
    idslb = p(1);
    xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);
    xlabel('Rescaled l_i');
    ylabel('Rescaled \delta_i_d');
    title(sprintf('slope = %.2f',idslb));
    set(gca,'fontsize',14);

    subplot(1,3,3)
    disc = discretize(scaleVi1,bin);
    for i = 1:length(bin)
        meanscaledTid(i) = nanmean(scaledTid(disc==i));
    end
    scatter(scaleVi1,scaledTid,10,'filled'); hold all
    plot(bin,meanscaledTid,'ro-','linewidth',2);
    p=polyfit(scaleVi1,scaledTid,1);
    idslc = p(1);
    xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);
    xlabel('Rescaled l_i');
    ylabel('Rescaled \tau_i_d');
    title(sprintf('slope = %.2f',idslc));
    set(gca,'fontsize',14);
    drawnow;

    %% plot auto correlations
    bin = 0.5:0.1:1.5;
    h4=figure('position',[900 400 900 300]);
    subplot(1,3,1)
    dVi1 = scaledVi(1:end-1);
    dVi2 = scaledVi(2:end);
    disc = discretize(dVi1,bin);
    for i = 1:length(bin)
        meandVi2(i) = nanmean(dVi2(disc==i));
    end
    scatter(dVi1,dVi2,10,'filled'); hold all
    plot(bin,meandVi2,'ro-','linewidth',2);
    p=polyfit(dVi1,dVi2,1);
    slautoii = p(1);
    title(sprintf('slope = %.2f',slautoii));
    xlabel('Rrescaled \delta_i_i^n^-^1');
    ylabel('Rrescaled \delta_i_i^n');
    set(gca,'fontsize',14);
    xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);

    subplot(1,3,2)
    dVid1 = scaledVid(1:end-1);
    dVid2 = scaledVid(2:end);
    disc = discretize(dVid1,bin);
    for i = 1:length(bin)
        meandVid2(i) = nanmean(dVid2(disc==i));
    end
    scatter(dVid1,dVid2,10,'filled'); hold all
    plot(bin,meandVid2,'ro-','linewidth',2);
    p=polyfit(dVid1,dVid2,1);
    slautoid = p(1);
    title(sprintf('slope = %.2f',slautoid));
    xlabel('Rrescaled \delta_i_d^n^-^1');
    ylabel('Rrescaled \delta_i_d^n');
    set(gca,'fontsize',14);
    xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);

    subplot(1,3,3)
    dVdd1 = scaledVd(1:end-1);
    dVdd2 = scaledVd(2:end);
    disc = discretize(dVdd1,bin);
    for i = 1:length(bin)
        meandVdd2(i) = nanmean(dVdd2(disc==i));
    end
    scatter(dVdd1,dVdd2,10,'filled'); hold all
    plot(bin,meandVdd2,'ro-','linewidth',2);
    p=polyfit(dVdd1,dVdd2,1);
    slautodd = p(1);
    title(sprintf('slope = %.2f',slautodd));
    xlabel('Rrescaled \delta_d_d^n^-^1');
    ylabel('Rrescaled \delta_d_d^n');
    set(gca,'fontsize',14);
    xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);
    drawnow;


    k_delta_dd(ir) = slautodd;

    %% save data

    ii_sizer_slope(ir) = isla;
    dd_sizer_slope(ir) = dsla;
    id_sizer_slope(ir) = idsla;
    ii_adder_slope(ir) = islb;
    dd_adder_slope(ir) = dslb;
    id_adder_slope(ir) = idslb;
    ii_timer_slope(ir) = islc;
    dd_timer_slope(ir) = dslc;
    id_timer_slope(ir) = idslc;
    k_delta_ii(ir) = slautoii;
    k_delta_dd(ir) = slautodd;



    lambdaCD(ir) = idsla;
    lambdaI(ir) = isla;
    lambdaH(ir) = dsla;

    meanVi(ir) = nanmean(rep_Vi1);
    meandVi(ir) = nanmean(rep_dVi);
    meandTi(ir) = nanmean(rep_dTi);
    meanVd(ir) = nanmean(rep_Vd1);
    meandVd(ir) = nanmean(rep_dVd);
    meandTd(ir) = nanmean(rep_dTd);
    meandVid(ir) = nanmean(rep_dVid);
    meandTid(ir) = nanmean(rep_dTid);

    disp([num2str(ir) '/' num2str(length(DTlist))]);
end
%%
figure;
subplot(1,2,1)
plot(alist,smooth(meanVi,5),'linewidth',2); hold all
xlim([0 20]);
ylim([0 10]);
xlabel('DnaA expression level, \alpha');
ylabel('initiation mass');
axis square

subplot(1,2,2)
plot(alist,smooth(ii_adder_slope,5),'linewidth',2); hold all
plot(alist,zeros(size(alist)),'--');
xlim([0 20]);
ylim([-1 1]);
xlabel('DnaA expression level, \alpha');
ylabel('II adder slope');
axis square