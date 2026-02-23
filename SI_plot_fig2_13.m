clear all
close all

lambdalist = [0.1:0.1:2];
DTlist = round(log(2)./lambdalist * 60);
alist(1:6) = [20 17 15 14 13 12];
alist(7:20) = 12:0.43:18;
alpha_m = 0.28;       beta_m=0.99;      m0 = 1;
mi=m0/log(2)*(alpha_m+beta_m*lambdalist).*exp(-alpha_m-beta_m*lambdalist);
md=m0/log(2)*(0.28+0.99*lambdalist);
dVidlist = md-mi;

exp_lambda = [1.60 1.27 1.49 1.14 0.8 0.66 0.68 0.59 0.25 0.44];
exp_kdeltaid = [0.49 0.47 0.36 0.37 0.34 0.30 0.29 0.20 0.13 0.06];
p_exp = polyfit(exp_lambda,exp_kdeltaid,2);
auto_corr_list = lambdalist.^2*p_exp(1)+lambdalist.^1*p_exp(2)+p_exp(3);
auto_corr_list(auto_corr_list<0)=0;

ir = [3 7 18];
for i = 1:length(ir)
    %% global parameters
    gamma0 = lambdalist(ir(i)) ; % cell mass growth rate (min^-1)
    DT = log(2)/gamma0 ; % doubling rate (hr)
    TT = DT*50 ; % total calculation time (hr)
    dt = 0.05 ; % calculation time step (hr)

    epsilon_A = 0.1;
    epsilon_V = 0.1;
    epsilon_id = 0.1;
    auto_corr = auto_corr_list(ir(i));

    gamma = gamma0*(epsilon_V * randn(1)+1);
    Cperiod = round((log(2)./gamma).*(0.8801*gamma*60+0.2375));
    alpha0 = alist(ir(i))*gamma*60;

    datAn = 250 ; % dnaA binding site on datA
    Chns = 150 ; % number of hns binding to each chromosome

    DNAlength = ceil(Cperiod/dt) ; % DNA length is set to fix V fork = 1 ;
    dVidc0 = dVidlist(ir(i));


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
    % Aboxlist =  NaN(1,ceil(TT/DT)+10) ; % threshold at initiation
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
        Abox = Nd*datAn + sum(DNAcopy(1:round(0.8*Cperiod/dt)))/numel(DNAcopy(1:round(0.8*Cperiod/dt)))*Chns;

        % check new replication initiation condition
        if A>Abox  
            nrep = nrep +1; % replication number +1
            Prep(nrep) = 1 ; % generate new replication position
            rpt(nrep) = t*dt; % replication time +dt
            Vi(nrep) = V; % note initiation mass
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

        Tlist(t) = t*dt;
        Vlist(i,t) = V;
        Numlist(i,t) = Num(nd);
        Alist(i,t) = A;
        Aboxlist(i,t) = Abox;
        Olist(i,t) = DNAcopy(1);

    end

end
%% initiation plots
figure;
ir = [3 7 18];
for i = 1:3
    subplot(4,3,i)
    plot(Tlist,Vlist(i,:)./Olist(i,:),'linewidth',2); 
    xlim([10 15]);
    xlabel('Time (hr)');
    ylabel('Mass/OriC');
    title(sprintf('Growth rate = %.1f hr^-^1',lambdalist(ir(i))));
    axis('square');
    set(gca,'fontsize',14);
end

%% case DT = 139, lambda = 0.3
load lambda=0.3.mat
subplot(4,3,4)
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
xlabel('Rescaled V_i^n');
ylabel('Rescaled V_i^n^+^1');
title(sprintf('sizer slope = %.2f',isla));
set(gca,'fontsize',14);
axis('square');


subplot(4,3,7)
meandVi = nan(size(bin));    stddV = nan(size(bin));
for i = 1:length(bin)
    meandVi(i) = nanmean(scaledVi(disc==i));
    stddVi(i) = nanstd(scaledVi(disc==i));
end
scatter(scaleVi1,scaledVi,10,'filled'); hold all
plot(bin,meandVi,'ro-','linewidth',2);
xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);
pb=polyfit(scaleVi1,scaledVi,1);
islb = pb(1);
xlabel('Rescaled V_i');
ylabel('Rescaled \DeltaV_i_i');
title(sprintf('adder slope = %.2f',islb));
set(gca,'fontsize',14);
axis('square');


subplot(4,3,10)
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
xlabel('Rescaled V_i');
ylabel('Rescaled \tau_i_i');
title(sprintf('timer slope = %.2f',islc));
set(gca,'fontsize',14);
axis('square');
drawnow;



%% case DT = 59, lambda = 0.7
load lambda=0.7.mat
subplot(4,3,5)
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
xlabel('Rescaled V_i^n');
ylabel('Rescaled V_i^n^+^1');
title(sprintf('sizer slope = %.2f',isla));
set(gca,'fontsize',14);
axis('square');



subplot(4,3,8)
meandVi = nan(size(bin));    stddV = nan(size(bin));
for i = 1:length(bin)
    meandVi(i) = nanmean(scaledVi(disc==i));
    stddV(i) = nanstd(scaledVi(disc==i));
end
scatter(scaleVi1,scaledVi,10,'filled'); hold all
plot(bin,meandVi,'ro-','linewidth',2);
xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);
pb=polyfit(scaleVi1,scaledVi,1);
islb = pb(1);
xlabel('Rescaled V_i');
ylabel('Rescaled \DeltaV_i_i');
title(sprintf('adder slope = %.2f',islb));
set(gca,'fontsize',14);
axis('square');


subplot(4,3,11)
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
xlabel('Rescaled V_i');
ylabel('Rescaled \tau_i_i');
title(sprintf('timer slope = %.2f',islc));
set(gca,'fontsize',14);
axis('square');
drawnow;



%% case DT = 23, lambda = 1.8
load lambda=1.8.mat
subplot(4,3,6)
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
xlabel('Rescaled V_i^n');
ylabel('Rescaled V_i^n^+^1');
title(sprintf('sizer slope = %.2f',isla));
set(gca,'fontsize',14);
axis('square');


subplot(4,3,9)
meandVi = nan(size(bin));    stddV = nan(size(bin));
for i = 1:length(bin)
    meandVi(i) = nanmean(scaledVi(disc==i));
    stddV(i) = nanstd(scaledVi(disc==i));
end
scatter(scaleVi1,scaledVi,10,'filled'); hold all
plot(bin,meandVi,'ro-','linewidth',2);
xlim([0 2]);    ylim([0 2]);    set(gca,'xtick',[0 1 2],'ytick',[0 1 2]);
pb=polyfit(scaleVi1,scaledVi,1);
islb = pb(1);
xlabel('Rescaled V_i');
ylabel('Rescaled \DeltaV_i');
title(sprintf('adder slope = %.2f',islb));
set(gca,'fontsize',14);
axis('square');


subplot(4,3,12)
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
xlabel('Rescaled V_i');
ylabel('Rescaled \tau_i');
title(sprintf('timer slope = %.2f',islc));
set(gca,'fontsize',14);
axis('square');
drawnow;



%% division plots
figure;
ir = [3 7 18];
for i = 1:3
    subplot(4,3,i)
    plot(Tlist,Vlist(i,:)./Numlist(i,:),'linewidth',2); 
    xlim([10 15]);
    xlabel('Time (hr)');
    ylabel('Cell Mass');
    title(sprintf('Growth rate = %.1f hr^-^1',lambdalist(ir(i))));
    set(gca,'fontsize',14);
    axis('square');
end

%% case DT = 139, lambda = 0.3
load lambda=0.3.mat
subplot(4,3,4)
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
xlabel('Rescaled V_d^n');
ylabel('Rescaled V_d^n^+^1');
title(sprintf('sizer slope = %.2f',dsla));
set(gca,'fontsize',14);
axis('square');


subplot(4,3,7)
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
xlabel('Rescaled V_d');
ylabel('Rescaled \DeltaV_i_d');
title(sprintf('adder slope = %.2f',dslb));
set(gca,'fontsize',14);
axis('square');


subplot(4,3,10)
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
xlabel('Rescaled V_d');
ylabel('Rescaled \tau_i_d');
title(sprintf('timer slope = %.2f',dslc));
set(gca,'fontsize',14);
axis('square');




%% case DT = 59, lambda = 0.7
load lambda=0.7.mat
subplot(4,3,5)
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
xlabel('Rescaled V_d^n');
ylabel('Rescaled V_d^n^+^1');
title(sprintf('sizer slope = %.2f',dsla));
set(gca,'fontsize',14);
axis('square');

subplot(4,3,8)
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
xlabel('Rescaled V_d');
ylabel('Rescaled \DeltaV_i_d');
title(sprintf('adder slope = %.2f',dslb));
set(gca,'fontsize',14);
axis('square');

subplot(4,3,11)
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
xlabel('Rescaled V_d');
ylabel('Rescaled \tau_i_d');
title(sprintf('timer slope = %.2f',dslc));
set(gca,'fontsize',14);
axis('square');

%% case DT = 23, lambda = 1.8
load lambda=1.8.mat
subplot(4,3,6)
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
xlabel('Rescaled V_d^n');
ylabel('Rescaled V_d^n^+^1');
title(sprintf('sizer slope = %.2f',dsla));
set(gca,'fontsize',14);
axis('square');

subplot(4,3,9)
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
xlabel('Rescaled V_d');
ylabel('Rescaled \DeltaV_i_d');
title(sprintf('adder slope = %.2f',dslb));
set(gca,'fontsize',14);
axis('square');

subplot(4,3,12)
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
ylabel('Rescaled \tau_i_d');
title(sprintf('timer slope = %.2f',dslc));
set(gca,'fontsize',14);
axis('square');

