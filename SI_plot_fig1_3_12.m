clear all
close all

load stochastic_titration_correlation_model.mat

%% figure 1
h1 = figure;
plot(lambdalist,meanVi*10,'o','linewidth',2); hold all
plot(lambdalist,mi*10,'-','linewidth',2); hold all
xlim([0 2.1]);
ylim([0 10]);
xlabel('growth rate (hr^-^1)');
ylabel('Initiation Volume (10^?^1^0 OD_6_0_0 ml)');
legend('V_i simulation','V_i exp fitted','location','best');
set(gca,'fontsize',14);


%% figure S12

h1 = figure;
plot(lambdalist,meanVi,'o','linewidth',2); hold all
plot(lambdalist,meanVd,'o','linewidth',2); hold all
plot(lambdalist,meandVid,'o','linewidth',2); hold all
plot(lambdalist,mi,'-','linewidth',2,'color',[0, 114, 189]/255); hold all
plot(lambdalist,md,'-','linewidth',2,'color',[217, 83, 25]/255); hold all
plot(lambdalist,dVidlist,'-','linewidth',2,'color',[237, 177, 32]/255); hold all
xlim([0 2.1]);
xlabel('growth rate (hr^-^1)');
ylabel('Mean Volume (10^-^9 OD_6_0_0 ml)');
legend('V_i simulation','V_d simulation','\Delta_i_d simulation','V_i exp fitted','V_d exp fitted','\Delta_i_d exp fitted','location','best');
set(gca,'fontsize',14);


%% figure 3 



h5 = figure('position',[100 100 1200 400]);
subplot(1,3,1)
plot(lambdalist, smooth(ii_sizer_slope,5),'o-','linewidth',2);  hold all
plot([0 2],[0 0],'r--');
ylim([-1 1]);xlim([0 2.1]);
ylabel('Sizer slope');
xlabel('growth rate (hr^-^1)');
% legend('k_i_i','k_d_d','k_i_d','location','best');
set(gca,'fontsize',14);

subplot(1,3,2)
plot(lambdalist, smooth(ii_adder_slope,5),'o-','linewidth',2);  hold all
plot([0 2],[0 0],'r--');
ylim([-1 1]);xlim([0 2.1]);
ylabel('Adder slope');
xlabel('growth rate (hr^-^1)');
% legend('k_i_i','k_d_d','k_i_d','location','best');
set(gca,'fontsize',14);


subplot(1,3,3)
plot(lambdalist, smooth(ii_timer_slope,5),'o-','linewidth',2);  hold all
plot([0 2],[0 0],'r--');
ylim([-1.5 1.5]);xlim([0 2.1]);
ylabel('Timer slope');
xlabel('growth rate (hr^-^1)');
% legend('k_i_i','k_d_d','k_i_d','location','best');
set(gca,'fontsize',14);


%% figure 5
