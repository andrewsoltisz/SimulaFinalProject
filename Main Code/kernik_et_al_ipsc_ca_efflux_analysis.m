function [SERCA_percent, NCX_percent, PMCA_percent] = kernik_et_al_ipsc_ca_efflux_analysis(Time, Iup, INaCa, IpCa, Cai)
% Analysis of cytosolic calcium removal mechanisms (B3 analysis)

%% Constants (copied from ipsc_function)
F = 96.4853415;   % coulomb_per_mmole (in model_parameters)

Cm = 60; %pF
V_tot = 3960; %um^3
Vc_tenT = 16404;
VSR_tenT = 1094;
V_tot_tenT = Vc_tenT+VSR_tenT; %V_total data from Hwang et al., V_c and V_SR  proportionally scaled from Ten Tusscher 2004 values
Vc = V_tot*(Vc_tenT/V_tot_tenT); %=3712.4 um^3 (93.7% total volume)
%V_SR = V_tot*(VSR_tenT/V_tot_tenT);%=247.6 um^3 (6.3% total volume)

%% Find first beat to analyze
inds_time_800=find(Time>800);inds_time_800=inds_time_800(1);
inds_time_1600=find(Time>1600); inds_time_1600=inds_time_1600(1);
[~,inds1]=min(Cai(1:inds_time_800));
[~,inds2]=min(Cai(inds_time_800:inds_time_1600)); inds2=inds_time_800+inds2;

%% Calculate Normalized Ca2+ flux 
% take integral
intJserca = cumtrapz(Time,Iup);
intIncx_ca = cumtrapz(Time,-INaCa*2*Cm/(2.0*Vc*F));
intIpca = cumtrapz(Time,IpCa*Cm/(2.0*Vc*F));

% integral for first beat
fluxJserca = intJserca(inds1:inds2)-intJserca(inds1);
fluxIncx_ca = intIncx_ca(inds1:inds2)-intIncx_ca(inds1);
fluxIpca = intIpca(inds1:inds2)-intIpca(inds1);

% normalize flux
flux_total=fluxJserca+fluxIncx_ca+fluxIpca;
ref=max(flux_total);
fluxJserca_norm=fluxJserca./ref;
fluxIncx_ca_norm=fluxIncx_ca./ref;
fluxIpca_norm=fluxIpca./ref;
Time_flux=Time(inds1:inds2);

% outputs
SERCA_percent=fluxJserca_norm(end).*100;
NCX_percent=fluxIncx_ca_norm(end).*100;
PMCA_percent =fluxIpca_norm(end).*100;

%% plot figure 10A-C
figure, set(gcf,'color','w')
subplot(2,1,1)
plot((Time(inds1:inds2)-Time(inds1)), Cai(inds1:inds2).*1e6, 'Color', [.8 0 .18])
set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('[Ca]i (nM)')
xlabel('Time (ms)');

subplot(2,1,2)
plot((Time_flux-Time(inds1)),fluxJserca_norm,(Time_flux-Time(inds1)),fluxIncx_ca_norm,(Time_flux-Time(inds1)),fluxIpca_norm);
set(gca,'box','off','tickdir','out','fontsize',12)
legend('SERCA', 'NCX', 'non-NCX')
ylabel('Int Ca efflux normalized')
xlabel('Time (ms)')

end