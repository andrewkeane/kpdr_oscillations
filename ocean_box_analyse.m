%% Set up rhs function
funcs=set_funcs(...
    'sys_rhs',@(x,p)ocean_box_rhs(x,p),...
    'sys_tau',@()[4],...
    'sys_ntau',@()1,...
    'x_vectorized',true);

%% Create initial steady-state solution with sigma=0
par=[-0.3 1 0 900 0];
x0=[34.2458 34.6701]'; % Initial solution (can be found by numerical simulation)
stst.kind='stst';
stst.parameter=par;
stst.x=x0;

%%
flag_newhheur=2;
method=df_mthod(funcs,'stst',flag_newhheur);
method.point.print_residual_info=1;
method.point.minimal_accuracy=1e-10;

%% Correct steady state
[stst,success]=p_correc(funcs,stst,[],[],method.point)

%% Perturb in parameter F1 and correct for 2nd steady-state solution
stst2=stst;
stst2.parameter(1)=stst2.parameter(1)+0.001;
[stst2,success]=p_correc(funcs,stst2,[],[],method.point)

%% Initialise branch of steady states in F1
branch_stst_sigma0=df_brnch(funcs,[1],'stst');
branch_stst_sigma0.parameter.max_step(1,:)=[0 0.001];
branch_stst_sigma0.point=stst;
branch_stst_sigma0.point(2)=stst2;
branch_stst_sigma0.method=method;
branch_stst_sigma0.method.continuation.plot=1;

%% Continue in one direction
figure(1); 
[branch_stst_sigma0,s,f,r]=br_contn(funcs,branch_stst_sigma0,100)

%% Calculate stability
branch_stst_sigma0=br_stabl(funcs,branch_stst_sigma0,0,0);

%% Plot real part of eigenvalues
[xm,ym]=df_measr(1,branch_stst_sigma0);
figure(2); 
br_plot(branch_stst_sigma0,[],ym);
xlabel('point number') ;
ylabel('Re($\lambda$)','Interpreter','latex');
% Hopf bifurcation occurs near point 79

%% Branch off Hopf point
[branch_per_sigma0,suc]=SetupPsol(funcs,branch_stst_sigma0,79,'degree',4,'intervals',200,'contpar',1,'radius',1e-3,'minimal_accuracy',1e-6,'print_residual_info',1,'newton_max_iterations',10,'newton_nmon_iterations',4);
branch_per_sigma0.method.point.print_residual_info=1;
branch_per_sigma0.parameter.max_step=[0 1];
figure(3);hold on;
[branch_per_sigma0,s,f,r]=br_contn(funcs,branch_per_sigma0,50)

%% Calculate stability
branch_per_sigma0=br_stabl(funcs,branch_per_sigma0,0,0);

%% Plot absolute value of Floquet multipliers
[xm,ym]=df_measr(1,branch_per_sigma0);
figure(4); 
br_plot(branch_per_sigma0,[],ym);
xlabel('point number') ;
ylabel('$| \Lambda |$','Interpreter','latex');

%% Plot branches together
figure(5);hold on;
[xm,ym]=df_measr(0,branch_stst_sigma0);
br_plot(branch_stst_sigma0,xm,ym);
[xm,ym]=df_measr(0,branch_per_sigma0);
ym.col='max';
br_plot(branch_per_sigma0,xm,ym);
axis([-0.258 -0.242 34.21 34.32]);
xlabel('$F_1$','Interpreter','latex') ;
ylabel('max($S_1$)','Interpreter','latex');

%% Create initial steady-state solution with sigma=11
par=[-0.3 1 11 900 0];
x0=[34.2458 34.6701]'; % Initial solution (can be found by numerical simulation)
stst.kind='stst';
stst.parameter=par;
stst.x=x0;

%% Correct steady state
[stst,success]=p_correc(funcs,stst,[],[],method.point)

%% Perturb in parameter F1 and correct for 2nd steady-state solution
stst2=stst;
stst2.parameter(1)=stst2.parameter(1)+0.001;
[stst2,success]=p_correc(funcs,stst2,[],[],method.point)

%% Initialise branch of steady states in F1
branch_stst_sigma11=df_brnch(funcs,[1],'stst');
branch_stst_sigma11.parameter.max_step(1,:)=[0 0.001];
branch_stst_sigma11.point=stst;
branch_stst_sigma11.point(2)=stst2;
branch_stst_sigma11.method=method;
branch_stst_sigma11.method.continuation.plot=1;

%% Continue in one direction
figure(1);clf;
[branch_stst_sigma11,s,f,r]=br_contn(funcs,branch_stst_sigma11,200)

%% Calculate stability
branch_stst_sigma11=br_stabl(funcs,branch_stst_sigma11,0,0);

%% Plot real part of eigenvalues
[xm,ym]=df_measr(1,branch_stst_sigma11);
figure(2); clf; 
br_plot(branch_stst_sigma11,[],ym);
xlabel('point number') ;
ylabel('Re($\lambda$)','Interpreter','latex');
% Hopf bifurcation occurs near point 155

%% Branch off Hopf point
[branch_per_sigma11,suc]=SetupPsol(funcs,branch_stst_sigma11,155,'degree',4,'intervals',200,'contpar',1,'radius',1e-3,'minimal_accuracy',1e-6,'print_residual_info',1,'newton_max_iterations',10,'newton_nmon_iterations',4);
branch_per_sigma11.method.point.print_residual_info=1;
branch_per_sigma11.parameter.max_step=[0 5];
figure(3); clf;
[branch_per_sigma11,s,f,r]=br_contn(funcs,branch_per_sigma11,100)

%% Calculate stability
branch_per_sigma11=br_stabl(funcs,branch_per_sigma11,0,0);

%% Plot absolute value of Floquet multipliers
[xm,ym]=df_measr(1,branch_per_sigma11);
figure(4); clf; 
br_plot(branch_per_sigma11,[],ym);
xlabel('point number') ;
ylabel('$| \Lambda |$','Interpreter','latex');

%% Plot branches together
figure(6);hold on;
[xm,ym]=df_measr(0,branch_stst_sigma11);
br_plot(branch_stst_sigma11,xm,ym);
[xm,ym]=df_measr(0,branch_per_sigma11);
ym.col='max';
br_plot(branch_per_sigma11,xm,ym);
axis([-0.213 -0.205 34.19 34.285]);
xlabel('$F_1$','Interpreter','latex') ;
ylabel('max($S_1$)','Interpreter','latex');

%% Branch of Hopf bifurcations for tau=900
[hbranch_tau900,suc]=SetupHopf(funcs,branch_stst_sigma0,79,'contpar',[1 3],'dir',3,'step',1e-3,'correc',1,'newton_max_iterations',10,'newton_nmon_iterations',4)
hbranch_tau900.parameter.max_step=[0 1e-1];
hbranch_tau900.parameter.min_bound(2,:)=[3 0];

figure(7);hold on;
[hbranch_tau900,s,f,r]=br_contn(funcs,hbranch_tau900,250)

%% Compute first Lyapunov coefficient to determine criticality
[L1_hbr_tau900,L1low_hbr_tau900,hbranch_tau900]=HopfLyapunovCoefficients(funcs,hbranch_tau900);
figure(8);
plot(1:size(L1_hbr_tau900,2),L1_hbr_tau900);

%% Branch of Hopf bifurcations for tau=1100
branch_stst_sigma0.point(79).parameter(4)=1100;
[hbranch_tau1100,suc]=SetupHopf(funcs,branch_stst_sigma0,79,'contpar',[1 3],'dir',3,'step',1e-3,'correc',1,'newton_max_iterations',10,'newton_nmon_iterations',4)
hbranch_tau1100.parameter.max_step=[0 1e-1];
hbranch_tau1100.parameter.max_step=[1 5e-4];
hbranch_tau1100.parameter.min_bound(2,:)=[3 0];

figure(7);hold on;
[hbranch_tau1100,s,f,r]=br_contn(funcs,hbranch_tau1100,300)

%% Compute first Lyapunov coefficient to determine criticality
[L1_hbr_tau1100,L1low_hbr_tau1100,hbranch_tau1100]=HopfLyapunovCoefficients(funcs,hbranch_tau1100);
figure(8); clf;
plot(1:size(L1_hbr_tau1100,2),L1_hbr_tau1100);