% Script for simulating the individual based model(IBM) and numerically solving the mean-field model in the article 
% Surendran et al. (2020), Small-scale spatial structure affects predator-prey dynamics and coexistence. 
%
%
% Author: Anudeep Surendran
%         anudeep.surendran@hdr.qut.edu.au
%         School of Mathematical Sciences
%         Queensland University of Technology
%         Brisbane, Australia
%
% Last update: 20 March 2020
%
% Instructions:
% 1- This script outputs four Matlab figure files ('.fig'):
%    a) snapshot of the IBM showing the locations of consumers and resources at the initial time (t=0)
%    b) snapshot of the IBM showing the locations of consumers and resources at the final time (t=200)
%    c) plot showing the density of consumers and resources as a function of time
%    d) Plot showing the pair correlation functions (PCFs) computed at t=200 as a function of separation distance.
% 2- Change the model parameters below to simulate various cases considered in the article and any other new senarios.

%.............................Model parameters.............................
% 
% Initial population sizes

Initial_size=300; % Total population size
Initial_size_c=100;Initial_size_r=Initial_size-Initial_size_c; % Population sizes of consumers(c) and resources(r)

% Constant neighbour independent rates

m_c=0.1;m_r=0.1; % Intrinsic movement rate of consumers(c) and resources(r)
p_c=0;p_r=0.2;   % Intrinsic proliferation rate of consumers(c) and resources(r)
d_c=0.1;d_r=0;   % Intrinsic death rate of consumers(c) and resources(r)

% Movement and dispersal distance parameters

mus_c=0.4;sigmas_c=0.1;sigmad_c=4.0; % Mean and standard deviation movement distances and dispersal range for consumers 					
mus_r=0.4;sigmas_r=0.1;sigmad_r=4.0; % Mean and standard deviation movement distances and dispersal range for resources

% Competition and predation strengths

gammac_cc=0.001;gammac_rr=0.001;gammap_cr=0.003;gammap_rc=0.004;

% Competition and predation ranges

sigmac_cc=4.0;sigmac_rr=4.0;sigmap_cr=4.0;sigmap_rc=4.0;

L=20; % Domain length

% Parameters and variables introduced for simulation of the IBM and calculate PCFs....
N_realizations=1; % Number of realizations
t_final=200; % Final time
xi_final=8.0;xi_initial=0.1;xi_increment=0.04;dxi=0.2; % Parameters for binning the distance($|\xi|$) to calculate PCFs 
steps=5000000; % Number of steps used in the Gillespie algorithm

%.................Preallocating arrays.....................................

t_ar=0:.1:t_final; % Time points at which density of consumers and resources are recorded for averaging.
xi_ar=xi_initial:xi_increment:xi_final;% For binning the distance
Z_c_ar_final=zeros(size(t_ar));Z_r_ar_final=zeros(size(t_ar)); % Array to store density of consumers(c) and resources(c) 
Z_c_ar=zeros(size(t_ar));Z_r_ar=zeros(size(t_ar));% Array to store density of consumers(c) and resources(c) (in each realizations)
C_cc_final=zeros(size(xi_ar));C_cr_final=zeros(size(xi_ar));C_rr_final=zeros(size(xi_ar));% Array to store auto and cross PCFs values.

%.....................The actual code begins from here.....................

% Loop for creating ensamble of realizations
for rlzn=1:N_realizations
           
    % Initializing the population and calculating initial rates 
    t=0; 
    num=Initial_size;num_c=Initial_size_c;num_r=Initial_size_r;
    x_initial=(rand(1,num)*L)+(-L/2.0);y_initial=(rand(1,num)*L)+(-L/2.0); % Individuals are randomely placed on LxL domain
    x=x_initial;y=y_initial;% x and y coordinates of the locations of individuals.
    % Here we consider first num_c individuals as consumers and remaining num_r=num-num_c as resources
    
    % Pre allocating arrays
    distance_ar=zeros(num); % Distance matrix storing the separation distance of all individuals from each other
    D_ar=zeros(1,num);P_ar=zeros(1,num);M_ar=zeros(1,num); % Array for storing the Movement, proliferation and death rates of each individuals    
    D_ar(1:num_c)=d_c*ones(1,num_c);D_ar(num_c+1:num)=d_r*ones(1,num_r);
    P_ar(1:num_c)=p_c*ones(1,num_c);P_ar(num_c+1:num)=p_r*ones(1,num_r);
    M_ar(1:num_c)=m_c*ones(1,num_c);M_ar(num_c+1:num)=m_r*ones(1,num_r);
    count=1;
    Z_c_ar(count)=num_c/(L^2);Z_r_ar(count)=num_r/(L^2); % Computing initial density of consumers and resources
    
    for i=1:num % In these loops, we calculate the neighbour dependent event rates  as a result of the interaction of individuals
        for j=1:num
            xi_x=periodic(x(j)-x(i),L);xi_y=periodic(y(j)-y(i),L);
            xi_square=xi_x^2+xi_y^2;  
            distance_ar(i,j)=sqrt(xi_square); % Distance between individuals
            
            if (j ~= i)  % Neglecting the self pairs                
                if (i<=num_c) && (j<=num_c) % Interaction between consumers
                    D_ar(i)=D_ar(i)+w(xi_square,gammac_cc,sigmac_cc);
                elseif (i>num_c) && (j>num_c)% Interaction between resources
                    D_ar(i)=D_ar(i)+w(xi_square,gammac_rr,sigmac_rr);
                elseif (i<=num_c) && (j>num_c)% Interaction of a consumer with a resource
                    P_ar(i)=P_ar(i)+w(xi_square,gammap_cr,sigmap_cr);
                else% Interaction of a resource with a consumer
                    D_ar(i)=D_ar(i)+w(xi_square,gammap_rc,sigmap_rc);
                end
            end
        end
    end
%................Gillespie simulation starts from here.....................

    rn=rand(1,steps); % uniform random numbers to use in Gillespie algorithm
    for ns=1:steps    
        % Calculating net event rates and time increment        
		rate_P = sum(P_ar);              % Total proliferation rate
		rate_D = sum(D_ar);              % Total death rate
        rate_M = sum(M_ar);              % Total movement rate
		net_rate = rate_P+rate_D+rate_M; % Net events rate		
		dt=exprnd(1/net_rate);  % Time step for the event (exponentially distributed) 
        
        if (rn(ns) <= (rate_M/net_rate)) % If the event is movement
            cum_sum=cumsum(M_ar);% Cumulative sum of movement rate
            R=rand*rate_M; % Uniform random number on interval [0,rate_M]
            choose=find(cum_sum>R, 1, 'first'); % find first indice of cum_sum that is greater than R. This individual is selected for next event
            old_loc_x=x(choose);old_loc_y=y(choose); % Previous location of the individual chosen for movement event
            if_consumer=(choose <= num_c);if_resource=(choose > num_c);% Checking whether chosen individual is conumer or resource
            mus=(if_consumer*mus_c)+(if_resource*mus_r);% Assigning the movement distance parameters according to which subpopulation the individual belongs
            sigmas=(if_consumer*sigmas_c)+(if_resource*sigmas_r); %            ""
            magnitude_xi=-1; % Computing the movement distance
            while (magnitude_xi < mus-(4.0*sigmas) || magnitude_xi > (mus+(4.0*sigmas))) % Truncation 
                magnitude_xi=normrnd(mus,sigmas);
            end
            angle_xi=2.0*pi*rand; % Uniformly distributed direction angle
            new_loc_x=x(choose)+(magnitude_xi*cos(angle_xi)); % Updating the new location after the movement
            new_loc_y=y(choose)+(magnitude_xi*sin(angle_xi)); %                      ""
            x(choose)=periodic(new_loc_x,L);y(choose)=periodic(new_loc_y,L);%        ""            
            M_ar(choose)=(if_consumer*m_c)+(if_resource*m_r); % Updating the event rates after the movement event
            P_ar(choose)=(if_consumer*p_c)+(if_resource*p_r); %                      ""
            D_ar(choose)=(if_consumer*d_c)+(if_resource*d_r); %                      ""
            for i=1:num % Updating the event rates after the movement event           
                xi_x=periodic(x(i)-x(choose),L);xi_y=periodic(y(i)-y(choose),L);
                xi_x_old=periodic(x(i)-old_loc_x,L);xi_y_old=periodic(y(i)-old_loc_y,L);
                xi_square=xi_x^2+xi_y^2;xi_square_old=xi_x_old^2+xi_y_old^2;
                distance_ar(choose,i)=sqrt(xi_square); % Distance between individuals
                distance_ar(i,choose)=sqrt(xi_square);                
                if (i ~= choose)
                    if (if_consumer && (i <= num_c))
                        D_ar(choose)=D_ar(choose)+w(xi_square,gammac_cc,sigmac_cc);
                        D_ar(i)=D_ar(i)+w(xi_square,gammac_cc,sigmac_cc)-w(xi_square_old,gammac_cc,sigmac_cc);
                    elseif (if_resource && (i > num_c))
                        D_ar(choose)=D_ar(choose)+w(xi_square,gammac_rr,sigmac_rr);
                        D_ar(i)=D_ar(i)+w(xi_square,gammac_rr,sigmac_rr)-w(xi_square_old,gammac_rr,sigmac_rr);
                    elseif (if_consumer && (i > num_c))
                        P_ar(choose)=P_ar(choose)+w(xi_square,gammap_cr,sigmap_cr);
                        D_ar(i)=D_ar(i)+w(xi_square,gammap_rc,sigmap_rc)-w(xi_square_old,gammap_rc,sigmap_rc);
                    else
                        D_ar(choose)=D_ar(choose)+w(xi_square,gammap_rc,sigmap_rc);
                        P_ar(i)=P_ar(i)+w(xi_square,gammap_cr,sigmap_cr)-w(xi_square_old,gammap_cr,sigmap_cr);
                    end
                end
            end
        elseif (rn(ns) > (rate_M/net_rate) && rn(ns) <= ((rate_M+rate_P)/net_rate)) % If the event is proliferation
            cum_sum=cumsum(P_ar);% Cumulative sum of proliferation rates
            R=rand*rate_P; % Uniform random number on interval [0,rate_P]
            choose=find(cum_sum>R, 1, 'first'); % find first indice of cum_sum that is greater than R. This individual is selected for next event
            if_consumer=(choose <= num_c);if_resource=(choose > num_c);% Checking whether chosen individual is conumer or resource
            sigmad=(if_consumer*sigmad_c)+(if_resource*sigmad_r);% Assiging the dispersal range according to which subpopulation the individual belong
            dispersal=mvnrnd(0,sigmad,2);% New location of individual from bivariate normal distribution
            new_loc_x=periodic(x(choose)+dispersal(1),L);new_loc_y=periodic(y(choose)+dispersal(2),L);
            x=[x(1:choose) new_loc_x x(choose+1:end)];y=[y(1:choose) new_loc_y y(choose+1:end)];            
            M_ar=[M_ar(1:choose) (if_consumer*m_c)+(if_resource*m_r) M_ar(choose+1:end)];
            P_ar=[P_ar(1:choose) (if_consumer*p_c)+(if_resource*p_r) P_ar(choose+1:end)];
            D_ar=[D_ar(1:choose) (if_consumer*d_c)+(if_resource*d_r) D_ar(choose+1:end)];
            row=zeros(1,num+1);coloumn=zeros(num+1,1);
            for i=1:num+1
                xi_x=periodic(x(i)-x(choose+1),L);xi_y=periodic(y(i)-y(choose+1),L);
                xi_square=xi_x^2+xi_y^2;
                row(i)=sqrt(xi_square);coloumn(i)=sqrt(xi_square);% To record distance from daughter individual                
                if (i ~= choose+1)
                    if (if_consumer && (i <= num_c))
                        D_ar(choose+1)=D_ar(choose+1)+w(xi_square,gammac_cc,sigmac_cc);
                        D_ar(i)=D_ar(i)+w(xi_square,gammac_cc,sigmac_cc);
                    elseif (if_resource && (i > num_c))
                        D_ar(choose+1)=D_ar(choose+1)+w(xi_square,gammac_rr,sigmac_rr);
                        D_ar(i)=D_ar(i)+w(xi_square,gammac_rr,sigmac_rr);
                    elseif (if_consumer && (i > num_c))
                        P_ar(choose+1)=P_ar(choose+1)+w(xi_square,gammap_cr,sigmap_cr);
                        D_ar(i)=D_ar(i)+w(xi_square,gammap_rc,sigmap_rc);
                    else
                        D_ar(choose+1)=D_ar(choose+1)+w(xi_square,gammap_rc,sigmap_rc);
                        P_ar(i)=P_ar(i)+w(xi_square,gammap_cr,sigmap_cr);
                    end
                end
            end
            dummy=zeros(num+1);% Tempor
            dummy(1:choose,1:choose)=distance_ar(1:choose,1:choose);
            dummy(1:choose,choose+2:end)=distance_ar(1:choose,choose+1:end);
            dummy(choose+2:end,1:choose)=distance_ar(choose+1:end,1:choose);
            dummy(choose+2:end,choose+2:end)=distance_ar(choose+1:end,choose+1:end);
            dummy(1:end,choose+1)=coloumn;dummy(choose+1,1:end)=row;
            distance_ar=dummy;
            num=num+1;num_c=num_c+if_consumer;num_r=num_r+if_resource; % Updating the number of individuals after proliferation event
        else
            cum_sum=cumsum(D_ar);% Cumulative sum of death rates
            R=rand*rate_D; % Uniform random number on interval [0,rate_D]
            choose=find(cum_sum>R, 1, 'first'); % find first indice of cum_sum that is greater than R. This individual is selected for next event
            if_consumer=(choose <= num_c);if_resource=(choose > num_c);% Checking whether chosen individual is conumer or resource            
            for i=1:num
                xi_x=periodic(x(choose)-x(i),L);xi_y=periodic(y(choose)-y(i),L);
                xi_square=xi_x^2+xi_y^2;
                if (i ~= choose)
                    if (if_consumer && (i <= num_c))
                        D_ar(i)=D_ar(i)-w(xi_square,gammac_cc,sigmac_cc);
                    elseif (if_resource && (i > num_c))
                        D_ar(i)=D_ar(i)-w(xi_square,gammac_rr,sigmac_rr);
                    elseif (if_consumer && (i > num_c))
                        D_ar(i)=D_ar(i)-w(xi_square,gammap_rc,sigmap_rc);
                    else
                        P_ar(i)=P_ar(i)-w(xi_square,gammap_cr,sigmap_cr);
                    end
                end
            end
            x=[x(1:choose-1) x(choose+1:end)];y=[y(1:choose-1) y(choose+1:end)];M_ar=[M_ar(1:choose-1) M_ar(choose+1:end)];
            P_ar=[P_ar(1:choose-1) P_ar(choose+1:end)];D_ar=[D_ar(1:choose-1) D_ar(choose+1:end)];
            distance_ar=[distance_ar(1:choose-1,1:choose-1) distance_ar(1:choose-1,choose+1:end);distance_ar(choose+1:end,1:choose-1) distance_ar(choose+1:end,choose+1:end)];
            num=num-1;num_c=num_c-if_consumer;num_r=num_r-if_resource; % Updating the number of individuals after death event
        end
        t=t+dt; % Updating the time
        
        if (t >= t_ar(count))
            Z_c_ar(count)=num_c/(L^2);Z_r_ar(count)=num_r/(L^2); % Storing the desnsity of consumers and resource at time t
            if (t >= t_final)% Computing the PCFs, C_{ij}(|\xi|, t), at t=t_final
                C_cc=zeros(size(xi_ar));C_cr=zeros(size(xi_ar));C_rr=zeros(size(xi_ar));% Pre-allocation
                for ab=1:length(xi_ar)
                    xi=xi_ar(ab);
                    count_Ccc=0;count_Ccr=0;count_Crr=0;
                    for i=1:num
                        for j=1:num
                            if (i ~= j)% Counting distances excluding the self distance (when i=j).
                                xival=distance_ar(j,i);
                                if (xival >= (xi-(dxi/2.0)) && xival <= (xi+(dxi/2.0))) % Binning the distances
                                    count_Ccc=count_Ccc+(1.0*(i <= num_c && j <= num_c));
                                    count_Crr=count_Crr+(1.0*(i > num_c && j > num_c));
                                    count_Ccr=count_Ccr+(1.0*(i <= num_c && j > num_c));
                                end
                            end
                        end
                    end
                    C_cc(ab)=count_Ccc/((num_c*num_c)*(2.0*pi*xi*dxi)/(L^2));% Normalizing the count of distance to plot PCFs
                    C_cr(ab)=count_Ccr/((num_c*num_r)*(2.0*pi*xi*dxi)/(L^2));%               ""
                    C_rr(ab)=count_Crr/((num_r*num_r)*(2.0*pi*xi*dxi)/(L^2));%               ""
                end
                C_cc_final=C_cc_final+C_cc;C_cr_final=C_cr_final+C_cr;C_rr_final=C_rr_final+C_rr;
                break
            end
            count=count+1;
        end
    end
    Z_c_ar_final=Z_c_ar_final+Z_c_ar;Z_r_ar_final=Z_r_ar_final+Z_r_ar;
end
% Averaged results from IBM
Z_c_ar_final=Z_c_ar_final/N_realizations;Z_r_ar_final=Z_r_ar_final/N_realizations;
C_cc_final=C_cc_final/N_realizations;C_cr_final=C_cr_final/N_realizations;C_rr_final=C_rr_final/N_realizations;

%...........................IBM simulation ends here.............................

% Numerical solution of the Mean-field model using ode45 routine

a=2.0*pi*gammap_cr*(sigmap_cr^2);
b=2.0*pi*gammac_cc*(sigmac_cc^2);
c=2.0*pi*gammap_rc*(sigmap_rc^2);
e=2.0*pi*gammac_rr*(sigmac_rr^2);
f = @(t,x) [x(1)*(a*x(2)-d_c-b*x(1));x(2)*(p_r-c*x(1)-e*x(2))];
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t_MF,Z_MF]=ode45(f,[0 t_final],[Initial_size_c/(L^2) Initial_size_r/(L^2)],opts);

%.................Generating plots.........................................

figure(1)
% Snapshot of the IBM showing the locations of consumers and resources at the initial time (t=0)
plot(x_initial(1:Initial_size_c),y_initial(1:Initial_size_c),'.r','markersize',16);hold on
plot(x_initial(Initial_size_c+1:end),y_initial(Initial_size_c+1:end),'.k','markersize',16);
xlabel('x','FontSize',15,'FontWeight','bold')
ylabel('y','FontSize',15,'FontWeight','bold')
set(gca,'FontSize',14)
savefig('Initial_population.fig')

figure(2)
% Snapshot of the IBM showing the locations of consumers and resources at the final time (t=200)
plot(x(1:num_c),y(1:num_c),'.r','markersize',16);hold on
plot(x(num_c+1:end),y(num_c+1:end),'.k','markersize',16);
xlabel('x','FontSize',15,'FontWeight','bold')
ylabel('y','FontSize',15,'FontWeight','bold')
set(gca,'FontSize',14)
savefig('Final_population.fig')

figure(3)
% plot showing the density of consumers and resources as a function of time
p1=plot(t_ar,Z_c_ar_final,'-r',t_MF,Z_MF(:,1),'--r');hold on
p2=plot(t_ar,Z_r_ar_final,'-k',t_MF,Z_MF(:,2),'--k');
legend({'Z_{c}(t)','Z_{c}(t)','Z_{r}(t)','Z_{r}(t)'},'FontSize',15);
axis([0 200 0 1.5])
xlabel('t','FontSize',15,'FontWeight','bold')
ylabel('Z_{i}(t)','FontSize',15,'FontWeight','bold')
set(p1,'LineWidth',2);set(p2,'LineWidth',2);
set(gca,'FontSize',14)
hold off
savefig('Z_Vs_t.fig')

figure(4)
% Plot showing the pair correlation functions (PCFs) computed at t=200 as a function of separation distance.
p3=plot(xi_ar,C_cc_final,'-r',xi_ar,C_rr_final,'-k',xi_ar,C_cr_final,'-b');
axis([0 8 0 2])
legend({'C_{cc}(|\xi|, t)','C_{rr}(|\xi|, t)','C_{cr}(|\xi|, t)'},'FontSize',15)
set(p3,'LineWidth',2);
xlabel('|\xi|','FontSize',15,'FontWeight','bold')
ylabel('C_{ij}(|\xi|, t)','FontSize',15,'FontWeight','bold')
set(gca,'FontSize',14)
hold off
savefig('C_Vs_xi.fig')
%.......................Main code ends here................................

 function a = periodic(x,L)
% Function to impliment periodic boundary condition.
% Input   :  x  =  x or y cordinate of an individual's location
%		     L	=  spatial domain length
%
% Output  :  impliment periodic boundary condition 
    if x<(-L/2.0)	
        a=x+L;
    elseif x>=(L/2.0)	
		a=x-L;
    else	
        a=x;
    end
end

function a = w(dsquare,gamma,sigma) 
% Function to calculate interaction kernel
% Input   :  dsquare =  distance square
%		     gamma   =  interaction strength
%            sigma   =  spatial extent of interaction
%
% Output  :  value of interaction kernel
a=gamma*exp(-(dsquare/(2.0*(sigma^2))));
end                