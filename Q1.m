%the yuvals of matlab
%Q1_a.......................................................
w=1;
d=0.1;
v_0=1;
delta=0.01;
N=(w./delta);
eps_0=8.85*10^(-12);
a=10000;
etta = calc_etta(delta,a,eps_0,N,d);

%Q1_a1.......................................................
figure(1)
%we plot etta as function of x
x = 1:N;
y_1 = etta(x,1);
y_2 = etta(x+N,1);
plot((x-N/2)/N,y_1);
%legend('etta(top)')
hold on
plot((x-N/2)/N,y_2);
grid on
xlabel('X[m]');
ylabel('etta[coulomb/m^2]');
title('Q1_a1: etta as a function of x)')
legend('etta(bottom)','etta(top)');

%Q1_a2.......................................................
q_Q1a2 = abs(sum(delta*etta(1:N,1)));    %sum across all of charge units
C_Q1a2 = q_Q1a2/(v_0*eps_0);            %capacitance is that formula bro

fprintf("Q1_a2 : Q = %d\n",q_Q1a2);        
fprintf("Q1_a2 : C = %d\n",C_Q1a2);
%Q1_b1.......................................................
d = 0.01;
delta = w/2;
N = (w/delta);
a = 100000;

C_Q1b1 = get_c(delta,a,eps_0,v_0,N,d);  %capacitance per meter

fprintf("Q1_b1: C =  %d\n",C_Q1b1);

%rest of Q1_b................................................
c_prev = C_Q1b1;
delta = delta/2;
c_cur = get_c(delta,a,eps_0,v_0,w/delta,d);
figure(2)
error = 1;
while((abs(error))>0.01)
    c_prev = c_cur;
    delta = delta/2;
    c_cur = get_c(delta,a,eps_0,v_0,w/delta,d);
    error = (abs(c_cur-c_prev)/(c_prev));
    plot(d/delta,c_cur,'xr');
    hold on
end
title("Q1_b: capacitance as a function of (d/delta)")
xlabel('d/delta[scalar]');
ylabel('C per length [scalar]');
grid on
legend('discrete capacitance values');


%Q1_c.......................................................
w = 1;
v_0 = 1;
N = (w/delta);
eps_0 = 8.85*10^(-12);
a = 10000;

d1 = 1;
d2 = 0.1;
d3 = 0.01;
eta1 = calc_etta(delta,a,eps_0,N,d1);
eta2 = calc_etta(delta,a,eps_0,N,d2);
eta3 = calc_etta(delta,a,eps_0,N,d3);
%distance decreases - capasitance increase - charges on plates increase
figure(3)
x = 1:N;
y_1 = eta1(x,1);
y_2 = eta1(x+N,1);
subplot(3,1,1);
plot((x-N/2)/N,y_1),title("Q1_c: d=1 meter, etta as a function of x");
hold on
grid on
plot((x-N/2)/N,y_2);
xlabel('X[m]');
ylabel('etta[coulomb/m^2]');


x = 1:N;
y_1 = eta2(x,1);
y_2 = eta2(x+N,1);
subplot(3,1,2);
plot((x-N/2)/N,y_1),title("Q1_c: d=0.1 meter, etta as a function of x");
hold on
grid on
plot((x-N/2)/N,y_2);
xlabel('X[m]');
ylabel('etta[coulomb/m^2]');


x = 1:N;
y_1 = eta3(x,1);
y_2 = eta3(x+N,1);
subplot(3,1,3);
plot((x-N/2)/N,y_1),title("Q1_c: d=0.01 meter, etta as a function of x");
hold on
grid on
plot((x-N/2)/N,y_2);
xlabel('X[m]');
ylabel('etta[coulomb/m^2]');


%Q1_d.......................................................
% we build a for loop and calc numeric and analytic
figure(4)
range = 0.01; 
j = 1;
for d = 0.01:range:1
    c_numer(j,1) = get_c(delta,a,eps_0,v_0,N,d);                
    c_anal(j,1) = w/d;                      
    j = j+1; 
end
x = 0.01:range:1;
plot(x,c_numer),title("Q1_d: capasitance as a function of distance");
hold on
plot(x,c_anal);
xlabel('distance[m]');
ylabel('capacitance[F]');
grid on
legend('numeric','analitic');

%now plotting the error graph

err = zeros(length(c_numer));
for i = 1:length(c_numer)
    err(i) = abs(c_anal(i)-c_numer(i))/(0.5*(c_anal(i)+c_numer(i)));
end
figure(5);
plot(x,err, 'b');
title("Q1_d: relative error between numeric and analytic");
legend('error');
xlabel('d[m]');
ylabel('error[unitless]');
grid on
hold off

%Q1_e.......................................................
d=0.1;
etta=calc_etta(delta,a,eps_0,w/delta,d);
R_xy=zeros(2*N,2);

pixel=0.1;
phi_s=zeros(2*3*w/pixel,2*3*w/pixel);     %initialing the potential matrix
d_s=zeros(2*N,1);                         %initialing the distance vector

for i=1:N                              
    %we build an array of x,y coordinates of each tiny pice of electrode
    R_xy(i,2)=(-d/2);
    R_xy(i,1)=(-w/2)+delta*i-0.5*delta;
    R_xy(i+N,2)=(d/2);
    R_xy(i+N,1)=(-w/2)+delta*i-0.5*delta;
   
end
 
%in this loop we pass upon every pixel in our field
%and calculate its potential, so we have a 2D matrix of every potential at
%every point in space
for j=1:(2*3*w/pixel)
    for k=1:(2*3*w/pixel)
        X_=-3*w+j*pixel-0.5*pixel;
        Y_=3*w-k*pixel+0.5*pixel;
        % x_ and y_ are reletive variables helping us calculate the
        % distance from each pixel
        for l=1:(2*N)
            d_s(l,1)=sqrt((Y_-R_xy(l,1)).^2 + (X_-R_xy(l,2)).^2);
        end
        %now calculate for every pixel its potential 
        for l=1:(2*N)
            phi_s(j,k)=phi_s(j,k)+log(d_s(l,1)/a)*(-etta(l,1)*delta)/(2*pi*eps_0);
        end
    end
end

x=-3*w+0.5*pixel:pixel:3*w-0.5*pixel;
space=-3*w+0.5*pixel:pixel:3*w-0.5*pixel;
figure(106);
contour(space,x,phi_s),title("Equal potential lines(x,y)");
xlabel('X[m]');
ylabel('Y[m]');

figure(107);
[Ex,Ey]= gradient(phi_s);
Ex = -Ex;
Ey = -Ey;
quiver(x,space,Ex,Ey),title("Eelctric vector field(x,y)");
xlabel('X[m]');
ylabel('Y[m]');
 
%Q1_f.......................................................
d = 0.1;
w = 1;
etta = calc_etta(delta,a,eps_0,N,d);
P_tot = d*abs(delta*sum(etta(1:N,1)));

fprintf("Q1_f total dipole: %d\n",P_tot);     

space = zeros(2*3*w/pixel,1);                     %the spce on which we calculate the potential
phi_nm = ones(6*w/pixel,1);                       %the angle matrix, consist of 1 and -1
pot_y = zeros(6*w/pixel,1);                       %the Y axis potential space
pot_x = zeros(6*w/pixel,1);                       %the X axis potential space
phi_s_x = zeros(6*w/pixel,1);                     %the x axis potential space (for numeric use)
axis = -3*w+0.5*pixel : pixel : 3*w-0.5*pixel;    %range for plotting the graph

for i=1:(2*3*w/pixel)
    space(i)=-3*w-(0.5-i)*pixel;
    if(i<=(3*w/pixel))
        phi_nm(i,1)=-1;
    end
    %X = 0 analytic potential calculations
    pot_y(i,1)=(P_tot*phi_nm(i,1))/(2*eps_0*pi*abs(space(i,1)));
    %Y = 0 analytic potential calculations
    pot_x(i,1)=(phi_s(3*w/pixel,i)+phi_s(3*w/pixel+1,i))/2;
end

figure(108);
plot(axis,pot_y);
title("potential at x=0");
hold on
plot(axis,phi_s(1:(6*w/pixel),3*w/pixel)),legend("analytic","numeric");
grid on
xlabel('y[m]');
ylabel('V[v]');
hold off

figure(109);
plot(axis,pot_x);
title("potential at y=0");
hold on
plot(axis,phi_s_x),legend("numeric","analytic");
grid on
xlabel('X[m]');
ylabel('V[v]');
hold off


%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
%&&&&&&&&&&&&&&&&&&&&&&&&&&----functions of Question 1 --------&&&&&&&&&&&&&%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%


% the R function is building us the correct division of the plates, and saves it to a vector of tupels 
% these tuples are our coordinates

function R = R_dist(n,m,N,delta,d)    
    if ((n>N && m>N)||(n<=N && m<=N)) %both n,m are on the same plate
        R = abs((m-n)*delta);
    else
        if(n>N)
            n_x=n-N;
        else
            n_x=n;
        end
        if(m>N)
            m_x=m-N;
        else
            m_x=m;
        end
        % n,m are on two different plate
        R = sqrt(((n_x-m_x)*delta).^2   + d.^2);
    end
end

% the etta function is building us the charge density vector

function etta=calc_etta(delta,a,eps_0,N,d)
    z_mat=zeros(2*N,2*N);
    for n=1:(2*N)
        for m=1:(2*N)
            %bulding the Z matrix accourding to the momentum method
            if(n==m)
                z_mat(m,n)=(log(delta./(2*a))-1)*(-0.5)*delta./(pi*eps_0);            
            else
                z_mat(m,n)=((-0.5)*delta./(pi*eps_0))*log(R_dist(n,m,N,delta,d)/a);
            end
        end
    end
    phi=zeros(2*N,1);
    for n=1:N
        phi(n,1)=-0.5;
        phi(n+N,1)=0.5;
    end
    
  %USING LINEAR SOLVING TO find the etta vector
    etta = (z_mat)^(-1)*phi;
end

function c = get_c(delta,a,eps_0,v_0,N,d)
    etta = calc_etta(delta,a,eps_0,N,d);
    c = abs(delta*sum(etta(1:N,1)))/(v_0*eps_0);
end

