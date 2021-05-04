%Q2_a............................................................
w=1;
d=.01;
delta=w/2;
N=(w./delta);
v_0=1;
eps_0=8.85*10^(-12);
a=10000;

%same principle as in Q1: we keep the previous C value and keep dividing
%delta by 2 untill the approximation is good enough

c_pre=calc_C_Q2(delta,N,d);
delta = delta/2;
c_cur=calc_C_Q2(delta,w/delta,d);
figure(1)
err = (c_cur-c_pre)/(c_pre);
while(abs(err) > 0.01)
    delta=delta/2;
    c_pre=c_cur;
    c_cur=calc_C_Q2(delta,w/delta,d);
    err=((c_cur-c_pre)/(c_pre));
    plot(d/delta,c_cur,'xr');
    hold on
end
hold on
title("Q2_a: Capacitance as a function of (d/delta)")
xlabel('d/delta [scalar]');
ylabel('C per length [scalar])');
grid on
legend('discrete capacitance values');
fprintf("Q2_a: delta =  %d\n",delta);

%Q2_b............................................................

%same logic us in Q1_c
N=w/delta;
d1=[0.01, 0.1, 1];
%d2=0.1;
%d3=0.01;
etta1=get_etta_Q2(delta,N,d1(1));
etta2=get_etta_Q2(delta,N,d1(2));
etta3=get_etta_Q2(delta,N,d1(3));
%distance decreases - capasitance increase - charges on plates increase
figure(2)
x= (1:N);

y1=etta3(x,1);
y2=etta3(x+N,1);
subplot(3,1,3);
plot((-d1(3))/2-w+x/N,y1),title("Q2_b: etta  as a function of X, d=1(m)");
hold on
plot((d1(3))/2+x/N,y2);
xlabel('X[m]');
ylabel('etta [coulomb/m^2]');
grid on

y1=etta2(x,1);
y2=etta2(x+N,1);
subplot(3,1,2);
plot((-d1(2))/2-w+x/N,y1),title("Q2_b: etta  as a function of X, d=0.1(m)");
hold on
plot((d1(2))/2+x/N,y2);
xlabel('X[m]');
ylabel('etta [coulomb/m^2]');
grid on

y1=etta1(x,1);
y2=etta1(x+N,1);
subplot(3,1,1);
plot((-d1(1))/2-w+x/N,y1),title("Q2_b: etta  as a function of X, d=0.01(m)");
hold on
plot((d1(1))/2+x/N,y2);
xlabel('X [m]');
ylabel('etta [coulomb/m^2]');
grid on


%Q2_c............................................................
numeri_C=zeros(101,1);
analy_C=zeros(101,1);
j=1; 
x= 0.01:0.01:1.01;
%after setting the desired resolution, we calc 2 vectors - numeric and
%analytic calculation for Capacitance

for d = x
    numeri_C(j,1) =calc_C_Q2(delta,N,d);
    analy_C(j,1)=calc_C_Q2_analytic(d,w);
    j = j+1;
end

figure(3);
plot(x,analy_C),title('Q2_c: Capacitance as a function of d');
hold on
plot(x,numeri_C);
grid on
legend('Analytic', 'Numeric');
xlabel('d [m]');
ylabel('capacitance per length [unitless]');

%now to compare how good was our error
figure(4);
plot(x,abs(numeri_C-analy_C));
title( 'Q2_c: error between numeric and analytic calculations');
xlabel('d [m]');
ylabel('error in C [F]');
grid on


%Q2_d............................................................
d=0.1;
pixel=0.05;
eta=get_etta_Q2(delta,w/delta,d);
phi_s=zeros(2*3*w/pixel,2*3*w/pixel);
d_s=zeros(2*N,1);
%first we build a vector of x,y coordinates, the size of 2N
R_xy=zeros(2*N,2);
for i=1:N                  
    R_xy(i,1)=-d/2-delta*(N-i);
    R_xy(i,2)=0;
    R_xy(i+N,1)=d/2+delta*(i-N);
    R_xy(i+N,2)=0;
end    
for j=1:(2*3*w/pixel)
    for k=1:(2*3*w/pixel)
        Xi=-3*w+j*pixel-0.5*pixel;
        Yi=3*w-k*pixel+0.5*pixel;
        for l=1:(2*N)
            d_s(l,1)=sqrt((Yi-R_xy(l,1)).^2 + (Xi-R_xy(l,2)).^2);
        end
        for l=1:(2*N)
            phi_s(j,k)=phi_s(j,k)+log(d_s(l,1)/a)*(-eta(l,1)*delta)/(2*pi*eps_0);
        end
    end
end

x=-3*w+0.5*pixel:pixel:3*w-0.5*pixel;
y=-3*w+0.5*pixel:pixel:3*w-0.5*pixel;
figure(5);
contour(x-0.5,y,phi_s),title("Equal potential lines (x,y)[V]");
xlabel('X[m]');
ylabel('Y[m]');

[Ex,Ey]= gradient(-phi_s);
figure(6);
quiver(x-0.5,y,Ex,Ey),title("Electric vector field (x,y)[V/m]");
xlabel('X[m]');
ylabel('Y[m]');

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
%&&&&&&&&&&&&&&&&&&&&&&&&&&---------functions of Question 2 --------&&&&&&&&&&&&&%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%


function c=calc_C_Q2_analytic(d,w)
    LS=3*10^8;
    v_0=1;
    eps_0=8.85*10^(-12);
    a=10000;
    WD =(sqrt(d/(d+2*w)));
    if(0<WD<(1/sqrt(2)))
        c=calc_1(WD);
    else
        c=calc_2(WD);
    end        
end

function c=calc_1(WD)
    LS=3*10^8;
    eps_0=8.85*10^(-12);
    c = log((-2*(nthroot(1 - WD.^2, 4) + 1))/(nthroot(1 - WD.^2, 4) - 1))/(377*pi*eps_0*LS); 
end

function c=calc_2(WD)
    LS=3*10^8;
    eps_0=8.85*10^(-12);
    c=1/(120*eps_0*LS*log(-2*(WD+1)/(WD-1)));
end


function r = r_dist_Q2(n,m,N,delta,d)
    if ((n<=N && m<=N)||(n>N && m>N)) %both n,m are on the same plate
        r = abs((m-n)*delta);
    else
        if(n>N)
            n_x=d/2+delta*(n-N);
        else
            n_x=-d/2-delta*(N-n);
        end
        if(m>N)
            m_x=d/2+delta*(m-N);
        else
            m_x=-d/2-delta*(N-m);
        end
        % n,m are on two different plate
        r = abs(n_x-m_x);
    end
end


function z=Z_value(delta,N,d,n,m)
%just like in Q1, we build the Z mat accourding to the momentum method
    eps_0=8.85*10^(-12);
    a=10000;
    if(n==m)
        z=(log(delta./(2*a))-1)*(-0.5)*delta./(pi*eps_0);
    else
        z=log(r_dist_Q2(n,m,N,delta,d)/a)*(-0.5)*delta./(pi*eps_0);
    end
end


function etta=get_etta_Q2(delta,N,d)
    z_mat=zeros(2*N,2*N);
    phi_mat=zeros(2*N,1);
    %now building potential (phi in greek) matrix
    for n=1:N
        phi_mat(n,1)=0.5;
        phi_mat(n+N,1)=-0.5;
    end
    %now building Z matrix
    for n=1:(2*N)
        for m=1:(2*N)
            z_mat(n,m)=Z_value(delta,N,d,n,m);
        end
    end
    etta = (z_mat)^(-1)*phi_mat;
end


function c = calc_C_Q2(delta,N,d)
    v_0=1;
    eps_0=8.85*10^(-12);
    a=10000;
    etta = get_etta_Q2(delta,N,d);
    c = abs(delta*sum(etta(1:N,1)))/(v_0*eps_0);
end
