
%%% Constants %%%
w=1;
d=0.01;
v_0=1;
N = 2;
eps_0=8.85*10^(-12);
a=10000;

A = 0;
B = 6;
C = 6;
fprintf("Q3: A B C =  %d %d %d\n", A, B, C);
%%%%%%% MAIN %%%%%%

%%% Q3_a %%%
figure(302);
cap3 = create_cap(N,w,A,B,C);
c_cur = get_c(cap3,N,v_0,eps_0);
plot(N,c_cur,'*b');
title("C_n as a function of W/delta (number of segments)");
hold on
c_prev = c_cur;
N = N*2;
cap3_new = create_cap(N,w,A,B,C);
c_cur = get_c(cap3_new,N,v_0,eps_0);
err = abs((c_cur-c_prev)/c_prev); 
plot(N,c_cur,'*b');
hold on

while(err>=0.01)
    c_prev = c_cur;
    N = N*2;
    cap3_new = create_cap(N,w,A,B,C);
    c_cur = get_c(cap3_new,N,v_0,eps_0);
    err = abs((c_cur-c_prev)/c_prev);
    plot(N,c_cur,'*b');
    hold on
end
xlabel('W/delta [scalar]');
ylabel('C per length [scalar]');
legend('discrete capacitance values');
grid on
fprintf("Q3_a: N =  %d\n",N);

%Plot segment for demonstration
cap3 = create_cap(N,w,A,B,C);
figure(301);
plot(cap3(1:N+1,1)', cap3(1:N+1,2)','g');
hold on
plot(cap3(N+2:2*N+2,1)', cap3(N+2:2*N+2,2)','r');
hold on

delta_mat = delta_cap(cap3,N);
plot(delta_mat(1:N,1)', delta_mat(1:N,2)','pb');
hold on
plot(delta_mat(N+1:2*N,1)', delta_mat(N+1:2*N,2)','pb');
grid on
xlabel('X[m]');
ylabel('Y[m]');

%%% Q3_b %%%

etta_tot = get_eta(cap3, N,v_0,eps_0);
X_vector = linspace(-0.5,0.5,N);
figure(3);
plot(X_vector,etta_tot(N+1:2*N,1));
grid on
title ('etta as a function of x (bottom plate)');
xlabel('X[m]');
ylabel('etta [coulomb/m^2]');

figure(304);
phi_vector = zeros(N,1);
for i = 1:N
    phi_vector(i,1) = (2*pi/N)*i;
end
plot(phi_vector,etta_tot(1:N,1));
grid on
title ('\eta as a function of \phi (top plate)');
xlabel('\phi[rad]');
set(gca,'XTick',0:pi/2:2*pi); 
set(gca,'XTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'});
ylabel('\eta [coulomb/m^2]');

%%% Q3_c %%%

delta = w/N;
top_plate_etta = sum(etta_tot(1:N,1));
bottom_plate_etta = sum(etta_tot(N+1:2*N,1));
Q_3_etta = get_eta(cap3,N,v_0, eps_0);

pixel = 0.05;
phi_s = zeros(2*3*w/pixel,2*3*w/pixel);     %initialing the potential matrix
d_s = zeros(2*N,1);                         %initialing the distance vector
 
%in this long end painstaiking loop we pass upon every pixel in our field
%and calculate its potential, so we have a 2D matrix of every potential at
%every point in space
for j=1:(2*3*w/pixel)
    for k=1:(2*3*w/pixel)
        X_=-3*w+j*pixel-0.5*pixel;
        Y_=3*w-k*pixel+0.5*pixel;
        % x_ and y_ are reletive variables helping us calculate the
        % distance from each pixel
        for l=1:(2*N)
            d_s(l,1)=sqrt((Y_-delta_mat(l,1)).^2 +(X_-delta_mat(l,2)).^2);
        end
        %now calculate for every pixel its potential 
        for l=1:(2*N)
            %pot_loc_mat_q3(y,x) = pot_loc_mat_q3(y,x) + -((etta_B_q3(n,1).*delta_tot_vec(n,1))./(2*pi*eps_0))*log((dist)./a); 
            
            phi_s(j,k)=phi_s(j,k)+log(d_s(l,1)/a)*(-1*Q_3_etta(l,1).*delta_mat(l,3)/(2*pi*eps_0));
        end
    end
end

x=-3*w+0.5*pixel:pixel:3*w-0.5*pixel;
y=-3*w+0.5*pixel:pixel:3*w-0.5*pixel;
figure(305);
contour(x,y,phi_s),title("Equal potential lines(x,y) [V]");
xlabel('X[m]');
ylabel('Y[m]');

[Ex,Ey]= gradient(-phi_s);
figure(306);
quiver(x,y,Ex,Ey),title("Electric vector field(x,y) [V/m]");
xlabel('X[m]');
ylabel('Y[m]');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function cap3 = create_cap(N,W,A,B,C)
% creates capatitor coordinates
    phi = linspace(0,2*pi(),N+1);
    rho = (W/3.6)*(1+0.5*sin(A*phi)+0.2*sin(B*phi)+0.1*sin(C*phi));  % for ID - 203244066 - A = 0, B = 6, C = 6
    cap3 = zeros(2*(length(phi)),2);
    for i = 1:N+1
    % top flower-electrode
        cap3(i,1) = rho(i)*cos(phi(i));
        cap3(i,2) = rho(i)*sin(phi(i));
    end
    v = linspace(-0.5,0.5,N+1);
    for i = N+2:2*N+2
    % bottom plate    
        cap3(i,1) = v(i-N-1);
        cap3(i,2) = -0.6;
    end
end


function delta_mat = delta_cap(cap,N) % for capacitor #3
    % returns delta matrix - 1st - X coor, 2nd - Y coor, 3rd - i'th delta
    delta_mat = zeros(2*N,3);
    for i = 1:N
    % Calc top plate middle segment point
        delta_mat(i,1) = (cap(i,1)+cap(i+1,1))/2; % X coor
        delta_mat(i,2) = (cap(i,2)+cap(i+1,2))/2; % Y coor
        delta_mat(i,3) = sqrt((cap(i,1)-cap(i+1,1)).^2+(cap(i,2)-cap(i+1,2)).^2); %delta
    end
    for i = N+1:2*N
    % Calc bottom plate middle segment point
        delta_mat(i,1) = (cap(i+1,1)+cap(i+2,1))/2; % X coor
        delta_mat(i,2) = (cap(i+1,2)+cap(i+2,2))/2; % Y coor
        delta_mat(i,3) = sqrt((cap(i+1,1)-cap(i+2,1)).^2+(cap(i+1,2)-cap(i+2,2)).^2); %delta
    end 
end


function z = z_nm(cap,a,eps_0,N,n,m)
%returns Z[n,m] coordinate in Z matrix   
    delta_mat = delta_cap(cap,N); % 2Nx2 matrix
    if(n==m)
        z = (log(delta_mat(n,3)./(2*a))-1)*(-0.5)*delta_mat(n,3)./(pi*eps_0);
    else
        r_nm = sqrt((delta_mat(n,1)-delta_mat(m,1))^2+(delta_mat(n,2)-delta_mat(m,2))^2);
        % r_nm is the abs(r_n-r_m)
        z = ((-0.5)*delta_mat(m,3)./(pi*eps_0))*log(r_nm/a);
    end
end


function eta = get_eta(cap,N,V_0,eps_0) % for CAPACITOR 3
% gets back the Etta vector - surface charge density on plates
    a = 10000;
    z_mat = zeros(2*N,2*N); %(2N,2N)
    for i=1:(2*N)
        for j=1:(2*N)
            z_mat(i,j)=z_nm(cap,a,eps_0,N,i,j);
        end
    end
    v0 = zeros(2*N,1);
    for i=1:N
        v0(i,1)=0.5*V_0;
        v0(i+N,1)=(-0.5)*V_0;
    end
    eta = z_mat\v0;
end

function c = get_c(cap,N,V_0,eps_0)
% Gets: Capacitor coordinates - cap, Voltage on electrodes - V_0, 
% number of segments - N, distance between plates - d, and gets back
% its CAPACITY
    eta = get_eta(cap,N,V_0,eps_0);
    c = abs(sum(eta(N+1:2*N), 1))/(N*V_0*eps_0);
end