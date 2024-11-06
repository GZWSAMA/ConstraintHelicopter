clear;clc;close all;

t = 0.005;
dt = 0.005;
i= 1;

e3 = [0,0,1]';
g = 9.8;
m = 8.2;
Jx = 0.18; Jy = 0.34; Jz = 0.28;
J = diag([Jx,Jy,Jz]);

p = [0.2,-0.2,2]';
v = [0,0,0]';
gamma = [0,0,0]';
omega = [0,0,0]';
Inte_bar_R3e_t = [0,0]';
Inte_psie_t = 0;
Inte_omegae_t = [0,0,0]';

xi = 0.07;
wn = 20;

z1 = zeros(2,1);
dz1 = zeros(2,1);

z2 = 0;
dz2 = 0;

z3 = zeros(3,1);
dz3 = zeros(3,1);
%-------------PID------------------------%
kz = 1;
kw = 0.5;
kp = 1.2;
kv = 0.4;
kgammap = 2.12;
kgammai = 2.25;
kpsip = 0.35;
kpsii = 0.06;
komegap = 5;
komegai = 12.96;
%--------------- Parameters of helicopter------------------%
rho = 1.2; Ma = 25.23; Lb = 25.23;
lm = 0.01; hm = 0.14; sm = 0.05; am = 5.4; Am = 2.06; Omegam = 167; 
Rm = 0.81;lt = 0.95; ht = 0.05; st = 0.15; at = 5; At = 0.06; Omegat = 776;
Rt = 0.14;Alon = 4.2; Blat = 4.2; tf = 0.0278;deltad = 0.012;
az = 10;aw = 10;ap = 3;av = 3;

while t<80
    %--------------------参考信号-----------------------------------%
    pr = [9.6*10^(-8)*t^5-1.12*10^(-5)*t^4+3.2*10^(-4)*t^3+0.2;
        -5.76*10^(-8)*t^5+6.4*10^(-6)*t^4-1.6*10^(-4)*t^3-0.2;
        1.152*10^(-7)*t^5-1.44*10^(-5)*t^4+4.8*10^(-4)*t^3];
    dpr = [5*9.6*10^(-8)*t^4-4*1.12*10^(-5)*t^3+3*3.2*10^(-4)*t^2;
        -5*5.76*10^(-8)*t^4+4*6.4*10^(-6)*t^3-3*1.6*10^(-4)*t^2;
        5*1.152*10^(-7)*t^4-4*1.44*10^(-5)*t^3+3*4.8*10^(-4)*t^2];
    ddpr = [4*5*9.6*10^(-8)*t^3-3*4*1.12*10^(-5)*t^2+2*3*3.2*10^(-4)*t;
        -4*5*5.76*10^(-8)*t^3+3*4*6.4*10^(-6)*t^2-2*3*1.6*10^(-4)*t;
        4*5*1.152*10^(-7)*t^3-3*4*1.44*10^(-5)*t^2+2*3*4.8*10^(-4)*t];
    psir = atan2(dpr(2),dpr(1));
    %--------------------信号转换----------------------------------%
    %高度信号
    zr = pr(3);dzr = dpr(3);ddzr = ddpr(3);
    w = v(3);z = p(3);
    %角度信号
    phi = gamma(1);
    theta = gamma(2);
    psi = gamma(3);
    R = [cos(psi)*cos(theta), ...
         cos(psi)*sin(theta)*sin(phi) - sin(psi)*cos(phi), ...
         cos(psi)*sin(theta)*cos(phi) + sin(psi)*sin(phi); ...
         
         sin(psi)*cos(theta), ...
         sin(psi)*sin(theta)*sin(phi) + cos(psi)*cos(phi), ...
         sin(psi)*sin(theta)*cos(phi) - cos(psi)*sin(phi); ...
         
         -sin(theta), ...
         cos(theta)*sin(phi), ...
         cos(theta)*cos(phi)];
    %水平运动信号
    bar_pr = pr(1:2);dbar_pr=dpr(1:2);ddbar_pr=ddpr(1:2);
    bar_p = p(1:2);
    bar_v = v(1:2);
    bar_R3 = [R(1,3),R(2,3)]';
    %姿态信号
    hat_R = [-R(1,2), R(1,1);
        -R(2,2),R(2,1)];
    pp = omega(1);q = omega(2);r = omega(3);
    bar_omega = [pp,q]';
    Ggamma = [hat_R, zeros(2,1);
        zeros(1,2),cos(phi)/cos(theta)];
    %--------------------高度速度回路----------------------------------%
    ze = z - zr;
    we = w - dzr;
    Tm = m*(g+ddzr-kz*tanh(az*ze+aw*we)-kw*tanh(aw*we));
    %--------------------水平速度回路----------------------------------%
    bar_pe = bar_p-bar_pr;
    bar_ve = bar_v - dbar_pr;
    bar_alphap = m/Tm*(ddbar_pr-kp*tanh(ap*bar_pe+av*bar_ve)-kv*tanh(av*bar_ve));
    %--------------------pq姿态回路----------------------------------%
    bar_R3e = bar_R3 - bar_alphap;

    ddz1 = -2*xi*wn*dz1-wn*(z1-bar_alphap);
    dbar_alphap = dz1;

    bar_alphaR = hat_R^(-1)*(-kgammap*bar_R3e-kgammai*Inte_bar_R3e_t+dbar_alphap);
    bar_omegae = bar_omega-bar_alphaR;
    %--------------------r姿态回路----------------------------------%
    psie = psi - psir;

    ddz2 = -2*xi*wn*dz2-wn*(z2-psir);
    dpsir_test = dz2;

    dpsir = -(ddpr(1)*dpr(2)-dpr(1)*ddpr(2))/(dpr(1)^2+dpr(2)^2);

    alphapsi = -sin(phi)/cos(phi)*q-cos(theta)/cos(phi)*(kpsip*psie+kpsii*Inte_psie_t-dpsir);
    re = r - alphapsi;
    %--------------------姿态跟踪----------------------------------%
    alphaR = [bar_alphaR',alphapsi]';
    omegae = [bar_omegae',re]';

    ddz3 = -2*xi*wn*dz3-wn*(z3-alphaR);
    dalphaR = dz3;
    bar_gammae = [bar_R3e',psie]';

    taur = skew(omega)*J*omega+J*dalphaR-komegap*omegae-komegai*Inte_omegae_t-Ggamma*bar_gammae;
    
    %-------------------推力和转矩---------------------------------%
    tcm = Tm/(rho*sm*Am*Omegam^2*Rm^2);
    qcm = deltad/8+1.13*tcm^1.5*sqrt(sm/2);
    Qm = qcm*rho*sm*Am*Omegam^2*Rm^3;
    QA = [ht Qm Tm*hm;
        0 Tm*hm -Qm;
        -lt 0 -Tm*lm];
    taurB = [0;Tm*lm;Qm];
    taurA = QA^(-1)*(taur - taurB);
    Tt = taurA(1);
    as = taurA(2);
    bs = taurA(3);

    f = [Tm*sin(as);
        -Tm*sin(bs)+Tt;
        Tm*(cos(bs)*cos(as))];
    %------------------ODE-----------------------------------%
    dp = v;
    dv = -g*e3+1/m*R*f;
    dR = R*skew(omega);
    domega = J^(-1)*(-skew(omega)*J*omega+taur);
    %-------------------update CF1 -----------------------%
    dz1 = dz1+ddz1*dt;
    z1 = z1+dz1*dt;
    %-------------------update CF2 -----------------------%
    dz2 = dz2+ddz2*dt;
    z2 = z2+dz2*dt;
    %-------------------update CF3 -----------------------%
    dz3 = dz3+ddz3*dt;
    z3 = z3+dz3*dt;
    %-------------------update inte -----------------------%
    Inte_psie_t = Inte_psie_t+psie*dt;
    Inte_bar_R3e_t = Inte_bar_R3e_t+bar_R3e*dt;
    Inte_omegae_t = Inte_omegae_t+omegae*dt;
    
    %--------------------ODE update -------------------------%
    p = p+dp*dt;
    v = v+dv*dt;
    omega = omega+domega*dt;
    R = R +dR*dt;
    gamma = rotm2eul(R,'ZYX');
    t = t+dt;
    i = i+1;
    disp("i =");
    disp(i);

    output(i,1) = t;
    output(i,2:4) = p;
    output(i,5:7) = pr;
    output(i,8) = Tm;
    output(i,9) = dpsir_test;
    output(i,10) = dpsir;

end
pic(output);