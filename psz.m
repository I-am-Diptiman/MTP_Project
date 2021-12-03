v=input('input speed=');		% v=input('cutting speed in m/min=')
F=input('feed rate=');		% feed in mm/rev
Flo=input('flow rate=');		% flow rate in ml/hr
P=input('pressure=');		% pressure in bar

Fx = 22.3+(1.83*v)+(2307*F)-(0.754*Flo)+(12.8*P)-(0.00646*v*v)-(2724*F*F)+(0.001087*Flo*Flo)-(1.54*P*P)-(13.99*v*F)+(0.00339*v*Flo)-(0.127*v*P)+(0.36*F*Flo)+(27.0*F*P)+(0.0342*Flo*P);
Fc = -232+(2.18*v)+(3716*F)+(0.021*Flo)+(70.6*P)-(0.00768*v*v)-(6280*F*F)-(0.00071*Flo*Flo)-(4.52*P*P)-(4.94*v*F)+(0.00551*v*Flo)-(0.459*v*P)-(1.24*F*Flo)-(57*F*P)-(0.0050*Flo*P);
Ft = -123+(2.63*v)+(2282*F)-(0.510*Flo)+(37.7*P)-(0.01323*v*v)-(4707*F*F)+(0.00057*Flo*Flo)-(2.75*P*P)-(8.42*v*F)+(0.00543*v*Flo)-(0.288*v*P)-(0.79*F*Flo)+(1.9*F*P)+(0.0103*Flo*P); 
t2mu = 376-(4.52*v)+(563*F)-(0.619*Flo)-(28.0*P)+(0.01810*v*v)+(1310*F*F)-(0.00038*Flo*Flo)+(0.53*P*P)-(13.11*v*F)+(0.00923*v*Flo)+(0.275*v*P)+(3.85*F*Flo)+(200*F*P)-(0.0839*Flo*P);    %chip thickness in micron
TCL= 1.549 - 0.01607*v - 0.68*F - 0.00335*Flo - 0.0456*P + 0.000051*v*v + 2.13*F*F + 0.000002*Flo*Flo - 0.00576*P*P + 0.0119*v*F + 0.000006*v*Flo + 0.000668*v*P+ 0.00753*F*Flo + 0.006*F*P + 0.000183*Flo*P;
MUsl = 0.042+(0.00218*v)+(0.00379*Flo)+(0.1357*P)-(0.000050*v*v)-(0.000003*Flo*Flo)-(0.00508*P*P)-(0.000001*v*Flo)+(0.000063*v*P)-(0.000725*Flo*P);
t2=t2mu*0.001;		%cut-chip thickness
%alpha=input('alpha=') %alpha
alpha = -6;
t1 = F;		%for orthogonal cutting

phi= atand(cosd(alpha)/((t2/t1)-sind(alpha))); %shearangle
E= cosd(alpha)/(sind(phi)*cosd(phi-alpha));    %shear strain
n=0.62;		%hardening coefficient
m=1.2;		%thermal softening coefficient

Beta= (1.73/m)*((n/E)+(0.02/(8180*466*(1+1.328*((11.4*E)/(v*t1))^0.5))))*(n+1-((0.664)*((11.4*E)/(v*t1))^0.5/(1+1.328*((11.4*E)/(v*t1))^0.5)));
delY=0.9*(m*t1*sind(phi))/(Beta*cosd(alpha));
E_rate= (cosd(alpha)/cosd(phi-alpha))*((v*10^3)/(60*delY));
E_refrate=1*10^3;
T_rt=25;                                         %Room temperature
T_melt=1300;                                     %Melting temperature

w= input('width of cut=');
k=11.4;															% k=input('thermal conductivity=');
rho=8.18;														% rho = input('enter the value of density in g/cm^3 =');
c=0.466;														% c = input('enter the value of specific heat of work material in J/gC=');
Cn = input('enter the value of Kronenberg constant=');			%0.7
Fu = Fc/(w*t1*10^(-6)); 										%unit cutting force in N/m^2
h = rho*c*10^(6); 												%volumetric specific heat
Ac = w*t1*10^(-6); 												%chip cross sectional area in m^2
T_psz = (Cn*Fu*(v/60)^0.44*Ac^0.22)/(10*k^0.44*h^0.56);			%tempr at primary shear zone

A=762;															% A=input('A=')
B=1812;															% B=input('B=')
C=0.013;															% C=input('C=')
D=log(E_rate/E_refrate);

Sig_psz = 0.57*[A+B*(E)^n]*[1+(C*D)]*[1-((T_psz-T_rt)/(T_melt-T_rt))^m];	%flow stress

Epsz = (Sig_psz*E)/(n+1)

%==============================================================================================================================

lambda = atand(Ft/Fc);
alpha_new = phi + lambda- alpha;
Lt = (t1*sind(alpha_new))/(cosd(alpha_new+alpha-phi)*sind(phi));

zetta = 2;
% Rho1=input('rho1 in kg/m3=');
Rho1=8180;
% c1=input('enter the value of specific heat of work material in J/kgC=');
c1=466;
% v1=input('cutting speed in m/s=');
v1=(v/60);
% t1m=input('uncut chip thickness in m=');
t1m=(t1/1000);
R=(Rho1*c1*v1*t1m)/k;
Tc=Ft*v1*(t1/t2);
x=10^(0.06-((0.02/t2)*(R*t2/TCL)^0.5)+(0.5*log10(R*t2/TCL)));
Tm=Tc*x;
T_ssz= 25+T_psz+ (0.8*Tm);
Sig_ssz=0.57*[A+B*(E)^n]*[1+(C*D)]*[1-((T_ssz-T_rt)/(T_melt-T_rt))^m];
Sig_not= 4*Sig_psz*((1+zetta)/(2+zetta))*(((cosd(lambda))^2)/sind(2*phi+2*lambda-2*alpha));
Lst=(-((Sig_ssz/(Sig_not*MUsl))^(1/zetta))+1)*TCL;

Vchip = v*(t1/t2);
syms x
Vfun = Vchip*(1-((Lt-x)/(Lt-Lst))^2);

Vreg2 = eval(int(Vfun,Lst,Lt));
Vreg2 = Vreg2/(Lt-Lst);
Vreg3 = Vchip;

fa = 0.56616;
fc = 0.16664;
fk = 0.92035;

Sig_not= 4*Sig_psz*((1+zetta)/(2+zetta))*(((cosd(lambda))^2)/sind(2*phi+2*lambda-2*alpha));

sigma_n = (Sig_not*((1 - (x/TCL))^zetta));

%test = eval(int(sigma_n,Lst,Lt))
power = 1-fk;
sigma_nreg2 = (sigma_n)^(power);

F_expo = w*fa*exp(-1*fc*Vfun);

%F_expo = fa*((Sig_not*((1 - (x/TCL))^zetta))^(1-fk))*exp(-1*fc*Vreg2);

%test = sigma_nreg2*F_expo;
Fric_reg2 = vpa(int(sigma_nreg2*F_expo,Lst,Lt));

%============================================================================================================================


Fric_fun_reg3 = w*fa*sigma_n*exp(-1*fc*Vchip);
Fric_reg3 = vpa(int(Fric_fun_reg3,Lt,TCL));

denom = (Sig_not*((1 - (Lt/TCL))^zetta))^fk;
z1 = Fric_reg2 * Vreg2
z2 = (Fric_reg3*Vchip/denom)
Essz = (Fric_reg2 * Vreg2 + (Fric_reg3*Vchip/denom))/(0.5*v*F)