%Meshless2D_accuracy_verification
clc
tic;
%****************************************Node Settings****************************************%
dl=200*ones(1,100);
dr=dl;
da=200*ones(1,100);
dg=da;
dx=[fliplr(dl),dr];
dz=[fliplr(da),dg];
xr=cumsum(dr);
xl=-fliplr(cumsum(dl));
x=[xl,0,xr];
za=-fliplr(cumsum(da));
zg=cumsum(dg);
z=[za,0,zg];
nx=numel(x);
nz=numel(z);
ndx=numel(dx);
ndz=numel(dz);
[xm,zm]=meshgrid(x,z);
an=nz*nx;
yuan=find(xm==0 & zm==0);
z0=find(zm==0);

%****************************************Electromagnetic parameter settings****************************************%
I=1;
f=1;
miu=4*pi*1e-7;
eps=1/36/pi*1e-9;
sig=1/100;
siga=1e-12;

%****************************************Gauss integral settings****************************************%
gs=4;
xk=[-0.8611363116,-0.3399810436,0.3399810436,0.8611363116];
Hk=[0.3478548451,0.65214515497,0.65214515497,0.3478548451];
xi=reshape((0.5*x(2:end)'*(xk+1)+0.5*x(1:end-1)'*(1-xk))',1,gs*ndx);
zi=reshape((0.5*z(2:end)'*(xk+1)+0.5*z(1:end-1)'*(1-xk))',1,gs*ndz);
nxi=numel(xi);
nzi=numel(zi);
[xim,zim]=meshgrid(xi,zi);
Him=repmat(Hk'*Hk,ndz,ndx);
bl=dz'*dx/4;
bli=kron(bl,ones(gs,gs));
sigma=sig*ones(nzi,nxi);
sigma(zim<0)=siga;

%****************************************Meshless parameter settings****************************************
ac=0.1; %MQ
q=0.5;
sx=ceil((1:numel(xi))/gs);
sz=ceil((1:numel(zi))/gs);
sn0=(sx-1)*nz+sz';
nsn=4;

%****************************************Frequency domain calculation****************************************%
for nf=1:length(f)
    omg=2*pi*f(nf);
    tau=1/(1i*omg*miu);
    lamda=(sigma+1i*omg*eps);
    K=zeros(an,an);
    P=zeros(an,1);
    P(yuan)=-I;
    for m=1:nzi
        for n=1:nxi
            xjs=xim(m,n);
            zjs=zim(m,n);
            sn1=sn0(m,n);
            sn=[sn1;sn1+1;sn1+nz;sn1+nz+1];
            xjd=xm(sn)';
            zjd=zm(sn)';
            dc2=40000;
            cm=ac^2*dc2;
            r=((xjs-xjd).^2+(zjs-zjd).^2+cm).^q;
            rx=2*q*(xjs-xjd).*((xjs-xjd).^2+(zjs-zjd).^2+cm).^(q-1);
            rz=2*q*(zjs-zjd).*((xjs-xjd).^2+(zjs-zjd).^2+cm).^(q-1);
            r0=((xjd'-xjd).^2+(zjd'-zjd).^2+cm).^q;
            p=[1 xjs zjs];
            px=[0 1 0];
            pz=[0 0 1];
            pm=[ones(1,nsn);xjd;zjd];
            g=[r0 pm'
                pm zeros(3,3)];
            fi=[r p]/g;
            fix=[rx px]/g;
            fiz=[rz pz]/g;
            fi(:,(nsn+1):end)=[];
            fix(:,(nsn+1):end)=[];
            fiz(:,(nsn+1):end)=[];
            fij=fi'*fi;
            fixx=fix'*fix;
            fizz=fiz'*fiz;
            Ka=bli(m,n)*Him(m,n)*tau*(fixx+fizz);
            Kb=bli(m,n)*Him(m,n)*lamda(m,n)*fij;
            K(sn,sn)=K(sn,sn)+Ka+Kb;
        end
    end
    K=sparse(K);
    P=sparse(P);
    eyxs=K\P;
    eyl=full(eyxs);
    ey=eyl(z0);
end
toc;


