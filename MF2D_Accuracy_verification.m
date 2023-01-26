%无网格法和李俊杰2015年的应用地球物理上的论文做对比
clc
tic;
%****************************************节点设置****************************************%
dl=200*ones(1,100);
dr=dl;
% da=[0.25*ones(1,4),0.5*ones(1,2),1*ones(1,2),2*ones(1,2),4*ones(1,2),8*ones(1,2),16*ones(1,2),32*ones(1,2)...
%     ,64*ones(1,2),128*ones(1,2),256*ones(1,2),512*ones(1,2),1024*ones(1,2),2048*ones(1,2)];
% dg=[0.25*ones(1,4),0.5*ones(1,2),1*ones(1,2),2*ones(1,2),4*ones(1,2),8*ones(1,2),16*ones(1,2),32*ones(1,2)...
%     ,64*ones(1,2),128*ones(1,2),256*ones(1,2),512*ones(1,2),1024*ones(1,2),2048*ones(1,2),4096*ones(1,2)...
%     ,8192*ones(1,2),16384*ones(1,2),32768*ones(1,2),65536*ones(1,2)];
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
[xm,zm]=meshgrid(x,z); %组成二维坐标矩阵
an=nz*nx; %节点总个数
yuan=find(xm==0 & zm==0);
delta=min(dz);
z0=find(zm==0); %地表的z坐标

dbc=find(zm==z(end) | zm==z(1) | xm==x(end) | xm==x(1));
% dbc0=find(zm==0 & (xm==x(end) | xm==x(1)));

%****************************************电磁参数设置****************************************%
I=1; %电流
f=1; %频点
miu=4*pi*1e-7;
eps=1/36/pi*1e-9;
sig=1/100; %大地电导率
siga=1e-8; %空气层电导率


%****************************************高斯积分设置****************************************%
gs=4;
xk=[-0.8611363116,-0.3399810436,0.3399810436,0.8611363116];
Hk=[0.3478548451,0.65214515497,0.65214515497,0.3478548451];
% gs=2;
% xk=[-0.5773502691896250 0.5773502691896250];
% Hk=[1 1];
%面积分高斯点及其系数
xi=reshape((0.5*x(2:end)'*(xk+1)+0.5*x(1:end-1)'*(1-xk))',1,gs*ndx);
zi=reshape((0.5*z(2:end)'*(xk+1)+0.5*z(1:end-1)'*(1-xk))',1,gs*ndz);
nxi=numel(xi);
nzi=numel(zi);
[xim,zim]=meshgrid(xi,zi); %组成二维高斯点矩阵
Him=repmat(Hk'*Hk,ndz,ndx);
bl=dz'*dx/4;
bli=kron(bl,ones(gs,gs)); %对应于每个面积内积分点的尺度变换系数
sigma=sig*ones(nzi,nxi);
sigma(zim<0)=siga;


%****************************************无网格参数设置****************************************
ac=0.1; %MQ形函数的参数
q=0.5;
sx=ceil((1:numel(xi))/gs);
sz=ceil((1:numel(zi))/gs);
sn0=(sx-1)*nz+sz'; %找到各个积分点所在背景网格的左上角节点
nsn=4;


%****************************************代入频点计算****************************************%
for nf=1:length(f)
    %和频率有关的参数
    omg=2*pi*f(nf); %角频率
    tau=1/(1i*omg*miu);
    lamda=(sigma+1i*omg*eps);
    %解析解设置+1i*omg*eps
    Eyjx=Function_LineS_Ey(I,f(nf),sig,x,nx);
    %矩阵初始化
    K=zeros(an,an);
    P=zeros(an,1);
    P(yuan)=-I;
    px=[0 1 0];
    pz=[0 0 1];
    %面积分
    for m=1:nzi
        for n=1:nxi
            %计算点坐标
            xjs=xim(m,n);
            zjs=zim(m,n);
            sn1=sn0(m,n); %找到各个积分点所在背景网格的左上角节点
            sn=[sn1;sn1+1;sn1+nz;sn1+nz+1]; %各个积分点所在背景网格的四个节点
            xjd=xm(sn)';
            zjd=zm(sn)';
%             dc2=(xjd(3)-xjd(1))^2+(zjd(2)-zjd(1))^2;
            dc2=40000;
            cm=ac^2*dc2;
            r=((xjs-xjd).^2+(zjs-zjd).^2+cm).^q; %径向基函数
            rx=2*q*(xjs-xjd).*((xjs-xjd).^2+(zjs-zjd).^2+cm).^(q-1);
            rz=2*q*(zjs-zjd).*((xjs-xjd).^2+(zjs-zjd).^2+cm).^(q-1);
            r0=((xjd'-xjd).^2+(zjd'-zjd).^2+cm).^q;
            p=[1 xjs zjs];
            pm=[ones(1,nsn);xjd;zjd];
            g=[r0 pm'
                pm zeros(3,3)];
            fi=[r p]/g; %形函数
            fix=[rx px]/g;
            fiz=[rz pz]/g;
            fi(:,(nsn+1):end)=[];
            fix(:,(nsn+1):end)=[];
            fiz(:,(nsn+1):end)=[];
            fij=fi'*fi;
            fixx=fix'*fix;
            fizz=fiz'*fiz;
            Ka=bli(m,n)*Him(m,n)*tau*(fixx+fizz); %刚度矩阵
            Kb=bli(m,n)*Him(m,n)*lamda(m,n)*fij;
            K(sn,sn)=K(sn,sn)+Ka+Kb; %把积分累加进去
        end
    end
%     K(dbc,:)=0;
%     for mm=1:length(dbc)
%         K(dbc(mm),dbc(mm))=1;
%     end
%     P(dbc0)=Eyjx(end);
    %求解
    K=sparse(K);
    P=sparse(P);
    eyxs=K\P;
    eyl=full(eyxs);
    % eym=reshape(eyl,nz,nx);
    ey=eyl(z0);
%     hx(nf,:)=tau*(-50/24*eyl(z0)+4*eyl(z1)-3*eyl(z2)+4/3*eyl(z3)-1/4*eyl(z4))/(zg(2)-zg(1));
%     Z(nf,:)=-ey(nf,:)./hx(nf,:);
%     rho(nf,:)=1/omg/miu*(abs(Z(nf,:))).^2;
%     phi(nf,:)=atan(imag(Z(nf,:))./real(Z(nf,:)))*180/pi;
end
ey=ey.';
%数值解
eyr=real(ey);
eyi=imag(ey);
%解析解
eyjxr=real(Eyjx);
eyjxi=imag(Eyjx);
%误差
wcr=abs((eyr-eyjxr)./eyjxr)*1e2;
wci=abs((eyi-eyjxi)./eyjxi)*1e2;
wc=(abs(ey)-abs(Eyjx))./abs(Eyjx)*1e2;
% xoffset=find(x==6e3);
% rho(:,xoffset)
% phi(:,xoffset)
toc;


