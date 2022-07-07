clear all 
close all
rng shuffle
clc

T=100;
time=1:T;

a=zeros(1,T);
o=zeros(1,T);
ml=zeros(1,T);
mr=zeros(1,T);
pl=zeros(1,T);
pr=zeros(1,T);
Gl=zeros(1,T);
Gr=zeros(1,T);
fl=zeros(1,T);
fr=zeros(1,T);


tt=1:T;

Pl=(tt<=T/2)*0.5+(tt>T/2)*0.5;
Pr=(tt<=T/2)*1+(tt>T/2)*1;

vw=0.4;
am=0.05;

po=0.1;
mo=0;

%initial value
a(1)=rand<0.5;

ml(1)=mo; 
pl(1)=po;

mr(1)=mo;
pr(1)=po;

dt=1;

Pofix=[0.5:0.01:0.99,0.999999999];%x-axis
Cfixtmp=[-10:0.5:10];
Cfix=sort(Cfixtmp,'descend');%y-axis


for l=1:length(Pofix)
for k=1:length(Cfix)
C=Cfix(k)*ones(1,T);
Po=Pofix(l);

for i=1:10000

for t=2:T
   
    if t==2
        a(t)=rand<0.5;
    else
       al_prob=1/(1+(exp(-(Gr(t-1)-Gl(t-1)))));
       a(t)=rand<al_prob;
    end
    
    fl(t)=(1/(1+exp(-ml(t-1))));   
    fr(t)=(1/(1+exp(-mr(t-1)))); 
    
    if a(t)==1 
        o(t) = rand < Pl(t);
        
        dmldt=((1/pl(t-1))+vw)*(o(t)-fl(t));
        ml(t)=ml(t-1)+am*dmldt*dt;
        Sl(t)=(1/(1+exp(-ml(t))));
        pl(t)=((1/vw*pl(t-1))/(pl(t-1)+1/vw))+(Sl(t)*(1-Sl(t)));
        
        mr(t)=mr(t-1);
        pr(t)=((1/vw*pr(t-1))/(pr(t-1)+1/vw));
    
    else
        o(t) = rand < Pr(t);
        
        dmrdt=((1/pr(t-1))+vw)*(o(t)-fr(t));
        mr(t)=mr(t-1)+am*dmrdt*dt;
        Sr(t)=(1/(1+exp(-mr(t))));
        pr(t)=((1/vw*pr(t-1))/(pr(t-1)+1/vw))+(Sr(t)*(1-Sr(t)));
        
        ml(t)=ml(t-1);
        pl(t)=((1/vw*pl(t-1))/(pl(t-1)+1/vw));
         
    end 
    lamdal(t)=1/(1+exp(-ml(t)));
    lamdar(t)=1/(1+exp(-mr(t)));
   
  %%%%%%%%%%%%%expected free energy%%%%%%%%%%%%%%
  
POAl(t)=lamdal(t)+0.5*lamdal(t)*(1-lamdal(t))*(1-2*lamdal(t))*((1/pl(t))+vw);  
Al(t)= -lamdal(t)*log(lamdal(t))-(1-lamdal(t))*log(1-lamdal(t));
Bl(t)=-0.5*(lamdal(t)*(1-lamdal(t))*(1+(1-2*lamdal(t))*(log(lamdal(t))-log(1-lamdal(t)))))*((1/pl(t))+vw);
Cl(t)=(1-POAl(t))*log(1-POAl(t))+POAl(t)*log(POAl(t));
Dl(t)=-POAl(t)*log(Po/(1-Po))-(1-POAl(t))*0;

POAr(t)=lamdar(t)+0.5*lamdar(t)*(1-lamdar(t))*(1-2*lamdar(t))*((1/pr(t))+2*vw);
Ar(t)= -lamdar(t)*log(lamdar(t))-(1-lamdar(t))*log(1-lamdar(t));
Br(t)=-0.5*(lamdar(t)*(1-lamdar(t))*(1+(1-2*lamdar(t))*(log(lamdar(t))-log(1-lamdar(t)))))*((1/pr(t))+vw);
Cr(t)=(1-POAr(t))*log(1-POAr(t))+POAr(t)*log(POAr(t));
Dr(t)=-POAr(t)*log(Po/(1-Po))-(1-POAr(t))*0;

Gl(t)=C(t)*(Al(t)+Bl(t)+Cl(t))+Dl(t);
Gr(t)=C(t)*(Ar(t)+Br(t)+Cr(t))+Dr(t);

Precisionl(t)=pl(t-1)/(fl(t)^2)./((1-fl(t))^2);
Precisionr(t)=pr(t-1)/(fr(t)^2)./((1-fr(t))^2); 

end

Times_of_al_tmp=sum(a);
Times_of_ar_tmp=T-Times_of_al_tmp;

Rate_of_al_tmp(i)=Times_of_al_tmp/T;
Rate_of_ar_tmp(i)=Times_of_ar_tmp/T;
end

Rate_of_al_ave(k,l)=mean(Rate_of_al_tmp);
Rate_of_ar_ave(k,l)=mean(Rate_of_ar_tmp);


end
end



f = figure;
f.Position(3:4) = [500 1500];


subplot(4,2,1:4)
h=imagesc(Pofix,Cfix,Rate_of_ar_ave)
axis xy
grid off
colormap(redblue)
caxis([0,1])

set(gca, 'XTick', [0.5,0.75,1], 'XTickLabel', [0.5,0.75,1]) % 10 ticks 
set(gca, 'YTick', [-10,0,10], 'YTickLabel', [-10,0,10]) % 20 ticks


hold on
plot(0.75,-10,'k .','MarkerSize',30)
hold on
plot(0.5,10,'k .','MarkerSize',30)
hold on
plot(1,0,'k .','MarkerSize',30)


%% fig3B'-a

clear all
T=100;
time=1:T;

a=zeros(1,T);
o=zeros(1,T);
ml=zeros(1,T);
mr=zeros(1,T);
pl=zeros(1,T);
pr=zeros(1,T);
Gl=zeros(1,T);
Gr=zeros(1,T);
fl=zeros(1,T);
fr=zeros(1,T);

I=1000;
Po=0.999999999;
C=1*ones(1,T);


tt=1:T;
Pl=(tt<=T/2)*0.5+(tt>T/2)*0.5;
Pr=(tt<=T/2)*1+(tt>T/2)*1;

  

vw=0.5;
am=0.3;

po=0.1;
mo=0;

%initial value
a(1)=rand<0.5;

ml(1)=mo; 
pl(1)=po;

mr(1)=mo;
pr(1)=po;

dt=1;

for i=1:I

for t=2:T
   
    if t==2
        a(t)=rand<0.5;
    else
       al_prob=1/(1+(exp(-(Gr(t-1)-Gl(t-1)))));
       a(t)=rand<al_prob;
    end
    
    fl(t)=(1/(1+exp(-ml(t-1))));   
    fr(t)=(1/(1+exp(-mr(t-1)))); 
    
    if a(t)==1 
        o(t) = rand < Pl(t);
        
        dmldt=((1/pl(t-1))+vw)*(o(t)-fl(t));
        ml(t)=ml(t-1)+am*dmldt*dt;
        Sl(t)=(1/(1+exp(-ml(t))));
        pl(t)=((1/vw*pl(t-1))/(pl(t-1)+1/vw))+(Sl(t)*(1-Sl(t)));
        
        
        mr(t)=mr(t-1);
        pr(t)=((1/vw*pr(t-1))/(pr(t-1)+1/vw));
    
    else
        o(t) = rand < Pr(t);
        
        dmrdt=((1/pr(t-1))+vw)*(o(t)-fr(t));
        mr(t)=mr(t-1)+am*dmrdt*dt;
        Sr(t)=(1/(1+exp(-mr(t))));
        pr(t)=((1/vw*pr(t-1))/(pr(t-1)+1/vw))+(Sr(t)*(1-Sr(t)));
        
        
        ml(t)=ml(t-1);
        pl(t)=((1/vw*pl(t-1))/(pl(t-1)+1/vw));
         
    end 
    lamdal(t)=1/(1+exp(-ml(t)));
    lamdar(t)=1/(1+exp(-mr(t)));
   
  %%%%%%%%%%%%%expected free energy%%%%%%%%%%%%%%
  
POAl(t)=lamdal(t)+0.5*lamdal(t)*(1-lamdal(t))*(1-2*lamdal(t))*((1/pl(t))+vw);  
Al(t)= -lamdal(t)*log(lamdal(t))-(1-lamdal(t))*log(1-lamdal(t));
Bl(t)=-0.5*(lamdal(t)*(1-lamdal(t))*(1+(1-2*lamdal(t))*(log(lamdal(t))-log(1-lamdal(t)))))*((1/pl(t))+vw);
Cl(t)=(1-POAl(t))*log(1-POAl(t))+POAl(t)*log(POAl(t));
Dl(t)=-POAl(t)*log(Po/(1-Po))-(1-POAl(t))*0;

POAr(t)=lamdar(t)+0.5*lamdar(t)*(1-lamdar(t))*(1-2*lamdar(t))*((1/pr(t))+2*vw);
Ar(t)= -lamdar(t)*log(lamdar(t))-(1-lamdar(t))*log(1-lamdar(t));
Br(t)=-0.5*(lamdar(t)*(1-lamdar(t))*(1+(1-2*lamdar(t))*(log(lamdar(t))-log(1-lamdar(t)))))*((1/pr(t))+vw);
Cr(t)=(1-POAr(t))*log(1-POAr(t))+POAr(t)*log(POAr(t));
Dr(t)=-POAr(t)*log(Po/(1-Po))-(1-POAr(t))*0;

Gl(t)=C(t)*(Al(t)+Bl(t)+Cl(t))+Dl(t);
Gr(t)=C(t)*(Ar(t)+Br(t)+Cr(t))+Dr(t);

Precisionl(t)=pl(t-1)/(fl(t)^2)./((1-fl(t))^2);
Precisionr(t)=pr(t-1)/(fr(t)^2)./((1-fr(t))^2); 

end

Times_of_al_tmp=sum(a);
Times_of_ar_tmp=T-Times_of_al_tmp;

Rate_of_al_tmp(i)=Times_of_al_tmp/T;
Rate_of_ar_tmp(i)=Times_of_ar_tmp/T;

end

Rate_of_al_ave=mean(Rate_of_al_tmp);
Rate_of_ar_ave=mean(Rate_of_ar_tmp)
Selections=[Rate_of_al_ave,Rate_of_ar_ave];

SD_al=sqrt(var(Rate_of_al_tmp));
SD_ar=sqrt(var(Rate_of_ar_tmp));
SD=[SD_al,SD_ar];


name=categorical({'Left';'Right'});
subplot(4,2,5)
b=bar(name,Selections,'stack')
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];
ylim([0,1])
yticks([0,0.5,1])
ylabel('Ratio','Fontsize',20,'Fontweight','bold')
hold on

%% fig3B'-b
clear all
T=100;
time=1:T;

a=zeros(1,T);
o=zeros(1,T);
ml=zeros(1,T);
mr=zeros(1,T);
pl=zeros(1,T);
pr=zeros(1,T);
Gl=zeros(1,T);
Gr=zeros(1,T);
fl=zeros(1,T);
fr=zeros(1,T);

I=1000;
Po=0.5;
C=10*ones(1,T);


tt=1:T;
Pl=(tt<=T/2)*0.5+(tt>T/2)*0.5;
Pr=(tt<=T/2)*1+(tt>T/2)*1;

  

vw=0.5;
am=0.3;
po=0.1;
mo=0;

%initial value
a(1)=rand<0.5;

ml(1)=mo; 
pl(1)=po;

mr(1)=mo;
pr(1)=po;

dt=1;

for i=1:I

for t=2:T
   
    if t==2
        a(t)=rand<0.5;
    else
       al_prob=1/(1+(exp(-(Gr(t-1)-Gl(t-1)))));
       a(t)=rand<al_prob;
    end
    
    fl(t)=(1/(1+exp(-ml(t-1))));   
    fr(t)=(1/(1+exp(-mr(t-1)))); 
    
    if a(t)==1 
        o(t) = rand < Pl(t);
        
        dmldt=((1/pl(t-1))+vw)*(o(t)-fl(t));
        ml(t)=ml(t-1)+am*dmldt*dt;
        Sl(t)=(1/(1+exp(-ml(t))));
        pl(t)=((1/vw*pl(t-1))/(pl(t-1)+1/vw))+(Sl(t)*(1-Sl(t)));
        
        mr(t)=mr(t-1);
        pr(t)=((1/vw*pr(t-1))/(pr(t-1)+1/vw));
    
    else
        o(t) = rand < Pr(t);
        
        dmrdt=((1/pr(t-1))+vw)*(o(t)-fr(t));
        mr(t)=mr(t-1)+am*dmrdt*dt;
        Sr(t)=(1/(1+exp(-mr(t))));
        pr(t)=((1/vw*pr(t-1))/(pr(t-1)+1/vw))+(Sr(t)*(1-Sr(t)));
        
        ml(t)=ml(t-1);
        pl(t)=((1/vw*pl(t-1))/(pl(t-1)+1/vw));
         
    end 
    lamdal(t)=1/(1+exp(-ml(t)));
    lamdar(t)=1/(1+exp(-mr(t)));
   
  %%%%%%%%%%%%%expected free energy%%%%%%%%%%%%%%
  
POAl(t)=lamdal(t)+0.5*lamdal(t)*(1-lamdal(t))*(1-2*lamdal(t))*((1/pl(t))+vw);  
Al(t)= -lamdal(t)*log(lamdal(t))-(1-lamdal(t))*log(1-lamdal(t));
Bl(t)=-0.5*(lamdal(t)*(1-lamdal(t))*(1+(1-2*lamdal(t))*(log(lamdal(t))-log(1-lamdal(t)))))*((1/pl(t))+vw);
Cl(t)=(1-POAl(t))*log(1-POAl(t))+POAl(t)*log(POAl(t));
Dl(t)=-POAl(t)*log(Po/(1-Po))-(1-POAl(t))*0;

POAr(t)=lamdar(t)+0.5*lamdar(t)*(1-lamdar(t))*(1-2*lamdar(t))*((1/pr(t))+2*vw);
Ar(t)= -lamdar(t)*log(lamdar(t))-(1-lamdar(t))*log(1-lamdar(t));
Br(t)=-0.5*(lamdar(t)*(1-lamdar(t))*(1+(1-2*lamdar(t))*(log(lamdar(t))-log(1-lamdar(t)))))*((1/pr(t))+vw);
Cr(t)=(1-POAr(t))*log(1-POAr(t))+POAr(t)*log(POAr(t));
Dr(t)=-POAr(t)*log(Po/(1-Po))-(1-POAr(t))*0;

Gl(t)=C(t)*(Al(t)+Bl(t)+Cl(t))+Dl(t);
Gr(t)=C(t)*(Ar(t)+Br(t)+Cr(t))+Dr(t);

Precisionl(t)=pl(t-1)/(fl(t)^2)./((1-fl(t))^2);
Precisionr(t)=pr(t-1)/(fr(t)^2)./((1-fr(t))^2); 

end

Times_of_al_tmp=sum(a);
Times_of_ar_tmp=T-Times_of_al_tmp;

Rate_of_al_tmp(i)=Times_of_al_tmp/T;
Rate_of_ar_tmp(i)=Times_of_ar_tmp/T;

end

Rate_of_al_ave=mean(Rate_of_al_tmp);
Rate_of_ar_ave=mean(Rate_of_ar_tmp)
Selections=[Rate_of_al_ave,Rate_of_ar_ave];

SD_al=sqrt(var(Rate_of_al_tmp));
SD_ar=sqrt(var(Rate_of_ar_tmp));
SD=[SD_al,SD_ar];


name=categorical({'Left';'Right'});
subplot(4,2,6)
b=bar(name,Selections,'stack')
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];
ylim([0,1])
yticks([0,0.5,1])
ylabel('Ratio','Fontsize',20,'Fontweight','bold')
hold on

%% fig3B'-c
clear all

a(1)=1;

T=100;
time=1:T;

a=zeros(1,T);
o=zeros(1,T);
ml=zeros(1,T);
mr=zeros(1,T);
pl=zeros(1,T);
pr=zeros(1,T);
Gl=zeros(1,T);
Gr=zeros(1,T);
fl=zeros(1,T);
fr=zeros(1,T);


I=1;
C=-10*ones(1,T);
Po=0.75;


tt=1:T;
Pl=(tt<=T/2)*0.5+(tt>T/2)*0.5;
Pr=(tt<=T/2)*1+(tt>T/2)*1;

  

vw=0.5;
am=0.3;
po=0.1;
mo=0;



ml(1)=mo; 
pl(1)=po;

mr(1)=mo;
pr(1)=po;

dt=1;

for i=1:I

for t=2:T
   
    if t==2
        a(t)=rand<0.5;
    else
       al_prob=1/(1+(exp(-(Gr(t-1)-Gl(t-1)))));
       a(t)=rand<al_prob;
    end
    
    fl(t)=(1/(1+exp(-ml(t-1))));   
    fr(t)=(1/(1+exp(-mr(t-1)))); 
    
    if a(t)==1 
        o(t) = rand < Pl(t);
        
        dmldt=((1/pl(t-1))+vw)*(o(t)-fl(t));
        ml(t)=ml(t-1)+am*dmldt*dt;
        Sl(t)=(1/(1+exp(-ml(t))));
        pl(t)=((1/vw*pl(t-1))/(pl(t-1)+1/vw))+(Sl(t)*(1-Sl(t)));
        
        mr(t)=mr(t-1);
        pr(t)=((1/vw*pr(t-1))/(pr(t-1)+1/vw));
    
    else
        o(t) = rand < Pr(t);
        
        dmrdt=((1/pr(t-1))+vw)*(o(t)-fr(t));
        mr(t)=mr(t-1)+am*dmrdt*dt;
        Sr(t)=(1/(1+exp(-mr(t))));
        pr(t)=((1/vw*pr(t-1))/(pr(t-1)+1/vw))+(Sr(t)*(1-Sr(t)));
        
        ml(t)=ml(t-1);
        pl(t)=((1/vw*pl(t-1))/(pl(t-1)+1/vw));
         
    end 
    lamdal(t)=1/(1+exp(-ml(t)));
    lamdar(t)=1/(1+exp(-mr(t)));
   
  %%%%%%%%%%%%%expected free energy%%%%%%%%%%%%%%
  
POAl(t)=lamdal(t)+0.5*lamdal(t)*(1-lamdal(t))*(1-2*lamdal(t))*((1/pl(t))+vw);  
Al(t)= -lamdal(t)*log(lamdal(t))-(1-lamdal(t))*log(1-lamdal(t));
Bl(t)=-0.5*(lamdal(t)*(1-lamdal(t))*(1+(1-2*lamdal(t))*(log(lamdal(t))-log(1-lamdal(t)))))*((1/pl(t))+vw);
Cl(t)=(1-POAl(t))*log(1-POAl(t))+POAl(t)*log(POAl(t));
Dl(t)=-POAl(t)*log(Po/(1-Po))-(1-POAl(t))*0;


POAr(t)=lamdar(t)+0.5*lamdar(t)*(1-lamdar(t))*(1-2*lamdar(t))*((1/pr(t))+2*vw);
Ar(t)= -lamdar(t)*log(lamdar(t))-(1-lamdar(t))*log(1-lamdar(t));
Br(t)=-0.5*(lamdar(t)*(1-lamdar(t))*(1+(1-2*lamdar(t))*(log(lamdar(t))-log(1-lamdar(t)))))*((1/pr(t))+vw);
Cr(t)=(1-POAr(t))*log(1-POAr(t))+POAr(t)*log(POAr(t));
Dr(t)=-POAr(t)*log(Po/(1-Po))-(1-POAr(t))*0;

Gl(t)=C(t)*(Al(t)+Bl(t)+Cl(t))+Dl(t);
Gr(t)=C(t)*(Ar(t)+Br(t)+Cr(t))+Dr(t);

Precisionl(t)=pl(t-1)/(fl(t)^2)./((1-fl(t))^2);
Precisionr(t)=pr(t-1)/(fr(t)^2)./((1-fr(t))^2); 

end

Times_of_al_tmp=sum(a);
Times_of_ar_tmp=T-Times_of_al_tmp;

Rate_of_al_tmp(i)=Times_of_al_tmp/T;
Rate_of_ar_tmp(i)=Times_of_ar_tmp/T;

end

Rate_of_al_ave=mean(Rate_of_al_tmp);
Rate_of_ar_ave=mean(Rate_of_ar_tmp)
Selections=[Rate_of_al_ave,Rate_of_ar_ave];

SD_al=sqrt(var(Rate_of_al_tmp));
SD_ar=sqrt(var(Rate_of_ar_tmp));
SD=[SD_al,SD_ar];


name=categorical({'Left';'Right'});
subplot(4,2,7)
b=bar(name,Selections,'stack')
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];
ylim([0,1])
yticks([0,0.5,1])
ylabel('Ratio','Fontsize',20,'Fontweight','bold')
hold on

er = errorbar(1:2,Selections,SD,SD);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off




%% fig3B'-d

clear all

%initial value
a(1)=0; %‰E

T=100;
time=1:T;

a=zeros(1,T);
o=zeros(1,T);
ml=zeros(1,T);
mr=zeros(1,T);
pl=zeros(1,T);
pr=zeros(1,T);
Gl=zeros(1,T);
Gr=zeros(1,T);
fl=zeros(1,T);
fr=zeros(1,T);


I=1;
C=-10*ones(1,T);
Po=0.75;


tt=1:T;
Pl=(tt<=T/2)*0.5+(tt>T/2)*0.5;
Pr=(tt<=T/2)*1+(tt>T/2)*1;

  

vw=0.5;
am=0.3;

po=0.1;
mo=0;



ml(1)=mo; 
pl(1)=po;

mr(1)=mo;
pr(1)=po;

dt=1;

for i=1:I

for t=2:T
   
    if t==2
        a(t)=rand<0.5;
    else
       al_prob=1/(1+(exp(-(Gr(t-1)-Gl(t-1)))));
       a(t)=rand<al_prob;
    end
    
    fl(t)=(1/(1+exp(-ml(t-1))));   
    fr(t)=(1/(1+exp(-mr(t-1)))); 
    
    if a(t)==1 
        o(t) = rand < Pl(t);
        
        dmldt=((1/pl(t-1))+vw)*(o(t)-fl(t));
        ml(t)=ml(t-1)+am*dmldt*dt;
        Sl(t)=(1/(1+exp(-ml(t))));
        pl(t)=((1/vw*pl(t-1))/(pl(t-1)+1/vw))+(Sl(t)*(1-Sl(t)));
        
        mr(t)=mr(t-1);
        pr(t)=((1/vw*pr(t-1))/(pr(t-1)+1/vw));
    
    else
        o(t) = rand < Pr(t);
        
        dmrdt=((1/pr(t-1))+vw)*(o(t)-fr(t));
        mr(t)=mr(t-1)+am*dmrdt*dt;
        Sr(t)=(1/(1+exp(-mr(t))));
        pr(t)=((1/vw*pr(t-1))/(pr(t-1)+1/vw))+(Sr(t)*(1-Sr(t)));
        
        
        ml(t)=ml(t-1);
        pl(t)=((1/vw*pl(t-1))/(pl(t-1)+1/vw));
         
    end 
    lamdal(t)=1/(1+exp(-ml(t)));
    lamdar(t)=1/(1+exp(-mr(t)));
   
  %%%%%%%%%%%%%expected free energy%%%%%%%%%%%%%%
  
POAl(t)=lamdal(t)+0.5*lamdal(t)*(1-lamdal(t))*(1-2*lamdal(t))*((1/pl(t))+vw);  
Al(t)= -lamdal(t)*log(lamdal(t))-(1-lamdal(t))*log(1-lamdal(t));
Bl(t)=-0.5*(lamdal(t)*(1-lamdal(t))*(1+(1-2*lamdal(t))*(log(lamdal(t))-log(1-lamdal(t)))))*((1/pl(t))+vw);
Cl(t)=(1-POAl(t))*log(1-POAl(t))+POAl(t)*log(POAl(t));
Dl(t)=-POAl(t)*log(Po/(1-Po))-(1-POAl(t))*0;


POAr(t)=lamdar(t)+0.5*lamdar(t)*(1-lamdar(t))*(1-2*lamdar(t))*((1/pr(t))+2*vw);
Ar(t)= -lamdar(t)*log(lamdar(t))-(1-lamdar(t))*log(1-lamdar(t));
Br(t)=-0.5*(lamdar(t)*(1-lamdar(t))*(1+(1-2*lamdar(t))*(log(lamdar(t))-log(1-lamdar(t)))))*((1/pr(t))+vw);
Cr(t)=(1-POAr(t))*log(1-POAr(t))+POAr(t)*log(POAr(t));
Dr(t)=-POAr(t)*log(Po/(1-Po))-(1-POAr(t))*0;

Gl(t)=C(t)*(Al(t)+Bl(t)+Cl(t))+Dl(t);
Gr(t)=C(t)*(Ar(t)+Br(t)+Cr(t))+Dr(t);

Precisionl(t)=pl(t-1)/(fl(t)^2)./((1-fl(t))^2);
Precisionr(t)=pr(t-1)/(fr(t)^2)./((1-fr(t))^2); 

end

Times_of_al_tmp=sum(a);
Times_of_ar_tmp=T-Times_of_al_tmp;

Rate_of_al_tmp(i)=Times_of_al_tmp/T;
Rate_of_ar_tmp(i)=Times_of_ar_tmp/T;

end

Rate_of_al_ave=mean(Rate_of_al_tmp);
Rate_of_ar_ave=mean(Rate_of_ar_tmp)
Selections=[Rate_of_al_ave,Rate_of_ar_ave];

SD_al=sqrt(var(Rate_of_al_tmp));
SD_ar=sqrt(var(Rate_of_ar_tmp));
SD=[SD_al,SD_ar];


name=categorical({'Left';'Right'});
subplot(4,2,8)
b=bar(name,Selections,'stack')
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];
ylim([0,1])
yticks([0,0.5,1])
ylabel('Ratio','Fontsize',20,'Fontweight','bold')
hold on

er = errorbar(1:2,Selections,SD,SD);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off











