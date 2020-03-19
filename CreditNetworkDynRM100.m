% CreditNetworkDyn

clear all
warning off

T=1000; %input('Number of simulation periods=   ');
Nd=500; %input('Number of downstream firms=   ');
Nu=250; %input('Number of upstream firms=   ');
Nb=100; %input('Number of banks=   ');
MC=100;

rand('state',15);

for mc=1:MC;
    
    disp(mc)

    %% parameters
    phi=1.5;
    beta=0.8;
    gammad=0.5;
    delta=0.5;
    gammau=1;
    par1(1:Nb)=0.05;
    par2=0.05;
    w=1;
    pmin=0;
    pmax=2-pmin;

    %     if mc==1
    %         epslon=1;
    %     elseif mc==2
    %         epslon=0.01;
    %     end

    %%%%%%%%%%%%%epslon=0.01; % --> "Prefered-Partner Choice"
    epslon=1; % --> "Random Matching"

    %% initial conditions: individual variables
    Ru(1:Nu)=0.1;
    Rbd(1:Nd)=0;
    Rbu(1:Nu)=0;
    Ad(1:Nd)=1;
    Au(1:Nu)=1;
    Ab(1:Nb)=1;
    Yu2(1:Nu)=0;
    Yd(1:Nd)=0;
    Ld(1:Nd)=0;
    Wd(1:Nd)=0;
    Prd(1:Nd)=0;
    Bd(1:Nd)=0;
    pd(1:Nd)=1;
    Yu(1:Nu)=0;
    Lu(1:Nu)=0;
    Wu(1:Nu)=0;
    Pru(1:Nu)=0;
    Bu(1:Nu)=0;
    pu(1:Nu)=1+Ru;
    Badu(1:Nu)=0;
    Badb(1:Nb)=0;
    Prb(1:Nb)=0;
    falld(1:Nd)=0;
    fallu(1:Nu)=0;
    fallb(1:Nb)=0;
    link_db(1:Nd)=ceil(rand(1,Nd)*Nb);
    link_ub(1:Nu)=ceil(rand(1,Nu)*Nb);
    link_du(1:Nd)=ceil(rand(1,Nd)*Nu);

    %% initial conditions: aggregate variables
    YD(1,mc)=sum(Yd);
    YU(1,mc)=sum(Yu);
    AB(1,mc)=sum(Ab);
    BAD(1,mc)=0;
    FALLD(1,mc)=0;
    FALLU(1,mc)=0;
    ERD(1,mc)=1;
    ERU(1,mc)=1;
    %%%LR(1,mc)=0;

    tent=5;

    for j=1:Nu
        linkU(j,mc)=length(find(link_du==j));
    end
    for j=1:Nb
        linkB(j,mc)=length(find(link_db==j))+length(find(link_ub==j));
    end


    %% MAIN PROGRAM
    for t=2:T

        %disp(t)

        %%%%% entry-exit %%%%%
        Ad=Ad+Prd; Ad(falld==1)=2*rand(length(find(falld==1)),1);
        Au=Au+Pru; Au(fallu==1)=2*rand(length(find(fallu==1)),1);
        Ab=Ab-Badb+Prb; Ab(fallb==1)=2*rand(length(find(fallb==1)),1);

%         %% SHOCK THE ECONOMY %%%%%%%%%%%%
%         xx=max(linkU(:,mc)); xxx=find(linkU(:,mc)==xx); %% TARGETED SHOCK
%         %xxx=ceil(Nu*rand(1,1)); %% RANDOM SHOCK
%         Au(find(linkU(:,mc)==xxx(1)))=0.001; %Au(find(linkU(:,mc)==xxx(1)))*0.01;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     numstoc=rand(1,1); shock=0;
        %     if numstoc<=1/3 & shock==1;
        %         xx=max(Ad); xxx=find(Ad==xx);
        %         %xxx=ceil(Nd*rand(1,1));
        %         Ad(find(Ad==xxx))=Ad(find(Ad==xxx))*0.01;
        %     elseif numstoc>1/3 & numstoc<=2/3 & shock==1
        %         xx=max(Au); xxx=find(Au==xx);
        %         %xxx=ceil(Nu*rand(1,1));
        %         Ad(find(Au==xxx))=Au(find(Au==xxx))*0.01;
        %     elseif numstoc>2/3 & shock==1;
        %         xx=max(Ab); xxx=find(Ab==xx);
        %         %xxx=ceil(Nb*rand(1,1));
        %         Ab(find(Ab==xxx))=Ab(find(Ab==xxx))*0.01;
        %     end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %     numstoc=rand(1,1); shock=1;
        %     if numstoc<=1/2 & shock==1;
        %         xx=max(linkU); xxx=find(linkU==xx);
        %         %xxx=ceil(Nu*rand(1,1));
        %         Au(find(linkU==xxx(1)))=Au(find(linkU==xxx(1)))*0.01;
        %     elseif numstoc>1/2 & shock==1
        %         xx=max(linkB); xxx=find(linkB==xx);
        %         %xxx=ceil(Nb*rand(1,1));
        %         Ab(find(linkB==xxx(1)))=Ab(find(linkB==xxx(1)))*0.01;
        %     end

        Add=Ad; Auu=Au; Abb=Ab;

        %%%%% downstream firms %%%%%%
        Yd=phi.*Ad.^beta; Qd=delta*Yd; Ld=gammad*Yd; Wd=w*Ld;
        Bd(Wd-Ad>0)=Wd(Wd-Ad>0)-Ad(Wd-Ad>0); Bd(Wd-Ad<=0)=0;
        Rbd=par1(link_db)+par2*(Bd./Ad).^par2;
        Rbd(Rbd<0.01)=0.01;
        Ad(Wd-Ad>0)=0; Ad(Wd-Ad<=0)=Ad(Wd-Ad<=0)-Wd(Wd-Ad<=0);
        pd=(pmin+(pmax-pmin)*rand(1,Nd));
        Prd=pd.*Yd-pu(link_du).*Qd-(1+Rbd).*Bd;
        falld(Ad+Prd>0)=0; falld(Ad+Prd<=0)=1;

        %%%%%% upstream firms %%%%%%
        Ru=0.1.*Au.^-0.1; pu=1+Ru;
        Yu(1:Nu)=0;
        Yu2(1:Nu)=0;
        Badu(1:Nu)=0;
        for j=1:Nd
            ind=link_du(j);
            Yu(ind)=Yu(ind)+Qd(j);
            if falld(j)==0
                Yu2(ind)=Yu(ind);
                Badu(ind)=0;
            else
                Yu2(ind)=0;
                Badu(ind)=Badu(ind)+Qd(j);
            end
        end
        Lu=gammau*Yu; Wu=w*Lu;
        Bu(Wu-Au>0)=Wu(Wu-Au>0)-Au(Wu-Au>0); Bu(Wu-Au<=0)=0;
        Rbu=par1(link_ub)+par2*(Bu./Au).^par2;
        Rbu(Rbu<0.01)=0.01;
        Au(Wu-Au>0)=0; Au(Wu-Au<=0)=Au(Wu-Au<=0)-Wu(Wu-Au<=0);
        Pru=pu.*(1+Ru).*Yu2-(1+Rbu).*Bu;
        fallu(Au+Pru>0)=0; fallu(Au+Pru<=0)=1;

        %%%%% aggregate variables %%%%%
        YD(t,mc)=sum(Yd); YU(t,mc)=sum(Yu);
        FALLD(t,mc)=length(find(falld==1)); FALLU(t,mc)=length(find(fallu==1));

        %%% preferential partner choice  %%%%%
        for j=1:Nd
            %%% ...downstream vs. upstream firms %%%%
            if rand(1,1) < epslon
                link_du(j) = ceil(rand(1,1)*Nu);
            else
                newdu=ceil(rand(tent,1)*Nu);
                if min(pu(newdu))<pu(link_du(j))
                    sfacc=find(pu(newdu)==min(pu(newdu)));
                    link_du(j)=newdu(sfacc(1));
                else
                    link_du(j)=link_du(j);
                end
                if falld(j)==1
                    link_du(j)=ceil(rand(1,1)*Nu);
                end
                while fallu(link_du(j))==1
                    link_du(j)=ceil(rand(1,1)*Nu);
                end
            end
            %%% ...downstream vs. banks %%%%
            if rand(1,1) < epslon
                link_db(j)=ceil(rand(1,1)*Nb);
            else
                newdb=ceil(rand(tent,1)*Nb);
                if min(par1(newdb)) < par1(link_db(j))
                    sfacc2=find(par1(newdb)==min(par1(newdb)));
                    link_db(j)=newdb(sfacc2(1));
                else
                    link_db(j)=link_db(j);
                end
                if falld(j)==1
                    link_db(j)=ceil(rand(1,1)*Nb);
                end
            end
        end
        %%% ... upstream firms vs. banks %%%%%
        for j=1:Nu
            if rand(1,1) < epslon
                link_ub(j)=ceil(rand(1,1)*Nb);
            else
                newub=ceil(rand(tent,1)*Nb);
                if min(par1(newub)) < par1(link_ub(j))
                    sfacc3=find(par1(newub)==min(par1(newub)));
                    link_ub(j)=newub(sfacc3(1));
                else
                    link_ub(j)=link_ub(j);
                end
                if fallu(j)==1
                    link_ub(j)=ceil(rand(1,1)*Nb);
                end
            end
        end

        %%%% banks %%%%%
        Prb(1:Nb)=0;
        for j=1:Nb
            Badb(j)=sum(Bd(falld==1&link_db==j))+sum(Bu(fallu==1&link_ub==j));
            Prb(j)=Bd(link_db==j)*Rbd(link_db==j)'+Bu(link_ub==j)*Rbu(link_ub==j)';
        end
        fallb(Ab-Badb+Prb>0)=0; fallb(Ab-Badb+Prb<=0)=1;
        par1=0.1*Ab.^-0.1;

        %% random connections for bankrupted agents
        for j=1:Nd
            while fallb(link_db(j))==1 & sum(fallb)<Nb
                link_db(j)=ceil(rand(1,1)*Nb);
            end
        end
        for j=1:Nu
            while fallb(link_ub(j))==1 & sum(fallb)<Nb
                link_ub(j)=ceil(rand(1,1)*Nb);
            end
        end

        %%% aggregate variables
        AB(t,mc)=sum(Ab); AD(t,mc)=sum(Ad); AU(t,mc)=sum(Au); BB(t,mc)=sum(Bd)+sum(Bu);
        ABB(t,mc)=sum(Abb); ADD(t,mc)=sum(Add); AUU(t,mc)=sum(Auu);
        RBD(t,mc)=mean(Rbd); RBU(t,mc)=mean(Rbu); RU(t,mc)=mean(Ru);
        BAD(t,mc)=sum(Badb)+sum(Badu); FALLB(t,mc)=length(find(fallb==1));
        PRU(t,mc)=mean(Pru); PRD(t,mc)=mean(Prd); PRB(t,mc)=mean(Prb);
        num=0;
        for j=1:Nd
            num=num+(Add(j)/(Add(j)+Bd(j)))*Add(j);
        end
        ERD(t,mc)=num/sum(Add);
        num=0;
        for j=1:Nu
            num=num+(Auu(j)/(Auu(j)+Bu(j)))*Auu(j);
        end
        ERU(t,mc)=num/sum(Auu);
        %%%ERD(t,mc)=sum(Add)/(sum(Add)+sum(Bd)); ERU(t,mc)=sum(Auu)/(sum(Auu)+sum(Bu));
        %%%LRD(t,mc)=sum(Bd)/sum(Add); LRU(t,mc)=sum(Bu)/sum(Auu);
        
        %%% degree distribution(s)
        for j=1:Nu
            linkU(j,mc)=length(find(link_du==j));
        end
        for j=1:Nb
            linkB(j,mc)=length(find(link_db==j))+length(find(link_ub==j));
        end

        %     figure(01)
        %     subplot(2,2,1)
        %     semilogy(YD)
        %     title('aggregate production')
        %     subplot(2,2,2)
        %     loglog(abs(sort(-[linkU linkB])),1:length([linkU linkB]),'b.')
        %     title('degree distribution of the network')
        %     subplot(2,2,3)
        %     loglog(abs(sort(-[Ad Au])),1:length([Ad Au]),'.')
        %     title('firm size distribution')
        %     subplot(2,2,4)
        %     plot(BAD)
        %     title('bad debt')
        %     pause(0.01)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         BDAb=Badb./Ab; BDAu=Badu./Au;
%         ProbBDAu(t,mc)=length(find(abs(BDAu-median(BDAu))>=4*std(BDAu)))/length(BDAu);
%         ProbBDAb(t,mc)=length(find(abs(BDAb-median(BDAb))>=4*std(BDAb)))/length(BDAb);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %% extreme event statistics (MC simulations) %%%%
    
%     % bad debt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     BADmean(mc)=mean(BAD); BADvar(mc)=var(BAD); BADskew(mc)=skewness(BAD); BADkurt(mc)=kurtosis(BAD); 
%     BADmedian(mc)=median(BAD); BADstd(mc)=std(BAD); BADiqr(mc)=iqr(BAD);
%     VEC=[]; VEC=abs(BAD-BADmedian(mc));
%     
%     ww=1;
%     stdstep(ww)=0;
%     VECEE=[]; VECEE=find(VEC>=stdstep(ww)*BADstd(mc)); 
%     ProbBAD(1,mc)=length(VECEE)/length(VEC);
%     while length(VECEE)>0
%         ww=ww+1;
%         stdstep(ww)=stdstep(ww-1)+0.1;
%         VECEE=find(VEC>=stdstep(ww)*BADstd(mc));
%         ProbBAD(ww,mc)=length(VECEE)/length(VEC);
%     end
% 
%     % bad debt to net worth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     BDAmean(mc)=mean(BAD./(AB+AU));
%     BDAvar(mc)=var(BAD./(AB+AU));
%     BDAskew(mc)=skewness(BAD./(AB+AU));
%     BDAkurt(mc)=kurtosis(BAD./(AB+AU)); 
%     
%     BDAmedian(mc)=median(BAD./(AB+AU));
%     BDAstd(mc)=std(BAD./(AB+AU));
%     BDAiqr(mc)=iqr(BAD./(AB+AU));
%     VEC=[]; VEC=abs(BAD./(AB+AU)-BDAmedian(mc));
%     
%     ww=1;
%     stdstep(ww)=0;
%     VECEE=[]; VECEE=find(VEC>=stdstep(ww)*BDAstd(mc)); 
%     ProbBDA(1,mc)=length(VECEE)/length(VEC);
%     while length(VECEE)>0
%         ww=ww+1;
%         stdstep(ww)=stdstep(ww-1)+0.1;
%         VECEE=find(VEC>=stdstep(ww)*BDAstd(mc));
%         ProbBDA(ww,mc)=length(VECEE)/length(VEC);
%     end
% 
%     % growth rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     gr=log(YD(3:T)./YD(2:T-1));
%     GRmean(mc)=mean(gr); GRvar(mc)=var(gr); GRskew(mc)=skewness(gr); GRkurt(mc)=kurtosis(gr);
%     GRmedian(mc)=median(gr); GRstd(mc)=std(gr); GRiqr(mc)=iqr(gr); 
%     VECC=[]; VECC=abs(gr-GRmedian(mc));
%     
%     ww=1;
%     stdstep1(ww)=0;
%     VECCEE=[]; VECCEE=find(VECC>=stdstep1(ww)*GRstd(mc));
%     ProbGR(1,mc)=length(VECCEE)/length(VECC);
%     while length(VECCEE)>0
%         ww=ww+1;
%         stdstep1(ww)=stdstep1(ww-1)+0.1;
%         VECCEE=find(VECC>=stdstep1(ww)*GRstd(mc));
%         ProbGR(ww,mc)=length(VECCEE)/length(VECC);
%     end
    
    % bankruptcy probability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PB(mc)=(sum(FALLD(1:T,mc))+sum(FALLU(1:T,mc))+sum(FALLB(1:T,mc)))/((Nd+Nu+Nb)*T);
    PBD(mc)=sum(FALLD(1:T,mc))/(Nd*T);
    PBU(mc)=sum(FALLU(1:T,mc))/(Nu*T);
    PBB(mc)=sum(FALLB(1:T,mc))/(Nb*T);
    
    % FSD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AAd(:,mc)=Ad; AAu(:,mc)=Au; AAb(:,mc)=Ab;
    YYd(:,mc)=Yd; YYu(:,mc)=Yu;
end

GR(2:T,:)=log(YD(2:T,:)./YD(1:T-1,:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Robust measure of SKEWNESS
for j=1:MC
    div=[]; div=mean(abs(GR(3:T,j)-median(GR(3:T,j))));
    SK(j)=(mean(GR(3:T,j))-median(GR(3:T,j)))/div;
end

% Robust measure of KURTOSIS
for j=1:MC
    num=0;
    for k=1:8
        num=num+12.5;
        E(k,j)=prctile(GR(3:T,j),num);
    end
end
for j=1:MC
    KR(j)=(E(7,j)-E(5,j)+E(3,j)-E(1,j))/(E(6,j)-E(2,j));
end

% Bankruptcy correlation across sectors
for j=1:MC
    vec=corrcoef(FALLD(:,j),FALLU(:,j)); CORR_D_U(j)=vec(1,2);
    vec=corrcoef(FALLD(:,j),FALLB(:,j)); CORR_D_B(j)=vec(1,2);
    vec=corrcoef(FALLU(:,j),FALLB(:,j)); CORR_U_B(j)=vec(1,2);
    vec=corrcoef(FALLD(:,j)+FALLU(:,j),FALLB(:,j)); CORR_DU_B(j)=vec(1,2);
end

% Correlation between aggregate debt-to-equity ratio and business cycle
for j=1:MC
    vec=corrcoef(BB(251:T,j)./(AD(251:T,j)+AU(251:T,j)),YD(251:T,j)); CCC(j)=vec(1,2);
end

save RM100
