clear all
clc
% Simulation of the geodetic network and calculation of Mean Succes Rates
% datas: network, network2a, network2b,
% model: fsme (forward search of model error)

aralik=3; % magnitude of outlier
delta=6-3; % magnitude of outlier

format long g
aradongu=0;
sayac=0;
dongu=0;
saybesdurumu=0;
saydortdurumu=0;
sayteklidurumu=0;
sayucdurumu=0;
sayikidurumu=0;
msr1=0;
msr2=0;
msr3=0;
msr4=0;
tii1 = tic ;

for basla=1:1  % number of outliers
    % top(basla)=0;
    birim=1000; % unit mm
    sigma0=1; % a priori standard deviation
    m=100;  % producing standard normal distribution
    mm=100; % producing outlier
    u=7;    % number of unknown

    network2b; % data: coefficient, heights, leveling line and weights for the network
    sayac=sayac+1;
    for tt=1:m
        tt;
        [n nu]=size(Olculer);
        f=n-u+1; % degrees of freedom
        % normally distributed error to be added to the "good" observations
        randn('seed',(12349.567*tt));
        eh11=randn(100,1);
        rand('seed',(12348.567*tt));
        ara=randperm(100);
        for i=1:n
            eh1(i,1)=eh11(ara(i));% normally distributed errors
        end

        % weight matrix
        for i=1:n
            P2(i)=1/ghetero(i);
        end

        % standard deviations for each observation
        for i=1:n
            stan(i)=sqrt(ghetero(i));
        end
        eh=eh1.*sqrt(ghetero);
        Agirlik=diag(P2);
        disp(tt)
        aradongu=aradongu+1;

        for jj=1:mm
            jj;
            dongu=dongu+1;

            rand('seed',12349.89*jj*tt)
            t2=rand(1,n);

            rand('seed',12348.89*jj*tt)
            t1=rand(1,n);

            rand('seed',12347.89*jj*tt)
            I=randperm(n);


            % random selection of observations to be attributed to gross error

            %             aaa=10; % if not random selected
            %             bbb=3;
            aaa=I(1);
            bbb=I(2);
            ccc=I(3);
            ddd=I(4);
            eee=I(5);
            fff=I(6);
            ggg=I(7); %until all possible number of potential outlier

            dy=zeros(1,n); % outlier vector

            if basla==1
                uy=1;
                dy(aaa)=aralik*stan(aaa)+delta*t1(aaa)*stan(aaa);
                %                 dy(aaa)=50; % gross error
                noktalar=aaa;
            elseif basla==2
                uy=2;
                dy(aaa)=aralik*stan(aaa)+delta*t1(aaa)*stan(aaa);
                dy(bbb)=(aralik*stan(bbb)+delta*t1(bbb)*stan(bbb));
                %                   dy(aaa)=1000; % extreme outliers
                %                   dy(bbb)=1000;
                noktalar=sort([aaa bbb]);
            elseif basla==3
                uy=3;
                dy(aaa)=aralik*stan(aaa)+delta*t1(aaa)*stan(aaa);
                dy(bbb)=aralik*stan(bbb)+delta*t1(bbb)*stan(bbb);
                dy(ccc)=aralik*stan(ccc)+delta*t1(ccc)*stan(ccc);
                %                 dy(bbb)=1000;
                noktalar=sort([aaa bbb ccc]);
            elseif basla==4
                uy=4;
                dy(aaa)=aralik*stan(aaa)+delta*t1(aaa)*stan(aaa);
                dy(bbb)=aralik*stan(bbb)+delta*t1(bbb)*stan(bbb);
                dy(ccc)=aralik*stan(ccc)+delta*t1(ccc)*stan(ccc);
                dy(ddd)=aralik*stan(ddd)+delta*t1(ddd)*stan(ddd);
                noktalar=sort([aaa bbb ccc ddd]);
            elseif basla==5
                uy=5;
                dy(aaa)=aralik*stan(aaa)+delta*t1(aaa)*stan(aaa);
                dy(bbb)=aralik*stan(bbb)+delta*t1(bbb)*stan(bbb);
                dy(ccc)=aralik*stan(ccc)+delta*t1(ccc)*stan(ccc);
                dy(ddd)=aralik*stan(ddd)+delta*t1(ddd)*stan(ddd);
                dy(eee)=aralik*stan(eee)+delta*t1(eee)*stan(eee);
                noktalar=sort([aaa bbb ccc ddd eee]);
            elseif basla==6
                uy=6;
                dy(aaa)=aralik*stan(aaa)+delta*t1(aaa)*stan(aaa);
                dy(bbb)=aralik*stan(bbb)+delta*t1(bbb)*stan(bbb);
                dy(ccc)=aralik*stan(ccc)+delta*t1(ccc)*stan(ccc);
                dy(ddd)=aralik*stan(ddd)+delta*t1(ddd)*stan(ddd);
                dy(eee)=aralik*stan(eee)+delta*t1(eee)*stan(eee);
                dy(fff)=aralik*stan(fff)+delta*t1(fff)*stan(fff);
                noktalar=sort([aaa bbb ccc ddd eee fff]);
            elseif basla==7
                uy=7;
                dy(aaa)=aralik*stan(aaa)+delta*t1(aaa)*stan(aaa);
                dy(bbb)=aralik*stan(bbb)+delta*t1(bbb)*stan(bbb);
                dy(ccc)=aralik*stan(ccc)+delta*t1(ccc)*stan(ccc);
                dy(ddd)=aralik*stan(ddd)+delta*t1(ddd)*stan(ddd);
                dy(eee)=aralik*stan(eee)+delta*t1(eee)*stan(eee);
                dy(fff)=aralik*stan(fff)+delta*t1(fff)*stan(fff);
                dy(ggg)=aralik*stan(ggg)+delta*t1(ggg)*stan(ggg);
                noktalar=sort([aaa bbb ccc ddd eee fff ggg]);

            end

            % Calculation of "good" and "bad" observations
            % sign function to select randomly
            %             dy(aaa)=tt;
            for i=1:n
                if t2(i)<0.5
                    dy(i)=-dy(i);
                else dy(i)=dy(i);
                end
            end

            %                 dy(aaa)=jj*tt;

            % adding random error and outlier vectors
            for ii=1:n
                if dy(ii)~=0
                    holcu(ii,1)=Olculer(ii)+dy(ii)/birim;
                else

                    holcu(ii,1)=Olculer(ii)+eh(ii)/birim;
                end
            end

            for ii=1:n
                holcu2(ii,1)=Olculer(ii)+eh(ii)/birim;
            end

            % calculating approximate heights
            H1(1)=H(1);
            H1(2)=H1(1)-holcu(3);
            H1(3)=H1(2)-holcu(6);
            H1(5)=H1(1)-holcu(1);
            H1(4)=H1(5)+holcu(11);
            H1(6)=H1(1)-holcu(2);
            H1(7)=H1(2)+holcu(5);

            % normal equations
            araolcu=Coeff*H;
            l1=(holcu-araolcu)*birim;
            A=Coeff;
            P=Agirlik;
            N=A'*P*A;
            x1=pinv(N)*A'*P*l1;

            %residuals
            v1=A*x1-l1;

            % % control check
            % denetim21=l1'*P*l1-(A'*P*l1)'*x1;
            % denetim22=v1'*P*v1;
            %
            %
            % denetim11=[];
            % denetim11=(v1/1000)+holcu;
            % HH=[];
            % for i=1:uu
            %     HH(i)=H(i)+(x1(i)/1000);
            % end
            % HH=HH';
            % holcuu=[];
            %
            %
            % holcuu(1)=HH(1)-HH(5);
            % holcuu(2)=HH(1)-HH(6);
            % holcuu(3)=HH(1)-HH(2);
            % holcuu(4)=HH(6)-HH(2);
            % holcuu(5)=HH(7)-HH(2);
            % holcuu(6)=HH(2)-HH(3);
            % holcuu(7)=HH(7)-HH(3);
            % holcuu(8)=HH(4)-HH(3);
            % holcuu(9)=HH(4)-HH(7);
            % holcuu(10)=HH(4)-HH(6);
            % holcuu(11)=HH(4)-HH(5);
            % holcuu(12)=HH(6)-HH(5);
            % holcuu(13)=HH(7)-HH(6);
            % holcuu(14)=HH(7)-HH(1);
            % holcuu(15)=HH(7)-HH(5);

            fsme  % Forward Search of Model Error

            % calculation of mean success rates
            if durum2==1
                teklikom=uyusumsuz1;
                gercektekli=[aaa];
                if artim2==1
                    if gercektekli==teklikom
                        sayteklidurumu=sayteklidurumu+1;
                    end
                end
            end

            if durum2==1
                if durum==1
                    ikilikom=uyusumsuz2;
                    gercekikili=[aaa bbb];
                    gercekikili=sort(gercekikili);
                    if artim==1
                        if gercekikili==ikilikom
                            sayikidurumu=sayikidurumu+1;
                        end
                    end
                end
            end

            if durum2==1
                if durum==1
                    if durum3==1
                        uclukom=uyusumsuz3;
                        gercekuclu=[aaa bbb ccc];
                        gercekuclu=sort(gercekuclu);
                        if artim3==1
                            if gercekuclu==uclukom
                                sayucdurumu=sayucdurumu+1;
                            end
                        end
                    end
                end
            end


            if durum2==1
                if durum==1
                    if durum3==1
                        if durum4==1
                            dortlukom =sort(uyusumsuz4);
                            gercekdortlu=[aaa bbb ccc ddd];
                            gercekdortlu=sort(gercekdortlu);
                            if artim4==1
                                if gercekdortlu==dortlukom
                                    saydortdurumu=saydortdurumu+1;
                                end
                            end
                        end
                    end
                end
            end

            %
            % if durum2==1
            % if durum==1
            % if durum3==1
            % if durum4==1
            % if durum5==1
            % %                     beslikom=[besinci1 besinci2 besinci3 besinci4 besinci5];      % comb
            %
            % %                     beslikom=sort(beslikom);   % comb
            %     beslikom=sort(uyusumsuz5);  % model error
            %     gercekbesli=[aaa bbb ccc ddd eee];
            %     gercekbesli=sort(gercekbesli);
            %     if artim5==1
            %         if gercekbesli==beslikom
            %             saybesdurumu=saybesdurumu+1;
            %         end
            %     end
            % end
            % end
            % end
            % end
            % end


            %             elseif durum3==3

            %             elseif durum==1

            % end
            %             elseif durum2==1


            %             uyusumsuztop=[];
            %             uyusumsuztop=[uyusumsuz1 uyusumsuz2 uyusumsuz3 uyusumsuz4];
            %             uyusumsuztop=sort(uyusumsuztop);
            %              [~,~,ix] = unique(uyusumsuztop);
            %              CC = accumarray(ix,1).'


            %             [aaa bbb]
            %             gercekikili
            %             pause

        end
    end
end


msr1=(sayteklidurumu*100)/(m*mm) % msr for one outlier when basla=1
msr2=(sayikidurumu*100)/(m*mm) % msr for two outliers when basla=2
msr3=(sayucdurumu*100)/(m*mm)
msr4=(saydortdurumu*100)/(m*mm)
msr5=(saybesdurumu*100)/(m*mm)




tii1 = toc(tii1)
