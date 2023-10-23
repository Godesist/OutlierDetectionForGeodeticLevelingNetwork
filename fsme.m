
% forward search of model error

alfa=0.05;
s0=1;
%% one model error
ek1=1;
nm2=n;
PP=P;
LL=l1;
A1=A;
nsay=u;
durum2=0;
top=0;
katsayiy=[];
uz=nsay+ek1;
artim2=0;
Ny=[];
Qxxy=[];
Bilinmeyenlery=[];
Duzeltmelery=[];
dvaryans=[];
Bil=[];
R1=[];
R2=[];
s00=[];
tta1=[];
% ttb=[];
uyusumsuz1=[];

for j1=1:nm2
    top=top+1;
    for ii=1:nm2
        for jj=1:u
            katsayiy(ii,jj)=A1(ii,jj);
        end
    end
    for ii=1:nm2
        if j1==ii
            katsayiy(ii,uz)=1;
        else
            katsayiy(ii,uz)=0;
        end
    end
    Ny=katsayiy'*PP*katsayiy;
    %Qxxy=pinv(Ny+GT'*GT)-GT'*inv(GT*GT'*GT*GT')*GT;
    Qxxy=pinv(Ny,0.001);

    Bilinmeyenlery=Qxxy*katsayiy'*PP*LL;
    Bil(j1)=Bilinmeyenlery(uz);
    Duzeltmelery=katsayiy*Bilinmeyenlery-LL;
    dvaryans(j1)=(Duzeltmelery'*PP*Duzeltmelery)/f;
end
[R1,R2]=min(dvaryans');
f1=n-u-1+1;
s00=sqrt(R1);
tta1=Bil(R2)/(s00*sqrt(Qxxy(8,8)));

if abs(tta1)>tinv(0.95,f1)
    uyusumsuz1=R2;
    durum2=1;
    artim2=1;
end

%% two model errors
me=2;
ek1=1;
nm2=n;
PP=P;
LL=l1;
A1=A;
nsay=u;
durum=0;
top2=0;
katsayit=[];
uz=nsay+ek1;
Nt=[];
Qxxt=[];
Bilinmeyenlert=[];
Duzeltmelert=[];
dvaryans1=[];
var1=[];
var2=[];
bil1=[];
bil2=[];
R3=[];
R4=[];
s00=[];
tta2=[];
ttb2=[];
uyusumsuz11=[];
uyusumsuz12=[];
uyusumsuz2=[];
artim=0;

for k1=1:nm2-1
    for k2=k1+1:nm2
        top2=top2+1;
        for ii=1:nm2
            for jj=1:nsay
                katsayit(ii,jj)=A1(ii,jj);
            end
        end
        for ii=1:nm2
            if k1==ii
                katsayit(ii,uz)=1;
            else
                katsayit(ii,uz)=0;
            end
            if k2==ii
                katsayit(ii,uz+1)=1;
            else
                katsayit(ii,uz+1)=0;
            end
        end
        Nt=katsayit'*PP*katsayit;
        %Qxx=inv(Nt+GT'*GT)-GT'*inv(GT*GT'*GT*GT')*GT;
        Qxxt=pinv(Nt,0.001);
        Bilinmeyenlert=Qxxt*katsayit'*PP*LL;
        Duzeltmelert=katsayit*Bilinmeyenlert-LL;
        dvaryans1(top2)=(Duzeltmelert'*PP*Duzeltmelert)/f;
        var1(top2)=k1;
        var2(top2)=k2;
        bil1(top2)=Bilinmeyenlert(uz);
        bil2(top2)=Bilinmeyenlert(uz+1);
    end
end
dvaryans1';
[R3,R4]=min(dvaryans1');
var1(R4);
var2(R4);
s00=sqrt(R3);
tta2=bil1(R4)/(s00*sqrt(Qxxt(8,8)));
ttb2=bil2(R4)/(s00*sqrt(Qxxt(9,9)));
f2=nm2-nsay-me+1;
if abs(tta2)>tinv(0.95,f2) && abs(ttb2)>tinv(0.95,f2)
    uyusumsuz11=var1(R4);
    uyusumsuz12=var2(R4);
    uyusumsuz2=[uyusumsuz11 uyusumsuz12];
    uyusumsuz2=sort(uyusumsuz2);
    durum=1;
    artim=1;
end

%% three model errors
me=3;
ek1=1;
nm2=n;
PP=P;
LL=l1;
A1=A;
nsay=u;
durum3=0;
top3=0;
katsayiz=[];
uz=nsay+ek1;
Nz=[];
Qxxz=[];
Bilinmeyenlerz=[];
Duzeltmelerz=[];
dvaryans3=[];
bil3=[];
bil4=[];
bil5=[];
R6=[];
R5=[];
var1=[];
var2=[];
var3=[];
s00=[];
tta3=[];
ttb3=[];
ttc=[];
% ttb=[];
uyusumsuz3=[];
uyusumsuz31=[];
uyusumsuz31=[];
uyusumsuz33=[];
artim3=0;

for t1=1:nm2-2
    for t2=t1+1:nm2-1
        for t3=t2+1:nm2
            top3=top3+1;
            for ii=1:nm2
                for jj=1:nsay
                    katsayiz(ii,jj)=A1(ii,jj);
                end
            end
            for ii=1:nm2
                if t1==ii
                    katsayiz(ii,uz)=1;
                else
                    katsayiz(ii,uz)=0;
                end
                if t2==ii
                    katsayiz(ii,uz+1)=1;
                else
                    katsayiz(ii,uz+1)=0;
                end
                if t3==ii
                    katsayiz(ii,uz+2)=1;
                else
                    katsayiz(ii,uz+2)=0;
                end
            end

            Nz=katsayiz'*PP*katsayiz;
            Qxxz=pinv(Nz);
            % Qxxt=inv(Nz+GT'*GT)-GT'*inv(GT*GT'*GT*GT')*GT;
            Bilinmeyenlerz=Qxxz*katsayiz'*PP*LL;
            Duzeltmelerz=katsayiz*Bilinmeyenlerz-LL;
            dvaryans3(top3)=(Duzeltmelerz'*PP*Duzeltmelerz)/f;
            var1(top3)=t1;
            var2(top3)=t2;
            var3(top3)=t3;
            bil3(top3)=Bilinmeyenlerz(uz);
            bil4(top3)=Bilinmeyenlerz(uz+1);
            bil5(top3)=Bilinmeyenlerz(uz+2);
        end
    end
end
top3;
[R5,R6]=min(dvaryans3');
var1(R6);
var2(R6);
var3(R6);
s00=sqrt(R5);
tta3=bil3(R6)/(s00*sqrt(Qxxz(8,8)));
ttb3=bil4(R6)/(s00*sqrt(Qxxz(9,9)));
ttc3=bil5(R6)/(s00*sqrt(Qxxz(10,10)));,
f3=n-u-3+1;
if abs(tta3)>tinv(0.95,f3) && abs(ttb3)>tinv(0.95,f3) && abs(ttc3)>tinv(0.95,f3)
    uyusumsuz31=var1(R6);
    uyusumsuz32=var2(R6);
    uyusumsuz33=var3(R6);
    uyusumsuz3=[uyusumsuz31 uyusumsuz32 uyusumsuz33];
    uyusumsuz3=sort(uyusumsuz3);
    uyusumsuz3;
    durum3=1;
    artim3=1;
end

%% four model errors
me=4;
ek1=1;
nm2=n;
PP=P;
LL=l1;
A1=A;
nsay=u;
durum4=0;
top4=0;
katsayiw=[];
uz=nsay+ek1;
Nw=[];
Qxxw=[];
Bilinmeyenlerw=[];
Duzeltmelerw=[];
dvaryans4=[];
bil6=[];
bil7=[];
bil8=[];
bil9=[];
R7=[];
R8=[];
var1=[];
var2=[];
var3=[];
var4=[];
s00=[];
tta4=[];
ttb4=[];
ttc4=[];
ttd4=[];
uyusumsuz4=[];
uyusumsuz41=[];
uyusumsuz42=[];
uyusumsuz43=[];
uyusumsuz44=[];
artim4=0;

for t1=1:nm2-3
    for t2=t1+1:nm2-2
        for t3=t2+1:nm2-1
            for t4=t3+1:nm2
                top4=top4+1;
                for ii=1:nm2
                    for jj=1:nsay
                        katsayiw(ii,jj)=A1(ii,jj);
                    end
                end
                for ii=1:nm2
                    if t1==ii
                        katsayiw(ii,uz)=1;
                    else
                        katsayiw(ii,uz)=0;
                    end
                    if t2==ii
                        katsayiw(ii,uz+1)=1;
                    else
                        katsayiw(ii,uz+1)=0;
                    end
                    if t3==ii
                        katsayiw(ii,uz+2)=1;
                    else
                        katsayiw(ii,uz+2)=0;
                    end
                    if t4==ii
                        katsayiw(ii,uz+3)=1;
                    else
                        katsayiw(ii,uz+3)=0;
                    end
                end
                Nw=katsayiw'*PP*katsayiw;
                %Qxx=inv(N+GT'*GT)-GT'*inv(GT*GT'*GT*GT')*GT;
                Qxxw=pinv(Nw,0.001);
                Bilinmeyenlerw=Qxxw*katsayiw'*PP*LL;
                Duzeltmelerw=katsayiw*Bilinmeyenlerw-LL;
                dvaryans4(top4)=(Duzeltmelerw'*PP*Duzeltmelerw)/f;
                var1(top4)=t1;
                var2(top4)=t2;
                var3(top4)=t3;
                var4(top4)=t4;
                bil6(top4)=Bilinmeyenlerw(uz);
                bil7(top4)=Bilinmeyenlerw(uz+1);
                bil8(top4)=Bilinmeyenlerw(uz+2);
                bil9(top4)=Bilinmeyenlerw(uz+3);
            end
        end
    end
end
top4;
[R7,R8]=min(dvaryans4');
var1(R8);
var2(R8);
var3(R8);
var4(R8);
s00=sqrt(R7);
tta4=bil6(R8)/(s00*sqrt(Qxxw(8,8)));
ttb4=bil7(R8)/(s00*sqrt(Qxxw(9,9)));
ttc4=bil8(R8)/(s00*sqrt(Qxxw(10,10)));
ttd4=bil9(R8)/(s00*sqrt(Qxxw(11,11)));
f4=n-u-me+1;

if abs(tta4)>tinv(0.95,f4) && abs(ttb4)>tinv(0.95,f4) && abs(ttc4)>tinv(0.95,f4) && abs(ttd4)>tinv(0.95,f4)
    uyusumsuz41=var1(R8);
    uyusumsuz42=var2(R8);
    uyusumsuz43=var3(R8);
    uyusumsuz44=var4(R8);
    uyusumsuz4=[uyusumsuz41 uyusumsuz42 uyusumsuz43 uyusumsuz44];
    uyusumsuz4=sort(uyusumsuz4);
    uyusumsuz4;
    durum4=1;
    artim4=1;
end

%% five model errors
me=5;
ek1=1;
nm2=n;
PP=P;
LL=l1;
A1=A;
nsay=u;
durum5=0;
top5=0;
katsayiw=[];
uz=nsay+ek1;
Np=[];
Qxxp=[];
Bilinmeyenlerp=[];
Duzeltmelerp=[];
dvaryans5=[];
bil6=[];
bil7=[];
bil8=[];
bil9=[];
bil10=[];
R9=[];
R10=[];
var1=[];
var2=[];
var3=[];
var4=[];
var5=[];
s00=[];
tta5=[];
ttb5=[];
ttc5=[];
ttd5=[];
tte5=[];
uyusumsuz5=[];
uyusumsuz51=[];
uyusumsuz52=[];
uyusumsuz53=[];
uyusumsuz54=[];
uyusumsuz55=[];
artim5=0;

for t1=1:nm2-4
    for t2=t1+1:nm2-3
        for t3=t2+1:nm2-2
            for t4=t3+1:nm2-1
                for t5=t4+1:nm2
                    top5=top5+1;
                    for ii=1:nm2
                        for jj=1:nsay
                            katsayip(ii,jj)=A1(ii,jj);
                        end
                    end
                    for ii=1:nm2
                        if t1==ii
                            katsayip(ii,uz)=1;
                        else
                            katsayip(ii,uz)=0;
                        end
                        if t2==ii
                            katsayip(ii,uz+1)=1;
                        else
                            katsayip(ii,uz+1)=0;
                        end
                        if t3==ii
                            katsayip(ii,uz+2)=1;
                        else
                            katsayip(ii,uz+2)=0;
                        end
                        if t4==ii
                            katsayip(ii,uz+3)=1;
                        else
                            katsayip(ii,uz+3)=0;
                        end
                        if t5==ii
                            katsayip(ii,uz+4)=1;
                        else
                            katsayip(ii,uz+4)=0;
                        end
                    end
                    Np=katsayip'*PP*katsayip;
                    %Qxx=inv(N+GT'*GT)-GT'*inv(GT*GT'*GT*GT')*GT;
                    Qxxp=pinv(Np,0.001);
                    Bilinmeyenlerp=Qxxp*katsayip'*PP*LL;
                    Duzeltmelerp=katsayip*Bilinmeyenlerp-LL;
                    dvaryans5(top5)=(Duzeltmelerp'*PP*Duzeltmelerp)/f;
                    var1(top5)=t1;
                    var2(top5)=t2;
                    var3(top5)=t3;
                    var4(top5)=t4;
                    var5(top5)=t5;
                    bil6(top5)=Bilinmeyenlerp(uz);
                    bil7(top5)=Bilinmeyenlerp(uz+1);
                    bil8(top5)=Bilinmeyenlerp(uz+2);
                    bil9(top5)=Bilinmeyenlerp(uz+3);
                    bil10(top5)=Bilinmeyenlerp(uz+4);
                end
            end
        end
    end
end
top5;
[R9,R10]=min(dvaryans5');
var1(R10);
var2(R10);
var3(R10);
var4(R10);
var5(R10);
s00=sqrt(R9);
tta5=bil6(R10)/(s00*sqrt(Qxxp(8,8)));
ttb5=bil7(R10)/(s00*sqrt(Qxxp(9,9)));
ttc5=bil8(R10)/(s00*sqrt(Qxxp(10,10)));
ttd5=bil9(R10)/(s00*sqrt(Qxxp(11,11)));
tte5=bil10(R10)/(s00*sqrt(Qxxp(12,12)));
f5=n-u-me+1;

if abs(tta5)>tinv(0.95,f5) && abs(ttb5)>tinv(0.95,f5) && abs(ttc5)>tinv(0.95,f5) && abs(ttd5)>tinv(0.95,f5) && abs(tte5)>tinv(0.95,f5)
    uyusumsuz51=var1(R10);
    uyusumsuz52=var2(R10);
    uyusumsuz53=var3(R10);
    uyusumsuz54=var4(R10);
    uyusumsuz55=var5(R10);
    uyusumsuz5=[uyusumsuz51 uyusumsuz52 uyusumsuz53 uyusumsuz54 uyusumsuz55];
    uyusumsuz5=sort(uyusumsuz5);
    uyusumsuz5;
    durum5=1;
    artim5=1;
end




