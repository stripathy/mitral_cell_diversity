

function [ac1,ac2,xc,ac1all,ac2all,xcall]=crossCorrSonya(ST1,ST2,tmax)

nsweeps=size(ST1,1);

span=50;
pad=0;

ac1=zeros(1,span*2+1);
ac2=zeros(1,span*2+1);
xc1=zeros(1,span*2+1);
xc2=zeros(1,span*2+1);
xc = zeros(1,span*2+1);;
ac1all=zeros(nsweeps,span*2+1);
ac2all=zeros(nsweeps,span*2+1);
xcall=zeros(nsweeps,span*2+1);


inst1=zeros(nsweeps,tmax);
inst2=zeros(nsweeps,tmax);

for ee=1:nsweeps
    k=round(ST1(ee,:));
    k=k(k>0);
    k2=round(ST2(ee,:));
    k2=k2(k2>0);
    for bb=0:pad*2
        inst1(ee,k+bb-pad)=1;
        inst2(ee,k2+bb-pad)=1;
    end
end

timeOffset = 0;
for aa=1:nsweeps
    
    if nnz(ST1(aa,:))>0
        
        div1=sum(ST1(aa,:)>(timeOffset+span) & ST1(aa,:)<(tmax-span));
        div2=sum(ST2(aa,:)>(timeOffset+span) & ST2(aa,:)<(tmax-span));
        for dd=1:nnz(ST1(aa,:))
            time=round(ST1(aa,dd));

            if time>(timeOffset+span) && time<(tmax-span)
                
                    ac1=ac1+(inst1(aa,time-span:time+span)-mean(inst1(aa,time-span:time+span)))/(div1*nsweeps);
                    ac1all(aa,:)=ac1all(aa,:)+(inst1(aa,time-span:time+span)-mean(inst1(aa,time-span:time+span)))/(div1);
                    if nnz(inst2(aa,:))>0
                    xc1=xc1+(inst2(aa,time-span:time+span)-(mean(inst2(aa,time-span:time+span))))/(div1*nsweeps);
                    xcall(aa,:)=xcall(aa,:)+(inst2(aa,time-span:time+span)-(mean(inst2(aa,time-span:time+span))))/(div2);
                    end
            end
        end
    end

    if nnz(ST2(aa,:))>0
        div1=sum(ST2(aa,:)>(timeOffset+span) & ST2(aa,:)<(tmax-span));
        div2=sum(ST1(aa,:)>(timeOffset+span) & ST1(aa,:)<(tmax-span));
        
        for dd=1:nnz(ST2(aa,:))
            time=round(ST2(aa,dd));

            if time>(timeOffset+span) && time<(tmax-span)
               
                ac2=ac2+(inst2(aa,time-span:time+span)-mean(inst2(aa,time-span:time+span)))/(div1*nsweeps);
                ac2all(aa,:)=ac2all(aa,:)+(inst2(aa,time-span:time+span)-mean(inst2(aa,time-span:time+span)))/(div1);
                if nnz(inst1(aa,:))>0
                xc2=xc2+(inst1(aa,time-span:time+span)-mean(inst1(aa,time-span:time+span)))/(div2*nsweeps);
                xcall(aa,:)=xcall(aa,:)+(inst1(aa,time-span:time+span)-mean(inst1(aa,time-span:time+span)))/(div1);
                end
            end
        end
    end
    
    
end

xc=(xc1+xc2)/2;
xcall=xcall/2;


end