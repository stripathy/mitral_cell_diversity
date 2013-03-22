

function [ac1,ac2,xc,ac1all,ac2all,xcall]=crossCorrSonyaTest(ST1,ST2,tmax)

nsweeps=size(ST1,1);

span=50;
pad=0;
corrFxnLen = span*2+1;

ac1=zeros(1,span*2+1);
ac2=zeros(1,span*2+1);
xc=zeros(1,span*2+1);
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
        inds1 = makeValidInds(k+bb-pad, tmax);
        inds2 = makeValidInds(k2+bb-pad, tmax);
        inst1(ee,inds1)=1;
        inst2(ee,inds2)=1;
    end
end

timeOffset = 0;
for aa=1:nsweeps
    
    if nnz(ST1(aa,:))>0
        
        div1=sum(ST1(aa,:)>(timeOffset) & ST1(aa,:)<(tmax));
        div2=sum(ST2(aa,:)>(timeOffset) & ST2(aa,:)<(tmax));
        for dd=1:nnz(ST1(aa,:))
            time=round(ST1(aa,dd));

            if time>=(timeOffset) && time<=(tmax)
                    windowInds = time-span:time+span;
                    inds1 = makeValidInds(windowInds,tmax);
                    corrInds = windowInds(windowInds>0 & windowInds <=tmax) - time + span+1;
                    
                    ac1(corrInds)=ac1(corrInds)+(inst1(aa,inds1)-mean(inst1(aa,inds1)))/(div1*nsweeps);
                    ac1all(aa,corrInds)=ac1all(aa,corrInds)+(inst1(aa,inds1)-mean(inst1(aa,inds1)))/(div1);
                    if nnz(inst2(aa,inds1))>0
                    xc(corrInds)=xc(corrInds)+(inst2(aa,inds1)-(mean(inst2(aa,inds1))))/(div2*nsweeps);
                    xcall(aa,corrInds)=xcall(aa,corrInds)+(inst2(aa,inds1)-(mean(inst2(aa,inds1))))/(div2);
                    end
            end
        end
    end

    if nnz(ST2(aa,:))>0
        div1=sum(ST2(aa,:)>(timeOffset) & ST2(aa,:)<(tmax));
        div2=sum(ST1(aa,:)>(timeOffset) & ST1(aa,:)<(tmax));
        
        for dd=1:nnz(ST2(aa,:))
            time=round(ST2(aa,dd));

            if time>=(timeOffset) && time<=(tmax)
                    windowInds = time-span:time+span;
                    inds1 = makeValidInds(windowInds,tmax);
                    corrInds = windowInds(windowInds>0 & windowInds <=tmax) - time + span+1;
               
                ac2(corrInds)=ac2(corrInds)+(inst2(aa,inds1)-mean(inst2(aa,inds1)))/(div1*nsweeps);
                ac2all(aa,corrInds)=ac2all(aa,corrInds)+(inst2(aa,inds1)-mean(inst2(aa,inds1)))/(div1);
                if nnz(inst1(aa,inds1))>0
                xc(corrInds)=xc(corrInds)+(inst1(aa,inds1)-mean(inst1(aa,inds1)))/(div2*nsweeps);
                xcall(aa,corrInds)=xcall(aa,corrInds)+(inst1(aa,inds1)-mean(inst1(aa,inds1)))/(div2);
                end
            end
        end
    end
    
    
end

xc=xc/2;
xcall=xcall/2;


end