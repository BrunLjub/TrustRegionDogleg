function [retcode,xp,fp,delta,maxtaken] = dogdriver(n,xc,fc,f,gc,L,sn,maxstep,steptol,delta)
%   Prona�i xp na dvostrukoj dogleg krivulji tako da zadovoljava
%   f(xp) <= f(xc) + alpha*gc(xp-xc) i vrati skaliran delta

%   n ... dimenzija domene funkcije f
%   xc ... po�etna to�ka
%   fc ... vrijednost po�etne to�ke
%   f ... funkcija
%   gc ... aproksimacija gradijenta od f u xc
%   delta ... radijus trust-regiona
%   L ... H = LL' (L iz Cholesky faktorizacije Hessijana)
%   sn ... Newtonov korak
%   maxstep ... vrijednost koji korak ne smije pre�i
%   steptol ... ako je korak premali, prekidamo


%   retcode ... 0 sve u redu, 
%           ... 1 rutina ne nalazi xp dovoljno razli�it od xc
%           ... 4 inicijalni poziv funkcije dogstep
%   xp ... prona�ena to�ka
%   fp ... vrijednost u prona�enoj to�ki
%   delta ... novi delta
%   maxtaken ... jesmo li uzeli newtonov korak ili ne

    retcode = 4;
    firstdog = 1;
    Newtlen =  norm(sn);
    
    %U provm koraku se ne�e koristiti sljede�e vrijednosti
    Cauchylen = NaN;
    eta = NaN;
    ssd = NaN;
    v = NaN;
    xpprev = NaN;
    fpprev = NaN;
    
    while 1
        [delta,firstdog,Cauchylen,eta,ssd,v,s,Newttaken] = dogstep(n,gc,L,sn,Newtlen,maxstep,delta,firstdog,Cauchylen,eta,ssd,v);
        [delta,retcode,xpprev,fpprev,xp,fp,maxtaken] = trustregup(n,xc,fc,f,gc,L,s,Newttaken,maxstep,steptol,delta,retcode,xpprev,fpprev);

        if retcode < 2
            break
        end
    end
end

