function [retcode,xp,fp,delta,maxtaken] = dogdriver(n,xc,fc,f,gc,L,sn,maxstep,steptol,delta)
%   Pronaði xp na dvostrukoj dogleg krivulji tako da zadovoljava
%   f(xp) <= f(xc) + alpha*gc(xp-xc) i vrati skaliran delta

%   n ... dimenzija domene funkcije f
%   xc ... poèetna toèka
%   fc ... vrijednost poèetne toèke
%   f ... funkcija
%   gc ... aproksimacija gradijenta od f u xc
%   delta ... radijus trust-regiona
%   L ... H = LL' (L iz Cholesky faktorizacije Hessijana)
%   sn ... Newtonov korak
%   maxstep ... vrijednost koji korak ne smije preæi
%   steptol ... ako je korak premali, prekidamo


%   retcode ... 0 sve u redu, 
%           ... 1 rutina ne nalazi xp dovoljno razlièit od xc
%           ... 4 inicijalni poziv funkcije dogstep
%   xp ... pronaðena toèka
%   fp ... vrijednost u pronaðenoj toèki
%   delta ... novi delta
%   maxtaken ... jesmo li uzeli newtonov korak ili ne

    retcode = 4;
    firstdog = 1;
    Newtlen =  norm(sn);
    
    %U provm koraku se neæe koristiti sljedeæe vrijednosti
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

