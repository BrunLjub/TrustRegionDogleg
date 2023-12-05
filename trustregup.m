function [delta,retcode,xpprev,fpprev,xp,fp,maxtaken] = trustregup(n,xc,fc,f,gc,L,s,Newttaken,maxstep,steptol,delta,retcode,xpprev,fpprev)
    %provjerava je li novi korak dovoljno dobar

    %n dimenzija domene funkcije
    %xc trenutna toèka
    %fc f(xc)
    %f funkcija
    %gc gradijent od f u xc
    %L cholesky od Hessijana od f u xc
    %s trenutni korak
    %Newttaken, bool - je li korak Newtonov
    %maxstep - najveæi dopušteni korak
    %steptol - 
    %delta - trenutni delta
    %xprev - prethodni xc
    %fprev - prethodni f(xc)
    
    %retcode ... 0 dovoljno dobra nova iteracija
    %        ... 1 premala razlika izmeðu toèaka - nema napretka (ovisi o steptol)
    %        ... 2 f(xp) prevelik, smanji delta
    %        ... 3 f(xp) dovoljno mali, ali potencijalno možemo uzeti bolji korak   
    
    
    
    maxtaken = 0; %bool ...
    alpha = 1e-4; %parametar
    steplen = norm(s); %velièina trenutnog koraka
    xp = xc + s; %nova toèka
    fp = f(xp); %nova vrijednost
    deltaf = fp-fc;
    initslope = gc'*s;
    if retcode ~= 3
        fpprev = 0;
    end
    if retcode == 3 && ((fp >= fpprev) || deltaf > alpha * initslope) 
        
        %Ako išta od ovog vrijedi, ne prihvaæamo novi korak, vraæamo se na
        %prethodni i smanjujemo delta
        
        retcode = 0;
        xp = xpprev;
        fp = fpprev;
        delta = delta/2; 
        return;
    else
        if deltaf >= alpha*initslope
            
            %Inaèe, fp je prevelik
            
            rellength = max( abs(s)./max(abs(xp),ones(n,1)) );
            if rellength < steptol 
                %rellength je premal, prihvaæamo xp (iako nismo presretni)
                retcode = 1;
                xp = xc;
            else
                %nastavljamo, smanjujemo delta
                retcode = 2;
                deltatemp = -initslope * steplen /(2 * (deltaf-initslope));
                if deltatemp < 0.1*delta
                    delta = 0.1*delta;
                else
                    if deltatemp > 0.5*delta
                        delta = 0.5*delta;
                    else
                        delta = deltatemp;
                    end
                end
                return;
            end
        else 
            %fp je dovoljno mali
            %izraèunaj deltafpred
            deltafpred = initslope;
            for i = 1:n 
                temp = L(i:n,i)'*s(i:n);
                deltafpred = deltafpred + (temp^2/2);   
            end
            if retcode ~= 2 && ((abs(deltafpred - deltaf) <= 0.1*abs(deltaf)) || (deltaf <= initslope)) && Newttaken == 0 && (delta <= (0.99 * maxstep))
                %dupliraj delta (ako ne prelazi maxstep)
                retcode = 3;
                xpprev = xp;
                fpprev = fp;
                delta = min(2*delta,maxstep);
                return
            else
                %prihvaæamo xp
                retcode = 0;
                if steplen > 0.99 * maxstep
                    %Uzeli smo najveæi moguæi korak
                    maxtaken = 1;
                end
                
                if deltaf >= 0.1*deltafpred
                    %smanjujemo delta
                    delta = delta/2;
                else
                    if deltaf <= 0.75*deltafpred
                        %poveæavamo delta
                        delta = min(2*delta,maxstep);
                    else
                        %delta je nepromjenjen
                    end
                end            
            end            
        end
    end
    
end

