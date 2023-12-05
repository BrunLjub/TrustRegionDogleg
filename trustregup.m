function [delta,retcode,xpprev,fpprev,xp,fp,maxtaken] = trustregup(n,xc,fc,f,gc,L,s,Newttaken,maxstep,steptol,delta,retcode,xpprev,fpprev)
    %provjerava je li novi korak dovoljno dobar

    %n dimenzija domene funkcije
    %xc trenutna to�ka
    %fc f(xc)
    %f funkcija
    %gc gradijent od f u xc
    %L cholesky od Hessijana od f u xc
    %s trenutni korak
    %Newttaken, bool - je li korak Newtonov
    %maxstep - najve�i dopu�teni korak
    %steptol - 
    %delta - trenutni delta
    %xprev - prethodni xc
    %fprev - prethodni f(xc)
    
    %retcode ... 0 dovoljno dobra nova iteracija
    %        ... 1 premala razlika izme�u to�aka - nema napretka (ovisi o steptol)
    %        ... 2 f(xp) prevelik, smanji delta
    %        ... 3 f(xp) dovoljno mali, ali potencijalno mo�emo uzeti bolji korak   
    
    
    
    maxtaken = 0; %bool ...
    alpha = 1e-4; %parametar
    steplen = norm(s); %veli�ina trenutnog koraka
    xp = xc + s; %nova to�ka
    fp = f(xp); %nova vrijednost
    deltaf = fp-fc;
    initslope = gc'*s;
    if retcode ~= 3
        fpprev = 0;
    end
    if retcode == 3 && ((fp >= fpprev) || deltaf > alpha * initslope) 
        
        %Ako i�ta od ovog vrijedi, ne prihva�amo novi korak, vra�amo se na
        %prethodni i smanjujemo delta
        
        retcode = 0;
        xp = xpprev;
        fp = fpprev;
        delta = delta/2; 
        return;
    else
        if deltaf >= alpha*initslope
            
            %Ina�e, fp je prevelik
            
            rellength = max( abs(s)./max(abs(xp),ones(n,1)) );
            if rellength < steptol 
                %rellength je premal, prihva�amo xp (iako nismo presretni)
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
            %izra�unaj deltafpred
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
                %prihva�amo xp
                retcode = 0;
                if steplen > 0.99 * maxstep
                    %Uzeli smo najve�i mogu�i korak
                    maxtaken = 1;
                end
                
                if deltaf >= 0.1*deltafpred
                    %smanjujemo delta
                    delta = delta/2;
                else
                    if deltaf <= 0.75*deltafpred
                        %pove�avamo delta
                        delta = min(2*delta,maxstep);
                    else
                        %delta je nepromjenjen
                    end
                end            
            end            
        end
    end
    
end

