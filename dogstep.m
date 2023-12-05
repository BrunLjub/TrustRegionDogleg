function [delta,firstdog,Cauchylen,eta,ssd,v,s,Newttaken] = dogstep(n,gc,L,sn,Newtlen,maxstep,delta,firstdog,Cauchylen,eta,ssd,v)
%   Radi double dogleg korak  

%   n dimenzija domene od f
%   gc, gradijent od f u xc
%   L, matrica iz Choleskyjeve dekompozicije Hessijana
%   sn, Newtonov korak
%   Newtlen, du�ina Newtonovog koraka
%   maxstep, najve�a dopu�tena du�ina koraka (zadaje korisnik)
%   delta, trenutni deltac
%   firstdog, bool - je li ovo prvi prolazak kroz petlju
%   Cauchylen, du�ina Cauchyjevog koraka
%   eta, du�ina Newtonovom koraku
%   ssd, Cauchyjev korak
%   v, konveksna kombinacija ssd i eta*sn

%   s, korak koji uzimamo
%   Newttaken, bool - jesmo li uzeli Newtonov korak
    

    if Newtlen <= delta
        %uzimamo Newtonov korak
        Newttaken = 1; 
        s = sn;
        %novi delta je veli�ina Newtonovog koraka
        delta = Newtlen;
        return;
    else
        %Newtonov korak je prevelik
        Newttaken = 0;
        if firstdog %Ako prvi put ulazimo u petlju
            firstdog = 0;
            
            
            alpha = norm(gc)^2;
            beta  = norm(L'*gc)^2;
            ssd = -(alpha/beta)*gc; %Cauchyjev korak
            Cauchylen = alpha * sqrt(alpha)/beta; %Du�ina Cauchyjevog koraka
            eta = 0.2 + (0.8 * (alpha^2)/(beta*abs(gc'*sn))); %gamma = (alpha^2)/(beta*abs(gc'*sn))
            
            v = eta*sn - ssd; %smjer od ssd do eta*sn
            if delta == -1 %Ako user ne postavi deltu (samo u prvom koraku)
                delta = min(Cauchylen,maxstep);
            end
        end
    end
    if eta * Newtlen <= delta %Uzmi parcijalni Newtonov korak
        s = (delta/Newtlen)*sn;
        return;
    else
        if Cauchylen >= delta %Uzmi parcijalni Cauchyjev korak
            s = (delta/Cauchylen)*ssd;
            return;
        else %Izra�unaj konveksnu kombinaciju ssd i eta*sn koji le�i na trust-region kru�nici
            temp = v'*ssd;
            tempv = v'*v;
            %Rje�enje kv. jednad�be ||ssd+lambda*v||^2=delta^2
            lambda = (-temp + sqrt(temp^2 - tempv*(Cauchylen^2-delta^2)))/tempv;
            s = ssd + lambda*v;
        end
    end
end

