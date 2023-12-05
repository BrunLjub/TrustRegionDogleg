function [x,k,X,FX,mi,Deltas] = dogleg(f,df,df2,x0,dim,epsilon,delta,N, deltac)
    x = x0;    
    gc = df(x0);
    fx1 = f(x0);
    X = [x0'];
    FX = [fx1];
    mi= [];
    k = 0;
    x1 = x0;
    Deltas = [];
    
    while norm(gc) >= delta && k < N
        k = k+1;
        H = df2(x1);
        %izraèunaj mik
        d = eig(H);
        d = sort(d,'descend');
        if d(1)==d(dim) && d(dim)<1e-5
            mik = 1e-5-d(dim);
        else
            if d(1) > 1e5*d(dim) || (d(dim) < 1e5*d(1) && d(1)<=0 )
                if abs(d(dim)) < d(1)
                    mik = 1e-5*d(1)-d(dim);
                else
                    mik = -d(dim)*(1+1e-5);
                end
            else
                mik = 0;
            end
        end
        if mik ~= 0 %ako mik != 0, H = H+mik*I
            for j = 1:dim
                H(j,j) = H(j,j) + mik;
            end
        end
        mi = [mi,mik];
        
        L = chol(H,'lower'); %LL's = -df(x1)        
        %Newtonov korak
        sn = -L'\(L\gc);
        
        [~,x1,fx1,deltac,~] = dogdriver(dim,x1,fx1,f,gc,L,sn,1,1e-10,deltac);
        Deltas = [Deltas;deltac];
        X = [X;x1'];
        FX = [FX,fx1];
        if max(abs((x1 - x0)./x0)) < epsilon
            break;
        end
        x0 = x1;
        gc = df(x1);
        fx1 = f(x1);
    end
    x = x1;
    
end
    
    


