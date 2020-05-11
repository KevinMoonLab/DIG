function PDX = P_distances(P_t, pot_method, Gamma)

switch pot_method
    case 'info_geometry'
        disp 'using information geometry distance (arccos)'
                Pot = P_t;
                PDX = P_t;
        [si, ~] = size(Pot);
        for o = 1:si
            for oo = o:si
                if o == oo 
                    PDX(o,oo)=0;
                else 
                pro = Pot(o,:).*Pot(oo,:);
                PDX(o,oo) = acos(sum(pro.^(1/2)));
                PDX(oo,o) = PDX(o,oo);
                end 
            end 
        end 
        Pot = PDX;
    case 'log'
        disp 'using -log(P) potential distance'
        pot_eps = 10e-7;
        Pot = -log(P_t + pot_eps);
        PDX = squareform(pdist(Pot, 'euclidean'));
    case 'sqrt'
        disp 'using sqrt(P) potential distance'
        Pot = sqrt(P_t);
        PDX = squareform(pdist(Pot, 'euclidean'));
    case 'none'
        Pot = P_t;
        PDX = squareform(pdist(Pot, 'euclidean'));
    case 'gamma'
        pot_eps = 10e-7;
        Pot = -log(P_t + pot_eps);
        PDX = squareform(pdist(Pot, 'euclidean'));
        if Gamma == -1
            Pot = P_t;
            PDX = squareform(pdist(Pot, 'euclidean'));
        elseif Gamma == 1
            pot_eps = 10e-7;
            Pot = -log(P_t + pot_eps);
            PDX = squareform(pdist(Pot, 'euclidean'));
        elseif abs(Gamma) < 1
            disp 'Pot = 2/(1-\gamma)*P^((1-\gamma)/2)'
            disp(['gamma = ' num2str(Gamma)]);
            Gamma = min(Gamma, 0.95);
            Pot = 2/(1-Gamma)*P_t.^((1-Gamma)/2); 
            PDX = squareform(pdist(Pot, 'euclidean'));
        else
            error 'not a valid gamma value'
        end 
    otherwise
        error 'potential method unknown'
end


end 