function pCO2 = pCO2_fun(DIC,alk,pCO2_dict)
%
% Calculate pCO2 through combination of a lookup table (which is passed
% to the routine as a dictionary) and bilinear interpolation. This is done
% because individual calls of f_csys_alk_DIC.m are relativley slow, while
% calls of this function are fast. The function is called for every
% gridcell at every timestep, and thus has to be fast.
%
% inputs: DIC (scalar of vector) [uM]
%         alk (scalar or vector; same size as DIC) [uM]
%         pCO2_dict dictionary of pCO2 values taylored to this case
%
% output: pCO2 (scalar or vector; same size as inputs) [uatm]
%

nx = length(DIC);
pCO2 = zeros(size(DIC));

for i=1:nx

    x = DIC(i);
    x1 = floor(x);
    x2 = ceil(x);
    y = alk(i);
    y1 = floor(y);
    y2 = ceil(y);

    if x1 ~= x2
        % need to interpolate in x direction

        if y1 ~= y2
            % need to interpolate in y direction too, i.e. bilinear
            
            key11 = x1*10000 + y1;
            key12 = x1*10000 + y2;
            key21 = x2*10000 + y1;
            key22 = x2*10000 + y2;
            
            pCO2_11 = lookup(pCO2_dict,key11);
            pCO2_12 = lookup(pCO2_dict,key12);
            pCO2_21 = lookup(pCO2_dict,key21);
            pCO2_22 = lookup(pCO2_dict,key22);
            
            % then perform bilinear interpolation from the neighboring pCO2 values
            pCO2(i) = 1./(x2-x1)./(y2-y1).*(pCO2_11.*(x2-x).*(y2-y)+pCO2_21.*(x-x1).*(y2-y) ...
                                        +pCO2_12.*(x2-x).*(y-y1)+pCO2_22.*(x-x1).*(y-y1));

        else
            % no interpolation in y direction, just x

            key1 = x1*10000 + y1;
            key2 = x2*10000 + y1;
            
            pCO2_1 = lookup(pCO2_dict,key1);
            pCO2_2 = lookup(pCO2_dict,key2);

            % linear interpolation in x direction
            pCO2(i) = (x2-x)/(x2-x1)*pCO2_1 + (x-x1)/(x2-x1)*pCO2_2;
        end
    else
        if y1 ~= y2
            % need to interpolate in y direction only
            key1 = x1*10000 + y1;
            key2 = x1*10000 + y2;
            
            pCO2_1 = lookup(pCO2_dict,key1);
            pCO2_2 = lookup(pCO2_dict,key2);

            % linear interpolation in y direction
            pCO2(i) = (y2-y)/(y2-y1)*pCO2_1 + (y-y1)/(y2-y1)*pCO2_2;
        else
            % no interpolation at all

            key = x*10000 + y;
            pCO2(i) = lookup(pCO2_dict,key);
        end
    end
end

end