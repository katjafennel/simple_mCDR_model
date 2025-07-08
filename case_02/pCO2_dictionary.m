function dict = pCO2_dictionary(T,S,DIC_min,DIC_max,alk_min,alk_max,del)

% ranges of DIC and alk
DIC_v = DIC_min:del:DIC_max;
alk_v = alk_min:del:alk_max;

% turn into matrices
DIC = DIC_v'*ones(1,length(alk_v));
alk = ones(length(DIC_v),1)*alk_v;

% initialize pCO2 matrix
pCO2 = zeros(size(DIC));

% calculate values
for i=1:length(pCO2(:))
    pCO2(i) = f_csys_alk_DIC(T,S,alk(i),DIC(i)); % [uatm]
end

% create dictionary
keys = DIC(:)*10000+alk(:);
values = pCO2(:);
dict = dictionary(keys,values);
end