function [index,findex]=surfpatch(Source,vertices,faces,d0)
%% Initialization
index  = Source; %Indices of patch's elements
findex = Source; %Indices of patch's frontier elements
d      = 0;
%% Cycle while the geodesic distance of the frontier to the center is not greater than 'd0'
while d<=d0
    findex_new  = [];
    for cont = 1:length(findex)
        fpoint        = findex(cont);             %Pick point at the frontier 'fpoint'
        %% Search the neighbors of 'fpoint' out of the patche 'nfpoint'
        [row,col]      = find(faces==fpoint);      
        neig_fpoint    = faces(row,:);
        neig_fpoint    = neig_fpoint(:);
        neig_fpoint    = setdiff(neig_fpoint,index);
        findex_new     = [findex_new; setdiff(neig_fpoint,findex_new)];
    end
    findex = findex_new;
    index  = [index; findex];
    %% Compute geodesic distance of the frontier points 'd'
    d      = vertices(findex,:) - repmat(vertices(Source,:),length(findex),1);
    d      = mean(sqrt(sum(d.^2,2)));
%     J      = zeros(length(vertices),1); 
%     J(findex) = 5;
%     Plot_sources_Haufe(J,vertices,faces,'simple')
end
end