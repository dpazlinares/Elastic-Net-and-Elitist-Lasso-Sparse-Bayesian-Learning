function Plot_sources_Haufe(solution,vertices,faces,plottype)
% Plot inverse solutions in a surface
% (or any data in the corresponding surface) using:
% Haufe's tools if plottype = 'haufe' (default)
% Simple matlab view if plottype = 'simple'
%
% INPUT
% solution -> data to plot. Each column is plotted in separate figures
%             No. of rows must match the no. of vertices
% vertices -> Nx3 matrix where each row is the coordinates of vertices of the surface mesh
% faces    -> Mx3 matrix where each row contains the indices of vertices conforming each face or triangle of the surface mesh
% plottype -> define type of plot. for using Haufe's plot the toolbox is needed
%
% Eduardo adapated from Mayrim, who adpated from Haufe... ;) 2017
% Ejemplo:
%   Plot_sources_Haufe(zeros(81924,1),Patch.vertices,Patch.faces,'simple');


if size(faces,2)==4; faces(:,1)=[];elseif size(faces,2)~=3; error('incorrect input "faces"');end
if nargin<4 || isempty(plottype) || isnumeric(plottype);plottype='haufe';end
switch lower(plottype)
    case 'simple'
        for i= 1:size(solution,2)
            colormap jet;
            patch('Vertices',vertices,'Faces',faces,'FaceVertexCData',solution(:,i),'EdgeColor','b','FaceAlpha',1,'FaceColor','interp'); %[202 189 38]/255);%
            title(['Column ' num2str(i)]);
            caxis([-5 5]);
        end
    otherwise   %'haufe'
        vc.tri = faces;
        vc.vc_ori = vertices;
        para.colormap = jet;
        for i= 1:size(solution,2)
            colormap jet;
            showsurface(vc,para,solution(:,i));
            title(['Column ' num2str(i)]);            
            caxis([-5 5]);
        end
end