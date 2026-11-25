% ==================================
% To compute strain on a quadtree mesh
% ----------------------------------
% switch MethodToComputeStrain
%   case 3: finite element Gauss points;
% ----------------------------------------------
% Author: Jin Yang.
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 2020.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch DICpara.MethodToComputeStrain

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0 % ALDIC directly solved deformation gradients
        FSubpb2=FLocal;
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1 % Central finite difference
        
        % Compute strain method I: Use Finite difference operator or FEM solver
        minCoordStep = min( [DICmesh.elementMinSize] );
        xList = [ceil(min(coordinatesFEM(:,1))) : minCoordStep : floor(max(coordinatesFEM(:,1)))]';
        yList = [ceil(min(coordinatesFEM(:,2))) : minCoordStep : floor(max(coordinatesFEM(:,2)))]';
        
        [xGrid,yGrid] = ndgrid(xList, yList);
        uGrid = gridfit( coordinatesFEM(:,1), coordinatesFEM(:,2), ULocal(1:2:end), xList, yList,'regularizer','springs' ); uGrid=uGrid';
        vGrid = gridfit( coordinatesFEM(:,1), coordinatesFEM(:,2), ULocal(2:2:end), xList, yList ,'regularizer','springs'); vGrid=vGrid';
         
        for tempi = 1:length(uGrid(:))
           tempx = xGrid(tempi); tempy = yGrid(tempi);
           if  Df.ImgRefMask(tempx,tempy) == 0
               uGrid(tempi) = nan; 
               vGrid(tempi) = nan;
           end
        end
        
        %%%%%%%%%%%%%%%%%%%%%
        % ====== Central finite difference ======
        % uvGrid = [uGrid(:),vGrid(:)]'; uvGrid=uvGrid(:);
        % 
        % D = funDerivativeOp(length(xList),length(yList),minCoordStep);
        % FStrainGrid = D*uvGrid(:);
        % 
        % figure, surf(reshape(FStrainGrid(1:4:end),length(xList),length(yList)),'edgecolor','none');
        % view(2);  axis equal; axis tight; colorbar; caxis([-0.2,0.2])
        
        % Try to use regularization  
        % tempF11 = regularizeNd([xGrid(:),yGrid(:)],FStrainGrid(1:4:end),{xList,yList},DICpara.smoothness);
        % tempF21 = regularizeNd([xGrid(:),yGrid(:)],FStrainGrid(2:4:end),{xList,yList},DICpara.smoothness);
        % tempF12 = regularizeNd([xGrid(:),yGrid(:)],FStrainGrid(3:4:end),{xList,yList},DICpara.smoothness);
        % tempF22 = regularizeNd([xGrid(:),yGrid(:)],FStrainGrid(4:4:end),{xList,yList},DICpara.smoothness);
        % figure, surf(tempF22,'edgecolor','none');  view(2);  axis equal; axis tight; colorbar; caxis([-0.2,0.2])
         
          
        Rad = 1*ceil(mean(DICpara.winstepsize)/minCoordStep/2);
        [UNew,dudx,dudy,iGrid,jGrid] = PlaneFit22(reshape(uGrid,length(xList),length(yList)), minCoordStep, minCoordStep, Rad);
        [VNew,dvdx,dvdy,iGrid,jGrid] = PlaneFit22(reshape(vGrid,length(xList),length(yList)), minCoordStep, minCoordStep, Rad);
    
        [row,col] = find(isnan(VNew)==1);
        nonNanIndtemp = sub2ind([length(xList)-2*Rad,length(yList)-2*Rad], row,col);
        nonNanInd = sub2ind([length(xList),length(yList)], iGrid(nonNanIndtemp),jGrid(nonNanIndtemp));
        F_F11 = scatteredInterpolant(xGrid(nonNanInd),yGrid(nonNanInd),dudx(nonNanIndtemp));
        F11 = F_F11(coordinatesFEM(:,1),coordinatesFEM(:,2));
        F_F21 = scatteredInterpolant(xGrid(nonNanInd),yGrid(nonNanInd),dvdx(nonNanIndtemp));
        F21 = F_F21(coordinatesFEM(:,1),coordinatesFEM(:,2));
        F_F12 = scatteredInterpolant(xGrid(nonNanInd),yGrid(nonNanInd),dudy(nonNanIndtemp));
        F12 = F_F12(coordinatesFEM(:,1),coordinatesFEM(:,2));
        F_F22 = scatteredInterpolant(xGrid(nonNanInd),yGrid(nonNanInd),dvdy(nonNanIndtemp));
        F22 = F_F22(coordinatesFEM(:,1),coordinatesFEM(:,2));
         
        FSubpb2 = [F11(:),F21(:),F12(:),F22(:)]'; FSubpb2 = FSubpb2(:);
        %for tempk=0:3
        %    FSubpb2(4*DICmesh.markCoordHoleEdge-tempk) = FSubpb1(4*DICmesh.markCoordHoleEdge-tempk);
        %end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2 % Plane fitting
        % Not implemented yet.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3 % Finite element method
        
        FSubpb2 = FLocal;
        try
            if DICpara.DoYouWantToSmoothOnceMore == 0
                FSubpb2 = funSmoothStrainQuadtree(FSubpb2,DICmesh,DICpara);
                for tempk=0:3
                    FSubpb2(4*DICmesh.markCoordHoleEdge-tempk) = FLocal(4*DICmesh.markCoordHoleEdge-tempk);
                end
            end
        catch
        end
          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
        disp('Wrong Input to compute strain field!')
        
end




%% Update infinitesimal strain to other finite strains
FStrain = FSubpb2;
FStrainFinite = FStrain;
for tempi = 1:4:length(FStrain)
    
    % Obtain each component of def grad tensor
    dudx = FStrain(tempi);
    dvdx = FStrain(tempi+1);
    dudy = FStrain(tempi+2);
    dvdy = FStrain(tempi+3);
    
    switch DICpara.StrainType
        case 0 % Infinitesimal stran
            % Do nothing
        case 1 % Eluerian strain
            FStrainFinite(tempi) = 1/(1-dudx)-1;
            FStrainFinite(tempi+3) = 1/(1-dvdy)-1;
            FStrainFinite(tempi+2) = dudy/(1-dvdy);
            FStrainFinite(tempi+1) = dvdx/(1-dudx);
        case 2 % Green-Lagrangian strain: E=(C-I)/2
            FStrainFinite(tempi) = 0.5*(dudx*2+dudx^2+dvdx^2);
            FStrainFinite(tempi+3) = 0.5*(dvdy*2+dudy^2+dvdy^2);
            FStrainFinite(tempi+2) = 0.5*(dudy+dvdx+dudx*dudy+dvdx*dvdy);
            FStrainFinite(tempi+1) = 0.5*(dvdx+dudy+dudy*dudx+dvdy*dvdx);
        case 3
            disp('Press "Ctrl+C" to modify by yourself.'); pause;
        otherwise
            disp('Wrong strain type!');
    end
    
end

FStraintemp = FStrainFinite;
FStrainWorld = FStraintemp; FStrainWorld(2:4:end) = -FStrainWorld(2:4:end); FStrainWorld(3:4:end) = -FStrainWorld(3:4:end); 


