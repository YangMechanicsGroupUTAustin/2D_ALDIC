% ==================================
% To compute strain on a uniform mesh
% ----------------------------------
% switch MethodToComputeStrain
%   case 0: direct solved results;
%   case 1: central finite difference;
%   case 2: plane fitting method;
%   case 3: finite element Gauss points;
% ==================================

switch DICpara.MethodToComputeStrain

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0 % ALDIC directly solved deformation gradients
        FStrain = FLocal;   
        FStraintemp = FStrain;
        Rad = 0;
        
        temp = 1:1:size(coordinatesFEM,1); temp = temp';
        temp = reshape(temp,M,N); temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad);
        temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad),1);
        
        temp3 = zeros(4*(M-2*Rad)*(N-2*Rad),1);
        for i = 1:(M-2*Rad)*(N-2*Rad)
            temp3(4*i-3:4*i) = [4*temp2(i)-3; 4*temp2(i)-2; 4*temp2(i)-1; 4*temp2(i)];
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1 % Central finite difference
         
        % Compute strain method I: Use Finite difference operator 
        D = funDerivativeOp(M,N,DICpara.winstepsize); % D = sparse(4*M*N, 2*M*N);
        FStrain = D*reshape(ULocal,length(ULocal),1);
        
        Rad = 1;  
        
        % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
        temp = 1:1:size(coordinatesFEM,1); temp = temp';
        temp = reshape(temp,M,N); temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad);
        temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad),1);
        
        temp3 = zeros(4*(M-2*Rad)*(N-2*Rad),1);
        for i = 1:(M-2*Rad)*(N-2*Rad)
            temp3(4*i-3:4*i) = [4*temp2(i)-3; 4*temp2(i)-2; 4*temp2(i)-1; 4*temp2(i)];
        end
        FStraintemp = FStrain(temp3);
         
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2 % Plane fitting
        D = funDerivativeOp(M,N,DICpara.winstepsize); % D = sparse(4*M*N, 2*M*N);
        FStrain = D*reshape(ULocal,length(ULocal),1);  
        
        % Compute strain method II: Use Plane Fitting method
        prompt = 'What is your half window size: ';
        Rad = input(prompt);
        
        [Uytemp, Uxtemp, UNewtemp] = PlaneFit2(reshape(ULocal(1:2:end),M,N), DICpara.winstepsize, DICpara.winstepsize,Rad);
        [Vytemp, Vxtemp, VNewtemp] = PlaneFit2(reshape(ULocal(2:2:end),M,N), DICpara.winstepsize, DICpara.winstepsize,Rad);
        
        FStraintemp = zeros(4*(M-2*Rad)*(N-2*Rad),1);
        FStraintemp(1:4:end) = reshape(Uxtemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        FStraintemp(2:4:end) = reshape(Vxtemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        FStraintemp(3:4:end) = reshape(Uytemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        FStraintemp(4:4:end) = reshape(Vytemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        
        % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
        temp = 1:1:size(coordinatesFEM,1); temp = temp';
        temp = reshape(temp,M,N); temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad);
        temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad),1);
         
        temp3 = zeros(4*(M-2*Rad)*(N-2*Rad),1);
        for i = 1:(M-2*Rad)*(N-2*Rad)
            temp3(4*i-3:4*i) = [4*temp2(i)-3; 4*temp2(i)-2; 4*temp2(i)-1; 4*temp2(i)];
        end
        
        FStrain(temp3) = FStraintemp;
        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3 % Finite element method
        
        GaussPtOrder = 2; [FStrain] = funGlobal_NodalStrainAvg(coordinatesFEM,elementsFEM,USubpb2,GaussPtOrder);
        Rad = 1;  
        
        % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
        temp = 1:1:size(coordinatesFEM,1); temp = temp';
        temp = reshape(temp,M,N); temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad);
        temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad),1);
        
        temp3 = zeros(4*(M-2*Rad)*(N-2*Rad),1);
        for i = 1:(M-2*Rad)*(N-2*Rad)
            temp3(4*i-3:4*i) = [4*temp2(i)-3; 4*temp2(i)-2; 4*temp2(i)-1; 4*temp2(i)];
        end
        FStraintemp = FStrain(temp3);
        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
        disp('Wrong Input to compute strain field!')
        
end


%% Update infinitesimal strain to other finite strains
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

FStraintemp = FStrainFinite(temp3);
FStrainWorld = FStraintemp; FStrainWorld(2:4:end) = -FStrainWorld(2:4:end); FStrainWorld(3:4:end) = -FStrainWorld(3:4:end); 
    