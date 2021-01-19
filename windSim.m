function [u,v,w,t,nodes] = windSim(filename,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [u,v,w,t,nodes] = windSim(filename) generates spatially correlated wind
% histories based on an input file "filename" which is a text file.
% The space is defined using a cartesian coordinate system (x,y,z).
% The x axis is the horizontal axis aligned with the wind direction.
% The y axis is the horizontal axis normal to the wind direction.
% The z axis is the veertical axis.
% In this simulation the flow is assumed to have only a mean value in the x
% direction. In other words, mean(v)=mean(w)=0 m/s. This is a common
% assumption in wind engineering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: text file .
% example: filename = 'INPUT.txt'
% example: [u,v,w,t,nodes] = windSim(filename,'geometry','circle.txt')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u: [Nyy x N] matrix of the along wind component
% v: [Nyy x N] matrix of the cross-wind component
% w: [Nyy x N] matrix of the vertical wind component
% t: time vector
% nodes: structure variables that contains informations about mean wind
% speed, and coordinates of each nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example: [u,w,t,nodes] = windSim('INPUT.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author. Etienne Cheynet  -- last modified: 01-05-2020
%  see also windSim2.m coherence.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('geometry',[]);
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
geometry = p.Results.geometry ;




%% TIME DEFINITION
% import data from input file
[data] = importfile(filename, 7, 8,['%*s%f%*[^\n]'],'\t');
data=cell2mat(data);
if data(1)<0, error('"fs" must be positiv'); end
if data(2)<0, error('The variable "Duration" must be positiv'); end
fs=data(1); % sampling frequency
tmax=data(2); % duration of time series
dt=1./fs; % time step
[t,f] = getSamplingPara(nextpow2(tmax*fs),fs);
Nfreq = numel(f);

%% WIND DATA
[data] = importfile(filename, 12, 23,['%*s%f%*[^\n]'],'\t');
data=cell2mat(data);
if any(data(1:3)<0), error('The standard deviations of wind velocity must be larger than 0'); end
if any(data(4:6)<=0), error('The turbulence length scales must be positiv'); end
if any(data(7:12)<0), error('The decay coefficients must be positiv'); end

stdU=data(1);
stdV=data(2);
stdW=data(3);

Lux=data(4);
Lvx=data(5);
Lwx=data(6);

Cuy=data(7);
Cuz=data(8);
Cvy=data(9);

Cvz=data(10);
Cwy=data(11);
Cwz=data(12);

[typeWind] = importfile(filename, 29, 30,['%*s%s%*[^\n]'],'\t');
type=char(typeWind{1});
windProfile=char(typeWind{2}); % power or log

[data] = importfile(filename, 32,36,['%*s%f%*[^\n]'],'\t');
data=cell2mat(data);
Uref=data(1);
zr = data(2);
a=data(3);
u_star=data(4);
z0=data(5);

%% GRID GENERATION and node structure

[nodes,Nm] = createGrid(geometry,filename);
M = 3*Nm; % u,v and w are generated together, requiring 3 times more points than Nm (u,v and w are not necessarily independant)
%% Input data for wind coherence
dy = abs(bsxfun(@minus,nodes.Y(:)',nodes.Y(:))); % Matrix distance along y
dz = abs(bsxfun(@minus,nodes.Z(:)',nodes.Z(:))); % Matrix distance along z
meanU = 0.5*abs(bsxfun(@plus,nodes.U(:)',nodes.U(:))); % Mean wind velocity between each nodes

%% Compute the core spectral matrix A
A = zeros(Nfreq,M);
for ii=1:Nfreq
    if ii==2,    tStart = tic; end
    randPhase = rand(M,1); % set random phase
    
    % compute the coherence at each frequency step
    [cohU] = cohDavenport(meanU,dy,dz,f(ii),Cuy,Cuz);
    [cohV] = cohDavenport(meanU,dy,dz,f(ii),Cvy,Cvz);
    [cohW] = cohDavenport(meanU,dy,dz,f(ii),Cwy,Cwz);

    
    % turbulence spectrum
    if strcmpi(type,'von karman')
        Su= VonKarmanSpectrum(f(ii),nodes.U,stdU,Lux,'u');
        Sv= VonKarmanSpectrum(f(ii),nodes.U,stdV,Lvx,'v');
        Sw= VonKarmanSpectrum(f(ii),nodes.U,stdW,Lwx,'w');
        Suw = zeros(size(Su));
        Svw = zeros(size(Su));
    elseif strcmpi(type,'kaimal')
        Su= KaimalSpectrum(f(ii),nodes.U,u_star,nodes.Z,'u');
        Sv= KaimalSpectrum(f(ii),nodes.U,u_star,nodes.Z,'v');
        Sw= KaimalSpectrum(f(ii),nodes.U,u_star,nodes.Z,'w');
        Suw= KaimalSpectrum(f(ii),nodes.U,u_star,nodes.Z,'uw');
        Svw = zeros(size(Su)); % Can be modified by the user
    else
        error(' spectrum type is unknown')
    end
    
    
    Suu = sqrt(Su*Su').*cohU;
    Svv = sqrt(Sv*Sv').*cohV;
    Sww = sqrt(Sw*Sw').*cohW;
    Suw2 = sqrt(Suw*Suw').*sqrt(cohU.*cohW); % The cross-coherence is here not very well defined, but this should be good enough at a first approximation
    Svw2 = sqrt(Svw*Svw').*sqrt(cohU.*cohW); % The cross-coherence is here not very well defined, but this should be good enough at a first approximation
    
    ZERO = zeros(size(Suu));
    S = [Suu,   ZERO,   Suw2;...
        ZERO,  Svv,    Svw2;...
        Suw2,  Svw2,   Sww];
    
    [L,D]=ldl(S,'lower'); % a LDL decomposition is applied this time
    G = L*sqrt(D);
    A(ii,:)= G*exp(1i*2*pi*randPhase);
    if ii==2,    fprintf(['Expected computation time: From ',num2str(round(min(toc(tStart))*Nfreq/2)),' to ',num2str(round(min(toc(tStart))*Nfreq)),' seconds \n']); end
end

%% Apply IFFT

Nu = [A(1:Nfreq,:) ; real(A(Nfreq,:)); conj(flipud(A(2:Nfreq,:)))];
speed=real(ifft(Nu).*sqrt(Nfreq./(dt)));

u = speed(:,1:Nm);
v = speed(:,Nm+1:2*Nm);
w = -speed(:,2*Nm+1:3*Nm);


[~,indmax] = min(abs(t-tmax));
u = u(1:indmax,:)'-mean(u(1:indmax,:)',2);
v = v(1:indmax,:)'-mean(v(1:indmax,:)',2);
w = w(1:indmax,:)'-mean(w(1:indmax,:)',2);
t = t(1:indmax);




%% Nested functions
    function [coh] = cohDavenport(meanU,dy,dz,f,Cy,Cz)
        
        % lateral separation
        ay = Cy(1).*dy;
        % vertical separation
        az = Cz(1).*dz;
        % Combine them into the coherence matrix for lateral and vertical
        % separations
        coh = exp(-sqrt(ay.^2+az.^2).*f./meanU);
    end
    function [data] = importfile(filename, startRow, endRow,formatSpec,delimiter)
        % GOAL
        % extract data from txt files or excel files to store them into matrix
        
        %                               INPUT
        %  filename:
        %           type: string with extension.
        %           definition: name of the file whose information are extracted
        %  startRow:
        %           type: integer
        %           definition: first row of the file to be read
        %  endRow:
        %           type: integer
        %           definition: last row of the file to be read
        %  formatSpec:
        %           type: string
        %           definition: format of the data to be read.
        %           e.g. : formatSpec = ['%f%f%f%s%s%*[^\n]'];
        %  delimiter:
        %           type: string
        %           definition: symbol used to delimite the columns in the file
        %           e.g. : delimiter = '\t'; or delimiter = ',';
        
        %                               OUTPUT
        %  data:
        %           type: matrix [N1 x N2] where N1 and N2 are integers.
        %           N1 = endRow-startRow or total number of rows if enRow is not
        %           specified
        %           N2 = defined by formatSpec.
        %           definition: extracted data from the file, stored as a matrix
        % Copyright (C) Etienne Cheynet, 2015.
        % last modification: 27/01/2015 10:56
        % Initialize variables.
        if nargin<=2
            startRow = 1;
            endRow = inf;
        end
        if endRow == [],
            endRow = inf;
        end
        %  open the text file.
        fileID = fopen(filename,'r');
        %Read columns of data according to format string.
        dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
        for block=2:length(startRow)
            frewind(fileID);
            dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
            for col=1:length(dataArray)
                dataArray{col} = [dataArray{col};dataArrayBlock{col}];
            end
        end
        % Close the text file.
        fclose(fileID);
        %% Create output variable
        dataArray = cellfun(@(x) num2cell(x), dataArray, 'UniformOutput', false);
        data = [dataArray{1:end}];
    end
    function [S] = VonKarmanSpectrum(f,meanU,stdV,L,component)
        % ---------------------------------------------
        % INPUT
        % f: float; frequency is [1 x 1]
        % V: float; Mean wind speed Normal to the deck is [1x1]
        % std_speed : float; std of speed is [1 x 1]
        % L =  float; turbulence length scales is [1x1]
        % stdV = ; float; std of wind velocity fluctuations [1x1]
        % component : string; is 'u' or 'w'
        % ---------------------------------------------
        % OUTPUT
        % Sv: float; [1x1] value of Spectrum for a given frequency
        % ---------------------------------------------
        %calculation of S/std^2
        fr = L.*f./meanU;
        %         f0 = logspace(-4,3,100);
        %         fr0 = L.*f0./meanU(1);
        
        if strcmpi(component,'u')
            S = 4*fr./(1+71.*fr.^2).^(5/6);
            %             S0 = 4*fr0./(1+71.*fr0.^2).^(5/6);
        elseif strcmpi(component,'v')
            S=  4*fr.*(1+755*fr.^2)./(1+283.*fr.^2).^(11/6);
            %             S0=  4*fr.*(1+755*fr0.^2)./(1+283.*fr0.^2).^(11/6);
            
        elseif strcmpi(component,'w')
            S=  4*fr.*(1+755*fr.^2)./(1+283.*fr.^2).^(11/6);
            %             S0=  4*fr.*(1+755*fr0.^2)./(1+283.*fr0.^2).^(11/6);
        else
            error('component unknown \n\n')
        end
        S = S./f.*stdV.^2;
        %         coeffcorr = trapz(f0,S0./f0);
        %         S = S./coeffcorr.*stdV.^2;
        
        
    end
    function [S] = KaimalSpectrum(f,meanU,u_star,z,component)
        % ---------------------------------------------
        % INPUT
        % f: double; frequency is [1 x 1]
        % meanU: double; Mean wind speed [Nyy*Nzzx1]
        % u_star : double; friction velocity
        % z =  double; turbulence length scales is [Nyy*Nzzx1]
        % component : string; is 'u', 'v' or 'w'
        % ---------------------------------------------
        % OUTPUT
        % Sv: double; [Nyy*Nzzx1] value of Spectrum for a given frequency
        % ---------------------------------------------
        % ---------------------------------------------
        %%
        fr = f.*z(:)./meanU(:);
        if strcmpi(component,'u')
            S = u_star.^2./f.*(102.*fr)./(1+33.*fr).^(5/3);
        elseif strcmpi(component,'v')
            S = u_star.^2./f.*(17.*fr)./(1+9.5.*fr).^(5/3);
        elseif strcmpi(component,'w')
            S = u_star.^2./f.*(2.*fr)./(1+5.*(fr).^(5/3));
        elseif strcmpi(component,'uw')
            S = u_star.^2./f.*-14.*fr./(1+10.5*fr).^(7/3); % corrected (by me) Kaimal cross-spectrum model (NOT normalized)
        else
            fprintf(' spectrum type is unknown \n')
            return
        end
    end
    function [nodes,Nm] = createGrid(geometry,filename)
        if isempty(geometry)
            [data0] = importfile(filename, 40,45,['%*s%f%*[^\n]'],'\t');
            data0=cell2mat(data0);
            if any(data0(1:2)<0), error('Nyy and Nzz must be positiv'); end
            Nyy = data0(1);
            Nzz	= data0(2);
            Zmin= data0(3);
            Zmax= data0(4);
            Ymin= data0(5);
            Ymax= data0(6);
            
            if Zmin>Zmax
                warning('Zmin > Zmax. Their values have been switched');
                dummy = Zmax;
                Zmax = Zmin;
                Zmin = dummy;
                clear dummy
            end
            
            if Ymin>Ymax
                warning('Ymin > Ymax. Their values have been switched');
                dummy = Ymax;
                Ymax = Ymin;
                Ymin = dummy;
                clear dummy
            end
            
            % Check compatibility between grid and node number
            if and(Ymin==Ymax,Nyy>1)
                warning('Ymin = Ymax but Nyy > 1, Nyy is set to 1')
                Nyy=1;
            end
            
            if and(Zmin==Zmax,Nzz>1)
                warning('Zmin = Zmax but Nzz > 1, Nzz is set to 1')
                Nzz=1;
            end
            
            %% Create the grid
            
            y = linspace(Ymin,Ymax,Nyy);
            z = linspace(Zmin,Zmax,Nzz);
            [Z,Y] = meshgrid(z,y);
            Nm = numel(Y(:)); % number of nodes in the grid, equal to Nyy*Nzz
            %% Wind profile
            if strcmpi(windProfile,'power')
                U= Uref.*(Z./zr).^(a);
            elseif strcmpi(windProfile,'log')
                k=0.4; % Von karman constant
                U = u_star./k.*log(z./z0);
            else
                error('wind profile selected is unknown\n')
            end
            %% Create the structure nodes
            nodes.U = U;
            nodes.Y = Y;
            nodes.Z = Z;
            % Affect one name and to each nodes
            for jj=1:Nm,    nodes.name{jj} = strcat('N',num2str(jj));end
            nodes.U = U(:);
        elseif exist(geometry,'file')
            [data0] = importfile(geometry,1,inf,['%f%f%*[^\n]'],'\t');
            data0=cell2mat(data0);
            if any(data0(1:2)<0), error('Nyy and Nzz must be positiv'); end
            nNodes = numel(data0(:,1));
            Zmin= min(data0(:,2));
            Zmax= max(data0(:,2));
            Ymin= min(data0(:,1));
            Ymax= max(data0(:,1));
            y = data0(:,1);
            z = data0(:,2);
            Nm = numel(y(:)); % number of nodes in the grid, equal to Nyy*Nzz
            % Check compatibility between grid and node number
            if and(Ymin==Ymax,nNodes>1),
                error('Ymin = Ymax but Nyy > 1, please correct your geometry')
            end
            if and(Zmin==Zmax,nNodes>1),
                error('Zmin = Zmax but Nzz > 1, please correct your geometry')
            end
            
            % Create the grid
            clear nodes
            nodes.Y = y;
            nodes.Z = z;
            
            % Wind profile
            if strcmp(windProfile,'power'),
                meanU= Uref.*(nodes.Z./zr).^(a);
            elseif strcmp(windProfile,'log'),
                kappa=0.4; % Von karman constant
                meanU = u_star./kappa.*log(z./z0);
            else
                error('wind profile selected is unknown\n')
            end
            nodes.U = meanU(:);
            for jj=1:Nm,    nodes.name{jj} = strcat('N',num2str(jj));end
            
        else
            error('Input "geometry" is unknown')
        end
        
        
    end
end

