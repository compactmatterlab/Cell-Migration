%Fiber based cell migration code - Single cell persistent random walk
%code written by Ben Yeoman
%created 01/21/2016
%edited 11/28/2017

clear all
c = parcluster;
c.NumWorkers = 20;
saveProfile(c);
parpool(20);

%File settings or data storage
filename = strcat('RT2_',date);
rng(datenum(date));%......................................For repeatability
tic

%Model Definitions
runTime = 12*3600;%...........................Total runtime for model (sec)
dt = 2;%.......................................Size of each time step (sec)
elmtSize = 0.01;%............Size of binding site element (um) - max is 0.3
samples = 10;%..................................Sample size to get averages
xPnts = 1;%.............................Number of parameter data points
yPnts = 20;%.............................Number of parameter data points



%Initiate matrices for MSD calculation
time = ((0:1+runTime/dt-1)');%...................Time at each cell position
tMax = (floor((1+runTime/dt)/8));%.............Maximum time for MSD fitting
tMsd = dt*((1:tMax)')/3600;%.....................Time steps for MSD fitting

%Min and max values for each parameter
pr = [2 6;%.......................................Gel concentration (mg/ml)
    0.0005 0.002;%..............................Fiber density (fibers/um^3)
    1/16 1/2;%.........................Pseudopod extension frequency (s^-1)
    0 0.8;%............................................Fiber alignment (AI)
    5 12];%........................Binding Site (RGD peptides/tropocollagen)

%Parameter ranges
CGel = (linspace(pr(1,1),pr(1,2),xPnts));%.......Gel concentration range
PFiber = (linspace(pr(2,1),pr(2,2),xPnts));%.........Fiber Density range
Tsearch = (linspace(pr(3,1),pr(3,2),xPnts));%..........Search time range
FiberAI = (linspace(pr(4,1),pr(4,2),xPnts));%......Fiber Alignment range
rgdPeptides = [5.2,5.75,7]; %(linspace(pr(5,1),pr(5,2),yPnts));%.....RGD Peptides range

%Constant Cell Attributes
Vpseudo = 0.45;%...................Velecity of extending pseudopod (um/sec)
Lsearch = 0.3;%...............Length of pseudopod tip in contact with fiber
Fmax = 10e-9;%........................Max force provided by single cell (N)
Dcell = 15;%...........................................Cell's diameter (um)
Vcell = (4/3)*pi*(Dcell/2)^3;%.........................Cell's volume (um^3)
dRecp = 0.5*10^6;%..............................Integrin receptors per cell 
dR = dRecp/saellipsoid(15/2);%...Integrin receptor density (receptors/um^2)
tC = Lsearch/Vpseudo;%............................Binding site contact time
Bhalf = 200;%....................Number of bound ligands at half max. force
Bmin = 100;%...............Min bound ligands needed for pseudopod outgrowth 
Bmax = 530;%.............Min bound ligands needed for pseudopod contraction 
kcell = 1e-8;%.........................................Cell spring constant

%Constant ECM Attributes
Ltropo = 0.3;%.................................Length of tropocollagen (um)
Dtropo = 0.0015;%............................Diameter of tropocollagen (um)
Lfiber = 75;%.............................................Fiber length (um)
otAngle = 0;%...................Average angle between ECM fibers and x-axis
vRef = [1,0,0];%.........................Fiber orientation referance vector
eta = 1e-11;%............................Extracellular viscosity (N*s/um^2)
kI = 0.25e-9;%........................................Bond stiffness (N/um)
koff = 5;%.............................Bond dissociation rate .1-100 (s^-1)
kon = 1.4e-3;%...........................Bond association rate  (um^2/s^-1)
kB = 1.38e-17;%.........................Boltzmann's contant (um^2*kg/s^2*K)
T = 311;%...................................................Temperature (K)

%Output Parameter Preallocations
avgVel = zeros(xPnts,yPnts,samples);%.......................Average cell speed
LpR = zeros(xPnts,yPnts,samples);%.................Persistence length by <R^2>
D = zeros(xPnts,yPnts,samples);%...................Random motility coefficient
alp = zeros(xPnts,yPnts,samples);%.......................................Alpha
R1 = zeros(xPnts,yPnts,samples);%.....................Velocity Goodness of fit
R2 = zeros(xPnts,yPnts,samples);%..........................MSD Goodness of fit
R3 = zeros(xPnts,yPnts,samples);%..........................LpR Goodness of fit

polarize = 0.2;
    
%=========================================================================%
%Parametric Sweep
%=========================================================================%
parfor xd=1:samples
    
    
    p = [3.7,0.001,1/16,0,rgdPeptides(xd)]; 
    
    RGD = p(5);%.......................................RGD Peptides/monomer
    fiberDens = p(2);%.......................................Fiber density
    C_gel =  p(1);%..................................Collegen concentration
    Fiber_D = sqrt((C_gel*6.022e20*(Ltropo+0.067)*Dtropo^2)/...
    (1e12*fiberDens*805000*0.7*Lfiber*0.9));%.......Average fiber diameter
    dL = (RGD*.7*.9)/((Ltropo+0.067)*Dtropo);%.............RGD ligand density
    
    avgSearchTime = 1/p(3);%.......................Pseudopod extension time
    
    AI = p(4);%.............................................Fiber Alignment
    sigma = sqrt(-1642*reallog(AI));%.Angle deviation from reference vector   
    if isinf(sigma)
        sigma = 1e10;
    end
    crslnkPerL = ((4.439*sigma)/(sigma+4.55))/8;%..........Crosslinks per fiber 
    crslnkPerF = Lfiber*crslnkPerL;
    crslnkDens = fiberDens*crslnkPerF;%..................Crosslink density
    
    %Equation from Lin et. al. 2015
    if crslnkPerF < 3.5
        gel_stiffness = 1039.9*crslnkPerF-1992.9;
    else
        gel_stiffness = 5247.9*crslnkPerF-16274;
    end%.................................................Gel stiffness (Pa)     
    kecm = (1/crslnkDens^(1/3))*...
        gel_stiffness/1e12;%.....................ECM spring constant (N/um)
    


%=========================================================================%
%Sample Sweep - repeats simulation for number of given samples
%=========================================================================%
    for yd=1:yPnts
        disp([yd xd])
        
        tic
        
        rt = linspace(12*3600,240*3600,yPnts);
        runTime = rt(yd);
        time = ((0:1+runTime/dt-1)');%...................Time at each cell position
        tMax = (floor((1+runTime/dt)/4));%.............Maximum time for MSD fitting
        tMsd = dt*((1:tMax)')/3600;%.....................Time steps for MSD fitting
        
        %=================================================================%
        %Position and Phase Matrices
        %=================================================================%
        cellPos = zeros(1,3);%.........................Cell position matrix
        phaseRT = zeros;%....................Pseudopod phase runtime matrix
        pseudoVect = zeros(1,3);.................Pseudopod directoin vector
        
        %=================================================================%
        %Tracking variables
        %=================================================================%
        Lpseudo = 0;%............................Tracks length of pseudopod
        newAngle = 1;%.....................Tracks when cell finds new fiber
        retracting = 1;%............................Tracks retracting phase
        outgrowth = 0;%..............................Tracks outgrowth phase
        contracting = 0;%..........................Tracks contracting phase
        t = 0;%........................Tracks pseudopod extension frequency
        pos = 1;%....................Tracks cell position at each time step
        
        %=================================================================%
        %Temporary Variables
        %=================================================================%
        fiberDir = 0;%.........................Unit vector of guiding fiber
        searchTime = 0;%...........................Pseudopod extension time
        numElmts = 0;%...........................Discrete elements of fiber
        Bsites = 0;%..............................Distributed binding sites
        searchLength = 0;%.............Number of elements searched per step
        bonds = 0;%.........................................Filopodia bonds
        Bpseudo = 0;%...................Number bonds along entire pseudopod
        dMove = 0;%......................Distance cell moves in contraction
        localFibers = 0;%....................Number of fibers touching cell
        cell = 0;%.............................................Axes of cell
        dR = 0;%...........................................Receptor density
        polAngle = 0;........................................Polarity angle
        dl = 0;
        
        %=================================================================%
        %Generates Cell's Initial Position (+/- 0.5um from origin)
        %=================================================================%
        cellPos(pos,:) = 2*rand(1,3)-1;%...................Initial position
        
        %Calculates Fiber's Angle from Reference Vector
        if sigma == 1e10
            fiberAngle = (2*rand-1)*360;%...........Randomly aligned fibers
        else
            fiberAngle = (180*(randi(2)-1))+...
                normrnd(otAngle,sigma);%.....................Aligned fibers
        end
        
        %Sets Initial Dinstance to Fiber Inersection
        Ltemp  = exprnd(1/crslnkPerL);
        if Ltemp <= Lsearch
            Ltemp = Lsearch;%............Limited to length of pseudopod tip
        end
        
        %Sets Initial Polarity
        vPrev = 2*rand(1,3)-1;
        vPrev = vPrev/norm(vPrev);
        
        %=================================================================%
        %Calculates path of cell
        %=================================================================%
        for i=dt:dt:runTime
            
            %Generates new fiber direction
            if newAngle

                %Generates local fibers in contact with cell
                localFibers = poissrnd(Vcell*fiberDens);
                numFibers = localFibers;
                if localFibers < 1 || outgrowth
                    localFibers = 1;
                end
                if numFibers < 1
                    numFibers = 1;
                end
                
                %Reset fiber, angle, and length matrices
                fiber = zeros(1,3);%.................Fiber direction matrix
                fb = zeros(1,3);
                angle = zeros(1,1);%.....................Fiber angle matrix
                acAng = zeros(1,1);%...............Fiber acute angle matrix
                searchLength = round(Lsearch/elmtSize);%.....Element length 
                L = zeros(localFibers,1);%.........Distance to intersection
                cs = zeros;
                sn = zeros;
                

                for j=1:localFibers
                    
                    %Calculates Distance to Next Crosslink
                    L(j) = exprnd(1/(crslnkDens)^(1/3));
                    if L(j) <= Lsearch
                        L(j)  = Lsearch;
                    end
                    
                    %Calculates random fiber angle
                    if sigma == 1e10
                        newDir = 2*rand(1,3)-1;%........New fiber direction
                        fiber(j,:) = newDir(1,:)/norm(newDir);%..New fibers
                    
                    %Calculates aligned fiber angle
                    else
                        fiberAngle = (180*(randi(2)-1))+...
                            normrnd(otAngle,sigma);
                        v1 = (L(j)*cosd(fiberAngle)*vRef)/norm(vRef);
                        vRand = 2*rand(1,3)-1;
                        xProd = (cross(vRef,vRand));
                        v2 =(L(j)*sind(fiberAngle)*xProd)/norm(xProd);
                        newDir = (v1 + v2);%............New fiber direction
                        fiber(j,:) = newDir(1,:)/norm(newDir);%..New fibers
                    end                    
                    
                    %Get fiber vectors in direction of previous direction
                    fb(j,:) = fiber(j,:);
                    an = acosd(dot(vPrev,fb(j,:))/(norm(vPrev)*norm(fb(j,:))));
                    if an > 90
                        fb = -fb;
                    end
                end
                
                %Elongation Vector
                ve = sum(fb,1)+vPrev/norm(vPrev);

                %Cell follows acute angle between ve and prev direction
                veAngle = acosd(dot(vPrev,ve)/(norm(vPrev)*norm(ve)));
                if veAngle > 90
                    ve = -ve;
                end


                %Get ratio of a0 to b0
                fbp = vertcat(vPrev,fb);%...........All fibers cell touches
                for j=1:size(fbp,1)
                    cs(j,1) = abs(dot(ve,fbp(j,:))/(norm(ve)*norm(fbp(j,:))));              
                    sn(j,1) = abs(sind(acosd(dot(ve,fbp(j,:))/(norm(ve)*...
                        norm(fbp(j,:)))))); 
                end
                a0 = sum(cs);   
                b0 = sum(sn);
                
                

                bondsTotal = localFibers*RGD*(0.7*(Dcell/2)/Ltropo)...
                    *(0.9*Fiber_D/Dtropo)*pi/2;%......Bonds at trailing end
                keq = (bondsTotal*kI*kecm)/(bondsTotal*kI+kecm);

                %Define axes of the cell
                b = ((3*Vcell*(kcell+keq))/(4*pi*((a0/b0)*keq+kcell)))^(1/3);
                c = b;
                a = (3*Vcell)/(4*pi*b^2);
                cell = [a,b,c];%........Cell axes 
                
                %Calculate cell sphericity
                A_cell = saellipsoid(cell);%..........Cell surface area 
                s = (pi^(1/3)*(6*Vcell)^(2/3))/A_cell;%.Cell sphericity
                dR = dRecp/A_cell;%...........Membrane integrin density
                polAngle = 180*s^2;......................Polarity angle

                %Fiber angles with respect to previous direction
                for j=1:localFibers
                    angle(j,1) = acosd(dot(fiber(j,:),ve)/...
                        (norm(fiber(j,:))*norm(ve)));  
                end
                
                %Gets acute angle for all fibers
                for j=1:size(angle,1)
                    if angle(j) > 90
                        acAng(j,1) = 180-angle(j);
                    else
                        acAng(j,1) = angle(j);
                    end
                end
                
                %Determines which fiber to follow based on cell's polarity
                cell_pol = polAngle*(2*rand-1);
                ang = abs(acAng-cell_pol);                                
                min_idx = ang == min(ang);
                fiberDir = fiber(min_idx,:);
                
                %Limits reversal of polarity if new pseudopod
                if retracting && angle(min_idx) > 90 
                    if rand >= polarize
                        fiberDir = -1*fiberDir;
                    end
                %Determines which way cell moves at crosslink
                elseif outgrowth && angle(min_idx) > 90
                    if rand >= 0.01
                        fiberDir = -1*fiberDir;
                    end
                end
                
                %Gets distance to next crosslink
                crossL = L(min_idx);%............Distance to next crosslink
                numElmts = floor(crossL/elmtSize);%......Elements to search 
                
                
                %Distributes # of adhesion sites in each element - assumes 
                %cell touches half of fiber surface
                Fdm = normrnd(Fiber_D,0.02);%......Randomize fiber diameter
                if (Fdm <= 0)
                    Fdm = 0.01;
                end
                dmax = max([dL,dR]);%.Greater of ligand or receptor density
                Pb = ((kon*dmax)/(kon*dmax+koff))*...
                    (1-exp(-(kon*dmax+koff)*tC));%......Binding probability 
                Bsites = poissrnd(Pb*RGD*0.7*0.9*(elmtSize/(Ltropo+.067))*(Fdm/Dtropo)*...
                    (pi/2),numElmts,1);%..........Distributed binding sites
                Bpseudo = Bpseudo + sum(Bsites(1:searchLength));
                
                %Sets time until new pseudopod extends
                if (retracting)
                    searchTime = exprnd(avgSearchTime);%.....Extension time        
                    if searchTime > 2*avgSearchTime
                        searchTime = 2*avgSearchTime;%..........Extension limit
                    end
                end
                
                %Reset values
                newAngle = 0;
                retracting = 0;
                outgrowth = 0;
            end
            
            %Calculates bonds formed for each time step
            if not(contracting)
                dl = 0;  
                while(searchLength < numElmts && dl < dt*Vpseudo)
                    searchLength = searchLength+1;
                    dl = dl+elmtSize;
                    bonds = sum(Bsites((searchLength-(29)):searchLength));
                    Bpseudo = Bpseudo+Bsites(searchLength);
                        
                    %Break if retracting or contracting phase
                    if bonds >= Bmax || bonds < Bmin
                        break
                    end                   
                end
            end
            
            %=============================================================%
            %RETRACTING PHASE - pseudopod will retract if weak bonds are 
            %formed, or if new pseudopod starts extending
            %=============================================================%
            if bonds < Bmin || (bonds < Bmax && t > searchTime)
                retracting = 1;
                phaseRT(pos,1) = 1;
                
                %Resets values
                Lpseudo = 0;.........................Reset pseudopod length
                t = 0;%................................Reset time searching
                bonds = 0;%.................Reset bonds at tip of pseudopod
                Bpseudo = 0;%...................Reset bonds along pseudopod
                pseudoVect = zeros(1,3);%.........Reset pseudopod direction
                newAngle = 1;%.............................Select new fiber
                cellPos(pos+1,:) = cellPos(pos,:);%.....Store cell position
                
            %=============================================================%
            %OUTGROWTH PHASE - pseudopod will continue growing as long as 
            %enough stable bonds form to allow for actin polymerization
            %=============================================================%
            elseif bonds >= Bmin && bonds < Bmax
                outgrowth = 1;
                phaseRT(pos,1) = 2;
                
                %Tracks distance and time pseudopod extends along 
                %current fiber
                pseudoVect = pseudoVect+dl*fiberDir;
                Lpseudo = norm(pseudoVect);%............Length of pseudopod
                t = t+dt;%...................Increment time spent searching 
                
                %Get new fiber if cell reaches cross-link
                if searchLength >= numElmts
                    newAngle = 1;
                end
                cellPos(pos+1,:) = cellPos(pos,:);%.....Store cell position
                
            %=============================================================%
            %CONTRACTING PHASE - pseudopod will contract when enough bonds 
            %are formed to withstand the acto-myosin contractile force
            %=============================================================%
            elseif bonds >= Bmax 
                phaseRT(pos,1) = 3;
                
                %Calculates final pseudopod direction, magnitude, 
                %and contractile force
                if not(contracting)
                    
                    contracting = 1;%.................Set contracting phase
                    Ltouching = mean(cell)/4;%Length of fibers touching cell
                    bondsTotal = localFibers*RGD*(0.7*Ltouching/Ltropo)...
                        *(0.9*Fiber_D/Dtropo)*pi/2;%Bonds at trailing end
                    
                    pseudoVect = pseudoVect+dl*fiberDir;
                    Lpseudo = norm(pseudoVect);%..Final length of pseudopod
                    
                    
                 
                    F0 = (Fmax*Bpseudo)/(Bhalf+Bpseudo);%....Adhesion force

                    kT_eff = (F0*kecm^2*Lpseudo^3)/...
                        (F0+kecm*Lpseudo)^2;%...Effective equivalent of kBT
                    fb = bondsTotal*((kecm*kI)/((kecm+kI)*koff))...
                        *exp(-(kT_eff)/(bondsTotal*kB*T));%Bond friction 
                    
                    a = mean([min(cell),median(cell)]);%Estimated minor axis
                    b = max(cell);%......................Major axis of cell
                    
                    if b < a
                        beta = b/a;%.........................Inverse sphericity
                        Kprime = ((4/3)*(beta^2-1))/(((2*beta^2-1)/...
                            sqrt(beta^2-1))*reallog(beta+sqrt(beta^2-1))-...
                            beta);%.....................Drag adjustment factor
                    else
                        beta = a/b;
                        Kprime = ((4/3)*(beta^2-1))/((beta*(beta^2-2)/...
                            sqrt(beta^2-1))*atan(sqrt(beta^2-1))+...
                            beta);%.....................Drag adjustment factor
                    end
                    
                    fv = 6*pi*eta*a*Kprime;%...............Viscous friction
                    if ~isreal(fv)
                        disp(fv)
                    end
                        
                    Vinst = (F0*(F0+kecm*Lpseudo)-F0^2*...
                       reallog(F0+kecm*Lpseudo))/((fb+fv)*kecm*Lpseudo)-...
                       (F0*(F0+kecm*0)-F0^2*reallog(F0+kecm*0))/...
                       ((fb+fv)*kecm*Lpseudo);%......Insantaneous velocity
                    dMove = Vinst*dt;%.................Distance moved in dt   
                    
                end
                
                
                
                dMove = min([Lpseudo,dMove]);%.........Eliminates overshoot
                Lpseudo = Lpseudo-dMove;%..Decrement Lspeudo each time step
                
                %Stores new cell position and direction
                cellPos(pos+1,:) = dMove*(pseudoVect./norm(pseudoVect))+...
                    cellPos(pos,:);
                vPrev = cellPos(pos+1,:)-cellPos(pos,:);
               
                
                %Resets values after cell is done contracting
                if Lpseudo <= 0
                    searchTime = exprnd(avgSearchTime);
                    if searchTime > 2*avgSearchTime
                        searchTime = 2*avgSearchTime;
                    end
                    Lpseudo = 0;
                    t = 0;
                    bonds = 0;
                    Bpseudo = 0;
                    pseudoVect = zeros(1,3);
                    contracting = 0;
                end
            end
            %=============================================================%  

    
            pos = pos+1;
        end%.................................................End simulation
        rw = cellPos;
        
%=========================================================================%
% Output Paramter Calculations
%=========================================================================%        
        
        %=================================================================%
        % Cell Velocity
        %=================================================================%
        w = 1;
        tm = (1:runTime/10:runTime)./3600;%..10 time points to get distance
        dst = zeros;
        pos = 1;
        for q=1:(size(rw,1))/10:size(rw,1)%........Gets distance at time tm
            Mag = zeros;
            for j=1:q
                if cellPos(j+1,:) ~= cellPos(j,:)
                    Mag(pos,1) = norm(cellPos(j+1,:)-cellPos(j,:));
                    pos = pos+1;
                end
            end
            dst(w,1) = sum(Mag);%.Gets distance traveled at each time point
            w = w+1;
        end
        try%....................Curve fit to get slope for average velocity
            [f,g] = fit(tm',dst,'poly1');
            c = coeffvalues(f);
            avgVel(xd,yd) = c(1);%............................Cell Speed
            R1(xd,yd) = g.rsquare;%......................Goodness of fit
        catch
            disp('Velocity Fitting Error')
            avgVel(xd,yd) = 0;
            R1(xd,yd) = 0;
        end
        
        
        
        %=================================================================%
        % Non-Overlapping MSD
        %=================================================================%
        msd = zeros;
        m = 1;
        for tau = 1:tMax
            
            %Gets nonoverlapping points for tau
            idx = mod(time,tau+1) == 0;
            rwIdx = rw(idx,:);
            
            dspl = rwIdx(2:end,1:3)-rwIdx(1:end-1,1:3);%...xyz displacement
            D_squared = sum(dspl.^2,2);%...............Squared displacement
            msd(m,1) = mean(D_squared);%................................MSD
            m=m+1;
        end
        
        try%.........Curve fit to get random motility coefficient and alpha
            [f,g] = fit(tMsd,msd,'power1');
            c = coeffvalues(f);
            D(xd,yd) = c(1);%................Random motility coefficient 
            alp(xd,yd) = c(2);%....................................Alpha
            R2(xd,yd) = g.rsquare;%......................Goodness of fit
        catch
            disp('MSD Fitting Error')
            D(xd,yd) = 0;
            alp(xd,yd) = 0;
            R2(xd,yd) = 0;
        end

        %=================================================================%
        %Persistence Length 
        %=================================================================%
        %Creates indices and variables
        nPos = zeros;
        nVect = zeros;
        nMag = zeros;
        avgCos = zeros;
        avgRsquared = zeros;
        contourLength = zeros;
        pos = 1;
        
        
        %Position, vector, and magnitude for each change in position
        for j=1:size(cellPos,1)-1
            if cellPos(j+1,:) ~= cellPos(j,:)
                nPos(pos,1:3) = cellPos(j,:);
                nVect(pos,1:3) = cellPos(j+1,:)-cellPos(j,:);
                nMag(pos,1) = norm(cellPos(j+1,:)-cellPos(j,:));
                pos = pos+1;
            end
        end
        nPos(pos,1:3) = cellPos(size(cellPos,1),1:3);
        
        xl = 60; %floor(sum(nMag)/4);%.....Contour length limit for Lp plots
        n = 0;%.........................................Tracks steps for L0
        
        %Loop for incresing L0
        for L0=1:xl/300:xl
            
            n = n + 1;%.............................Increments each L0 step
            endPos = 1;%.............................First direction vector
            m = 0;%....................Tracks cos(theta)/R for each L0 step
            stop = 0;%Breaks while loop if pos equals size of vector matrix
            totalCos = 0;
            totalR = 0;
            
            %Calculates average cos(theta)/R^2 for given L0
            while endPos < size(nVect,1)
                
                Length = 0;%........Resets L after moving step length of L0
                m = m + 1;%...............Increments steps for cos(theta)/R
                strt_pos = endPos;%....Gets starting position of L0 contour
                
                %Gets the contour position at end of length L0
                while Length <= L0
                    %Breaks loop at end of contour 
                    if endPos == size(nVect,1)
                        stop = 1;
                        break;
                    end
                    
                    Length = Length + nMag(endPos);%......Length of contour
                    endPos = endPos + 1;%...Gets end position of L0 contour
                end
                
                if stop == 0
                    %cos(theta) between vectors at start and end of contour
%                     cos_theta = dot(nVect(endPos-1,:),nVect(strt_pos,:))/...
%                         (nMag(endPos-1)*nMag(strt_pos));
%                     totalCos = totalCos + cos_theta;
                    
                    %Distance R between start and end points of contour
                    Rsquared = (norm(nPos(endPos,:) - nPos(strt_pos,:)))^2;
                    totalR = totalR + Rsquared;
                end
            end
            
            %avgCosStep = totalCos/m;
            avgRstep = totalR/m;
            
            %avgCos(n,1) = avgCosStep;%...................Average cos(theta) 
            avgRsquared(n,1) = avgRstep;%.................Average R squared 
            
            contourLength(n,1) = L0;%.......................Contour lengths
        end
        
%         ft = fittype('exp(-x/b1)');
%         try%............................Curve fit to get persistence length
%             [f,g] = fit(contourLength,avgCos,ft,'StartPoint',1);
%             c = coeffvalues(f);
%             LpC(par,samp) = c(1);%.......................Persistence length
%             R(par,samp) = g.rsquare;%......................Goodness of fit
%         catch
%             LpC(par,samp) = 0;
%             R(par,samp) = 0;
%         end
        
        ft = fittype('2*(b1^2)*(exp(-(x/b1))-1+x/b1)');
        try%............................Curve fit to get persistence length
            [f,g] = fit(contourLength,avgRsquared,ft,'StartPoint',1);
            c = coeffvalues(f);
            LpR(xd,yd) = c(1);%.......................Persistence length
            R3(xd,yd) = g.rsquare;%......................Goodness of fit
        catch
            LpR(xd,yd) = 0;
            R3(xd,yd) = 0;
        end
        
        %RandomWalk(:,:,par,samp) = rw;
%         retr(xd,yd,samp) = sum(phaseRT == 1);
%         outg(xd,yd,samp) = sum(phaseRT == 2);
%         cont(xd,yd,samp) = sum(phaseRT == 3);

    tclock(xd,yd) = toc; %[floor(toc/3600), mod(floor(toc/60),60), mod(toc,60)];

    end%....................................................End sample loop
    
end%........................................................End parfor Loop

v1 = avgVel(:,1:samples/2);
v2 = avgVel(:,samples/2+1:samples);

p1 = LpR(:,1:samples/2);
p2 = LpR(:,samples/2+1:samples);

d1 = D(:,1:samples/2);
d2 = D(:,samples/2+1:samples);

a1 = alp(:,1:samples/2);
a2 = alp(:,samples/2+1:samples);

save(filename,'avgVel','v1','v2','LpR','p1','p2','D','d1','d2','alp','a1','a2','R1','R2','R3')
%delete(gcp('nocreate'))
