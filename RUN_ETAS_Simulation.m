
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This simulator is simplified temporal-only version of AFTSimulator.m 
% Felzer, K. R., T. W. Becker, R. E. Abercrombie, G. Ekstrom, and J. R.
% Rice, Triggering of the 1999 Mw 7.1 Hector Mine earthquake by aftershocks
% of the 1992 Mw 7.3 Landers earthquake, J. Geophys. Res., 107, 2190,
% doi:10.1029/2001JB000911, 2002. 
%
% Everything except temporal ETAS simulation is removed from original code 
% by Kyungjae Im.
% The simulation result will be saved as "CatalogETAS". 
% The parameter inversion of saved catalog can be done by 
% RUN_MLE_Inversion.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
rng shuffle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Prameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mu0 = 0.01; % expected events per day
C=0.001;
P=1.20;
K0=0.004;
Alpha=log(10);

LowerLimitMagnitude = 2.5; % magnitud lower limit 
UpperLimitMagnitude = 7;  % magnitud upperlimit limit 
SimulationDays = 50*365;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AFT Simulator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Felzer et al., 2002 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Simplified and Modified by Im %%%%%%%%%%%%%%%%%%%%%

AveNumOfBkgndQuakes=Mu0*SimulationDays;
NumOfBkgndQuakesInThisSim = poissinv(rand(1),AveNumOfBkgndQuakes); 
TimesOfBkgndQuakes = rand(NumOfBkgndQuakesInThisSim,1) *(SimulationDays) ;
TimesOfBkgndQuakes = sort(TimesOfBkgndQuakes);
MagsOfBkgndQuakes ... 
   = GetMags(NumOfBkgndQuakesInThisSim,LowerLimitMagnitude, UpperLimitMagnitude);
BkgndQuakes = [TimesOfBkgndQuakes, MagsOfBkgndQuakes];
Catalog = BkgndQuakes;

if(P ~= 1)
    AftProductivity = (K0/(1-P))*((SimulationDays+C)^(1-P) - C^(1-P));
else
    AftProductivity = K0*(log(SimulationDays+C) - log(C));
end  
        
NewQuakes = BkgndQuakes;
NumOfNewQuakes = length(BkgndQuakes);
NumOfIterations = 0;


while NumOfNewQuakes>0
    TimeOfNewQuakes = NewQuakes(:,1);
    MagsOfNewQuakes = NewQuakes(:,2);
    
    % Calculate how many aftershock expected for each new quakes
    AveNumOfAftsOfEachQuake ... 
       = AftProductivity.*exp(Alpha*(MagsOfNewQuakes - LowerLimitMagnitude));
   
   % Pick actual number of new quakes as a random poisson process
    NumOfAftsOfEachQuake ... 
        = poissinv(rand(length(AveNumOfAftsOfEachQuake),1), ... 
            AveNumOfAftsOfEachQuake);
        
    DaughterCount=0;
    DaughterTime = 0;
    DaughterMagnitude = 0;
    for ParentIndex=1:length(MagsOfNewQuakes)
        for DaughterIndex=1:NumOfAftsOfEachQuake(ParentIndex)

            TimesBetweenAftsAndParents ... 
                = GetTimesBetweenAftsAndParents(P,C, 0,SimulationDays);
            
            if TimeOfNewQuakes(ParentIndex) + TimesBetweenAftsAndParents < SimulationDays
                DaughterCount=DaughterCount+1;
                DaughterTime(DaughterCount) = TimeOfNewQuakes(ParentIndex) + TimesBetweenAftsAndParents;
                DaughterMagnitude(DaughterCount) = GetMags(1,LowerLimitMagnitude, UpperLimitMagnitude);
            end
        end
    end
    DaughterCount
    if DaughterCount>0
    NewQuakes = [DaughterTime', DaughterMagnitude'];
    Catalog = [Catalog; NewQuakes];
    end
    NumOfNewQuakes = DaughterCount;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save the Catalog %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Day = Catalog(:,1);
Magnitude = Catalog(:,2);

figure(1)
hold on
scatter(Day, Magnitude)
xlabel("Time (Day)")
ylabel("Magnitude")
xlim([0, SimulationDays])
    

ThetaOriginal=[Mu0, C, P, K0, Alpha]
M0_Original = LowerLimitMagnitude
clearvars -except Day Magnitude ThetaOriginal M0_Original
save CatalogETAS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AFT Simulator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Felzer et al., 2002 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------------%
%----------------------------Function GetMags---------------------------%
%-----------------------------------------------------------------------%

function Mags = GetMags(NumOfAfts,MinMagOfActiveQuakes,MaxMagOfSimAft)

Mags = MinMagOfActiveQuakes - log10(rand(NumOfAfts,1));

MagsTooLarge = find(Mags > MaxMagOfSimAft);
NumOfMagsTooLarge = length(MagsTooLarge);

for i=1:NumOfMagsTooLarge
    
    while Mags(MagsTooLarge(i))>MaxMagOfSimAft
        Mags(MagsTooLarge(i)) = MinMagOfActiveQuakes - log10(rand(1));
    end
end
end   %GetMags


%-----------------------------------------------------------------------%
%---------------Function GetTimesBetweenAftsAndParents------------------%
%-----------------------------------------------------------------------%

function TimesBetweenAftsAndParents ... 
        = GetTimesBetweenAftsAndParents(P,CX,tmin,tmax);
%This function generates random points from the function (t+c)^-p, between
%tmin and tmax.  0 may be entered for tmin.

r = rand(1);

if(P~=1)
    
    a1 = (tmax + CX)^(1-P);
    a2 = (tmin + CX)^(1-P);
    a3 = r*a1 + (1-r)*a2;
    TimesBetweenAftsAndParents = a3.^(1/(1-P)) - CX;
    
else
    
    a1 = log(tmax+CX);
    a2 = log(tmin + CX);
    a3 = r*a1 + (1-r)*a2;
    TimesBetweenAftsAndParents = exp(a3) - CX;
     
end

end   %TimesBetweenAftsAndParents


