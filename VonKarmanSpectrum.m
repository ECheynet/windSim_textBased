function [S] = VonKarmanSpectrum(f,meanU,L,component)
        % ---------------------------------------------
        % INPUT
        % f: float; frequency is [1 x 1]
        % meanU: float; Mean wind speed Normal to the deck is [1x1]
        % stdVel : float; std of speed is [1 x 1]
        % L =  float; turbulence length scales is [1x1]
%         stdV = ; float; std of wind velocity fluctuations [1x1]
        % component : string; is 'u' or 'w'
        % ---------------------------------------------
        % OUTPUT
        % Sv: float; [1x1] value of Spectrum for a given frequency
        % ---------------------------------------------
        % Von Karman coefficent
        % dimension of output
        S = zeros(size(meanU)); % is [Nzz,1]
        %calculation of S/std^2
        fr = L.*meanU.^(-1).*f;
        if strcmp(component,'u')
            S =  (4.*fr)./(1+70.7.*fr.^2).^(5/6);
%             S = S.*2.4.^2;
        elseif strcmp(component,'v') || strcmp(component,'w')
            S=  (4.*fr).*(1+754*fr.^2)./(1+283.*fr.^2).^(11/6);
        else
            fprintf('error: component unknown \n\n')
            return
        end
        
%         S = S.*stdVel.^2./f;
        
        
        
end
    

 
