   function [S] = NPDSpectrum(f,U10,z)
        % ---------------------------------------------
        % INPUT
        % f: float; frequency is [1 x 1]
        % U10: float; Hourly mean wind speed at 10 m asl [1x1]
        % z =  double; height about the surface is [Nyy*Nzzx1]
        % ---------------------------------------------
        % OUTPUT
        % S: float; [1x1] value of Spectrum for a given frequency
        % ---------------------------------------------

        
        fr = 172*f.*(z(:)/10).^(2/3)*(U10./10).^(-3/4);
        n = 0.468;
        S = 3.2.*U10.^2 .* (z(:)/10).^(0.45) .* (1+fr.^n).^(-5/(3*n));
        
        
   end