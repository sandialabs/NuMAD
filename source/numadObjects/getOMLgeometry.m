function coords=getOMLgeometry(profile,naturaloffset,maxThicknessCordLocation, chordoffset,chordLength,rotorspin,degreestwist,sweep,prebend,spanLocation)
				
				x=profile(:,1);
				y=profile(:,2);
				
                if naturaloffset
                    x = x - maxThicknessCordLocation; 
                end
                x = x - chordoffset;   % apply chordwise offset
                x = x * chordLength * -1*rotorspin;  % scale by chord
                y = y * chordLength;                     % scale by chord
                twist = -1*rotorspin * degreestwist;
                % prepare for hgtransform rotate & translate
                coords(:,1) = cosd(twist) * x - sind(twist) * y;
                coords(:,2) = sind(twist) * x + cosd(twist) * y;
                coords(:,3) = zeros(size(x));
                coords(:,4) = ones(size(x));
                
                % use the generating line to translate and rotate the coordinates
                [sweep_rot, prebend_rot] = deal(0);
% jcb: This code, copied from NuMAD 2.0, causes each section to rotate out
% of plane so that its normal follows the generating line direction. Need
% to replace 'twistFlag' with '-1*obj.rotorspin' and calculate the slopes
% based on the available data. For now, default to parallel sections.
%                 if isequal(blade.PresweepRef.method,'normal')
%                     sweep_slope = ppval(blade.PresweepRef.dpp,sta.LocationZ);
%                     sweep_rot = atan(sweep_slope*twistFlag);
%                 end
%                 if isequal(blade.PrecurveRef.method,'normal')
%                     prebend_slope = ppval(blade.PrecurveRef.dpp,sta.LocationZ);
%                     prebend_rot = atan(-prebend_slope);
%                 end
                transX = -1*rotorspin*sweep;
                transY = prebend;
                transZ = spanLocation;
                R = makehgtform('yrotate',sweep_rot,'xrotate',prebend_rot);
                T = makehgtform('translate',transX,transY,transZ);
                coords = coords * R' * T';
				coords(:,4)=[];
end