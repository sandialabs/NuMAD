function gM=setIEC5(analysis,env,temp,mfg,calca,loada,varargin)
    if ~isempty(varargin)
        calcb = varargin{1};
        loadb = varargin{2};

    end
    gM0=1.2;
    
    switch analysis
		case 'LU'
			%**************************************************************************
			%Laminate Ultimate
			switch env
				case 'noEnv'; gM1=1.2;
				case 'env'; gM1=1.0;
				otherwise warning('Laminate ultimate environmental degradation not specified correctly');
			end
			
			switch temp
				case 'noTemp'; gM2=1.1;
				case 'temp'; gM2=1.0;
				otherwise warning('Laminate ultimate temperature effects not specified correctly');
			end
			
			switch mfg
				case 'noFlaw'; gM3=1.3;
				case 'basicFlaw'; gM3=1.1;
				case 'fullFlaw'; gM3=1.0;
				otherwise warning('Laminate ultimate manufacturing effects not specified correctly');
			end
			
			switch calca
				case 'noTest'; gM4a=1.2;
				case 'test'; gM4a=1.0;
				otherwise warning('Laminate ultimate calculation method and/or validation not specified correctly');
			end
			
			gM4b=1;
			
			switch loada
				case '4_directions'; gM5a=1.2;
				case '12_directions'; gM5a=1.0;
				otherwise warning('Laminate ultimate load characterization not specified correctly');
			end
			
			gM5b=1;
						
		case 'LF'
			%**************************************************************************
			%Laminate Fatigue
			switch env
                case 'noEnv'; gM1=1.1;
                case 'env'; gM1=1.0;
				otherwise warning('Laminate fatigue environmental degradation not specified correctly');
			end
			
			gM2=1.0;
			
			switch mfg
				case 'noFlaw'; gM3=1.3;
				case 'basicFlaw'; gM3=1.1;
				case 'fullFlaw'; gM3=1.0;
				otherwise warning('Laminate fatigue manufacturing effects not specified correctly');
			end
			
			switch calca
				case 'noTest'; gM4a=1.2;
				case 'test'; gM4a=1.0;
				otherwise warning('Laminate fatigue model validation not specified correctly');
			end
			
			switch calcb
				case 'assumedSlope'; gM4b=1.2;
				case 'measuredSlope'; gM4b=1.1;
				case 'full'; gM4b=1.0;
				otherwise warning('Laminate fatigue calculation method and/or validation not specified correctly');
			end
			
			switch loada
				case '2_directions'; gM5a=1.2;
				case '6_directions'; gM5a=1.0;
				otherwise warning('Laminate fatigue calculation method and/or validation not specified correctly');
			end
			
			switch loadb
				case 'DEL'; gM5b=1.2;
				case 'Markov'; gM5b=1.0;
				otherwise warning('Laminate fatigue load characterization not specified correctly');
			end

		case 'IFF'
			%**************************************************************************
			%Inter Fiber Failure
			gM1=1.0;
			
			switch temp	
				case 'noTemp'; gM2=1.1;
				case 'temp'; gM2=1.0;
			end
			
			gM3=1.1;

			switch calca	
				case 'FEA'; gM4a=1.1;
				case 'analytical'; gM4a=1.0;
			end
			
			gM4b=1.0;
			
			switch loada
				case '4_directions'; gM5a=1.1;
				case '12_directions'; gM5a=1.0;
			end

			gM5b=1.0;
			
		case 'SC'
			%**************************************************************************
			%Sandwhich Core Ultimate
			switch env
				case 'open'; gM1=1.3;
				case 'closed'; gM1=1.1;
				case 'env'; gM1=1.0;
			end
			
			switch temp
				case 'litNoTemp'; gM2=1.2;
				case 'noTemp'; gM2=1.1;
				case 'temp'; gM2=1.0;
			end
			
			switch mfg
				case 'noFlaw'; gM3=1.3;
				case 'basicFlaw'; gM3=1.1;
				case 'fullFlaw'; gM3=1.0;
			end

			switch calca
				case 'analyticalNoTest'; gM4a=1.35;
				case 'FEAShellNoTest'; gM4a=1.20;
				case 'FEABrickNoTest'; gM4a=1.20;
				case 'FEAShellTest'; gM4a=1.0;
				case 'FEABrick3DTest'; gM4a=1.0;
			end
			
			gM4b=1.0;
			
			switch loada
				case '4_directions'; gM5a=1.2;
				case '12_directions'; gM5a=1.0;
			end

			gM5b=1.0;
			
		case 'GS'
			%**************************************************************************
			%Global Stability
			gM1=1.0;
				
			switch temp
% 				case 'litNoTemp'; gM2=1.2;
				case 'noTemp'; gM2=1.1;
				case 'temp'; gM2=1.0;
            end
			
            
            gM3 = 1.0;
% 			switch mfg
% 				case 'noFlaw'; gM3=1.3;
% 				case 'basicFlaw'; gM3=1.1;
% 				case 'fullFlaw'; gM3=1.0;
% 			end

            switch calca
                case 'analytical'; gM4a=1.4;
                case 'linearFEA'; gM4a=1.2;
                case 'nonlinearFEA'; gM4a=1.0;
            end
            switch calcb
				case 'not_validated'; gM4b=1.0;
				case 'validated'; gM4b=1.0;
			end
			
			gM4b=1.0;
			
			switch loada
				case '4_directions'; gM5a=1.2;
				case '12_directions'; gM5a=1.0;
			end

			gM5b=1.0;
		
		case 'LS'
			%**************************************************************************
			%Local Stability
			gM1=1.0;
				
			switch temp
				case 'noTemp'; gM2=1.1;
				case 'temp'; gM2=1.0;
			end
			
			gM3=1.0;

			switch calca
				case 'analyticalNoTest'; gM4a=1.35;
				case 'FEAShellNoTest'; gM4a=1.20;
				case 'FEABrickNoTest'; gM4a=1.20;
				case 'FEAShellTest'; gM4a=1.0;
				case 'FEABrick3DTest'; gM4a=1.0;
                otherwise warning('Laminate ultimate calculation method and/or validation not specified correctly');
			end
			
			gM4b=1.0;
			
			switch loada
				case '4_directions'; gM5a=1.2;
				case '12_directions'; gM5a=1.0;
			end

			gM5b=1.0;
	
		case 'BU'
			%**************************************************************************
			%Bonded Ultimate
			switch env
				case 'noEnv'; gM1=1.2;
				case 'env'; gM1=1.0;
				otherwise warning('Laminate ultimate environmental degradation not specified correctly');
			end
			
			switch temp
				case 'noTemp'; gM2=1.1;
				case 'temp'; gM2=1.0;
				otherwise warning('Laminate ultimate temperature effects not specified correctly');
			end
			
			switch mfg
				case 'noFlaw'; gM3=1.3;
				case 'basicFlaw'; gM3=1.1;
				case 'fullFlaw'; gM3=1.0;
				otherwise warning('Laminate ultimate manufacturing effects not specified correctly');
			end
			
			switch calca
				case 'analytical'; gM4a=2.0;
				case 'FEstressBased'; gM4a=1.3;
                case 'FE_FMbased'; gM4a=1.1;    
                case 'test'; gM4a=1.0;
				otherwise warning('Laminate ultimate calculation method and/or validation not specified correctly');
			end
			
			gM4b=1;
			
			gM5a=1;
			gM5b=1;

		case 'BF'
			%**************************************************************************
			%Bonded Fatigue
            switch env
				case 'noEnv'; gM1=1.1;
				case 'env'; gM1=1.0;
				otherwise warning('Laminate fatigue environmental degradation not specified correctly');
			end
			
			gM2=1.0;
			
			switch mfg
				case 'noFlaw'; gM3=1.3;
				case 'basicFlaw'; gM3=1.1;
				case 'fullFlaw'; gM3=1.0;
				otherwise warning('Laminate fatigue manufacturing effects not specified correctly');
			end
			
			switch calca
				case 'analytical'; gM4a=2.0;
				case 'FEstressBased'; gM4a=1.3;
                case 'FE_FMbased'; gM4a=1.1;    
                case 'test'; gM4a=1.0;
				otherwise warning('Laminate ultimate calculation method and/or validation not specified correctly');
			end
			
			
			switch loada
				case 'DEL'; gM5a=1.05;
				case 'Markov'; gM5a=1.0;
				otherwise warning('Laminate fatigue load characterization not specified correctly');
            end
            gM5b=1.0; 
            			
    end
	
	gM=gM0*gM1*gM2*gM3*gM4a*gM4b*gM5a*gM5b;

end