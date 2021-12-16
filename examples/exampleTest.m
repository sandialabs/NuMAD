classdef exampleTest < matlab.unittest.TestCase
    
    properties
        numadDir = ''
        exampleDir = ''
        exceldDir = '' 
        OriginalDefault
    end
    
    methods (Access = 'public')        
        function obj = exampleTest            
            % Set NuMAD and Example directories
            obj.numadDir = fullfile(obj.exampleDir,'..');
            obj.exampleDir = fileparts(mfilename('fullpath'));
            obj.exceldDir = fullfile(obj.exampleDir,'ExcelToObject');            
            % Save the visibility state at construction
            obj.OriginalDefault = get(0,'DefaultFigureVisible');
            set(0,'DefaultFigureVisible','off')            
        end        
    end
    
    methods(TestMethodTeardown)        
        function closePlots(testCase)
            close all
            set(0,'DefaultFigureVisible',testCase.OriginalDefault);
        end        
    end
    
    methods(TestClassTeardown)        
        function cdToexampleDir(testCase)
            cd(testCase.exampleDir);
        end        
    end
    
    methods(Test)        
        function bladeObject(testCase)
            % Kelley: for now, I'm just testing that blade is a BladeDef object, but
            % the test will fail if any of these lines fail. 
            cd(fullfile(testCase.exceldDir))                        
            designFile = 'Excel2ObjectExample.xlsx';
            blade = xlsBlade(designFile);
            assert(isa(blade,'BladeDef'),'blade is not a BladeDef object.')
        end        
        function excelToObject(testCase)
            % % Write ExcelToObject test here
            cd(fullfile(testCase.exceldDir))                        
            % excel2NumadObject            
        end              
    end
end
