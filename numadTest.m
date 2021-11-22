%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Part of the SNL NuMAD Toolbox                    
%  Developed by Sandia National Laboratories Wind Energy Technologies 
%              See license.txt for disclaimer information             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = numadTest(options)
    % ``numadTest`` runs NuMAD continuous integration testing suite.
    %
    %   results = numadTest returns a matlab.unittest.TestResult object
    %
    %   results = numadTest(..., NAME1, VALUE1, NAME2, VALUE2, ...) returns
    %   a matlab.unittest.TestResult object, depending on the values of the 
    %   optional parameter name/value pairs. See Parameters below.
    %
    %   Parameters
    %   ----------
    %   'exampleTest'     Run tests for NuMAD examples. Default is true.
    %

    arguments
        options.exampleTest = true
    end
    
    import matlab.unittest.TestSuite
    import matlab.unittest.Test
    import matlab.unittest.TestRunner
    import matlab.unittest.plugins.DiagnosticsRecordingPlugin
    import matlab.unittest.plugins.CodeCoveragePlugin
    
    suites = Test.empty();
    
    if options.exampleTest
        suites = [suites TestSuite.fromFile('examples/exampleTest.m')];
    end
        
    % Create TestRunner
    runner = TestRunner.withTextOutput; % Contains TestRunProgressPlugin, DiagnosticsOutputPlugin
    runner.addPlugin(DiagnosticsRecordingPlugin);
    runner.addPlugin(CodeCoveragePlugin.forFolder('./source','IncludingSubfolders',true));
    
    % Run the tests
    results = runner.run(suites);
    results.table    
end
