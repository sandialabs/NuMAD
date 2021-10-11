classdef fileReader < handle
    properties
        fileID
%         fmt_num = '  %-8g';
%         fmt_str = '  %-8s';
%         fmt_lst = '%-8g ';
    end
    
    properties (Hidden)
        lineBuffer = '';
    end
    
    methods
        function self = fileReader(filename)
            self.fileID = fopen(filename,'rt');
            assert(self.fileID ~= -1,'fileReader:cannotOpenFile',...
                'Could not open file %s.',filename);
        end
        
        function delete(self)
            if (self.fileID ~= -1)
                fclose(self.fileID);
            end
        end
        
        function discardLine(self)
            if isempty(self.lineBuffer)
                fgetl(self.fileID);
            else
                self.lineBuffer = '';
            end
        end
        
        function s = getLine(self)
            if isempty(self.lineBuffer)
                s = fgetl(self.fileID);
            else
                s = self.lineBuffer;
                self.lineBuffer = '';
            end
        end
        
        function pushLine(self,s)
            n = length(s)+2;
            fseek(self.fileID,-n,'cof');
%             if isempty(self.lineBuffer)
%                 self.lineBuffer = s;
%             else
%                 error('cannot pushLine: lineBuffer not empty')
%             end
        end
        
        function s = str(self)
            s = sscanf(self.getLine,'%s',[1 1]);
        end
        
        function val = num(self)
            val = sscanf(self.getLine,'%g',[1 1]);
        end
        
%         function row = fmt(self,fmtstr)
%             row = fscanf(self.fileID,'%g',[1 1]); 
%         end
    end
    
end
    