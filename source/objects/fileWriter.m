classdef fileWriter < handle
    properties
        fileID
        fmt_num = '  %-8g';
        fmt_str = '  %-8s';
        fmt_lst = '%-8g ';
    end
    
    methods
        function self = fileWriter(filename)
            if exist('filename','var') && ~isempty(filename)
                self.fileID = fopen(filename,'wt');
                assert(self.fileID ~= -1,'fileWriter:cannotOpenFile',...
                    'Could not open file %s.',filename);
            else
                self.fileID = 1;
            end
        end
        
        function delete(self)
            if self.fileID ~= 1
                fclose(self.fileID);
            end
        end
        
        function text(self,str)
            fprintf(self.fileID,'%s\n',str);
        end
        
        function num(self,value,descrip)
            fprintf(self.fileID,self.fmt_num,value);
            fprintf(self.fileID,' %s\n',descrip);
        end
        
        function str(self,value,descrip)
            fprintf(self.fileID,self.fmt_str,value);
            fprintf(self.fileID,' %s\n',descrip);
        end
        
        function csv(self,values,descrip)
            str = sprintf('%g, ',values); %create list of numbers
            fprintf(self.fileID,self.fmt_str,str(1:end-2));
            fprintf(self.fileID,' %s\n',descrip);
        end
        
        function list(self,values,descrip)
            str = sprintf(self.fmt_lst,values); %create list of numbers
            fprintf(self.fileID,'%s',str);
            fprintf(self.fileID,' %s\n',descrip);
        end
        
        function qlist(self,values,descrip)
            str = sprintf(self.fmt_lst,values); %create list of numbers
            fprintf(self.fileID,'"%s"',str);
            fprintf(self.fileID,' %s\n',descrip);
        end
        
    end
    
end