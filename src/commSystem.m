classdef commSystem
    % Digital Communication System
    
    properties
        KwonModTypes = ['4-pam', '4-qam', 'qpsk'];
        ModulationType = strings
    end
    
    methods
        function obj = commSystem(ModulationType)
            %commSystem Construct an instance of this class
            if (strfind(obj.KwonModTypes, ModulationType) >= 1)
                obj.ModulationType = ModulationType;
            else
                warning("Invalid modulation type. Default is 4-PAM");
                obj.ModulationType = "4-pam";
            end
        end
        
        function out = type(obj)
            %TYPE
            disp(obj.ModulationType)
            out = obj.ModulationType;
        end
        
        function out = tx(obj)
            %TYPE
            disp(obj.ModulationType)
            out = obj.ModulationType;
        end
        
        function out = channel(obj)
            %TYPE
            disp(obj.ModulationType)
            out = obj.ModulationType;
        end
        
        function out = rx(obj)
            %TYPE
            disp(obj.ModulationType)
            out = obj.ModulationType;
        end
    end
end

