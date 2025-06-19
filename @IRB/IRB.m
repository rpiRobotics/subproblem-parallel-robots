classdef IRB
    properties
        kin
    end

    methods
        function obj = IRB()
            obj.kin = hardcoded_IK_setups.IRB_6640.get_kin;
        end
    end
end