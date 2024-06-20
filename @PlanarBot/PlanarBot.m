classdef PlanarBot    
    properties
        kin
    end
    
    methods
        function obj = PlanarBot()
            ex = [1;0;0];
            ey = [0;1;0];
            ez = [0;0;1];
            zv = [0;0;0];

            obj.kin.H = [ez ez];
            obj.kin.P = [zv ex ex];
            obj.kin.joint_type = [0 0];
        end
    end
end

