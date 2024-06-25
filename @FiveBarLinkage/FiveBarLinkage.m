classdef FiveBarLinkage    
    properties
        armA
        armB
    end
    
    methods
        function obj = FiveBarLinkage()
            obj.armA = PlanarBot;
            obj.armB = PlanarBot;

            ex = [1;0;0];
            ey = [0;1;0];
            ez = [0;0;1];
            zv = [0;0;0];

            obj.armA.kin.P = [zv ex ex];
            obj.armB.kin.P = [ex ex ex];
            obj.armA.kin.H = [ez ez];
            obj.armB.kin.H = [ez ez];
            obj.armA.kin.joint_type = [0 0];
            obj.armB.kin.joint_type = [0 0];
        end
    end
end

