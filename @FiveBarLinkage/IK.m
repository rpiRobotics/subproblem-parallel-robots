function [Q_A, Q_B, is_LS_A, is_LS_B] =  IK(obj, p_0T)

[Q_A, is_LS_A] = obj.armA.IK(p_0T);
[Q_B, is_LS_B] = obj.armB.IK(p_0T);
end