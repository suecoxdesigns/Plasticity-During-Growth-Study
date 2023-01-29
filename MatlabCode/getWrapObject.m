function [wrap2,properties] = getWrapObject(model,bodyName,wrapObjectName)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
import org.opensim.modeling.*
        wrap = model.getBodySet().get(bodyName).getWrapObjectSet().get(wrapObjectName);
        wrap2 = WrapCylinder.safeDownCast(wrap);
        translation = wrap2.getPropertyByName('translation');
        rotation = wrap2.getPropertyByName('xyz_body_rotation');
        radius = wrap2.getPropertyByName('radius');
        length = wrap2.getPropertyByName('length');
        properties = struct('radius',radius,'rotation',rotation,'translation',translation,'length',length)
end

