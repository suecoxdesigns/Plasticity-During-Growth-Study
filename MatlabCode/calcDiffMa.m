function [diffMa]= calcDiffMa(d3,model,wrapObject,x0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import org.opensim.modeling.*

radius = x0(1);
translation = x0(2:4);
rotation = x0(5:end);

state=model.initSystem;
jointCoords= model.updCoordinateSet().get('ankle_flexion');
LG = model.getMuscles().get('LG');

%first adjust the wrap object
wrapObject.set_radius(radius);
wrapObject.upd_radius();
wrapObject.set_xyz_body_rotation(Vec3(rotation(1),rotation(2),rotation(3)));
wrapObject.upd_xyz_body_rotation();
wrapObject.set_xyz_body_rotation(Vec3(translation(1),translation(2),translation(3)));
wrapObject.upd_translation();

%then calculate the moment arm and difference from experimental
count=1;
for ang = [d3.maAngles]
    model.updCoordinateSet().get('ankle_flexion').setValue(state,ang);
    LgMa(count)=LG.computeMomentArm(state,jointCoords);
    count=count+1;
end

diffMa = sum((LgMa+[d3.ma]).^2);
end

