function [wrapObject] = adjustWrapSurface(radius,translation,rotation,wrapObject)
%This function sets and updates the properties of a wrap surface

wrapObject.set_radius(radius);
wrapObject.upd_radius()
wrapObject.set_rotation(Vec3(rotation));
wrapObject.upd_rotation()
wrapObject.set_translation(Vec3(translation));
wrapObject.upd_translation()
end

