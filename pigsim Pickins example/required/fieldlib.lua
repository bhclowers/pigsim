--[[
 fieldlib.lua - Utility functions for using multiple PA's to
 represent a vector field.
 D.Manura, 2017-01-14
 (c) 2011-2017 Scientific Instrument Services, Inc. (Licensed SIMION 8.1)
--]]

local FL = {}

-- Make field function from PA instances containing
-- x, y, and z components of vector field.
function FL.make_field(vxinst, vyinst, vzinst)
  if vxinst and not vyinst and not vzinst then -- scalar
    return function(x,y,z)
      local v = vxinst:potential_wc(x,y,z)
      return v
    end
  end
  local is_cyl = (vxinst.pa.symmetry_type == '2dcylindrical')
  local function field(x,y,z)
    local vx = vxinst:potential_wc(x,y,z)
    local vy = vyinst:potential_wc(x,y,z)
    local vz = vzinst and vzinst:potential_wc(x,y,z) or 0
    if vx and vy and vz then
      x,y,z = vxinst:wb_to_pa_coords(x,y,z)
      if is_cyl then
        -- rotate vector around axis.
        local theta = math.atan2(z,y)
        local cost, sint = math.cos(theta), math.sin(theta)
        vy,vz = cost*vy - sint*vz,
                sint*vy + cost*vz
      end
      return vx,vy,vz
    else
      return 0,0,0
    end
  end
  return field
end

return FL
