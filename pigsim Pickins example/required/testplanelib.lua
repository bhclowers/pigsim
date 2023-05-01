--[[
 testplanelib.lua
 Utility routines for defining and recording data on additional test
 planes.

 This can be more flexible that the "Crossing Plane" option in
 SIMION Data Recording.  Here you can define multiple test planes,
 even test planes non-orthogonal to the axes.

 By default, this generates "markers" when particles reach test
 planes.  SIMION Data Recording can be configured to record data
 when markers occur ("All Markers").  Alternately you can define
 your own "action" function to execute when particles hit the test
 plane.

 SIMION segments only get executed when particles are inside
 potential array instances.  Therefore, these test planes will not
 get triggered when particles are outside array instances.  If that
 is a problem, add to your workbench a large PA instance of low
 priority in the PA instance order (PAs tab) to cover the region in
 which test planes should be triggered.

 v2014-04-21/2008-07-24 D.Manura
 (c) 2008-2017 Scientific Instrument Services, Inc. (Licensed under SIMION 8.1)
--]]

local abs = math.abs

-- vector cross product of A x B.
local function vcross(ax, ay, az, bx, by, bz)
  return ay*bz - az*by,
         az*bx - ax*bz,
         ax*by - ay*bx
end

-- magnitude of vector A
local function vmag(ax, ay, az)
  return math.sqrt(ax^2 + ay^2 + az^2)
end

-- Creates an arbitrary coordinate system in an plane normal to given unit vector.
-- returns u and v unit vectors for axes.
local function get_coords(wx,wy,wz)
  local vx,vy,vz = vcross(wx,wy,wz, 1,0,0)
  if vmag(vx,vy,vz) < 0.1 then -- try another
    vx,vy,vz = vcross(wx,wy,wz, 0,1,0)
  end
  local ux,uy,uz = vcross(vx,vy,vz, wx,wy,wz)
  return ux,uy,uz, vx,vy,vz
end

-- Calculate time for particle at point u=<x,y,z> traveling at velocity
-- v=<vx,vy,vz> to reach test plane defined by some point p0=<p0_x,p0_y,p0_z>
-- on the test plane and some (normalized) normal vector pn=<pn_x,pn_y,pn_z>.
local function time_to_plane_helper(p0_x,p0_y,p0_z, pn_x,pn_y,pn_z,
                                    x,y,z, vx,vy,vz)
  -- dt = pn * (p0 - u) / (pn * v)
  local dt = (pn_x*(p0_x - x) + pn_y*(p0_y - y) + pn_z*(p0_z - z))
           / (pn_x*vx + pn_y*vy + pn_z*vz)
  return dt
end

-- Creates test plane definition.
-- Test plane contains point (p0_x,p0_y,p0_z) and normal vector
-- (pn_x,pn_y,pn_z).  Normal vector need not be normalized.
--
-- Returns a table containing these fields:
--
--   time_to_plane
--     This is a function (x,y,z, vx,vy,vz) --> t that returns
--     time t (microseconds) for particle at position (x,y,z) (mm)
--     traveling with velocity (vx,vy,vz) to reach the test plane,
--     assuming constant velocity.  t will be negative if particle
--     is traveling away from test plane or +-infinity if traveling
--     exactly parallel to test plane.
--
--   tstep_adjust
--     SIMION tstep_adjust segment that can be used to reduce
--     time-steps when particles approach test planes (binary boundary
--     approach).  This ensures that test plane measurements are precisely
--     on the test plane.
--
--   other_actions
--     SIMION other_actions segment that calls action() when particles
--     reach the test plane.  If action is nil, then the mark() (generate
--     marker) function is used.
--
-- tstep_adjust and other_actions must both be used since they
-- depend on each other in this implementation.
--
-- To use, markers, enable record when "All Markers" in Data Recording.
--
-- This supports both Grouped ("Particles|Grouped" enabled) and
-- non-Grouped flying.  However, markers might not work properly for
-- grouped flying since the mark() command marks *all* currently
-- flying particles, even if just one of those is currently at the
-- test plane.
--
local function create(p0_x,p0_y,p0_z, pn_x,pn_y,pn_z, action)
  action = action or mark
  local pnr = 1 / math.sqrt(pn_x*pn_x + pn_y*pn_y + pn_z*pn_z)
  pn_x,pn_y,pn_z = pn_x*pnr,pn_y*pnr,pn_z*pnr
  local function time_to_plane(x,y,z, vx,vy,vz)
    return time_to_plane_helper(p0_x,p0_y,p0_z, pn_x,pn_y,pn_z,
                                x,y,z, vx,vy,vz)
  end

  -- hit_state[n] is hit state of partice n.  Values:
  --   'approaching' - likely to hit plane soon, so do binary
  --                   boundary approach.
  --   'willhit'     - hit scheduled for next time-step
  --   'hitting'     - hitting in this time-step
  --   'hitted'      - recent hit complete, so prevent further hits
  --                   for some a short time
  --   nil           - particle has not yet approached a test plane
  local hit_state = {}

  -- number of times other_actions segment has been called.
  local oa_count = 0

  -- previous value of oa_count seen by tstep_adjust segment.
  local oa_count_last

  -- previous value of ion_number seend by tstep_adjust segment.
  local _ion_number_last

  -- Apply binary boundary approach so that time-step sizes decrease
  -- exponentially with decreasing distance to test plane.  This make
  -- time-steps stop end very close to test plane (for better
  -- accuracy).
  local function tstep_adjust()
    -- Note: special care is required to support these conditions:
    --
    --   * Multiple particles -- solved by making hit_state a table.
    --   * Grouped flying enabled or disabled -- solved by making hit_state
    --     a table.
    --   * Overlapping electrostatic and magnetic potential array instances,
    --     in which case segments are called twice (once for each array
    --     instance) -- solved by the oa_count test.

    -- When using overlapping electric and magnetic PA instances, this
    -- segment is called twice (for each particle): first for the
    -- electric instance and then for the magnetic instance.  To avoid
    -- complications, the following lines ensure this code is only
    -- called once.
    if ion_number == _ion_number_last and oa_count == oa_count_last then
      return
    end
    oa_count_last = oa_count
    _ion_number_last = ion_number

    -- estimated time to reach test plane.
    local dt  = time_to_plane(ion_px_mm,ion_py_mm,ion_pz_mm,
                              ion_vx_mm,ion_vy_mm,ion_vz_mm)

    -- smallest approach time for binary boundary approach.  Note:
    -- ion_time_step normally originally corresponds to 1 grid unit or
    -- a fraction of that, though it could be much smaller somtimes
    -- (e.g. when approaching an electrode).
    local m = ion_time_step * 1E-5

    if ion_time_step * 2 > abs(dt) then
      -- particle within hitting distance of test plane.
      -- We allow a factor of 2 in case the estimated time-step length
      -- is underpredicted.

      -- print('DEBUG:', ion_instance, ion_number, dt, hit_state[ion_number])

      local hs = hit_state[ion_number]
      if hs == 'hitted' and dt < m then  -- already hit
        -- do nothing
      elseif dt <= 0 or hs == 'willhit' then  -- hitting
        hit_state[ion_number] = 'hitting'          

        -- Set time-step to 0 when hitting to avoid confusion over
        -- whether variables represent the beginning or end of a
        -- time-step.
        ion_time_step = 0
      else  -- dt > 0 -- not yet hit
        if dt <= m then  -- (within smallest approach time)
          hit_state[ion_number] = 'willhit'  -- schedule forced hit
          ion_time_step = dt  -- attempt to hit plane exactly
        else
          hit_state[ion_number] = 'approaching'
          ion_time_step = dt * 0.5  -- half time to plane
        end
      end
    end
  end

  local function other_actions()
    oa_count = oa_count + 1

    -- Trigger action when particles reach test plane.
    local hs = hit_state[ion_number]
    if hs == 'hitting' then
      hit_state[ion_number] = 'hitted'
      action()
    end
    -- unused: local dt = time_to_plane(ion_px_mm,ion_py_mm,ion_pz_mm,
    --                                  ion_vx_mm,ion_vy_mm,ion_vz_mm)
  end

  -- Draws test plane on screen with two lines intersecting origin point.
  local function draw(color)
    color = color or 0
    local ux,uy,uz, vx,vy,vz = get_coords(pn_x, pn_y, pn_z)
    local v = 0 -- voltage
    local l = 10000 -- some large length (mm)
    simion.experimental.plot_line_segment(
      p0_x - ux*l, p0_y - uy*l, p0_z - uz*l,
      p0_x + ux*l, p0_y + uy*l, p0_z + uz*l,
      v,v, color, color, false
    )
    simion.experimental.plot_line_segment(
      p0_x - vx*l, p0_y - vy*l, p0_z - vz*l,
      p0_x + vx*l, p0_y + vy*l, p0_z + vz*l,
      v,v, color, color, false
    )
  end
  
  return {
    time_to_plane = time_to_plane,
    tstep_adjust  = tstep_adjust,
    other_actions = other_actions,
    draw = draw
  }
end

return create
