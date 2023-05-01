--[[
LUA MODULE

  contourlib81 - Field line plotting library.  (Requires SIMION 8.1)

EXAMPLE USAGE
 
  In the View screen, enter this on the SIMION command bar:
  
    dofile'contourlib81.lua'.plot(function(x,y,z) return x,y^2,0 end)
    
  Press the 'Del' button on the Particles tab to clear the plot.
    
  If using the collision_sds example with a gas flow, you can do things like
  
    dofile'contourlib81.lua'.plot(SDS.velocity)
    
  See also magnetic field plotting examples in the README.

Source: D.Manura
(c) 2011-2016 Scientific Instrument Services, Inc. (SIMION 8.1 Licensed)
--]]

local CON = {_VERSION='0.1.20161209'}

-- helper: e.g. `check_exclusive_opts(t, 'a,b;c,d')` ensures
--   `not(t.a and t.b) and not(t.b and t.c)`.
local function check_exclusive_opts(t, specs)
  for spec in specs:gmatch'[^;]+' do local found
    for name in spec:gmatch'[^,]+' do
      if t[name] then
        if found then error(found..' and '..name..' cannot both be set', 3) end
        found = name end end end
end

local IS51 = (_VERSION == 'Lua 5.1')

local myload = IS51 and function(str, source, env)
  local f = assert(loadstring('return '..str, source))
  setfenv(f, env)
  return f
end or function(str, source, env) return load(str, source, 't', env) end

-- Gets magnitude of vector
local function mag(x,y,z) return math.sqrt(x*x + y*y + z*z) end


--[[
  local CON = simion.import 'contourlib81.lua'

  CON.plot {func=func,
           npoints=npoints, npointsx=npointsx,npointsy=npointsy,npointsz=npointsz,
           xl=xl, yl=yl, zl=zl, xr=xr, yr=yr, zr=zr, x=x,y=y,z=z,
           vmin=vmin,vmax=vmax, v0=v0,v1=v1, vautoexpand=vautoexpand, vscale=vscale,
           color=color, mark=mark
  }
 
  Plots function `func` over all PA instances on the workbench View screen.

  The function `plot` must be passed a table with some of these fields
  (many of which are optional):

    `func` - a function to plot.  It must be in one of these forms:
        vx,vy,vz = f(x,y,z)
        vx,vy = f(x,y,z)
        vx = f(x,y,z)
      (x,y,z) is in workbench coordinates (mm).
      The first two forms return a vector; the last form returns a scalar.
    `npoints`
      number of points to plot in each dimension (default 20)
    `npointsx,npointsy,npointsz`
      number of points to plot in X, Y, and Z
      dimensions.  These default to `npoints` if omitted.
    `xl,yl,zl,xr,yr,zr`
      boundary box (xl,yl,zl) to (xr,yr,zr) in workbench
      coordinates (mm) to plot over.  Each value defaults to that of
      entire workbench volume.
    `x,y,z`
      `x` sets both `xl` and `xr`, and likewise for `y` and `z`.
      These may be used instead of above, for convenience, particularly
      if you want to plot over just a single cross-sectional plane.
    `vmin,vmax`
      minimum and maximum function values to plot.
      When plotting scalars, only values in this range are plotted.
      When plotting vectors, only vectors with magnitudes in this range are plotted.
      `vmin` defaults to the minimum observed value (for scalars) or 0 (for vectors).
      `vmax` defaults to the maximum observed value (for scalars) or maximum magnitude
      (for vectors).  These values may also be strings of mathematical
      expressions involving 'minimum', 'maximum', 'span'
      (i.e. maximum-minimum), 'percentile(d)' (for some integer 0..100)
      and 'median' (which is equivalent to 'percentile(50)').  Examples:
        '2*median' - twice the median value
        'percentile(10)' -- the 10% percentile
        'maximum-5' -- the maximum value minus 5.
      Note: things like vmax='percentile(98)' can be useful if certain regions have
      very high fields you want to ignore (e.g. magnetic fields near very thin
      current carrying wires).
    `v0,v1`
      function values corresponding to screen drawing of zero length and
      "reference" length (i.e. approximately the distance between points drawn)
      vectors or scalar crosses.  By default, v0 will be
      vmin (for scalars) or 0 (for vectors) and v1 will be
      vmax (for scalars and vectors).  [v0,v1] may be a superset of the
      interval [vmin,vmax] or vice-versa.
      You can also pass strings of mathematical expressions, like with `vmin/vmax`.
    `vscale` - scale factor for drawing vectors or scalar crosses.
      A scalar value or vector magnitude of v1 will be plotted on the screen
      with a size of about the distance between points in the contour plot
      times `vscale`.  `vscale` defaults to 0.9 to provide some spacing.
      (This was field was previously named `scalef`, which
      is still supported for compatibility but deprecated.)
    `vautoexpand`
      This only affects vector (not scalar) plots.
      If true, the "reference length" for drawing vectors may be automatically
      increased (for easier viewing) if it will not cause the vectors currently
      being drawn to overlap.  This is by default true.  Set it to false
      if you want more reproducible vector lengths, such as if you're
      simultaneously plotting two fields on the screen for comparison.
    `color` - (0..15) - index of color to draw with (defaults to 0 - black)
    `mark` - (true or false) - whether to draw marker at start of vector.

  The function can also accept multiple plots at the same time.  For example,

      CON.plot{npoints=41, mark=true, z=0, xl=70,xr=130,yl=70,yr=130,
        {func=bfield, color=2},
        {func=btheo, color=2},
      }
   
  will do two plots defined by the two child tables provided.  Parameters from
  the parent table are inherited in this child table, so the above is equivalent
  to (but shorter than) this:	
	
      CON.plot{
        {npoints=41, mark=true, z=0, xl=70,xr=130,yl=70,yr=130, func=bfield, color=2},
        {npoints=41, mark=true, z=0, xl=70,xr=130,yl=70,yr=130, func=btheo, color=2},
      }
	  
  The above is almost identical to doing two separate calls:
	
      CON.plot {npoints=41, mark=true, z=0, xl=70,xr=130,yl=70,yr=130, func=bfield, color=2}
      CON.plot {npoints=41, mark=true, z=0, xl=70,xr=130,yl=70,yr=130, func=btheo, color=2}
	  
  except that the two separate plots given above may have different scales
  (i.e. v0/v1/vautoexpand), whereas the combined plot ensures both are on the same scale.
--]]  
function CON.plot(...)
  local t = type((...)) == 'table' and (...) or {...}
  
  local plots = type(t[1]) == 'table' and t or {t}
  
  -- Get statistics on field.
  local vals = {}
  local vxmin,vymin,vzmin =  math.huge, math.huge, math.huge
  local vxmax,vymax,vzmax = -math.huge,-math.huge,-math.huge
  local minimum =  math.huge
  local maximum = -math.huge
  local median
  local is_vector
  local function percentile(percent)
    assert(percent >= 0 and percent <= 100, "percent not in range 0..100")
    local ir = 1 + (#vals-1)*percent/100
    local im,ip = math.floor(ir), math.ceil(ir)
    local f = ir-im
    local v = vals[im]*(1-f) + vals[ip]*f
    return v
  end
  local function summarize_statistics()
    median = percentile(50)
    print('Plotting over volume')
    if is_vector then
      print('|v|(min/max/median)=', minimum, maximum, median)
      print('vx(min/max)=', vxmin, vxmax)
      print('vy(min/max)=', vymin, vymax)
      print('vz(min/max)=', vzmin, vzmax)
    else
      print('v(min/max/median)=', minimum, maximum, median)
    end
  end

  local maxmagx,maxmagy,maxmagz = 0,0,0

  local function get_plot_params(plot)
    local t = {}
    for k,v in pairs(plots) do
      if type(k) ~= 'number' then t[k] = v end
    end -- inherit
    for k,v in pairs(plot) do t[k] = v end  -- override
  
    -- Get arguments.
    check_exclusive_opts(t, 'x,xl;x,xr;y,yl;y,yr;z,zl;z,zr')
    local func = t.func or t[1]
      if func == nil then error("parameter 'func' undefined.", 3) end
    local npoints  = t.npoints or 20
    local npointsx = t.npointsx or npoints
    local npointsy = t.npointsy or npoints
    local npointsz = t.npointsz or npoints
    local vscale = t.vscale or t.scalef or 0.9
    local vmin = t.vmin or 'minimum'
    local vmax = t.vmax or 'maximum'
    local v0 = t.v0 or nil
    local v1 = t.v1 or nil
    local vautoexpand = t.vautoexpand; if vautoexpand == nil then vautoexpand = true end
    local color = t.color or 0
    local mark = t.mark or false
    local bounds = simion.wb.bounds
    local xl = t.xl or t.x or bounds.xl
    local yl = t.yl or t.y or bounds.yl
    local zl = t.zl or t.z or bounds.zl
    local xr = t.xr or t.x or bounds.xr
    local yr = t.yr or t.y or bounds.yr
    local zr = t.zr or t.z or bounds.zr
    if xl > xr then xl, xr = xr, xl end -- swap if out of order
    if yl > yr then yl, yr = yr, yl end
    if zl > zr then zl, zr = zr, zl end


    -- Determine spacing between points to plot.
    local skip_wx = (xr - xl) / npointsx
    local skip_wy = (yr - yl) / npointsy
    local skip_wz = (zr - zl) / npointsz
  
    -- Evaluate points, caching in case field computation is expensive (e.g. Biot-Savart).
    local xs,ys,zs = {},{},{}
    local vxs,vys,vzs = {},{},{}
    local is_vector_current = false
    for wz=zl+skip_wz*0.5,zr,skip_wz == 0 and 1 or skip_wz do
    for wy=yl+skip_wy*0.5,yr,skip_wy == 0 and 1 or skip_wy do
    for wx=xl+skip_wx*0.5,xr,skip_wx == 0 and 1 or skip_wx do
      local vx,vy,vz = func(wx, wy, wz)
      if vy then is_vector_current = true end
      vx,vy,vz = vx or 0, vy or 0, vz or 0
      local n = #xs+1
      xs[n],ys[n],zs[n], vxs[n],vys[n],vzs[n] = wx,wy,wz, vx,vy,vz
    end end end
  
    if #xs == 0 then error('no points plotted', 2) end
  
    -- Update statistics.
    assert(is_vector == nil or is_vector == is_vector_current,
           'scalar and vectors not supported in same plot')
    is_vector = is_vector_current
    for i=1,#xs do
      local vx,vy,vz = vxs[i],vys[i],vzs[i]
      local v = is_vector and mag(vx,vy,vz) or vx
      minimum = math.min(minimum, v)
      maximum = math.max(maximum, v)
      vxmin = math.min(vxmin, vx); vymin = math.min(vymin, vy); vzmin = math.min(vzmin, vz)
      vxmax = math.max(vxmax, vx); vymax = math.max(vymax, vy); vzmax = math.max(vzmax, vz)
      vals[#vals+1] = v
    end
    table.sort(vals)

    function t.normalize()
      -- Expand vmin/vmax/v0/v1 value limits.
      local function expandstr(s)
        local env = {minimum=minimum, maximum=maximum, median=median, span=maximum-minimum,
                     percentile=percentile}
        local f = myload(s, nil, env)
        return f()
      end
      if type(vmin) == 'string' then vmin = expandstr(vmin) end
      if type(vmax) == 'string' then vmax = expandstr(vmax) end
      if type(v0)   == 'string' then v0   = expandstr(v0) end
      if type(v1)   == 'string' then v1   = expandstr(v1) end
      v0 = v0 or (is_vector and 0 or vmin)
      v1 = v1 or vmax
  
      -- Filter values, normalize, and get vector x,y,z component ranges.
      for i=1,#xs do
        local vx,vy,vz = vxs[i],vys[i],vzs[i]
        local v = is_vector and mag(vx,vy,vz) or vx
        if v >= vmin and v <= vmax then
          if is_vector then
            if v ~= 0 then
              local rescale = (v0==v1) and 1 or (v - v0)/(v1-v0)
              vx = (vx/v)*rescale
              vy = (vy/v)*rescale
              vz = (vz/v)*rescale
            end -- else vx,vy,vz=0,0,0 already
          else -- scalar
            vx = (v0==v1) and 1 or (vx - v0)/(v1-v0)
            if vx < 0 then vx = 0 end -- won't plot scalars under v0
          end
          vxs[i],vys[i],vzs[i] = vx,vy,vz

          if true then maxmagx = math.max(maxmagx, math.abs(vx)) end
          if vy   then maxmagy = math.max(maxmagy, math.abs(vy)) end
          if vz   then maxmagz = math.max(maxmagz, math.abs(vz)) end
        else
          vxs[i],vys[i],vzs[i] = nil,nil,nil -- clear
        end
      end
      if not vautoexpand or not is_vector then
        maxmagx,maxmagy,maxmagz = 1,1,1
      end
    end

    function t.plot_single()
      assert(simion.pas, 'This function requires SIMION 8.1.')

      local plot_segment = assert(simion.experimental.plot_line_segment)
      
  
      -- Get scaling factor to screen.
      local cell_wx = (skip_wx == 0) and math.max(skip_wy, skip_wz) or skip_wx 
      local cell_wy = (skip_wy == 0) and math.max(skip_wx, skip_wz) or skip_wy
      local cell_wz = (skip_wz == 0) and math.max(skip_wx, skip_wy) or skip_wz
      local scalex = (cell_wx == 0) and 1 or cell_wx / maxmagx
      local scaley = (cell_wy == 0) and 1 or cell_wy / maxmagy
      local scalez = (cell_wz == 0) and 1 or cell_wz / maxmagz
      local scale = vscale * math.min(scalex, scaley, scalez)
      if not is_vector then scale=scale*0.5 end
  
      -- Rescale to screen.
      for i=1,#xs do
        if vxs[i] then
          vxs[i] = vxs[i]*scale
          vys[i] = vys[i]*scale
          vzs[i] = vzs[i]*scale
        end
      end

      if math.abs(scale) == math.huge then
        print 'warning: zero values over entire plot'
        return
      end

      -- Plot field.
      for i=1,#xs do
        local wx,wy,wz, vx,vy,vz = xs[i],ys[i],zs[i], vxs[i],vys[i],vzs[i]
        if vx then
          if mag(vx,vy,vz) > 0 then -- non-zero
            if is_vector then
              plot_segment(wx+vx, wy+vy, wz+vz, wx, wy, wz, nil,nil,color,nil, mark)
            else -- scalar
              -- Draw line in all three directions (for easy viewing in all orientations)
              plot_segment(wx,wy-vx,wz, wx,wy+vx,wz, nil,nil,color,nil, false)
              plot_segment(wx,wy,wz-vx, wx,wy,wz+vx, nil,nil,color,nil, false)
              plot_segment(wx,wy,wz, wx-vx,wy,wz, nil,nil,color,nil, false)
              plot_segment(wx,wy,wz, wx+vx,wy,wz, nil,nil,color,nil, mark)  -- mark last (top)
            end
          end
        end
      end -- for

    end -- function
   
    return t
  end -- function

  
  for i=1,#plots do
    plots[i] = get_plot_params(plots[i])
  end
  summarize_statistics()
  
  for i=1,#plots do
    plots[i].normalize()    
  end

  for _, plot in ipairs(plots) do
    plot.plot_single()
  end
end

local function vmag(ax,ay,az)
  return math.sqrt(ax^2 + ay^2 + az^2)
end
local function vnorm(ax,ay,az)
  local m = vmag(ax,ay,az)
  if m == 0 then
    return 0,0,0
  else
    return ax/m, ay/m, az/m
  end
end
-- Returns vector cross product of A x B.
local function vcross(ax,ay,az, bx,by,bz)
  return ay*bz - az*by,
         az*bx - ax*bz,
         ax*by - ay*bx
end

-- Gets minimum length of all grid cells in all PA's the workbench.
function CON.min_grid_size()
  local d = math.huge
  for i=1, #simion.wb.instances do
    local inst = simion.wb.instances[i]
    local di = inst.scale * math.min(inst.pa.dx_mm, inst.pa.dy_mm, inst.pa.dz_mm)
    d = math.min(d, di)
  end
  return d
end

--[[
 Plots streamlines for given field, starting at seed point (x,y,z) mm
 with given step distance d in mm (defaults to 1/10 grid cell size if nil),
 in both directions, and color number (defaults to 2 green if nil).
 Field field(x,y,z) is a function given point (x,y,z) mm and returning
 vector (Ex,Ey,Ez) in any unit.
 field defaults to workbech E field if nil.
 Will terminates if cycle detected in path to avoid infinite or long loop.
--]]
function CON.streamline(field, x,y,z, d, color)
  field = field or function(x,y,z) return simion.wb:efield(x,y,z) end
  d = d or CON.min_grid_size() / 10
  color = color or 2
  local plot_segment = assert(simion.experimental.plot_line_segment)
  local x0,y0,z0 = x,y,z
  for dir = -1,1, 2 do
    local x1,y1,z1 = x0,y0,z0
    local xt,yt,zt = x0,y0,z0
    local ds = d * dir
    local len = 0
    for i=1,math.huge do
      -- Avoid cycles
      local dist = vmag(x1-x0, y1-y0, z1-z0)
      if len > dist * 5 then break end -- large
      if i % 10 == 0 then xt,yt,zt = x1,y1,z1 end
      if i % 10 == 9 and vmag(x1-xt,y1-yt,z1-zt) < 2*d then break end -- small

      local ex,ey,ez = field(x1,y1,z1)
      if ex and vmag(ex,ey,ez) ~= 0 then
        ex,ey,ez = vnorm(ex,ey,ez)
        local x2,y2,z2 = x1 + ex*ds, y1 + ey*ds, z1 + ez*ds
        plot_segment(x1,y1,z1, x2,y2,z2, nil,nil,color,nil, false)
        --print('DEBUG', i,x1,y1,z1, len, dist) simion.sleep(0.001)
        x1,y1,z1 = x2,y2,z2
        len = len + d
      else
        break
      end
    end
  end
end

--[[
 Plots n streamlines for given field (defaulting to workbech E field if nil),
 starting at seed points uniformly spaced
 between (x1,y1,z1) and (x2,y2,z2) with given step.
 distance d in mm (defaults to 1/10 grid cell size if nil)
 and color number (defaults to 2 green if nil).
--]]
function CON.streamlines(field, x1,y1,z1, x2,y2,z2, n, d, color)
  n = n or 10
  for i=1,n do
    local f = n == 0 and 0 or (i-1)/(n-1)
    local x,y,z = x1 + (x2-x1)*f, y1 + (y2-y1)*f, z1 + (z2-z1)*f
    CON.streamline(field, x,y,z, d, color)
  end
end

--[[
 Plots n random streamlines given field
 (defaulting to workbech E field if nil),
 using given color.
--]]
function CON.streamlines_random(field, n, d, color)
  for i=1,n do
    local bounds = simion.wb.bounds
    local x = bounds.xl + (bounds.xr - bounds.xl) * rand()
    local y = bounds.yl + (bounds.yr - bounds.yl) * rand()
    local z = bounds.zl + (bounds.zr - bounds.zl) * rand()
    CON.streamline(field, x,y,z, d, color)
  end
end

--[[
 Walk along streamline for field (defaulting to workbech E field if nil)
 starting at point (x0,y0,0) mm, within plane perpendicular
 to vector (ux,uz,uz), with step distance d mm.
 At each point (x,y,z) mm, calls function f(x,y,z).
--]]
function CON.walk(field, x0,y0,z0, ux,uy,uz, d, f)
  field = field or function(x,y,z) return simion.wb:efield(x,y,z) end
  d = d or CON.min_grid_size() / 10
  for dir=1,-1,-2 do
    local x,y,z = x0,y0,z0
    local xt,yt,zt = x0,y0,z0
    local len = 0
    for i=1,math.huge do
      -- Avoid cycles
      local dist = vmag(x-x0, y-y0, z-z0)
      if len > dist * 5 then break end -- large
      if i % 10 == 0 then xt,yt,zt = x,y,z end
      if i % 10 == 9 and vmag(x-xt,y-yt,z-zt) < 2*d then break end -- small

      local ex,ey,ez = f(x,y,z)
      if ex and vmag(ex,ey,ez) ~= 0 then
        f(x,y,z)
        local sx,sy,sz = vnorm(vcross(ex,ey,ez, ux,uy,uz))
        x = x + sx * d * dir
        y = y + sy * d * dir
        z = z + sz * d * dir
        len = len + d
      else
        break
      end
    end
  end
end

--[[
 Plots streamlines for given field (defaulting to workbech E field if nil)
 starting at seedpoint along given walk (CON.walk), spaced by
 about space mm.
--]]
function CON.streamlines_walk(field, x0,y0,z0, ux,uy,uz, space, d, color)
  d = d or CON.min_grid_size() / 10
  space = space or CON.min_grid_size() * 10
  local skip = math.max(math.floor(space / d + 0.5), 1)
  local i = 0
  CON.walk(field, x0,y0,z0,ux,uy,uz, d, function(x,y,z)
    if i % skip == 0 then
      CON.streamline(x,y,z, d, color)
    end
    i = i + 1
  end)
end

-- Also dofile'contourlib.lua'(...) as shorthand for dofile'contourlib.lua'.plot(...)
setmetatable(CON, {__call=function(_, ...) return CON.plot(...) end})

-- Whether drawing is supported in this version of SIMION.
CON.can_draw = (simion.experimental.plot_line_segment ~= nil)


return CON
