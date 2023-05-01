--[[
 comsol_to_pa.lua
 Converts COMSOL(R) calculated gas flows to SIMION PA files.
 
 This will create one SIMION PA for each flow parameter.
 
 Some adjustment to this program may be necessary for your own purposes.
 
 WARNING: What are units for distance, velocity, and pressure?
 Assume SI?
 See SIMION SDS and HS1 for flow units
 they use.  This program makes some attempt to convert these, but
 you may need to rescale units.
 
 WARNING: For best precision, take care about whether parameters
 represent the corners or centers of grid cells.  Flows may be
 off by one grid unit if you don't accounts for this when
 creating PA's and interpreting them in SIMION.  So, please check
 if maximum precision is desired.

 A typical way to plot velocity vectors (using "contourlib81" in the "contour"
 example), assuming the vx,vy,vz arrays are PA instances #1-3 on the PA
 Instances list on the View screen PAs tab, is as follows
 (in SIMION >= 8.1.1.0):
 
   local vxinst = simion.wb.instances[1]
   local vyinst = simion.wb.instances[2]
   local vzinst = simion.wb.instances[3]
   local function velocity(x,y,z)
     return vxinst:potential_wc(x,y,z),
            vyinst:potential_wc(x,y,z),
            vzinst:potential_wc(x,y,z)
   end
   simion.import'../contour/contourlib81.lua'.plot{func=velocity, mark=true}
 
 Any NaN (not-a-number) values in the COMSOL file are converted to 0.
 
 D.Manura, v2015-05-25 (based on fluent_ip_to_simion_pas.lua)
 (c) 2012-2016 Scientific Instrument Services, Inc. (Licensed SIMION 8.1)
--]]

--[[
  Convert COMSOL file (filename) to
  SIMION PA (pa_filename).
  Optionally rescale values by factor `rescale_factor` (defaults to 1).
--]]
local function convert_file(filename, pa_filename, rescale_factor)
  rescale_factor = rescale_factor or 1

  -- Open input file.
  local fh = assert(io.open(filename, 'r'))

  -- Attempts to read line of numbers.
  local function read_line_of_numbers(fh)
    local nums = {}
    local line = assert(fh:read'*l')
    for word in line:gmatch'%S+' do
      if word == 'NaN' then word = 0 end  -- convert NaN values to 0.
      if not tonumber(word) then return nil, line end
      local val = tonumber(word)
      table.insert(nums, val)
    end
    return nums
  end

  -- Attempts to read line matching string pattern `pattern`.
  local function expect_line(fh, pattern, line)
    local line = line or assert(fh:read'*l')
    if not line:match(pattern) then
      error 'format not recognized'
    end
  end

  -- Read header of coordinates for each dimension of
  -- rectilinear grid.
  -- There can be 1 to 3 coordinates.
  local line
  expect_line(fh, '%% Grid')
  local xs = read_line_of_numbers(fh)
  local ys, zs
  ys, line = read_line_of_numbers(fh)
  if ys then
    zs, line = read_line_of_numbers(fh)
    if zs then  -- 3D array
      line = assert(fh:read'*l')
    else -- 2D array
      zs = {0}
    end
  else -- 1D array
    ys = {0}
    error('1D array not allowed')
  end
  expect_line(fh, '%% Data', line)

  -- Calculates statitics (average, min, and max) for
  -- step sizes in numbers in array `vals`.
  local function calculate_delta(vals)
    if #vals < 2 then return end

    local dvals = {}
    for i=1,#vals-1 do
      dvals[i] = vals[i+1] - vals[i]
    end
    local dmin
    local dmax
    local dsum = 0
    for i=1,#dvals do
      local d = dvals[i]
      dmin = math.min(dmin or d, d)
      dmax = math.min(dmax or d, d)
      dsum = dsum + d
    end
    local dave = dsum / #dvals

    local derror = math.abs((dmax - dmin) / dave)
    print('delta percent span=', derror)
    if derror > 0.01 then
      error('grid points not equally spaced')
    end
  
    return dave, dmin, dmax
  end

  -- Compute statistics on array dimensions.
  local dx_ave = calculate_delta(xs)
  local dy_ave = calculate_delta(ys)
  local dz_ave = calculate_delta(zs)

  -- Conversion factors.
  local MM_PER_M = 1000 -- assuming positions in original file in units of m.

  -- Create PA object.
  local pa = simion.pas:open()
  pa:size(#xs,#ys,#zs)
  pa.filename = pa_filename
  pa.refinable = false  -- prevent SIMION prompting to refine this PA.
  pa.dx_mm = MM_PER_M * dx_ave
  pa.dy_mm = MM_PER_M * dy_ave
  if #zs > 1 then pa.dz_mm = MM_PER_M * dz_ave end
  --pa.mirror_x = true  -- optionally enable mirroring like this

  -- Copy values from input array to SIMION PA.
  for zg=0,#zs-1 do
  for yg=0,#ys-1 do
    local vals = read_line_of_numbers(fh)
    if #vals ~= #xs then
      error('unexpected number of values found on line')
    end
  
    for xg=0,#vals-1 do
      local val = vals[xg+1]
      if val ~= val then -- not-a-number (NaN)
        val = 0  -- treat as zero since NaN's should not be stored in PA's.
      end
      pa:potential(xg,yg,zg, val * rescale_factor)
    end
  end
  end

  -- Any extra data in input file?
  local more = fh:read'*a'
  if not (more or ''):match'^%s*$' then
    error('unexpected extra data found in file: ['..more:sub(1,100)..'...]')
  end

  -- Close input file.
  fh:close()

  -- Save PA.
  pa:save()
end

convert_file('FullVelocity_x_3.txt', 'comsol_vx3.pa', 1)
convert_file('FullVelocity_y_3.txt', 'comsol_vy3.pa', 1)
convert_file('FullVelocity_z_3.txt', 'comsol_vz3.pa', 1)
convert_file('FullPressure_3.txt', 'pressure3.pa', 1)
--convert_file('FullPressure_2_Torr.txt', 'pressure2_Torr.pa', 1)

print 'DONE'
