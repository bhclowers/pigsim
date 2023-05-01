--[[
  SIMION Lua workbench user program for
  Statistical Diffusion Simulation (SDS) Model.
  
  Runs under SIMION version >= 8.1.1.16 as is.
  
  == DESCRIPTION ==
  
  This is a refined model of using gas kinetic numbers to simulate
  ion mobility and diffusion via Stokes' Law, with diffusion derived
  from kinetic jump statistics data, efficiently simulating pressures
  in the atmospheric range.  This version supports multiple
  ion definitions as well as temperature, pressure, and bulk gas
  velocity fields (incorporated as array files or external functions).
  
  See the accompanying README.html, collision_sds_documentation.pdf
  files for full details on this program (as well as the papers
  cited).
  
  This program requires the existence of the following file:
  
    mbmr.dat --  collision statistics file
                 (see load_diffusion_statistics below)
  
  Each file below is optional and ignored if omitted:
  
    m_defs.dat -- ion parameter definition file
                  holds ion diameter(nm), Ko(10-4m2V-1s-1) by ion mass
                  (see load_mass_data below)
  
    field definition files...Note: each file must have the same
    dimensions as the associated potential array (see load_array below):
  
      p_defs.dat    holds pressure field data (torr)
      t_defs.dat    holds Temperature field data (K)
      vx_defs.dat   holds Vx bulk gas velocity (m/s)
      vy_defs.dat   holds Vy bulk gas velocity (m/s)
      vz_defs.dat   holds Vz bulk gas velocity (m/s)
  
  WARNING: If magnetic and electrostatic PAs overlap, then this
  program must be active only in the magnetic one.  You can disable
  SDS on a certain instance (e.g. #2) by doing this:
  
    local SDS = simion.import 'collision_sds.lua'
    SDS.instances[2] = false
  
  Trajectory quality of 0 or negative is recommended for speed.
  
  == SOURCE ==
  
  This code and corresponding documentation are based on the original
  SDS model by Dave Dahl, 2004-09-27, included in SIMION PRG format in
  the supplementary material in the paper:
  
    Anthony D. Appelhans and David A. Dahl.
    SIMION ion optics simulations at atmospheric pressure.
    International Journal of Mass Spectrometry.
    Volume 244, Issue 1, 15 June 2005, Pages 1-14.
    http://dx.doi.org/10.1016/j.ijms.2005.03.010
  
  The SIMION PRG was converted to SIMION Lua by D.Manura, 2007-05,
  and incorporated into SIMION 8 by permission of Appelhans and Dahl.
  See the README.html for changes made.

  The current form is
    (c) 2007-2014 Scientific Instrument Services, Inc. (Licensed SIMION 8.1)
--]]
simion.workbench_program()

local M = {}
M._VERSION = '8.1.2.7.20140407'

-- Options passed as second argument to simion.import().
local opt = ... or ''
local opt_noinstall = opt:match'noinstall'
local opt_version   = opt:match'v([%d%.]+)'

-- Version check.  (Added in 8.1.1.16.)
local function compare_version(ver1, ver2)
  local function normalize(ver)
    return (ver:gsub('%d+', function(s) return ('%010d'):format(s) end))
  end
  local s1, s2 = normalize(ver1), normalize(ver2)
  return s1 == s2 and 0 or s1 >= s2 and 1 or -1
end
if opt_version then
  if compare_version(opt_version, M._VERSION) > 0 then
    error(('Please upgrade collision_sds.lua.\n'..
           'Current Version=%s\nRequired Version=%s')
          :format(M._VERSION, opt_version), 2)
  end
end


M.segment = {}

-- detect program errors ( http://simion.com/issue/490 )
if checkglobals then checkglobals() end

-- check SIMION version.
assert(simion.VERSION >= simion.VERSION '8.1.1.16', "SIMION 8.1.1.16 or above required")


-- Gets path of file relative to current working directory or current file.
local function find_file(filename)
  local function file_exists(filename)
    local fh = io.open(filename)
    if fh then fh:close(); return true else return false end
  end
  if file_exists(filename) then
    return filename
  else
    -- Relative directory path to this file.
    local PATH = debug.getinfo(2).source:gsub('^@(.-)[^/\\]*$', '%1'):gsub('^$', './')
    if file_exists(PATH .. filename) then
      return PATH .. filename
    end
  end
end

local TF = simion.import(find_file('textfilelib.lua'))

--##
--## SECTION: Common Constants
--##


-- Copy of function references (can increase performance and conciseness).
local abs   = math.abs
local min   = math.min
local max   = math.max
local sqrt  = math.sqrt
local log10 = math.log10
local exp   = math.exp
local modf  = math.modf
local floor = math.floor
local rand  = rand
local print = print
local assert= assert
local type  = type
local error = error
local pairs = pairs
local ipairs= ipairs


-- Physical constants used.
local ELEMENTARY_CHARGE = 1.602176462e-19  -- elementary charge (C/e)
local K_BOLTZMANN  = 1.3806503e-23         -- Boltzmann constant (JK-1)
local N_AVOGADRO   = 6.02214199e23         -- Avogadro's number
local MOL_VOLUME   = 22.413996e-3          -- Volume (m3) of one mol
                                           -- of ideal gas at 0 C, 101.325 kPa
local AMU_TO_KG    = 1.66053873e-27        -- mass of one u in kg
local PI           = math.pi               -- value of PI
local MM_HG_TO_PA  = 133.322               -- conversion of mmHg to Pa
local STP_TEMP     = 273.15                -- standard temperature (K)
local M_AIR        = 28.94515              -- Effective mass of air (u)
local D_AIR        = 0.366                 -- Effective diameter of air (nm)

--##
--## SECTION: Util
--##


-- Tests whether string `s` ends with string `s2` (case insensitive).
local function endswithi(s, s2)
  s  =  s:lower()
  s2 = s2:lower()
  return #s >= #s2 and s:sub(-#s2) == s2
end


--##
--## SECTION: I/O Utilities - general functions for input/output
--##


local function load_array(filename)
  if TF.file_exists(filename) then
    return simion.import(find_file('arraylib.lua')).load_array(filename)
  end
end


-- Formatted print.
local function printf(...) print(string.format(...)) end


--##
--## SECTION: PA Instance
--##


-- rotate vector CCW around X (looking down positive X axis)
local function x_rotate(theta, x, y, z)
  local sint, cost = math.sin(theta), math.cos(theta)
  local yb = cost * y - sint * z
  local zb = sint * y + cost * z
  return x, yb, zb
end

local function objtype(o)
  return (type(o) == 'userdata' and o.az) and 'painstance' or
         (type(o) == 'userdata' and o.nx) and 'pa' or
         (type(o) == 'table' and getmetatable(o) == nil) and 'rawtable' or
         (type(o) == 'table'    and o.nx) and 'luaarray' or
         type(o)
end


--[[
 Rotates vector (vx,vy,vz) from PA array coordinates to PA volume coordinates,
 This is the inverse rotation of a point (px,py,pz) in PA volume coordinates
 to PA array coordinates (via instance:pa_to_array_coords).
 We assume antisymmetry for vector components perpendicular to mirror planes,
 e.g. vx(x,y,z) = -vx(-x,+-y,+-z) in 3D.
--]]
local array_to_pa_orient = {} -- keyed by PA symmetry_type.
array_to_pa_orient['2dcylindrical'] = function(px,py,pz, vx,vy,vz)
  vx = vx * (px >= 0 and 1 or -1)
  vx,vy,vz = x_rotate(math.atan2(pz,py), vx,vy,vz)
  return vx,vy,vz
end
array_to_pa_orient['2dplanar'] = function(px,py,pz, vx,vy,vz)
  vx = vx * (px >= 0 and 1 or -1)
  vy = vy * (py >= 0 and 1 or -1)
  -- note: no antisymmetry in vz.
  return vx,vy,vz
end
array_to_pa_orient['3dplanar'] = function(px,py,pz, vx,vy,vz)
  vx = vx * (px >= 0 and 1 or -1)
  vy = vy * (py >= 0 and 1 or -1)
  vz = vz * (pz >= 0 and 1 or -1)
  return vx,vy,vz
end


--[[
  Returns a function `s = f(x,y,z)` that maps the workbench coordinates
  `(x,y,z)` in mm to the value `s` of some scalar field.  The array-like
  object `o` is used to build `f`.  `o` may be `f` itself, a PA instance,
  a PA, or a Lua array returned by `load_array`.
--]]
local function wrap_scalar_field(o)
  if objtype(o) == 'function' then
    local f = o
    return f
  elseif objtype(o) == 'painstance' then
    local painst = o
    local pa = painst.pa
    local function f(x,y,z)
      local xp,yp,zp = painst:wb_to_pa_coords(x,y,z)
      return pa:inside_vc(xp,yp,zp) and pa:potential_vc(xp,yp,zp) or 0
    end
    return f
  elseif objtype(o) == 'pa' then  -- [*1]
    local pa = o
    local function f(x,y,z)
      local xp,yp,zp = wb_coords_to_pa_coords(x,y,z)
      return pa:inside_vc(xp,yp,zp) and pa:potential_vc(xp,yp,zp) or 0
    end
    return f
  elseif objtype(o) == 'luaarray' then  --[*1]
    local function f(x,y,z)
      local px,py,pz = wb_coords_to_pa_coords(x,y,z)
      local xa,ya,za = pa_coords_to_array_coords(px,py,pz)
      return o(xa,ya,za)
    end
    return f
  else
    error('unrecognized array type: ' .. tostring(o))
  end
end
-- [*1] Assumes points correspond in `simion.wb.instances[ion_instance].pa`
--      and `o`.


--[[
  Returns a function `vx,vy,vz = f(x,y,z)` that maps the workbench coordinates
  `(x,y,z)` in mm to the value `(vx,vy,vz)` of some vector field.  The array-like
  object `o` is used to build `f`.  `o` may be `nil`, `f` itself, or a table
  of three objects (three PA instances or three Lua arrays returned by `load_array`).
  Some of those three object may be `nil` and default to `defaults.vx|vy|vz` respectively,
  which are functions that return a value.
--]]
local function wrap_vector_field(o, defaults)
  if type(o) == 'function' then
    local f = o
    return f
  elseif type(o) == 'table' and 
        (objtype(o[1])=='painstance' or objtype(o[2])=='painstance' or objtype(o[3])=='painstance')
  then
    local xpa, ypa, zpa = unpack(o, 1, 3)
    local function ok(o) return objtype(o) == 'painstance' or o == nil end
    assert(ok(xpa) and ok(ypa) and ok(zpa), 'type mismatch')
    local anyinstance = xpa or ypa or zpa
    local anypa = anyinstance.pa
    local vxpa = xpa and xpa.pa or nil
    local vypa = ypa and ypa.pa or nil
    local vzpa = zpa and zpa.pa or nil
    local fvx = vxpa and vxpa.potential_vc or defaults.vx
    local fvy = vypa and vypa.potential_vc or defaults.vy
    local fvz = vzpa and vzpa.potential_vc or defaults.vz
    local array_to_pa_orient_ = array_to_pa_orient[anypa.symmetry_type]
    local function f(x,y,z)
      local px,py,pz = anyinstance:wb_to_pa_coords(x,y,z)
      if anypa:inside_vc(px,py,pz) then
        local vx = fvx(vxpa, px,py,pz)
        local vy = fvy(vypa, px,py,pz)
        local vz = fvz(vzpa, px,py,pz)
        vx,vy,vz = array_to_pa_orient_(px,py,pz, vx,vy,vz)
        vx,vy,vz = anyinstance:pa_to_wb_orient(vx,vy,vz)  --rotate
        return vx,vy,vz
      else
        return 0, 0, 0
      end
    end
    return f
  elseif type(o) == 'table' and
         (objtype(o[1])=='luaarray' or objtype(o[2])=='luaarray' or objtype(o[3])=='luaarray')
  then
    local xobj, yobj, zobj = unpack(o)
    local function ok(o) return objtype(o) == 'luaarray' or o == nil end
    assert(ok(xobj) and ok(yobj) and ok(zobj), 'type mismatch')
    local anyobj = xobj or yobj or zobj
    local fvx = xobj or defaults.vx
    local fvy = yobj or defaults.vy
    local fvz = zobj or defaults.vz
    local symmetry_type = anyobj.symmetry:gsub('%b[]', '')
    local array_to_pa_orient_ = array_to_pa_orient[symmetry_type]
    local function f(x,y,z)
      local px,py,pz = wb_coords_to_pa_coords(x,y,z)
      local xa,ya,za = pa_coords_to_array_coords(px,py,pz)
      local vx = fvx(xa,ya,za)
      local vy = fvy(xa,ya,za)
      local vz = fvz(xa,ya,za)
      vx,vy,vz = array_to_pa_orient_(px,py,pz, vx,vy,vz)
      vx,vy,vz = pa_orient_to_wb_orient(vx,vy,vz)
      return vx,vy,vz
    end
    return f
  end
end


--[[
 Attempts to locate flow data in workbench PA instance list.
 PA instance file name must end with '_name.pa' for some
 value of 'name' (ignoring case, but preferrably lowercase).
 Prior to SIMION 8.1, always returns `nil`.
--]]
local function find_instance(name)
  if not simion.wb then return nil end  -- skip if not SIMION 8.1
  local painst
  for i=1,#simion.wb.instances do
    local _painst = simion.wb.instances[i]
    local filename = _painst.filename
    if endswithi(filename, '_' .. name .. '.pa') then
      if i == #simion.wb.instances then   
        error('Gas flow array ' .. filename .. ' has the highest priority ' ..
              'number ('..i..') in the PA instances list.  Please use L-/L+ ' ..
              'buttons to give gas flow arrays lower priority than ' ..
              'electric and magnetic arrays.\n')
      end
      painst = _painst; break
    end
  end
  return painst
end
M.find_instance = find_instance



--##
--## SECTION: Diffusion Statistics
--##


-- Number of collisions represented in distribution data (see
-- load_diffusion_statistics).
local N_DIST_COLLISIONS = 100000

-- Number of ICDFs in represented in distribution data (see
-- load_diffusion_statistics)
local N_DIST = 5

-- Number of data points per ICDF.
-- (Note: the extra point, 1002, might be unnecesary.)
local N_DIST_POINTS = 1002


--[[
  Loads and returns diffusion statistics data from file.
  
  This data represents a number of inverse cumulative probability
  density functions (ICDFs), one for each of a series of mass ratios
  (mass_ion/mass_gas), which are powers of 10 from 1 to 10,000
  (referred to as functions #1 - #5).
  
  Each function represents the distance that a particle of the given
  mass ratio travels (between starting and ending points) upon doing
  N_DIST_COLLISIONS collisions, assuming a mean-free-path of 1 unit.
  The function is represented by N_DIST_POINTS data points to cover
  percentiles 0% to 100% in increments of 0.1%.
  
  The data are fully randomized Maxwell-Boltzmann (MB) statistics:
    MB randomized initial energy (normalized to average velocity of 1).
    Randomized initial direction.
    Poisson randomized collision distance (normalize to average of 1).
    MB randomized collision gas energies (normalized to have same mean
      energy of the ions)
    Randomized equally probable impact points.
  
   N_DIST_COLLISIONS is also stored in the returned table.
  
   Raises on error.
--]]
local function load_diffusion_statistics(filename)
  local raw = TF.read_file_numbers(filename)  -- raises
  assert(#raw == N_DIST_POINTS * N_DIST)

  -- Split ICDFs into separate arrays.
  local stats = {}
  for i=1,N_DIST do
    local icdf = {}
    for j=1,N_DIST_POINTS do
      icdf[j] = raw[(i-1) * N_DIST_POINTS + (j-1) + 1]
    end
    stats[i] = icdf
  end
  stats.N_DIST_COLLISIONS = N_DIST_COLLISIONS

  return stats
end


--[[
  Computes a random walk jump distance for diffusion effect.  This is
  the distance traveled (between starting and ending points) by a
  particle after having done N_DIST_COLLISIONS collisions, assuming
  log10(mass_ion/mass_gas) = log_mr_ratio and normalized conditions
  described in load_diffusion_statistics.
  
  This value is interpolated from statistics in stats, which
  represents inverse cumulative density function (ICDF).  It involves
  a log-log scale interpolation of ICDFs at the nearest mass ratios
  (the value of each ICDFs is itself interpolated as well).
--]]
local function diff_dist_steps(stats, log_mr_ratio)
  -- Select uniformly random percentile * 10 (0, 1000) for the random
  -- jump distance.
  local n = rand() * (N_DIST_POINTS - 2)

  -- Select the two ICDFs to interpolate between.  We'll interpolate
  -- between ICDFs icdf1 and icdf2.  These are selected based on mass
  -- ratio.
  local iicdf = (log_mr_ratio <= 1) and 1 or
                (log_mr_ratio <= 2) and 2 or
                (log_mr_ratio <= 3) and 3 or 4
  local icdf1,icdf2 = stats[iicdf], stats[iicdf+1]

  -- Select the indices in each ICDF to interpolate between.
  local ilow = floor(n) + 1; local ihigh = ilow + 1

  -- Interpolate inside both ICDFs.
  local weight = n % 1
  local dist1_steps = (icdf1[ihigh] - icdf1[ilow]) * weight + icdf1[ilow]
  local dist2_steps = (icdf2[ihigh] - icdf2[ilow]) * weight + icdf2[ilow]

  -- Interpolate between both ICDFs (log-log interpolation).
  dist1_steps, dist2_steps = log10(dist1_steps), log10(dist2_steps)
  local weight = log_mr_ratio - (iicdf - 1)
  local dist_steps = 10^((dist2_steps - dist1_steps) * weight + dist1_steps)

  return dist_steps
end


--##
--## SECTION: Mass Data (massdata)
--##


--[[
  Gets scaling constant air_to_gas such that Kogas = Koair *
  air_to_gas.  Given diameter d_ion (nm) and mass mass_ion (u)
  of ion and given diameter d_gas (nm) and mass mass_gas (u) of gas.
--]]
local function get_air_to_gas(d_ion,mass_ion, d_gas,mass_gas)
  -- Reduced mass with respect to air and gas
  local reduced_mass_air = mass_ion * M_AIR    / (mass_ion + M_AIR)
  local reduced_mass_gas = mass_ion * mass_gas / (mass_ion + mass_gas)

  -- Scaling constant to convert from Koair to Kogas.
  local air_to_gas =
    ((d_ion + D_AIR) / (d_ion + d_gas) )^2 *
    sqrt(reduced_mass_air / reduced_mass_gas)

  return air_to_gas
end


--[[
  Estimates ion diameter d_ion (nm) from mass mass_ion (u) and
  (optionally) mobility (ko).  Given diameter d_gas (nm) and mass
  mass_gas (u) of gas.
--]]
local function estimate_d_ion(mass_ion,ko, d_gas,mass_gas)
  -- Get rough estimate of ion diameter (in nm) from mass_ion (u).
  -- Emperical formula noted in paper.
  local d_ion = 0.120415405 * mass_ion^(1/3)

  if ko then -- refinement (if ko given)
    -- Compute Koair by scaling Ko(gas).
    local Koair = ko / get_air_to_gas(d_ion,mass_ion, d_gas,mass_gas)

    -- Estimate d (nm) from Ko (10-4 m2 V-1 s-1) from ko.
    local C = 1.0e5  -- converts to (10-9 m2 V-1 s-1)
    local logkm = log10(Koair * C)
    -- Emperical formula noted in the paper.
    d_ion = 10^(3.0367 - 0.8504 * logkm + 0.1137 * logkm^2 - 0.0135 * logkm^3)
  end

  return d_ion
end


--[[
  Estimates ion mobility (ko) from mass mass_ion (u) and diameter
  d_ion (nm) of ion.  Given diameter d_gas (nm) and mass mass_gas
  (u) of gas.
--]]
local function estimate_ko(mass_ion,d_ion, d_gas,mass_gas)
  -- Estimate Koair (10-4 m2 V-1 s-1) from d_ion
  -- Emperical formula noted in paper.
  local logdm = log10(d_ion)
  local Koair =
    1.0e-5 * 10^(4.9137 - 1.4491*logdm - 0.2772*logdm^2 + 0.0717*logdm^3)

  -- Compute Kogas by scaling Koair (10-4 m2 V-1 s-1)
  local ko = Koair * get_air_to_gas(d_ion,mass_ion, d_gas,mass_gas)

  return ko
end


--[[
  Adds ions STP average velocity Vo and mean free path MFPo data to
  given mass record mdata.  The data is computed from the other data
  in the mass record.  Given diameter d_gas (nm) and mass mass_gas
  (u) of gas.
  returns mfpmass, vmass (the Vo and MFPo values added).
--]]
local function massdata_add_Vo_MFPo(mdata, d_gas,mass_gas)
  local mass_ion, d_ion, ko = mdata.mass, mdata.d, mdata.ko

  -- Mean ion speed at STP (mm/usec)
  local MM_USEC__M_S = 1.0e-3  -- (mm/usec)/(m/s)
  local Vk = sqrt(8 * K_BOLTZMANN * STP_TEMP / PI / AMU_TO_KG) * MM_USEC__M_S
  -- print('DEBUG:check: ', Vk, '=', 2.404850043)

  -- Thermal velocity of ion (mm/usec)
  local Vio = Vk * sqrt(1 / mass_ion)

  -- Thermal velocity of gas molecule (mm/usec)
  local Vgo = Vk * sqrt(1 / mass_gas)

  -- Particles per volume, No, (n/mm3) at STP.
  local MM3_PER_M3 = 1.0e9     -- (mm3/m3)
  local No = N_AVOGADRO        -- (n/mol)
           / MOL_VOLUME        -- (m3/mol)
           / MM3_PER_M3

  -- Collision frequency (collisions/usec)
  -- Collision frequency constant Fk 
  local MM2_PER_NM2 = 1.0e-12  -- (mm2/nm2)
  local Fio = MM2_PER_NM2 *  No * PI *
    ((sqrt(2)-1/4) * ((d_gas + d_ion)/2)^2 * Vio + (1/4) * d_ion^2 * Vgo)

  -- Ion's mean free path (mm) at STP
  local Lio = Vio / Fio

  local vmass   = Vio      -- ion's STP mean velocity (mm/usec), Vo
  local mfpmass = Lio      -- ion's STP mean free path, MFPo

  -- Store Vo and MFPo in mass record.
  mdata.vo   = vmass
  mdata.mfpo = mfpmass

  return mfpmass, vmass
end


--[[
  Gets data (MFPo, Vo, and Ko) for given mass mass_ion
  (u) from massdata.  Fills in mass record with estimated values if
  no match.
--]]
local function massdata_for_mass(massdata, mass_ion)
  -- Return any record that matches.
  for i,mdata in ipairs(massdata) do
    if mass_ion == mdata.mass then
      return mdata.mfpo, mdata.vo, mdata.ko
    end
  end
  -- Otherwise, estimate values...

  local d_gas, mass_gas = massdata.d_gas, massdata.mass_gas

  -- Estimate ion diameter (nm) from mass_ion
  local d_ion = estimate_d_ion(mass_ion,nil, d_gas,mass_gas)

  -- Estimate Ko (10-4 m2 V-1 s-1) from mass_ion and d_ion.
  local ko = estimate_ko(mass_ion,d_ion, d_gas,mass_gas)

  -- Create new mass record.
  local mdata = {
    mass       = mass_ion,
    d          = d_ion,
    ko         = ko,
    estimation = 'dk'
  }
  -- add Vo and MFPo also to mdata
  local mfpmass, vmass = massdata_add_Vo_MFPo(mdata, d_gas,mass_gas)
  massdata[#massdata+1] = mdata

  return mfpmass, vmass, ko
end


--[[
  Fills in missing data in massdata mass records, after diameter d_gas
  (nm) and mass mass_gas (u) of background gas is known.
--]]
local function massdata_complete_records(massdata, d_gas,mass_gas)
  -- Preserve.
  massdata.d_gas, massdata.mass_gas = d_gas,mass_gas

  -- Update each record.
  for i,mdata in ipairs(massdata) do
    local mass_ion, d_ion, ko = mdata.mass, mdata.d, mdata.ko

    -- Compute d or Ko if necessary.
    -- Update mass record.
    if ko * d_ion ~= 0 then     -- both defined
      mdata.estimation = ''
    elseif d_ion == 0 and ko == 0 then -- both undefined
      mdata.d = estimate_d_ion(mass_ion,nil, d_gas,mass_gas)  -- ko = nil
      mdata.ko = estimate_ko(mass_ion,mdata.d, d_gas,mass_gas)
      mdata.estimation = 'dk'
    elseif d_ion == 0 then      -- d_ion undefined
      mdata.d = estimate_d_ion(mass_ion,ko, d_gas,mass_gas)
      mdata.estimation = 'd'
    elseif ko == 0 then         -- Ko undefined
      mdata.ko = estimate_ko(mass_ion,d_ion, d_gas,mass_gas)
      mdata.estimation = 'k'
    end

    -- Add Vo and Mean free path to record.
    massdata_add_Vo_MFPo(mdata, d_gas,mass_gas)
  end
end


--[[
  Prints data for all mass records.
--]]
local function massdata_print(massdata, tlocal,plocal)
  -- Print mass data.

  -- sort data by increasing mass
  local mass_idxs = {}
  for i=1,#massdata do mass_idxs[i] = i end
  table.sort(mass_idxs,
    function(a,b) return massdata[a].mass < massdata[b].mass end)

  -- Print mass data.
  print()
  print("SDS Mass definitions loaded or created (sorted by mass)")
  print("  N , mass, dia,  Ko Mobility  ,MFPo, Avg Vo")
  print(" (n),(amu),(nm),(10-4 m2V-1s-1),(mm),(mm/usec)")
  local missing = false
  for _,i in ipairs(mass_idxs) do
    local mdata = massdata[i]
    -- add line to new mass data output
    printf(
      " n=%d, m=%g, d=%g(%s), Ko=%g(%s), MFPo=%g, Vo=%g %s", 
      i,
      mdata.mass,
      mdata.d,
      (mdata.estimation:find'd' and "est" or "def"),
      mdata.ko,
      (mdata.estimation:find'k' and "est" or "def"),
      mdata.mfpo,
      mdata.vo,
      mdata.estimation ~= '' and '[*]' or ''
    )
    missing = missing or mdata.estimation ~= ''
  end
  if missing then
    print(" [*] WARNING: Estimations made due to absent data in m_defs.dat")
  end

  -- Print mass data at local conditions.
  if plocal ~= 0 then
    printf("")
    printf("SDS Local values for ions (sorted by mass) at user defined")
    printf("   Pressure = %g Torr, Temperature = %g K", plocal, tlocal)
    printf(" mass ,  K Mobility   , MFP,  Avg V  ,B1 MFP/rB1=1")
    printf(" (amu),(10-4 m2V-1s-1),(mm),(mm/usec),(gauss)")
    for _,i in ipairs(mass_idxs) do
      local mdata = massdata[i]

      -- Correct for local temperature and pressure.
      local vlocal   = mdata.vo * sqrt(tlocal / STP_TEMP)
      local mfplocal = mdata.mfpo * (760 / plocal) * (tlocal / STP_TEMP)
      local klocal   = mdata.ko   * (760 / plocal) * (tlocal / STP_TEMP)

      -- Magnetic field B when MFP = r_cyclotron (thermal).
      local B1 = 1439.74 *
                 sqrt(mdata.mass * speed_to_ke(vlocal, mdata.mass)) / mfplocal

      printf(" m=%g, K=%g, MFP=%g, v=%g, B1=%g",
             mdata.mass, klocal, mfplocal, vlocal, B1)
    end
  end

  print()
end


--[[
  Loads mass definitions (ion diameters and mobility constants) keyed to mass
  from file with name filename.
  
  The file must have the format as described in read_file_numbers,
  where each row contains a comma-separated list of these numbers:
  
    mass: ion mass (u)
    d:    ion diameter d (nm)
    ko:   ion reduced mobility Ko (10-4 m2 V-1 s-1) (at STP).

  This routine creates an array of mass records.  Records are tables
  with the above fields.
  
  Records can also contain these fields (after calling
  massdata_complete_records):
  
    vo       -- ion velocity Vo (mm/usec)
    mfpo     -- ion MFPo (mm)
    estimation -- which of ion d and Ko are estimated rather than defined:
                  '', 'd', 'k', 'dk'
  
  The containing array can also contain these fields (after calling
  massdata_complete_records):
  
    d_gas    -- background gas diameter (nm)
    mass_gas -- background gas mass (u)
  
  returns empty list if file not exist.
--]]
local function load_massdata(filename)
  local raw = TF.opt_read_file_numbers(filename)
  local massdata = {}
  -- Check mass records.
  if raw then
    local ncols = raw.ncols
    if (ncols ~= 3 or ncols ~= 5) and #raw % ncols ~= 0 then
      error(filename .. " does not contain three or five numbers per line.")
    end
    local ir, i = 1, 1
    while raw[ir] and raw[ir] ~= 0 do
      massdata[i] = {
        mass  = abs(raw[ir]),
        d     = abs(raw[ir + 1]),
        ko    = abs(raw[ir + 2]),
        alpha = ncols >= 5 and raw[ir + 3] or nil,
        beta  = ncols >= 5 and raw[ir + 4] or nil
      }
      ir = ir + ncols
      i = i + 1
    end
    -- check for duplicated masses
    local found = {}
    for _,mdata in ipairs(massdata) do
      local mass = mdata.mass
      if found[mass] then error("SDS ERROR: Mass duplicated: #", mass) end
      found[mass] = true
    end
  else
    printf("SDS Warning: mass definitions \"%s\" not found", filename)
  end
  return massdata
end


--##
--## SECTION: More Locals
--##


-- Collision statistics data.
local s_dist = load_diffusion_statistics(find_file("mbmr.dat"))

-- Data for diffusion and mobility parameters keyed to mass.
local massdata = load_massdata(find_file("m_defs.dat"))

-- User supplied initialization function.
-- If not nil, this function is called at the beginning of SDS initialization.
-- A typical application of this is to intialize SDS adjustable variables
-- based on other adjustable variables. Example:
--
--   adjustable pressure_pascals = 101325
--   function SDS.init()
--     adjustable SDS_pressure_torr = pressure_pascals * (760/101325)
--   end
local init = nil


-- will be set to true if any local gas parameters are defined.
local is_local_params = false

-- will be set to one of the arrays to use a template for dimensions.
local array_template


---- Data for each ion.
---- For each variable t below, t[i] is the data for ion number i.
---- Note that STP standard for "standard temperature and pressure".
--
-- Ion's Stokes' law damping (usec-1) at STP.
-- This is converted from ion mobility.
--
local ions_STP_damping = {}
--
-- Ion's mean free path (mm) at STP.
--
local ions_STP_mfp_mm = {}
--
-- Ion's average thermal velocity (mm/usec) at STP.
--
local ions_STP_Vo_mm_per_usec = {}
--
-- Ion's Stoke's law damping (usec-1) at the local temperature and
-- pressure.  This is calculated in part from ions_STP_damping.
--
local ions_local_damping = {}
--
-- Ion's local pressure, temp corrected mean free path (mm).
-- This is calculated in part from ions_STP_mfp_mm.
--
local ions_local_mfp_mm = {}
M.ions_local_mfp_mm = ions_local_mfp_mm  -- used by RS
--
-- Ion's average thermal velocity rate (mm/usec) at local temperature.
-- This is calculated in part from ions_STP_Vo_mm_per_usec.
--
local ions_local_V_mm_per_usec = {}
--
-- Ion's log10(mass_ion/mass_gas), where mass_ion/mass_gas is ratio of
-- ion mass to gas mass.
--
local ions_log_mr_ratio = {}



-- System parameters at current ion position.
-- These apply to the current ion.
local local_pressure_torr   -- local pressure (torr)
local local_temperature_K   -- local temperature (K)
local local_vx_mm_per_usec  -- local bulk gas velocity components (mm/usec)
local local_vy_mm_per_usec  -- ..
local local_vz_mm_per_usec  -- ..


-- The user program is enabled when instances[ion_instance] = true.
-- Set to false to disable SDS on certain PA instance numbers
-- (e.g. when electrostatic and magnetic PAs overlap, you only want
-- SDS to be active only on the magnetic PA).
local instances = {}
for i=1,200 do instances[i] = true end -- all instances enabled by default.


--##
--## SECTION: SDS Util
--##


--[[
  Updates local data of ion number `nion`:
    ions_local_damping[nion]
    ions_local_mfp_mm[nion]
    ions_local_V_mm_per_usec[nion]
  for local state (e.g. temperature and pressure):
    local_temperature_K
    local_pressure_torr
  as well as ion's STP variables:
    ions_STP_mfp_mm[nion]
    ions_STP_Vo_mm_per_usec[nion]
--]]
local function update_ions_local(nion)
  local mfpmass = ions_STP_mfp_mm[nion]
  local vmass   = ions_STP_Vo_mm_per_usec[nion]

  -- Correct for local temperature and pressure.
  local t_ratio = local_temperature_K / STP_TEMP          -- local T
  local pt_ratio = t_ratio * (760 / local_pressure_torr)  -- local P
  ions_local_damping[nion] = ions_STP_damping[nion] / pt_ratio
  ions_local_mfp_mm[nion] = mfpmass * pt_ratio
  ions_local_V_mm_per_usec[nion] = vmass * sqrt(t_ratio)
end



--##
--## SECTION: SIMION Adjustable Variables (used in segments).
--##


-- Whether SDS effects are enabled (1=yes, 0=no).
-- Allows quickly comparing with and without collisional effects.
adjustable SDS_enable = 1

-- Whether diffusion effects enabled.
-- 0=no, 1=yes (default yes).
-- It can sometimes be helpful to disable this to observe the behavior
-- without diffusion effects, such as to observe the average expected
-- behavior using a single ion.
adjustable SDS_diffusion = 1

-- Mass of background gas particle (u).  Default assumes normal air
-- mixture.
adjustable SDS_collision_gas_mass_amu = 28.94515

-- Effective diameter of background gas molecules (nm).  Default
-- assumes air.
adjustable SDS_collision_gas_diameter_nm = 0.366

-- Background gas pressure (torr).  Default is 760 Torr = 101,325 Pa
-- (atmospheric pressure).
-- This will be overriden by pres_defs (if defined)
adjustable SDS_pressure_torr = 2  --760.00

-- Background gas temperature (K).  Default is standard temp (25 C).
-- This will be overriden by temp_defs (if defined)
adjustable SDS_temperature_K = 298.15

-- Background gas mean velocity in x (m s-1).  Default is stationary.
-- This will be overriden by vel*_defs (if any defined)
adjustable SDS_vx_m_per_sec  = 0.0

-- Background gas mean velocity in y (m s-1).  Default is stationary.
-- This will be overriden by vel*_defs (if any defined)
adjustable SDS_vy_m_per_sec  = 0.0

-- Background gas mean velocity in z (m s-1).  Default is stationary.
-- This will be overriden by vel*_defs (if any defined)
adjustable SDS_vz_m_per_sec  = 0.0

-- Scaling factor to apply on pressure array fields (pres_defs) if
-- defined.  Default is 1 (no scaling).
adjustable SDS_P_field_scale_factor = 1.0

-- Scaling factor to apply on temperature array fields (temp_defs) if
-- defined.  Default is 1 (no scaling).
adjustable SDS_T_field_scale_factor = 1.0

-- Scaling factor to apply on velocity (Vx, Vy, Vz) array fields
-- (vel*_defs) if defined.  Default is 1 (no scaling).
adjustable SDS_v_field_scale_factor = 1.0

-- Force minimum time step in usecs.
-- 0 (the default) disables this. >0 forces code speedup.
adjustable SDS_min_time_step_usec = 0.0



-- default functions
local defaults = {}
function defaults.vx() return SDS_vx_m_per_sec end
function defaults.vy() return SDS_vy_m_per_sec end
function defaults.vz() return SDS_vz_m_per_sec end
function defaults.v()  return SDS_vx_m_per_sec,
                              SDS_vy_m_per_sec,
                              SDS_vz_m_per_sec  end
function defaults.p()  return SDS_pressure_torr end
function defaults.t()  return SDS_temperature_K end


-- Definitions for pressure, temperature, and background gas velocity.
-- These consist of an object and a function derived from the object.
-- Object may be `nil` if global state value is used everywhere.
--
-- pressure function (torr)
--
local pres_obj  = find_instance'p' or load_array("p_defs.dat") or nil
local pres_defs = wrap_scalar_field(pres_obj or defaults.p)
--
-- temperature function (deg K)
--
local temp_obj  = find_instance't' or load_array("t_defs.dat") or nil
local temp_defs = wrap_scalar_field(temp_obj or defaults.t)
    
-- bulk gas velocity function (m/s) for all three components
-- (vx,vy,vz).
-- Note: vx,vy,vz, if defined via a function, are oriented with respect
-- to the workbench (not the PA), which is different from the above.
--
local vel_obj, vel_defs do
  local vxobj, vyobj, vzobj = find_instance 'vx', find_instance 'vy', find_instance 'vz'
  if not(vxobj or vyobj or vzobj) then
    vxobj,vyobj,vzobj = load_array("vx_defs.dat"),load_array("vy_defs.dat"),load_array("vz_defs.dat")
  end
  if vxobj or vyobj or vzobj then
    vel_obj = {vxobj,vyobj,vzobj}
  end
  vel_defs = wrap_vector_field(vel_obj or defaults.v, defaults)
end
local function velx_defs(...) local vx,_,_ = vel_defs(...); return vx end
local function vely_defs(...) local _,vy,_ = vel_defs(...); return vy end
local function velz_defs(...) local _,_,vz = vel_defs(...); return vz end


-- (See section L.1.6 "Units and Coordinate Conventions" of
-- the printed SIMION manual for details.)
--  "mm"     -- workbench coordinates
--  "gu"     -- array volume coordinates
--  "abs_gu" -- absolute array coordinates (e.g removing symmetry)


--##
--## SECTION: SIMION Segments (or Code Designed to Run in them)
--##


--[[
  Prints coordinates of ion and local state (e.g. pressure, temperature,
  and velocity of gas).  Useful for debugging.
 
  This is designed to be called in certain SIMION segments, such as
  initialize, other_actions, terminate, or accel_adjust.
--]]
local function print_coordinates()
  printf("-----")
  printf("  Ion at array coords (gu): x=%f, y=%f, z=%f",
         ion_px_abs_gu, ion_py_abs_gu, ion_pz_abs_gu)
  printf("  Ion at PA coords (gu): x=%f, y=%f, z=%f",
         ion_px_gu, ion_py_gu, ion_pz_gu)
  printf("  Ion at workbench coords (mm): x=%f, y=%f, z=%f",
         ion_px_mm, ion_py_mm, ion_pz_mm)
  printf("  Bulk gas vel (m/s): vx=%f, vy=%f, vz=%f",
         local_vx_mm_per_usec * 1.0e3, local_vy_mm_per_usec * 1.0e3,
         local_vz_mm_per_usec * 1.0e3)  -- scale to m/s
  printf("  Local P(torr)=%f, T(K)=%f",
         local_pressure_torr, local_temperature_K)
end


--[[
  Print SDS parameters to log.
  (useful to ensure desired parameters are active)
--]]
local function print_sds_parameters()
  print("SDS collision gas diameter (nm)=", SDS_collision_gas_diameter_nm)
  print("SDS collision gas mass (u)=", SDS_collision_gas_mass_amu)
  
  local function mytostring(o)
    local s = (type(o) == 'function') and 'Lua function' or objtype(o)
    if (type(o)=='table' or type(o)=='userdata') and o.filename then s = s..' '..o.filename end
    return s
  end

  local s= "SDS pressure (torr)= "
  if pres_obj then
    s = s .. mytostring(pres_obj)
    if SDS_P_field_scale_factor ~= 1 then
      s = s .. ", scaled by factor=" .. SDS_P_field_scale_factor
    end
  else
    s = s .. SDS_pressure_torr
  end
  print(s)

  local s = "SDS temperature (K)= "
  if temp_obj then
    s = s .. mytostring(temp_obj)
    if SDS_T_field_scale_factor ~= 1 then
      s = s .. ", scaled by factor=" .. SDS_T_field_scale_factor
    end
  else
    s = s .. SDS_temperature_K
  end
  print(s)

  local s = "SDS velocity (m/s)= "
  if vel_obj then
    local vel_table = (type(vel_obj) == 'table') and vel_obj or {non=true}
    s = s ..
      (vel_table.non  and "\n vxyz=" .. mytostring(vel_obj)  or '') ..
      (vel_table[1] and "\n vx="   .. mytostring(vel_table[1]) or '') ..
      (vel_table[2] and "\n vy="   .. mytostring(vel_table[2]) or '') ..
      (vel_table[3] and "\n vz="   .. mytostring(vel_table[3]) or '')
    if SDS_v_field_scale_factor ~= 1 then
      s = s .. "\n scaled by factor=" .. SDS_v_field_scale_factor
    end
  else
    s = s .. SDS_vx_m_per_sec .. ', ' ..  SDS_vy_m_per_sec .. ', '
          .. SDS_vz_m_per_sec
  end
  print(s)

  if not is_local_params then
    print("SDS Note: No local P, T, or v field defined.")
  end

  if SDS_min_time_step_usec ~= 0 then
    print("SDS min time step (usec)=", SDS_min_time_step_usec)
  end
end




-- Gets current ion position in given coordinate system.
local getcoords = {}
function getcoords.abs_gu() return ion_px_abs_gu, ion_py_abs_gu, ion_pz_abs_gu end
function getcoords.gu() return ion_px_gu, ion_py_gu, ion_pz_gu end
function getcoords.mm() return ion_px_mm, ion_py_mm, ion_pz_mm end


--[[
  Checks data arrays: check flow arrays.
  Has side-effects.
--]]
local function check_all_arrays()
  -- Check all arrays (these can be nil)
  local more = (objtype(vel_obj)=='rawtable') and vel_obj or {vel_obj}
  local objs = {pres_obj, temp_obj, unpack(more, 1, 3)}
  for _,obj in pairs(objs) do
    if obj then is_local_params = true end
    if objtype(obj) == 'luaarray' or objtype(obj) == 'pa' then
      if not array_template then array_template = obj end
      if obj.nx ~= array_template.nx then
         error("ERROR: nx dimension of arrays don't match") end
      if obj.ny ~= array_template.ny then
         error("ERROR: ny dimension of arrays don't match") end
      if obj.nz ~= array_template.nz then
         error("ERROR: nz dimension of arrays don't match") end
    end
  end
end


--[[
  Updates local parameters (pressure, temperature, velocity, damping,
  MFP, and average speed) with data at current ion position.
 
  This is designed to be called from a SIMION accel_adjust segment.
 
  Writes local_* variables and also variables written by update_ions_local.
--]]
local function update_local_parameters()
  -- Compute local pressure.
  if pres_obj then
    local px,py,pz = getcoords.mm()
    local_pressure_torr = pres_defs(px,py,pz) * SDS_P_field_scale_factor
  end
  if not(local_pressure_torr >= 0) then error("local pressure < 0") end

  -- Compute local temperature.
  if temp_obj then
    local px,py,pz = getcoords.mm()
    local_temperature_K = temp_defs(px,py,pz) * SDS_T_field_scale_factor
  end
  if not(local_temperature_K > 0) then error("local temperature <= 0") end

  -- Need to recalculate ion's local parameters.
  if (pres_obj or temp_obj) and local_pressure_torr ~= 0 then
    update_ions_local(ion_number)
  end

  -- Compute local velocity.
  if vel_obj then
    local vx,vy,vz = vel_defs(getcoords.mm())
    local MM_USEC__M_SEC = 1.0e-3  -- (mm/usec)/(m/s)
    local f = SDS_v_field_scale_factor * MM_USEC__M_SEC
    local_vx_mm_per_usec = vx * f
    local_vy_mm_per_usec = vy * f
    local_vz_mm_per_usec = vz * f
  end
end


-- Generates vector of length r and uniformly random direction.
-- r defaults to 1 if omitted.
-- From simionx.Statistics sphere_rand.
local function sphere_rand(r)
  r = r or 1
  -- Marsaglia, Ann. Math. Stat. 43,645-646 (1972) doi:10.1214/aoms/1177692644 .
  -- Knop, CACM 13,5,326 (1970) doi:10.1145/362349.362377 .
  -- Sample z and azimuthal angle (around z) each from uniform distributions.
  -- First select point (x',y') uniformly within unit circle.
  -- x'^2+y'^2 has a uniform distribution used to define z.
  -- (x',y') is uniformly distributed in azimuthal angle and is rescaled
  -- to define (x,y) such that x^2+y^2+z^2 = r^2.
  local xp,yp
  local S
  repeat
    xp = 2*rand()-1
    yp = 2*rand()-1
    S = xp*xp + yp*yp
  until S <= 1  -- rejection method
  local z = (2*S - 1)*r
  local f = 2*r*sqrt(1-S)  -- rescaling factor for x,y
  local x, y = xp*f, yp*f
  return x,y,z
end


--[[
  Applies (random-walk) diffusion effect to current ion's position
  (ion_p[xyz]_mm) in this time-step, given time step delta_time (usec),
  mass ratio ions_mr_ratio (= mass_ion/mass_gas), mean-free-path
  ions_MFP (mm), and mean thermal velocity ions_V (mm/usec) of ion
  in workbench orientation.
  Uses diffusion statistics (s_dist).
 
  This is designed to be called inside a SIMION other_actions segment.
--]]
local function apply_diffusion(ions_mr_ratio, ions_MFP, ions_V, delta_time)
  -- Unitless distance traveled in N_DIST_COLLISIONS collisions,
  -- assuming normalized conditions.
  local dist_steps = diff_dist_steps(s_dist, ions_mr_ratio)

  -- Estimate number of collisions in current time step (delta_time),
  -- assuming local T and P (i.e. current MFP).
  local ncollisions =
      ions_V             -- (mm/usec) ion's current average thermal velocity
      / ions_MFP         -- (1/mm)    ion's current MFP
      * delta_time       -- (usec)    next time step

  -- Estimate jump distance (r) obtained over ncollisions collisions,
  -- assuming local T and P.  This is a simple scaling of dist_steps
  -- to actual conditions.
  local r =
    sqrt(ncollisions / s_dist.N_DIST_COLLISIONS) -- square root scaling law
    * dist_steps         -- r assuming MFP = 1.0 mm
    * ions_MFP           -- scale r with ion's MFP at local T and P

  -- Jump displacement vector (dx,dy,dz) in uniformly random direction.
  local dx,dy,dz = sphere_rand(r)

  -- Add jump vector to current ion position.
  ion_px_mm = ion_px_mm + dx
  ion_py_mm = ion_py_mm + dy
  ion_pz_mm = ion_pz_mm + dz
end


--[[
  Applies Stokes' Law viscous damping in this time step
  (ion_time_step) by damping acceleration (ion_a[xyz]_mm), given
  damping factor (usec-1) and mean bulk velocity of background gas
  (vx,vy,vz) (mm/usec in workbench orientation).
  
  This is designed to be called inside a SIMION accel_adjust segment.
--]]
local function apply_stokes_damping(damping, vx,vy,vz)
  -- See examples\drag\drag.lua for futher details on this implementation.
  if damping ~= 0 and ion_time_step ~= 0 then
    assert(damping > 0)
    local tterm = damping * ion_time_step  -- time constant
    local factor = (1 - exp(-tterm)) / tterm

    -- Store as new acceleration components.
    ion_ax_mm = factor*(ion_ax_mm - (ion_vx_mm - vx)*damping)
    ion_ay_mm = factor*(ion_ay_mm - (ion_vy_mm - vy)*damping)
    ion_az_mm = factor*(ion_az_mm - (ion_vz_mm - vz)*damping)
  end
end

--[[
 Updates mass and charge dependent variables.
 This is designed to be called inside initialize or other_actions segments.
--]]
function M.update_ion()
  -- Get or estimate parameters for current ion.
  local mfpmass, vmass, ko = massdata_for_mass(massdata, ion_mass)

  -- Estimate Stokes' Law damping (usec-1)
  local emu = ELEMENTARY_CHARGE/AMU_TO_KG    -- (C e^-1 kg-1 u)
  local damping =
      emu
      * 0.01       -- 10-4 * (10+6 usec sec-1)
      / ko         -- (10-4 m2 V-1 s-1) = 10-4 s C kg-1
      * (ion_charge / ion_mass)   -- (e/u)

  -- Compute log of mass ratio.
  local logmrratio = log10(ion_mass / SDS_collision_gas_mass_amu)

  -- Set current ion's parameters at standard temperature pressure (STP)
  -- or independent of position.
  local nion = ion_number
  ions_log_mr_ratio[nion]       = logmrratio
  ions_STP_damping[nion]        = damping
  ions_STP_mfp_mm[nion]         = mfpmass
  ions_STP_Vo_mm_per_usec[nion] = vmass

  -- Compute initial values of current ion's parameters at current position.
  update_ions_local(ion_number)
end
--[[
 Updates ion mass and mass dependent variables.
 This is designed to be called inside initialize or other_actions segments.
--]]
function M.update_ion_mass(new_mass)
  ion_mass = new_mass
  M.update_ion()
end
--[[
 Updates ion charge and charge dependent variables.
 This is designed to be called inside initialize or other_actions segments.
--]]
function M.update_ion_charge(new_charge)
  ion_charge = new_charge
  M.update_ion()
end

--[[
  SIMION segment called on IOB loading.
--]]
function M.segment.load()
  if sim_trajectory_quality > 0 then
    print("SDS changing Trajectory Quality (TQual) to zero (recommended for speed).")
    sim_trajectory_quality = 0
  end
end

--[[
  SIMION segment called once on beginning of run.
--]]
function M.segment.initialize_run()
  if SDS_enable == 0 then
    if ion_run == 1 then print "SDS disabled" end
    return
  end

  -- user initialization code
  if init then init() end

  -- checks
  -- Check adjustable variables.
  assert(SDS_pressure_torr >= 0, "SDS_pressure_torr not >= 0")
  assert(SDS_temperature_K > 0, "SDS_temperature_K not >= 0")
  assert(SDS_collision_gas_mass_amu~=0,"SDS_collision_gas_mass_amu is zero")
  assert(SDS_collision_gas_diameter_nm ~= 0,
         "SDS_collision_gas_diameter_nm is zero")
  check_all_arrays()

  -- Fill in mass data records, given adjusted vars on background gas.
  massdata_complete_records(massdata,
    SDS_collision_gas_diameter_nm, SDS_collision_gas_mass_amu)

  -- printing
  if ion_run == 1 then
    if SDS_diffusion ~= 0 then
      print "SDS enabled (diffusion enabled)"
    else
      print "SDS enabled (WARNING: diffusion disabled: SDS_diffusion==0)"
    end
    print_sds_parameters()
  end

  --seed(1) --DEBUG (disable randomization between runs)
  
  -- Set local temperature, pressure, and background gas velocity to
  -- global ones.  This is the default unless arrays/functions are
  -- defined.
  local_temperature_K = SDS_temperature_K
  local_pressure_torr = SDS_pressure_torr
  local MM_USEC__M_S = 1.0e-3  -- (mm/usec)/(m/s)
  local_vx_mm_per_usec = SDS_vx_m_per_sec * MM_USEC__M_S
  local_vy_mm_per_usec = SDS_vy_m_per_sec * MM_USEC__M_S
  local_vz_mm_per_usec = SDS_vz_m_per_sec * MM_USEC__M_S
end

--[[
  SIMION segment called for each particle creation.
 
  This is used to initialize some parameters.
--]]
function M.segment.initialize()
  if SDS_enable == 0 then return end

  M.update_ion()
end


--[[
  SIMION segment called to override time step sizes.
--]]
function M.segment.tstep_adjust()
  if SDS_enable == 0 then return end

  -- Increase time step to minimum.
  ion_time_step = max(ion_time_step, SDS_min_time_step_usec)
end


--[[
  SIMION segment called to override acceleration vector.
  
  This is used to apply the viscous mobility effect.
--]]
function M.segment.accel_adjust()
  if SDS_enable == 0 then return end

  -- WARNING: This segment can be limited to only certain PA
  -- instances.  It's particularly important that if magnetic and
  -- electrostatic arrays overlap then this segment should only
  -- be called in the magnetic one.
  if not instances[ion_instance] then return end

  -- Update local_* parameters.
  -- (affect the Stokes' Law damping constant).
  if is_local_params then update_local_parameters() end

  -- stokes' law viscous mobility effect.
  if local_pressure_torr ~= 0 then
    -- Obtain Stokes' Law damping constant for current ion.
    local damping = ions_local_damping[ion_number]
    if not damping then
      error(string.format(
        "Particle #%d was not initialized.  Ensure that the particle " ..
        "originates strictly inside a potential array instance in which " ..
        "the initialize segment is called.", ion_number))
    end
    -- stokes' law viscous mobility effect.
    apply_stokes_damping(
      damping, local_vx_mm_per_usec,local_vy_mm_per_usec,local_vz_mm_per_usec)
  end
end


--[[
  SIMION segment called on every time-step.
  This is used to apply the random-walk diffusion effect.
--]]
local first = true
function M.segment.other_actions()
  if SDS_enable == 0 then return end

  -- WARNING: This segment can be limited to only certain PA instances.
  -- It's particularly important that if magnetic and electrostatic arrays
  -- overlap, then this segment should only be called in the magnetic one.
  if not instances[ion_instance] then return end

  -- Code executed on first call.
  if first then first = false
    massdata_print(massdata, SDS_temperature_K, SDS_pressure_torr)
  end
  
  -- random-walk diffusion effect.
  if SDS_diffusion ~= 0 and local_pressure_torr ~= 0 then
    apply_diffusion(ions_log_mr_ratio[ion_number],
                    ions_local_mfp_mm[ion_number],
                    ions_local_V_mm_per_usec[ion_number],
                    ion_time_step)
  end
  --print('time constant (usec)=', 1 / ions_local_damping[ion_number])
  -- print_coordinates()
end


local function merge_segments(t)
  for name,newseg in pairs(t) do
    local oldseg = segment[name]
    segment[name] =
      oldseg and function() oldseg(); newseg() end
             or  newseg
  end
end


--[[
  Installs SIMION segments (copy functions from M.segment table into
  segment table).
--]]
function M.install()
  merge_segments(M.segment)
end


-- Install segments on load.
if not opt_noinstall then
  M.install()
end



-- Catch getting and setting on this module table.
-- For example, support things like
--   local SDS = simion.import("collision_sds.lua")
--   SDS.pressure = function(x,y,z) return x * 2 end
--   SDS.pressure = function() return ion_px_mm * 2 end
--   SDS.pressure = simion.wb.instances[1]
--   SDS.velocity = function(x,y,z) return x*2,0,0 end
--   SDS.velocity = {simion.wb.instances[2], simion.wb.instances[3], simion.wb.instances[4])
local mt = {}
function mt.__index(t,k)  -- GET t[k]
  if     k == 'pressure'    then return pres_defs
  elseif k == 'temperature' then return temp_defs
  elseif k == 'velocity'    then return vel_defs
  elseif k == 'velocity_x'  then return velx_defs
  elseif k == 'velocity_y'  then return vely_defs
  elseif k == 'velocity_z'  then return velz_defs
  elseif k == 'velocity_coordinates'    then error('no longer implemented: ' .. k)
  elseif k == 'pressure_coords'         then error('no longer implemented: ' .. k)
  elseif k == 'temperature_coordinates' then error('no longer implemented: ' .. k)
  elseif k == 'instances'   then return instances
  elseif k == 'init'        then return init
  else error(k) end
end
function mt.__newindex(t,k,v)  -- SET t[k]=v
  print("SDS Defining " .. k .. ":",
      type(v) == 'function' and 'as Lua function' or v)
  if     k == 'pressure'    then pres_defs = wrap_scalar_field(v); pres_obj = v
  elseif k == 'temperature' then temp_defs = wrap_scalar_field(v); temp_obj = v
  elseif k == 'velocity'    then vel_defs  = wrap_vector_field(v); vel_obj  = v
  elseif k == 'velocity_x'              then error('no longer implemented: ' .. k)
  elseif k == 'velocity_y'              then error('no longer implemented: ' .. k)
  elseif k == 'velocity_z'              then error('no longer implemented: ' .. k)
  elseif k == 'velocity_coordinates'    then error('no longer implemented: ' .. k)
  elseif k == 'pressure_coordinates'    then error('no longer implemented: ' .. k)
  elseif k == 'temperature_coordinates' then error('no longer implemented: ' .. k)
  elseif k == 'instances'   then instances = v
  elseif k == 'init'        then init = v
  else error(k) end
end
setmetatable(M, mt)

_G.SDS = M -- allows plotting with contourlib81.lua.

return M
