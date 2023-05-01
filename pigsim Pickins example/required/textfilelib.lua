--[[
 testfilelib.lua - utilities for text files containing numerical data.
 
 D.Manura 2011-11-12.
 (c) 2007-2011 Scientific Instrument Services, Inc. (SIMION 8.1/8.0 Licensed)
--]]

local TF = {}


--[[
 Returns (true/false) whether the file of the given name exists
 (i.e. is readable at least).
--]]
function TF.file_exists(filename)
  local fh = io.open(filename)
  local is_exist = (fh ~= nil)
  if fh then fh:close() end
  return is_exist
end


--[[
 Reads numeric data file with name `filename` into array that is
 returned.  Commas and white-space are ignored.  Comments (which
 begin with a semicolon and continue to the end of the line) are
 ignored.  The table has an additional field `ncols` indicating
 the maximum number of columns on a line.
 Raises on error.
--]]
function TF.read_file_numbers(filename)
  local ncols = 0
  local array = {}
  local f = assert(io.open(filename, 'rb'))
  for line in f:lines() do
    line = line:match'^[^;\r]*' -- remove comments and \r
    local ncols_cur = 0
    local pos
    while 1 do
      local vstr, pos0 = line:match('^[ \t,]*([0-9%.eE+-]+)()', pos)
      if not vstr then
        if not line:match('^[ \t,]*$', pos) then
          error('unexpected ' .. line:sub(pos or 1), 2)
        end
        break
      end
      local v = tonumber(vstr)
      if v == nil then
        error(vstr .. ' not a number', 2)
      end
      array[#array+1] = v
      ncols_cur = ncols_cur + 1
      pos = pos0
    end
    ncols = math.max(ncols, ncols_cur)
  end
  array.ncols = ncols
  f:close()
  return array
end


--[[
 Same as `read_file_numbers`, except it just returns `nil`
 (not raises an error) if file doesn't exist.
--]]
function TF.opt_read_file_numbers(filename)
  return TF.file_exists(filename) and TF.read_file_numbers(filename) or nil
end


return TF
