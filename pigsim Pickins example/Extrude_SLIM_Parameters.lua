wElectrode = 0.4064    --Width of the DC electrode (mm)
lElectrode = 1.016    --Length of the DC electrode (mm)
wRF = 0.4064    --Width of the RF (mm)
wGuard = 3.1496    --Width of the Guard (mm)
sp = 0.127    --Spacing between electrodes (mm)
numDC = 8    --Number of DC electrodes per TW Segment
numRowsDC = 5    --Number of rows of TWAVE segments
numRF = 6    --Number of RF Segments
numGuards = 2    --Number of Guard Segments
totalLen = 39.852599999999995    --Total Length of the monomer (mm)
totalWidth = 67.1576    --Total Width of the monomer (mm)
trackWidth = 5.9944    --Track Width (mm)
d = 1.5875    -- half distance between the layers/boards (mm)
elecDict = {["DC1"] = 1, ["DC2"] = 2, ["DC3"] = 3, ["DC4"] = 4, ["DC5"] = 5, ["DC6"] = 6, ["DC7"] = 7, ["DC8"] = 8, ["RF1"] = 9, ["RF2"] = 10, ["G"] = 11}    -- elec keys : SIMION adj_elec dexes
design_title = 'pickins_test'    --  Name for associated files
splat_line = {38.77, 57.15, 39.79, 64.87}        ----- rectangle box defining a fly finish line
