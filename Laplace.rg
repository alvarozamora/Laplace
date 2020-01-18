import "regent"

-- Helper modules to handle PNG files and command line arguments
local Config = require("config")
local coloring   = require("coloring_util")

-- Some C APIs
local c     = regentlib.c
local sqrt  = regentlib.sqrt(double)
local cmath = terralib.includec("math.h")
local PI = cmath.M_PI

-- Field space for reduced distribution functions and their gradients
fspace grid{

  g : double,
  x : double,
  y : double,
  z : double
}

  reads(r_sig, r_mesh, plx_mesh, prx_mesh, plx_sig, prx_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 0  --change when copy
  var Dim2 : int32 = 0

  -- Cell Centers
  var xC : double
  var xL : double
  var xR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  -- Indices
  var bc : int32[6]
  var IL : int32 
  var JL : int32 
  var KL : int32 
  var IR : int32 
  var JR : int32 
  var KR : int32 
  
  var e3 : int8d 
  var eL3 : int8d 
  var eR3 : int8d 

  var e7 : int8d 
  var e8 : int8d 
  var eL7 : int8d 
  var eR7  : int8d

  -- Left/Right values of g/b
  var gsigL : double
  var gsigR : double
  var bsigL : double
  var bsigR : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    xC = r_mesh[e3].x -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IL = bc[0]
    JL = bc[1]
    KL = bc[2]
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]

    eL3 = {IL, JL, KL, 0, 0, 0, 0, 0}
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eL7 = {IL, JL, KL, Dim, 0, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
     
      if r_sig.bounds.lo.x == s.x then
        gsigL = plx_sig[eL7].g 
        bsigL = plx_sig[eL7].b 
        xL = plx_mesh[eL3].x
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        xL = r_mesh[eL3].x
      end

      if r_sig.bounds.hi.x == s.x then
        gsigR = prx_sig[eR7].g
        bsigR = prx_sig[eR7].b
        xR = prx_mesh[eR3].x
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        xR = r_mesh[eR3].x
      end

      r_sig2[e8].g = VanLeer(gsigL, r_sig[e7].g, gsigR, xL, xC, xR)
      r_sig2[e8].b = VanLeer(bsigL, r_sig[e7].b, bsigR, xL, xC, xR)
      --if s.y == 32 and v.x == 32 then c.printf("sigx_x[%d].g = %f\n", Dim2, r_sig2[e8].g) end
    end
  end
end

task Step1b_sigy_x(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            plx_mesh : region(ispace(int8d), mesh),
            prx_mesh : region(ispace(int8d), mesh),
            plx_sig : region(ispace(int8d), grid),
            prx_sig : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, plx_mesh, prx_mesh, plx_sig, prx_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 1  --change when copy
  var Dim2 : int32 = 0

  -- Cell Centers
  var xC : double
  var xL : double
  var xR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  -- Indices
  var bc : int32[6]
  var IL : int32 
  var JL : int32 
  var KL : int32 
  var IR : int32 
  var JR : int32 
  var KR : int32 
  
  var e3 : int8d 
  var eL3 : int8d 
  var eR3 : int8d 

  var e7 : int8d 
  var e8 : int8d 
  var eL7 : int8d 
  var eR7  : int8d

  -- Left/Right values of g/b
  var gsigL : double
  var gsigR : double
  var bsigL : double
  var bsigR : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    xC = r_mesh[e3].x -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IL = bc[0]
    JL = bc[1]
    KL = bc[2]
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]

    eL3 = {IL, JL, KL, 0, 0, 0, 0, 0}
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eL7 = {IL, JL, KL, Dim, 0, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
     
      if r_sig.bounds.lo.x == s.x then
        gsigL = plx_sig[eL7].g 
        bsigL = plx_sig[eL7].b 
        xL = plx_mesh[eL3].x
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        xL = r_mesh[eL3].x
      end

      if r_sig.bounds.hi.x == s.x then
        gsigR = prx_sig[eR7].g
        bsigR = prx_sig[eR7].b
        xR = prx_mesh[eR3].x
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        xR = r_mesh[eR3].x
      end

      r_sig2[e8].g = VanLeer(gsigL, r_sig[e7].g, gsigR, xL, xC, xR)
      r_sig2[e8].b = VanLeer(bsigL, r_sig[e7].b, bsigR, xL, xC, xR)

      --if s.x == 32 and v.x == 32 then c.printf("sigy_x[%d].g = %f\n", Dim2, r_sig2[e8].g) end
      if s.x == 32 and v.x == 32 and v.y == 32 then
        --c.printf("r_sigy_x[{%d, %d}, {%d, %d}]: gsigL = %f, r_gridbarp.g = %f, gsigR = %f\n", s.x, s.y, v.x, v.y, gsigL, r_sig[e7].g, gsigR)
      end

    end
  end
end

task Step1b_sigz_x(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            plx_mesh : region(ispace(int8d), mesh),
            prx_mesh : region(ispace(int8d), mesh),
            plx_sig : region(ispace(int8d), grid),
            prx_sig : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, plx_mesh, prx_mesh, plx_sig, prx_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 2  --change when copy
  var Dim2 : int32 = 0

  -- Cell Centers
  var xC : double
  var xL : double
  var xR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  -- Indices
  var bc : int32[6]
  var IL : int32 
  var JL : int32 
  var KL : int32 
  var IR : int32 
  var JR : int32 
  var KR : int32 
  
  var e3 : int8d 
  var eL3 : int8d 
  var eR3 : int8d 

  var e7 : int8d 
  var e8 : int8d 
  var eL7 : int8d 
  var eR7  : int8d

  -- Left/Right values of g/b
  var gsigL : double
  var gsigR : double
  var bsigL : double
  var bsigR : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    xC = r_mesh[e3].x -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IL = bc[0]
    JL = bc[1]
    KL = bc[2]
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]

    eL3 = {IL, JL, KL, 0, 0, 0, 0, 0}
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eL7 = {IL, JL, KL, Dim, 0, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
     
      if r_sig.bounds.lo.x == s.x then
        gsigL = plx_sig[eL7].g 
        bsigL = plx_sig[eL7].b 
        xL = plx_mesh[eL3].x
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        xL = r_mesh[eL3].x
      end

      if r_sig.bounds.hi.x == s.x then
        gsigR = prx_sig[eR7].g
        bsigR = prx_sig[eR7].b
        xR = prx_mesh[eR3].x
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        xR = r_mesh[eR3].x
      end

      r_sig2[e8].g = VanLeer(gsigL, r_sig[e7].g, gsigR, xL, xC, xR)
      r_sig2[e8].b = VanLeer(bsigL, r_sig[e7].b, bsigR, xL, xC, xR)
    end
  end
end

task Step1b_sigx_y(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            ply_mesh : region(ispace(int8d), mesh),
            pry_mesh : region(ispace(int8d), mesh),
            ply_sig : region(ispace(int8d), grid),
            pry_sig : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, ply_mesh, pry_mesh, ply_sig, pry_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 0  --change when copy
  var Dim2 : int32 = 1

  -- Cell Centers
  var yC : double
  var yL : double
  var yR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  -- Indices
  var bc : int32[6]
  var IL : int32 
  var JL : int32 
  var KL : int32 
  var IR : int32 
  var JR : int32 
  var KR : int32 
  
  var e3 : int8d 
  var eL3 : int8d 
  var eR3 : int8d 

  var e7 : int8d 
  var e8 : int8d 
  var eL7 : int8d 
  var eR7  : int8d

  -- Left/Right values of g/b
  var gsigL : double
  var gsigR : double
  var bsigL : double
  var bsigR : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    yC = r_mesh[e3].y -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IL = bc[0]
    JL = bc[1]
    KL = bc[2]
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]

    eL3  = {IL, JL, KL, 0, 0, 0, 0, 0}
    eR3  = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eL7 = {IL, JL, KL, Dim, 0, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
     
      if r_sig.bounds.lo.y == s.y then
        gsigL = ply_sig[eL7].g 
        bsigL = ply_sig[eL7].b 
        yL = ply_mesh[eL3].y
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        yL = r_mesh[eL3].y
      end

      if r_sig.bounds.hi.y == s.y then
        gsigR = pry_sig[eR7].g
        bsigR = pry_sig[eR7].b
        yR = pry_mesh[eR3].y
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        yR = r_mesh[eR3].y
      end

      r_sig2[e8].g = VanLeer(gsigL, r_sig[e7].g, gsigR, yL, yC, yR)
      r_sig2[e8].b = VanLeer(bsigL, r_sig[e7].b, bsigR, yL, yC, yR)

      --if s.x == 32 and v.x == 32 then c.printf("sigx_y[%d].g = %f\n", Dim2, r_sig2[e8].g) end
      --if s.x == 32 and v.x == 32 and v.y == 32 then
      --  c.printf("r_sigx_y[{%d, %d}, {%d, %d}]: gsigL = %f, r_sig.g = %f, gsigR = %f\n", s.x, s.y, v.x, v.y, gsigL, r_sig[e7].g, gsigR)
      --end
    end
  end
end

task Step1b_sigy_y(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            ply_mesh : region(ispace(int8d), mesh),
            pry_mesh : region(ispace(int8d), mesh),
            ply_sig : region(ispace(int8d), grid),
            pry_sig : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, ply_mesh, pry_mesh, ply_sig, pry_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 1  --change when copy
  var Dim2 : int32 = 1

  -- Cell Centers
  var yC : double
  var yL : double
  var yR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  -- Indices
  var bc : int32[6]
  var IL : int32 
  var JL : int32 
  var KL : int32 
  var IR : int32 
  var JR : int32 
  var KR : int32 
  
  var e3 : int8d 
  var eL3 : int8d 
  var eR3 : int8d 

  var e7 : int8d 
  var e8 : int8d 
  var eL7 : int8d 
  var eR7  : int8d

  -- Left/Right values of g/b
  var gsigL : double
  var gsigR : double
  var bsigL : double
  var bsigR : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    yC = r_mesh[e3].y -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IL = bc[0]
    JL = bc[1]
    KL = bc[2]
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]

    eL3 = {IL, JL, KL, 0, 0, 0, 0, 0}
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eL7 = {IL, JL, KL, Dim, 0, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
     
      if r_sig.bounds.lo.y == s.y then
        gsigL = ply_sig[eL7].g 
        bsigL = ply_sig[eL7].b 
        yL = ply_mesh[eL3].y
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        yL = r_mesh[eL3].y
      end

      if r_sig.bounds.hi.y == s.y then
        gsigR = pry_sig[eR7].g
        bsigR = pry_sig[eR7].b
        yR = pry_mesh[eR3].y
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        yR = r_mesh[eR3].y
      end

      r_sig2[e8].g = VanLeer(gsigL, r_sig[e7].g, gsigR, yL, yC, yR)
      r_sig2[e8].b = VanLeer(bsigL, r_sig[e7].b, bsigR, yL, yC, yR)

      --if s.x == 32 and v.x == 32 then c.printf("sigy_y[%d].g = %f\n", Dim2, r_sig2[e8].g) end

    end
  end
end

task Step1b_sigz_y(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            ply_mesh : region(ispace(int8d), mesh),
            pry_mesh : region(ispace(int8d), mesh),
            ply_sig : region(ispace(int8d), grid),
            pry_sig : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, ply_mesh, pry_mesh, ply_sig, pry_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 2  --change when copy
  var Dim2 : int32 = 1

  -- Cell Centers
  var yC : double
  var yL : double
  var yR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  -- Indices
  var bc : int32[6]
  var IL : int32 
  var JL : int32 
  var KL : int32 
  var IR : int32 
  var JR : int32 
  var KR : int32 
  
  var e3 : int8d 
  var eL3 : int8d 
  var eR3 : int8d 

  var e7 : int8d 
  var e8 : int8d 
  var eL7 : int8d 
  var eR7  : int8d

  -- Left/Right values of g/b
  var gsigL : double
  var gsigR : double
  var bsigL : double
  var bsigR : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    yC = r_mesh[e3].y -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IL = bc[0]
    JL = bc[1]
    KL = bc[2]
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]

    eL3 = {IL, JL, KL, 0, 0, 0, 0, 0}
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eL7 = {IL, JL, KL, Dim, 0, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
     
      if r_sig.bounds.lo.y == s.y then
        gsigL = ply_sig[eL7].g 
        bsigL = ply_sig[eL7].b 
        yL = ply_mesh[eL3].y
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        yL = r_mesh[eL3].y
      end

      if r_sig.bounds.hi.y == s.y then
        gsigR = pry_sig[eR7].g
        bsigR = pry_sig[eR7].b
        yR = pry_mesh[eR3].y
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        yR = r_mesh[eR3].y
      end

      r_sig2[e8].g = VanLeer(gsigL, r_sig[e7].g, gsigR, yL, yC, yR)
      r_sig2[e8].b = VanLeer(bsigL, r_sig[e7].b, bsigR, yL, yC, yR)
    end
  end
end

task Step1b_sigx_z(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            plz_mesh : region(ispace(int8d), mesh),
            prz_mesh : region(ispace(int8d), mesh),
            plz_sig : region(ispace(int8d), grid),
            prz_sig : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, plz_mesh, prz_mesh, plz_sig, prz_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 0  --change when copy
  var Dim2 : int32 = 2

  -- Cell Centers
  var zC : double
  var zL : double
  var zR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  -- Indices
  var bc : int32[6]
  var IL : int32 
  var JL : int32 
  var KL : int32 
  var IR : int32 
  var JR : int32 
  var KR : int32 
  
  var e3 : int8d 
  var eL3 : int8d 
  var eR3 : int8d 

  var e7 : int8d 
  var e8 : int8d 
  var eL7 : int8d 
  var eR7  : int8d

  -- Left/Right values of g/b
  var gsigL : double
  var gsigR : double
  var bsigL : double
  var bsigR : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    zC = r_mesh[e3].z -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IL = bc[0]
    JL = bc[1]
    KL = bc[2]
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]

    eL3 = {IL, JL, KL, 0, 0, 0, 0, 0}
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eL7 = {IL, JL, KL, Dim, 0, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
     
      if r_sig.bounds.lo.z == s.z then
        gsigL = plz_sig[eL7].g 
        bsigL = plz_sig[eL7].b 
        zL = plz_mesh[eL3].z
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        zL = r_mesh[eL3].z
      end

      if r_sig.bounds.hi.z == s.z then
        gsigR = prz_sig[eR7].g
        bsigR = prz_sig[eR7].b
        zR = prz_mesh[eR3].z
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        zR = r_mesh[eR3].z
      end

      r_sig2[e8].g = VanLeer(gsigL, r_sig[e7].g, gsigR, zL, zC, zR)
      r_sig2[e8].b = VanLeer(bsigL, r_sig[e7].b, bsigR, zL, zC, zR)
    end
  end
end

task Step1b_sigy_z(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            plz_mesh : region(ispace(int8d), mesh),
            prz_mesh : region(ispace(int8d), mesh),
            plz_sig : region(ispace(int8d), grid),
            prz_sig : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, plz_mesh, prz_mesh, plz_sig, prz_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 1  --change when copy
  var Dim2 : int32 = 2

  -- Cell Centers
  var zC : double
  var zL : double
  var zR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  -- Indices
  var bc : int32[6]
  var IL : int32 
  var JL : int32 
  var KL : int32 
  var IR : int32 
  var JR : int32 
  var KR : int32 
  
  var e3 : int8d 
  var eL3 : int8d 
  var eR3 : int8d 

  var e7 : int8d 
  var e8 : int8d 
  var eL7 : int8d 
  var eR7  : int8d

  -- Left/Right values of g/b
  var gsigL : double
  var gsigR : double
  var bsigL : double
  var bsigR : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    zC = r_mesh[e3].z -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IL = bc[0]
    JL = bc[1]
    KL = bc[2]
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]

    eL3 = {IL, JL, KL, 0, 0, 0, 0, 0}
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eL7 = {IL, JL, KL, Dim, 0, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
     
      if r_sig.bounds.lo.z == s.z then
        gsigL = plz_sig[eL7].g 
        bsigL = plz_sig[eL7].b 
        zL = plz_mesh[eL3].z
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        zL = r_mesh[eL3].z
      end

      if r_sig.bounds.hi.z == s.z then
        gsigR = prz_sig[eR7].g
        bsigR = prz_sig[eR7].b
        zR = prz_mesh[eR3].z
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        zR = r_mesh[eR3].z
      end

      r_sig2[e8].g = VanLeer(gsigL, r_sig[e7].g, gsigR, zL, zC, zR)
      r_sig2[e8].b = VanLeer(bsigL, r_sig[e7].b, bsigR, zL, zC, zR)
    end
  end
end

task Step1b_sigz_z(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            plz_mesh : region(ispace(int8d), mesh),
            prz_mesh : region(ispace(int8d), mesh),
            plz_sig : region(ispace(int8d), grid),
            prz_sig : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_mesh, plz_mesh, prz_mesh, plz_sig, prz_sig),
  reads writes(r_sig2)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 2  --change when copy
  var Dim2 : int32 = 2

  -- Cell Centers
  var zC : double
  var zL : double
  var zR : double

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  -- Indices
  var bc : int32[6]
  var IL : int32 
  var JL : int32 
  var KL : int32 
  var IR : int32 
  var JR : int32 
  var KR : int32 
  
  var e3 : int8d 
  var eL3 : int8d 
  var eR3 : int8d 

  var e7 : int8d 
  var e8 : int8d 
  var eL7 : int8d 
  var eR7  : int8d

  -- Left/Right values of g/b
  var gsigL : double
  var gsigR : double
  var bsigL : double
  var bsigR : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    zC = r_mesh[e3].z -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IL = bc[0]
    JL = bc[1]
    KL = bc[2]
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]

    eL3 = {IL, JL, KL, 0, 0, 0, 0, 0}
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eL7 = {IL, JL, KL, Dim, 0, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
     
      if r_sig.bounds.lo.z == s.z then
        gsigL = plz_sig[eL7].g 
        bsigL = plz_sig[eL7].b 
        zL = plz_mesh[eL3].z
      else
        gsigL = r_sig[eL7].g
        bsigL = r_sig[eL7].b
        zL = r_mesh[eL3].z
      end

      if r_sig.bounds.hi.z == s.z then
        gsigR = prz_sig[eR7].g
        bsigR = prz_sig[eR7].b
        zR = prz_mesh[eR3].z
      else
        gsigR = r_sig[eR7].g
        bsigR = r_sig[eR7].b
        zR = r_mesh[eR3].z
      end

      r_sig2[e8].g = VanLeer(gsigL, r_sig[e7].g, gsigR, zL, zC, zR)
      r_sig2[e8].b = VanLeer(bsigL, r_sig[e7].b, bsigR, zL, zC, zR)
    end
  end
end

task Step1b_sigx_x2(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_sigb : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            prx_mesh : region(ispace(int8d), mesh),
            prx_sig : region(ispace(int8d), grid),
            prx_sig2 : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_sig2, r_mesh, vxmesh, vymesh, vzmesh, prx_mesh, prx_sig, prx_sig2),
  reads writes(r_sigb)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 0  --change when copy
  var Dim2 : int32 = 0

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  var sC : double
 
  var bc : int32[6]
  var IR : int32
  var JR : int32
  var KR : int32
  
  var e3 : int8d
  var e7 : int8d
  var e8 : int8d
  var eR3 : int8d
  var eR7 : int8d
  var eR8 : int8d
 
  var gsig : double
  var bsig : double
  var gsig2 : double
  var bsig2 : double
  var swap : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    sC = r_mesh[e3].dx -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
      eR8 = {IR, JR, KR, Dim, Dim2, v.x, v.y, v.z}
  
      swap = 1.0
      if (vxmesh[v.x].v < 0 and Dim == 0) then
        if s.x == r_sigb.bounds.hi.x then
          gsig = prx_sig[eR7].g
          bsig = prx_sig[eR7].b
          gsig2 = prx_sig2[eR8].g
          bsig2 = prx_sig2[eR8].b
          sC = prx_mesh[eR3].dx
        else
          gsig = r_sig[eR7].g
          bsig = r_sig[eR7].b
          gsig2 = r_sig2[eR8].g
          bsig2 = r_sig2[eR8].b
          sC = r_mesh[eR3].dx
        end
        swap = -1
      else
        gsig = r_sig[e7].g
        bsig = r_sig[e7].b
        gsig2 = r_sig2[e8].g
        bsig2 = r_sig2[e8].b
        sC = r_mesh[e3].dx
      end
 
      r_sigb[e8].g = gsig + swap*(sC/2.0)*gsig2
      r_sigb[e8].b = bsig + swap*(sC/2.0)*bsig2
      
      --if v.x == 32 and s.x == 32 then c.printf("rsigb[%d] = {%f, %f}\n", s.x, r_sigb[e8].g, r_sigb[e8].b) end
    end
  end
end

task Step1b_sigy_x2(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_sigb : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            prx_mesh : region(ispace(int8d), mesh),
            prx_sig : region(ispace(int8d), grid),
            prx_sig2 : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_sig2, r_mesh, vxmesh, vymesh, vzmesh, prx_mesh, prx_sig, prx_sig2),
  reads writes(r_sigb)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 1  --change when copy
  var Dim2 : int32 = 0

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  var sC : double
 
  var bc : int32[6]
  var IR : int32
  var JR : int32
  var KR : int32
  
  var e3 : int8d
  var e7 : int8d
  var e8 : int8d
  var eR3 : int8d
  var eR7 : int8d
  var eR8 : int8d
 
  var gsig : double
  var bsig : double
  var gsig2 : double
  var bsig2 : double
  var swap : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    sC = r_mesh[e3].dx -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
      eR8 = {IR, JR, KR, Dim, Dim2, v.x, v.y, v.z}
  
      swap = 1.0
      if (vxmesh[v.x].v < 0 and Dim == 0) then
        if s.x == r_sigb.bounds.hi.x then
          gsig = prx_sig[eR7].g
          bsig = prx_sig[eR7].b
          gsig2 = prx_sig2[eR8].g
          bsig2 = prx_sig2[eR8].b
          sC = prx_mesh[eR3].dx
        else
          gsig = r_sig[eR7].g
          bsig = r_sig[eR7].b
          gsig2 = r_sig2[eR8].g
          bsig2 = r_sig2[eR8].b
          sC = r_mesh[eR3].dx
        end
        swap = -1
      else
        gsig = r_sig[e7].g
        bsig = r_sig[e7].b
        gsig2 = r_sig2[e8].g
        bsig2 = r_sig2[e8].b
        sC = r_mesh[e3].dx
      end
 
      r_sigb[e8].g = gsig + swap*(sC/2.0)*gsig2
      r_sigb[e8].b = bsig + swap*(sC/2.0)*bsig2
      
      --if v.x == 128 then c.printf("rsigb[%d] = {%f, %f}\n", s.x, r_sigb[e8].g, r_sigb[e8].b) end
    end
  end
end

task Step1b_sigz_x2(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_sigb : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            prx_mesh : region(ispace(int8d), mesh),
            prx_sig : region(ispace(int8d), grid),
            prx_sig2 : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_sig2, r_mesh, vxmesh, vymesh, vzmesh, prx_mesh, prx_sig, prx_sig2),
  reads writes(r_sigb)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 2  --change when copy
  var Dim2 : int32 = 0

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  var sC : double
 
  var bc : int32[6]
  var IR : int32
  var JR : int32
  var KR : int32
  
  var e3 : int8d
  var e7 : int8d
  var e8 : int8d
  var eR3 : int8d
  var eR7 : int8d
  var eR8 : int8d
 
  var gsig : double
  var bsig : double
  var gsig2 : double
  var bsig2 : double
  var swap : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    sC = r_mesh[e3].dx -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
      eR8 = {IR, JR, KR, Dim, Dim2, v.x, v.y, v.z}
  
      swap = 1.0
      if (vxmesh[v.x].v < 0 and Dim == 0) then
        if s.x == r_sigb.bounds.hi.x then
          gsig = prx_sig[eR7].g
          bsig = prx_sig[eR7].b
          gsig2 = prx_sig2[eR8].g
          bsig2 = prx_sig2[eR8].b
          sC = prx_mesh[eR3].dx
        else
          gsig = r_sig[eR7].g
          bsig = r_sig[eR7].b
          gsig2 = r_sig2[eR8].g
          bsig2 = r_sig2[eR8].b
          sC = r_mesh[eR3].dx
        end
        swap = -1
      else
        gsig = r_sig[e7].g
        bsig = r_sig[e7].b
        gsig2 = r_sig2[e8].g
        bsig2 = r_sig2[e8].b
        sC = r_mesh[e3].dx
      end
 
      r_sigb[e8].g = gsig + swap*(sC/2.0)*gsig2
      r_sigb[e8].b = bsig + swap*(sC/2.0)*bsig2
      
      --if v.x == 128 then c.printf("rsigb[%d] = {%f, %f}\n", s.x, r_sigb[e8].g, r_sigb[e8].b) end
    end
  end
end

task Step1b_sigx_y2(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_sigb : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            pry_mesh : region(ispace(int8d), mesh),
            pry_sig : region(ispace(int8d), grid),
            pry_sig2 : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_sig2, r_mesh, vxmesh, vymesh, vzmesh, pry_mesh, pry_sig, pry_sig2),
  reads writes(r_sigb)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 0  --change when copy
  var Dim2 : int32 = 1

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  var sC : double
 
  var bc : int32[6]
  var IR : int32
  var JR : int32
  var KR : int32
  
  var e3 : int8d
  var e7 : int8d
  var e8 : int8d
  var eR3 : int8d
  var eR7 : int8d
  var eR8 : int8d
 
  var gsig : double
  var bsig : double
  var gsig2 : double
  var bsig2 : double
  var swap : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    sC = r_mesh[e3].dy -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
      eR8 = {IR, JR, KR, Dim, Dim2, v.x, v.y, v.z}
  
      swap = 1.0
      if (vymesh[v.y].v < 0 and Dim == 1) then
        if s.y == r_sigb.bounds.hi.y then
          gsig = pry_sig[eR7].g
          bsig = pry_sig[eR7].b
          gsig2 = pry_sig2[eR8].g
          bsig2 = pry_sig2[eR8].b
          sC = pry_mesh[eR3].dy
        else
          gsig = r_sig[eR7].g
          bsig = r_sig[eR7].b
          gsig2 = r_sig2[eR8].g
          bsig2 = r_sig2[eR8].b
          sC = r_mesh[eR3].dy
        end
        swap = -1
      else
        gsig = r_sig[e7].g
        bsig = r_sig[e7].b
        gsig2 = r_sig2[e8].g
        bsig2 = r_sig2[e8].b
        sC = r_mesh[e3].dy
      end
 
      r_sigb[e8].g = gsig + swap*(sC/2.0)*gsig2
      r_sigb[e8].b = bsig + swap*(sC/2.0)*bsig2
      
      --if v.y == 128 then c.printf("rsigb[%d] = {%f, %f}\n", s.y, r_sigb[e8].g, r_sigb[e8].b) end
    end
  end
end

task Step1b_sigy_y2(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_sigb : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            pry_mesh : region(ispace(int8d), mesh),
            pry_sig : region(ispace(int8d), grid),
            pry_sig2 : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_sig2, r_mesh, vxmesh, vymesh, vzmesh, pry_mesh, pry_sig, pry_sig2),
  reads writes(r_sigb)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 1  --change when copy
  var Dim2 : int32 = 1

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  var sC : double
 
  var bc : int32[6]
  var IR : int32
  var JR : int32
  var KR : int32
  
  var e3 : int8d
  var e7 : int8d
  var e8 : int8d
  var eR3 : int8d
  var eR7 : int8d
  var eR8 : int8d
 
  var gsig : double
  var bsig : double
  var gsig2 : double
  var bsig2 : double
  var swap : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    sC = r_mesh[e3].dy -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
      eR8 = {IR, JR, KR, Dim, Dim2, v.x, v.y, v.z}
  
      swap = 1.0
      if (vymesh[v.y].v < 0 and Dim == 1) then
        if s.y == r_sigb.bounds.hi.y then
          gsig = pry_sig[eR7].g
          bsig = pry_sig[eR7].b
          gsig2 = pry_sig2[eR8].g
          bsig2 = pry_sig2[eR8].b
          sC = pry_mesh[eR3].dy
        else
          gsig = r_sig[eR7].g
          bsig = r_sig[eR7].b
          gsig2 = r_sig2[eR8].g
          bsig2 = r_sig2[eR8].b
          sC = r_mesh[eR3].dy
        end
        swap = -1
      else
        gsig = r_sig[e7].g
        bsig = r_sig[e7].b
        gsig2 = r_sig2[e8].g
        bsig2 = r_sig2[e8].b
        sC = r_mesh[e3].dy
      end
 
      r_sigb[e8].g = gsig + swap*(sC/2.0)*gsig2
      r_sigb[e8].b = bsig + swap*(sC/2.0)*bsig2
      
      --if v.y == 32 and s.x == 32 then c.printf("rsigb[%d] = {%f, %f}\n", s.y, r_sigb[e8].g, r_sigb[e8].b) end
    end
  end
end

task Step1b_sigz_y2(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_sigb : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            pry_mesh : region(ispace(int8d), mesh),
            pry_sig : region(ispace(int8d), grid),
            pry_sig2 : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_sig2, r_mesh, vxmesh, vymesh, vzmesh, pry_mesh, pry_sig, pry_sig2),
  reads writes(r_sigb)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 2  --change when copy
  var Dim2 : int32 = 1

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  var sC : double
 
  var bc : int32[6]
  var IR : int32
  var JR : int32
  var KR : int32
  
  var e3 : int8d
  var e7 : int8d
  var e8 : int8d
  var eR3 : int8d
  var eR7 : int8d
  var eR8 : int8d
 
  var gsig : double
  var bsig : double
  var gsig2 : double
  var bsig2 : double
  var swap : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    sC = r_mesh[e3].dy -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
      eR8 = {IR, JR, KR, Dim, Dim2, v.x, v.y, v.z}
  
      swap = 1.0
      if (vymesh[v.y].v < 0 and Dim == 1) then
        if s.y == r_sigb.bounds.hi.y then
          gsig = pry_sig[eR7].g
          bsig = pry_sig[eR7].b
          gsig2 = pry_sig2[eR8].g
          bsig2 = pry_sig2[eR8].b
          sC = pry_mesh[eR3].dy
        else
          gsig = r_sig[eR7].g
          bsig = r_sig[eR7].b
          gsig2 = r_sig2[eR8].g
          bsig2 = r_sig2[eR8].b
          sC = r_mesh[eR3].dy
        end
        swap = -1
      else
        gsig = r_sig[e7].g
        bsig = r_sig[e7].b
        gsig2 = r_sig2[e8].g
        bsig2 = r_sig2[e8].b
        sC = r_mesh[e3].dy
      end
 
      r_sigb[e8].g = gsig + swap*(sC/2.0)*gsig2
      r_sigb[e8].b = bsig + swap*(sC/2.0)*bsig2
      
      --if v.y == 128 then c.printf("rsigb[%d] = {%f, %f}\n", s.x, r_sigb[e8].g, r_sigb[e8].b) end
    end
  end
end

task Step1b_sigx_z2(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_sigb : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            prz_mesh : region(ispace(int8d), mesh),
            prz_sig : region(ispace(int8d), grid),
            prz_sig2 : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_sig2, r_mesh, vxmesh, vymesh, vzmesh, prz_mesh, prz_sig, prz_sig2),
  reads writes(r_sigb)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 0  --change when copy
  var Dim2 : int32 = 2

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  var sC : double
 
  var bc : int32[6]
  var IR : int32
  var JR : int32
  var KR : int32
  
  var e3 : int8d
  var e7 : int8d
  var e8 : int8d
  var eR3 : int8d
  var eR7 : int8d
  var eR8 : int8d
 
  var gsig : double
  var bsig : double
  var gsig2 : double
  var bsig2 : double
  var swap : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    sC = r_mesh[e3].dz -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
      eR8 = {IR, JR, KR, Dim, Dim2, v.x, v.y, v.z}
  
      swap = 1.0
      if (vzmesh[v.z].v < 0 and Dim == 2) then
        if s.z == r_sigb.bounds.hi.z then
          gsig = prz_sig[eR7].g
          bsig = prz_sig[eR7].b
          gsig2 = prz_sig2[eR8].g
          bsig2 = prz_sig2[eR8].b
          sC = prz_mesh[eR3].dz
        else
          gsig = r_sig[eR7].g
          bsig = r_sig[eR7].b
          gsig2 = r_sig2[eR8].g
          bsig2 = r_sig2[eR8].b
          sC = r_mesh[eR3].dz
        end
        swap = -1
      else
        gsig = r_sig[e7].g
        bsig = r_sig[e7].b
        gsig2 = r_sig2[e8].g
        bsig2 = r_sig2[e8].b
        sC = r_mesh[e3].dz
      end
 
      r_sigb[e8].g = gsig + swap*(sC/2.0)*gsig2
      r_sigb[e8].b = bsig + swap*(sC/2.0)*bsig2
      
      --if v.z == 128 then c.printf("rsigb[%d] = {%f, %f}\n", s.z, r_sigb[e8].g, r_sigb[e8].b) end
    end
  end
end

task Step1b_sigy_z2(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_sigb : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            prz_mesh : region(ispace(int8d), mesh),
            prz_sig : region(ispace(int8d), grid),
            prz_sig2 : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_sig2, r_mesh, vxmesh, vymesh, vzmesh, prz_mesh, prz_sig, prz_sig2),
  reads writes(r_sigb)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 1  --change when copy
  var Dim2 : int32 = 2

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  var sC : double
 
  var bc : int32[6]
  var IR : int32
  var JR : int32
  var KR : int32
  
  var e3 : int8d
  var e7 : int8d
  var e8 : int8d
  var eR3 : int8d
  var eR7 : int8d
  var eR8 : int8d
 
  var gsig : double
  var bsig : double
  var gsig2 : double
  var bsig2 : double
  var swap : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    sC = r_mesh[e3].dz -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
      eR8 = {IR, JR, KR, Dim, Dim2, v.x, v.y, v.z}
  
      swap = 1.0
      if (vzmesh[v.z].v < 0 and Dim == 2) then
        if s.z == r_sigb.bounds.hi.z then
          gsig = prz_sig[eR7].g
          bsig = prz_sig[eR7].b
          gsig2 = prz_sig2[eR8].g
          bsig2 = prz_sig2[eR8].b
          sC = prz_mesh[eR3].dz
        else
          gsig = r_sig[eR7].g
          bsig = r_sig[eR7].b
          gsig2 = r_sig2[eR8].g
          bsig2 = r_sig2[eR8].b
          sC = r_mesh[eR3].dz
        end
        swap = -1
      else
        gsig = r_sig[e7].g
        bsig = r_sig[e7].b
        gsig2 = r_sig2[e8].g
        bsig2 = r_sig2[e8].b
        sC = r_mesh[e3].dz
      end
 
      r_sigb[e8].g = gsig + swap*(sC/2.0)*gsig2
      r_sigb[e8].b = bsig + swap*(sC/2.0)*bsig2
      
      --if v.z == 128 then c.printf("rsigb[%d] = {%f, %f}\n", s.z, r_sigb[e8].g, r_sigb[e8].b) end
    end
  end
end

task Step1b_sigz_z2(r_sig : region(ispace(int8d), grid),
            r_sig2 : region(ispace(int8d), grid),
            r_sigb : region(ispace(int8d), grid),
            r_mesh : region(ispace(int8d), mesh),
            prz_mesh : region(ispace(int8d), mesh),
            prz_sig : region(ispace(int8d), grid),
            prz_sig2 : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            BCs : int32[3], N : int32[3])
where
  reads(r_sig, r_sig2, r_mesh, vxmesh, vymesh, vzmesh, prz_mesh, prz_sig, prz_sig2),
  reads writes(r_sigb)
do
  -- Dim  is vector component that is being interpolated.
  -- Dim2 is direction of interpolation.
  var Dim : int32 = 2  --change when copy
  var Dim2 : int32 = 2

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  var sC : double
 
  var bc : int32[6]
  var IR : int32
  var JR : int32
  var KR : int32
  
  var e3 : int8d
  var e7 : int8d
  var e8 : int8d
  var eR3 : int8d
  var eR7 : int8d
  var eR8 : int8d
 
  var gsig : double
  var bsig : double
  var gsig2 : double
  var bsig2 : double
  var swap : double

  for s in s3 do
    
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    sC = r_mesh[e3].dz -- change when copy

    -- Gather Left and Right Indices
    bc = BC(s.x, s.y, s.z, Dim2, BCs, N)
    IR = bc[3]
    JR = bc[4]
    KR = bc[5]
    eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

    for v in v3 do

      e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
      e8 = {s.x, s.y, s.z, Dim, Dim2, v.x, v.y, v.z}
      eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}
      eR8 = {IR, JR, KR, Dim, Dim2, v.x, v.y, v.z}
  
      swap = 1.0
      if (vzmesh[v.z].v < 0 and Dim == 2) then
        if s.z == r_sigb.bounds.hi.z then
          gsig = prz_sig[eR7].g
          bsig = prz_sig[eR7].b
          gsig2 = prz_sig2[eR8].g
          bsig2 = prz_sig2[eR8].b
          sC = prz_mesh[eR3].dz
        else
          gsig = r_sig[eR7].g
          bsig = r_sig[eR7].b
          gsig2 = r_sig2[eR8].g
          bsig2 = r_sig2[eR8].b
          sC = r_mesh[eR3].dz
        end
        swap = -1
      else
        gsig = r_sig[e7].g
        bsig = r_sig[e7].b
        gsig2 = r_sig2[e8].g
        bsig2 = r_sig2[e8].b
        sC = r_mesh[e3].dz
      end
 
      r_sigb[e8].g = gsig + swap*(sC/2.0)*gsig2
      r_sigb[e8].b = bsig + swap*(sC/2.0)*bsig2
      
      --if v.z == 128 then c.printf("rsigb[%d] = {%f, %f}\n", s.z, r_sigb[e8].g, r_sigb[e8].b) end
    end
  end
end

task Step1b_b(r_sig : region(ispace(int8d), grid),
              r_mesh : region(ispace(int8d), mesh),
              r_gridbarp : region(ispace(int8d), grid),
              r_gridbarpb : region(ispace(int8d), grid),
              prx_gridbarp : region(ispace(int8d), grid),
              pry_gridbarp : region(ispace(int8d), grid),
              prz_gridbarp : region(ispace(int8d), grid),
              prx_sig : region(ispace(int8d), grid),
              pry_sig : region(ispace(int8d), grid),
              prz_sig : region(ispace(int8d), grid),
	      prx_mesh : region(ispace(int8d), mesh),
	      pry_mesh : region(ispace(int8d), mesh),
	      prz_mesh : region(ispace(int8d), mesh),
              vxmesh : region(ispace(int1d), vmesh),
              vymesh : region(ispace(int1d), vmesh),
              vzmesh : region(ispace(int1d), vmesh),
              BCs : int32[3], N : int32[3], effD : int32)
where
  reads writes(r_gridbarpb),
  reads(r_sig, r_gridbarp, r_mesh.{dx,dy,dz}, vxmesh, vymesh, vzmesh),
  reads(prx_gridbarp, pry_gridbarp, prz_gridbarp, prx_sig, pry_sig, prz_sig, prx_mesh.dx, pry_mesh.dy, prz_mesh.dz)
do

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  var sC : double[3]
  
  var bc : int32[6]
  var IR : int32
  var JR : int32
  var KR : int32

  var e3 : int8d
  var e6 : int8d
  var e7 : int8d
  var eR3 : int8d
  var eR6 : int8d
  var eR7 : int8d

  var swap : double
  var gb : double
  var bb : double
  var gsig : double
  var bsig : double

  for s in s3 do

    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}

    for Dim = 0, effD do

      -- Gathering Right Indices
      bc = BC(s.x, s.y, s.z, Dim, BCs, N)
      IR = bc[3]
      JR = bc[4]
      KR = bc[5]
      eR3 = {IR, JR, KR, 0, 0, 0, 0, 0}

      for v in v3 do      

        e6 = {s.x, s.y, s.z, 0, 0, v.x, v.y, v.z}
        e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}

        -- Gathering Right Indices 
        eR6 = {IR, JR, KR, 0, 0, v.x, v.y, v.z}
        eR7 = {IR, JR, KR, Dim, 0, v.x, v.y, v.z}

        -- Dot Product is just a single product when using rectangular mesh
        swap = 1.0
        gb = r_gridbarp[e6].g
        bb = r_gridbarp[e6].b
        gsig = r_sig[e7].g
        bsig = r_sig[e7].b
      
        sC[0] = r_mesh[e3].dx
        sC[1] = r_mesh[e3].dy
        sC[2] = r_mesh[e3].dz

        if (vxmesh[v.x].v < 0 and Dim == 0) then
          swap = -1
          if s.x == r_mesh.bounds.hi.x then
            gsig = prx_sig[eR7].g
            bsig = prx_sig[eR7].b
            gb = prx_gridbarp[eR6].g
            bb = prx_gridbarp[eR6].b
            sC[Dim] = prx_mesh[eR3].dx
           else
            gsig = r_sig[eR7].g
            bsig = r_sig[eR7].b
            gb = r_gridbarp[eR6].g
            bb = r_gridbarp[eR6].b
            sC[Dim] = r_mesh[eR3].dx
           end
        elseif (vymesh[v.y].v < 0 and Dim == 1) then
          swap = -1
          if s.y == r_mesh.bounds.hi.y then
            gsig = pry_sig[eR7].g
            bsig = pry_sig[eR7].b
            gb = pry_gridbarp[eR6].g
            bb = pry_gridbarp[eR6].b
            sC[Dim] = pry_mesh[eR3].dy
          else
            gsig = r_sig[eR7].g
            bsig = r_sig[eR7].b
            gb = r_gridbarp[eR6].g
            bb = r_gridbarp[eR6].b
            sC[Dim] = r_mesh[eR3].dy
          end
        elseif (vzmesh[v.z].v < 0 and Dim == 2) then
          swap = -1
          if s.z == r_mesh.bounds.hi.z then
            gsig = prz_sig[eR7].g
            bsig = prz_sig[eR7].b
            gb = prz_gridbarp[eR6].g
            bb = prz_gridbarp[eR6].b
            sC[Dim] = prz_mesh[eR3].dz
          else
            gsig = r_sig[eR7].g
            bsig = r_sig[eR7].b
            gb = r_gridbarp[eR6].g
            bb = r_gridbarp[eR6].b
            sC[Dim] = r_mesh[eR3].dz
          end
        end
     
        r_gridbarpb[e7].g = gb + swap*sC[Dim]/2.0*gsig 
        r_gridbarpb[e7].b = bb + swap*sC[Dim]/2.0*bsig

        -- NAN checker
        if (isnan(r_gridbarpb[e7].g) == 1 or isnan(r_gridbarpb[e7].b) == 1) then
 
          c.printf("Step 1b: r_gridbarp.g = %f, r_gridbarp.b = %f, r_sig.g = %f, r_sig.b = %f\n", gb, bb, gsig, bsig)

          regentlib.assert(not [bool](isnan(r_gridbarpb[e7].g)), "Step 1b\n")
          regentlib.assert(not [bool](isnan(r_gridbarpb[e7].b)), "Step 1b\n")
    
        end

      end
    end 
  end
  --c.printf("Step1b rhotest[32] = %f\n", rhotest)
end

-- Step 1c: Compute phibar at interface by interpolating w/ phisigma2, x-Xi*dt/2
task Step1c(r_gridbarpb : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),        
            vymesh : region(ispace(int1d), vmesh),        
            vzmesh : region(ispace(int1d), vmesh),        
            r_sigb : region(ispace(int8d), grid),
            dt : double, BCs : int32[3],  N : int32[3], effD : int32)
where
  reads(vxmesh, vymesh, vzmesh, r_sigb), 
  reads writes(r_gridbarpb)
do     
  -- Compute gbar/bbar @ t=n+1/2  with interface sigma
  var Xi : double[3]
  
  var slo : int3d = {r_gridbarpb.bounds.lo.x, r_gridbarpb.bounds.lo.y, r_gridbarpb.bounds.lo.z}
  var shi : int3d = {r_gridbarpb.bounds.hi.x, r_gridbarpb.bounds.hi.y, r_gridbarpb.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  var e7 : int8d
  var e8 : int8d

  for v in v3 do
    Xi[0] = vxmesh[v.x].v
    Xi[1] = vymesh[v.y].v
    Xi[2] = vzmesh[v.z].v

    for s in s3 do
      for Dim = 0, effD do
        for Dim2 = 0, effD do

          -- Gather Indices
          e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
          e8 = {s.x, s.y, s.z, Dim2, Dim, v.x, v.y, v.z}

          r_gridbarpb[e7].g = r_gridbarpb[e7].g - dt/2.0*Xi[Dim2]*r_sigb[e8].g
          r_gridbarpb[e7].b = r_gridbarpb[e7].b - dt/2.0*Xi[Dim2]*r_sigb[e8].b
    
          regentlib.assert(not [bool](isnan(r_gridbarpb[e7].g)), "Step 1c\n")
          regentlib.assert(not [bool](isnan(r_gridbarpb[e7].b)), "Step 1c\n")
        
        end 
      end 
    end
  end 
end

--Step 2: Microflux
--Step 2a: Compute W at interface.
task Step2a(r_gridbarpb : region(ispace(int8d), grid),
            vxmesh : region(ispace(int1d), vmesh),
            vymesh : region(ispace(int1d), vmesh),
            vzmesh : region(ispace(int1d), vmesh),
            r_Wb   : region(ispace(int8d), W),
            dt : double, effD : int32)
where
  reads(r_gridbarpb, vxmesh, vymesh, vzmesh),
  reads writes(r_Wb)
do    
  var slo : int3d = {r_Wb.bounds.lo.x, r_Wb.bounds.lo.y, r_Wb.bounds.lo.z}
  var shi : int3d = {r_Wb.bounds.hi.x, r_Wb.bounds.hi.y, r_Wb.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  var e4 : int8d
  var e7 : int8d

  -- First do density at boundary, density is needed for others.
  for s in s3 do
    for Dim = 0, effD do
      e4 = {s.x, s.y, s.z, Dim, 0, 0, 0, 0}
      r_Wb[e4].rho = 0

      for v in v3 do
        e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
        r_Wb[e4].rho = r_Wb[e4].rho + vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w*r_gridbarpb[e7].g
      end
    end
  end
      
  -- Then do momentum and energy
  -- Initialize, then iterate over contributions in velocity space
  var U : double[3]
  for s in s3 do
    for Dim = 0, effD do

      e4 = {s.x, s.y, s.z, Dim, 0, 0, 0, 0}  
      for d = 0, effD do
        r_Wb[e4].rhov[d] = dt/2.0*r_Wb[e4].rho*0 -- TODO: In future, replace 0 with acceleration field 
      end
      r_Wb[e4].rhoE = dt/2.*r_Wb[e4].rho*0 -- TODO: In future replace 0 with u.dot(a), vel dot acc

      for v in v3 do

        U[0] = vxmesh[v.x].v
        U[1] = vymesh[v.y].v
        U[2] = vzmesh[v.z].v
 

        e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}    

        for d = 0, effD do
          r_Wb[e4].rhov[d] = r_Wb[e4].rhov[d] + vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w*U[d]*r_gridbarpb[e7].g 
        end

        r_Wb[e4].rhoE = r_Wb[e4].rhoE + vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w*r_gridbarpb[e7].b
      end
    end
  end

  -- NAN checker
  for e in r_Wb do
    
    --if e.x == 32 then c.printf("r_Wb[%d, %d, Dim = %d] = {%f, {%f, %f}, %f}\n", e.x, e.y, e.w, r_Wb[e].rho, r_Wb[e].rhov[0], r_Wb[e].rhov[1], r_Wb[e].rhoE) end
    regentlib.assert(not [bool](isnan(r_Wb[e].rho)), "Step 2a rho\n")
    regentlib.assert(not [bool](isnan(r_Wb[e].rhov[0])), "Step 2a rhov0\n")
    regentlib.assert(not [bool](isnan(r_Wb[e].rhov[1])), "Step 2a rhov1\n")
    regentlib.assert(not [bool](isnan(r_Wb[e].rhov[2])), "Step 2a rhov2\n")
    regentlib.assert(not [bool](isnan(r_Wb[e].rhoE)), "Step 2a rhoE\n")
    
  end

end

-- Step 2b: compute original phi at interface using gbar, W at interface
-- Memory Recycling: phibar @ interface is used to store phi @ interface.
task Step2b(r_gridbarpb : region(ispace(int8d), grid), 
            r_Wb      : region(ispace(int8d), W),
            vxmesh    : region(ispace(int1d), vmesh),
            vymesh    : region(ispace(int1d), vmesh),
            vzmesh    : region(ispace(int1d), vmesh),
            dt : double, R : double, K : double, Cv : double, g : double,
            w : double, ur : double, Tr : double, Pr : double, effD : int32)
where
  reads writes(r_gridbarpb),
  reads(r_Wb, vxmesh, vymesh, vzmesh)
do
  var slo : int3d = {r_Wb.bounds.lo.x, r_Wb.bounds.lo.y, r_Wb.bounds.lo.z}
  var shi : int3d = {r_Wb.bounds.hi.x, r_Wb.bounds.hi.y, r_Wb.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)
  
  var tg : double
  var tb : double
  var u : double
  var T : double
  var Xi : double[3]
  var c2 : double
  var g_eq : double
  var b_eq : double

  var e3 : int8d
  var e4 : int8d
  var e7 : int8d

  for s in s3 do

    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}

    for Dim = 0, effD do

      e4 = {s.x, s.y, s.z, Dim, 0, 0, 0, 0}

      u = 0 
      for d = 0, effD do
        u = u + r_Wb[e4].rhov[d]/r_Wb[e4].rho*r_Wb[e4].rhov[d]/r_Wb[e4].rho
      end
      u = sqrt(u)
      
      T = Temperature(r_Wb[e4].rhoE/r_Wb[e4].rho, u, g, R)
      if T < 0 then
        c.printf("T < 0, r_Wb[{%d, %d, %d, %d}].rhoE = %f, r_Wb[e4].rho = %f, u = %f, g = %f, R = %f\n", e4.x, e4.y, e4.z, e4.w, r_Wb[e4].rhoE, r_Wb[e4].rho, u, g, R)
        regentlib.assert(T >= 0, "Negative Temperature\n")
      end

      tg = visc(T, ur, Tr, w)/r_Wb[e4].rho/R/T
      tb = tg/Pr
      --c.printf("tg = %f, tb = %f\n", tg, tb)

      for v in v3 do

        e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
        c2 = 0

        Xi[0] = vxmesh[v.x].v
        Xi[1] = vymesh[v.y].v
        Xi[2] = vzmesh[v.z].v

        for d = 0, effD do
          c2 =  c2 + (Xi[d] - r_Wb[e4].rhov[d]/r_Wb[e4].rho)*(Xi[d] - r_Wb[e4].rhov[d]/r_Wb[e4].rho)
        end

        g_eq = geq(c2, r_Wb[e4].rho, T, R, effD)
        b_eq = g_eq*(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2] + (3.0-effD+K)*R*T)/2.0

        if (isnan(g_eq) == 1 or isnan(b_eq) == 1) then

          c.printf("c2 = %f, r_Wb[e4].rho = %f, T = %f, R = %f, effD = %d\n", c2, r_Wb[e4].rho, T, R, effD) 
      
          regentlib.assert(not  [bool](isnan(g_eq)), "Step 2b\n")
          regentlib.assert(not  [bool](isnan(b_eq)), "Step 2b\n")
    
        end

        -- this is actually the original distribution function, recycling memory from gbar
        r_gridbarpb[e7].g = 2*tg/(2*tg + dt/2.)*r_gridbarpb[e7].g + dt/(4*tg + dt)*g_eq + dt*tg/(4*tg + dt)*0 -- TODO replace this last *0 with source term 
        r_gridbarpb[e7].b = 2*tb/(2*tb + dt/2.)*r_gridbarpb[e7].b + dt/(4*tb + dt)*b_eq + dt*tb/(4*tb + dt)*0 -- TODO replace this last *0 with source term

        if (isnan(r_gridbarpb[e7].g) == 1 or isnan(r_gridbarpb[e7].b) == 1) then

          c.printf("gbar = %.10f, bbar = %.10f, g_eq = %.10f, tg = %.10f, tb = %.10f\n", r_gridbarpb[e7].g, r_gridbarpb[e7].b, g_eq, tg, tb)
      
          regentlib.assert(not [bool](isnan(r_gridbarpb[e7].g)), "Step 2b\n")
          regentlib.assert(not [bool](isnan(r_gridbarpb[e7].b)), "Step 2b\n")
    
        end

      end
    end 
  end
end

-- Step 2c: Compute Microflux F at interface at half timestep using W/phi at interface.
task Step2c(r_gridbarpb : region(ispace(int8d), grid),
            r_F       : region(ispace(int8d), grid),
            r_mesh    : region(ispace(int8d), mesh),
            vxmesh    : region(ispace(int1d), vmesh),
            vymesh    : region(ispace(int1d), vmesh),
            vzmesh    : region(ispace(int1d), vmesh),
            plx_gridbarpb : region(ispace(int8d), grid),
            ply_gridbarpb : region(ispace(int8d), grid),
            plz_gridbarpb : region(ispace(int8d), grid),
            BCs : int32[3], R : double, K : double, Cv : double, g : double,
            w : double, Tr : double, Pr : double, effD : int32, N : int32[3])
where
  reads(r_gridbarpb, vxmesh, vymesh, vzmesh, r_mesh, plx_gridbarpb, ply_gridbarpb, plz_gridbarpb),
  reads writes(r_F)
do    
  var A : double[3]
  var Xi : double[3]

  fill(r_F.g, 0)
  fill(r_F.b, 0)

  var slo : int3d = {r_mesh.bounds.lo.x, r_mesh.bounds.lo.y, r_mesh.bounds.lo.z}
  var shi : int3d = {r_mesh.bounds.hi.x, r_mesh.bounds.hi.y, r_mesh.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  var e3 : int8d

  var bc : int32[6]
  var IL : int32
  var JL : int32
  var KL : int32
  var IR : int32
  var JK : int32
  var KR : int32

  var right : double = 1.0 
  var left : double = 1.0

  var e7 : int8d 
  var e6 : int8d 
  var eL7 : int8d

  var gL : double
  var bL : double

  for s in s3 do
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}

    A[0] = r_mesh[e3].dy*r_mesh[e3].dz
    A[1] = r_mesh[e3].dx*r_mesh[e3].dz
    A[2] = r_mesh[e3].dx*r_mesh[e3].dy


    for Dim = 0, effD do
      bc = BC(s.x, s.y, s.z, Dim, BCs, N)
      IL = bc[0]
      JL = bc[1]
      KL = bc[2]
      IR = bc[3]
      JK = bc[4]
      KR = bc[5]

      -- Periodic Boundary Conditions
      right = 1.0 
      left = 1.0

      -- Dirichlet Boundary Conditions
      if (Dim == 0 and BCs[0] == 1 and (s.x == 0 or s.x == N[0] - 1)) then
        left = 0
        right = 0
      end 
      if (Dim == 1 and BCs[1] == 1 and (s.y == 0 or s.y == N[1] - 1)) then
        left = 0
        right = 0
      end
      if (Dim == 2 and BCs[2] == 1 and (s.z == 0 or s.z == N[2] - 1)) then
        left = 0
        right = 0
      end

      -- TODO: Neumann Boundary Conditions
      if (Dim == 0 and BCs[0] == 1) then
        if s.x == 0 then 
          left = 0
        elseif s.x == N[0] - 1 then 
          right = 0
        end
      end 
      if (Dim == 1 and BCs[1] == 1) then
        if s.y == 0 then 
          left = 0
        elseif s.y == N[1] - 1 then 
          right = 0
        end
      end
      if (Dim == 2 and BCs[2] == 1) then
        if s.z == 0 then 
          left = 0
        elseif s.z == N[2] - 1 then 
          right = 0
        end
      end

      for v in v3 do 
        -- Gather Left Indices
        e7 = {s.x, s.y, s.z, Dim, 0, v.x, v.y, v.z}
        e6 = {s.x, s.y, s.z, 0, 0, v.x, v.y, v.z}
        eL7 = {IL, JL, KL, Dim, 0, v.x, v.y, v.z}

        Xi[0] = vxmesh[v.x].v
        Xi[1] = vymesh[v.y].v
        Xi[2] = vzmesh[v.z].v
  
        if Dim == 0 then
          if s.x == s3.bounds.lo.x then
            gL = plx_gridbarpb[eL7].g
            bL = plx_gridbarpb[eL7].b
          else
            gL = r_gridbarpb[eL7].g
            bL = r_gridbarpb[eL7].b
          end
        end

        if Dim == 1 then
          if s.y == s3.bounds.lo.y then
            gL = ply_gridbarpb[eL7].g
            bL = ply_gridbarpb[eL7].b
          else
            gL = r_gridbarpb[eL7].g
            bL = r_gridbarpb[eL7].b
          end
        end

        if Dim == 2 then
          if s.z == s3.bounds.lo.z then
            gL = plz_gridbarpb[eL7].g
            bL = plz_gridbarpb[eL7].b
          else
            gL = r_gridbarpb[eL7].g
            bL = r_gridbarpb[eL7].b
          end
        end

        if (v.x == 32 and v.y == 32 and s.x == 32) then
          --c.printf("pre r_F[{%d, %d, Dim = %d}}]= {%f, %f}, right/left = {%f, %f}, gL/gR = {%f, %f}, bL/bR = {%f, %f}, A[%d] = %f\n", s.x, s.y, Dim, r_F[e6].g, r_F[e6].b, right, left, gL, r_gridbarpb[e7].g, bL, r_gridbarpb[e7].b, Dim, A[Dim])
        end
           
        r_F[e6].g = r_F[e6].g + Xi[Dim]*A[Dim]*(right*r_gridbarpb[e7].g - left*gL)
        r_F[e6].b = r_F[e6].b + Xi[Dim]*A[Dim]*(right*r_gridbarpb[e7].b - left*bL)

        if (v.x == 32 and v.y == 32 and s.x == 32) then
          --c.printf("post r_F[{%d, %d, Dim = %d}}]= {%f, %f}, right/left = {%f, %f}, gL/gR = {%f, %f}, bL/bR = {%f, %f}, A[%d] = %f\n", s.x, s.y, Dim, r_F[e6].g, r_F[e6].b, right, left, gL, r_gridbarpb[e7].g, bL, r_gridbarpb[e7].b, Dim, A[Dim])
        end
        regentlib.assert(not [bool](isnan(r_F[e6].g)), "Step 2c\n")
        regentlib.assert(not [bool](isnan(r_F[e6].b)), "Step 2c\n")
      end 
    end
  end 
end

--Step 3: Source Terms
task Step3()
  --TODO
end

--Step 4: Update Conservative Variables W at cell center at next timestep
--Step 5: Update Phi at cell center at next time step
task Step4and5(r_grid : region(ispace(int8d), grid),
               r_W    : region(ispace(int8d), W),
               r_mesh : region(ispace(int8d), mesh),
               r_F    : region(ispace(int8d), grid),
               vxmesh : region(ispace(int1d), vmesh),
               vymesh : region(ispace(int1d), vmesh),
               vzmesh : region(ispace(int1d), vmesh),
               dt : double, BCs : int32[3], R : double, K : double, Cv : double, N : int32[3],
               g : double, w : double, ur : double, Tr : double, Pr : double, effD : int32)
where
  reads(vxmesh, vymesh, vzmesh, r_mesh, r_F),
  reads writes(r_W, r_grid)
do
  var V : double      -- Volume of Cell
  var Xi : double[3]  -- Discrete Velocity 
  var uo : double     -- Old Flow Velocity
  var To : double     -- Old Temp
  var tgo : double    -- Old visc
  var tbo : double    -- Old visc
  var u : double      -- Flow Velocity
  var T : double      -- Temperature 
  var tg : double     -- visc
  var tb : double     -- visc 
  var c2 : double     -- peculiar velocity sq
  var g_eq : double   -- Equils
  var b_eq : double
  var g_eqo : double
  var b_eqo : double

  var slo : int3d = {r_F.bounds.lo.x, r_F.bounds.lo.y, r_F.bounds.lo.z}
  var shi : int3d = {r_F.bounds.hi.x, r_F.bounds.hi.y, r_F.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  
  var e3 : int8d 
  var e6 : int8d 
  var i : int32
  var j : int32
  var k : int32

  for s in s3 do

    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    
    V = r_mesh[e3].dx*r_mesh[e3].dy*r_mesh[e3].dz

    -- Compute old flow velocity
    uo = 0 -- old flow velocity
    for d = 0, effD do
      uo += r_W[e3].rhov[d]/r_W[e3].rho*r_W[e3].rhov[d]/r_W[e3].rho
    end
    uo = sqrt(uo)
    regentlib.assert(bool(uo>=0), "uo")

    -- Compute old temperature
    To = Temperature(r_W[e3].rhoE/r_W[e3].rho, uo, g, R)
    regentlib.assert(bool(To>=0), "To")
  
    -- Compute old taus
    tgo = visc(To, ur, Tr, w)/r_W[e3].rho/R/To
    tbo = tgo/Pr 
  
    for v in v3 do 
      Xi[0] = vxmesh[v.x].v
      Xi[1] = vymesh[v.y].v
      Xi[2] = vzmesh[v.z].v
      e6 = {s.x, s.y, s.z, 0, 0, v.x, v.y, v.z}

      -- Compute old eq's
      c2 = 0
      for d = 0, effD do
        c2 += (Xi[d]-r_W[e3].rhov[d]/r_W[e3].rho)*(Xi[d]-r_W[e3].rhov[d]/r_W[e3].rho)
      end
      g_eqo = geq(c2, r_W[e3].rho, To, R, effD)
      b_eqo = g_eqo*(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2] + (3-effD+K)*R*To)/2

      -- First Update phi
      i = s.x
      j = s.y
      k = s.z

      -- Update First Step (terms involving oldW)
      if ((BCs[0] == 1 and i == 0) or (BCs[0] == 1 and i == N[0] - 1) or
          (BCs[1] == 1 and j == 0 and effD > 1) or (BCs[1] == 1 and j == N[1] - 1 and effD > 1) or
          (BCs[2] == 1 and k == 0 and effD > 2) or (BCs[2] == 1 and k == N[2] - 1 and effD > 2)) then
       
        r_grid[e6].g = r_grid[e6].g
        r_grid[e6].b = r_grid[e6].b

        --c.printf("Not updating : s = {%d, %d, %d}\n", s.x, s.y, s.z)
      else
        r_grid[e6].g = r_grid[e6].g + dt/2.0*(g_eqo-r_grid[e6].g)/tgo - dt/V*r_F[e6].g + dt*0 -- TODO replace 0 with source term
        r_grid[e6].b = r_grid[e6].b + dt/2.0*(b_eqo-r_grid[e6].b)/tbo - dt/V*r_F[e6].b + dt*0 -- TODO replace 0 with source term
      end

    end
  end


  for s in s3 do
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    
    for v in v3 do 
      Xi[0] = vxmesh[v.x].v
      Xi[1] = vymesh[v.y].v
      Xi[2] = vzmesh[v.z].v

      e6 = {s.x, s.y, s.z, 0, 0, v.x, v.y, v.z}
      -- Step 4: Update W at cell center
      r_W[e3].rho = r_W[e3].rho - dt*(r_F[e6].g/V - 0)*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w -- TODO replace 0 with source term.
      for d = 0, effD do
        r_W[e3].rhov[d] = r_W[e3].rhov[d] - dt*Xi[d]*(r_F[e6].g/V - 0)*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w -- TODO replace 0 with source term
      end
  
      r_W[e3].rhoE = r_W[e3].rhoE - dt*(r_F[e6].b/V - 0)*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w -- TODO replace 0 with source term      
    end
  end     

  -- Second Update Phi at cell center using new tau/W
  -- Compute flow velocity u

  for s in s3 do
    e3 = {s.x, s.y, s.z, 0, 0, 0, 0, 0}

    u = 0 
    for d = 0, effD do 
      u += r_W[e3].rhov[d]/r_W[e3].rho*r_W[e3].rhov[d]/r_W[e3].rho
    end
    u = sqrt(u)
    regentlib.assert(bool(u>=0), "u")

    -- Compute T
    T = Temperature(r_W[e3].rhoE/r_W[e3].rho, u, g, R)
    if T < 0 then
      c.printf("r_W[%d].rhoE = %f, r_W[e3].rho = %f, u = %f, g = %f, R = %f\n", e3.x, r_W[e3].rhoE, r_W[e3].rho, u, g, R)
    end
    regentlib.assert(bool(T>=0), "T")
  
    -- Compute new taus
    tg = visc(T, ur, Tr, w)/r_W[e3].rho/R/T 
    tb = tg/Pr 
    -- printf("taus are = {%f, %f}\n", tg, tb)

    for v in v3 do
      Xi[0] = vxmesh[v.x].v
      Xi[1] = vymesh[v.y].v
      Xi[2] = vzmesh[v.z].v

      -- Compute new eq's
      c2 = 0  -- reset c2 from before
      for d = 0, effD do
        c2 += (Xi[d]-r_W[e3].rhov[d]/r_W[e3].rho)*(Xi[d]-r_W[e3].rhov[d]/r_W[e3].rho)
      end
      g_eq = geq(c2, r_W[e3].rho, T, R, effD)
      b_eq = g_eq*(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2] + (3-effD+K)*R*T)/2


      -- Second Update phi
      i = s.x
      j = s.y
      k = s.z
      e6 = {s.x, s.y, s.z, 0, 0, v.x, v.y, v.z} 

      if ((BCs[0] == 1 and i == 0) or (BCs[0] == 1 and i == N[0] - 1) or
          (BCs[1] == 1 and j == 0 and effD > 1) or (BCs[1] == 1 and j == N[1] - 1 and effD > 1) or
          (BCs[2] == 1 and k == 0 and effD > 2) or (BCs[2] == 1 and k == N[2] - 1 and effD > 2)) then
                
        r_grid[e6].g = r_grid[e6].g
        r_grid[e6].b = r_grid[e6].b
        --c.printf("Not updating : s = {%d, %d, %d}\n", s.x, s.y, s.z)
      else
        r_grid[e6].g = (r_grid[e6].g + dt/2.0*g_eq/tg)/(1+dt/2.0/tg)
        r_grid[e6].b = (r_grid[e6].b + dt/2.0*b_eq/tb)/(1+dt/2.0/tb)

        --r_grid[e6].g = r_grid[e6].g/(1+dt/2.0/tg)
        --r_grid[e6].b = r_grid[e6].b/(1+dt/2.0/tb)
      end

      if isnan(r_grid[e6].g) == 1 then
        c.printf("Step4and5: g_eq = %f, tg = %f, tgo = %f, r_F[e6].g = %f\n", g_eq, tg, tgo, r_F[e6].g)
      end 
      if isnan(r_grid[e6].b) == 1 then
        c.printf("Step4and5: b_eq = %f, tb = %f, tbo = %f, r_F[e6].b = %f\n", b_eq, tb, tbo, r_F[e6].b)
      end 
      regentlib.assert(not [bool](isnan(r_grid[e6].g)), "Step4and5\n")
      regentlib.assert(not [bool](isnan(r_grid[e6].b)), "Step4and5\n")
    end
  end
  --c.printf("Step4and5 Complete\n")
end



task MaxwellianInitialization(r_grid  : region(ispace(int8d), grid),
         r_mesh : region(ispace(int8d), mesh),
         r_W    : region(ispace(int8d), W),
         vxmesh : region(ispace(int1d), vmesh),
         vymesh : region(ispace(int1d), vmesh),
         vzmesh : region(ispace(int1d), vmesh),
         testProblem: int32, R : double, K : double, Cv : double, 
         g : double, w : double, ur : double, Tr : double, Pr : double,
         N : int32[3], NV : int32[3], effD : int32)
where 
  reads writes(r_grid), 
  reads(r_mesh, r_W, vxmesh, vymesh, vzmesh)
do
  var T : double
  var u : double 

  var rhotest : double = 0
  var rhovxtest : double = 0
  var rhovytest : double = 0
  var Etest : double = 0

  var slo : int3d = {r_grid.bounds.lo.x, r_grid.bounds.lo.y, r_grid.bounds.lo.z}
  var shi : int3d = {r_grid.bounds.hi.x, r_grid.bounds.hi.y, r_grid.bounds.hi.z}
  var vlo : int3d = {vxmesh.bounds.lo, vymesh.bounds.lo, vzmesh.bounds.lo}
  var vhi : int3d = {vxmesh.bounds.hi, vymesh.bounds.hi, vzmesh.bounds.hi}
  var s3 = ispace(int3d, shi - slo + {1,1,1}, slo)
  var v3 = ispace(int3d, vhi - vlo + {1,1,1}, vlo)

  for s in s3 do
    var e3 : int8d = {s.x, s.y, s.z, 0, 0, 0, 0, 0}
    u = 0
    for dim = 0, effD do
      u += r_W[e3].rhov[dim]/r_W[e3].rho*r_W[e3].rhov[dim]/r_W[e3].rho
    end
    u = sqrt(u)
      
    T = Temperature(r_W[e3].rhoE/r_W[e3].rho, u, g, R)
    
    for v in v3 do
    
      var e6 : int8d = {s.x, s.y, s.z, 0, 0, v.x, v.y, v.z}

      var c2 : double = 0
      var Xi : double[3]  
      Xi[0] = vxmesh[v.x].v
      Xi[1] = vymesh[v.y].v
      Xi[2] = vzmesh[v.z].v

      for d = 0, effD do
        c2 += (Xi[d] - r_W[e3].rhov[d]/r_W[e3].rho)*(Xi[d] - r_W[e3].rhov[d]/r_W[e3].rho) --TODO could be bugged
      end

      r_grid[e6].g = geq(c2, r_W[e3].rho, T, R, effD)
      r_grid[e6].b = r_grid[e6].g*(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2] + (3.0-effD+K)*R*T)/2.0
  
      if s.x == r_grid.bounds.lo.x and s.y == r_grid.bounds.lo.y and s.z == 0 then
        rhotest += r_grid[e6].g*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w
        rhovxtest += r_grid[e6].g*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w*Xi[0]
        rhovytest += r_grid[e6].g*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w*Xi[1]
        Etest += r_grid[e6].b*vxmesh[v.x].w*vymesh[v.y].w*vzmesh[v.z].w
      end
    end
  end
  --c.printf("rhotest = %f, rhovxtest = %f, rhovytest = %f, Etest = %f\n", rhotest, rhovxtest, rhovytest, Etest)
end

task InitializeGrid(r_grid  : region(ispace(int8d), grid),
                r_mesh : region(ispace(int8d), mesh),
		r_W    : region(ispace(int8d), W),
		vxmesh : region(ispace(int1d), vmesh),
		vymesh : region(ispace(int1d), vmesh),
		vzmesh : region(ispace(int1d), vmesh),
                testProblem: int32, R : double, K : double, Cv : double, g : double, w : double, ur : double, Tr : double, Pr : double,
		N : int32[3], NV : int32[3], effD : int32)
where
  reads writes(r_grid),
  reads(r_mesh, r_W, vxmesh, vymesh, vzmesh)
do

  -- InitializeTestProblem
  if testProblem == 0 then
    -- TODO : User Defined
  elseif testProblem > 0 or testProblem == -1 then
    c.printf("Maxwellian Initialization\n")
    MaxwellianInitialization(r_grid, r_mesh, r_W, vxmesh, vymesh, vzmesh, testProblem, R, K, Cv, g, w, ur, Tr, Pr, N, NV, effD)
  end
end

terra dumpdouble(f : &c.FILE, val : double)
  var a : double[1]
  a[0] = val
  c.fwrite(&a, 8, 1, f)
end

terra dumpint32(f : &c.FILE, val : int32)
  var a : int32[1]
  a[0] = val
  c.fwrite(&a, 4, 1, f)
end

terra dumpbool(f : &c.FILE, val : bool)
  var a : int32[1]
  if val==true then
    a[0] = 1
  else 
    a[0] = 0
  end
  c.fwrite(&a, 4, 1, f)
end

task factorize1d(parallelism : int) : int3d
  var sizex = parallelism
  return int3d {sizex, 1, 1} 
end

task factorize2d(parallelism : int) : int3d
  var limit = [int](cmath.sqrt([double](parallelism)))
  var size_x : int32 = 1
  var size_y : int32 = parallelism
  for i = 1, limit + 1 do
    if parallelism % i == 0 then
      size_x, size_y =  i, parallelism / i
      if size_x > size_y then
        size_x, size_y = size_y, size_x
      end
    end
  end
  return int3d { size_x, size_y, 1 }
end

task factorize3d(parallelism : int) : int3d
  var limit = [int](cmath.pow([double](parallelism), 1/3))
  var size_x = 1
  var size_y = 1
  var size_z = parallelism
  for i = 1, limit + 1 do
    if parallelism % i == 0 then
      size_x, size_z = i, parallelism / i
      if size_x > size_z then
        size_x, size_z = size_z, size_x
      end
    end

    if parallelism % i*i == 0 then
      size_x, size_y, size_z = i, i, parallelism/i*i
    end
  end
  return int3d { size_x, size_y, size_z }
end

task factorize(parallelism: int, effD : int32)

  var f6 : int6d = {1, 1, 1, 1, 1, 1}
  if effD == 1 then
    var f3 = factorize1d(parallelism)
    f6.x, f6.y, f6.z = f3.x, f3.y, f3.z    
  elseif effD == 2 then
    var f3 = factorize2d(parallelism)
    f6.x, f6.y, f6.z = f3.x, f3.y, f3.z    
  elseif effD == 3 then
    var f3 = factorize3d(parallelism)
    f6.x, f6.y, f6.z = f3.x, f3.y, f3.z    
  end

  return f6
end

terra wait_for(x : int) return 1 end
task block_task(r_W : region(ispace(int8d), W))
where
  reads writes(r_W)
do
  return 1
end


task Dump(r_W : region(ispace(int8d), W), iter : int32)
where
  reads (r_W)
do
  var filename : int8[1000]
  c.sprintf([&int8](filename), './Data/rho_%04d',iter)
  var g = c.fopen(filename,'wb')

  for e in r_W do
    dumpdouble(g, r_W[e].rho)
  end
  __fence(__execution, __block)
  c.fclose(g)
  return 1
end



task toplevel()
  var config : Config
  config:initialize_from_command() -- TODO : CPUs, output bool

  -- Simulation Parameters
  var testProblem : int32 = config.testproblem
  var r_params = region(ispace(int1d, 1), params)
  TestProblem(r_params, testProblem)
  var N  : int32[3] = r_params[0].N
  var NV : int32[3] = r_params[0].NV
  var Nc : int64 = r_params[0].Nc
  var Nv : int64 = r_params[0].Nv
  var R : double = r_params[0].R
  var K : double = r_params[0].K
  var g : double = r_params[0].g
  var w : double = r_params[0].w
  var Cv : double = r_params[0].Cv
  var ur : double = r_params[0].ur
  var Tr : double = r_params[0].Tr
  var Pr : double = r_params[0].Pr
  var effD : int32 = r_params[0].effD
  var BCs : int32[3] = r_params[0].BCs
  var Vmin : double[3] = r_params[0].Vmin
  var Vmax : double[3] = r_params[0].Vmax
  var Tf : double = r_params[0].Tf
  var dtdump : double = r_params[0].dtdump

  __fence(__execution, __block) 
  c.printf("Simulation Parameters\n")
  if testProblem > 0 then c.printf("testProblem = %d\n", testProblem) end
  c.printf("N = {%d, %d, %d}, NV = {%d, %d, %d}, effD = %d\n", N[0], N[1], N[2], NV[0], NV[1], NV[2], effD)
  c.printf("BCs = {%d, %d, %d}, Vmin = {%f, %f, %f}, Vmax = {%f, %f, %f}\n", BCs[0], BCs[1], BCs[2], Vmin[0], Vmin[1], Vmin[2], Vmax[0], Vmax[1], Vmax[2])
  c.printf("R = %f, K = %f, g = %f, Cv = %f\n", R, K, g, Cv)
  c.printf("w = %f, ur = %f, Tr = %f, Pr = %f\n", w, ur, Tr, Pr)
  c.printf("End Time = %f, dtdump %f\n", Tf, dtdump)

  -- Create regions for distribution functions and gradients
  var r_grid      = region(ispace(int8d, {N[0], N[1], N[2], 1, 1, NV[0], NV[1], NV[2]}), grid)
  var r_gridbarp  = region(ispace(int8d, {N[0], N[1], N[2], 1, 1, NV[0], NV[1], NV[2]}), grid)
  var r_gridbarpb = region(ispace(int8d, {N[0], N[1], N[2], effD, 1, NV[0], NV[1], NV[2]}), grid)
  var r_sig       = region(ispace(int8d, {N[0], N[1], N[2], effD, 1, NV[0], NV[1], NV[2]}), grid)
  var r_sig2      = region(ispace(int8d, {N[0], N[1], N[2], effD, effD, NV[0], NV[1], NV[2]}), grid)
  var r_sigb      = region(ispace(int8d, {N[0], N[1], N[2], effD, effD, NV[0], NV[1], NV[2]}), grid)
 
  -- Create regions for mesh and conserved variables (cell center and interface)
  var r_mesh = region(ispace(int8d, {N[0], N[1], N[2], 1, 1, 1, 1, 1}), mesh)
  var r_W    = region(ispace(int8d, {N[0], N[1], N[2], 1, 1, 1, 1, 1}), W)
  var r_Wb   = region(ispace(int8d, {N[0], N[1], N[2], effD, 1, 1, 1, 1}), W)
 
  -- Create regions for velocity space and initialize
  var vxmesh = region(ispace(int1d, NV[0]), vmesh) 
  var vymesh = region(ispace(int1d, NV[1]), vmesh) 
  var vzmesh = region(ispace(int1d, NV[2]), vmesh) 
  NewtonCotes(vxmesh, vymesh, vzmesh, NV, Vmin, Vmax)

  -- Create regions for source terms and flux
  var r_S = region(ispace(int8d, {N[0], N[1], N[2], 1, 1, NV[0], NV[1], NV[2]}), grid)
  var r_F = region(ispace(int8d, {N[0], N[1], N[2], 1, 1, NV[0], NV[1], NV[2]}), grid)


  -- Create partitions for regions
  var f6 : int6d = factorize(config.cpus, effD)  
  var f8 : int8d = {f6.x, f6.y, f6.z, 1, 1, f6.w, f6.v, f6.u}
  c.printf("Partitioning as {%d, %d, %d, %d, %d, %d, %d, %d}\n", f6.x, f6.y, f6.z, 1, 1, f6.w, f6.v, f6.u)
  var p8 = ispace(int8d, f8)
  var p_grid = partition(equal, r_grid, p8)
  var p_gridbarp = partition(equal, r_gridbarp, p8)
  var p_gridbarpb = partition(equal, r_gridbarpb, p8)
  var p_sig = partition(equal, r_sig, p8)
  var p_sig2 = partition(equal, r_sig2, p8)
  var p_sigb = partition(equal, r_sigb, p8)
  var p_mesh = partition(equal, r_mesh, p8)
  var p_W = partition(equal, r_W, p8)
  var p_Wb = partition(equal, r_Wb, p8)
  var p_S = partition(equal, r_S, p8)
  var p_F = partition(equal, r_F, p8)
  if config.debug == true then
    __fence(__execution, __block)
    c.printf("Equal Partitions Done\n")
  end

  -- Create coloring for partitions for left/right ghost regions
  var c3Lx = coloring.create()
  var c3Ly = coloring.create()
  var c3Lz = coloring.create()
  var c4Lx = coloring.create()
  var c4Ly = coloring.create()
  var c4Lz = coloring.create()
  var c6Lx = coloring.create()
  var c6Ly = coloring.create()
  var c6Lz = coloring.create()
  var c7Lx = coloring.create()
  var c7Ly = coloring.create()
  var c7Lz = coloring.create()
  var c8Lx = coloring.create()
  var c8Ly = coloring.create()
  var c8Lz = coloring.create()

  var c3Rx = coloring.create()
  var c3Ry = coloring.create()
  var c3Rz = coloring.create()
  var c4Rx = coloring.create()
  var c4Ry = coloring.create()
  var c4Rz = coloring.create()
  var c6Rx = coloring.create()
  var c6Ry = coloring.create()
  var c6Rz = coloring.create()
  var c7Rx = coloring.create()
  var c7Ry = coloring.create()
  var c7Rz = coloring.create()
  var c8Rx = coloring.create()
  var c8Ry = coloring.create()
  var c8Rz = coloring.create()
 
  var cvxmesh = coloring.create()
  var cvymesh = coloring.create()
  var cvzmesh = coloring.create()
  -- Create Rects for colorings for partitions
  for col8 in p_gridbarpb.colors do
    var bounds = p_gridbarpb[col8].bounds
    
    -- Leftmost and Rightmost indices
    var il : int32 = bounds.lo.x
    var jl : int32 = bounds.lo.y
    var kl : int32 = bounds.lo.z
    var ir : int32 = bounds.hi.x
    var jr : int32 = bounds.hi.y
    var kr : int32 = bounds.hi.z

    var lx = BC(il, jl, kl, 0, BCs, N)[0]
    var rx = BC(ir, jr, kr, 0, BCs, N)[3]
    var ly = BC(il, jl, kl, 1, BCs, N)[1]
    var ry = BC(ir, jr, kr, 1, BCs, N)[4]
    var lz = BC(il, jl, kl, 2, BCs, N)[2]
    var rz = BC(ir, jr, kr, 2, BCs, N)[5]

    var rleftx3 : rect8d = { {lx, bounds.lo.y, bounds.lo.z, 0, 0, 0, 0, 0}, 
                         {lx, bounds.hi.y, bounds.hi.z, 0, 0, 0, 0, 0}}
    var rlefty3 : rect8d = { {bounds.lo.x, ly, bounds.lo.z, 0, 0, 0, 0, 0},
                         {bounds.hi.x, ly, bounds.hi.z, 0, 0, 0, 0, 0}}
    var rleftz3 : rect8d = { {bounds.lo.x, bounds.lo.y, lz, 0, 0, 0, 0, 0},
                         {bounds.hi.x, bounds.hi.y, lz, 0, 0, 0, 0, 0}}

    var rrightx3 : rect8d = { {rx, bounds.lo.y, bounds.lo.z, 0, 0, 0, 0, 0}, 
                         {rx, bounds.hi.y, bounds.hi.z, 0, 0, 0, 0, 0}}
    var rrighty3 : rect8d = { {bounds.lo.x, ry, bounds.lo.z, 0, 0, 0, 0, 0},
                         {bounds.hi.x, ry, bounds.hi.z, 0, 0, 0, 0, 0}}
    var rrightz3 : rect8d = { {bounds.lo.x, bounds.lo.y, rz, 0, 0, 0, 0, 0},
                         {bounds.hi.x, bounds.hi.y, rz, 0, 0, 0, 0, 0}}

    var rleftx4 : rect8d = { {lx, bounds.lo.y, bounds.lo.z, bounds.lo.w, 0, 0, 0, 0}, 
                         {lx, bounds.hi.y, bounds.hi.z, bounds.hi.w, 0, 0, 0, 0}}
    var rlefty4 : rect8d = { {bounds.lo.x, ly, bounds.lo.z, bounds.lo.w, 0, 0, 0, 0},
                         {bounds.hi.x, ly, bounds.hi.z, bounds.hi.w, 0, 0, 0, 0}}
    var rleftz4 : rect8d = { {bounds.lo.x, bounds.lo.y, lz, bounds.lo.w, 0, 0, 0, 0},
                         {bounds.hi.x, bounds.hi.y, lz, bounds.hi.w, 0, 0, 0, 0}}

    var rrightx4 : rect8d = { {rx, bounds.lo.y, bounds.lo.z, bounds.lo.w, 0, 0, 0, 0}, 
                         {rx, bounds.hi.y, bounds.hi.z, bounds.hi.w, 0, 0, 0, 0}}
    var rrighty4 : rect8d = { {bounds.lo.x, ry, bounds.lo.z, bounds.lo.w, 0, 0, 0, 0},
                         {bounds.hi.x, ry, bounds.hi.z, bounds.hi.w, 0, 0, 0, 0}}
    var rrightz4 : rect8d = { {bounds.lo.x, bounds.lo.y, rz, bounds.lo.w, 0, 0, 0, 0},
                         {bounds.hi.x, bounds.hi.y, rz, bounds.hi.w, 0, 0, 0, 0}}

    var rleftx6 : rect8d = { {lx, bounds.lo.y, bounds.lo.z, 0, 0, bounds.lo.u, bounds.lo.t, bounds.lo.s}, 
                         {lx, bounds.hi.y, bounds.hi.z, 0, 0, bounds.hi.u, bounds.hi.t, bounds.hi.s}}
    var rlefty6 : rect8d = { {bounds.lo.x, ly, bounds.lo.z, 0, 0, bounds.lo.u, bounds.lo.t, bounds.lo.s},
                         {bounds.hi.x, ly, bounds.hi.z, 0, 0, bounds.hi.u, bounds.hi.t, bounds.hi.s}}
    var rleftz6 : rect8d = { {bounds.lo.x, bounds.lo.y, lz, 0, 0, bounds.lo.u, bounds.lo.t, bounds.lo.s},
                         {bounds.hi.x, bounds.hi.y, lz, 0, 0, bounds.hi.u, bounds.hi.t, bounds.hi.s}}

    var rrightx6 : rect8d = { {rx, bounds.lo.y, bounds.lo.z, 0, 0, bounds.lo.u, bounds.lo.t, bounds.lo.s}, 
                         {rx, bounds.hi.y, bounds.hi.z, 0, 0, bounds.hi.u, bounds.hi.t, bounds.hi.s}}
    var rrighty6 : rect8d = { {bounds.lo.x, ry, bounds.lo.z, 0, 0, bounds.lo.u, bounds.lo.t, bounds.lo.s},
                         {bounds.hi.x, ry, bounds.hi.z, 0, 0, bounds.hi.u, bounds.hi.t, bounds.hi.s}}
    var rrightz6 : rect8d = { {bounds.lo.x, bounds.lo.y, rz, 0, 0, bounds.lo.u, bounds.lo.t, bounds.lo.s},
                         {bounds.hi.x, bounds.hi.y, rz, 0, 0, bounds.hi.u, bounds.hi.t, bounds.hi.s}}

    var rleftx7 : rect8d = { {lx, bounds.lo.y, bounds.lo.z, bounds.lo.w, 0, bounds.lo.u, bounds.lo.t, bounds.lo.s}, 
                         {lx, bounds.hi.y, bounds.hi.z, bounds.hi.w, 0, bounds.hi.u, bounds.hi.t, bounds.hi.s}}
    var rlefty7 : rect8d = { {bounds.lo.x, ly, bounds.lo.z, bounds.lo.w, 0, bounds.lo.u, bounds.lo.t, bounds.lo.s},
                         {bounds.hi.x, ly, bounds.hi.z, bounds.hi.w, 0, bounds.hi.u, bounds.hi.t, bounds.hi.s}}
    var rleftz7 : rect8d = { {bounds.lo.x, bounds.lo.y, lz, bounds.lo.w, 0, bounds.lo.u, bounds.lo.t, bounds.lo.s},
                         {bounds.hi.x, bounds.hi.y, lz, bounds.hi.w, 0, bounds.hi.u, bounds.hi.t, bounds.hi.s}}

    var rrightx7 : rect8d = { {rx, bounds.lo.y, bounds.lo.z, bounds.lo.w, 0, bounds.lo.u, bounds.lo.t, bounds.lo.s}, 
                         {rx, bounds.hi.y, bounds.hi.z, bounds.hi.w, 0, bounds.hi.u, bounds.hi.t, bounds.hi.s}}
    var rrighty7 : rect8d = { {bounds.lo.x, ry, bounds.lo.z, bounds.lo.w, 0, bounds.lo.u, bounds.lo.t, bounds.lo.s},
                         {bounds.hi.x, ry, bounds.hi.z, bounds.hi.w, 0, bounds.hi.u, bounds.hi.t, bounds.hi.s}}
    var rrightz7 : rect8d = { {bounds.lo.x, bounds.lo.y, rz, bounds.lo.w, 0, bounds.lo.u, bounds.lo.t, bounds.lo.s},
                         {bounds.hi.x, bounds.hi.y, rz, bounds.hi.w, 0, bounds.hi.u, bounds.hi.t, bounds.hi.s}}

    var rleftx8 : rect8d = { {lx, bounds.lo.y, bounds.lo.z, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t, bounds.lo.s}, 
                         {lx, bounds.hi.y, bounds.hi.z, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t, bounds.hi.s}}
    var rlefty8 : rect8d = { {bounds.lo.x, ly, bounds.lo.z, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t, bounds.lo.s},
                         {bounds.hi.x, ly, bounds.hi.z, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t, bounds.hi.s}}
    var rleftz8 : rect8d = { {bounds.lo.x, bounds.lo.y, lz, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t, bounds.lo.s},
                         {bounds.hi.x, bounds.hi.y, lz, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t, bounds.hi.s}}

    var rrightx8 : rect8d = { {rx, bounds.lo.y, bounds.lo.z, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t, bounds.lo.s}, 
                         {rx, bounds.hi.y, bounds.hi.z, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t, bounds.hi.s}}
    var rrighty8 : rect8d = { {bounds.lo.x, ry, bounds.lo.z, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t, bounds.lo.s},
                         {bounds.hi.x, ry, bounds.hi.z, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t, bounds.hi.s}}
    var rrightz8 : rect8d = { {bounds.lo.x, bounds.lo.y, rz, bounds.lo.w, bounds.lo.v, bounds.lo.u, bounds.lo.t, bounds.lo.s},
                         {bounds.hi.x, bounds.hi.y, rz, bounds.hi.w, bounds.hi.v, bounds.hi.u, bounds.hi.t, bounds.hi.s}}

    var rvx : rect1d = {vxmesh.bounds.lo, vxmesh.bounds.hi}
    var rvy : rect1d = {vymesh.bounds.lo, vymesh.bounds.hi}
    var rvz : rect1d = {vzmesh.bounds.lo, vzmesh.bounds.hi}
    
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("rect8d Done\n")
    end

    -- Color in left strips
    coloring.color_domain(c3Lx, col8, rleftx3)
    coloring.color_domain(c3Ly, col8, rlefty3)
    coloring.color_domain(c3Lz, col8, rleftz3)
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("3d Coloring Done\n")
    end

    coloring.color_domain(c6Lx, col8, rleftx6)
    coloring.color_domain(c6Ly, col8, rlefty6)
    coloring.color_domain(c6Lz, col8, rleftz6)
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("6d Coloring Done\n")
    end


    coloring.color_domain(c7Lx, col8, rleftx7)
    coloring.color_domain(c7Ly, col8, rlefty7)
    coloring.color_domain(c7Lz, col8, rleftz7)
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("7d Coloring Done\n")
    end

    coloring.color_domain(c8Lx, col8, rleftx8)
    coloring.color_domain(c8Ly, col8, rlefty8)
    coloring.color_domain(c8Lz, col8, rleftz8)
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("8d Coloring Done\n")
    end

    -- Color in right strips
    coloring.color_domain(c3Rx, col8, rrightx3)
    coloring.color_domain(c3Ry, col8, rrighty3)
    coloring.color_domain(c3Rz, col8, rrightz3)

    coloring.color_domain(c6Rx, col8, rrightx6)
    coloring.color_domain(c6Ry, col8, rrighty6)
    coloring.color_domain(c6Rz, col8, rrightz6)

    coloring.color_domain(c7Rx, col8, rrightx7)
    coloring.color_domain(c7Ry, col8, rrighty7)
    coloring.color_domain(c7Rz, col8, rrightz7)

    coloring.color_domain(c8Rx, col8, rrightx8)
    coloring.color_domain(c8Ry, col8, rrighty8)
    coloring.color_domain(c8Rz, col8, rrightz8)

    coloring.color_domain(cvxmesh, col8, rvx) 
    coloring.color_domain(cvymesh, col8, rvy) 
    coloring.color_domain(cvzmesh, col8, rvz) 

    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Coloring Done\n")
    end
  end

  -- Create Partitions
  var plx_mesh = partition(disjoint, r_mesh, c3Lx, p8)
  var ply_mesh = partition(disjoint, r_mesh, c3Ly, p8)
  var plz_mesh = partition(disjoint, r_mesh, c3Lz, p8)
  var prx_mesh = partition(disjoint, r_mesh, c3Rx, p8)
  var pry_mesh = partition(disjoint, r_mesh, c3Ry, p8)
  var prz_mesh = partition(disjoint, r_mesh, c3Rz, p8)
  if config.debug == true then
    __fence(__execution, __block)
    c.printf("Mesh Strips Done\n")
  end

  var plx_gridbarp = partition(disjoint, r_gridbarp, c6Lx, p8)
  var ply_gridbarp = partition(disjoint, r_gridbarp, c6Ly, p8)
  var plz_gridbarp = partition(disjoint, r_gridbarp, c6Lz, p8)
  var prx_gridbarp = partition(disjoint, r_gridbarp, c6Rx, p8)
  var pry_gridbarp = partition(disjoint, r_gridbarp, c6Ry, p8)
  var prz_gridbarp = partition(disjoint, r_gridbarp, c6Rz, p8)
  if config.debug == true then
    __fence(__execution, __block)
    c.printf("gridbarp Strips Done\n")
  end

  var plx_gridbarpb = partition(disjoint, r_gridbarpb, c7Lx, p8)
  var ply_gridbarpb = partition(disjoint, r_gridbarpb, c7Ly, p8)
  var plz_gridbarpb = partition(disjoint, r_gridbarpb, c7Lz, p8)
  var prx_gridbarpb = partition(disjoint, r_gridbarpb, c7Rx, p8)
  var pry_gridbarpb = partition(disjoint, r_gridbarpb, c7Ry, p8)
  var prz_gridbarpb = partition(disjoint, r_gridbarpb, c7Rz, p8)
  if config.debug == true then
    __fence(__execution, __block)
    c.printf("gridbarpbb Strips Done\n")
  end

  var plx_sig = partition(disjoint, r_sig, c7Lx, p8)
  var ply_sig = partition(disjoint, r_sig, c7Ly, p8)
  var plz_sig = partition(disjoint, r_sig, c7Lz, p8)
  var prx_sig = partition(disjoint, r_sig, c7Rx, p8)
  var pry_sig = partition(disjoint, r_sig, c7Ry, p8)
  var prz_sig = partition(disjoint, r_sig, c7Rz, p8)
  if config.debug == true then
    __fence(__execution, __block)
    c.printf("sig Strips Done\n")
  end

  var plx_sig2 = partition(disjoint, r_sig2, c8Lx, p8)
  var ply_sig2 = partition(disjoint, r_sig2, c8Ly, p8)
  var plz_sig2 = partition(disjoint, r_sig2, c8Lz, p8)
  var prx_sig2 = partition(disjoint, r_sig2, c8Rx, p8)
  var pry_sig2 = partition(disjoint, r_sig2, c8Ry, p8)
  var prz_sig2 = partition(disjoint, r_sig2, c8Rz, p8)
  if config.debug == true then
    __fence(__execution, __block)
    c.printf("sig2 Strips Done\n")
  end
 
  var plx_sigb = partition(disjoint, r_sigb, c8Lx, p8)
  var ply_sigb = partition(disjoint, r_sigb, c8Ly, p8)
  var plz_sigb = partition(disjoint, r_sigb, c8Lz, p8)
  var prx_sigb = partition(disjoint, r_sigb, c8Rx, p8)
  var pry_sigb = partition(disjoint, r_sigb, c8Ry, p8)
  var prz_sigb = partition(disjoint, r_sigb, c8Rz, p8)
  if config.debug == true then
    __fence(__execution, __block)
    c.printf("sigb Strips Done\n")
    c.printf("Left/Right Strip Partitioning Done\n")
  end

  var pxmesh = partition(aliased, vxmesh, cvxmesh, p8)
  var pymesh = partition(aliased, vymesh, cvymesh, p8)
  var pzmesh = partition(aliased, vzmesh, cvzmesh, p8)

  --Initialize r_mesh
  var MeshType : int32 = 1
  InitializeMesh(r_mesh, N, MeshType) --TODO Needs more input for nested, user-def etc.
  if config.debug == true then
    __fence(__execution, __block)
    c.printf("Mesh Initialized\n")
  end

  --Initialize r_W
  __demand(__index_launch)
  for col8 in p_W.colors do
    InitializeW(p_W[col8], p_mesh[col8], N, NV, testProblem, R, Cv)
  end
  if config.debug == true then
    __fence(__execution, __block)
    c.printf("W Initialized\n")
  end
 
  --Initialize r_grid
  __demand(__index_launch)
  for col8 in p_grid.colors do
    InitializeGrid(p_grid[col8], p_mesh[col8], p_W[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], testProblem, R, K, Cv, g, w, ur, Tr, Pr, N, NV, effD)
  end
  if config.debug == true then
    __fence(__execution, __block)
    c.printf("Grid Initialized\n")
  end

  --Timestep
  var CFL : double = 0.8 -- Safety Factor
  var dxmin : double = 1.0/cmath.fmax(cmath.fmax(N[0],N[1]),N[2]) -- Smallest Cell Width (TODO : Non-Uniform Meshes)
  var umax : double  = 4.0 -- Estimated maximum flow velocity, TODO calculate at each iteration for stronger problems
  var calcdt : double = CFL*dxmin/(umax + sqrt(Vmax[0]*Vmax[0] + Vmax[1]*Vmax[1] + Vmax[2]*Vmax[2]))
  
  var Tsim : double = 0.0  -- Sim time
  var Tdump : double = 0.0 -- Time since last dump 
  
  var iter : int32 = 0
  var dumpiter : int32 = 0
  if testProblem > 0 then 
    --__fence(__execution, __block)
    --c.printf("Dumping %d\n", dumpiter)

    --__fence(__execution, __block)
    Dump(r_W, dumpiter) -- Initial Conditions

    --__fence(__execution, __block)
    c.printf("Dump %d\n", dumpiter)
  end
  
  var Start : double = c.legion_get_current_time_in_nanos()
  var End : double 
  
  while Tsim < Tf do --and iter < 10 do
    iter += 1

    var dt = TimeStep(calcdt, dtdump-Tdump, Tf-Tsim)

    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Starting Step1a\n")
      c.fflush(c.stdout)
    end
    -- Step 1a
    __demand(__index_launch)
    for col8 in p_grid.colors do 
      Step1a(p_grid[col8], p_gridbarp[col8], p_S[col8], p_W[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], dt, R, K, Cv, g, w, ur, Tr, Pr, effD)
    end
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Step1a Complete\n")
      c.fflush(c.stdout)
    end

    -- Step 1b: Compute Gradient Sigma
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Computing Gradients\n")
      c.fflush(c.stdout)
    end
    __demand(__index_launch)
    for col8 in p_grid.colors do 
      Step1b_sigx(p_gridbarp[col8], p_sig[col8], p_mesh[col8], plx_mesh[col8], prx_mesh[col8], plx_gridbarp[col8], prx_gridbarp[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N, effD)
    end
    if effD > 1 then
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigy(p_gridbarp[col8], p_sig[col8], p_mesh[col8], ply_mesh[col8], pry_mesh[col8], ply_gridbarp[col8], pry_gridbarp[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N, effD)
      end  
    end
    if effD > 2 then
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigz(p_gridbarp[col8], p_sig[col8], p_mesh[col8], plz_mesh[col8], prz_mesh[col8], plz_gridbarp[col8], prz_gridbarp[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N, effD)
      end  
    end
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Gradients Computed\n")
      c.fflush(c.stdout)
    end

    -- Step 1b: Compute Gradient of Gradient Sigma, Sigma2
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Computing Gradients of Gradients\n")
      c.fflush(c.stdout)
    end
    __demand(__index_launch)
    for col8 in p_grid.colors do 
      Step1b_sigx_x(p_sig[col8], p_sig2[col8], p_mesh[col8], plx_mesh[col8], prx_mesh[col8], plx_sig[col8], prx_sig[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
    end  
    if effD > 1 then
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigx_y(p_sig[col8], p_sig2[col8], p_mesh[col8], ply_mesh[col8], pry_mesh[col8], ply_sig[col8], pry_sig[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigy_y(p_sig[col8], p_sig2[col8], p_mesh[col8], ply_mesh[col8], pry_mesh[col8], ply_sig[col8], pry_sig[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigy_x(p_sig[col8], p_sig2[col8], p_mesh[col8], plx_mesh[col8], prx_mesh[col8], plx_sig[col8], prx_sig[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
    end
    if effD > 2 then
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigz_x(p_sig[col8], p_sig2[col8], p_mesh[col8], plx_mesh[col8], prx_mesh[col8], plx_sig[col8], prx_sig[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigz_y(p_sig[col8], p_sig2[col8], p_mesh[col8], ply_mesh[col8], pry_mesh[col8], ply_sig[col8], pry_sig[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigz_z(p_sig[col8], p_sig2[col8], p_mesh[col8], plz_mesh[col8], prz_mesh[col8], plz_sig[col8], prz_sig[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigx_z(p_sig[col8], p_sig2[col8], p_mesh[col8], plz_mesh[col8], prz_mesh[col8], plz_sig[col8], prz_sig[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigy_z(p_sig[col8], p_sig2[col8], p_mesh[col8], plz_mesh[col8], prz_mesh[col8], plz_sig[col8], prz_sig[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
    end
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Gradients of Gradients Computed\n")
      c.fflush(c.stdout)
    end

    --Step1b: Compute Sig at Boundary, Sigb
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Computing Gradient at Interface\n")
      c.fflush(c.stdout)
    end
    __demand(__index_launch)
    for col8 in p_grid.colors do 
      Step1b_sigx_x2(p_sig[col8], p_sig2[col8], p_sigb[col8], p_mesh[col8], prx_mesh[col8], prx_sig[col8], prx_sig2[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
    end  

    if effD > 1 then
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigy_x2(p_sig[col8], p_sig2[col8], p_sigb[col8], p_mesh[col8], prx_mesh[col8], prx_sig[col8], prx_sig2[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigx_y2(p_sig[col8], p_sig2[col8], p_sigb[col8], p_mesh[col8], pry_mesh[col8], pry_sig[col8], pry_sig2[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigy_y2(p_sig[col8], p_sig2[col8], p_sigb[col8], p_mesh[col8], pry_mesh[col8], pry_sig[col8], pry_sig2[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
    end

    if effD > 2 then
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigz_x2(p_sig[col8], p_sig2[col8], p_sigb[col8], p_mesh[col8], prx_mesh[col8], prx_sig[col8], prx_sig2[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigz_y2(p_sig[col8], p_sig2[col8], p_sigb[col8], p_mesh[col8], pry_mesh[col8], pry_sig[col8], pry_sig2[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigx_z2(p_sig[col8], p_sig2[col8], p_sigb[col8], p_mesh[col8], prz_mesh[col8], prz_sig[col8], prz_sig2[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigy_z2(p_sig[col8], p_sig2[col8], p_sigb[col8], p_mesh[col8], prz_mesh[col8], prz_sig[col8], prz_sig2[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
      __demand(__index_launch)
      for col8 in p_grid.colors do 
        Step1b_sigz_z2(p_sig[col8], p_sig2[col8], p_sigb[col8], p_mesh[col8], prz_mesh[col8], prz_sig[col8], prz_sig2[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N)
      end  
    end
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Gradients at boundary computed\n")
      c.fflush(c.stdout)
    end

    -- Step 1b_b: Interpolate boundary
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Interpolating to boundary\n")
      c.fflush(c.stdout)
    end
    __demand(__index_launch)
    for col8 in p_gridbarpb.colors do
      Step1b_b(p_sig[col8], p_mesh[col8], p_gridbarp[col8], p_gridbarpb[col8], prx_gridbarp[col8], pry_gridbarp[col8], prz_gridbarp[col8], prx_sig[col8], pry_sig[col8], prz_sig[col8], prx_mesh[col8], pry_mesh[col8], prz_mesh[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], BCs, N, effD)
    end
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Interpolating to boundary Complete\n")
      c.fflush(c.stdout)
    end

    -- Step 1c
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Interpolating to velocity-dependent past at boundary\n")
      c.fflush(c.stdout)
    end
    __demand(__index_launch)
    for col8 in p_grid.colors do
      Step1c(p_gridbarpb[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], p_sigb[col8], dt, BCs, N, effD)
    end
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Interpolating to Past Complete\n")
      c.fflush(c.stdout)
    end

    -- Step 2a
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Computing Wb\n")
      c.fflush(c.stdout)
    end
    __demand(__index_launch)
    for col8 in p_gridbarpb.colors do
      Step2a(p_gridbarpb[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], p_Wb[col8], dt, effD)
    end
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Computed Wb\n")
      c.fflush(c.stdout)
    end

    -- Step 2b
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Computing original phi at boundary\n")
      c.fflush(c.stdout)
    end
    __demand(__index_launch)
    for col8 in p_gridbarpb.colors do
      Step2b(p_gridbarpb[col8], p_Wb[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], dt, R, K, Cv, g, w, ur, Tr, Pr, effD)
    end
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Computed original phi at boundary\n")
      c.fflush(c.stdout)
    end

    -- Step 2c
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Computing Microflux\n")
      c.fflush(c.stdout)
    end
    __demand(__index_launch)
    for col8 in p_gridbarpb.colors do
      Step2c(p_gridbarpb[col8], p_F[col8], p_mesh[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], plx_gridbarpb[col8], ply_gridbarpb[col8], plz_gridbarpb[col8], BCs, R, K, Cv, g, w, Tr, Pr, effD, N)
    end
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Computed Microflux\n")
      c.fflush(c.stdout)
    end

    -- Step 3
    Step3() -- TODO
  
    -- Step 4 and 5
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Updating W and Phi\n")
      c.fflush(c.stdout)
    end
    __demand(__index_launch)
    for col8 in p_grid.colors do
      Step4and5(p_grid[col8], p_W[col8], p_mesh[col8], p_F[col8], pxmesh[col8], pymesh[col8], pzmesh[col8], dt, BCs, R, K, Cv, N, g, w, ur, Tr, Pr, effD)
    end
    if config.debug == true then
      __fence(__execution, __block)
      c.printf("Updated W and Phi\n")
      c.fflush(c.stdout)
    end
    if dt < calcdt then
      dumpiter += 1
      Tdump = 0

      --__fence(__execution, __block)
      --c.printf("Dumping %d\n", dumpiter)
      --c.fflush(c.stdout)

      --__fence(__execution, __block)
      Dump(r_W, dumpiter)

      --__fence(__execution, __block)
      c.printf("Dump %d\n", dumpiter)
      c.fflush(c.stdout)
    else
      Tdump += dt
    end

    Tsim += dt

    __fence(__execution, __block)
    End = c.legion_get_current_time_in_nanos()
    c.printf("Iteration = %d, Tsim = %f, Realtime = %f\n", iter, Tsim, (End-Start)*1e-9)
    c.fflush(c.stdout)
  end

  __fence(__execution, __block)
  End = c.legion_get_current_time_in_nanos()
  c.printf("Finished simulation in %.4f seconds.\n", (End-Start)*1e-9)
  c.fflush(c.stdout)
end

regentlib.start(toplevel)
