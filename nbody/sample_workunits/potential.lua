arg = { ... }

assert(#arg == 6, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed")

prng = DSFMT.create(argSeed)

evolveTime       = arg[1]
reverseOrbitTime = arg[1] / arg[2]

r0  = arg[3]
light_r_ratio = arg[4]

dwarfMass  = arg[5]
light_mass_ratio = arg[6]

model1Bodies = 10000
totalBodies = model1Bodies

nbodyLikelihoodMethod = "EMD"
nbodyMinVersion = "1.32"

function mw_lua_mulvs(x, y, z, s)
  return s * x, s * y, s * z
end

function mw_lua_absv(x,y,z)
  return sqrt(x*x + y*y + z*z)
end

function mw_lua_incaddv(x_1, y_1, z_1, x_2, y_2, z_2)
  return (x_1 + x_2), (y_1 + y_2), (z_1 + z_2)
end

function sphericalAccel(mass, scale, x, y, z, r)
  local tmp = scale + r

  return mw_lua_mulvs(x, y, z, (-mass / (r * tmp^2)))
end

function miyamotoNagaiDiskAccel(mass, scaleLength, scaleHeight, x, y, z, r)
  local zp  = sqrt(z^2 + scaleHeight^2)
  local azp = scaleLength + zp

  local rp  = x^2 + y^2 + azp^2
  local rth = sqrt(rp^3)

  local a_x = -mass * x / rth
  local a_y = -mass * y / rth
  local a_z = -mass * z * azp / (zp * rth)

  return a_x, a_y, a_z
end

function exponentialDiskAccel(mass, scaleLength, scaleHeight, x, y, z, r)
  local expPiece = exp(-r / scaleLength) * (r +scaleLength) / scaleLength
  local factor   = mass * (expPiece - 1.0) / r^3

  return mw_lua_mulvs(x, y, z, factor)
end

function logHaloAccel(vhalo, scaleLength, flattenZ, x, y, z, r)
  local tvsqr = -2.0 * vhalo^2
  local qsqr  = flattenZ^2
  local zsqr  = z^2

  local arst  = scaleLength^2 + x^2 + y^2
  local denom = (zsqr / qsqr) + arst

  local a_x = tvsqr * x / denom
  local a_y = tvsqr * y / denom
  local a_z = tvsqr * z / ((qsqr * arst) + zsqr)

  return a_x, a_y, a_z
end

function nfwHaloAccel(vhalo, scaleLength, flattenZ, x, y, z, r)
  local ar = scaleLength + r
  local c = a * vhalo^2 * (r - ar * log((r + a) / a)) / (0.2162165954 * r^3 * ar)

  return mw_lua_mulvs(x, y, z, c)
end

function getPotential()
   return  Potential.create{
      spherical = Spherical.spherical{ mass  = 1.52954402e5, scale = 0.7 },
      disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
      halo      = Halo.logarithmic{ vhalo = 73, scaleLength = 12.0, flattenZ = 1.0 }
   }
end

function makePotential()

  local self = {
    spherical = { mass  = 1.52954402e5, scale = 0.7,
                  accFunc = sphericalAccel },
    disk      = { mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26,
                  accFunc = miyamotoNagaiDiskAccel },
    halo      = { vhalo = 73, scaleLength = 12.0, flattenZ = 1.0,
                  accFunc = logHaloAccel }
  }


  local nbExtAcceleration = function(x,y,z)
    -- Return variables that contain the acceleration
    -- In all 3 Dimensions
    local a_x, a_y, a_z

    -- These will contain the acceleration due to the Halo
    local tmp_a_x, tmp_a_y, tmp_a_z

    -- Radial distance
    local r = mw_lua_absv(x,y,z)

    a_x, a_y, a_z = self.disk.accFunc(self.disk.mass, self.disk.scaleLength, self.disk.scaleHeight, x, y, z, r)

   tmp_a_x, tmp_a_y, tmp_a_z = self.halo.accFunc(self.halo.vhalo, self.halo.scaleLength, self.halo.flattenZ, x, y, z, r);

    -- Add the accelerations
    a_x, a_y, a_z = mw_lua_incaddv(a_x, a_y, a_z, tmp_a_x, tmp_a_y, tmp_a_z)

    -- Calculate the Bulge Accelerations
    tmp_a_x, tmp_a_y, tmp_a_z = self.spherical.accFunc(self.spherical.mass, self.spherical.scale, x, y, z, r);

    -- Add the accelerations
    a_x, a_y, a_z = mw_lua_incaddv(a_x, a_y, a_z, tmp_a_x, tmp_a_y, tmp_a_z)

    return a_x, a_y, a_z
  end

  return nbExtAcceleration
  -- return getPotential()
end

encMass = plummerTimestepIntegral(r0*light_r_ratio, sqr(r0) + sqr(r0/arg[4]) , dwarfMass, 1e-7)

function makeContext()
   return NBodyCtx.create{
      timeEvolve = evolveTime,
      timestep   = sqr(1/10.0) * sqrt((pi_4_3 * cube(r0)) / (encMass + dwarfMass)),
      eps2       = calculateEps2(totalBodies, r0),
      criterion  = "NewCriterion",
      useQuad    = true,
      theta      = 1.0
   }
end

-- Also required
function makeBodies(ctx, potential)
    local firstModel
    local finalPosition, finalVelocity = reverseOrbit{
        potential = potential or getPotential(),
        position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6)),
        velocity  = Vector.create(-156, 79, 107),
        tstop     = reverseOrbitTime,
        dt        = ctx.timestep / 10.0
    }

    firstModel = predefinedModels.isotropic{
        nbody       = model1Bodies,
        prng        = prng,
        position    = finalPosition,
        velocity    = finalVelocity,
        mass1        = dwarfMass * arg[6],
	mass2       = dwarfMass - (dwarfMass * arg[6]),
        scaleRadius1 = r0,
	scaleRadius2 = r0/arg[4],
        ignore      = true
    }

    if (tonumber(arg[6]) == 1.0) then
      for i,v in ipairs(firstModel)
      do
     	v.ignore = false
      end


    else
      count = 0
      while (count < arg[6] * totalBodies)
      do

      for i,v in ipairs(firstModel)
      do
         --Figure out properties
         r_2 = (finalPosition.x - v.position.x)^2 + (finalPosition.y - v.position.y)^2 + (finalPosition.z - v.position.z)^2
         r = r_2 ^ (0.5)
         scale1 = r0
         lightDensityToMaxRatio = 1/((1 + r^2/scale1^2)^(5/2))

         --Get light sphere properties

         --Chance object is light_mass = light/dwarf
         chance = lightDensityToMaxRatio

         ---Do random calculation
         if(prng:random() > chance and v.ignore==true and count < arg[6] * totalBodies)
          then
              v.ignore=false
              count = count + 1
          end
       end
     end
   end

return firstModel
end

function makeHistogram()
    return HistogramParams.create{
     phi = 128.79,
     theta = 54.39,
     psi = 90.70,
     lambdaStart = -300,
     lambdaEnd = 300,
     lambdaBins = 50,
     betaStart = -20,
     betaEnd = 20,
     betaBins = 1
}
end
