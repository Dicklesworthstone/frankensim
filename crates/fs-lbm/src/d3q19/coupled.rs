//! WS1 3-D feature-preservation ports (bead eg62): thermal double
//! population and non-Newtonian local rheology on D3Q19.
//!
//! The 2-D crate VALIDATED thermal Rayleigh–Bénard (onset bracket, Nu > 1)
//! and power-law rheology; this module ports both to the 3-D path without
//! forking constitutive code: the apparent-viscosity law, τ floor/cap and
//! ledger come from [`crate::rheology`], and the plane-channel analytic
//! profile is the SAME [`crate::rheology::powerlaw_poiseuille_analytic`]
//! (the plate geometry is one-dimensional in the wall-normal axis).
//!
//! Flow cells live on a plain scalar `Vec` grid ([`PlatesGrid3`]) with
//! per-cell relaxation time and per-cell external force — the 3-D sibling
//! of `core2::Grid`, deliberately NOT the tiled SoA `Duct`/`BoundaryGrid3`
//! fast path: those carry frozen bit-semantics goldens and stay untouched;
//! this is the correctness-first reference rung, mirroring how the 2-D
//! extensions run on the plain grid. Collision reuses the shared
//! [`collide_cell3`](super::collide_cell3) kernel (BGK + Guo forcing), so the collision algebra
//! is written exactly once.
//!
//! Temperature rides a D3Q7 advection-diffusion population (rest + six
//! axis directions, c_s² = 1/4) advected by the flow and coupled back
//! through a Boussinesq body force in the Guo term — the double-
//! distribution scheme the 2-D module uses, with anti-bounce-back
//! fixed-temperature plates.
//!
//! No-claims: single-threaded reference implementation; no SIMD/tile
//! claim, no turbulence model, no golden is frozen here (the seeded
//! determinism hash stays a replay-checked CANDIDATE until the
//! GOLDEN_POLICY four-quadrant ceremony, as in the 40p2 precedent).

use super::{CollisionModel3, E3, OPP3, Q3, collision_kernel3, equilibrium3};
use crate::rheology::{Rheology, RheologyStats};
use crate::thermal::{ThermalForcingError, forcing::validate_external_forces};

/// Flow-lattice sound speed squared (D3Q19).
const CS2_F: f64 = 1.0 / 3.0;

/// D3Q7 temperature-lattice velocities: rest + six axis directions.
const E7: [(i32, i32, i32); 7] = [
    (0, 0, 0),
    (1, 0, 0),
    (-1, 0, 0),
    (0, 1, 0),
    (0, -1, 0),
    (0, 0, 1),
    (0, 0, -1),
];
/// D3Q7 weights (rest 1/4, axes 1/8): c_s² = 1/4.
const W7: [f64; 7] = [0.25, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125];
/// D3Q7 opposite table.
const OPP7: [usize; 7] = [0, 2, 1, 4, 3, 6, 5];
/// D3Q7 sound speed squared.
const CS2_T7: f64 = 0.25;

/// D3Q7 advection-diffusion equilibrium.
fn geq3(t: f64, u: [f64; 3]) -> [f64; 7] {
    let mut g = [0.0f64; 7];
    for q in 0..7 {
        let eu = f64::from(E7[q].0) * u[0] + f64::from(E7[q].1) * u[1] + f64::from(E7[q].2) * u[2];
        g[q] = W7[q] * t * (1.0 + eu / CS2_T7);
    }
    g
}

/// A plain-storage D3Q19 grid between two rigid z-plates: wall layers at
/// `z = 0` and `z = nz − 1` (halfway bounce-back planes at z = ½ and
/// z = nz − 3⁄2), periodic in x and y, per-cell relaxation time and
/// per-cell external force through the shared Guo kernel.
pub struct PlatesGrid3 {
    /// Cells along x (periodic).
    pub nx: usize,
    /// Cells along y (periodic).
    pub ny: usize,
    /// Cells along z INCLUDING the two wall layers.
    pub nz: usize,
    /// Distributions, index `(z·ny + y)·nx + x`.
    pub f: Vec<[f64; Q3]>,
    /// Per-cell relaxation time (rheology writes here).
    pub tau: Vec<f64>,
    /// Per-cell external force (Boussinesq / body force writes here).
    pub fext: Vec<[f64; 3]>,
    post: Vec<[f64; Q3]>,
}

impl PlatesGrid3 {
    /// A resting unit-density plate channel: `nz_fluid` fluid layers
    /// between two wall layers, uniform relaxation time `tau`.
    ///
    /// # Panics
    /// If any dimension is zero, `nz_fluid < 2`, or `tau` is not finite
    /// and greater than one half.
    #[must_use]
    pub fn plates(nx: usize, ny: usize, nz_fluid: usize, tau: f64) -> PlatesGrid3 {
        assert!(
            nx > 0 && ny > 0 && nz_fluid >= 2,
            "plate channel needs positive x/y extents and at least two fluid layers"
        );
        assert!(
            tau.is_finite() && tau > 0.5,
            "flow relaxation time tau must be finite and greater than 0.5"
        );
        let nz = nz_fluid + 2;
        let cells = nx * ny * nz;
        let f0 = equilibrium3(1.0, [0.0; 3]);
        PlatesGrid3 {
            nx,
            ny,
            nz,
            f: vec![f0; cells],
            tau: vec![tau; cells],
            fext: vec![[0.0; 3]; cells],
            post: vec![f0; cells],
        }
    }

    /// Linear cell index.
    #[inline]
    #[must_use]
    pub fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        (z * self.ny + y) * self.nx + x
    }

    /// Whether layer `z` is a wall layer.
    #[inline]
    #[must_use]
    pub fn is_wall_layer(&self, z: usize) -> bool {
        z == 0 || z == self.nz - 1
    }

    /// Density, momentum-velocity (with the Guo half-force correction),
    /// at cell index `i`.
    #[must_use]
    pub fn moments3(&self, i: usize) -> (f64, [f64; 3]) {
        let mut rho = 0.0f64;
        let mut m = [0.0f64; 3];
        for (q, e) in E3.iter().enumerate() {
            let fq = self.f[i][q];
            rho += fq;
            m[0] += f64::from(e.0) * fq;
            m[1] += f64::from(e.1) * fq;
            m[2] += f64::from(e.2) * fq;
        }
        let u = [
            (m[0] + 0.5 * self.fext[i][0]) / rho,
            (m[1] + 0.5 * self.fext[i][1]) / rho,
            (m[2] + 0.5 * self.fext[i][2]) / rho,
        ];
        (rho, u)
    }

    /// Total mass over all cells (conserved up to roundoff: collision and
    /// Guo forcing are mass-neutral, walls bounce back, x/y wrap).
    #[must_use]
    pub fn total_mass(&self) -> f64 {
        self.f.iter().flat_map(|f| f.iter()).sum()
    }

    /// Deterministically perturb fluid-cell densities with the same
    /// integer-hash schedule the tiled `Duct` uses (no RNG).
    pub fn perturb(&mut self, seed: u64, amplitude: f64) {
        for z in 1..self.nz - 1 {
            for y in 0..self.ny {
                for x in 0..self.nx {
                    let mut h = seed
                        ^ (x as u64)
                            .wrapping_mul(0x9e37_79b9_7f4a_7c15)
                            .wrapping_add((y as u64).wrapping_mul(0xbf58_476d_1ce4_e5b9))
                            .wrapping_add((z as u64).wrapping_mul(0x94d0_49bb_1331_11eb));
                    h ^= h >> 30;
                    h = h.wrapping_mul(0xbf58_476d_1ce4_e5b9);
                    h ^= h >> 27;
                    let unit = (h >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0;
                    let rho = 1.0 + amplitude * unit;
                    let i = self.idx(x, y, z);
                    self.f[i] = equilibrium3(rho, [0.0; 3]);
                }
            }
        }
    }

    /// One collide-force-stream step: shared BGK + Guo collision per
    /// fluid cell (per-cell τ and force), then pull streaming with
    /// halfway bounce-back at the z-plates and periodic x/y wrap.
    /// Traversal is x-fastest, then y, then z — pinned and
    /// single-threaded.
    pub fn step(&mut self) {
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        let collide = collision_kernel3();
        for z in 1..nz - 1 {
            for y in 0..ny {
                for x in 0..nx {
                    let i = self.idx(x, y, z);
                    self.post[i] = collide(
                        self.f[i],
                        CollisionModel3::Bgk { tau: self.tau[i] },
                        self.fext[i],
                    )
                    .expect("plate-channel state stays admissible for BGK/Guo collision");
                }
            }
        }
        for z in 1..nz - 1 {
            for y in 0..ny {
                for x in 0..nx {
                    let i = self.idx(x, y, z);
                    for q in 0..Q3 {
                        let (ex, ey, ez) = E3[q];
                        let sx = (x as i64 - i64::from(ex)).rem_euclid(nx as i64) as usize;
                        let sy = (y as i64 - i64::from(ey)).rem_euclid(ny as i64) as usize;
                        let sz = z as i64 - i64::from(ez);
                        let sz = sz as usize; // 1..=nz-2 ± 1 stays in 0..nz
                        self.f[i][q] = if self.is_wall_layer(sz) {
                            // Halfway bounce-back off a resting plate.
                            self.post[i][OPP3[q]]
                        } else {
                            self.post[self.idx(sx, sy, sz)][q]
                        };
                    }
                }
            }
        }
    }
}

/// The local shear-rate magnitude from the D3Q19 non-equilibrium second
/// moment — the exact 3-D extension of the 2-D convention:
/// `S_ab = −(3/(2ρτ)) Σ_q e_qa e_qb (f_q − feq_q)`, `γ̇ = √(2 S:S)`.
/// This raw estimator assumes unforced populations. Guo-forced consumers
/// must include the half-source stress correction.
#[must_use]
pub fn shear_rate3(f: &[f64; Q3], feq: &[f64; Q3], rho: f64, tau: f64) -> f64 {
    shear_rate_from_stress3(nonequilibrium_stress3(f, feq), rho, tau)
}

fn nonequilibrium_stress3(f: &[f64; Q3], feq: &[f64; Q3]) -> [f64; 6] {
    let mut sxx = 0.0f64;
    let mut syy = 0.0f64;
    let mut szz = 0.0f64;
    let mut sxy = 0.0f64;
    let mut sxz = 0.0f64;
    let mut syz = 0.0f64;
    for q in 0..Q3 {
        let neq = f[q] - feq[q];
        let (ex, ey, ez) = (f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2));
        sxx += ex * ex * neq;
        syy += ey * ey * neq;
        szz += ez * ez * neq;
        sxy += ex * ey * neq;
        sxz += ex * ez * neq;
        syz += ey * ez * neq;
    }
    [sxx, syy, szz, sxy, sxz, syz]
}

fn shear_rate_from_stress3([sxx, syy, szz, sxy, sxz, syz]: [f64; 6], rho: f64, tau: f64) -> f64 {
    let c = -3.0 / (2.0 * rho * tau);
    let (sxx, syy, szz) = (c * sxx, c * syy, c * szz);
    let (sxy, sxz, syz) = (c * sxy, c * sxz, c * syz);
    let ss = sxx * sxx + syy * syy + szz * szz + 2.0 * (sxy * sxy + sxz * sxz + syz * syz);
    (2.0 * ss).sqrt()
}

/// Molecular strain rate in the physical Guo frame, shared by fixed-domain
/// and moving-liquid thermal consumers. The half-source correction restores
/// `(u F^T + F u^T)/2` before extracting the non-equilibrium stress. `u` is
/// already force-corrected; `force` is the force density in this population
/// epoch, without the collision prefactor `1 - 1/(2 tau)`.
pub(crate) fn force_corrected_shear_rate3(
    f: &[f64; Q3],
    rho: f64,
    u: [f64; 3],
    force: [f64; 3],
    tau: f64,
) -> f64 {
    let [sxx, syy, szz, sxy, sxz, syz] = nonequilibrium_stress3(f, &equilibrium3(rho, u));
    // Apply the exact source second moment directly: six scalar corrections
    // replace a full D3Q19 source construction and second population pass.
    shear_rate_from_stress3(
        [
            sxx + u[0] * force[0],
            syy + u[1] * force[1],
            szz + u[2] * force[2],
            sxy + 0.5 * (u[0] * force[1] + u[1] * force[0]),
            sxz + 0.5 * (u[0] * force[2] + u[2] * force[0]),
            syz + 0.5 * (u[1] * force[2] + u[2] * force[1]),
        ],
        rho,
        tau,
    )
}

/// Update every fluid cell's τ from its local shear rate under `law` —
/// the SAME self-consistent constitutive solve, floor, cap, and ledger as the 2-D
/// [`crate::rheology::update_tau`] (shared, not forked).
#[must_use]
pub fn update_tau3(grid: &mut PlatesGrid3, law: Rheology) -> RheologyStats {
    let law = law.prepare();
    let mut stats = RheologyStats {
        floored: 0,
        capped: 0,
        max_shear: 0.0,
        tau_range: (f64::INFINITY, 0.0),
    };
    for z in 1..grid.nz - 1 {
        for y in 0..grid.ny {
            for x in 0..grid.nx {
                let i = grid.idx(x, y, z);
                let (rho, u) = grid.moments3(i);
                let numerator = force_corrected_shear_rate3(&grid.f[i], rho, u, grid.fext[i], 1.0);
                let local = law.solve(numerator, grid.tau[i]);
                grid.tau[i] = local.tau;
                stats.floored += u32::from(local.floored);
                stats.capped += u32::from(local.capped);
                stats.max_shear = stats.max_shear.max(local.shear_rate);
                stats.tau_range.0 = stats.tau_range.0.min(local.tau);
                stats.tau_range.1 = stats.tau_range.1.max(local.tau);
            }
        }
    }
    if stats.tau_range.0 == f64::INFINITY {
        stats.tau_range = (0.0, 0.0);
    }
    stats
}

/// Run a force-driven plate channel (force along x, walls at the
/// z-plates, periodic x/y) with local rheology updates until the
/// center velocity stabilizes; returns (grid, steps, final stats).
/// The steady profile is one-dimensional in z, so it is certified
/// against the SHARED [`crate::rheology::powerlaw_poiseuille_analytic`].
///
/// # Panics
/// On inadmissible dimensions or a non-finite body force.
#[must_use]
pub fn plate_channel_flow3(
    nx: usize,
    ny: usize,
    nz_fluid: usize,
    gx: f64,
    law: Rheology,
    max_steps: usize,
) -> (PlatesGrid3, usize, RheologyStats) {
    assert!(gx.is_finite(), "body force must be finite");
    let mut grid = PlatesGrid3::plates(nx, ny, nz_fluid, 1.0);
    for i in 0..grid.fext.len() {
        grid.fext[i] = [gx, 0.0, 0.0];
    }
    let mut stats = RheologyStats::default();
    let mut last_peak = 0.0f64;
    let mut steps = 0;
    let mid = grid.idx(0, 0, (nz_fluid + 2) / 2);
    for s in 0..max_steps {
        stats = update_tau3(&mut grid, law);
        grid.step();
        steps = s + 1;
        if s % 200 == 199 {
            let peak = grid.moments3(mid).1[0];
            if (peak - last_peak).abs() < 1e-9 * peak.abs().max(1e-30) {
                break;
            }
            last_peak = peak;
        }
    }
    (grid, steps, stats)
}

/// Double-population thermal D3Q19 lattice: rigid plates at fixed
/// temperatures (bottom hot), periodic x/y, Boussinesq buoyancy along z
/// through the shared Guo term — the 3-D port of the 2-D `ThermalLbm`.
pub struct ThermalLbm3 {
    /// The flow grid (wall layers at z = 0 and z = nz − 1).
    pub grid: PlatesGrid3,
    /// D3Q7 temperature populations.
    pub gpop: Vec<[f64; 7]>,
    /// Thermal relaxation time (α = c_s²(τ_g − ½), c_s² = ¼).
    pub tau_g: f64,
    /// Bottom-plate temperature.
    pub t_bottom: f64,
    /// Top-plate temperature.
    pub t_top: f64,
    /// Buoyancy coefficient g·β (force = gβ(T − T_ref) ẑ).
    pub gbeta: f64,
    /// Reference temperature.
    pub t_ref: f64,
}

impl ThermalLbm3 {
    /// A quiescent conducting slab: `nx × ny × nz_fluid` fluid cells
    /// between the plates, linear conduction profile as the initial
    /// state.
    ///
    /// # Panics
    /// On inadmissible dimensions, relaxation times, or buoyancy.
    #[must_use]
    pub fn slab(
        nx: usize,
        ny: usize,
        nz_fluid: usize,
        tau_f: f64,
        tau_g: f64,
        gbeta: f64,
    ) -> ThermalLbm3 {
        assert!(
            tau_g.is_finite() && tau_g > 0.5,
            "thermal relaxation time tau_g must be finite and greater than 0.5"
        );
        assert!(gbeta.is_finite(), "buoyancy coefficient must be finite");
        let grid = PlatesGrid3::plates(nx, ny, nz_fluid, tau_f);
        let nz = grid.nz;
        let (t_bottom, t_top) = (1.0f64, 0.0f64);
        let mut gpop = vec![[0.0f64; 7]; nx * ny * nz];
        for y in 0..ny {
            for x in 0..nx {
                gpop[grid.idx(x, y, 0)] = geq3(t_bottom, [0.0; 3]);
                gpop[grid.idx(x, y, nz - 1)] = geq3(t_top, [0.0; 3]);
            }
        }
        for z in 1..nz - 1 {
            // Linear conduction profile between the halfway plate planes.
            let t = t_bottom + (t_top - t_bottom) * ((z as f64 - 0.5) / nz_fluid as f64);
            for y in 0..ny {
                for x in 0..nx {
                    gpop[grid.idx(x, y, z)] = geq3(t, [0.0; 3]);
                }
            }
        }
        ThermalLbm3 {
            grid,
            gpop,
            tau_g,
            t_bottom,
            t_top,
            gbeta,
            t_ref: 0.5,
        }
    }

    /// Temperature of cell (x, y, z).
    #[must_use]
    pub fn temperature(&self, x: usize, y: usize, z: usize) -> f64 {
        self.gpop[self.grid.idx(x, y, z)].iter().sum()
    }

    /// Thermal diffusivity α = c_s²(τ_g − ½) with the D3Q7 c_s² = ¼.
    #[must_use]
    pub fn diffusivity(&self) -> f64 {
        CS2_T7 * (self.tau_g - 0.5)
    }

    /// Seed a sinusoidal vertical-velocity roll (onset mode) spanning x,
    /// preserving density, transverse velocity and force-corrected viscous stress.
    pub fn perturb(&mut self, amplitude: f64) {
        assert!(
            amplitude.is_finite(),
            "perturbation amplitude must be finite"
        );
        if amplitude == 0.0 {
            return;
        }
        let (nx, ny, nz) = (self.grid.nx, self.grid.ny, self.grid.nz);
        for z in 1..nz - 1 {
            for y in 0..ny {
                for x in 0..nx {
                    let i = self.grid.idx(x, y, z);
                    let s = fs_math::det::sin(std::f64::consts::TAU * x as f64 / nx as f64)
                        * fs_math::det::sin(
                            std::f64::consts::PI * (z as f64 - 0.5) / (nz as f64 - 2.0),
                        );
                    let (rho, u) = self.grid.moments3(i);
                    let next = [u[0], u[1], amplitude.mul_add(s, u[2])];
                    let before = equilibrium3(rho, u);
                    let after = equilibrium3(rho, next);
                    let force = self.grid.fext[i];
                    crate::thermal::forcing::reframe3(&mut self.grid.f[i], u, force, [0.0; 3]);
                    for q in 0..Q3 {
                        self.grid.f[i][q] += after[q] - before[q];
                    }
                    crate::thermal::forcing::reframe3(&mut self.grid.f[i], next, [0.0; 3], force);
                }
            }
        }
    }

    /// One coupled step: Boussinesq force from T, flow step, then
    /// temperature collide + stream (anti-bounce-back at the fixed-T
    /// plates, periodic x/y).
    pub fn step(&mut self) {
        self.step_with_forces_inner(None, None);
    }

    /// Advance with external lattice force densities added to thermal buoyancy.
    /// The borrowed field has one vector per cell including the two wall
    /// planes, whose loads are ignored. It contains EXTERNAL forces only, not
    /// the previous step's total `fext`. Admission failures leave state intact;
    /// ordinary solver-state and stability requirements still apply.
    /// Total forces are held over this interval. Source-frame changes preserve
    /// physical velocity and force-corrected stress before the actual step, so
    /// changing a magnetic load cannot add a second half-force impulse.
    pub fn step_with_external_forces(
        &mut self,
        forces: &[[f64; 3]],
    ) -> Result<(), ThermalForcingError> {
        let plane = self
            .grid
            .nx
            .checked_mul(self.grid.ny)
            .ok_or(ThermalForcingError::Layout)?;
        validate_external_forces(forces, &self.gpop, plane, self.gbeta, self.t_ref)?;
        self.step_with_forces_inner(Some(forces), None);
        Ok(())
    }

    // The shared transactional source owner admits the complete rate field.
    pub(crate) fn step_with_local_thermal_relaxation(
        &mut self,
        forces: &[[f64; 3]],
        thermal_tau: &[f64],
    ) -> Result<(), ThermalForcingError> {
        let plane = self
            .grid
            .nx
            .checked_mul(self.grid.ny)
            .ok_or(ThermalForcingError::Layout)?;
        validate_external_forces(forces, &self.gpop, plane, self.gbeta, self.t_ref)?;
        self.step_with_forces_inner(Some(forces), Some(thermal_tau));
        Ok(())
    }

    fn step_with_forces_inner(
        &mut self,
        external: Option<&[[f64; 3]]>,
        thermal_tau: Option<&[f64]>,
    ) {
        let (nx, ny, nz) = (self.grid.nx, self.grid.ny, self.grid.nz);
        for z in 1..nz - 1 {
            for y in 0..ny {
                for x in 0..nx {
                    let i = self.grid.idx(x, y, z);
                    let t: f64 = self.gpop[i].iter().sum();
                    let mut total = [0.0, 0.0, self.gbeta * (t - self.t_ref)];
                    if let Some(forces) = external {
                        for (total, &load) in total.iter_mut().zip(&forces[i]) {
                            *total += load;
                        }
                    }
                    let previous = self.grid.fext[i];
                    if previous != total {
                        let velocity = self.grid.moments3(i).1;
                        crate::thermal::forcing::reframe3(&mut self.grid.f[i], velocity, previous, total);
                        self.grid.fext[i] = total;
                    }
                }
            }
        }
        self.grid.step();
        // Flow streaming has consumed grid.post. Its first 7*N scalars now
        // hold tightly packed thermal collision; the next flow collision
        // overwrites them before use. No full-grid allocation per step.
        for z in 1..nz - 1 {
            for y in 0..ny {
                for x in 0..nx {
                    let i = self.grid.idx(x, y, z);
                    let (_, u) = self.grid.moments3(i);
                    let t: f64 = self.gpop[i].iter().sum();
                    let eq = geq3(t, u);
                    let tau_g = thermal_tau.map_or(self.tau_g, |field| field[i]);
                    for q in 0..7 {
                        self.grid.post.as_flattened_mut()[7 * i + q] =
                            self.gpop[i][q] + (eq[q] - self.gpop[i][q]) / tau_g;
                    }
                }
            }
        }
        // Temperature pull-stream, anti-bounce-back at the plates.
        for z in 1..nz - 1 {
            for y in 0..ny {
                for x in 0..nx {
                    let i = self.grid.idx(x, y, z);
                    for q in 0..7 {
                        let (ex, ey, ez) = E7[q];
                        let sx = (x as i64 - i64::from(ex)).rem_euclid(nx as i64) as usize;
                        let sy = (y as i64 - i64::from(ey)).rem_euclid(ny as i64) as usize;
                        let sz = (z as i64 - i64::from(ez)) as usize;
                        self.gpop[i][q] = if self.grid.is_wall_layer(sz) {
                            let tw = if sz == 0 { self.t_bottom } else { self.t_top };
                            // Anti-bounce-back: fixed halfway temperature.
                            2.0f64.mul_add(W7[q] * tw, -self.grid.post.as_flattened()[7 * i + OPP7[q]])
                        } else {
                            self.grid.post.as_flattened()[7 * self.grid.idx(sx, sy, sz) + q]
                        };
                    }
                }
            }
        }
    }

    /// Total kinetic energy of the fluid.
    #[must_use]
    pub fn kinetic_energy(&self) -> f64 {
        let (nx, ny, nz) = (self.grid.nx, self.grid.ny, self.grid.nz);
        let mut ke = 0.0f64;
        for z in 1..nz - 1 {
            for y in 0..ny {
                for x in 0..nx {
                    let (rho, u) = self.grid.moments3(self.grid.idx(x, y, z));
                    ke += 0.5 * rho * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
                }
            }
        }
        ke
    }

    /// Nusselt number: 1 + ⟨u_z·T⟩·H / (α·ΔT).
    #[must_use]
    pub fn nusselt(&self) -> f64 {
        let (nx, ny, nz) = (self.grid.nx, self.grid.ny, self.grid.nz);
        let h = (nz - 2) as f64;
        let mut adv = 0.0f64;
        let mut count = 0usize;
        for z in 1..nz - 1 {
            for y in 0..ny {
                for x in 0..nx {
                    let i = self.grid.idx(x, y, z);
                    let t: f64 = self.gpop[i].iter().sum();
                    adv += self.grid.moments3(i).1[2] * t;
                    count += 1;
                }
            }
        }
        adv /= count as f64;
        1.0 + adv * h / (self.diffusivity() * (self.t_bottom - self.t_top))
    }

    /// The Rayleigh number of the current configuration.
    #[must_use]
    pub fn rayleigh(&self) -> f64 {
        let h = (self.grid.nz - 2) as f64;
        let i = self.grid.idx(0, 0, 1);
        let nu = CS2_F * (self.grid.tau[i] - 0.5);
        self.gbeta * (self.t_bottom - self.t_top) * h * h * h / (nu * self.diffusivity())
    }
}

/// gβ needed for a target Rayleigh number at the given 3-D lattice setup
/// (ν from the D3Q19 flow lattice, α from the D3Q7 temperature lattice).
#[must_use]
pub fn gbeta_for_rayleigh3(ra: f64, nz_fluid: usize, tau_f: f64, tau_g: f64) -> f64 {
    assert!(ra.is_finite(), "Rayleigh number must be finite");
    assert!(nz_fluid > 0, "Rayleigh height must be positive");
    assert!(
        tau_f.is_finite() && tau_f > 0.5,
        "flow relaxation time tau_f must be finite and greater than 0.5"
    );
    assert!(
        tau_g.is_finite() && tau_g > 0.5,
        "thermal relaxation time tau_g must be finite and greater than 0.5"
    );
    let h = nz_fluid as f64;
    let nu = CS2_F * (tau_f - 0.5);
    let alpha = CS2_T7 * (tau_g - 0.5);
    ra * nu * alpha / (h * h * h)
}
