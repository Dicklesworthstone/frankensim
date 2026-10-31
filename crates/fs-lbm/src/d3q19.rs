//! The D3Q19 core (bead 84hv): 3-D BGK stream-and-collide with Guo body
//! forcing and halfway bounce-back walls, over 128-byte-aligned SoA
//! distributions in tile-major layout (plan §5) with the PULL scheme —
//! the 3-D sibling of the D2Q9 module in `lib.rs`, which it deliberately
//! mirrors (and leaves untouched).
//!
//! DETERMINISTIC TRAVERSAL (pinned): collision visits tiles in ascending
//! tile index and cells in ascending local index (x-fastest, then y, then
//! z). Duct pull-streaming visits directions, z rows, y rows, then x tiles,
//! moving each four-cell x row as one block. Single-threaded in this bead;
//! when WS1-D parallelizes the sweep, per-cell writes stay slot-exclusive
//! (collision is pointwise; pull-streaming writes only the destination
//! population), so these orders document CANONICAL schedules rather than
//! hidden result dependencies — results must stay bit-identical across
//! worker counts by construction.
//!
//! Layout: the domain is split into 4×4×4 tiles (64 cells, 512 B per
//! field per tile, 128-byte aligned). Each of the 19 distribution
//! fields is its own `Vec<Tile>` (structure of arrays); cell `(x,y,z)`
//! lives in tile `(x/4, y/4, z/4)` at local `(x%4, y%4, z%4)`.
//! Dimensions must be multiples of 4 (asserted): the fixture scales in
//! WS1 are, and padding is a later concern, not silent slop.

mod boundary;
mod coupled;
mod fma;
pub mod freesurface3;
mod simd;
pub mod sparse;

pub use boundary::{
    BodyFrameEllipsoid3, BodyFrameSampler3, BodyFrameShape3, BoundaryGrid3, BoundaryLink3,
    BoundaryProfileError3, BoundarySgsError3, BoundarySgsReport3, BoundarySpec3,
    BoundaryStepError3, D3Q19_BOUNDARY_BIT_SEMANTICS_VERSION, Face3, FaceBoundary3,
    FiniteMassWall3, FiniteMassWallReceipt3, InletFlux3, LinkMaskTile3, LubricatedPlanesExchange3,
    LubricationPlane3, MovingLubricatedPlanesExchange3, MovingPartialSaturationBodiesExchange3,
    MovingPartialSaturationBody3, MovingPartialSaturationExchange3,
    PartialSaturationBodiesExchange3, PartialSaturationBody3, PartialSaturationBodyExchange3,
    PartialSaturationBodyResponse3, PartialSaturationError3, PartialSaturationExchange3,
    PartialSaturationForcing3, PeriodicMomentumTransport3, PlaneLubricationExchange3,
    SphereLubrication3, SphereLubricationExchange3, VelocityProfileReport3, WallMomentumExchange3,
};
pub(crate) use coupled::force_corrected_shear_rate3;
pub use coupled::{
    PlatesGrid3, ThermalLbm3, gbeta_for_rayleigh3, plate_channel_flow3, shear_rate3, update_tau3,
};
pub use simd::{
    D3q19BgkSimdTier, D3q19StreamSimdTier, d3q19_bgk_simd_tier, d3q19_stream_simd_tier,
};

/// Bit-semantics version of the D3Q19 surface (golden-couplings.json):
/// covers the velocity/weight/opposite tables and ordering, the
/// equilibrium form, the Guo forcing form, the pull-stream + halfway
/// bounce-back rules, and the pinned traversal order. Bump on ANY
/// change that can move result bits.
pub const D3Q19_BIT_SEMANTICS_VERSION: u32 = 1;

/// Bit-semantics version for the optional moment-space collision surface.
///
/// This is intentionally independent of [`D3Q19_BIT_SEMANTICS_VERSION`]: the
/// frozen BGK grids keep their established golden, while changes to the
/// centered basis, relaxation grouping, cumulant projection, `mul_add` chains,
/// or deterministic solve order bump this version and require their own
/// golden evidence. Version 3 replaces the per-cell dense back-transform with
/// an algebraic central-to-raw shift and fixed D3Q19 inverse. This changes
/// moment-collision roundoff, including the unforced path; BGK is unchanged.
pub const D3Q19_MOMENT_COLLISION_SEMANTICS_VERSION: u32 = 3;

/// The D3Q19 population count.
pub const Q3: usize = 19;

/// The D3Q19 lattice velocities: rest, 6 face neighbors, 12 edge
/// neighbors — opposites adjacent (`2k-1` ↔ `2k`), so the opposite
/// table is verifiable at a glance.
pub const E3: [(i32, i32, i32); Q3] = [
    (0, 0, 0),
    (1, 0, 0),
    (-1, 0, 0),
    (0, 1, 0),
    (0, -1, 0),
    (0, 0, 1),
    (0, 0, -1),
    (1, 1, 0),
    (-1, -1, 0),
    (1, -1, 0),
    (-1, 1, 0),
    (1, 0, 1),
    (-1, 0, -1),
    (1, 0, -1),
    (-1, 0, 1),
    (0, 1, 1),
    (0, -1, -1),
    (0, 1, -1),
    (0, -1, 1),
];

/// The D3Q19 weights: 1/3 rest, 1/18 per face, 1/36 per edge.
pub const W3: [f64; Q3] = [
    1.0 / 3.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// Integer weights ×36 — the EXACT arithmetic the lattice-invariant
/// tests use (Σ = 36, moments in integers, no float tolerance).
pub const W36: [i64; Q3] = [12, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

/// Opposite-direction indices (for bounce-back): rest is self-opposite,
/// then adjacent pairs.
pub const OPP3: [usize; Q3] = [
    0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17,
];

/// Tile edge in cells (4×4×4 = 64 cells per tile).
pub const TILE: usize = 4;
const TILE_CELLS: usize = TILE * TILE * TILE;

/// One 4×4×4 tile of one scalar field: 512 B, 128-byte aligned — the
/// plan §5 SoA tile-major unit.
#[derive(Clone)]
#[repr(align(128))]
struct Tile([f64; TILE_CELLS]);

impl Tile {
    fn filled(value: f64) -> Tile {
        Tile([value; TILE_CELLS])
    }
}

/// The D3Q19 equilibrium distribution at density `rho` and velocity `u`.
#[must_use]
pub fn equilibrium3(rho: f64, u: [f64; 3]) -> [f64; Q3] {
    let usq = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    let mut f = [0.0; Q3];
    for i in 0..Q3 {
        let (ex, ey, ez) = (f64::from(E3[i].0), f64::from(E3[i].1), f64::from(E3[i].2));
        let eu = ex * u[0] + ey * u[1] + ez * u[2];
        f[i] = W3[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * usq);
    }
    f
}

/// Full-rank monomial basis on D3Q19, ordered by total degree.
///
/// Every D3Q19 velocity has at least one zero component, so monomials that
/// contain all three axes vanish on the lattice. The remaining tensor-product
/// monomials through exponent two give exactly 19 independent rows. Centering
/// these rows at the local velocity is a triangular change of basis and
/// therefore preserves full rank.
const CENTRAL_MOMENT_EXPONENTS3: [[u8; 3]; Q3] = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    [2, 0, 0],
    [0, 2, 0],
    [0, 0, 2],
    [1, 1, 0],
    [1, 0, 1],
    [0, 1, 1],
    [2, 1, 0],
    [2, 0, 1],
    [1, 2, 0],
    [0, 2, 1],
    [1, 0, 2],
    [0, 1, 2],
    [2, 2, 0],
    [2, 0, 2],
    [0, 2, 2],
];

/// Selectable D3Q19 collision law for the shared per-cell kernel.
///
/// The frozen BGK law remains the grid default. The central-moment and reduced
/// cumulant rungs are correctness-first reference operators; neither
/// claims a production high-Reynolds-number model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CollisionModel3 {
    /// Single-relaxation-time BGK collision.
    Bgk {
        /// Relaxation time; must be finite and greater than one half.
        tau: f64,
    },
    /// Relax central moments around the measured local velocity.
    ///
    /// Rates must be finite and strictly between zero and two. Rows of total
    /// degree zero and one are conserved; degree-two rows use
    /// `second_order_rate`, and all remaining rows use `higher_order_rate`.
    /// Body forcing uses the physical-population transform `f + S/2` before
    /// relaxation and restores `S/2` afterward, with the original Guo source.
    CentralMoment {
        /// Relaxation rate for the six degree-two moments.
        second_order_rate: f64,
        /// Relaxation rate for the nine degree-three/four moments.
        higher_order_rate: f64,
    },
    /// Reduced cumulant projection onto the 19 independent D3Q19 moments.
    ///
    /// Second- and third-order cumulants equal their central moments. The
    /// three represented fourth-order cumulants (`C220`, `C202`, and `C022`)
    /// use the nonlinear product corrections from Geier et al. This is not the
    /// paper's D3Q27 operator: D3Q19 omits eight velocities and the associated
    /// independent moments, corrections, and isotropy.
    /// Guo forcing uses the same physical-population transform as central
    /// moments, so the cumulants are evaluated in their zero-first-moment frame.
    ReducedCumulant {
        /// Relaxation rate for the six degree-two cumulants.
        second_order_rate: f64,
        /// Relaxation rate for the six represented degree-three cumulants.
        third_order_rate: f64,
        /// Relaxation rate for the three represented degree-four cumulants.
        fourth_order_rate: f64,
    },
}

impl CollisionModel3 {
    /// Validate all declared relaxation parameters.
    ///
    /// # Errors
    /// [`CollisionError3::InvalidRelaxationTime`] or
    /// [`CollisionError3::InvalidMomentRelaxationRate`] when a parameter is
    /// outside its finite physical window.
    pub fn validate(self) -> Result<(), CollisionError3> {
        match self {
            Self::Bgk { tau } if tau.is_finite() && tau > 0.5 => Ok(()),
            Self::Bgk { tau } => Err(CollisionError3::InvalidRelaxationTime { tau }),
            Self::CentralMoment {
                second_order_rate,
                higher_order_rate,
            } => {
                validate_moment_rate(2, second_order_rate)?;
                validate_moment_rate(3, higher_order_rate)
            }
            Self::ReducedCumulant {
                second_order_rate,
                third_order_rate,
                fourth_order_rate,
            } => {
                validate_moment_rate(2, second_order_rate)?;
                validate_moment_rate(3, third_order_rate)?;
                validate_moment_rate(4, fourth_order_rate)
            }
        }
    }

    /// Kinematic viscosity implied by the BGK relaxation time or the
    /// moment-space second-order relaxation rate.
    ///
    /// Call [`CollisionModel3::validate`] first when the model did not come
    /// from a checked constructor.
    #[must_use]
    pub const fn kinematic_viscosity(self) -> f64 {
        match self {
            Self::Bgk { tau } => (tau - 0.5) / 3.0,
            Self::CentralMoment {
                second_order_rate, ..
            }
            | Self::ReducedCumulant {
                second_order_rate, ..
            } => (1.0 / second_order_rate - 0.5) / 3.0,
        }
    }

    /// Whether the current collision contract admits a nonzero body force.
    #[must_use]
    pub const fn supports_body_force(self) -> bool {
        match self {
            Self::Bgk { .. } | Self::CentralMoment { .. } | Self::ReducedCumulant { .. } => true,
        }
    }
}

fn validate_moment_rate(moment_order: u8, rate: f64) -> Result<(), CollisionError3> {
    if rate.is_finite() && rate > 0.0 && rate < 2.0 {
        Ok(())
    } else {
        Err(CollisionError3::InvalidMomentRelaxationRate { moment_order, rate })
    }
}

/// Fail-closed D3Q19 cell-collision diagnostics.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CollisionError3 {
    /// Relaxation time is non-finite or does not yield positive viscosity.
    InvalidRelaxationTime {
        /// Rejected relaxation time.
        tau: f64,
    },
    /// A moment-space relaxation rate is outside `(0, 2)`.
    InvalidMomentRelaxationRate {
        /// Lowest total moment order governed by this rate (`2`, `3`, or `4`).
        moment_order: u8,
        /// Rejected relaxation rate.
        rate: f64,
    },
    /// A body-force component is non-finite.
    NonFiniteForce {
        /// Cartesian axis, `0..3`.
        axis: usize,
        /// Rejected component.
        value: f64,
    },
    /// An incoming population is non-finite.
    NonFinitePopulation {
        /// D3Q19 direction index.
        direction: usize,
        /// Rejected population.
        value: f64,
    },
    /// Incoming populations do not define positive finite mass.
    NonPositiveDensity {
        /// Computed density.
        rho: f64,
    },
    /// A force-corrected velocity component is non-finite.
    NonFiniteVelocity {
        /// Cartesian axis, `0..3`.
        axis: usize,
        /// Computed component.
        value: f64,
    },
    /// Legacy central-moment force refusal, retained for diagnostic compatibility.
    /// Version 2's physical-population forcing no longer emits this variant.
    CentralMomentForceUnsupported {
        /// Rejected body-force vector.
        force: [f64; 3],
    },
    /// Legacy reduced-cumulant force refusal, retained for diagnostic compatibility.
    /// Version 2's physical-population forcing no longer emits this variant.
    ReducedCumulantForceUnsupported {
        /// Rejected body-force vector.
        force: [f64; 3],
    },
    /// The dense reference transform could not be solved reliably.
    /// The version-3 algebraic inverse does not require pivoting.
    SingularCentralMomentTransform {
        /// Pivot column that lost rank.
        column: usize,
        /// Absolute value of the rejected pivot.
        pivot_abs: f64,
    },
    /// Collision produced a non-finite outgoing population.
    NonFiniteOutput {
        /// D3Q19 direction index.
        direction: usize,
        /// Computed population.
        value: f64,
    },
}

impl core::fmt::Display for CollisionError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidRelaxationTime { tau } => {
                write!(f, "D3Q19 relaxation time {tau} must be finite and > 0.5")
            }
            Self::InvalidMomentRelaxationRate { moment_order, rate } => write!(
                f,
                "D3Q19 moment-order family {moment_order} relaxation rate {rate} must be finite and in (0, 2)"
            ),
            Self::NonFiniteForce { axis, value } => {
                write!(f, "D3Q19 force axis {axis} is non-finite ({value})")
            }
            Self::NonFinitePopulation { direction, value } => write!(
                f,
                "D3Q19 incoming population {direction} is non-finite ({value})"
            ),
            Self::NonPositiveDensity { rho } => {
                write!(f, "D3Q19 density must be positive and finite (got {rho})")
            }
            Self::NonFiniteVelocity { axis, value } => {
                write!(f, "D3Q19 velocity axis {axis} is non-finite ({value})")
            }
            Self::CentralMomentForceUnsupported { force } => write!(
                f,
                "D3Q19 central-moment collision does not yet admit body force {force:?}"
            ),
            Self::ReducedCumulantForceUnsupported { force } => write!(
                f,
                "D3Q19 reduced-cumulant collision does not yet admit body force {force:?}"
            ),
            Self::SingularCentralMomentTransform { column, pivot_abs } => write!(
                f,
                "D3Q19 central-moment transform lost rank at column {column} (pivot {pivot_abs})"
            ),
            Self::NonFiniteOutput { direction, value } => write!(
                f,
                "D3Q19 outgoing population {direction} is non-finite ({value})"
            ),
        }
    }
}

impl core::error::Error for CollisionError3 {}

/// Collide one D3Q19 population vector under an explicit checked model and
/// Guo body force.
///
/// The BGK arithmetic deliberately retains the two existing evaluation paths:
/// this public entry point uses the frozen boundary-grid expression, while the
/// [`Duct`] tile kernel retains its separately frozen axial-force expression.
/// A test-only axial cell oracle locks the tile scalar/SIMD twins to that
/// expression without silently refreezing either golden surface.
///
/// On x86-64 CPUs with FMA the same body runs compiled with native fused
/// multiply-add (the `fma` capsule), bit-identical to the portable codegen.
///
/// # Errors
/// [`CollisionError3`] for inadmissible model parameters, force, incoming
/// state, macroscopic state, or outgoing populations.
pub fn collide_cell3(
    populations: [f64; Q3],
    model: CollisionModel3,
    force: [f64; 3],
) -> Result<[f64; Q3], CollisionError3> {
    fma::collide_cell3(populations, model, force)
}

// Raw Guo source projection, before its weight and collision prefactor. The
// free-surface force-frame transform uses this SAME polynomial, so changing a
// source does not change physical momentum or the force-corrected stress.
pub(crate) fn guo_projection3(e: [f64; 3], velocity: [f64; 3], eu: f64, force: [f64; 3]) -> f64 {
    let cs2 = 1.0 / 3.0;
    let cs4 = cs2 * cs2;
    (0..3)
        .map(|axis| ((e[axis] - velocity[axis]) / cs2 + eu * e[axis] / cs4) * force[axis])
        .sum()
}

#[cfg(test)]
fn collide_axial_z_cell3(
    populations: [f64; Q3],
    model: CollisionModel3,
    gz: f64,
) -> Result<[f64; Q3], CollisionError3> {
    collide_cell3_with_projection(populations, model, [0.0, 0.0, gz], true)
}

// This body and every moment-space helper below that uses `mul_add` MUST stay
// `inline(always)`: the `fma` capsule recompiles them by inlining, and a
// non-inlined call keeps baseline codegen (one libm `fma` call per mul_add).
#[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
#[inline(always)]
fn collide_cell3_with_projection(
    populations: [f64; Q3],
    model: CollisionModel3,
    force: [f64; 3],
    frozen_axial_z_projection: bool,
) -> Result<[f64; Q3], CollisionError3> {
    model.validate()?;
    for (axis, value) in force.into_iter().enumerate() {
        if !value.is_finite() {
            return Err(CollisionError3::NonFiniteForce { axis, value });
        }
    }

    let mut rho = 0.0;
    let mut momentum = [0.0; 3];
    for (direction, value) in populations.into_iter().enumerate() {
        if !value.is_finite() {
            return Err(CollisionError3::NonFinitePopulation { direction, value });
        }
        rho += value;
        momentum[0] += f64::from(E3[direction].0) * value;
        momentum[1] += f64::from(E3[direction].1) * value;
        momentum[2] += f64::from(E3[direction].2) * value;
    }
    if !(rho.is_finite() && rho > 0.0) {
        return Err(CollisionError3::NonPositiveDensity { rho });
    }
    let velocity: [f64; 3] =
        core::array::from_fn(|axis| (momentum[axis] + 0.5 * force[axis]) / rho);
    for (axis, value) in velocity.into_iter().enumerate() {
        if !value.is_finite() {
            return Err(CollisionError3::NonFiniteVelocity { axis, value });
        }
    }

    let equilibrium = equilibrium3(rho, velocity);
    match model {
        CollisionModel3::Bgk { tau } => {
            let coefficient = 1.0 - 0.5 / tau;
            let mut post = [0.0; Q3];
            for direction in 0..Q3 {
                let e = [
                    f64::from(E3[direction].0),
                    f64::from(E3[direction].1),
                    f64::from(E3[direction].2),
                ];
                let eu = if frozen_axial_z_projection {
                    e[0] * velocity[0] + e[1] * velocity[1] + e[2] * velocity[2]
                } else {
                    e.iter()
                        .zip(velocity)
                        .map(|(component, u)| *component * u)
                        .sum::<f64>()
                };
                let forcing = if frozen_axial_z_projection {
                    coefficient
                        * W3[direction]
                        * (3.0 * (e[2] - velocity[2]) + 9.0 * eu * e[2])
                        * force[2]
                } else {
                    let projection = guo_projection3(e, velocity, eu, force);
                    coefficient * W3[direction] * projection
                };
                let value = populations[direction]
                    + (equilibrium[direction] - populations[direction]) / tau
                    + forcing;
                if !value.is_finite() {
                    return Err(CollisionError3::NonFiniteOutput { direction, value });
                }
                post[direction] = value;
            }
            Ok(post)
        }
        CollisionModel3::CentralMoment { .. } | CollisionModel3::ReducedCumulant { .. } => {
            collide_moment_with_force3(populations, equilibrium, velocity, force, model)
        }
    }
}

/// Symmetric source transform around the existing collision operator.
/// For the linear central-moment operator this is `(I - Lambda/2) S`
/// in every represented mode (Premnath et al., 2009, Eq. (1),
/// https://arxiv.org/abs/0901.2786). For reduced cumulants the nonlinear
/// projection acts on physical populations, whose first central moments
/// vanish; applying it to the stored `-F/2` momentum would violate that
/// projection's assumptions. No collision rate or source is fitted here.
#[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
#[inline(always)]
fn collide_moment_with_force3(
    populations: [f64; Q3],
    equilibrium: [f64; Q3],
    velocity: [f64; 3],
    force: [f64; 3],
    model: CollisionModel3,
) -> Result<[f64; Q3], CollisionError3> {
    if force.iter().all(|&component| component == 0.0) {
        // Avoid source arithmetic in unforced cells.
        return collide_moments3(populations, equilibrium, velocity, model);
    }
    let half_source: [f64; Q3] = core::array::from_fn(|q| {
        let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
        let eu = e.iter().zip(velocity).map(|(e, u)| e * u).sum();
        0.5 * W3[q] * guo_projection3(e, velocity, eu, force)
    });
    let physical = core::array::from_fn(|q| populations[q] + half_source[q]);
    let mut post = collide_moments3(physical, equilibrium, velocity, model)?;
    for q in 0..Q3 {
        post[q] += half_source[q];
        if !post[q].is_finite() {
            return Err(CollisionError3::NonFiniteOutput {
                direction: q,
                value: post[q],
            });
        }
    }
    Ok(post)
}

/// The moment-space collision of physical populations for a moment model.
/// A function rather than a closure, so it inlines into the `fma` capsule.
#[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
#[inline(always)]
fn collide_moments3(
    physical: [f64; Q3],
    equilibrium: [f64; Q3],
    velocity: [f64; 3],
    model: CollisionModel3,
) -> Result<[f64; Q3], CollisionError3> {
    match model {
        CollisionModel3::CentralMoment {
            second_order_rate,
            higher_order_rate,
        } => collide_central_moments3(
            physical,
            equilibrium,
            velocity,
            second_order_rate,
            higher_order_rate,
        ),
        CollisionModel3::ReducedCumulant {
            second_order_rate,
            third_order_rate,
            fourth_order_rate,
        } => collide_reduced_cumulants3(
            physical,
            equilibrium,
            velocity,
            second_order_rate,
            third_order_rate,
            fourth_order_rate,
        ),
        CollisionModel3::Bgk { .. } => {
            unreachable!("BGK relaxes populations directly, never through moments")
        }
    }
}

fn centered_power3(value: f64, exponent: u8) -> f64 {
    match exponent {
        0 => 1.0,
        1 => value,
        2 => value * value,
        _ => unreachable!("D3Q19 central-moment exponents are at most two"),
    }
}

fn central_moment_matrix3(velocity: [f64; 3]) -> [[f64; Q3]; Q3] {
    core::array::from_fn(|moment| {
        let exponent = CENTRAL_MOMENT_EXPONENTS3[moment];
        core::array::from_fn(|direction| {
            let centered = [
                f64::from(E3[direction].0) - velocity[0],
                f64::from(E3[direction].1) - velocity[1],
                f64::from(E3[direction].2) - velocity[2],
            ];
            centered_power3(centered[0], exponent[0])
                * centered_power3(centered[1], exponent[1])
                * centered_power3(centered[2], exponent[2])
        })
    })
}

// Plain loops, not an `array::from_fn` closure, so the whole transform
// inlines into the `fma` capsule; each moment accumulates in direction order.
#[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
#[inline(always)]
fn apply_moment_matrix3(matrix: &[[f64; Q3]; Q3], populations: [f64; Q3]) -> [f64; Q3] {
    let mut moments = [0.0; Q3];
    for (moment, row) in moments.iter_mut().zip(matrix) {
        let mut value = 0.0;
        for (coefficient, population) in row.iter().zip(populations) {
            value = coefficient.mul_add(population, value);
        }
        *moment = value;
    }
    moments
}

#[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
#[inline(always)]
fn collide_central_moments3(
    populations: [f64; Q3],
    equilibrium: [f64; Q3],
    velocity: [f64; 3],
    second_order_rate: f64,
    higher_order_rate: f64,
) -> Result<[f64; Q3], CollisionError3> {
    let matrix = central_moment_matrix3(velocity);
    let moments = apply_moment_matrix3(&matrix, populations);
    let equilibrium_moments = apply_moment_matrix3(&matrix, equilibrium);
    let mut relaxed = moments;
    for moment in 4..Q3 {
        let rate = if moment < 10 {
            second_order_rate
        } else {
            higher_order_rate
        };
        relaxed[moment] = rate.mul_add(
            equilibrium_moments[moment] - moments[moment],
            moments[moment],
        );
    }
    let post = populations_from_central_moments3(relaxed, velocity);
    for (direction, value) in post.into_iter().enumerate() {
        if !value.is_finite() {
            return Err(CollisionError3::NonFiniteOutput { direction, value });
        }
    }
    Ok(post)
}

#[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
#[inline(always)]
fn fourth_order_pair_product3(first: f64, second: f64, cross: f64) -> f64 {
    first.mul_add(second, 2.0 * cross * cross)
}

// Geier et al. (2015), Eqs. 46-52, use unnormalized cumulants C = rho*c.
// On this reduced D3Q19 basis the only independent order-four instances are
// C220, C202, and C022; each subtracts its paired normal stresses and twice
// the squared cross stress divided by density.
#[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
#[inline(always)]
fn reduced_cumulants_from_moments3(moments: [f64; Q3]) -> [f64; Q3] {
    let inverse_density = 1.0 / moments[0];
    let mut cumulants = moments;
    let xy = fourth_order_pair_product3(moments[4], moments[5], moments[7]);
    let xz = fourth_order_pair_product3(moments[4], moments[6], moments[8]);
    let yz = fourth_order_pair_product3(moments[5], moments[6], moments[9]);
    cumulants[16] = (-xy).mul_add(inverse_density, moments[16]);
    cumulants[17] = (-xz).mul_add(inverse_density, moments[17]);
    cumulants[18] = (-yz).mul_add(inverse_density, moments[18]);
    cumulants
}

#[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
#[inline(always)]
fn reduced_moments_from_cumulants3(cumulants: [f64; Q3]) -> [f64; Q3] {
    let inverse_density = 1.0 / cumulants[0];
    let mut moments = cumulants;
    let xy = fourth_order_pair_product3(cumulants[4], cumulants[5], cumulants[7]);
    let xz = fourth_order_pair_product3(cumulants[4], cumulants[6], cumulants[8]);
    let yz = fourth_order_pair_product3(cumulants[5], cumulants[6], cumulants[9]);
    moments[16] = xy.mul_add(inverse_density, cumulants[16]);
    moments[17] = xz.mul_add(inverse_density, cumulants[17]);
    moments[18] = yz.mul_add(inverse_density, cumulants[18]);
    moments
}

#[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
#[inline(always)]
fn collide_reduced_cumulants3(
    populations: [f64; Q3],
    equilibrium: [f64; Q3],
    velocity: [f64; 3],
    second_order_rate: f64,
    third_order_rate: f64,
    fourth_order_rate: f64,
) -> Result<[f64; Q3], CollisionError3> {
    let matrix = central_moment_matrix3(velocity);
    let moments = apply_moment_matrix3(&matrix, populations);
    let equilibrium_moments = apply_moment_matrix3(&matrix, equilibrium);
    let cumulants = reduced_cumulants_from_moments3(moments);
    let equilibrium_cumulants = reduced_cumulants_from_moments3(equilibrium_moments);
    let mut relaxed = cumulants;
    for moment in 4..Q3 {
        let rate = if moment < 10 {
            second_order_rate
        } else if moment < 16 {
            third_order_rate
        } else {
            fourth_order_rate
        };
        relaxed[moment] = rate.mul_add(
            equilibrium_cumulants[moment] - cumulants[moment],
            cumulants[moment],
        );
    }
    let relaxed_moments = reduced_moments_from_cumulants3(relaxed);
    let post = populations_from_central_moments3(relaxed_moments, velocity);
    for (direction, value) in post.into_iter().enumerate() {
        if !value.is_finite() {
            return Err(CollisionError3::NonFiniteOutput { direction, value });
        }
    }
    Ok(post)
}

// Lowering any nonzero exponent stays in the D3Q19 basis. Build these indices
// at compile time; no moment lookup or matrix factorization occurs per cell.
const LOWER_MOMENT3: [[usize; 3]; Q3] = {
    let mut lower = [[0; 3]; Q3];
    let mut moment = 0;
    while moment < Q3 {
        let mut axis = 0;
        while axis < 3 {
            let mut exponent = CENTRAL_MOMENT_EXPONENTS3[moment];
            if exponent[axis] > 0 {
                exponent[axis] -= 1;
                let mut index = 0;
                while index < moment
                    && (CENTRAL_MOMENT_EXPONENTS3[index][0] != exponent[0]
                        || CENTRAL_MOMENT_EXPONENTS3[index][1] != exponent[1]
                        || CENTRAL_MOMENT_EXPONENTS3[index][2] != exponent[2])
                {
                    index += 1;
                }
                assert!(
                    index < moment,
                    "D3Q19 basis must contain every lower moment"
                );
                lower[moment][axis] = index;
            }
            axis += 1;
        }
        moment += 1;
    }
    lower
};

/// Invert the central-moment transform with O(Q) work and O(Q) scratch.
/// First expand `c = (c-u) + u` by the binomial theorem, one axis at a time.
/// Descending total degree keeps each lower moment unmodified within a shift.
/// The raw D3Q19 inverse then follows from the three four-velocity coordinate
/// planes, six axial velocities and rest population; its coefficients are
/// exactly 1, 1/2 and 1/4. No velocity-dependent linear solve is needed.
#[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
#[inline(always)]
fn populations_from_central_moments3(mut raw: [f64; Q3], velocity: [f64; 3]) -> [f64; Q3] {
    for (axis, shift) in velocity.into_iter().enumerate() {
        let twice = 2.0 * shift;
        let square = shift * shift;
        for moment in (1..Q3).rev() {
            let lower = LOWER_MOMENT3[moment][axis];
            raw[moment] = match CENTRAL_MOMENT_EXPONENTS3[moment][axis] {
                0 => raw[moment],
                1 => shift.mul_add(raw[lower], raw[moment]),
                2 => square.mul_add(
                    raw[LOWER_MOMENT3[lower][axis]],
                    twice.mul_add(raw[lower], raw[moment]),
                ),
                _ => unreachable!("D3Q19 exponents are at most two"),
            };
        }
    }
    let mut populations = [0.0; Q3];
    populations[0] = raw[0] - raw[4] - raw[5] - raw[6] + raw[16] + raw[17] + raw[18];
    for (positive, normal, mixed, first, odd) in [
        (1, 4, [16, 17], 1, [12, 14]),
        (3, 5, [16, 18], 2, [10, 15]),
        (5, 6, [17, 18], 3, [11, 13]),
    ] {
        let even = raw[normal] - raw[mixed[0]] - raw[mixed[1]];
        let odd = raw[first] - raw[odd[0]] - raw[odd[1]];
        populations[positive] = 0.5 * (even + odd);
        populations[positive + 1] = 0.5 * (even - odd);
    }
    for (first, even, odd1, odd2, cross) in
        [(7, 16, 12, 10, 7), (11, 17, 14, 11, 8), (15, 18, 15, 13, 9)]
    {
        let same = raw[even] + raw[cross];
        let opposite = raw[even] - raw[cross];
        let sum = raw[odd1] + raw[odd2];
        let difference = raw[odd1] - raw[odd2];
        populations[first] = 0.25 * (same + sum);
        populations[first + 1] = 0.25 * (same - sum);
        populations[first + 2] = 0.25 * (opposite + difference);
        populations[first + 3] = 0.25 * (opposite - difference);
    }
    populations
}

// Keep the original pivoted solve as an independent test oracle. Production
// moment collisions use the fixed lattice inverse above.
#[cfg(test)]
fn solve_moment_system3(
    mut matrix: [[f64; Q3]; Q3],
    mut rhs: [f64; Q3],
) -> Result<[f64; Q3], CollisionError3> {
    for column in 0..Q3 {
        let mut pivot_row = column;
        let mut pivot_abs = matrix[column][column].abs();
        for (row, coefficients) in matrix.iter().enumerate().skip(column + 1) {
            let candidate = coefficients[column].abs();
            if candidate > pivot_abs {
                pivot_row = row;
                pivot_abs = candidate;
            }
        }
        if !pivot_abs.is_finite() || pivot_abs <= 256.0 * f64::EPSILON {
            return Err(CollisionError3::SingularCentralMomentTransform { column, pivot_abs });
        }
        if pivot_row != column {
            matrix.swap(column, pivot_row);
            rhs.swap(column, pivot_row);
        }

        let pivot_coefficients = matrix[column];
        let pivot = pivot_coefficients[column];
        let pivot_rhs = rhs[column];
        for (row, coefficients) in matrix.iter_mut().enumerate().skip(column + 1) {
            let factor = coefficients[column] / pivot;
            coefficients[column] = 0.0;
            for (coefficient, value) in coefficients.iter_mut().enumerate().skip(column + 1) {
                *value = (-factor).mul_add(pivot_coefficients[coefficient], *value);
            }
            rhs[row] = (-factor).mul_add(pivot_rhs, rhs[row]);
        }
    }

    let mut solution = [0.0; Q3];
    for row in (0..Q3).rev() {
        let mut value = rhs[row];
        for (column, solved) in solution.iter().enumerate().skip(row + 1) {
            value = (-matrix[row][column]).mul_add(*solved, value);
        }
        solution[row] = value / matrix[row][row];
    }
    Ok(solution)
}

#[cfg(test)]
mod forcing_tests {
    use super::*;

    #[test]
    fn algebraic_moment_inverse_matches_dense_reference() {
        // G0: every population basis vector exercises the complete inverse,
        // including moments that are absent from equilibrium/shear fixtures.
        for velocity in [
            [0.0; 3],
            [0.07, -0.11, 0.05],
            [-0.13, 0.09, -0.02],
            [0.1, 0.08, -0.07],
        ] {
            let matrix = central_moment_matrix3(velocity);
            for direction in 0..Q3 {
                let mut population = [0.0; Q3];
                population[direction] = 1.0;
                let moments = apply_moment_matrix3(&matrix, population);
                let actual = populations_from_central_moments3(moments, velocity);
                let reference = solve_moment_system3(matrix, moments).unwrap();
                for q in 0..Q3 {
                    assert!(
                        (actual[q] - population[q]).abs() < 2e-13,
                        "u={velocity:?}, basis={direction}, q={q}"
                    );
                    assert!(
                        (actual[q] - reference[q]).abs() < 2e-13,
                        "dense inverse: u={velocity:?}, basis={direction}, q={q}"
                    );
                }
            }
            // An arbitrary moment vector also checks mixed signs and scales,
            // independent of any equilibrium or collision relaxation law.
            let moments = core::array::from_fn(|i| {
                if i == 0 {
                    1.03
                } else {
                    (i as f64 - 8.0) * 0.017
                }
            });
            let actual = populations_from_central_moments3(moments, velocity);
            let reference = solve_moment_system3(matrix, moments).unwrap();
            for q in 0..Q3 {
                assert!(
                    (actual[q] - reference[q]).abs() < 2e-13,
                    "arbitrary moments: u={velocity:?}, q={q}"
                );
            }
        }
    }

    const MODELS: [CollisionModel3; 2] = [
        CollisionModel3::CentralMoment {
            second_order_rate: 1.23,
            higher_order_rate: 1.71,
        },
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.23,
            third_order_rate: 1.35,
            fourth_order_rate: 1.72,
        },
    ];

    fn source(u: [f64; 3], force: [f64; 3]) -> [f64; Q3] {
        let uf: f64 = u.iter().zip(force).map(|(u, f)| u * f).sum();
        core::array::from_fn(|q| {
            let c = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
            let cu: f64 = c.iter().zip(u).map(|(c, u)| c * u).sum();
            let cf: f64 = c.iter().zip(force).map(|(c, f)| c * f).sum();
            W3[q] * (3.0 * cf + 9.0 * cu * cf - 3.0 * uf)
        })
    }

    #[test]
    fn forced_moment_modes_follow_exact_uniform_acceleration() {
        // G0/G3: with no gradients the quadratic equilibrium evolves exactly
        // under trapezoidal forcing. Shear/bulk relaxation must not alter
        // momentum delivery or create stress over repeated accelerated steps.
        let rho = 1.07;
        let initial_u = [0.03, -0.02, 0.01];
        let force = [3e-5, -4e-5, 2e-5];
        let eq = equilibrium3(rho, initial_u);
        let initial_source = source(initial_u, force);
        for model in MODELS {
            assert!(model.supports_body_force());
            let mut f = core::array::from_fn(|q| eq[q] - 0.5 * initial_source[q]);
            for step in 1..=32 {
                f = collide_cell3(f, model, force).expect("supported physical-frame forcing");
                let velocity =
                    core::array::from_fn(|d| initial_u[d] + step as f64 * force[d] / rho);
                let expected_eq = equilibrium3(rho, velocity);
                let expected_source = source(velocity, force);
                for q in 0..Q3 {
                    let expected = expected_eq[q] - 0.5 * expected_source[q];
                    assert!(
                        (f[q] - expected).abs() < 1e-12,
                        "{model:?}, step {step}, direction {q}: {} vs {expected}",
                        f[q]
                    );
                }
            }
        }
    }

    #[test]
    fn forced_equal_rate_central_moments_recover_bgk() {
        // G0: independent production BGK is an oracle when every central
        // relaxation rate equals 1/tau, including non-equilibrium ghost modes.
        let mut populations = equilibrium3(0.97, [-0.025, 0.04, 0.013]);
        for (q, population) in populations.iter_mut().enumerate() {
            *population += (q as f64 - 9.0) * 2e-5;
        }
        let force = [6e-4, -3e-4, 2e-4];
        for tau in [0.56, 0.85, 1.3] {
            let bgk = collide_cell3(populations, CollisionModel3::Bgk { tau }, force).unwrap();
            let central = collide_cell3(
                populations,
                CollisionModel3::CentralMoment {
                    second_order_rate: 1.0 / tau,
                    higher_order_rate: 1.0 / tau,
                },
                force,
            )
            .unwrap();
            for q in 0..Q3 {
                assert!((central[q] - bgk[q]).abs() < 5e-13, "tau={tau}, q={q}");
            }
        }
    }

    #[test]
    fn forced_moment_modes_relax_physical_stress_at_the_viscous_rate() {
        // G0: physical Pi_neq is independent of the force-frame shift; both
        // the nonlinear fourth-order projection and independent higher rates
        // must leave the prescribed second-moment relaxation law intact.
        let (rho, u, force) = (1.02, [0.05, -0.03, 0.02], [7e-4, -2e-4, 5e-4]);
        let a = [
            [3e-4, -2e-4, 1e-4],
            [-2e-4, -1e-4, 2.5e-4],
            [1e-4, 2.5e-4, 4e-4],
        ];
        let eq = equilibrium3(rho, u);
        let source = source(u, force);
        let trace = a[0][0] + a[1][1] + a[2][2];
        let initial = core::array::from_fn(|q| {
            let c = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
            let mut cac = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    cac += c[i] * c[j] * a[i][j];
                }
            }
            eq[q] + 4.5 * W3[q] * (cac - trace / 3.0) - 0.5 * source[q]
        });
        for model in MODELS {
            let post = collide_cell3(initial, model, force).unwrap();
            let mut observed = [[0.0; 3]; 3];
            for q in 0..Q3 {
                let c = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
                let physical_neq = post[q] - eq[q] - 0.5 * source[q];
                for i in 0..3 {
                    for j in 0..3 {
                        observed[i][j] += c[i] * c[j] * physical_neq;
                    }
                }
            }
            for i in 0..3 {
                for j in 0..3 {
                    let expected = (1.0 - 1.23) * a[i][j];
                    assert!(
                        (observed[i][j] - expected).abs() < 5e-13,
                        "{model:?}, stress {i}{j}"
                    );
                }
            }
        }
    }
}

/// A D3Q19 duct: halfway bounce-back walls on the x and y boundaries,
/// periodic in z, driven by a body force `gz` along z — the 3-D
/// Poiseuille fixture (plan §14.1). Densities start at 1, velocities at
/// rest, unless seeded through [`Duct::perturb`].
pub struct Duct {
    nx: usize,
    ny: usize,
    nz: usize,
    tau: f64,
    gz: f64,
    /// SoA distributions: one tile-major field per population.
    f: [Vec<Tile>; Q3],
    /// Post-collision scratch (pull streaming reads this).
    post: [Vec<Tile>; Q3],
}

impl Duct {
    /// A duct at rest (unit density) with relaxation time `tau` and body
    /// force `gz`. Every dimension must be a positive multiple of
    /// [`TILE`].
    ///
    /// # Panics
    /// If any dimension is zero or not a multiple of [`TILE`], if `tau` is not
    /// finite and greater than one half, or if `gz` is not finite.
    #[must_use]
    pub fn new(nx: usize, ny: usize, nz: usize, tau: f64, gz: f64) -> Duct {
        assert!(
            nx > 0
                && ny > 0
                && nz > 0
                && nx.is_multiple_of(TILE)
                && ny.is_multiple_of(TILE)
                && nz.is_multiple_of(TILE),
            "duct dimensions must be positive multiples of {TILE} (got {nx}x{ny}x{nz})"
        );
        CollisionModel3::Bgk { tau }
            .validate()
            .expect("duct relaxation time must be finite and greater than 0.5");
        assert!(gz.is_finite(), "duct body force must be finite");
        let tiles = (nx / TILE) * (ny / TILE) * (nz / TILE);
        let f0 = equilibrium3(1.0, [0.0; 3]);
        let f = core::array::from_fn(|i| vec![Tile::filled(f0[i]); tiles]);
        let post = core::array::from_fn(|i| vec![Tile::filled(f0[i]); tiles]);
        Duct {
            nx,
            ny,
            nz,
            tau,
            gz,
            f,
            post,
        }
    }

    /// Tile index and local lane of cell `(x, y, z)` — the tile-major
    /// address map (x-fastest at both levels; the pinned order).
    #[inline]
    fn addr(&self, x: usize, y: usize, z: usize) -> (usize, usize) {
        let (ntx, nty) = (self.nx / TILE, self.ny / TILE);
        let tile = (z / TILE * nty + y / TILE) * ntx + x / TILE;
        let lane = (z % TILE * TILE + y % TILE) * TILE + x % TILE;
        (tile, lane)
    }

    /// The kinematic viscosity `ν = (τ − ½)/3`.
    #[must_use]
    pub fn viscosity(&self) -> f64 {
        (self.tau - 0.5) / 3.0
    }

    /// The macroscopic density at `(x, y, z)`.
    #[must_use]
    pub fn density(&self, x: usize, y: usize, z: usize) -> f64 {
        let (tile, lane) = self.addr(x, y, z);
        (0..Q3).map(|i| self.f[i][tile].0[lane]).sum()
    }

    /// The macroscopic velocity at `(x, y, z)` (with the Guo half-force
    /// momentum correction).
    #[must_use]
    pub fn velocity(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        let (tile, lane) = self.addr(x, y, z);
        let mut rho = 0.0;
        let mut m = [0.0; 3];
        for (i, e) in E3.iter().enumerate() {
            let fi = self.f[i][tile].0[lane];
            rho += fi;
            m[0] += f64::from(e.0) * fi;
            m[1] += f64::from(e.1) * fi;
            m[2] += f64::from(e.2) * fi;
        }
        // Guo: half the force is added to the momentum.
        [m[0] / rho, m[1] / rho, (m[2] + 0.5 * self.gz) / rho]
    }

    /// Total mass (conserved by construction; drift is roundoff only).
    #[must_use]
    pub fn total_mass(&self) -> f64 {
        self.f
            .iter()
            .flat_map(|field| field.iter())
            .map(|tile| tile.0.iter().sum::<f64>())
            .sum()
    }

    /// Deterministically perturb the resting state: cell `(x, y, z)`
    /// gets its density shifted by `amplitude · h(x, y, z)` where `h`
    /// is a fixed integer hash mapped to `[-1, 1)` — the seeded golden
    /// fixture's initial condition (no RNG dependency; the hash IS the
    /// seed schedule).
    pub fn perturb(&mut self, seed: u64, amplitude: f64) {
        for z in 0..self.nz {
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
                    // map to [-1, 1)
                    let unit = (h >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0;
                    let rho = 1.0 + amplitude * unit;
                    let feq = equilibrium3(rho, [0.0; 3]);
                    let (tile, lane) = self.addr(x, y, z);
                    for (field, feq_i) in self.f.iter_mut().zip(feq) {
                        field[tile].0[lane] = feq_i;
                    }
                }
            }
        }
    }

    /// One collide-force-stream step (BGK + Guo forcing, pull scheme,
    /// halfway bounce-back x/y walls, periodic z). Traversal follows the
    /// pinned tile-major order documented at module level.
    pub fn step(&mut self) {
        // Collide + Guo forcing (pointwise): write post-collision
        // populations into `post`, visiting tiles then lanes ascending.
        let tiles = self.f[0].len();
        for tile in 0..tiles {
            let input = self.f.each_ref().map(|field| &field[tile].0);
            let mut output = self.post.each_mut().map(|field| &mut field[tile].0);
            simd::collide_bgk_axial_z_tile(&input, &mut output, self.tau, self.gz)
                .expect("Duct constructor and prior collision state admit BGK/Guo collision");
        }
        // Pull streaming: SIMD rows are pure moves and retain the scalar
        // source map bit-for-bit (x/y halfway bounce, periodic z).
        simd::stream_duct(&self.post, &mut self.f, self.nx, self.ny, self.nz);
    }

    /// Run `steps` steps.
    pub fn run(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
        }
    }

    /// The `z`-velocity over the duct cross-section at `z = 0`, row-major
    /// in `(y, x)` — the profile the analytic duct series certifies.
    #[must_use]
    pub fn z_velocity_section(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.nx * self.ny);
        for y in 0..self.ny {
            for x in 0..self.nx {
                out.push(self.velocity(x, y, 0)[2]);
            }
        }
        out
    }
}

impl core::fmt::Debug for Duct {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("Duct")
            .field("nx", &self.nx)
            .field("ny", &self.ny)
            .field("nz", &self.nz)
            .field("tau", &self.tau)
            .field("gz", &self.gz)
            .finish_non_exhaustive()
    }
}

/// The analytic steady rectangular-duct `z`-velocity at lattice cell
/// `(x, y)` for a duct of `nx × ny` cells under body acceleration `gz`,
/// with halfway bounce-back walls at `x = −½`, `x = nx − ½`, `y = −½`,
/// `y = ny − ½` (so the full widths are exactly `nx` and `ny`):
///
/// `u(X,Y) = (16 gz a²)/(ν π³) Σ_{n odd} (−1)^((n−1)/2) n⁻³
///           [1 − cosh(nπY/2a)/cosh(nπb/2a)] cos(nπX/2a)`
///
/// with `a = nx/2`, `b = ny/2`, and `(X, Y)` the cell center relative to
/// the duct center. The series is truncated at `n = 99` (the 1/n³ decay
/// puts the tail below 1e-6 relative — far under the 3% acceptance bar);
/// aspect ratios up to ~4 stay clear of `cosh` overflow.
#[must_use]
pub fn duct_analytic(gz: f64, viscosity: f64, nx: usize, ny: usize, x: usize, y: usize) -> f64 {
    let a = nx as f64 / 2.0;
    let b = ny as f64 / 2.0;
    let cx = x as f64 - (nx as f64 - 1.0) / 2.0;
    let cy = y as f64 - (ny as f64 - 1.0) / 2.0;
    let mut sum = 0.0;
    let mut sign = 1.0;
    let mut n = 1u32;
    while n <= 99 {
        let nf = f64::from(n);
        // n³ and π³ as explicit products — `powi` is the build-mode
        // determinism hazard class (det:: doctrine / check-powi lint).
        let k = nf * core::f64::consts::PI / (2.0 * a);
        let term = (1.0 - (k * cy).cosh() / (k * b).cosh()) * (k * cx).cos() / (nf * nf * nf);
        sum += sign * term;
        sign = -sign;
        n += 2;
    }
    let pi = core::f64::consts::PI;
    16.0 * gz * a * a / (viscosity * (pi * pi * pi)) * sum
}
