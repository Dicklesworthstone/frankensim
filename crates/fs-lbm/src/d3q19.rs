//! The D3Q19 core (bead 84hv): 3-D BGK stream-and-collide with Guo body
//! forcing and halfway bounce-back walls, over 128-byte-aligned SoA
//! distributions in tile-major layout (plan §5) with the PULL scheme —
//! the 3-D sibling of the D2Q9 module in `lib.rs`, which it deliberately
//! mirrors (and leaves untouched).
//!
//! DETERMINISTIC TRAVERSAL (pinned): both the collide and stream loops
//! visit tiles in ascending tile index (x-fastest, then y, then z over
//! the tile grid) and cells within a tile in ascending local index
//! (x-fastest, then y, then z). Single-threaded in this bead; when WS1-D
//! parallelizes the sweep, per-cell writes stay slot-exclusive (collide
//! is pointwise; pull-streaming writes only the destination cell), so
//! this order is documentation of the CANONICAL schedule, not a hidden
//! result dependence — results must stay bit-identical across worker
//! counts by construction.
//!
//! Layout: the domain is split into 4×4×4 tiles (64 cells, 512 B per
//! field per tile, 128-byte aligned). Each of the 19 distribution
//! fields is its own `Vec<Tile>` (structure of arrays); cell `(x,y,z)`
//! lives in tile `(x/4, y/4, z/4)` at local `(x%4, y%4, z%4)`.
//! Dimensions must be multiples of 4 (asserted): the fixture scales in
//! WS1 are, and padding is a later concern, not silent slop.

mod boundary;

pub use boundary::{
    BoundaryGrid3, BoundaryLink3, BoundarySpec3, D3Q19_BOUNDARY_BIT_SEMANTICS_VERSION, Face3,
    FaceBoundary3, LinkMaskTile3,
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
/// golden evidence.
pub const D3Q19_MOMENT_COLLISION_SEMANTICS_VERSION: u32 = 1;

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
/// cumulant rungs are unforced, correctness-first reference operators; neither
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
    /// Body forcing is deliberately refused until a moment-space forcing
    /// contract is separately verified.
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
    /// Central-moment forcing is not yet admitted by its contract.
    CentralMomentForceUnsupported {
        /// Rejected body-force vector.
        force: [f64; 3],
    },
    /// Reduced-cumulant forcing is not yet admitted by its contract.
    ReducedCumulantForceUnsupported {
        /// Rejected body-force vector.
        force: [f64; 3],
    },
    /// The centered monomial transform could not be solved reliably.
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
/// this public entry point uses the frozen boundary-grid expression, while
/// [`Duct`] delegates through its private frozen axial-force expression. This
/// extraction therefore centralizes authority without silently refreezing
/// either golden surface.
///
/// # Errors
/// [`CollisionError3`] for inadmissible model parameters, force, incoming
/// state, macroscopic state, or outgoing populations.
pub fn collide_cell3(
    populations: [f64; Q3],
    model: CollisionModel3,
    force: [f64; 3],
) -> Result<[f64; Q3], CollisionError3> {
    collide_cell3_with_projection(populations, model, force, false)
}

fn collide_axial_z_cell3(
    populations: [f64; Q3],
    model: CollisionModel3,
    gz: f64,
) -> Result<[f64; Q3], CollisionError3> {
    collide_cell3_with_projection(populations, model, [0.0, 0.0, gz], true)
}

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
            let cs2 = 1.0 / 3.0;
            let cs4 = cs2 * cs2;
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
                    let projection = (0..3)
                        .map(|axis| {
                            ((e[axis] - velocity[axis]) / cs2 + eu * e[axis] / cs4) * force[axis]
                        })
                        .sum::<f64>();
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
        CollisionModel3::CentralMoment {
            second_order_rate,
            higher_order_rate,
        } => {
            if has_nonzero_force3(force) {
                return Err(CollisionError3::CentralMomentForceUnsupported { force });
            }
            collide_central_moments3(
                populations,
                equilibrium,
                velocity,
                second_order_rate,
                higher_order_rate,
            )
        }
        CollisionModel3::ReducedCumulant {
            second_order_rate,
            third_order_rate,
            fourth_order_rate,
        } => {
            if has_nonzero_force3(force) {
                return Err(CollisionError3::ReducedCumulantForceUnsupported { force });
            }
            collide_reduced_cumulants3(
                populations,
                equilibrium,
                velocity,
                second_order_rate,
                third_order_rate,
                fourth_order_rate,
            )
        }
    }
}

fn has_nonzero_force3(force: [f64; 3]) -> bool {
    force.iter().any(|component| component.abs() > 0.0)
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

fn apply_moment_matrix3(matrix: &[[f64; Q3]; Q3], populations: [f64; Q3]) -> [f64; Q3] {
    core::array::from_fn(|moment| {
        let mut value = 0.0;
        for (direction, population) in populations.into_iter().enumerate() {
            value = matrix[moment][direction].mul_add(population, value);
        }
        value
    })
}

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
    let post = solve_moment_system3(matrix, relaxed)?;
    for (direction, value) in post.into_iter().enumerate() {
        if !value.is_finite() {
            return Err(CollisionError3::NonFiniteOutput { direction, value });
        }
    }
    Ok(post)
}

fn fourth_order_pair_product3(first: f64, second: f64, cross: f64) -> f64 {
    first.mul_add(second, 2.0 * cross * cross)
}

// Geier et al. (2015), Eqs. 46-52, use unnormalized cumulants C = rho*c.
// On this reduced D3Q19 basis the only independent order-four instances are
// C220, C202, and C022; each subtracts its paired normal stresses and twice
// the squared cross stress divided by density.
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
    let post = solve_moment_system3(matrix, relaxed_moments)?;
    for (direction, value) in post.into_iter().enumerate() {
        if !value.is_finite() {
            return Err(CollisionError3::NonFiniteOutput { direction, value });
        }
    }
    Ok(post)
}

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
            for lane in 0..TILE_CELLS {
                let populations = core::array::from_fn(|direction| self.f[direction][tile].0[lane]);
                let post = collide_axial_z_cell3(
                    populations,
                    CollisionModel3::Bgk { tau: self.tau },
                    self.gz,
                )
                .expect("Duct constructor and prior collision state admit BGK/Guo collision");
                for (field, value) in self.post.iter_mut().zip(post) {
                    field[tile].0[lane] = value;
                }
            }
        }
        // Pull streaming: destination (x,y,z) reads post[source] where
        // source = destination − e_i; crossing an x or y boundary means
        // the population came off a wall — halfway bounce-back reflects
        // the OPPOSITE post-collision population of the destination cell
        // itself; z wraps (periodic). Offsets are ±1 only.
        for z in 0..self.nz {
            for y in 0..self.ny {
                for x in 0..self.nx {
                    let (dtile, dlane) = self.addr(x, y, z);
                    for i in 0..Q3 {
                        let (ex, ey, ez) = E3[i];
                        let sx = match ex {
                            1 if x == 0 => None,
                            1 => Some(x - 1),
                            -1 if x + 1 == self.nx => None,
                            -1 => Some(x + 1),
                            _ => Some(x),
                        };
                        let sy = match ey {
                            1 if y == 0 => None,
                            1 => Some(y - 1),
                            -1 if y + 1 == self.ny => None,
                            -1 => Some(y + 1),
                            _ => Some(y),
                        };
                        let sz = match ez {
                            1 => (z + self.nz - 1) % self.nz,
                            -1 => (z + 1) % self.nz,
                            _ => z,
                        };
                        self.f[i][dtile].0[dlane] = match (sx, sy) {
                            (Some(sx), Some(sy)) => {
                                let (stile, slane) = self.addr(sx, sy, sz);
                                self.post[i][stile].0[slane]
                            }
                            // wall crossing: bounce back in place.
                            _ => self.post[OPP3[i]][dtile].0[dlane],
                        };
                    }
                }
            }
        }
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
