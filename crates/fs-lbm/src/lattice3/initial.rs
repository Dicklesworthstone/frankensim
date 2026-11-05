//! Physical initial density/velocity in the native half-source frame.

use super::{D3q19, Grid3, GridBudget3, GridError3, sources::ForcedLattice3};
use crate::d3q27::{CollisionModel27, D3q27, Equilibrium27};

mod sealed {
    use super::ForcedLattice3;

    pub trait NativeEquilibrium<const Q: usize>: ForcedLattice3<Q> {
        fn equilibrium(collision: Self::Collision, rho: f64, velocity: [f64; 3]) -> [f64; Q];
        fn maximum_mach(collision: Self::Collision) -> f64;
        fn attach_force(
            populations: &mut [f64; Q],
            collision: Self::Collision,
            velocity: [f64; 3],
            force: [f64; 3],
        );
    }
}

/// Native lattices with a known equilibrium and half-source representation.
/// Sealed so arbitrary quadratures cannot inherit native moment identities.
pub trait EquilibriumLattice3<const Q: usize>: sealed::NativeEquilibrium<Q> {}
impl EquilibriumLattice3<19> for D3q19 {}
impl EquilibriumLattice3<27> for D3q27 {}

impl sealed::NativeEquilibrium<19> for D3q19 {
    fn equilibrium(_: Self::Collision, rho: f64, velocity: [f64; 3]) -> [f64; 19] {
        crate::d3q19::equilibrium3(rho, velocity)
    }

    fn maximum_mach(_: Self::Collision) -> f64 {
        crate::MACH_LIMIT
    }

    fn attach_force(
        populations: &mut [f64; 19],
        _: Self::Collision,
        velocity: [f64; 3],
        force: [f64; 3],
    ) {
        crate::thermal::forcing::reframe3(populations, velocity, [0.0; 3], force);
    }
}

fn equilibrium_form(model: CollisionModel27) -> Equilibrium27 {
    match model {
        CollisionModel27::Bgk { equilibrium, .. } | CollisionModel27::Trt { equilibrium, .. } => {
            equilibrium
        }
        CollisionModel27::CentralMoment { .. } | CollisionModel27::Cumulant { .. } => {
            Equilibrium27::Factorized
        }
    }
}

impl sealed::NativeEquilibrium<27> for D3q27 {
    fn equilibrium(collision: Self::Collision, rho: f64, velocity: [f64; 3]) -> [f64; 27] {
        crate::d3q27::equilibrium27(rho, velocity, equilibrium_form(collision.model))
    }

    fn maximum_mach(collision: Self::Collision) -> f64 {
        collision.maximum_mach
    }

    fn attach_force(
        populations: &mut [f64; 27],
        collision: Self::Collision,
        velocity: [f64; 3],
        force: [f64; 3],
    ) {
        use crate::d3q27::{CS2_27, E27, W27, factorized_source27};
        let half_force = force.map(|f| -0.5 * f);
        let source = match equilibrium_form(collision.model) {
            Equilibrium27::Factorized => factorized_source27(velocity, half_force),
            Equilibrium27::Quadratic => core::array::from_fn(|q| {
                let e = E27[q].map(f64::from);
                let eu = (0..3).map(|a| e[a] * velocity[a]).sum::<f64>();
                W27[q]
                    * (0..3)
                        .map(|a| {
                            ((e[a] - velocity[a]) / CS2_27 + eu * e[a] / (CS2_27 * CS2_27))
                                * half_force[a]
                        })
                        .sum::<f64>()
            }),
        };
        let rho = populations.iter().sum::<f64>();
        for (population, change) in populations.iter_mut().zip(source) {
            *population += change;
        }
        // The source has zero mass. Keep its summation roundoff in the rest
        // population, which changes neither momentum nor the stress tensor.
        populations[0] += rho - populations.iter().sum::<f64>();
    }
}

pub(super) fn admit_velocity<L: EquilibriumLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    velocity: [f64; 3],
    cell: usize,
    maximum: f64,
) -> Result<(), GridError3<L::Error>> {
    let maximum = maximum.min(L::maximum_mach(collision));
    let mach = velocity[0].hypot(velocity[1]).hypot(velocity[2]) / L::CS2.sqrt();
    if !mach.is_finite() || mach > maximum {
        return Err(GridError3::Mach {
            cell,
            observed: mach,
            maximum,
        });
    }
    Ok(())
}

pub(super) fn attach_force<L: EquilibriumLattice3<Q>, const Q: usize>(
    populations: &mut [f64; Q],
    collision: L::Collision,
    velocity: [f64; 3],
    force: [f64; 3],
) {
    // Preserve all original population bits, including signed zeros, when
    // there is no coordinate change to make.
    if force != [0.0; 3] {
        L::attach_force(populations, collision, velocity, force);
    }
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> Grid3<L, Q> {
    /// Initialize physical `(density, velocity)` in lattice units, with the
    /// requested uniform force density already represented by `-S(u,F)/2`.
    /// The initial physical stress is `rho * (cs² I + u u)`, without the false
    /// stress introduced by evaluating equilibrium at `u - F/(2 rho)`.
    ///
    /// D3Q19 uses its quadratic equilibrium and full Guo source. D3Q27 BGK/TRT
    /// retain the chosen equilibrium and its complete derivative source; its
    /// moment models start from the product equilibrium. This specifies
    /// hydrodynamic initial moments, not an additional finite-force accuracy
    /// claim for the moment collision models. No time, impulse or heat advances.
    ///
    /// The same topology, allocation and cell budgets as [`Self::try_new`]
    /// apply. Requested and represented physical velocities must satisfy the
    /// grid envelope and the collision's own Mach limit. No population clipping
    /// or inferred equilibrium choice is applied to BGK/TRT.
    pub fn try_new_equilibrium(
        dims: [usize; 3],
        collision: L::Collision,
        force: [f64; 3],
        periodic: [bool; 3],
        budget: GridBudget3,
        solid: impl Fn([usize; 3]) -> bool,
        initial: impl Fn([usize; 3]) -> (f64, [f64; 3]),
    ) -> Result<Self, GridError3<L::Error>> {
        Self::try_new_with_initial_cells(
            dims,
            collision,
            force,
            periodic,
            budget,
            solid,
            |cell, p| {
                let (rho, velocity) = initial(p);
                if !rho.is_finite() || rho <= 0.0 || velocity.iter().any(|u| !u.is_finite()) {
                    return Err(GridError3::NonFiniteCell { cell });
                }
                admit_velocity::<L, Q>(collision, velocity, cell, budget.maximum_mach)?;
                if force != [0.0; 3] && !L::supports_body_force(collision) {
                    return Err(GridError3::MaterialField {
                        cell,
                        reason: "collision model does not admit body force",
                    });
                }
                let mut populations = L::equilibrium(collision, rho, velocity);
                attach_force::<L, Q>(&mut populations, collision, velocity, force);
                let (_, actual) =
                    super::macro_cell::<L, Q>(&populations, force, cell, budget.maximum_mach)?;
                admit_velocity::<L, Q>(collision, actual, cell, budget.maximum_mach)?;
                Ok(populations)
            },
        )
    }
}
