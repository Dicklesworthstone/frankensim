//! Backward-Euler inertia on the original P1 inelastic Newton solve.
//!
//! The same endpoint constitutive stress, tangent, committed history and
//! backtracking owner solve `M (v_next-v_old)/dt + f_internal = f_external`,
//! with `v_next = (u_next-u_old)/dt`. Reference density supplies consistent
//! mass through the elastic owner's shared P1 element kernel. No material
//! response, geometry, constraints or return mapping is replaced here.
//!
//! Endpoint forces and support reactions do work through the actual increment.
//! Their inertial part equals the kinetic change plus `1/2 delta_v M delta_v`.
//! The latter is time-discretization loss, never physical material heat.
//! Constitutive storage/dissipation and thermal feedback remain with their
//! existing owners. Refine the physical timestep to control this first-order
//! approximation; a work identity is not a temporal-accuracy certificate.

pub mod fatigue;

use super::*;

pub(super) use crate::inertia::BackwardEulerInertia as Inertia;
/// Original public kinetic-work receipt, shared by both deformation regimes.
pub use crate::inertia::InertialEnergy;
use crate::inertia::{InertiaError, P1ReferenceMass};

/// Complete proposed motion; the original state/history is never mutated.
#[derive(Clone, Debug)]
pub struct BackwardEulerStep<S> {
    /// Native stress, displacement, plastic history and full support reactions.
    pub solution: SmallStrainTetSolution<S>,
    pub velocity_mps: Vec<[f64; 3]>,
    pub energy: InertialEnergy,
}

fn invalid(what: &str) -> SmallStrainTetError {
    SmallStrainTetError::InvalidInput { what: what.into() }
}
impl From<InertiaError> for SmallStrainTetError {
    fn from(error: InertiaError) -> Self {
        match error {
            InertiaError::InvalidInput { what } => Self::InvalidInput { what },
            InertiaError::Cancelled => Self::Cancelled,
        }
    }
}

impl<L: SmallStrainLaw> SmallStrainTetProblem<'_, L> {
    /// One physical inertial interval from the accepted full displacement,
    /// velocity and element history. `load_path` must contain only `1.0`:
    /// forces and prescribed displacements are actual endpoint values.
    /// Support reactions include consistent-mass coupling to constrained DOFs.
    /// Free bodies need no artificial restraints.
    ///
    /// Densities belong to the original reference inventory, in tet order.
    /// A material/thermal caller must retain that inventory across intervals;
    /// changing it requires a separate momentum/matter-transfer law. Every
    /// Newton proposal uses the same old constitutive history. This method
    /// reports inertial work but does not relabel it as constitutive heat.
    /// A motion increment too small to resolve against accumulated absolute
    /// displacement can refuse convergence; force tolerances are not relaxed.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_backward_euler(
        &self,
        displacement_m: &[f64],
        velocity_mps: &[[f64; 3]],
        element_states: &[L::State],
        reference_density_kg_m3: &[f64],
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<BackwardEulerStep<L::State>, SmallStrainTetError> {
        checkpoint(cx)?;
        if self.load_path != [1.0] {
            return Err(invalid(
                "one endpoint load at factor one is required per inertial interval",
            ));
        }
        let elements = self.validate()?;
        self.check(displacement_m, element_states)?;
        let mass = P1ReferenceMass::new(
            elements
                .iter()
                .map(|element| (element.nodes, element.volume)),
            reference_density_kg_m3,
            displacement_m.len(),
            cx,
        )?;
        let inertia = Inertia::new(&mass, displacement_m, velocity_mps.as_flattened(), dt_s, cx)?;
        let solution =
            self.solve_candidate(Some((displacement_m, element_states)), Some(&inertia), cx)?;
        let (velocity_mps, energy) = inertia.motion(solution.displacement_m.as_flattened(), cx)?;
        checkpoint(cx)?;
        Ok(BackwardEulerStep {
            solution,
            velocity_mps,
            energy,
        })
    }
}
