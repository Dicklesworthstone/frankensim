//! Zero-duration momentum transfer with the original mixed viscoelastic state.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::hyper::mixed::VolumetricTreatment;

use super::{
    MixedImpulseControls, MixedImpulseReport, MixedViscoTetError, MixedViscoTetMotion,
    MixedViscoTetSolid, SolidImpulseControls, SolidImpulseReport, SolidPointKinematics,
    SolidPointMobility, SolidPointMobilityMatrix, Tet10ImpulsePoint, checkpoint, inertia_error,
};
use crate::inertia::point::tet10;

mod constrained;
pub use constrained::MixedViscoConstrainedResponse;

/// Complete candidate retaining continuous pressure and every Maxwell history.
#[derive(Debug)]
pub struct MixedViscoMotionImpulse<'a> {
    pub motion: MixedViscoTetMotion<'a>,
    /// Actual momentum, kinetic change, applied impulses and signed support work.
    pub impulse: SolidImpulseReport,
}

/// Complete exact-volume candidate with separate pressure impulse accounting.
#[derive(Debug)]
pub struct MixedViscoConstrainedImpulse<'a> {
    pub motion: MixedViscoTetMotion<'a>,
    /// Pressure impulses [Pa s] never overwrite retained physical pressure [Pa].
    pub impulse: MixedImpulseReport,
}

impl<'a> MixedViscoTetMotion<'a> {
    /// Current material-point position and work-conjugate quadratic velocity.
    /// Uses all ten signed shape weights and the actual deformed geometry.
    /// Sampling is available for both original finite-compliance and exact
    /// pressure formulations; it changes no material state.
    pub fn point_kinematics(
        &self,
        point: Tet10ImpulsePoint,
        cx: &Cx<'_>,
    ) -> Result<SolidPointKinematics, MixedViscoTetError> {
        tet10::kinematics(
            point,
            self.solid.elements,
            self.solid.nodes_m,
            &self.solid.displacement_m,
            &self.velocity_mps,
            cx,
        )
        .map_err(|error| MixedViscoTetError::from(inertia_error(error)))
    }

    /// Supported directional inverse mass at an original quadratic material point.
    /// The direction must be unit length in world axes. The original consistent
    /// P2 mass supplies the response without a history clone or material trial.
    /// Exact incompressibility refuses: its pressure constraint must participate
    /// in mobility, and the mass-only value would give a different impact law.
    pub fn point_directional_inverse_mass_per_kg(
        &self,
        held: &[usize],
        point: Tet10ImpulsePoint,
        direction: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<SolidPointMobility, MixedViscoTetError> {
        checkpoint(cx)?;
        if self.solid.volumetric == VolumetricTreatment::ExactIncompressible {
            return Err(MixedViscoTetError::ExactIncompressibleImpulse);
        }
        tet10::directional_mobility(
            &self.mass,
            self.solid.elements,
            held,
            point,
            direction,
            controls,
            cx,
        )
        .map_err(|error| MixedViscoTetError::from(inertia_error(error)))
    }

    /// Complete simultaneous cross-point/direction mobility [1/kg].
    /// Each row uses a unit world direction and all ten signed quadratic weights.
    /// Original consistent mass supplies every cross term with bounded workspace,
    /// keeping one full response at a time and copying no material history.
    /// Exact incompressibility refuses this mass-only entry.
    pub fn point_mobility_matrix_per_kg(
        &self,
        held: &[usize],
        rows: &[(Tet10ImpulsePoint, [f64; 3])],
        maximum_rows: usize,
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<SolidPointMobilityMatrix, MixedViscoTetError> {
        checkpoint(cx)?;
        if self.solid.volumetric == VolumetricTreatment::ExactIncompressible {
            return Err(MixedViscoTetError::ExactIncompressibleImpulse);
        }
        tet10::mobility_matrix(
            &self.mass,
            self.solid.elements,
            held,
            rows,
            maximum_rows,
            controls,
            cx,
        )
        .map_err(|error| MixedViscoTetError::from(inertia_error(error)))
    }

    /// Apply all quadratic material-point impulses in one native velocity update.
    /// The work-conjugate traces retain curved geometry, signed corner weights,
    /// and small net nodal loads through original compensated accumulation.
    /// One complete candidate retains physical pressure and every Maxwell state.
    /// Exact incompressibility refuses this mass-only entry before distribution.
    pub fn point_impulses_candidate(
        &self,
        held: &[usize],
        impulses: &[(Tet10ImpulsePoint, [f64; 3])],
        origin: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<MixedViscoMotionImpulse<'a>, MixedViscoTetError> {
        checkpoint(cx)?;
        if self.solid.volumetric == VolumetricTreatment::ExactIncompressible {
            return Err(MixedViscoTetError::ExactIncompressibleImpulse);
        }
        let entries = tet10::distributions(
            self.solid.elements,
            self.velocity_mps.len(),
            impulses,
            controls,
            cx,
        )
        .map_err(inertia_error)?;
        self.impulse_candidate(held, &entries, origin, controls, cx)
    }

    /// Apply a world impulse through the transpose of the same quadratic trace.
    /// Original signed weights preserve point torque and work at the actual
    /// curved geometry; the native nodal owner retains pressure and all memory.
    /// Exact incompressibility refuses this mass-only entry before distribution.
    pub fn point_impulse_candidate(
        &self,
        held: &[usize],
        point: Tet10ImpulsePoint,
        impulse_world: [f64; 3],
        origin: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<MixedViscoMotionImpulse<'a>, MixedViscoTetError> {
        checkpoint(cx)?;
        if self.solid.volumetric == VolumetricTreatment::ExactIncompressible {
            return Err(MixedViscoTetError::ExactIncompressibleImpulse);
        }
        let entries = tet10::distribution(
            point,
            self.solid.elements,
            self.velocity_mps.len(),
            impulse_world,
            cx,
        )
        .map_err(inertia_error)?;
        self.impulse_candidate(held, &entries, origin, controls, cx)
    }

    /// Apply instantaneous nodal impulses through the original consistent P2 mass.
    ///
    /// Only the finite-compliance `OriginalMaterial` formulation is admitted.
    /// Exact incompressibility returns `ExactIncompressibleImpulse`; use
    /// [`Self::constrained_impulse_candidate`] for its pressure-constrained solve.
    ///
    /// Held velocity DOFs keep their accepted values, including moving supports.
    /// The native mass owner retains signed corner/edge entries, repeated-load
    /// accumulation, bounded CG and actual free-row/momentum/work gates. Angular
    /// momentum uses current deformed positions about the supplied origin.
    ///
    /// Only velocity changes. Pressure, deformation, all branch metrics and
    /// dissipated histories, stresses, storage, time and reference density are
    /// unchanged. Subsequent `advance` calls continue those same histories.
    /// The caller owns the contact/restitution law and opposite body's impulse;
    /// this event creates no material relaxation loss or heat.
    ///
    /// # Errors
    /// Exact incompressibility, invalid inputs, native mass-solve/work/momentum
    /// refusal, or cancellation. No accepted field changes on any refusal.
    pub fn impulse_candidate(
        &self,
        held_velocity_dofs: &[usize],
        nodal_impulses_kg_m_s: &[(usize, f64)],
        origin_m: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<MixedViscoMotionImpulse<'a>, MixedViscoTetError> {
        checkpoint(cx)?;
        if self.solid.volumetric == VolumetricTreatment::ExactIncompressible {
            return Err(MixedViscoTetError::ExactIncompressibleImpulse);
        }
        let positions_m = self.impulse_positions(cx)?;
        let (velocity_mps, impulse) = self
            .mass
            .impulse_candidate(
                &positions_m,
                &self.velocity_mps,
                held_velocity_dofs,
                nodal_impulses_kg_m_s,
                origin_m,
                controls,
                cx,
            )
            .map_err(inertia_error)?;
        Ok(MixedViscoMotionImpulse {
            motion: self.with_impulse_velocity(velocity_mps, cx)?,
            impulse,
        })
    }

    /// Apply an instantaneous load with the original exact P2/P1 volume constraint.
    ///
    /// The shared mixed impulse owner solves the original mass/pressure KKT
    /// system using cof(F) at each accepted material point. It checks actual
    /// momentum, pressure zero-work and weak volume rates before publication.
    /// Stored pressure, every Maxwell history, deformation, energy and time
    /// remain unchanged. Held components retain their original velocity.
    ///
    /// Pre-event velocity must already satisfy the supplied weak-rate allowance.
    /// A backward-Euler displacement chord is not automatically compatible:
    /// this event cannot conceal a large rate defect by projecting it away.
    /// A zero pressure-impulse gauge is allowed only for a true connected
    /// constant null mode; it never replaces accepted physical pressure.
    ///
    /// # Errors
    /// Finite-compliance formulation, incompatible velocity/gauge, invalid input,
    /// native factor/balance/budget refusal or cancellation; no state is mutated.
    #[allow(clippy::too_many_arguments)]
    pub fn constrained_impulse_candidate(
        &self,
        held_velocity_dofs: &[usize],
        nodal_impulses_kg_m_s: &[(usize, f64)],
        pressure_gauge_node: Option<usize>,
        origin_m: [f64; 3],
        controls: MixedImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<MixedViscoConstrainedImpulse<'a>, MixedViscoTetError> {
        self.constrained_impulse_response(held_velocity_dofs, pressure_gauge_node, controls, cx)?
            .nodal_impulse_candidate(nodal_impulses_kg_m_s, origin_m, cx)
    }

    fn impulse_positions(&self, cx: &Cx<'_>) -> Result<Vec<[f64; 3]>, MixedViscoTetError> {
        let mut positions = Vec::with_capacity(self.velocity_mps.len());
        for (node, (x, u)) in self
            .solid
            .nodes_m
            .iter()
            .zip(self.solid.displacement_m.iter())
            .enumerate()
        {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            positions.push(core::array::from_fn(|c| x[c] + u[c]));
        }
        Ok(positions)
    }

    fn with_impulse_velocity(
        &self,
        velocity_mps: Vec<[f64; 3]>,
        cx: &Cx<'_>,
    ) -> Result<Self, MixedViscoTetError> {
        let solid = MixedViscoTetSolid {
            nodes_m: self.solid.nodes_m,
            elements: self.solid.elements,
            material: self.solid.material,
            volumetric: self.solid.volumetric,
            mesh: Arc::clone(&self.solid.mesh),
            budget: self.solid.budget,
            quadrature_points_per_dir: self.solid.quadrature_points_per_dir,
            maximum_quadrature_points: self.solid.maximum_quadrature_points,
            displacement_m: Arc::clone(&self.solid.displacement_m),
            pressure_pa: Arc::clone(&self.solid.pressure_pa),
            states: Arc::clone(&self.solid.states),
            first_piola_pa: Arc::clone(&self.solid.first_piola_pa),
            remainder_energy_j: self.solid.remainder_energy_j,
            bulk_energy_j: self.solid.bulk_energy_j,
            time_s: self.solid.time_s,
            accepted_identity: Arc::clone(&self.solid.accepted_identity),
        };
        checkpoint(cx)?;
        Ok(Self {
            solid,
            velocity_mps,
            reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
            mass: Arc::clone(&self.mass),
        })
    }
}
