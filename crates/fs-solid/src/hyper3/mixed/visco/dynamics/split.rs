//! Original quadratic Maxwell force/contact/drift with pressure equilibration.
//!
//! Forces come from accepted Piola stresses and original reference quadrature.
//! The native consistent mass applies the kick at fixed geometry. Contact may
//! change velocity on that same immutable state. Kinetic drift then fixes every
//! displacement coordinate and the original mixed solver determines only its
//! continuous pressure unknowns while retaining the original frozen histories.
//! Reactions from those prescribed drift coordinates are not physical supports.
//! Physical storage remains remainder plus U(J(p)), never U(det F).

use fs_math::dd::Dd;

use super::{
    Arc, Cx, FiniteViscoContactDrift, FiniteViscoContactDriftControls, MixedUpSettings,
    MixedViscoTetError, MixedViscoTetMotion, ReferenceMass, SolidImpulseControls,
    SolidImpulseReport, SolidMomentum, VolumetricTreatment, checkpoint, gate, inertia_error,
    invalid,
};
use crate::hyper3::visco::FiniteViscoStepControls;
use crate::inertia::impulse_moments;

/// Native force account bound privately to one accepted force phase.
/// Editing public diagnostics cannot change continuation or its physical work.
#[derive(Clone, Debug)]
pub struct MixedViscoForceKick<'a> {
    pub duration_s: f64,
    pub impulse: SolidImpulseReport,
    pub source_impulse: SolidMomentum,
    /// Minus the actual accepted mixed Piola force times the interval.
    pub internal_impulse: SolidMomentum,
    pub source_impulse_work_j: f64,
    pub internal_impulse_work_j: f64,
    pub work_decomposition_residual_j: f64,
    pub linear_decomposition_residual_kg_m_s: [f64; 3],
    pub angular_decomposition_residual_kg_m2_s: [f64; 3],
    /// Original physical remainder plus pressure-reconstructed bulk storage.
    pub stored_before_j: f64,
    binding: ForceKickBinding<'a>,
}

#[derive(Clone, Debug)]
struct ForceKickBinding<'a> {
    accepted_identity: Arc<()>,
    mass: Arc<ReferenceMass>,
    reference_nodes_m: &'a [[f64; 3]],
    duration_s: f64,
    origin_m: [f64; 3],
    held: Vec<usize>,
    external_force_n: Vec<f64>,
    controls: SolidImpulseControls,
    source_work_j: f64,
    internal_work_j: f64,
    kick_work_residual_j: f64,
    decomposition_residual_j: f64,
}

impl<'a> MixedViscoTetMotion<'a> {
    /// Apply h(F_external-F_internal) through the original supported P2 mass.
    ///
    /// All ten original signed shape weights project body loads; every retained
    /// mixed Piola stress enters the original reference gradient quadrature.
    /// Source and internal resultants accumulate with compensation. Held DOFs
    /// must be stationary. Geometry, pressure, time and history remain fixed,
    /// sharing immutable buffers across this kick and subsequent contact impulses.
    /// Exact incompressibility requires a different constraint-compatible split.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn force_kick_candidate(
        &self,
        duration_s: f64,
        held: &[usize],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        origin_m: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, MixedViscoForceKick<'a>), MixedViscoTetError> {
        checkpoint(cx)?;
        if self.solid.volumetric == VolumetricTreatment::ExactIncompressible {
            return Err(MixedViscoTetError::ExactIncompressibleSplit);
        }
        let n = self.mass.ndofs();
        let next_time = self.solid.time_s + duration_s;
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !next_time.is_finite()
            || next_time <= self.solid.time_s
            || body_force_n_m3
                .iter()
                .chain(&origin_m)
                .any(|v| !v.is_finite())
            || held.len() > n
            || nodal_forces_n.len() > controls.maximum_workspace_components
        {
            return Err(
                invalid("finite progressing force interval and bounded loads required").into(),
            );
        }
        self.stationary_split_supports(held, cx)?;
        // Dd accumulation, rounded source/internal vectors, sparse entries,
        // positions and retained support indices coexist with native CG/output.
        let extra = n
            .checked_mul(10)
            .and_then(|count| count.checked_add(held.len()))
            .ok_or_else(|| invalid("quadratic force split workspace overflow"))?;
        let solve_controls = SolidImpulseControls {
            maximum_workspace_components: controls
                .maximum_workspace_components
                .checked_sub(extra)
                .ok_or_else(|| invalid("quadratic force split exceeds workspace cap"))?,
            ..controls
        };
        let mut external_sum = vec![Dd::ZERO; n];
        let mut internal_sum = vec![Dd::ZERO; n];
        for (index, &(dof, value)) in nodal_forces_n.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            if dof >= n || !value.is_finite() {
                return Err(invalid("finite original nodal force components required").into());
            }
            external_sum[dof] = external_sum[dof] + Dd::from_f64(value);
        }
        let points = self.solid.mesh.points_per_element;
        for (element_index, element) in self.solid.mesh.elements.iter().enumerate() {
            for (point_index, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let piola = &self.solid.first_piola_pa[element_index * points + point_index];
                let weight = Dd::from_f64(point.weight_m3);
                for (a, &node) in element.connectivity.iter().enumerate() {
                    for axis in 0..3 {
                        let dof = 3 * node + axis;
                        external_sum[dof] = external_sum[dof]
                            + weight
                                * Dd::from_f64(point.shapes[a])
                                * Dd::from_f64(body_force_n_m3[axis]);
                        let force = (0..3).fold(Dd::ZERO, |sum, j| {
                            sum + Dd::from_f64(piola[3 * axis + j])
                                * Dd::from_f64(point.gradients[a][j])
                        });
                        internal_sum[dof] = internal_sum[dof] + weight * force;
                    }
                }
            }
        }
        let mut external = Vec::with_capacity(n);
        let mut source = Vec::with_capacity(n);
        let mut internal = Vec::with_capacity(n);
        let mut impulses = Vec::with_capacity(n);
        for (dof, (applied, resisting)) in external_sum.into_iter().zip(internal_sum).enumerate() {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            let force = finite(applied.to_f64())?;
            let source_impulse = finite(duration_s * force)?;
            let internal_impulse = finite(-duration_s * finite(resisting.to_f64())?)?;
            external.push(force);
            source.push(source_impulse);
            internal.push(internal_impulse);
            impulses.push((dof, finite(source_impulse + internal_impulse)?));
        }
        let mut event = self.impulse_candidate(held, &impulses, origin_m, solve_controls, cx)?;
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
            positions.push([
                finite(x[0] + u[0])?,
                finite(x[1] + u[1])?,
                finite(x[2] + u[2])?,
            ]);
        }
        let source_impulse =
            impulse_moments(&positions, &source, origin_m, cx).map_err(inertia_error)?;
        let internal_impulse =
            impulse_moments(&positions, &internal, origin_m, cx).map_err(inertia_error)?;
        if magnitude(internal_impulse.linear_kg_m_s) > controls.linear_tolerance_kg_m_s
            || magnitude(internal_impulse.angular_kg_m2_s) > controls.angular_tolerance_kg_m2_s
        {
            return Err(
                invalid("accepted mixed Piola force misses internal momentum objectivity").into(),
            );
        }
        let source_work = mean_work(&self.velocity_mps, &event.motion.velocity_mps, &source, cx)?;
        let internal_work = mean_work(
            &self.velocity_mps,
            &event.motion.velocity_mps,
            &internal,
            cx,
        )?;
        let decomposition = finite(event.impulse.applied_work_j - source_work - internal_work)?;
        gate(
            "quadratic force kick work decomposition",
            decomposition,
            controls.work_tolerance_j,
        )?;
        let linear = core::array::from_fn(|axis| {
            event.impulse.applied.linear_kg_m_s[axis]
                - source_impulse.linear_kg_m_s[axis]
                - internal_impulse.linear_kg_m_s[axis]
        });
        let angular = core::array::from_fn(|axis| {
            event.impulse.applied.angular_kg_m2_s[axis]
                - source_impulse.angular_kg_m2_s[axis]
                - internal_impulse.angular_kg_m2_s[axis]
        });
        if !linear.iter().chain(&angular).all(|v| v.is_finite())
            || magnitude(linear) > controls.linear_tolerance_kg_m_s
            || magnitude(angular) > controls.angular_tolerance_kg_m2_s
        {
            return Err(invalid(
                "quadratic force kick misses source/internal momentum decomposition",
            )
            .into());
        }
        let continuation = Arc::new(());
        event.motion.solid.accepted_identity = Arc::clone(&continuation);
        let binding = ForceKickBinding {
            accepted_identity: continuation,
            mass: Arc::clone(&self.mass),
            reference_nodes_m: self.solid.nodes_m,
            duration_s,
            origin_m,
            held: held.to_vec(),
            external_force_n: external,
            controls,
            source_work_j: source_work,
            internal_work_j: internal_work,
            kick_work_residual_j: event.impulse.work_residual_j,
            decomposition_residual_j: decomposition,
        };
        let receipt = MixedViscoForceKick {
            duration_s,
            impulse: event.impulse,
            source_impulse,
            internal_impulse,
            source_impulse_work_j: source_work,
            internal_impulse_work_j: internal_work,
            work_decomposition_residual_j: decomposition,
            linear_decomposition_residual_kg_m_s: linear,
            angular_decomposition_residual_kg_m2_s: angular,
            stored_before_j: self.solid.stored_energy_j(),
            binding,
        };
        checkpoint(cx)?;
        Ok((event.motion, receipt))
    }

    /// Drift with unchanged post-contact velocity and solve the original P1 pressure.
    ///
    /// Every endpoint displacement is prescribed to the represented kinetic
    /// drift, leaving ONLY continuous pressure as a native nonlinear unknown.
    /// The original native solve freezes each accepted Maxwell history once,
    /// advances it at the actual endpoint and retains U(J(p)) bulk storage.
    /// The prescribed-coordinate reactions serve only that kinematic material
    /// evaluation; none is a physical support impulse or work contribution.
    /// Finite compliance fixes pressure uniquely, so no pressure gauge is inferred.
    ///
    /// The shared drift receipt's native_loading_work_j includes native remainder
    /// loading plus the reversible bulk-storage change. Its material balance is
    /// the original remainder identity; the pressure work is independently gated
    /// by the original mixed owner. Split errors are numerical, never heat.
    /// Callers must still admit the complete geometric path and refine this split.
    #[allow(clippy::too_many_lines)]
    pub fn contact_drift_candidate(
        &self,
        kick: &MixedViscoForceKick<'a>,
        settings: MixedUpSettings,
        controls: FiniteViscoContactDriftControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, FiniteViscoContactDrift), MixedViscoTetError> {
        checkpoint(cx)?;
        if self.solid.volumetric == VolumetricTreatment::ExactIncompressible {
            return Err(MixedViscoTetError::ExactIncompressibleSplit);
        }
        let binding = &kick.binding;
        if !Arc::ptr_eq(&self.solid.accepted_identity, &binding.accepted_identity)
            || !Arc::ptr_eq(&self.mass, &binding.mass)
            || !std::ptr::eq(self.solid.nodes_m, binding.reference_nodes_m)
        {
            return Err(invalid("quadratic contact drift requires its original kick configuration, mass and history").into());
        }
        if [
            controls.material_balance_tolerance_j,
            controls.maximum_integration_defect_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(
                invalid("finite nonnegative material and split allowances required").into(),
            );
        }
        self.stationary_split_supports(&binding.held, cx)?;
        let time_s = finite(self.solid.time_s + binding.duration_s)?;
        if time_s <= self.solid.time_s {
            return Err(invalid("quadratic contact drift cannot advance its clock").into());
        }
        let mut prescribed = Vec::with_capacity(self.mass.ndofs());
        for (node, (u, v)) in self
            .solid
            .displacement_m
            .iter()
            .zip(&self.velocity_mps)
            .enumerate()
        {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            for axis in 0..3 {
                let value = if v[axis] == 0.0 {
                    u[axis]
                } else {
                    finite(binding.duration_s.mul_add(v[axis], u[axis]))?
                };
                prescribed.push((3 * node + axis, value));
            }
        }
        // This is the original pressure-only problem: all displacement entries
        // are constrained, no inertia is supplied, and finite bulk retains all
        // pressure unknowns. Native KKT/volume residual gates remain unchanged.
        let endpoint = self.solid.advance(
            binding.duration_s,
            &prescribed,
            &[],
            [0.0; 3],
            None,
            settings,
            FiniteViscoStepControls {
                equilibrium_work_tolerance_j: controls.material_balance_tolerance_j,
                maximum_endpoint_work_defect_j: controls.maximum_integration_defect_j,
            },
            cx,
        )?;
        let energy = &endpoint.energy;
        let momentum = self
            .mass
            .constant_velocity_drift(
                self.solid.nodes_m,
                &self.solid.displacement_m,
                &endpoint.solid.displacement_m,
                &self.velocity_mps,
                binding.origin_m,
                cx,
            )
            .map_err(inertia_error)?;
        if magnitude(momentum.change.linear_kg_m_s) > binding.controls.linear_tolerance_kg_m_s
            || magnitude(momentum.change.angular_kg_m2_s)
                > binding.controls.angular_tolerance_kg_m2_s
        {
            return Err(invalid(
                "quadratic constant-velocity drift misses original-mass momentum conservation",
            )
            .into());
        }
        let mut source_displacement_work = Dd::ZERO;
        for (dof, ((new, old), force)) in endpoint
            .solid
            .displacement_m
            .as_flattened()
            .iter()
            .zip(self.solid.displacement_m.as_flattened())
            .zip(&binding.external_force_n)
            .enumerate()
        {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            source_displacement_work = source_displacement_work
                + (Dd::from_f64(*new) - Dd::from_f64(*old)) * Dd::from_f64(*force);
        }
        let source_displacement_work_j = finite(source_displacement_work.to_f64())?;
        let source_split = finite(binding.source_work_j - source_displacement_work_j)?;
        let force_split = finite(binding.internal_work_j + energy.internal_endpoint_work_j)?;
        let parts = [
            source_split,
            force_split,
            -energy.endpoint_work_defect_j,
            binding.kick_work_residual_j,
            binding.decomposition_residual_j,
        ];
        let integration = finite(parts.iter().sum())?;
        let absolute = finite(parts.iter().map(|v| v.abs()).sum())?;
        gate(
            "absolute quadratic force/material splitting work",
            absolute,
            controls.maximum_integration_defect_j,
        )?;
        let receipt = FiniteViscoContactDrift {
            duration_s: binding.duration_s,
            time_s,
            before: momentum.before,
            after: momentum.after,
            linear_residual_kg_m_s: momentum.change.linear_kg_m_s,
            angular_residual_kg_m2_s: momentum.change.angular_kg_m2_s,
            kinetic_before_j: momentum.kinetic_j,
            kinetic_after_j: momentum.kinetic_j,
            kinetic_change_j: 0.0,
            stored_before_j: energy.stored_before_j,
            stored_after_j: energy.stored_after_j,
            dissipated_j: energy.dissipated_j,
            native_loading_work_j: finite(
                (Dd::from_f64(energy.remainder_loading_work_j) + Dd::from_f64(energy.bulk_after_j)
                    - Dd::from_f64(energy.bulk_before_j))
                .to_f64(),
            )?,
            material_balance_error_j: energy.material_balance_error_j,
            internal_endpoint_work_j: energy.internal_endpoint_work_j,
            material_endpoint_defect_j: energy.endpoint_work_defect_j,
            source_displacement_work_j,
            source_splitting_work_j: source_split,
            force_splitting_work_j: force_split,
            integration_defect_j: integration,
            absolute_integration_defect_j: absolute,
        };
        let candidate = Self {
            solid: endpoint.solid,
            velocity_mps: self.velocity_mps.clone(),
            reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
            mass: Arc::clone(&self.mass),
        };
        checkpoint(cx)?;
        Ok((candidate, receipt))
    }

    fn stationary_split_supports(
        &self,
        held: &[usize],
        cx: &Cx<'_>,
    ) -> Result<(), MixedViscoTetError> {
        for (index, &dof) in held.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            if self
                .velocity_mps
                .get(dof / 3)
                .is_none_or(|v| v[dof % 3] != 0.0)
            {
                return Err(invalid(
                    "quadratic force/contact/drift requires stationary held velocities",
                )
                .into());
            }
        }
        Ok(())
    }
}

fn mean_work(
    before: &[[f64; 3]],
    after: &[[f64; 3]],
    impulse: &[f64],
    cx: &Cx<'_>,
) -> Result<f64, MixedViscoTetError> {
    let mut work = Dd::ZERO;
    for (dof, ((old, new), impulse)) in before
        .as_flattened()
        .iter()
        .zip(after.as_flattened())
        .zip(impulse)
        .enumerate()
    {
        if dof % 256 == 0 {
            checkpoint(cx)?;
        }
        work =
            work + (Dd::from_f64(0.5 * *old) + Dd::from_f64(0.5 * *new)) * Dd::from_f64(*impulse);
    }
    finite(work.to_f64())
}

fn finite(value: f64) -> Result<f64, MixedViscoTetError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("quadratic force/contact/drift arithmetic is unrepresentable").into())
    }
}

fn magnitude(value: [f64; 3]) -> f64 {
    value[0].hypot(value[1]).hypot(value[2])
}
