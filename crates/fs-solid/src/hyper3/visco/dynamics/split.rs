//! First-order native force/contact/drift with original Maxwell memory.
//!
//! The accepted Piola field supplies the force kick; contact acts at the same
//! configuration. Drift retains the post-contact velocity and advances the
//! original material return map at its actual endpoint. No equilibrium solve,
//! artificial support reaction, endpoint projection or heat partition is made.
//! A caller must admit the complete geometric path and refine the explicit split.

use fs_math::dd::Dd;

use super::{
    Arc, Cx, FiniteViscoTetError, FiniteViscoTetMotion, FiniteViscoTetSolid, P1ReferenceMass,
    SolidImpulseControls, SolidImpulseReport, SolidMomentum, checkpoint, energy_gate,
    inertia_error, invalid,
};
use crate::inertia::impulse_moments;

/// Native force kick and a private binding to the original accepted state.
/// Public diagnostics may be edited; continuation uses only the private binding.
#[derive(Clone, Debug)]
pub struct FiniteViscoForceKick<'a> {
    pub duration_s: f64,
    pub impulse: SolidImpulseReport,
    pub source_impulse: SolidMomentum,
    /// Minus the actual accepted Piola force times the interval.
    pub internal_impulse: SolidMomentum,
    pub source_impulse_work_j: f64,
    pub internal_impulse_work_j: f64,
    pub work_decomposition_residual_j: f64,
    pub linear_decomposition_residual_kg_m_s: [f64; 3],
    pub angular_decomposition_residual_kg_m2_s: [f64; 3],
    pub stored_before_j: f64,
    binding: ForceKickBinding<'a>,
}

#[derive(Clone, Debug)]
struct ForceKickBinding<'a> {
    accepted_identity: Arc<()>,
    mass: Arc<P1ReferenceMass>,
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

/// Independent physical material identity and numerical split allowances [J].
#[derive(Clone, Copy, Debug)]
pub struct FiniteViscoContactDriftControls {
    pub material_balance_tolerance_j: f64,
    /// Sum of absolute source/force/endpoint and native kick residuals.
    pub maximum_integration_defect_j: f64,
}

/// Original material update and constant-velocity drift, without fake equilibrium.
#[derive(Clone, Debug)]
pub struct FiniteViscoContactDrift {
    pub duration_s: f64,
    pub time_s: f64,
    pub before: SolidMomentum,
    pub after: SolidMomentum,
    pub linear_residual_kg_m_s: [f64; 3],
    /// Actual change from represented current positions; not external torque.
    pub angular_residual_kg_m2_s: [f64; 3],
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub kinetic_change_j: f64,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Native relaxation loss for this interval; never the split discrepancy.
    pub dissipated_j: f64,
    pub native_loading_work_j: f64,
    pub material_balance_error_j: f64,
    pub internal_endpoint_work_j: f64,
    /// Endpoint Piola work minus actual storage change and native relaxation.
    pub material_endpoint_defect_j: f64,
    pub source_displacement_work_j: f64,
    /// Source kick work minus held-source work on the actual displacement.
    pub source_splitting_work_j: f64,
    /// Minus-old-Piola kick work plus the actual new-Piola endpoint work.
    pub force_splitting_work_j: f64,
    /// Source split + force split - material endpoint defect + native kick
    /// work/decomposition residuals. Contact loss remains a separate account.
    pub integration_defect_j: f64,
    pub absolute_integration_defect_j: f64,
}

impl<'a> FiniteViscoTetMotion<'a> {
    /// Stage h(F_external-F_internal) with the original supported consistent mass.
    /// Held DOFs must be stationary. Positions, time and every material state
    /// remain fixed; the existing impulse owner retains their support reactions.
    #[allow(clippy::too_many_arguments)]
    pub fn force_kick_candidate(
        &self,
        duration_s: f64,
        held: &[usize],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        origin_m: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, FiniteViscoForceKick<'a>), FiniteViscoTetError> {
        checkpoint(cx)?;
        let n = self.mass.ndofs();
        let next_time = self.solid.time_s + duration_s;
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !next_time.is_finite()
            || next_time <= self.solid.time_s
            || body_force_n_m3
                .iter()
                .chain(&origin_m)
                .any(|x| !x.is_finite())
            || held.len() > n
            || nodal_forces_n.len() > controls.maximum_workspace_components
        {
            return Err(invalid(
                "finite progressing force interval and bounded loads required",
            ));
        }
        self.stationary(held, cx)?;
        // Source/internal impulses, force, positions, sparse entries and retained
        // support indices coexist with the original solver's own workspace.
        let extra = n
            .checked_mul(6)
            .and_then(|x| x.checked_add(held.len()))
            .ok_or_else(|| invalid("native split workspace overflow"))?;
        let solve_controls = SolidImpulseControls {
            maximum_workspace_components: controls
                .maximum_workspace_components
                .checked_sub(extra)
                .ok_or_else(|| invalid("native force split exceeds workspace cap"))?,
            ..controls
        };
        let mut accumulated = vec![Dd::ZERO; n];
        for (index, &(dof, value)) in nodal_forces_n.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            if dof >= n || !value.is_finite() {
                return Err(invalid("finite original nodal force components required"));
            }
            accumulated[dof] = accumulated[dof] + Dd::from_f64(value);
        }
        let mut internal = vec![0.0; n];
        for (index, element) in self.solid.elements.iter().enumerate() {
            checkpoint(cx)?;
            let piola = &self.solid.first_piola_pa[index];
            for (local, &node) in element.nodes.iter().enumerate() {
                for axis in 0..3 {
                    accumulated[3 * node + axis] = accumulated[3 * node + axis]
                        + Dd::from_f64(finite(body_force_n_m3[axis] * element.volume / 4.0)?);
                    let force = (0..3).fold(0.0, |sum, j| {
                        piola[3 * axis + j].mul_add(element.gradients[local][j], sum)
                    });
                    internal[3 * node + axis] =
                        finite(element.volume.mul_add(force, internal[3 * node + axis]))?;
                }
            }
        }
        // Round each authored nodal resultant once, before forming impulses.
        // Large opposing entries cannot erase a small net force.
        let mut external = Vec::with_capacity(n);
        for (dof, value) in accumulated.into_iter().enumerate() {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            external.push(finite(value.to_f64())?);
        }
        let mut source = Vec::with_capacity(n);
        let mut impulses = Vec::with_capacity(n);
        for dof in 0..n {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            source.push(finite(duration_s * external[dof])?);
            internal[dof] = finite(-duration_s * internal[dof])?;
            impulses.push((dof, finite(source[dof] + internal[dof])?));
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
            return Err(invalid(
                "accepted Piola force misses internal momentum objectivity",
            ));
        }
        let source_work = mean_work(&self.velocity_mps, &event.motion.velocity_mps, &source, cx)?;
        let internal_work = mean_work(
            &self.velocity_mps,
            &event.motion.velocity_mps,
            &internal,
            cx,
        )?;
        let decomposition = finite(event.impulse.applied_work_j - source_work - internal_work)?;
        energy_gate(
            "force kick work decomposition",
            decomposition,
            controls.work_tolerance_j,
        )?;
        let linear = std::array::from_fn(|a| {
            event.impulse.applied.linear_kg_m_s[a]
                - source_impulse.linear_kg_m_s[a]
                - internal_impulse.linear_kg_m_s[a]
        });
        let angular = std::array::from_fn(|a| {
            event.impulse.applied.angular_kg_m2_s[a]
                - source_impulse.angular_kg_m2_s[a]
                - internal_impulse.angular_kg_m2_s[a]
        });
        if !linear.iter().chain(&angular).all(|x| x.is_finite())
            || magnitude(linear) > controls.linear_tolerance_kg_m_s
            || magnitude(angular) > controls.angular_tolerance_kg_m2_s
        {
            return Err(invalid(
                "force kick misses source/internal momentum decomposition",
            ));
        }
        // A force phase starts a new continuation lineage. Contact impulses
        // retain it; the pre-kick input and later force phases cannot reuse it.
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
        let receipt = FiniteViscoForceKick {
            duration_s,
            impulse: event.impulse,
            source_impulse,
            internal_impulse,
            source_impulse_work_j: source_work,
            internal_impulse_work_j: internal_work,
            work_decomposition_residual_j: decomposition,
            linear_decomposition_residual_kg_m_s: linear,
            angular_decomposition_residual_kg_m2_s: angular,
            stored_before_j: self.solid.stored_energy_j,
            binding,
        };
        checkpoint(cx)?;
        Ok((event.motion, receipt))
    }

    /// Drift with unchanged post-contact velocity, then accept the original
    /// material endpoint and elapsed time. The private kick binding permits only
    /// velocity changes since the kick; other states or mass bindings refuse.
    /// Held support displacements stay fixed. No equilibrium solve is claimed.
    pub fn contact_drift_candidate(
        &self,
        kick: &FiniteViscoForceKick<'a>,
        controls: FiniteViscoContactDriftControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, FiniteViscoContactDrift), FiniteViscoTetError> {
        checkpoint(cx)?;
        let binding = &kick.binding;
        if !Arc::ptr_eq(&self.solid.accepted_identity, &binding.accepted_identity)
            || !Arc::ptr_eq(&self.mass, &binding.mass)
            || !std::ptr::eq(self.solid.nodes_m, binding.reference_nodes_m)
        {
            return Err(invalid(
                "contact drift requires its original kick configuration, mass and material history",
            ));
        }
        if [
            controls.material_balance_tolerance_j,
            controls.maximum_integration_defect_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(invalid(
                "finite nonnegative material and split allowances required",
            ));
        }
        self.stationary(&binding.held, cx)?;
        let time_s = finite(self.solid.time_s + binding.duration_s)?;
        if time_s <= self.solid.time_s {
            return Err(invalid("contact drift cannot advance its clock"));
        }
        let mut displacement = self.solid.displacement_m.to_vec();
        for (node, (u, v)) in displacement.iter_mut().zip(&self.velocity_mps).enumerate() {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            for axis in 0..3 {
                if v[axis] != 0.0 {
                    u[axis] = finite(binding.duration_s.mul_add(v[axis], u[axis]))?;
                }
            }
        }
        let endpoint =
            self.solid
                .material_endpoint(displacement.as_flattened(), binding.duration_s, cx)?;
        let energy = &endpoint.energy;
        energy_gate(
            "native drift material balance",
            energy.material_balance_error_j,
            controls.material_balance_tolerance_j,
        )?;
        if !energy.stored_after_j.is_finite()
            || !energy.dissipated_j.is_finite()
            || energy.dissipated_j < 0.0
        {
            return Err(invalid(
                "finite drift storage and nonnegative native relaxation required",
            ));
        }
        let momentum = self
            .mass
            .constant_velocity_drift(
                self.solid.nodes_m,
                &self.solid.displacement_m,
                &displacement,
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
                "constant-velocity drift misses original-mass momentum conservation",
            ));
        }
        let mut source_displacement_work = Dd::ZERO;
        for (dof, ((new, old), force)) in displacement
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
        energy_gate(
            "absolute force/material splitting work",
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
            native_loading_work_j: energy.native_loading_work_j,
            material_balance_error_j: energy.material_balance_error_j,
            internal_endpoint_work_j: energy.internal_endpoint_work_j,
            material_endpoint_defect_j: energy.endpoint_work_defect_j,
            source_displacement_work_j,
            source_splitting_work_j: source_split,
            force_splitting_work_j: force_split,
            integration_defect_j: integration,
            absolute_integration_defect_j: absolute,
        };
        let solid = FiniteViscoTetSolid {
            nodes_m: self.solid.nodes_m,
            tetrahedra: self.solid.tetrahedra,
            material: self.solid.material,
            elements: Arc::clone(&self.solid.elements),
            budget: self.solid.budget,
            displacement_m: displacement.into(),
            states: endpoint.states.into(),
            first_piola_pa: endpoint.first_piola_pa.into(),
            stored_energy_j: energy.stored_after_j,
            time_s,
            accepted_identity: Arc::new(()),
        };
        let candidate = Self {
            solid,
            velocity_mps: self.velocity_mps.clone(),
            reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
            reference_nodal_density_kg_m3: self
                .reference_nodal_density_kg_m3
                .as_ref()
                .map(Arc::clone),
            mass: Arc::clone(&self.mass),
        };
        checkpoint(cx)?;
        Ok((candidate, receipt))
    }

    fn stationary(&self, held: &[usize], cx: &Cx<'_>) -> Result<(), FiniteViscoTetError> {
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
                    "force/contact/drift requires stationary held velocities",
                ));
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
) -> Result<f64, FiniteViscoTetError> {
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
fn finite(value: f64) -> Result<f64, FiniteViscoTetError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "native force/contact/drift arithmetic is unrepresentable",
        ))
    }
}
fn magnitude(value: [f64; 3]) -> f64 {
    value[0].hypot(value[1]).hypot(value[2])
}
