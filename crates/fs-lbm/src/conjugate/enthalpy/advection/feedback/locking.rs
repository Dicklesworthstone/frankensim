//! Supported casting: native momentum becomes frame reaction and local heat.
//!
//! Freezing by itself does not remove momentum. This opt-in law declares an
//! infinitely massive stationary support and perfectly inelastic capture of
//! newly frozen material. Surviving liquid faces are never projected or clipped.
use std::collections::{BTreeMap, BTreeSet};

use fs_conduction::transient::deposition::{EnergyAudit, StepConfig};
use fs_material::phase::transport::PhaseTransport;

use super::{
    CellId, DepositedCell, DepositionError, DepositionSimulation, MaterialError, MaterialFlowError,
    MaterialFlowFeedback, PhaseFlowTransport, admit_density, check_gate, invalid,
};

/// A declared immobile support and its accepted impulse/heat history.
///
/// The frame receives momentum while its prescribed zero velocity does no
/// mechanical work. Lost fluid kinetic energy heats the captured material;
/// it is not also deposited in a second support-temperature reservoir. Angular
/// impulse is measured about the native grid origin in the same world frame.
#[derive(Debug, Clone, PartialEq)]
pub struct StationaryPhaseSupport {
    cells: BTreeSet<CellId>,
    time_s: f64,
    impulse_n_s: [f64; 3],
    angular_impulse_n_m_s: [f64; 3],
    capture_heat_j: f64,
}

impl StationaryPhaseSupport {
    /// Declare which original material cells can be mechanically captured by
    /// this stationary support. Missing/duplicate/inactive IDs refuse.
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        supported_cells: &[CellId],
    ) -> Result<Self, DepositionError> {
        let cells: BTreeSet<_> = supported_cells.iter().copied().collect();
        if cells.is_empty()
            || cells.len() != supported_cells.len()
            || cells.iter().any(|id| !thermal.cells().contains_key(id))
        {
            return Err(invalid(
                "stationary support requires distinct active material cells",
            ));
        }
        Ok(Self {
            cells,
            time_s: thermal.time_s(),
            impulse_n_s: [0.0; 3],
            angular_impulse_n_m_s: [0.0; 3],
            capture_heat_j: 0.0,
        })
    }

    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn impulse_n_s(&self) -> [f64; 3] {
        self.impulse_n_s
    }
    pub const fn angular_impulse_n_m_s(&self) -> [f64; 3] {
        self.angular_impulse_n_m_s
    }
    pub const fn capture_heat_j(&self) -> f64 {
        self.capture_heat_j
    }
}

/// Explicit original-enthalpy storage allowances for the capture impulse.
#[derive(Debug, Clone, Copy)]
pub struct SupportedFreezingControls {
    pub cell_energy_tolerance_j: f64,
    pub energy_tolerance_j: f64,
}

/// One accepted thermal advance followed by its same-time physical capture.
#[derive(Debug, Clone, PartialEq)]
pub struct SupportedFreezingReport {
    /// Includes native viscous heat AND capture heat exactly once. Both are
    /// external to the thermal subsystem, internal to the full coupled system.
    pub thermal: EnergyAudit,
    /// Native face-mass kinetic energies immediately before/after capture.
    pub initial_capture_kinetic_energy_j: f64,
    pub final_capture_kinetic_energy_j: f64,
    pub capture_heat_j: f64,
    pub captured_faces: usize,
    pub cell_capture_heat_j: BTreeMap<CellId, f64>,
    /// Impulse delivered to the stationary frame, opposite the fluid change.
    pub cell_support_impulse_n_s: BTreeMap<CellId, [f64; 3]>,
    pub support_impulse_n_s: [f64; 3],
    pub support_angular_impulse_n_m_s: [f64; 3],
}

impl<F: PhaseTransport + Clone, R> MaterialFlowFeedback<'_, F, R, (), PhaseFlowTransport>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Advance supported material through a bounded moving-to-solid endpoint.
    ///
    /// Original conduction/advection, latent storage, source material sampling
    /// and unsteady Stokes all run first on their existing accepted intervals.
    /// Only the requested final endpoint may defer its exact-solid constraint.
    /// Earlier moving-face freezing still refuses and requires a smaller caller
    /// interval. No hidden event search, new flow solve or thermal clock is used.
    ///
    /// Each captured native face must join TWO newly frozen supported cells.
    /// A nonzero face touching surviving liquid refuses instead of deleting its
    /// flux or treating pressure-projection loss as physical viscosity. This
    /// admits complete moving components and already zero-flux phase boundaries.
    /// The original half-cell inertial metric partitions the removed momentum
    /// and kinetic energy. The SAME thermal owner deposits the energy, then all
    /// endpoint material queries run again. Any changed phase mask, including
    /// capture-induced remelting, refuses; refine the interval/model explicitly.
    ///
    /// This entry requires native unsteady physical heating and a no-op inner
    /// observer. Other coupled history cannot observe an unheated intermediate
    /// endpoint. Authored thermal boundaries, source-dependent conduction,
    /// buoyancy and the existing external actuator still use their native paths.
    /// Refusal leaves thermal state, flow, support reaction and source histories
    /// unchanged. Ordinary `advance_observed` retains its original strict policy.
    /// Equal phase densities and a fixed skeleton remain mandatory. This is
    /// first-order endpoint splitting, not a front-localization, moving-solid,
    /// shrinkage, elastic-reference or arbitrary partial-freezing model.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_with_supported_freezing<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        support: &StationaryPhaseSupport,
        target_s: f64,
        step: StepConfig,
        controls: SupportedFreezingControls,
        mut check: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (
            DepositionSimulation<'t>,
            Self,
            StationaryPhaseSupport,
            SupportedFreezingReport,
        ),
        MaterialFlowError<DepositionError>,
    > {
        check_gate(self.gate, &mut check)?;
        if !self.exact_phase
            || !target_s.is_finite()
            || target_s <= thermal.time_s()
            || support.time_s.to_bits() != thermal.time_s().to_bits()
            || self.support_capture_endpoint_s.is_some()
            || !self
                .unsteady
                .as_ref()
                .is_some_and(|owner| owner.has_heating())
            || [
                controls.cell_energy_tolerance_j,
                controls.energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid("supported freezing needs synchronized phase, support and native physical-heating states").into());
        }
        self.check_map(thermal, &mut check)?;
        self.check_unsteady_start(thermal, &mut check)?;
        let mapped: BTreeSet<_> = self.input.cell_ids.iter().flatten().copied().collect();
        if !support.cells.is_subset(&mapped) {
            return Err(invalid("support cells must belong to the original phase-flow map").into());
        }
        let supported: Vec<_> = self
            .input
            .cell_ids
            .iter()
            .map(|id| id.is_some_and(|id| support.cells.contains(&id)))
            .collect();
        let mut prepared = self.clone();
        prepared.support_capture_endpoint_s = Some(target_s);
        let (endpoint, mut next, mut thermal_audit) =
            thermal.advance_observed(&prepared, target_s, step, &mut check)?;
        next.support_capture_endpoint_s = None;
        let old = next
            .applied
            .as_ref()
            .ok_or(invalid("capture lacks its actual source momentum epoch"))?;
        let source = next
            .momentum_domain(&old.samples, &mut check)?
            .ok_or(invalid("capture requires an exact phase domain"))?;
        let target = next
            .momentum_domain(&next.endpoint, &mut check)?
            .ok_or(invalid("capture requires an exact phase endpoint"))?;
        let fluid = next.reference.map(|_| next.reference_fluid()).transpose()?;
        let mut viscosity = next
            .reference
            .map(|s| vec![s.dynamic_viscosity_pa_s(); next.domain.cell_count()])
            .unwrap_or_default();
        for sample in &old.samples {
            if let Some(transport) = sample.transport.mobile_transport() {
                *viscosity.get_mut(sample.voxel).ok_or(invalid(
                    "capture lacks the actual mobile material reference",
                ))? = transport.dynamic_viscosity_pa_s();
            }
        }
        let native = next
            .unsteady
            .as_mut()
            .expect("admitted native heating")
            .capture_supported_phase(
                &source,
                &target,
                fluid.as_ref(),
                next.config,
                &viscosity,
                &supported,
                next.gate,
            )?;
        let mut heat = BTreeMap::new();
        let mut impulses = BTreeMap::new();
        let mut impulse = [0.0; 3];
        let mut angular = [0.0; 3];
        for (voxel, id) in next.input.cell_ids.iter().enumerate() {
            check_gate(next.gate, &mut check)?;
            let Some(id) = id else {
                continue;
            };
            if native.cell_heat_j[voxel] != 0.0 {
                heat.insert(*id, native.cell_heat_j[voxel]);
            }
            let p = native.cell_support_impulse_n_s[voxel];
            if p != [0.0; 3] {
                impulses.insert(*id, p);
            }
            let [x, y, z] = next.domain.coords(voxel);
            let r = next.domain.center(x, y, z);
            for axis in 0..3 {
                impulse[axis] += p[axis];
                let (b, c) = ((axis + 1) % 3, (axis + 2) % 3);
                angular[axis] += r[b] * p[c] - r[c] * p[b];
            }
        }
        let (final_thermal, capture_audit) = endpoint.apply_heat_j(
            &heat,
            controls.cell_energy_tolerance_j,
            controls.energy_tolerance_j,
            || check_gate(next.gate, &mut check),
        )?;
        if (capture_audit.external_heat_j - native.heat_j).abs() > controls.energy_tolerance_j {
            return Err(
                invalid("native capture heat differs from the original thermal deposit").into(),
            );
        }
        next.endpoint = next.sample_materials(&final_thermal, &mut check)?;
        admit_density(&next.endpoint, next.density_kg_m3)?;
        next.mushy_drag(&next.endpoint, &mut check)?;
        next.check_phase_conduction(&next.endpoint, &mut check)?;
        let final_domain = next
            .momentum_domain(&next.endpoint, &mut check)?
            .ok_or(invalid("capture lost its exact phase topology"))?;
        for cell in 0..target.cell_count() {
            check_gate(next.gate, &mut check)?;
            if target.is_fluid(cell) != final_domain.is_fluid(cell) {
                return Err(invalid("capture heat changes the final phase topology; refine the supported-freezing interval").into());
            }
        }
        next.check_phase_drives(&final_domain, &mut check)?;
        next.check_unsteady_start(&final_thermal, &mut check)?;
        let mut support = support.clone();
        support.time_s = final_thermal.time_s();
        support.capture_heat_j += native.heat_j;
        for axis in 0..3 {
            support.impulse_n_s[axis] += impulse[axis];
            support.angular_impulse_n_m_s[axis] += angular[axis];
        }
        thermal_audit.external_heat_j += capture_audit.external_heat_j;
        thermal_audit.stored_energy_change_j += capture_audit.stored_energy_change_j;
        thermal_audit.residual_j += capture_audit.residual_j;
        if impulse
            .iter()
            .chain(&angular)
            .chain(&support.impulse_n_s)
            .chain(&support.angular_impulse_n_m_s)
            .any(|v| !v.is_finite())
            || ![
                support.capture_heat_j,
                thermal_audit.external_heat_j,
                thermal_audit.stored_energy_change_j,
                thermal_audit.residual_j,
            ]
            .iter()
            .all(|v| v.is_finite())
            || thermal_audit.residual_j.abs()
                > step.energy_tolerance_j + controls.energy_tolerance_j
        {
            return Err(DepositionError::Arithmetic.into());
        }
        check_gate(next.gate, &mut check)?;
        let report = SupportedFreezingReport {
            thermal: thermal_audit,
            initial_capture_kinetic_energy_j: native.initial_kinetic_energy_j,
            final_capture_kinetic_energy_j: native.final_kinetic_energy_j,
            capture_heat_j: native.heat_j,
            captured_faces: native.captured_faces,
            cell_capture_heat_j: heat,
            cell_support_impulse_n_s: impulses,
            support_impulse_n_s: impulse,
            support_angular_impulse_n_m_s: angular,
        };
        Ok((final_thermal, next, support, report))
    }
}
