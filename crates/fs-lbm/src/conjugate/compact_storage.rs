//! One physical sensible-energy owner per collapsed compact junction.

use fs_exec::CancelGate;
use fs_material::thermal::transport::storage::{PreparedHeatStorage, TabulatedHeatStorage};

use super::{ChtError, EnergySolution, ThermalSetup, VoxelDomain, finite, finite_positive, poll};

/// Temperature-dependent total storage assigned to one compact component.
#[derive(Debug, Clone, PartialEq)]
pub struct CompactHeatStorage {
    /// Index in [`ThermalSetup::compact_components`].
    pub component: usize,
    /// Total junction storage; occupied voxel aliases add no heat capacity.
    pub law: CompactHeatStorageLaw,
}

/// Explicit total sensible capacity, separate from the component's resistors.
///
/// Both variants integrate the original breakpoints exactly. Neither infers
/// latent heat, extrapolates, nor borrows storage from an occupied voxel.
#[derive(Debug, Clone, PartialEq)]
pub enum CompactHeatStorageLaw {
    /// Caller-authored total capacity C(T), in J/K.
    Tabulated(TabulatedHeatStorage),
    /// Prepared source-qualified constituents with fixed masses or volumes.
    Material(PreparedHeatStorage),
}

impl CompactHeatStorageLaw {
    /// Inclusive temperature support, in absolute kelvin.
    #[must_use]
    pub fn temperature_span_k(&self) -> [f64; 2] {
        match self {
            Self::Tabulated(law) => law.temperature_span_k(),
            Self::Material(law) => law.temperature_span_k(),
        }
    }

    /// Total capacity at a supported junction temperature, J/K.
    pub fn heat_capacity_j_k(&self, temperature_k: f64) -> Result<f64, ChtError> {
        match self {
            Self::Tabulated(law) => law.heat_capacity_j_k(temperature_k),
            Self::Material(law) => law.heat_capacity_j_k(temperature_k),
        }
        .map_err(|error| invalid(error.to_string()))
    }

    /// Signed stored-energy change across all crossed breakpoints, J.
    pub fn energy_change_j(&self, initial_k: f64, final_k: f64) -> Result<f64, ChtError> {
        match self {
            Self::Tabulated(law) => law.energy_change_j(initial_k, final_k),
            Self::Material(law) => law.energy_change_j(initial_k, final_k),
        }
        .map_err(|error| invalid(error.to_string()))
    }

    fn is_constant(&self, gate: &CancelGate) -> Result<bool, ChtError> {
        let first = self.heat_capacity_j_k(self.temperature_span_k()[0])?;
        match self {
            Self::Tabulated(law) => {
                for (i, knot) in law.knots().iter().enumerate() {
                    if i.is_multiple_of(256) {
                        poll(gate)?;
                    }
                    if knot[1] != first {
                        return Ok(false);
                    }
                }
            }
            Self::Material(law) => {
                for (i, t) in law.temperatures_k().enumerate() {
                    if i.is_multiple_of(256) {
                        poll(gate)?;
                    }
                    if law
                        .heat_capacity_j_k(t)
                        .map_err(|e| invalid(e.to_string()))?
                        != first
                    {
                        return Ok(false);
                    }
                }
            }
        }
        Ok(true)
    }
}

/// Endpoint convergence controls for physical compact-junction storage.
///
/// The original time level stays fixed throughout every trial. Acceptance
/// requires each junction's exact energy balance; unrelated junction errors
/// cannot cancel. These controls do not change the physical time step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CompactStorageControls {
    /// Maximum shared storage/radiation endpoint passes per physical step.
    pub maximum_iterations: usize,
    /// Maximum unrelaxed temperature change for nonconstant laws, K.
    /// Constant tables need only their exact energy and support checks.
    pub temperature_tolerance_k: f64,
    /// Absolute per-junction sensible-energy tolerance, J.
    pub energy_absolute_tolerance_j: f64,
    /// Relative per-junction energy tolerance on physical transferred heat.
    pub energy_relative_tolerance: f64,
    /// Upper bound on trial relaxation, in (0, 1]. Capacity and support
    /// safeguards may reduce it; returned solutions are never clamped.
    pub temperature_relaxation: f64,
}

impl Default for CompactStorageControls {
    fn default() -> Self {
        Self {
            maximum_iterations: 100,
            temperature_tolerance_k: 1e-9,
            energy_absolute_tolerance_j: 1e-9,
            energy_relative_tolerance: 1e-10,
            temperature_relaxation: 1.0,
        }
    }
}

impl CompactStorageControls {
    fn validate(self) -> Result<(), ChtError> {
        if self.maximum_iterations == 0
            || !self.temperature_tolerance_k.is_finite()
            || self.temperature_tolerance_k < 0.0
            || !self.energy_absolute_tolerance_j.is_finite()
            || self.energy_absolute_tolerance_j < 0.0
            || !self.energy_relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.energy_relative_tolerance)
            || (self.energy_absolute_tolerance_j == 0.0 && self.energy_relative_tolerance == 0.0)
            || !self.temperature_relaxation.is_finite()
            || self.temperature_relaxation <= 0.0
            || self.temperature_relaxation > 1.0
        {
            return Err(ChtError::InvalidInput {
                field: "thermal.compact_storage_controls",
                reason: "needs a positive iteration budget, finite nonnegative tolerances, a nonzero energy tolerance, relative tolerance below one and relaxation in (0, 1]".into(),
            });
        }
        Ok(())
    }
}

fn invalid(reason: impl Into<String>) -> ChtError {
    ChtError::InvalidInput {
        field: "thermal.compact_heat_storage",
        reason: reason.into(),
    }
}

/// Validate ownership even for a steady solve, without querying its temperature.
/// Empty input avoids allocating an index for the original constant-C path.
pub(super) fn bindings<'a>(
    setup: &'a ThermalSetup,
    gate: &CancelGate,
) -> Result<Vec<Option<&'a CompactHeatStorageLaw>>, ChtError> {
    if setup.compact_heat_storage.is_empty() {
        return Ok(Vec::new());
    }
    setup.compact_storage_controls.validate()?;
    let parts = setup.compact_components.len();
    if setup.compact_heat_storage.len() > parts {
        return Err(invalid("more storage bindings than compact components"));
    }
    let mut laws = vec![None; parts];
    for (i, binding) in setup.compact_heat_storage.iter().enumerate() {
        if i.is_multiple_of(256) {
            poll(gate)?;
        }
        let j = binding.component;
        if j >= parts {
            return Err(invalid(format!("component {j} does not exist")));
        }
        if setup.compact_components[j].heat_capacity_j_k.is_some() {
            return Err(invalid(format!(
                "component {j}: constant capacity and a storage law are competing owners"
            )));
        }
        if laws[j].replace(&binding.law).is_some() {
            return Err(invalid(format!("component {j}: duplicate storage law")));
        }
    }
    Ok(laws)
}

struct JunctionStorage<'a> {
    component: usize,
    cell: usize,
    law: &'a CompactHeatStorageLaw,
    previous_k: f64,
    trial_k: f64,
    lo: [usize; 3],
    hi: [usize; 3],
    secant_j_k: f64,
    constant: bool,
}

/// Positive enthalpy secants preserve the existing conservative M-matrix.
pub(super) struct CompactStorageIteration<'a> {
    entries: Vec<JunctionStorage<'a>>,
    controls: CompactStorageControls,
    time_step_s: f64,
}

#[derive(Clone, Copy)]
pub(super) struct StorageResidual {
    pub normalized: f64,
    pub constitutive: f64,
    pub physical: f64,
}

impl<'a> CompactStorageIteration<'a> {
    pub(super) fn new(
        domain: &VoxelDomain,
        setup: &'a ThermalSetup,
        previous: &[f64],
        time_step_s: f64,
        gate: &CancelGate,
    ) -> Result<Option<Self>, ChtError> {
        let laws = bindings(setup, gate)?;
        if laws.is_empty() {
            return Ok(None);
        }
        let mut entries = Vec::with_capacity(setup.compact_heat_storage.len());
        for (component, law) in laws.into_iter().enumerate() {
            let Some(law) = law else {
                continue;
            };
            poll(gate)?;
            let part = &setup.compact_components[component];
            if (0..3).any(|axis| part.lo[axis] >= domain.dims()[axis]) {
                return Err(invalid(format!(
                    "component {component}: junction lies outside the domain"
                )));
            }
            let cell = domain.index(part.lo[0], part.lo[1], part.lo[2]);
            let previous_k = previous[cell];
            entries.push(JunctionStorage {
                component,
                cell,
                law,
                previous_k,
                trial_k: previous_k,
                lo: part.lo,
                hi: part.hi,
                secant_j_k: law.heat_capacity_j_k(previous_k)?,
                constant: law.is_constant(gate)?,
            });
        }
        Ok(Some(Self {
            entries,
            controls: setup.compact_storage_controls,
            time_step_s,
        }))
    }

    pub(super) fn maximum_iterations(&self) -> usize {
        self.controls.maximum_iterations
    }

    pub(super) fn write_rates(
        &mut self,
        trial: &[f64],
        junction_rates: &mut [f64],
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        for entry in &mut self.entries {
            poll(gate)?;
            let t = trial[entry.cell];
            let change = t - entry.previous_k;
            let capacity = if change == 0.0 {
                entry.law.heat_capacity_j_k(t)?
            } else {
                entry.law.energy_change_j(entry.previous_k, t)? / change
            };
            let rate = capacity / self.time_step_s;
            finite_positive("thermal.compact_heat_storage.secant_j_k", capacity)?;
            finite_positive("thermal.compact_heat_storage.rate_w_k", rate)?;
            finite(
                "thermal.compact_heat_storage.history_w",
                rate * entry.previous_k,
            )?;
            entry.trial_k = t;
            entry.secant_j_k = capacity;
            junction_rates[entry.component] = rate;
        }
        Ok(())
    }

    pub(super) fn supports(
        &self,
        temperature: &[f64],
        gate: &CancelGate,
    ) -> Result<bool, ChtError> {
        for entry in &self.entries {
            poll(gate)?;
            let [lo, hi] = entry.law.temperature_span_k();
            let t = temperature[entry.cell];
            finite("thermal.compact_heat_storage.candidate_k", t)?;
            if !(lo..=hi).contains(&t) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    pub(super) fn residual(
        &self,
        setup: &ThermalSetup,
        trial: &[f64],
        next: &EnergySolution,
        gate: &CancelGate,
    ) -> Result<StorageResidual, ChtError> {
        let mut result = StorageResidual {
            normalized: 0.0,
            constitutive: 0.0,
            physical: 0.0,
        };
        for entry in &self.entries {
            poll(gate)?;
            let j = entry.component;
            let junction = &next.junctions[j];
            let t = junction.temperature_k;
            let stored = entry.law.energy_change_j(entry.previous_k, t)?;
            let linearized = entry.secant_j_k * (t - entry.previous_k);
            let source = setup.compact_components[j].power_w;
            let supplied = self.time_step_s * (source - junction.board_w - junction.case_w);
            let constitutive = stored - linearized;
            let physical = stored - supplied;
            finite(
                "thermal.compact_heat_storage.constitutive_defect_j",
                constitutive,
            )?;
            finite("thermal.compact_heat_storage.physical_defect_j", physical)?;
            let mut flux_scale = source.abs();
            let mut flux_roundoff = source.abs();
            for &(cell, g) in &junction.links {
                flux_scale += g * (t - next.temperature[cell]).abs();
                flux_roundoff += g * (t.abs() + next.temperature[cell].abs());
            }
            let capacity = entry.law.heat_capacity_j_k(t)?.max(entry.secant_j_k);
            let scale = stored
                .abs()
                .max(linearized.abs())
                .max(self.time_step_s * flux_scale);
            // Absolute kelvin enters only the explicit floating-point
            // resolution allowance, never the physical relative tolerance.
            let roundoff = 64.0
                * f64::EPSILON
                * (capacity * t.abs().max(entry.previous_k.abs())
                    + self.time_step_s * flux_roundoff);
            let tolerance = self.controls.energy_absolute_tolerance_j
                + self.controls.energy_relative_tolerance * scale
                + roundoff;
            finite("thermal.compact_heat_storage.energy_scale_j", scale)?;
            finite("thermal.compact_heat_storage.energy_tolerance_j", tolerance)?;
            let c = normalized(constitutive.abs(), tolerance);
            let p = normalized(physical.abs(), tolerance);
            result.constitutive = result.constitutive.max(c);
            result.physical = result.physical.max(p);
            result.normalized = result.normalized.max(c).max(p);
            if !entry.constant {
                let change = (t - trial[entry.cell]).abs();
                let allowance = self.controls.temperature_tolerance_k
                    + 8.0 * f64::EPSILON * t.abs().max(trial[entry.cell].abs());
                result.normalized = result.normalized.max(normalized(change, allowance));
            }
        }
        Ok(result)
    }

    /// Advance only the constitutive trial. Return whether it is exactly the
    /// algebraic candidate, so radiation can reuse its endpoint evaluation.
    pub(super) fn advance_trial(
        &self,
        domain: &VoxelDomain,
        trial: &mut [f64],
        candidate: &[f64],
        gate: &CancelGate,
    ) -> Result<bool, ChtError> {
        let mut fraction = self.controls.temperature_relaxation;
        for entry in &self.entries {
            poll(gate)?;
            let target = candidate[entry.cell];
            let [lo, hi] = entry.law.temperature_span_k();
            let current_capacity = entry.law.heat_capacity_j_k(entry.trial_k)?;
            // Both ends matter: descending from a sharp high-C trial must
            // not leap into a low-C region without derivative damping.
            fraction = fraction.min(entry.secant_j_k / current_capacity);
            if (lo..=hi).contains(&target) {
                fraction = fraction.min(entry.secant_j_k / entry.law.heat_capacity_j_k(target)?);
            }
        }
        for (i, (&old, &target)) in trial.iter().zip(candidate).enumerate() {
            if i.is_multiple_of(4096) {
                poll(gate)?;
            }
            finite("thermal.compact_heat_storage.candidate_k", target)?;
            if target <= 0.0 {
                fraction = fraction.min(0.5 * (old / (old - target)));
            }
        }
        finite_positive("thermal.compact_heat_storage.trial_fraction", fraction)?;
        for (i, (value, &target)) in trial.iter_mut().zip(candidate).enumerate() {
            if i.is_multiple_of(4096) {
                poll(gate)?;
            }
            let next = if fraction == 1.0 {
                target
            } else {
                fraction.mul_add(target - *value, *value)
            };
            finite_positive("thermal.compact_heat_storage.trial_k", next)?;
            *value = next;
        }
        let mut exact_candidate = fraction == 1.0;
        for entry in &self.entries {
            poll(gate)?;
            let target = candidate[entry.cell];
            let [lo, hi] = entry.law.temperature_span_k();
            if (lo..=hi).contains(&target) {
                continue;
            }
            let old = entry.trial_k;
            let bound = if target > hi { hi } else { lo };
            // Limit this constitutive junction only. A junction already at
            // its bound may stay there while coupled neighbors find a valid
            // endpoint; a zero global step would incorrectly refuse it.
            let local_fraction = fraction.min(0.5 * ((bound - old) / (target - old)));
            if !local_fraction.is_finite() || local_fraction < 0.0 {
                return Err(invalid(
                    "no representable support-preserving junction trial",
                ));
            }
            let temperature = local_fraction.mul_add(target - old, old);
            if !(lo..=hi).contains(&temperature) {
                return Err(invalid(
                    "no representable trial inside the junction storage span",
                ));
            }
            for z in entry.lo[2]..entry.hi[2] {
                for y in entry.lo[1]..entry.hi[1] {
                    for x in entry.lo[0]..entry.hi[0] {
                        let cell = domain.index(x, y, z);
                        if cell.is_multiple_of(4096) {
                            poll(gate)?;
                        }
                        trial[cell] = temperature;
                    }
                }
            }
            exact_candidate = false;
        }
        if !self.supports(trial, gate)? {
            return Err(invalid(
                "no representable trial inside all junction storage spans",
            ));
        }
        // Stationary junction trials are permitted: other junctions or a
        // tighter linear solve can still resolve their physical heat balance.
        Ok(exact_candidate)
    }
}

fn normalized(error: f64, tolerance: f64) -> f64 {
    if error == 0.0 { 0.0 } else { error / tolerance }
}
