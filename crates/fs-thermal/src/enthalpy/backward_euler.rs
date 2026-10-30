//! Single-step backward-Euler residual over the existing thermal network.
//!
//! Unknowns are changes in CELL energy [J], in stable CellId order. For each
//! cell the public residual is m*(h_new-h_old) - dt*Q(h_new). Storage uses the
//! original mass and enthalpy law, including latent plateaus. Q is evaluated by
//! the SAME heat-rate routine used for explicit stepping. The caller supplies
//! a nonlinear solver; this module neither guesses derivatives nor accepts a
//! small iteration update in place of an actual local energy residual.
//!
//! Fixed topology/geometry, immutable caloric laws, endpoint-evaluated sources
//! and conductances. Optional entropy accounting uses the original caloric laws.
//! No moving mass, kinetics or unconditional
//! nonlinear-convergence claim. Refine time steps independently of tolerances.

use super::*;

impl DepositionSimulation<'_> {
    /// Evaluate the ORIGINAL thermal rate owner's external power at this state.
    /// Trial sources add to authored sources; link replacements keep their
    /// original coverage. Internal link transfers are excluded from the result.
    /// No time, enthalpy, event cursor or accepted coefficient is changed.
    ///
    /// Coupled consumers can pass only their independent external sources here,
    /// omitting internal Joule/gas transfers, to measure supply without deriving
    /// it by subtraction from an already combined energy account. This is one
    /// bounded rate evaluation, not a thermal step or an energy certificate.
    pub fn external_step_power_w(
        &self,
        sources_w: &BTreeMap<CellId, f64>,
        conductances_w_k: &BTreeMap<u64, f64>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, DepositionError> {
        checkpoint()?;
        let mut trial = self.clone();
        trial.set_step_conductances(conductances_w_k, &mut checkpoint)?;
        trial.add_step_sources(sources_w, &mut checkpoint)?;
        let power = trial.heat_rates(f64::MAX, &mut checkpoint)?.external_w;
        checkpoint()?;
        finite(power)
    }

    /// Check a coupled material binding against the timeline's actual borrowed
    /// caloric object, not merely one matching temperature sample or numeric ID.
    /// This compares in-process data addresses (not unstable trait vtable
    /// addresses); it is not portable content identity or material validation.
    /// Cloned timeline candidates keep the same original borrowed laws.
    pub fn uses_material_law(&self, id: MaterialId, law: &dyn ThermalLaw) -> bool {
        self.materials.get(&id).is_some_and(|actual| {
            std::ptr::addr_eq(
                *actual as *const dyn ThermalLaw,
                law as *const dyn ThermalLaw,
            )
        })
    }
}

/// Endpoint-only additions and replacements, not persistent authored inputs.
#[derive(Debug)]
pub struct StepCoefficients<T> {
    /// Additional signed cell powers [W]. Omitted cells receive zero.
    pub sources_w: BTreeMap<CellId, f64>,
    /// Replacement total link conductances [W/K]. Omitted links retain theirs.
    pub conductances_w_k: BTreeMap<u64, f64>,
    /// The exact constitutive/current evaluation that supplied these fields.
    pub payload: T,
}

/// One frozen initial state and time interval for an externally solved residual.
pub struct BackwardEulerStep<'s, 'm> {
    initial: &'s DepositionSimulation<'m>,
    target_s: f64,
    duration_s: f64,
    ids: Vec<CellId>,
}

impl<'s, 'm> BackwardEulerStep<'s, 'm> {
    /// Prepare an interval without consuming events or changing initial state.
    /// Process start-time events first with `advance_to(current_time, ...)`.
    /// An event strictly inside the interval refuses. Events exactly at target
    /// remain pending in an accepted left-limit state; process them explicitly
    /// AFTER all coupled owners have accepted this interval. The next interval
    /// refuses until those events have been processed.
    pub fn new(
        initial: &'s DepositionSimulation<'m>,
        target_s: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        let duration_s = target_s - initial.time_s;
        if !target_s.is_finite() || !positive(duration_s) || initial.cells.is_empty() {
            return Err(DepositionError::InvalidInput(
                "implicit step needs active cells and a later finite target",
            ));
        }
        if initial.budget.max_substeps == 0 {
            return Err(DepositionError::Budget("thermal substeps"));
        }
        if initial
            .schedule
            .get(initial.cursor)
            .is_some_and(|e| e.time_s < target_s)
        {
            return Err(DepositionError::InvalidInput(
                "process pending events or split implicit step at the next event",
            ));
        }
        let mut ids = Vec::with_capacity(initial.cells.len());
        for (i, &id) in initial.cells.keys().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            ids.push(id);
        }
        checkpoint()?;
        Ok(Self {
            initial,
            target_s,
            duration_s,
            ids,
        })
    }

    /// Stable ordering of energy-increment unknowns and residual components.
    pub fn cell_ids(&self) -> &[CellId] {
        &self.ids
    }
    /// Actual physical duration, not an explicit CFL-limited approximation.
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    /// Physical endpoint at which coefficients and temperatures are evaluated.
    pub const fn target_s(&self) -> f64 {
        self.target_s
    }

    /// Evaluate a complete trial without publishing ANY physical history.
    /// The callback sees one endpoint state and must not mutate external physical
    /// history. It may retain immutable source payloads; no Clone bound is needed.
    /// Numerical solver work must have its own evaluation/iteration budget.
    /// Unknowns are joules added to each initial cell, NOT Kelvin or J/kg.
    /// Every material query and all signed source/link admission remain active.
    pub fn evaluate<T, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(&DepositionSimulation<'m>) -> Result<StepCoefficients<T>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, T>, E> {
        self.evaluate_entropy_option(
            energy_increments_j,
            |state| resolve(state).map(|coefficients| (coefficients, BTreeMap::new())),
            false,
            checkpoint,
        )
    }

    /// Evaluate the SAME endpoint residual with exact caloric entropy storage
    /// and endpoint-quadrature conductive/boundary entropy. Every active law
    /// must provide its entropy integral. Numerical entropy remains a separate
    /// defect; applied sources remain unclassified until their owner supplies
    /// the missing heat/work and boundary-entropy semantics.
    /// Use this for a final candidate when iterative entropy queries are unwanted.
    pub fn evaluate_with_entropy<T, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(&DepositionSimulation<'m>) -> Result<StepCoefficients<T>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, T>, E> {
        self.evaluate_entropy_option(
            energy_increments_j,
            |state| resolve(state).map(|coefficients| (coefficients, BTreeMap::new())),
            true,
            checkpoint,
        )
    }

    /// Resolve endpoint link/source coefficients and surface-reservoir
    /// conductances together. Reservoir keys name active cells with an existing
    /// authored bath; its actual temperature and covered area are retained.
    /// The returned map replaces total W/K only for this endpoint evaluation,
    /// allowing a nonlinear half-cell or film resistance to use the SAME
    /// original residual as internal conduction. Accepted states retain authored
    /// boundaries. Missing/nonfinite/negative coefficients refuse transactionally.
    pub fn evaluate_with_reservoirs<T, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
        ) -> Result<(StepCoefficients<T>, BTreeMap<CellId, f64>), E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, T>, E> {
        self.evaluate_entropy_option(energy_increments_j, resolve, false, checkpoint)
    }

    /// The same endpoint evaluation with exact caloric entropy storage and
    /// reservoir entropy from the ACTUAL applied conductance and bath
    /// temperature. Numerical temporal entropy remains separately reported.
    pub fn evaluate_with_entropy_and_reservoirs<T, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
        ) -> Result<(StepCoefficients<T>, BTreeMap<CellId, f64>), E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, T>, E> {
        self.evaluate_entropy_option(energy_increments_j, resolve, true, checkpoint)
    }

    fn evaluate_entropy_option<T, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
        ) -> Result<(StepCoefficients<T>, BTreeMap<CellId, f64>), E>,
        with_entropy: bool,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, T>, E> {
        checkpoint()?;
        if energy_increments_j.len() != self.ids.len() {
            return Err(DepositionError::InvalidInput(
                "implicit energy vector length differs from active cells",
            )
            .into());
        }
        let mut next = self.initial.clone();
        next.time_s = self.target_s;
        let mut storage = Vec::with_capacity(self.ids.len());
        let mut stored_change = 0.0;
        let mut initial_h = with_entropy.then(BTreeMap::new);
        for (i, (&id, &joules)) in self.ids.iter().zip(energy_increments_j).enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            finite(joules)?;
            let old = &self.initial.cells[&id];
            if let Some(initial_h) = &mut initial_h {
                initial_h.insert(id, old.specific_h_j_kg);
            }
            let h = finite(old.specific_h_j_kg + finite(joules / old.mass_kg)?)?;
            // Refuse a source too small to change its actual stored coordinate.
            if joules != 0.0 && h == old.specific_h_j_kg {
                return Err(DepositionError::Arithmetic.into());
            }
            next.cells
                .get_mut(&id)
                .expect("frozen cell map")
                .specific_h_j_kg = h;
            next.temperature_k(id)?;
            let actual = finite(old.mass_kg * finite(h - old.specific_h_j_kg)?)?;
            stored_change = finite(stored_change + actual)?;
            storage.push(actual);
        }
        checkpoint()?;
        let (
            StepCoefficients {
                sources_w,
                conductances_w_k,
                payload,
            },
            reservoir_conductances,
        ) = resolve(&next)?;
        checkpoint()?;
        let links = next.set_step_conductances(&conductances_w_k, &mut checkpoint)?;
        let reservoirs =
            next.set_step_reservoir_conductances(&reservoir_conductances, &mut checkpoint)?;
        let sources = next.add_step_sources(&sources_w, &mut checkpoint)?;
        let rates = next.heat_rates(f64::MAX, &mut checkpoint)?;
        let entropy = if let Some(initial_h) = initial_h {
            Some(next.entropy_step_audit(
                &initial_h,
                &rates.temperatures,
                self.duration_s,
                &mut checkpoint,
            )?)
        } else {
            None
        };
        // The endpoint carries authored topology/boundaries, never trial inputs.
        for (id, value) in links {
            next.links
                .get_mut(&id)
                .expect("validated link")
                .conductance_w_k = value;
        }
        for (id, value) in reservoirs {
            next.cells
                .get_mut(&id)
                .expect("validated reservoir cell")
                .boundary
                .reservoir
                .as_mut()
                .expect("validated reservoir")
                .conductance_w_k = value;
        }
        for (id, value) in sources {
            next.cells
                .get_mut(&id)
                .expect("validated source")
                .boundary
                .source_w = value;
        }
        let mut residual = Vec::with_capacity(self.ids.len());
        for (i, (&id, &joules)) in self.ids.iter().zip(&storage).enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            residual.push(finite(
                joules - finite(self.duration_s * rates.watts[&id])?,
            )?);
        }
        let external_heat_j = finite(self.duration_s * rates.external_w)?;
        let audit = EnergyAudit {
            incoming_mass_kg: 0.0,
            incoming_energy_j: 0.0,
            external_heat_j,
            stored_energy_change_j: stored_change,
            residual_j: finite(stored_change - external_heat_j)?,
            substeps: 1,
            events: 0,
        };
        checkpoint()?;
        Ok(BackwardEulerTrial {
            next,
            payload,
            residual,
            audit,
            entropy,
        })
    }
}

/// A tentative physical endpoint together with the residual of THAT endpoint.
/// Private fields prevent unrelated state, source and audit substitution.
pub struct BackwardEulerTrial<'m, T> {
    next: DepositionSimulation<'m>,
    payload: T,
    residual: Vec<f64>,
    audit: EnergyAudit,
    entropy: Option<super::entropy::EntropyAudit>,
}
impl<'m, T> BackwardEulerTrial<'m, T> {
    /// Per-cell energy-balance residual [J] in the step's CellId order.
    pub fn residual_j(&self) -> &[f64] {
        &self.residual
    }
    /// Largest absolute local energy defect; total cancellation cannot hide it.
    pub fn maximum_cell_residual_j(&self) -> f64 {
        self.residual
            .iter()
            .copied()
            .map(f64::abs)
            .fold(0.0, f64::max)
    }
    /// Original endpoint evaluation, still tentative until acceptance succeeds.
    pub const fn payload(&self) -> &T {
        &self.payload
    }
    /// Total energy account from the actual material-backed storage differences.
    pub const fn audit(&self) -> &EnergyAudit {
        &self.audit
    }

    /// Optional entropy account of this exact tentative endpoint and its
    /// applied coefficients. It becomes an accepted account only if `accept`
    /// succeeds. Ordinary `evaluate` performs no entropy integration.
    pub const fn entropy_audit(&self) -> Option<&super::entropy::EntropyAudit> {
        self.entropy.as_ref()
    }

    /// Publish a candidate only after independent LOCAL and TOTAL joule gates.
    /// This is an algebraic acceptance gate, not a timestep-accuracy certificate.
    /// Accepted endpoint events are still pending, as described by `new`.
    pub fn accept(
        self,
        cell_tolerance_j: f64,
        total_tolerance_j: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, T, EnergyAudit), DepositionError> {
        checkpoint()?;
        if !nonnegative(cell_tolerance_j) || !nonnegative(total_tolerance_j) {
            return Err(DepositionError::InvalidInput(
                "implicit tolerances must be finite nonnegative joules",
            ));
        }
        for (i, &residual) in self.residual.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            if residual.abs() > cell_tolerance_j {
                return Err(DepositionError::EnergyBalance(residual));
            }
        }
        if self.audit.residual_j.abs() > total_tolerance_j {
            return Err(DepositionError::EnergyBalance(self.audit.residual_j));
        }
        checkpoint()?;
        Ok((self.next, self.payload, self.audit))
    }
}

#[cfg(test)]
#[path = "backward_euler/tests.rs"]
mod tests;
