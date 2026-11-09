//! Conservative caloric transport and self-heating on the native periodic jet.
//!
//! The original isochoric Cv(T) law supplies both integration and inversion at
//! the fixed material density. Eulerian cell masses rho*pi*A*dz evolve with the
//! SAME face q as native mass conservation. One coupled explicit SSPRK3 evolves
//! area, momentum and extensive internal energy; every actual stage resolves
//! temperature and the original Carreau stress again. This is third order in
//! time for smooth admitted semidiscrete solutions, not a held-temperature split.
//!
//! Energy advection is donor cell (first order in space). Optional axial
//! Fourier conduction uses the original source k(T) and the two half-cell
//! resistances dz/(2*pi*k*A). The free surface is insulated; density, composition
//! and surface tension stay fixed. Incompressibility includes the radial strain,
//! so material pressure work vanishes (div u=0); capillary work remains in the
//! native mechanical account. There is no thermal expansion, Marangoni stress,
//! phase change, pressure trajectory, radial temperature profile or implicit
//! thermal integration. The native radius/slope/resolution limits still apply.
//!
//! A fixed numerical reference temperature centers the transported internal
//! energy. Original caloric integrals and inverses retain that reference without
//! subtracting a large source datum. Changing the source's energy zero cannot
//! change transport or heating. Numerical momentum loss is never deposited.

use fs_conduction::enthalpy::ThermalLaw;
use fs_material::thermal::transport::ThermalTransportError;
use fs_material::thermal::transport::isochoric::IsochoricHeatCapacityLaw;
use fs_material::thermal::transport::isochoric::conduction::IsochoricConductionLaw;

use super::*;

/// Explicit flow bounds plus an absolute closed-domain caloric acceptance limit.
#[derive(Clone, Copy, Debug)]
pub struct ThermalJetStepBudget {
    pub flow: SlenderJetStepBudget,
    /// |change in stored internal energy - physical viscous heat| [J].
    /// This is an algebraic conservation allowance, not a time-error estimate.
    pub energy_tolerance_j: f64,
}

/// Source or numerical refusal leaves the entire accepted coupled state intact.
#[derive(Debug)]
pub enum ThermalJetError {
    InvalidInput(&'static str),
    Flow(NonlinearJetError),
    Thermal {
        cell: usize,
        source: ThermalTransportError,
    },
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl From<NonlinearJetError> for ThermalJetError {
    fn from(error: NonlinearJetError) -> Self {
        Self::Flow(error)
    }
}
impl core::fmt::Display for ThermalJetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "conservative thermal jet: {self:?}")
    }
}
impl std::error::Error for ThermalJetError {}

/// Native mechanical history and the same interval's caloric account.
#[derive(Clone, Copy, Debug)]
pub struct ThermalJetAdvance {
    pub flow: SlenderJetAdvance,
    /// Change of sum(m*[u(T)-u(T_reference)]) [J]. The domain is closed;
    /// using this fixed reference removes dependence on the source energy zero.
    pub internal_energy_change_j: f64,
    /// Stored change minus the native physical viscous heat [J].
    pub energy_balance_residual_j: f64,
    /// Sum of actual cell mass changes [kg], retained separately from energy.
    pub mass_balance_residual_kg: f64,
    /// Delta(kinetic + surface + internal) plus numerical momentum loss [J].
    /// Neither this truncation residual nor upwind loss becomes caloric heat.
    pub combined_energy_residual_j: f64,
}

struct Source<'m> {
    caloric: &'m IsochoricHeatCapacityLaw,
    conduction: Option<&'m IsochoricConductionLaw<'m>>,
    reference_temperature_k: f64,
    maximum_temperature_slope: f64,
}

struct Workspace {
    // Actual-stage material only. The accepted flow profile is untouched until
    // publication, including after source refusal or a rejected geometry stage.
    profile: CarreauProfile,
    specific: Vec<f64>,
    conductivity: Vec<f64>,
    conductance: Vec<f64>,
}

/// Original source-backed calorics and generalized-Newtonian periodic flow.
///
/// The source density must match the native incompressible fluid. The caller
/// binds its Cv and rheology; neither viscosity nor a material name supplies
/// calorimetry. An optional conductivity owner must retain this EXACT caloric
/// identity. All numerical arrays are prepared at construction; warm advances
/// reuse O(N) storage. Source preparation and caloric inversion run per stage.
pub struct ThermalCarreauJet<'m> {
    flow: PeriodicSlenderJet,
    source: Source<'m>,
    accepted: Vec<f64>,
    candidate: Vec<f64>,
    stage: Vec<f64>,
    rhs: Vec<f64>,
    workspace: Workspace,
}

impl<'m> ThermalCarreauJet<'m> {
    /// Bind evolving cell temperatures to an original explicit Carreau jet.
    ///
    /// Temperatures replace any previously prescribed field. The native flow
    /// and its clock are retained. The rate ceiling is a positive physical
    /// invariant-rate bound; every evolving temperature must support the whole
    /// original interval from rest through this ceiling. No source is clipped.
    ///
    /// None conduction means thermally insulated axial faces, not a zero
    /// conductivity inferred from viscosity. An implicit-viscosity flow policy
    /// refuses: this owner currently provides fully coupled explicit SSPRK3.
    pub fn new(
        mut flow: PeriodicSlenderJet,
        caloric: &'m IsochoricHeatCapacityLaw,
        conduction: Option<&'m IsochoricConductionLaw<'m>>,
        initial_temperatures_k: &[f64],
        maximum_invariant_rate_per_s: f64,
        gate: &CancelGate,
    ) -> Result<Self, ThermalJetError> {
        poll(gate)?;
        let n = flow.accepted.a.len();
        let rho = flow.parameters.density_kg_m3;
        if initial_temperatures_k.len() != n
            || flow.carreau.is_none()
            || flow.integration_method != SlenderJetIntegrationMethod::ExplicitSspRk3
            || (rho - caloric.density_kg_m3()).abs()
                > 64.0 * f64::EPSILON * rho.max(caloric.density_kg_m3())
            || conduction.is_some_and(|law| law.caloric().identity() != caloric.identity())
        {
            return Err(ThermalJetError::InvalidInput(
                "complete temperatures, explicit Carreau flow and identical caloric density/identity required",
            ));
        }
        let reference_temperature_k = initial_temperatures_k[0];
        let mut accepted = zeros(n, gate)?;
        let mass_factor = positive(core::f64::consts::PI * rho * flow.spacing_m)?;
        for (i, &temperature) in initial_temperatures_k.iter().enumerate() {
            poll_at(gate, i)?;
            let _original_u = caloric
                .specific_internal_energy_j_kg(temperature)
                .map_err(|source| ThermalJetError::Thermal { cell: i, source })?;
            let specific = caloric
                .specific_internal_energy_change_j_kg(reference_temperature_k, temperature)
                .map_err(|source| ThermalJetError::Thermal { cell: i, source })?;
            accepted[i] = finite(positive(mass_factor * flow.accepted.a[i])? * specific)?;
            if specific != 0.0 && accepted[i] == 0.0 {
                return Err(NonlinearJetError::NonRepresentable.into());
            }
            if let Some(law) = conduction {
                law.conductivity_at_temperature_w_m_k(temperature)
                    .map_err(|source| ThermalJetError::Thermal { cell: i, source })?;
            }
        }
        flow.set_carreau_cell_temperatures(
            initial_temperatures_k,
            maximum_invariant_rate_per_s,
            gate,
        )?;
        flow.enable_cell_viscous_heat(gate)?;
        let profile = flow
            .carreau
            .as_ref()
            .ok_or(NonlinearJetError::InvalidState(
                "thermal jet material absent",
            ))?
            .try_clone_prepared(gate)?;
        let workspace = Workspace {
            profile,
            specific: zeros(n, gate)?,
            conductivity: if conduction.is_some() {
                zeros(n, gate)?
            } else {
                Vec::new()
            },
            conductance: if conduction.is_some() {
                zeros(n, gate)?
            } else {
                Vec::new()
            },
        };
        let source = Source {
            caloric,
            conduction,
            reference_temperature_k,
            maximum_temperature_slope: caloric.maximum_temperature_slope(),
        };
        let candidate = zeros(n, gate)?;
        let stage = zeros(n, gate)?;
        let rhs = zeros(n, gate)?;
        poll(gate)?;
        Ok(Self {
            flow,
            source,
            accepted,
            candidate,
            stage,
            rhs,
            workspace,
        })
    }

    /// Accepted native flow and source field. No uncoupled mutable view exists.
    pub const fn flow(&self) -> &PeriodicSlenderJet {
        &self.flow
    }
    pub const fn caloric(&self) -> &'m IsochoricHeatCapacityLaw {
        self.source.caloric
    }
    pub const fn conduction(&self) -> Option<&'m IsochoricConductionLaw<'m>> {
        self.source.conduction
    }
    /// Fixed numerical origin, independent of the source's declared energy zero.
    pub const fn reference_temperature_k(&self) -> f64 {
        self.source.reference_temperature_k
    }
    /// Accepted m*[u(T)-u(T_reference)] per cell [J], including transported heat.
    pub fn cell_internal_energy_from_reference_j(&self) -> &[f64] {
        &self.accepted
    }
    pub fn temperatures_k(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.flow
            .prepared_carreau_cells()
            .expect("thermal jet retains prepared cells")
            .iter()
            .map(|cell| cell.temperature_k())
    }
    pub fn cell_masses_kg(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        let factor =
            core::f64::consts::PI * self.flow.parameters.density_kg_m3 * self.flow.spacing_m;
        self.flow.accepted.a.iter().map(move |&a| factor * a)
    }

    fn evaluate(&mut self, stage: bool, gate: &CancelGate) -> Result<Evaluation, ThermalJetError> {
        let fields = if stage {
            &self.flow.stage
        } else {
            &self.flow.candidate
        };
        let energy = if stage { &self.stage } else { &self.candidate };
        evaluate(
            fields,
            energy,
            &self.source,
            &mut self.workspace,
            self.flow.parameters,
            self.flow.limits,
            self.flow.spacing_m,
            &mut self.flow.spatial,
            &mut self.flow.rhs,
            &mut self.rhs,
            gate,
        )
    }

    /// Advance all conserved fields together. Every stage uses actual T and rate.
    ///
    /// A source/rate, budget, arithmetic, conservation or cancellation refusal
    /// rolls the entire call back. A native geometry stop publishes only complete
    /// admitted coupled substeps. The caller's maximum step controls accuracy;
    /// no global error certificate or unconditional thermal stability is claimed.
    pub fn advance(
        &mut self,
        duration_s: f64,
        budget: ThermalJetStepBudget,
        gate: &CancelGate,
    ) -> Result<ThermalJetAdvance, ThermalJetError> {
        poll(gate)?;
        let b = budget.flow;
        if !duration_s.is_finite()
            || duration_s < 0.0
            || !b.maximum_step_s.is_finite()
            || b.maximum_step_s <= 0.0
            || !b.cfl.is_finite()
            || b.cfl <= 0.0
            || b.cfl > 0.25
            || b.maximum_substeps == 0
            || !budget.energy_tolerance_j.is_finite()
            || budget.energy_tolerance_j < 0.0
        {
            return Err(ThermalJetError::InvalidInput(
                "duration or coupled step budget",
            ));
        }
        let target = self.flow.time_s + duration_s;
        if !target.is_finite() || (duration_s > 0.0 && target <= self.flow.time_s) {
            return Err(NonlinearJetError::NonRepresentable.into());
        }
        let n = self.accepted.len();
        for i in 0..n {
            poll_at(gate, i)?;
            self.flow.candidate.a[i] = self.flow.accepted.a[i];
            self.flow.candidate.q[i] = self.flow.accepted.q[i];
            self.candidate[i] = self.accepted[i];
        }
        self.flow.spatial.reset_viscous_heat(gate)?;
        let before = self.flow.diagnostics;
        let mut base = self.evaluate(false, gate)?;
        let (mut elapsed, mut step_cap) = (0.0, b.maximum_step_s);
        let (mut accepted_substeps, mut attempted_substeps) = (0, 0);
        let (mut viscous_loss, mut upwind_loss) = (0.0, 0.0);
        let mut outcome = SlenderJetOutcome::Completed;
        while elapsed < duration_s {
            poll(gate)?;
            if attempted_substeps >= b.maximum_substeps {
                return Err(NonlinearJetError::SubstepBudget.into());
            }
            let remaining = duration_s - elapsed;
            let dt = positive(
                remaining
                    .min(step_cap)
                    .min(stable_step(base.rate_per_s, b.cfl)),
            )?;
            let next_elapsed = if dt == remaining {
                duration_s
            } else {
                elapsed + dt
            };
            if next_elapsed <= elapsed
                || self.flow.time_s + next_elapsed <= self.flow.time_s + elapsed
            {
                return Err(NonlinearJetError::NonRepresentable.into());
            }
            attempted_substeps += 1;
            for i in 0..n {
                poll_at(gate, i)?;
                self.flow.stage.a[i] = self.flow.candidate.a[i] + dt * self.flow.rhs.a[i];
                self.flow.stage.q[i] = self.flow.candidate.q[i] + dt * self.flow.rhs.q[i];
                self.stage[i] = finite(self.candidate[i] + dt * self.rhs[i])?;
            }
            let (end, physical, numerical) = match self.finish_stages(dt, b.cfl, base, gate) {
                Err(ThermalJetError::Flow(NonlinearJetError::InitialDomain(violation))) => {
                    outcome = SlenderJetOutcome::DomainLimit(violation);
                    break;
                }
                Err(error) => return Err(error),
                Ok(None) => {
                    step_cap = positive(0.5 * dt)?;
                    base = self.evaluate(false, gate)?;
                    continue;
                }
                Ok(Some(result)) => result,
            };
            self.flow.spatial.accept_viscous_heat_step(gate)?;
            viscous_loss = finite(viscous_loss + physical)?;
            upwind_loss = finite(upwind_loss + numerical)?;
            core::mem::swap(&mut self.flow.candidate, &mut self.flow.stage);
            core::mem::swap(&mut self.candidate, &mut self.stage);
            base = end;
            elapsed = next_elapsed;
            accepted_substeps += 1;
            step_cap = b.maximum_step_s;
        }
        // A failed geometry stage may have changed only scratch material.
        // Rebuild precisely the last admitted temperature/source before publish.
        if matches!(outcome, SlenderJetOutcome::DomainLimit(_)) {
            base = self.evaluate(false, gate)?;
        }
        let mut storage = Sum::default();
        let mut mass = Sum::default();
        let mass_factor = positive(
            core::f64::consts::PI * self.flow.parameters.density_kg_m3 * self.flow.spacing_m,
        )?;
        for i in 0..n {
            poll_at(gate, i)?;
            storage.add(self.candidate[i] - self.accepted[i]);
            mass.add(mass_factor * (self.flow.candidate.a[i] - self.flow.accepted.a[i]));
        }
        let internal_energy_change_j = finite(storage.total())?;
        let thermal_residual = finite(internal_energy_change_j - viscous_loss)?;
        if thermal_residual.abs() > budget.energy_tolerance_j {
            return Err(ThermalJetError::EnergyBalance {
                residual_j: thermal_residual,
                tolerance_j: budget.energy_tolerance_j,
            });
        }
        let after = base.diagnostics;
        let mechanical_residual = finite(
            (after.kinetic_energy_j - before.kinetic_energy_j)
                + (after.surface_energy_j - before.surface_energy_j)
                + viscous_loss
                + upwind_loss,
        )?;
        let end_time = if elapsed == duration_s {
            target
        } else {
            self.flow.time_s + elapsed
        };
        let report = ThermalJetAdvance {
            flow: SlenderJetAdvance {
                start_time_s: self.flow.time_s,
                end_time_s: end_time,
                requested_duration_s: duration_s,
                accepted_substeps,
                attempted_substeps,
                outcome,
                before,
                after,
                integration_method: SlenderJetIntegrationMethod::ExplicitSspRk3,
                temporal_error: None,
                viscous_dissipation_j: viscous_loss,
                upwind_dissipation_j: upwind_loss,
                energy_balance_residual_j: mechanical_residual,
            },
            internal_energy_change_j,
            energy_balance_residual_j: thermal_residual,
            mass_balance_residual_kg: finite(mass.total())?,
            combined_energy_residual_j: finite(mechanical_residual + thermal_residual)?,
        };
        let partition = self
            .flow
            .spatial
            .viscous_heat_partition(viscous_loss, gate)?;
        let reference = positive(
            self.workspace.profile.coefficients(0.0, 0)?[0] / self.flow.parameters.density_kg_m3,
        )?;
        poll(gate)?;
        self.flow
            .spatial
            .publish_viscous_heat(viscous_loss, partition);
        core::mem::swap(&mut self.flow.accepted, &mut self.flow.candidate);
        core::mem::swap(&mut self.accepted, &mut self.candidate);
        // The constructor admitted this variant; the optional owner cannot change.
        if let Some(profile) = &mut self.flow.carreau {
            core::mem::swap(profile, &mut self.workspace.profile);
        }
        self.flow.parameters.kinematic_viscosity_m2_s = reference;
        self.flow.time_s = end_time;
        self.flow.diagnostics = after;
        Ok(report)
    }

    fn finish_stages(
        &mut self,
        dt: f64,
        cfl: f64,
        base: Evaluation,
        gate: &CancelGate,
    ) -> Result<Option<(Evaluation, f64, f64)>, ThermalJetError> {
        self.flow.spatial.start_viscous_heat(dt / 6.0, gate)?;
        let first = self.evaluate(true, gate)?;
        if dt > stable_step(first.rate_per_s, cfl) {
            return Ok(None);
        }
        self.flow.spatial.add_viscous_heat(dt / 6.0, gate)?;
        self.blend(dt, 0.25, gate)?;
        let second = self.evaluate(true, gate)?;
        if dt > stable_step(second.rate_per_s, cfl) {
            return Ok(None);
        }
        self.flow.spatial.add_viscous_heat((2.0 / 3.0) * dt, gate)?;
        self.blend(dt, 2.0 / 3.0, gate)?;
        let end = self.evaluate(true, gate)?;
        if dt > stable_step(end.rate_per_s, cfl) {
            return Ok(None);
        }
        let physical = finite(
            dt * (base.viscous_power_w / 6.0
                + first.viscous_power_w / 6.0
                + (2.0 / 3.0) * second.viscous_power_w),
        )?;
        let numerical = finite(
            dt * (base.upwind_power_w / 6.0
                + first.upwind_power_w / 6.0
                + (2.0 / 3.0) * second.upwind_power_w),
        )?;
        Ok(Some((end, physical, numerical)))
    }

    fn blend(&mut self, dt: f64, factor: f64, gate: &CancelGate) -> Result<(), ThermalJetError> {
        blend(
            &mut self.flow.stage,
            &self.flow.candidate,
            &self.flow.rhs,
            dt,
            factor,
            gate,
        )?;
        for i in 0..self.stage.len() {
            poll_at(gate, i)?;
            self.stage[i] = finite(
                self.candidate[i]
                    + factor * ((self.stage[i] - self.candidate[i]) + dt * self.rhs[i]),
            )?;
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn evaluate(
    fields: &Fields,
    energy: &[f64],
    source: &Source<'_>,
    workspace: &mut Workspace,
    parameters: SlenderJetParameters,
    limits: SlenderJetLimits,
    dx: f64,
    spatial: &mut Spatial,
    rhs: &mut Fields,
    energy_rhs: &mut [f64],
    gate: &CancelGate,
) -> Result<Evaluation, ThermalJetError> {
    let n = fields.a.len();
    let rho_pi = positive(parameters.density_kg_m3 * core::f64::consts::PI)?;
    let mass_factor = positive(rho_pi * dx)?;
    for i in 0..n {
        poll_at(gate, i)?;
        let mass = positive(mass_factor * fields.a[i])?;
        let specific = finite(energy[i] / mass)?;
        if energy[i] != 0.0 && specific == 0.0 {
            return Err(NonlinearJetError::NonRepresentable.into());
        }
        workspace.specific[i] = specific;
        let temperature = source
            .caloric
            .temperature_at_relative_internal_energy(source.reference_temperature_k, specific)
            .map_err(|source| ThermalJetError::Thermal { cell: i, source })?;
        workspace.profile.reprepare_temperature(i, temperature)?;
        if let Some(law) = source.conduction {
            workspace.conductivity[i] = positive(
                law.conductivity_at_temperature_w_m_k(temperature)
                    .map_err(|source| ThermalJetError::Thermal { cell: i, source })?,
            )?;
        }
    }
    let mut evaluation = evaluate_constitutive(
        fields,
        parameters,
        Some(&workspace.profile),
        limits,
        dx,
        spatial,
        rhs,
        true,
        gate,
    )?;
    // The native momentum flux has already supplied rhs.q. Reuse its scratch
    // for the distinct energy flux, without changing the area/momentum kernels.
    for i in 0..n {
        poll_at(gate, i)?;
        let r = next(i, n);
        let donor = if fields.q[i] >= 0.0 { i } else { r };
        let mut flux = finite(rho_pi * fields.q[i] * workspace.specific[donor])?;
        if source.conduction.is_some() {
            let left = positive(workspace.conductivity[i] * fields.a[i])?;
            let right = positive(workspace.conductivity[r] * fields.a[r])?;
            let low = left.min(right);
            let high = left.max(right);
            let harmonic = positive(low / (0.5 + 0.5 * (low / high)))?;
            let conductance = positive(core::f64::consts::PI * harmonic / dx)?;
            workspace.conductance[i] = conductance;
            flux = finite(
                flux + conductance
                    * (workspace.profile.temperature_k(i) - workspace.profile.temperature_k(r)),
            )?;
        }
        spatial.flux[i] = flux;
    }
    let powers = spatial.cell_viscous_power_w()?;
    let mut thermal_rate = 0.0_f64;
    for i in 0..n {
        poll_at(gate, i)?;
        let l = previous(i, n);
        energy_rhs[i] = finite(-(spatial.flux[i] - spatial.flux[l]) + powers[i])?;
        // Gross outgoing mass, not net drainage. Balanced flow through a thin
        // cell may have zero drainage and still require a small donor step.
        let advective =
            finite((fields.q[i].max(0.0) + (-fields.q[l]).max(0.0)) / fields.a[i] / dx)?;
        let conductive = if source.conduction.is_some() {
            finite(
                source.maximum_temperature_slope
                    * (workspace.conductance[l] + workspace.conductance[i])
                    / (mass_factor * fields.a[i]),
            )?
        } else {
            0.0
        };
        thermal_rate = thermal_rate.max(finite(advective + conductive)?);
    }
    evaluation.rate_per_s = finite(evaluation.rate_per_s + thermal_rate)?;
    Ok(evaluation)
}

#[cfg(test)]
mod tests;
