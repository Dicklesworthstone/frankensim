//! Original real-gas calorics in the shared conservative wall source.

use super::{RealGas, RealGasEuler1d, RealGasFlowError};
use crate::dynamics::{CaloricDomain, Conserved, Mesh1d, Primitive, allocated, poll, positive};
use crate::multicomponent::wall_exchange::{
    MaterialWallCaloric, WallCellReport, WallContact, WallControl, WallExchangeError,
    WallGasCaloric, WallGasState, WallThermal, exchange_caloric,
    exchange_caloric_with_heat_storage,
};
use fs_qty::Time;

mod shared;
pub use shared::{
    SharedWall, SharedWallControl, SharedWallGas, SharedWallLink, SharedWallLinkReport,
    SharedWallPort, SharedWallReport,
};

/// One cell's contact and optional original finite-wall caloric law.
///
/// Without a source, `Lumped` retains its explicitly authored constant capacity.
/// With a source, that capacity is only a tangent cache; its original integral
/// or h/u inventory owns storage, including latent heat. No Cp-to-Cv conversion
/// or unrepresented wall pressure/deformation work is inferred.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasWallContact<'law> {
    pub contact: WallContact,
    pub heat_storage: Option<MaterialWallCaloric<'law>>,
}

/// Accepted local sources and the complete transport-plus-wall storage balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasWallReport {
    pub elapsed: Time,
    /// Opposite axial impulse received by the walls/drives [N s].
    pub wall_impulse_ns: f64,
    pub heat_to_gas_j: f64,
    /// Finite-wall change in each original source's declared h/u coordinate [J].
    pub wall_storage_change_j: f64,
    pub external_heat_to_gas_j: f64,
    /// Prescribed sliding-wall work entering gas, separate from heat [J].
    pub mechanical_work_to_gas_j: f64,
    /// Relative kinetic energy already retained as gas internal energy [J].
    pub drag_heating_j: f64,
    pub source_storage_residual_j: f64,
    /// Gas plus declared wall storage minus boundary transport, reservoir heat
    /// sliding-drive work and external body-force work [J]. An enthalpy wall retains its work convention;
    /// this mixed-coordinate account is not a claim that its H is rigid-body U.
    pub gas_wall_storage_residual_j: f64,
    pub maximum_thermal_residual_j: f64,
    pub maximum_iterations: usize,
}

/// Prepared, per-cell real-gas wall coupling with preallocated transactional
/// contact storage. The original model, energy gauge and physical mesh are bound
/// once. Gas states, finite walls, reports and flow clock publish together.
/// Each contact owns an independent finite amount; a wall shared across cells
/// is bound explicitly with [`Self::with_shared_walls`]. Reusing a constitutive
/// source alone does not share its finite heat inventory.
///
/// Local heat is backward Euler. The stationary entry applies it after waves;
/// the piston entry surrounds pressure work with two source halves. Both
/// transport/source compositions are first order. Conductance and radial drag are
/// explicitly authored engineering inputs, not inferred dense-gas correlations.
pub struct RealGasWallOperator<'law> {
    gas: RealGas,
    mesh: Mesh1d,
    contacts: Vec<WallContact>,
    candidate_contacts: Vec<WallContact>,
    heat_storage: Vec<Option<MaterialWallCaloric<'law>>>,
    shared: Option<shared::SharedWalls<'law>>,
    last_report: Option<RealGasWallReport>,
    split: Option<SplitWorkspace>,
}

// One reusable intermediate wall inventory and unpublished moved geometry.
// The gas owner's existing source buffers retain the first thermal endpoint.
struct SplitWorkspace {
    contacts: Vec<WallContact>,
    mesh: Mesh1d,
}

impl<'law> RealGasWallOperator<'law> {
    pub fn new(
        flow: &RealGasEuler1d,
        inputs: Vec<RealGasWallContact<'law>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RealGasFlowError> {
        poll(&mut cancelled)?;
        let n = flow.states().len();
        if inputs.len() != n {
            return Err(RealGasFlowError::Invalid(
                "one wall contact per real-gas cell required",
            ));
        }
        let empty = WallContact {
            drag: crate::multicomponent::wall_exchange::WallDrag::None,
            thermal: WallThermal::Adiabatic,
        };
        let mut contacts = allocated(n, empty)?;
        let candidate_contacts = allocated(n, empty)?;
        let mut heat_storage = allocated(n, None)?;
        for (i, input) in inputs.into_iter().enumerate() {
            poll(&mut cancelled)?;
            input.contact.validate()?;
            contacts[i] = if let Some(source) = &input.heat_storage {
                let contact = source.bind_contact(input.contact)?;
                source.temperature_k(contact.thermal)?;
                contact
            } else {
                if matches!(input.contact.thermal, WallThermal::Caloric { .. }) {
                    return Err(RealGasFlowError::Invalid(
                        "caloric wall requires its original source",
                    ));
                }
                input.contact
            };
            heat_storage[i] = input.heat_storage;
        }
        poll(&mut cancelled)?;
        Ok(Self {
            gas: flow.gas().clone(),
            mesh: flow.mesh().clone(),
            contacts,
            candidate_contacts,
            heat_storage,
            shared: None,
            last_report: None,
            split: None,
        })
    }

    pub fn contacts(&self) -> &[WallContact] {
        &self.contacts
    }
    pub fn heat_storage(&self) -> &[Option<MaterialWallCaloric<'law>>] {
        &self.heat_storage
    }
    /// Failed or cancelled complete steps retain the preceding accepted report.
    pub const fn last_report(&self) -> Option<RealGasWallReport> {
        self.last_report
    }

    pub(super) fn validate_binding(&self, flow: &RealGasEuler1d) -> Result<(), RealGasFlowError> {
        if self.gas != *flow.gas() || self.mesh != *flow.mesh() {
            return Err(RealGasFlowError::Invalid(
                "real-gas wall source, gauge or mesh binding differs",
            ));
        }
        Ok(())
    }

    pub(super) fn prepare(
        &mut self,
        states: &mut [Conserved],
        primitives: &mut [Primitive],
        duration: Time,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RealGasWallReport, RealGasFlowError> {
        let mut total = HeatSum::from_report(prepare_contacts(
            &self.gas,
            self.mesh.volumes(),
            &self.contacts,
            &self.heat_storage,
            &mut self.candidate_contacts,
            states,
            primitives,
            duration,
            controls,
            cancelled,
        )?);
        if let Some(shared) = &mut self.shared {
            total.add(shared.prepare(
                &self.gas, &self.mesh, states, primitives, duration, controls, cancelled,
            )?)?;
        }
        Ok(total.report(duration))
    }

    pub(super) fn prepare_piston_first(
        &mut self,
        states: &mut [Conserved],
        primitives: &mut [Primitive],
        duration: Time,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RealGasWallReport, RealGasFlowError> {
        if self.split.is_none() {
            self.split = Some(SplitWorkspace {
                contacts: allocated(self.contacts.len(), self.contacts[0])?,
                mesh: self.mesh.try_copy()?,
            });
        }
        let split = self
            .split
            .as_mut()
            .expect("prepared real-gas thermal split");
        let mut total = HeatSum::from_report(prepare_contacts(
            &self.gas,
            self.mesh.volumes(),
            &self.contacts,
            &self.heat_storage,
            &mut split.contacts,
            states,
            primitives,
            duration,
            controls,
            cancelled,
        )?);
        if let Some(shared) = &mut self.shared {
            total.add(shared.prepare_piston_first(
                &self.gas, &self.mesh, states, primitives, duration, controls, cancelled,
            )?)?;
        }
        Ok(total.report(duration))
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare_piston_second(
        &mut self,
        states: &mut [Conserved],
        primitives: &mut [Primitive],
        mesh: &Mesh1d,
        duration: Time,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RealGasWallReport, RealGasFlowError> {
        let split = self
            .split
            .as_mut()
            .expect("prepared real-gas thermal split");
        // The accepted binding stays unchanged until the entire piston step
        // succeeds. This copy preserves the original GCL volumes without a
        // second geometry integration or any per-step allocation.
        split.mesh.reset_motion_from(mesh, cancelled)?;
        let mut total = HeatSum::from_report(prepare_contacts(
            &self.gas,
            split.mesh.volumes(),
            &split.contacts,
            &self.heat_storage,
            &mut self.candidate_contacts,
            states,
            primitives,
            duration,
            controls,
            cancelled,
        )?);
        if let Some(shared) = &mut self.shared {
            total.add(shared.prepare_piston_second(
                &self.gas,
                &split.mesh,
                states,
                primitives,
                duration,
                controls,
                cancelled,
            )?)?;
        }
        Ok(total.report(duration))
    }

    pub(super) fn commit(&mut self, report: RealGasWallReport) {
        core::mem::swap(&mut self.contacts, &mut self.candidate_contacts);
        if let Some(shared) = &mut self.shared {
            shared.commit();
        }
        self.last_report = Some(report);
    }

    pub(super) fn commit_piston(&mut self, report: RealGasWallReport) {
        let split = self
            .split
            .as_mut()
            .expect("prepared real-gas thermal split");
        core::mem::swap(&mut self.mesh, &mut split.mesh);
        self.commit(report);
    }
}

#[allow(clippy::too_many_arguments)]
fn prepare_contacts(
    gas: &RealGas,
    volumes: &[f64],
    contacts: &[WallContact],
    heat_storage: &[Option<MaterialWallCaloric<'_>>],
    candidate_contacts: &mut [WallContact],
    states: &mut [Conserved],
    primitives: &mut [Primitive],
    duration: Time,
    controls: WallControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RealGasWallReport, RealGasFlowError> {
    let mut total = HeatSum::default();
    for (i, state) in states.iter_mut().enumerate() {
        poll(cancelled)?;
        let source = heat_storage[i].as_ref();
        let contact = contacts[i];
        let hint = match source {
            Some(source) => Some(source.temperature_k(contact.thermal)?),
            None => match contact.thermal {
                WallThermal::Reservoir { temperature_k, .. }
                | WallThermal::Lumped { temperature_k, .. } => Some(temperature_k),
                _ => None,
            },
        };
        let caloric = RealGasCaloric {
            gas,
            density: state.density,
            hint,
        };
        let update = match source {
            Some(source) => exchange_caloric_with_heat_storage(
                &caloric,
                state.momentum,
                state.energy,
                volumes[i],
                duration,
                contact,
                source,
                None,
                controls,
                &mut *cancelled,
            )?,
            None => exchange_caloric(
                &caloric,
                state.momentum,
                state.energy,
                volumes[i],
                duration,
                contact,
                controls,
                &mut *cancelled,
            )?,
        };
        let candidate = Conserved {
            density: state.density,
            momentum: update.momentum,
            energy: update.energy,
        };
        // Final pressure, transport and acoustic admission always returns
        // to the same original density/energy owner, including unchanged cells.
        primitives[i] = gas.decode(candidate, cancelled)?;
        total.add(update.report)?;
        candidate_contacts[i] = update.contact;
        *state = candidate;
    }
    poll(cancelled)?;
    Ok(total.report(duration))
}

/// Borrow the SAME density-dependent original calorics and common energy gauge.
/// The shared heat kernel needs rho*U(T) and rho*Cv(T) at this frozen density;
/// neither ideal-gas Cp-R nor a constant tangent substitutes for real-gas Cv.
struct RealGasCaloric<'a> {
    gas: &'a RealGas,
    density: f64,
    hint: Option<f64>,
}

impl RealGasCaloric<'_> {
    fn specific_caloric(&self, temperature: f64) -> Result<(f64, f64), WallExchangeError> {
        if !self.gas.domain.admits_density(self.density)
            || !self.gas.domain.admits_temperature(temperature)
        {
            return Err(
                RealGasFlowError::Invalid("wall trial outside original real-gas domain").into(),
            );
        }
        let reference = self.gas.reference_state();
        // One finite increment checks both original EOS/caloric endpoints and
        // supplies complete Cv. Transport and the PT root inventory belong to
        // the unconditional accepted-state decode in local/shared prepare().
        let caloric = self
            .gas
            .model()
            .specific_internal_energy_change_and_cv_between_densities(
                reference.density_kg_m3(),
                reference.temperature_k(),
                self.density,
                temperature,
            )
            .map_err(RealGasFlowError::from)?;
        Ok((caloric.internal_energy_change_j_kg(), caloric.final_cv_j_kg_k()))
    }
}

impl WallGasCaloric for RealGasCaloric<'_> {
    fn domain(&self) -> CaloricDomain {
        let energy = self.gas.domain().energy;
        CaloricDomain {
            minimum_temperature: energy.minimum_temperature_k,
            maximum_temperature: energy.maximum_temperature_k,
            absolute_energy_tolerance: energy.energy_tolerance_j_kg,
            relative_energy_tolerance: 0.0,
            maximum_iterations: energy.maximum_iterations,
        }
    }

    fn primitive(
        &self,
        momentum: f64,
        energy: f64,
        mut cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<WallGasState, WallExchangeError> {
        poll(&mut cancelled)?;
        let state = Conserved { density: self.density, momentum, energy };
        if !state.finite() || !self.gas.domain.admits_density(self.density) {
            return Err(RealGasFlowError::Invalid(
                "conserved density outside real-gas domain or nonfinite state",
            ).into());
        }
        // Match RealGas::resolve's specific total/kinetic subtraction. The
        // source-relative caloric query never subtracts formation energies.
        let velocity = momentum / self.density;
        let specific_total = energy / self.density;
        let kinetic = (0.5 * velocity) * velocity;
        let relative_energy = specific_total - kinetic;
        if !velocity.is_finite() || !relative_energy.is_finite() {
            return Err(RealGasFlowError::Invalid(
                "unrepresentable real-gas specific energy",
            ).into());
        }
        if let Some(temperature) = self.hint.filter(|&t| self.gas.domain.admits_temperature(t)) {
            let (thermal, _) = self.caloric(temperature)?;
            let kinetic_density = 0.5 * momentum * velocity;
            if energy == thermal + kinetic_density {
                poll(&mut cancelled)?;
                return Ok(WallGasState {
                    density: self.density,
                    velocity,
                    temperature,
                });
            }
        }
        let reference = self.gas.reference_state();
        let energy_controls = self.gas.domain.energy;
        // Preserve the same bounded reference-density subtraction-roundoff
        // shortcut as RealGas::resolve. Conserved energy is never altered.
        let roundoff = (8.0 * f64::EPSILON * specific_total.abs()
            + 8.0 * f64::EPSILON * kinetic.abs())
            .min(energy_controls.energy_tolerance_j_kg);
        if self.density.to_bits() == reference.density_kg_m3().to_bits()
            && relative_energy.abs() <= roundoff
        {
            poll(&mut cancelled)?;
            return Ok(WallGasState {
                density: self.density, velocity, temperature: reference.temperature_k(),
            });
        }
        // The original increment inversion does not absorb nonzero relative
        // energy merely because it is smaller than the absolute allowance.
        let tolerance = if relative_energy == 0.0 {
            energy_controls.energy_tolerance_j_kg
        } else {
            energy_controls.energy_tolerance_j_kg.min(0.5 * relative_energy.abs())
        };
        // Solve a one-kilogram caloric balance: the numerical Delta u and Cv
        // values are joules and J/K for that unit mass. A zero equation scale
        // keeps the original absolute J/kg tolerance; WallControl's required
        // positive relative allowance contributes exactly zero. No rho*u/rho
        // round trip or invented physical conductance enters this inversion.
        let control = WallControl {
            absolute_energy_tolerance_j: tolerance,
            maximum_iterations: energy_controls.maximum_iterations,
            ..WallControl::default()
        };
        control.validate()?;
        let (minimum, maximum) = if self.density != reference.density_kg_m3() {
            (energy_controls.minimum_temperature_k, energy_controls.maximum_temperature_k)
        } else if relative_energy > 0.0 {
            (reference.temperature_k(), energy_controls.maximum_temperature_k)
        } else {
            (energy_controls.minimum_temperature_k, reference.temperature_k())
        };
        let initial = self.hint.filter(|&t| self.gas.domain.admits_temperature(t))
            .unwrap_or(reference.temperature_k());
        let sample = |temperature| -> Result<_, WallExchangeError> {
            let (internal, cv) = self.specific_caloric(temperature)?;
            let residual = internal - relative_energy;
            if !residual.is_finite() || !positive(cv) {
                return Err(RealGasFlowError::Invalid(
                    "unrepresentable original real-gas caloric inversion",
                ).into());
            }
            Ok((residual, Some(cv), 0.0))
        };
        // Preserve the original caloric inversion's endpoint acceptance. The
        // shared scalar owner otherwise takes strictly interior Newton steps,
        // which cannot reach a boundary at sub-ULP energy tolerances. Cache
        // these samples so an interior solve does not evaluate either twice.
        poll(&mut cancelled)?;
        let lower = sample(minimum)?;
        poll(&mut cancelled)?;
        let upper = sample(maximum)?;
        let temperature = if lower.0.abs() <= tolerance {
            minimum
        } else if upper.0.abs() <= tolerance {
            maximum
        } else {
            crate::multicomponent::wall_exchange::solve_scalar(
                minimum,
                maximum,
                initial,
                control,
                |temperature, mut checkpoint| {
                    poll(&mut checkpoint)?;
                    if temperature == minimum {
                        Ok(lower)
                    } else if temperature == maximum {
                        Ok(upper)
                    } else {
                        sample(temperature)
                    }
                },
                &mut cancelled,
            )?.0
        };
        poll(&mut cancelled)?;
        Ok(WallGasState {
            density: self.density,
            velocity,
            temperature,
        })
    }

    fn caloric(&self, temperature: f64) -> Result<(f64, f64), WallExchangeError> {
        let (internal, cv) = self.specific_caloric(temperature)?;
        let energy = self.density * internal;
        let capacity = self.density * cv;
        if !energy.is_finite() || !positive(capacity) {
            return Err(
                RealGasFlowError::Invalid("unrepresentable original real-gas calorics").into(),
            );
        }
        Ok((energy, capacity))
    }
}

#[derive(Default)]
struct HeatSum {
    values: [f64; 7],
    correction: [f64; 7],
    maximum_thermal_residual_j: f64,
    maximum_iterations: usize,
}

impl HeatSum {
    fn from_report(report: RealGasWallReport) -> Self {
        Self {
            values: [
                report.wall_impulse_ns,
                report.heat_to_gas_j,
                report.wall_storage_change_j,
                report.external_heat_to_gas_j,
                report.mechanical_work_to_gas_j,
                report.drag_heating_j,
                report.source_storage_residual_j,
            ],
            correction: [0.0; 7],
            maximum_thermal_residual_j: report.maximum_thermal_residual_j,
            maximum_iterations: report.maximum_iterations,
        }
    }

    fn add(&mut self, report: WallCellReport) -> Result<(), RealGasFlowError> {
        for (i, value) in [
            report.wall_impulse_ns,
            report.heat_to_gas_j,
            report.wall_energy_change_j,
            report.external_heat_to_gas_j,
            report.mechanical_work_to_gas_j,
            report.drag_heating_j,
            report.energy_residual_j,
        ]
        .into_iter()
        .enumerate()
        {
            let increment = value - self.correction[i];
            let next = self.values[i] + increment;
            if !next.is_finite() {
                return Err(RealGasFlowError::Invalid(
                    "unrepresentable combined real-gas wall source",
                ));
            }
            self.correction[i] = (next - self.values[i]) - increment;
            self.values[i] = next;
        }
        self.maximum_thermal_residual_j = self
            .maximum_thermal_residual_j
            .max(report.thermal_residual_j.abs());
        self.maximum_iterations = self.maximum_iterations.max(report.iterations);
        Ok(())
    }

    fn report(self, elapsed: Time) -> RealGasWallReport {
        RealGasWallReport {
            elapsed,
            wall_impulse_ns: self.values[0],
            heat_to_gas_j: self.values[1],
            wall_storage_change_j: self.values[2],
            external_heat_to_gas_j: self.values[3],
            mechanical_work_to_gas_j: self.values[4],
            drag_heating_j: self.values[5],
            source_storage_residual_j: self.values[6],
            gas_wall_storage_residual_j: 0.0,
            maximum_thermal_residual_j: self.maximum_thermal_residual_j,
            maximum_iterations: self.maximum_iterations,
        }
    }
}
