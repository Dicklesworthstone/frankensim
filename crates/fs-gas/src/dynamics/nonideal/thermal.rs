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
pub use shared::{SharedWall, SharedWallControl, SharedWallGas, SharedWallPort, SharedWallReport};

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
/// Local heat is backward Euler after the accepted wave interval: the overall
/// transport/source splitting is first order. Conductance and radial drag are
/// explicitly authored engineering inputs, not inferred dense-gas correlations.
pub struct RealGasWallOperator<'law> {
    gas: RealGas,
    mesh: Mesh1d,
    contacts: Vec<WallContact>,
    candidate_contacts: Vec<WallContact>,
    heat_storage: Vec<Option<MaterialWallCaloric<'law>>>,
    shared: Option<shared::SharedWalls<'law>>,
    last_report: Option<RealGasWallReport>,
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
        let mut total = HeatSum::default();
        for (i, state) in states.iter_mut().enumerate() {
            poll(cancelled)?;
            let source = self.heat_storage[i].as_ref();
            let contact = self.contacts[i];
            let hint = match source {
                Some(source) => Some(source.temperature_k(contact.thermal)?),
                None => match contact.thermal {
                    WallThermal::Reservoir { temperature_k, .. }
                    | WallThermal::Lumped { temperature_k, .. } => Some(temperature_k),
                    _ => None,
                },
            };
            let caloric = RealGasCaloric {
                gas: &self.gas,
                density: state.density,
                hint,
            };
            let update = match source {
                Some(source) => exchange_caloric_with_heat_storage(
                    &caloric,
                    state.momentum,
                    state.energy,
                    self.mesh.volumes()[i],
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
                    self.mesh.volumes()[i],
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
            primitives[i] = self.gas.decode(candidate, cancelled)?;
            total.add(update.report)?;
            self.candidate_contacts[i] = update.contact;
            *state = candidate;
        }
        if let Some(shared) = &mut self.shared {
            total.add(shared.prepare(&self.gas, states, primitives, duration, controls, cancelled)?)?;
        }
        poll(cancelled)?;
        Ok(total.report(duration))
    }

    pub(super) fn commit(&mut self, report: RealGasWallReport) {
        core::mem::swap(&mut self.contacts, &mut self.candidate_contacts);
        if let Some(shared) = &mut self.shared {
            shared.commit();
        }
        self.last_report = Some(report);
    }
}

/// Borrow the SAME density-dependent original calorics and common energy gauge.
/// The shared heat kernel needs rho*U(T) and rho*Cv(T) at this frozen density;
/// neither ideal-gas Cp-R nor a constant tangent substitutes for real-gas Cv.
struct RealGasCaloric<'a> {
    gas: &'a RealGas,
    density: f64,
    hint: Option<f64>,
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
        if let Some(temperature) = self.hint.filter(|&t| self.gas.domain.admits_temperature(t)) {
            let (thermal, _) = self.caloric(temperature)?;
            let kinetic = 0.5 * momentum * (momentum / self.density);
            if energy == thermal + kinetic {
                return Ok(WallGasState {
                    density: self.density,
                    velocity: momentum / self.density,
                    temperature,
                });
            }
        }
        let state = self.gas.decode(
            Conserved {
                density: self.density,
                momentum,
                energy,
            },
            &mut cancelled,
        )?;
        Ok(WallGasState {
            density: state.density,
            velocity: state.velocity,
            temperature: state.temperature,
        })
    }

    fn caloric(&self, temperature: f64) -> Result<(f64, f64), WallExchangeError> {
        if !self.gas.domain.admits_density(self.density)
            || !self.gas.domain.admits_temperature(temperature)
        {
            return Err(
                RealGasFlowError::Invalid("wall trial outside original real-gas domain").into(),
            );
        }
        let reference = self.gas.reference_state();
        // One original finite energy increment already resolves the endpoint
        // Cv and checks both EOS/caloric states. Speculative heat trials need
        // neither a PT cubic inventory nor dilute-gas transport. prepare()
        // still decodes every accepted candidate before publishing any state,
        // admitting its complete source axes, transport and density branch.
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
        let energy = self.density * caloric.internal_energy_change_j_kg();
        let capacity = self.density * caloric.final_cv_j_kg_k();
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
