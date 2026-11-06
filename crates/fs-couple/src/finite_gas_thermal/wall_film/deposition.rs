//! A separate dry-capable liquid film on an original fixed-mass thermal cell.
//!
//! The native heat graph owns every wall caloric law, contact and source. One
//! exposed patch replaces its authored ambient exchange with the actual finite
//! gas. Film sensible/latent heat and reciprocal gas convection have explicit
//! opposite ports; neither changes the wall material mass or its phase chart.

use super::{
    WallFilmThermalBoundary, WallFilmThermalError,
    implicit::{WallFilmImplicitControls, WallFilmImplicitError, WorkMeter},
};
use fs_material::evaporation::gas::{
    GasEvaporationError, GasEvaporationSolve,
    implicit::{ImplicitGasControls, ImplicitGasError},
    wall_film::{
        GasWallFilmInventory, GasWallFilmStep, GasWallFilmSurface,
        endpoint::{GasWallFilmEndpoint, GasWallFilmEndpointTrial},
    },
};
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::parcel::{ClosedGasParcel, transfer::species::boundary::GasSpeciesBoundaryCandidate},
};
use fs_thermal::enthalpy::backward_euler::{
    BackwardEulerStep, BackwardEulerTrial, BoundaryExchangeStepCoefficients,
    BoundaryStepCoefficients, StepCoefficients, TransportStepCoefficients,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit};
use std::cell::RefCell;

mod solve;
pub use solve::advance_deposition_with_wall_film;
pub(super) use solve::advance_deposition_with_wall_film_metered;

/// Original thermal events around the same coupled heat/film/gas endpoint.
pub mod timeline;

/// One fixed native wall and an independent original liquid inventory.
#[derive(Clone, Copy)]
pub struct DepositionWallFilmBinding<'a> {
    pub cell_id: CellId,
    pub film: GasWallFilmInventory,
    /// Area is an occupied patch of the native cell's currently exposed area.
    pub boundary: WallFilmThermalBoundary<'a>,
    /// Reciprocal sensible convection on this same patch [W/(m² K)].
    pub heat_transfer_coefficient_w_m2_k: f64,
}

/// Reuses the native heat/gas numerical and physical allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DepositionWallFilmControls {
    pub implicit: ImplicitGasControls,
    /// Original liquid h(T)/T(h) round trip at the actual wall temperature.
    pub temperature_tolerance_k: f64,
}
impl DepositionWallFilmControls {
    pub fn new(phase: GasEvaporationSolve) -> Self {
        Self {
            implicit: ImplicitGasControls::new(phase),
            temperature_tolerance_k: 1e-7,
        }
    }
    pub(super) fn numerics(self) -> WallFilmImplicitControls {
        let native = self.implicit.thermal;
        WallFilmImplicitControls {
            newton: native.newton,
            difference_step_j: native.difference_step_j,
            mass_coordinate_scale_j_kg: self.implicit.mass_coordinate_scale_j_kg,
            maximum_iterations: native.maximum_iterations,
            maximum_evaluations: native.maximum_evaluations,
        }
    }
    fn validate(self) -> Result<(), ImplicitGasError> {
        let native = self.implicit.thermal;
        native.newton.validate()?;
        if !positive(self.temperature_tolerance_k)
            || !positive(native.difference_step_j)
            || !positive(self.implicit.mass_coordinate_scale_j_kg)
            || !nonnegative(native.cell_energy_tolerance_j)
            || !nonnegative(native.total_energy_tolerance_j)
            || !nonnegative(native.mass_tolerance_kg)
            || native.maximum_iterations == 0
            || native.maximum_evaluations == 0
        {
            return Err(ImplicitGasError::Invalid(
                "invalid spatial wall-film controls",
            ));
        }
        Ok(())
    }
}

/// Native physical accounts; vapor energy already includes its latent portion.
#[derive(Clone, Debug, PartialEq)]
pub struct DepositionWallFilmReport {
    /// Native heat audit INCLUDING internal film/gas debits. Use the separately
    /// evaluated external_thermal_heat_j for the enclosing system's supply.
    pub thermal: EnergyAudit,
    pub phase: GasWallFilmStep,
    pub external_thermal_heat_j: f64,
    pub external_gas_heat_j: f64,
    pub sensible_heat_into_gas_j: f64,
    pub wall_heat_into_phase_system_j: f64,
    pub coupled_energy_residual_j: f64,
    pub sum_absolute_energy_residual_j: f64,
    pub mass_residual_kg: f64,
    pub maximum_mass_residual_kg: f64,
    pub sum_absolute_mass_residual_kg: f64,
    pub maximum_cell_residual_j: f64,
    pub source_combination_residual_j: f64,
    pub gas_source_combination_residual_j: f64,
    /// Native heat-graph source application versus the independently observed
    /// external heat and complete physical patch debit.
    pub port_application_residual_j: f64,
    pub newton_iterations: usize,
    pub residual_evaluations: usize,
    pub gas_iterations: usize,
    pub maximum_accepted_residual_norm_j: f64,
}

/// Publish the original thermal graph, finite gas and film together.
pub struct DepositionWallFilmCandidate<'m, 'catalog, T> {
    thermal: DepositionSimulation<'m>,
    film: GasWallFilmInventory,
    gas: GasSpeciesBoundaryCandidate<'catalog>,
    payload: T,
    report: DepositionWallFilmReport,
}
impl<'m, 'catalog, T> DepositionWallFilmCandidate<'m, 'catalog, T> {
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.gas.gas()
    }
    pub const fn film(&self) -> GasWallFilmInventory {
        self.film
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.gas.catalog()
    }
    pub const fn payload(&self) -> &T {
        &self.payload
    }
    pub const fn report(&self) -> &DepositionWallFilmReport {
        &self.report
    }
    pub const fn gas_exchange(&self) -> &GasSpeciesBoundaryCandidate<'catalog> {
        &self.gas
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'m>,
        GasWallFilmInventory,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        T,
        DepositionWallFilmReport,
    ) {
        let (gas, catalog, _) = self.gas.into_parts();
        (
            self.thermal,
            self.film,
            gas,
            catalog,
            self.payload,
            self.report,
        )
    }
}

struct Prepared<'s, 'm, 'catalog> {
    heat: BackwardEulerStep<'s, 'm>,
    binding: DepositionWallFilmBinding<'s>,
    phase: GasWallFilmEndpoint<'s, 'catalog>,
    gas: &'s ClosedGasParcel,
    controls: DepositionWallFilmControls,
    external_gas_heat_j: f64,
    ambient_conductance_w_k: f64,
    gas_conductance_w_k: f64,
    gas_capacity_j_k: f64,
    vapor_mass_kg: f64,
    initial_wall_temperature_k: f64,
}
struct Ports<T> {
    payload: T,
    phase: GasWallFilmEndpointTrial,
    external: Option<EnergyAudit>,
    sensible_heat_into_gas_j: f64,
    source_roundoff_j: f64,
    gas_roundoff_j: f64,
}
struct Trial<'m, T> {
    heat: BackwardEulerTrial<'m, Ports<T>>,
    residual_j: Vec<f64>,
}

impl<'s, 'm, 'catalog> Prepared<'s, 'm, 'catalog> {
    #[allow(clippy::too_many_arguments)]
    fn new(
        initial: &'s DepositionSimulation<'m>,
        binding: DepositionWallFilmBinding<'s>,
        gas: &'s ClosedGasParcel,
        catalog: &'catalog ThermochemicalGasModel,
        target_s: f64,
        external_gas_heat_j: f64,
        controls: DepositionWallFilmControls,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitGasError> {
        checkpoint()?;
        controls.validate()?;
        let heat = BackwardEulerStep::new(initial, target_s, &mut *checkpoint)?;
        let area = initial.exposed_area_m2(binding.cell_id)?;
        let patch = binding.boundary.area_m2;
        if !nonnegative(patch)
            || patch > area
            || !nonnegative(binding.heat_transfer_coefficient_w_m2_k)
            || !external_gas_heat_j.is_finite()
        {
            return Err(ImplicitGasError::Invalid(
                "wall-film patch must fit the native exposed area with finite sources",
            ));
        }
        let remaining = finite(area - patch)?;
        if patch > 0.0 && remaining == area {
            return Err(ImplicitGasError::Invalid(
                "nonzero wall-film occupancy rounds away",
            ));
        }
        let ambient_conductance_w_k =
            finite(remaining * initial.cells()[&binding.cell_id].boundary.htc_w_m2_k)?;
        let gas_conductance_w_k = finite(patch * binding.heat_transfer_coefficient_w_m2_k)?;
        if patch > 0.0
            && binding.heat_transfer_coefficient_w_m2_k > 0.0
            && (gas_conductance_w_k == 0.0 || gas_conductance_w_k * heat.duration_s() == 0.0)
        {
            return Err(ImplicitGasError::Invalid(
                "nonzero spatial gas heat exposure underflows",
            ));
        }
        let initial_wall_temperature_k = initial.temperature_k(binding.cell_id)?;
        let h = liquid_h(
            binding,
            initial_wall_temperature_k,
            controls.temperature_tolerance_k,
        )?;
        let phase = GasWallFilmEndpoint::new(
            GasWallFilmSurface {
                law: binding.boundary.law,
                binding: binding.boundary.binding,
                film: binding.film,
                area_m2: patch,
                interface_liquid_enthalpy_j_kg: h,
                temperature_from_enthalpy: binding.boundary.liquid_temperature_from_enthalpy,
                mass_transfer_coefficient_m_s: binding.boundary.mass_transfer_coefficient_m_s,
            },
            gas,
            catalog,
            heat.duration_s(),
            controls.implicit.phase,
        )?;
        let state = gas.state().state();
        let gas_capacity_j_k = finite(gas.mass_kg() * state.specific_heat_cp / state.gamma)?;
        if !positive(gas_capacity_j_k) {
            return Err(ImplicitGasError::Invalid(
                "invalid original gas energy scaling",
            ));
        }
        let vapor_mass_kg = gas
            .species_mass_kg(&binding.boundary.binding.species)
            .map_err(GasEvaporationError::from)?;
        checkpoint()?;
        Ok(Self {
            heat,
            binding,
            phase,
            gas,
            controls,
            external_gas_heat_j,
            ambient_conductance_w_k,
            gas_conductance_w_k,
            gas_capacity_j_k,
            vapor_mass_kg,
            initial_wall_temperature_k,
        })
    }
    fn dimension(&self) -> usize {
        self.heat.cell_ids().len() + 2
    }
    fn project(&self, coordinates: &[f64]) -> Vec<f64> {
        let n = self.heat.cell_ids().len();
        let mut point = self.heat.resolvable_increments_j(&coordinates[..n]);
        for (value, initial, scale) in [
            (
                coordinates[n],
                self.vapor_mass_kg,
                self.controls.implicit.mass_coordinate_scale_j_kg,
            ),
            (
                coordinates[n + 1],
                self.gas.state().state().temperature,
                self.gas_capacity_j_k,
            ),
        ] {
            point.push(if value.is_finite() && initial + value / scale == initial {
                0.0
            } else {
                value
            });
        }
        point
    }
    fn evaluate<T>(
        &self,
        coordinates: &[f64],
        actual_gas: Option<&ClosedGasParcel>,
        resolve_heat: &mut impl FnMut(
            &DepositionSimulation<'m>,
            &ClosedGasParcel,
        )
            -> Result<BoundaryExchangeStepCoefficients<T>, ImplicitGasError>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Trial<'m, T>, ImplicitGasError> {
        let checkpoint = RefCell::new(checkpoint);
        let n = self.heat.cell_ids().len();
        let dt = self.heat.duration_s();
        let gas_t = actual_gas.map_or_else(
            || {
                finite(
                    self.gas.state().state().temperature
                        + coordinates[n + 1] / self.gas_capacity_j_k,
                )
            },
            |gas| Ok(gas.state().state().temperature),
        )?;
        let heat = self.heat.evaluate_with_boundary_exchanges::<_, ImplicitGasError>(&coordinates[..n], |thermal| {
            let mut poll = || (checkpoint.borrow_mut())();
            poll()?;
            let wall_t = thermal.temperature_k(self.binding.cell_id)?;
            if (wall_t - self.initial_wall_temperature_k).abs() > self.controls.implicit.phase.maximum_temperature_increment_k {
                return Err(GasEvaporationError::StepSizeLimit { maximum_duration_s: (0.5 * dt).next_down().max(0.0),
                    reason: "complete spatial wall-film temperature departure" }.into());
            }
            let h = liquid_h(self.binding, wall_t, self.controls.temperature_tolerance_k)?;
            let sensible = finite(dt * finite(self.gas_conductance_w_k * finite(wall_t - gas_t)?)?)?;
            let (gas_heat, gas_roundoff) = combine(self.external_gas_heat_j, sensible)?;
            let phase = phase_call(&mut poll, |cancel| match actual_gas {
                Some(gas) => self.phase.evaluate_at_gas(h, gas, gas_heat, cancel),
                None => self.phase.evaluate(h, coordinates[n] / self.controls.implicit.mass_coordinate_scale_j_kg, gas_t, gas_heat, cancel),
            })?;
            let mut coefficients = resolve_heat(thermal, phase.gas())?;
            if coefficients.transport.heat.boundaries.ambient_w_k.contains_key(&self.binding.cell_id)
                || coefficients.exchanges.ambient.contains_key(&self.binding.cell_id) {
                return Err(ImplicitGasError::Invalid("bound wall-film ambient channel cannot also be replaced by the heat callback"));
            }
            coefficients.transport.heat.boundaries.ambient_w_k.insert(self.binding.cell_id, self.ambient_conductance_w_k);
            // Only final physical acceptance needs the independent external
            // account. Newton/FD trials retain one native heat-graph traversal.
            // The fresh actual-gas replay uses the SAME callback result and
            // occupied ambient area without internal gas/phase source debits.
            let external = if actual_gas.is_some() {
                let physical = self.heat.evaluate_with_boundary_exchanges(&coordinates[..n],
                    |_| Ok::<_, ImplicitGasError>(without_payload(&coefficients)), &mut poll)?;
                Some(*physical.audit())
            } else { None };
            let (port, port_roundoff) = combine(phase.step().wall_heat_into_system_j, sensible)?;
            let power = finite(-port / dt)?;
            if port != 0.0 && power == 0.0 { return Err(ImplicitGasError::Invalid("nonzero wall-film heat cannot be represented as native power")); }
            let source = coefficients.transport.heat.bulk.sources_w.entry(self.binding.cell_id).or_default();
            let (updated, source_roundoff) = combine(*source, power)?;
            *source = updated;
            let source_roundoff_j = sum(&[-source_roundoff * dt, port_roundoff, power * dt + port])?;
            let BoundaryExchangeStepCoefficients { transport: TransportStepCoefficients {
                heat: BoundaryStepCoefficients { bulk: StepCoefficients { sources_w, conductances_w_k, payload }, boundaries }, transfers }, exchanges } = coefficients;
            Ok(BoundaryExchangeStepCoefficients { transport: TransportStepCoefficients {
                heat: BoundaryStepCoefficients { bulk: StepCoefficients { sources_w, conductances_w_k,
                    payload: Ports { payload, phase, external, sensible_heat_into_gas_j: sensible, source_roundoff_j, gas_roundoff_j: gas_roundoff } }, boundaries }, transfers }, exchanges })
        }, || (checkpoint.borrow_mut())())?;
        let mut residual_j = heat.residual_j().to_vec();
        residual_j.push(finite(
            heat.payload().phase.gas_mass_residual_kg()
                * self.controls.implicit.mass_coordinate_scale_j_kg,
        )?);
        residual_j.push(heat.payload().phase.gas_energy_residual_j());
        Ok(Trial { heat, residual_j })
    }
    fn check<T>(
        &self,
        trial: &Trial<'m, T>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitGasError> {
        checkpoint()?;
        self.phase.check_balance(&trial.heat.payload().phase)?;
        gate(
            trial.heat.maximum_cell_residual_j(),
            self.controls.implicit.thermal.cell_energy_tolerance_j,
        )?;
        gate(
            trial.heat.audit().residual_j,
            self.controls.implicit.thermal.total_energy_tolerance_j,
        )?;
        Ok(())
    }
    fn accept<T>(
        &self,
        trial: Trial<'m, T>,
        gas: GasSpeciesBoundaryCandidate<'catalog>,
        iterations: usize,
        evaluations: usize,
        norm: f64,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<DepositionWallFilmCandidate<'m, 'catalog, T>, ImplicitGasError> {
        self.check(&trial, checkpoint)?;
        let ports = trial.heat.payload();
        let external = ports.external.ok_or(ImplicitGasError::Invalid(
            "final spatial wall-film acceptance needs its independent external heat account",
        ))?;
        let phase = ports.phase.step();
        let physical_thermal_residual = sum(&[
            external.residual_j,
            phase.wall_heat_into_system_j,
            ports.sensible_heat_into_gas_j,
        ])?;
        let port_application_residual_j =
            sum(&[physical_thermal_residual, -trial.heat.audit().residual_j])?;
        let coupled_energy_residual_j = sum(&[
            physical_thermal_residual,
            phase.energy_residual_j,
            -ports.gas_roundoff_j,
        ])?;
        let gas_original_residual = gas.step().energy_residual_j.abs();
        gate(
            sum(&[phase.sum_absolute_energy_residual_j, gas_original_residual])?,
            self.controls.implicit.phase.gas.energy_tolerance_j,
        )?;
        let native_absolute = sum(&trial
            .heat
            .residual_j()
            .iter()
            .map(|value| value.abs())
            .collect::<Vec<_>>())?;
        let sum_absolute_energy_residual_j = sum(&[
            native_absolute.max(trial.heat.audit().residual_j.abs()),
            phase.sum_absolute_energy_residual_j,
            gas_original_residual,
            port_application_residual_j.abs(),
            ports.source_roundoff_j.abs(),
            ports.gas_roundoff_j.abs(),
        ])?;
        gate(
            sum_absolute_energy_residual_j.max(coupled_energy_residual_j.abs()),
            self.controls.implicit.thermal.total_energy_tolerance_j,
        )?;
        let mass_residual_kg = sum(&[
            gas.gas().mass_kg() - self.gas.mass_kg(),
            ports.phase.film().mass_kg() - self.binding.film.mass_kg(),
        ])?;
        let maximum_mass_residual_kg = ports
            .phase
            .maximum_mass_residual_kg()
            .max(mass_residual_kg.abs());
        let sum_absolute_mass_residual_kg = sum(&[
            ports.phase.maximum_mass_residual_kg(),
            gas.step().mass_residual_kg.abs(),
            gas.step().maximum_species_mass_residual_kg,
        ])?;
        let mass_tolerance = self
            .controls
            .implicit
            .thermal
            .mass_tolerance_kg
            .min(self.controls.implicit.phase.gas.mass_tolerance_kg);
        if maximum_mass_residual_kg.max(sum_absolute_mass_residual_kg) > mass_tolerance {
            return Err(ImplicitGasError::GasMassBalance {
                residual_kg: maximum_mass_residual_kg.max(sum_absolute_mass_residual_kg),
                tolerance_kg: mass_tolerance,
            });
        }
        let maximum_cell_residual_j = trial.heat.maximum_cell_residual_j();
        let (thermal, ports, thermal_audit) = trial.heat.accept(
            self.controls.implicit.thermal.cell_energy_tolerance_j,
            self.controls.implicit.thermal.total_energy_tolerance_j,
            &mut *checkpoint,
        )?;
        let report = DepositionWallFilmReport {
            thermal: thermal_audit,
            phase: ports.phase.step().clone(),
            external_thermal_heat_j: external.external_heat_j,
            external_gas_heat_j: self.external_gas_heat_j,
            sensible_heat_into_gas_j: ports.sensible_heat_into_gas_j,
            wall_heat_into_phase_system_j: ports.phase.step().wall_heat_into_system_j,
            coupled_energy_residual_j,
            sum_absolute_energy_residual_j,
            mass_residual_kg,
            maximum_mass_residual_kg,
            sum_absolute_mass_residual_kg,
            maximum_cell_residual_j,
            source_combination_residual_j: ports.source_roundoff_j,
            gas_source_combination_residual_j: -ports.gas_roundoff_j,
            port_application_residual_j,
            newton_iterations: iterations,
            residual_evaluations: evaluations,
            gas_iterations: gas.step().iterations,
            maximum_accepted_residual_norm_j: norm,
        };
        checkpoint()?;
        Ok(DepositionWallFilmCandidate {
            thermal,
            film: ports.phase.film(),
            gas,
            payload: ports.payload,
            report,
        })
    }
}

fn without_payload<T>(
    value: &BoundaryExchangeStepCoefficients<T>,
) -> BoundaryExchangeStepCoefficients<()> {
    BoundaryExchangeStepCoefficients {
        transport: TransportStepCoefficients {
            heat: BoundaryStepCoefficients {
                bulk: StepCoefficients {
                    sources_w: value.transport.heat.bulk.sources_w.clone(),
                    conductances_w_k: value.transport.heat.bulk.conductances_w_k.clone(),
                    payload: (),
                },
                boundaries: value.transport.heat.boundaries.clone(),
            },
            transfers: value.transport.transfers.clone(),
        },
        exchanges: value.exchanges.clone(),
    }
}
fn liquid_h(
    binding: DepositionWallFilmBinding<'_>,
    temperature: f64,
    tolerance: f64,
) -> Result<f64, ImplicitGasError> {
    let h = (binding.boundary.liquid_enthalpy_at_temperature)(temperature)
        .map_err(GasEvaporationError::from)?;
    let t = (binding.boundary.liquid_temperature_from_enthalpy)(h)
        .map_err(GasEvaporationError::from)?;
    if !h.is_finite() || !t.is_finite() || (t - temperature).abs() > tolerance {
        return Err(ImplicitGasError::Invalid(
            "original liquid h(T)/T(h) disagrees with spatial wall temperature",
        ));
    }
    Ok(h)
}
fn phase_call<T>(
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    operation: impl FnOnce(&mut dyn FnMut() -> bool) -> Result<T, GasEvaporationError>,
) -> Result<T, ImplicitGasError> {
    let mut interrupted = None;
    let result = operation(&mut || match checkpoint() {
        Ok(()) => false,
        Err(error) => {
            interrupted = Some(error);
            true
        }
    });
    if let Some(error) = interrupted {
        return Err(error.into());
    }
    Ok(result?)
}
fn finite(value: f64) -> Result<f64, ImplicitGasError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ImplicitGasError::Invalid(
            "spatial wall-film arithmetic is unrepresentable",
        ))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn sum(values: &[f64]) -> Result<f64, ImplicitGasError> {
    let (mut total, mut correction) = (0.0, 0.0);
    for &value in values {
        let (next, remainder) = fs_math::eft::two_sum(total, value);
        total = finite(next)?;
        correction = finite(correction + remainder)?;
    }
    finite(total + correction)
}
fn combine(a: f64, b: f64) -> Result<(f64, f64), ImplicitGasError> {
    let (total, roundoff) = fs_math::eft::two_sum(a, b);
    finite(total)?;
    finite(roundoff)?;
    if (a != 0.0 && total == b) || (b != 0.0 && total == a) {
        return Err(ImplicitGasError::Invalid(
            "nonzero spatial wall-film source rounds away",
        ));
    }
    Ok((total, roundoff))
}
fn gate(residual_j: f64, tolerance_j: f64) -> Result<(), ImplicitGasError> {
    if !residual_j.is_finite() || residual_j.abs() > tolerance_j {
        Err(ImplicitGasError::CoupledEnergyBalance {
            residual_j,
            tolerance_j,
        })
    } else {
        Ok(())
    }
}

/// The shared meter produces budget or bookkeeping errors, never a constitutive
/// refusal. Original phase/thermal errors bypass this numerical conversion.
pub(super) fn work_error(error: WallFilmImplicitError) -> ImplicitGasError {
    match error {
        WallFilmImplicitError::EvaluationBudget => {
            DepositionError::Budget("spatial wall-film residual evaluations").into()
        }
        WallFilmImplicitError::IterationBudget => {
            DepositionError::Budget("spatial wall-film Newton iterations").into()
        }
        WallFilmImplicitError::GasIterationBudget => {
            DepositionError::Budget("spatial wall-film gas iterations").into()
        }
        WallFilmImplicitError::Physical(WallFilmThermalError::Invalid(message)) => {
            ImplicitGasError::Invalid(message)
        }
        WallFilmImplicitError::Numerical(error) => ImplicitGasError::Numerical(error),
        WallFilmImplicitError::NonConverged(report) => ImplicitGasError::NonConverged(report),
        WallFilmImplicitError::Interrupted(progress) => ImplicitGasError::Interrupted(progress),
        _ => ImplicitGasError::Invalid("shared implicit meter returned a non-work refusal"),
    }
}
