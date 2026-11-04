//! Implicit heat exchange with original constrained mixture phase equilibrium.
//!
//! The heat coordinate closes `Q = dt G [Tw(-Q) - Tg(Q)]`. Every gas trial
//! uses the original U,V or prescribed-pressure H owner; no equilibrium Cp/Cv
//! or separate latent source is constructed. Original finite wall storage
//! remains in its declared h/u coordinate. This backward-Euler lumped interface needs time refinement and
//! a uniform-inventory justification; it is not a spatial or kinetic model.

use crate::multicomponent::wall_exchange::{
    MaterialWallCaloric, WallCellReport, WallContact, WallControl, WallDrag,
    WallExchangeError as Error, WallThermal, solve_scalar,
};
use fs_material::thermochemical_gas::nonideal::flash::{
    ThermochemicalFlashError,
    equilibrium::{PhaseEnthalpyControls, volume::PhaseEnergyVolumeControls},
    parcel::{
        IsochoricMixtureEquilibriumParcel, MixtureEquilibriumHeatStep,
        isobaric::{IsobaricMixtureEquilibriumParcel, IsobaricMixtureHeatStep},
    },
};
use fs_qty::Time;

mod pair;
pub use pair::{
    IsobaricMixturePairThermalControl, IsobaricMixturePairThermalStep,
    MixedMixturePairThermalControl, MixedMixturePairThermalStep,
    MixturePairThermalControl, MixturePairThermalStep,
    exchange_isobaric_mixture_pair, exchange_isochoric_isobaric_mixture_pair,
    exchange_isochoric_mixture_pair,
};

/// Explicit heat bracket and source work for one original-law thermal interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureThermalControl<C = PhaseEnergyVolumeControls> {
    /// Signed heat into the gas [J]; finite increasing endpoints must contain
    /// zero. Intersected with the passive initial-flux and wall-source limits.
    /// Every remaining endpoint must be supported by the original caloric solver.
    pub heat_bracket_j: [f64; 2],
    /// Original species/phase policy, physical bounds and per-trial work limits.
    pub equilibrium: C,
    /// Physical joule residuals and safeguarded scalar iteration limit.
    pub thermal: WallControl,
    /// Total heat-coordinate queries, including zero and bracket endpoints.
    pub maximum_heat_evaluations: usize,
    /// Total original PT selections across all heat trials.
    pub maximum_pt_evaluations: usize,
}

impl<C> MixtureThermalControl<C> {
    pub fn new(heat_bracket_j: [f64; 2], equilibrium: C) -> Self {
        Self {
            heat_bracket_j,
            equilibrium,
            thermal: WallControl {
                absolute_energy_tolerance_j: 1e-5,
                relative_energy_tolerance: 1e-9,
                maximum_iterations: 96,
            },
            maximum_heat_evaluations: 132,
            maximum_pt_evaluations: 65_536,
        }
    }

    fn validate_common(&self) -> Result<(), Error> {
        self.thermal.validate()?;
        let [lo, hi] = self.heat_bracket_j;
        if !lo.is_finite()
            || !hi.is_finite()
            || lo > 0.0
            || hi < 0.0
            || lo >= hi
            || !(1..=132).contains(&self.maximum_heat_evaluations)
            || !(1..=1_048_576).contains(&self.maximum_pt_evaluations)
        {
            return Err(Error::Invalid("mixture thermal bracket or work limits"));
        }
        Ok(())
    }
}

impl MixtureThermalControl {
    pub fn validate(self) -> Result<(), Error> {
        self.equilibrium.validate()?;
        self.validate_common()
    }
}

/// The same bounded heat-transfer solve at a prescribed constant pressure.
/// The original PH owner retains latent enthalpy and the resulting expansion.
pub type IsobaricMixtureThermalControl = MixtureThermalControl<PhaseEnthalpyControls>;

impl IsobaricMixtureThermalControl {
    pub fn validate(self) -> Result<(), Error> {
        self.equilibrium.validate()?;
        self.validate_common()
    }
}

/// Accepted thermal account and original mixture endpoint, without a second
/// inventory. Species stay fixed; the original constraint determines volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureThermalStep<R = MixtureEquilibriumHeatStep> {
    /// Existing gas/wall heat accounts; drag and support impulse are zero.
    pub heat: WallCellReport,
    /// Actual source heat/storage/work accounts at the accepted endpoint.
    /// Its work counts describe the accepted trial only; whole-call counts follow.
    pub equilibrium: R,
    /// Initial gas and wall/bath temperatures [K].
    pub initial_temperature_k: [f64; 2],
    /// Accepted gas and wall/bath temperatures [K].
    pub endpoint_temperature_k: [f64; 2],
    pub heat_evaluations: usize,
    pub pt_evaluations: usize,
    /// Original pressure-root samples; zero for the prescribed-pressure owner.
    pub pressure_evaluations: usize,
    /// Sum of the original source EOS ceilings, not an observed evaluation count.
    pub eos_evaluation_upper_bound: usize,
}

/// Accepted constant-pressure heat exchange, including original Delta H,
/// Delta U and outward p Delta V in `equilibrium`. The common heat report's
/// mechanical work is negative when the mixture exports expansion work.
pub type IsobaricMixtureThermalStep = MixtureThermalStep<IsobaricMixtureHeatStep>;

// Only the native caloric constraint varies. Heat transfer, wall response,
// bounded scalar iteration, physical acceptance and publication stay shared.
trait ThermalParcel: Clone {
    type Control: Copy;
    type Report: Copy;
    fn temperature(&self) -> f64;
    fn mass(&self) -> f64;
    fn validate(control: Self::Control) -> Result<(), Error>;
    fn refine(control: &mut Self::Control, tolerance_j_kg: f64);
    fn resolve(
        &mut self,
        heat: f64,
        control: Self::Control,
        remaining_pt: usize,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<CaloricReport<Self::Report>, ThermochemicalFlashError>;
}

struct CaloricReport<R> {
    source: R,
    heat_change_j: f64,
    internal_energy_change_j: f64,
    boundary_work_j: f64,
    entropy_change_j_k: f64,
    pt: usize,
    pressure: usize,
    eos: usize,
}

impl ThermalParcel for IsochoricMixtureEquilibriumParcel {
    type Control = PhaseEnergyVolumeControls;
    type Report = MixtureEquilibriumHeatStep;
    fn temperature(&self) -> f64 {
        self.state().temperature_k()
    }
    fn mass(&self) -> f64 {
        self.mass_kg()
    }
    fn validate(control: Self::Control) -> Result<(), Error> {
        Ok(control.validate()?)
    }
    fn refine(control: &mut Self::Control, tolerance_j_kg: f64) {
        control.internal_energy_tolerance_j_kg =
            control.internal_energy_tolerance_j_kg.min(tolerance_j_kg);
    }
    fn resolve(
        &mut self,
        heat: f64,
        mut control: Self::Control,
        remaining_pt: usize,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<CaloricReport<Self::Report>, ThermochemicalFlashError> {
        if heat != 0.0 {
            control.maximum_pt_evaluations = control.maximum_pt_evaluations.min(remaining_pt);
        }
        let source = self.add_heat(heat, control, cancelled)?;
        Ok(CaloricReport {
            heat_change_j: source.stored_energy_change_j,
            internal_energy_change_j: source.stored_energy_change_j,
            boundary_work_j: 0.0,
            entropy_change_j_k: source.entropy_change_j_k,
            pt: source.pt_evaluations,
            pressure: source.pressure_evaluations,
            eos: source.eos_evaluation_upper_bound,
            source,
        })
    }
}

impl ThermalParcel for IsobaricMixtureEquilibriumParcel {
    type Control = PhaseEnthalpyControls;
    type Report = IsobaricMixtureHeatStep;
    fn temperature(&self) -> f64 {
        self.state().temperature_k()
    }
    fn mass(&self) -> f64 {
        self.mass_kg()
    }
    fn validate(control: Self::Control) -> Result<(), Error> {
        Ok(control.validate()?)
    }
    fn refine(control: &mut Self::Control, tolerance_j_kg: f64) {
        control.enthalpy.enthalpy_tolerance_j_kg =
            control.enthalpy.enthalpy_tolerance_j_kg.min(tolerance_j_kg);
    }
    fn resolve(
        &mut self,
        heat: f64,
        mut control: Self::Control,
        remaining_pt: usize,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<CaloricReport<Self::Report>, ThermochemicalFlashError> {
        if heat != 0.0 {
            control.enthalpy.maximum_flash_evaluations =
                control.enthalpy.maximum_flash_evaluations.min(remaining_pt);
        }
        let source = self.add_heat(heat, control, cancelled)?;
        Ok(CaloricReport {
            heat_change_j: source.enthalpy_change_j,
            internal_energy_change_j: source.stored_energy_change_j,
            boundary_work_j: source.boundary_work_j,
            entropy_change_j_k: source.entropy_change_j_k,
            pt: source.pt_evaluations,
            pressure: 0,
            eos: source.eos_evaluation_upper_bound,
            source,
        })
    }
}

#[derive(Default)]
struct Work {
    heat: usize,
    pt: usize,
    pressure: usize,
    eos: usize,
}
impl Work {
    fn error(&self) -> Error {
        Error::MixtureThermalWork {
            heat_evaluations: self.heat,
            pt_evaluations: self.pt,
        }
    }
    fn begin<C>(&mut self, controls: &MixtureThermalControl<C>, heat: f64) -> Result<usize, Error> {
        if self.heat == controls.maximum_heat_evaluations
            || (heat != 0.0 && self.pt == controls.maximum_pt_evaluations)
        {
            return Err(self.error());
        }
        self.heat += 1;
        Ok(controls.maximum_pt_evaluations - self.pt)
    }
    fn add(&mut self, pt: usize, pressure: usize, eos: usize) -> Result<(), Error> {
        self.pt = self.pt.checked_add(pt).ok_or_else(|| self.error())?;
        self.pressure = self
            .pressure
            .checked_add(pressure)
            .ok_or_else(|| self.error())?;
        self.eos = self.eos.checked_add(eos).ok_or_else(|| self.error())?;
        Ok(())
    }
    fn source_error(&mut self, error: ThermochemicalFlashError) -> Error {
        if let ThermochemicalFlashError::WorkLimit { flash_evaluations } = error {
            self.pt += flash_evaluations;
            self.error()
        } else {
            error.into()
        }
    }
}

struct WallSample {
    state: WallThermal,
    temperature: f64,
    change: f64,
}

struct Trial<P: ThermalParcel> {
    heat: f64,
    parcel: P,
    equilibrium: CaloricReport<P::Report>,
    wall: WallSample,
}

fn trial<P: ThermalParcel>(
    parcel: &P,
    wall: WallThermal,
    source: Option<&MaterialWallCaloric<'_>>,
    heat: f64,
    controls: MixtureThermalControl<P::Control>,
    work: &mut Work,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<Trial<P>, Error> {
    let remaining = work.begin(&controls, heat)?;
    let mut candidate = parcel.clone();
    let equilibrium = candidate
        .resolve(heat, controls.equilibrium, remaining, cancelled)
        .map_err(|e| work.source_error(e))?;
    work.add(equilibrium.pt, equilibrium.pressure, equilibrium.eos)?;
    Ok(Trial {
        heat,
        parcel: candidate,
        equilibrium,
        wall: wall_after(wall, source, heat)?,
    })
}

fn wall_after(
    state: WallThermal,
    source: Option<&MaterialWallCaloric<'_>>,
    heat_into_gas: f64,
) -> Result<WallSample, Error> {
    if let Some(source) = source {
        let next = source.after(
            WallContact {
                drag: WallDrag::None,
                thermal: state,
            },
            -heat_into_gas,
        )?;
        return Ok(WallSample {
            state: if heat_into_gas == 0.0 {
                state
            } else {
                next.contact.thermal
            },
            temperature: next.temperature_k,
            change: if heat_into_gas == 0.0 {
                0.0
            } else {
                next.energy_change_j
            },
        });
    }
    match state {
        WallThermal::Reservoir { temperature_k, .. } => Ok(WallSample {
            state,
            temperature: temperature_k,
            change: 0.0,
        }),
        WallThermal::Lumped {
            temperature_k,
            heat_capacity_j_per_k,
            conductance_w_per_k,
        } => {
            let next = temperature_k - heat_into_gas / heat_capacity_j_per_k;
            if !positive(next) {
                return Err(Error::Invalid("mixture finite-wall temperature"));
            }
            Ok(WallSample {
                state: WallThermal::Lumped {
                    temperature_k: next,
                    heat_capacity_j_per_k,
                    conductance_w_per_k,
                },
                temperature: next,
                change: heat_capacity_j_per_k * (next - temperature_k),
            })
        }
        _ => Err(Error::Invalid(
            "mixture thermal wall requires a source or reservoir",
        )),
    }
}

/// Exchange heat with an existing rigid phase-equilibrium mixture and publish
/// gas plus wall atomically. `None` admits an authored constant-capacity wall
/// or a held-temperature reservoir; a caloric wall requires its original source.
/// Conductance is the complete caller-declared W/K interface, held this interval.
///
/// Both endpoint heat/storage equations are checked in actual represented
/// joules. Original enthalpy walls retain their pressure-work convention. The
/// gas entropy increment is available, but an arbitrary wall need not provide
/// an entropy law, so this API makes no combined entropy-production theorem.
/// Cancellation, source/phase refusal and exhausted work leave both inputs intact.
pub fn exchange_mixture_parcel(
    parcel: &mut IsochoricMixtureEquilibriumParcel,
    wall: &mut WallThermal,
    wall_source: Option<&MaterialWallCaloric<'_>>,
    duration: Time,
    controls: MixtureThermalControl,
    cancelled: impl FnMut() -> bool,
) -> Result<MixtureThermalStep, Error> {
    exchange(parcel, wall, wall_source, duration, controls, cancelled)
}

/// Exchange heat while the closed phase-equilibrium mixture expands against a
/// prescribed constant pressure. The original PH owner determines T, phase
/// amounts and volume at every trial; heat enters Delta H exactly once.
///
/// Uses the same conservative implicit wall/contact solve as the rigid-vessel
/// path. Actual Delta U, outward p Delta V and wall storage close the complete
/// first law before gas and wall publish together. The mechanical load supplies
/// or receives boundary work; it is not a second heat source. Heat-transfer
/// geometry and conductance are held over this interval. Piston inertia, phase kinetics and
/// resolved spatial heat transport remain outside this lumped equilibrium model.
pub fn exchange_isobaric_mixture_parcel(
    parcel: &mut IsobaricMixtureEquilibriumParcel,
    wall: &mut WallThermal,
    wall_source: Option<&MaterialWallCaloric<'_>>,
    duration: Time,
    controls: IsobaricMixtureThermalControl,
    cancelled: impl FnMut() -> bool,
) -> Result<IsobaricMixtureThermalStep, Error> {
    exchange(parcel, wall, wall_source, duration, controls, cancelled)
}

fn exchange<P: ThermalParcel>(
    parcel: &mut P,
    wall: &mut WallThermal,
    wall_source: Option<&MaterialWallCaloric<'_>>,
    duration: Time,
    mut controls: MixtureThermalControl<P::Control>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixtureThermalStep<P::Report>, Error> {
    poll(&mut cancelled)?;
    P::validate(controls.equilibrium)?;
    controls.validate_common()?;
    if !duration.value().is_finite() || duration.value() < 0.0 {
        return Err(Error::Invalid(
            "nonnegative finite mixture thermal duration required",
        ));
    }
    let contact = WallContact {
        drag: WallDrag::None,
        thermal: *wall,
    };
    contact.validate()?;
    let conductance = match *wall {
        WallThermal::Reservoir {
            conductance_w_per_k,
            ..
        }
        | WallThermal::Lumped {
            conductance_w_per_k,
            ..
        }
        | WallThermal::Caloric {
            conductance_w_per_k,
            ..
        } => conductance_w_per_k,
        WallThermal::Adiabatic => {
            return Err(Error::Invalid(
                "use zero conductance for an insulated mixture wall",
            ));
        }
    };
    if let Some(source) = wall_source {
        source.bind_contact(contact)?;
    }
    let initial = wall_after(*wall, wall_source, 0.0)?;
    let tg0 = parcel.temperature();
    let coefficient = duration.value() * conductance;
    let initial_flux = coefficient * (initial.temperature - tg0);
    if !coefficient.is_finite()
        || !initial_flux.is_finite()
        || (coefficient == 0.0 && duration.value() != 0.0 && conductance != 0.0)
        || (initial_flux == 0.0 && coefficient != 0.0 && initial.temperature != tg0)
    {
        return Err(Error::Invalid("mixture thermal flux is unrepresentable"));
    }
    let mut work = Work::default();
    let mut accepted = None;
    let (heat, iterations) = if initial_flux == 0.0 {
        accepted = Some(trial(
            parcel,
            *wall,
            wall_source,
            0.0,
            controls,
            &mut work,
            &mut cancelled,
        )?);
        (0.0, 0)
    } else {
        controls.thermal.absolute_energy_tolerance_j = controls
            .thermal
            .absolute_energy_tolerance_j
            .min(0.125 * initial_flux.abs());
        P::refine(
            &mut controls.equilibrium,
            0.125 * controls.thermal.absolute_energy_tolerance_j / parcel.mass(),
        );
        P::validate(controls.equilibrium)?;
        controls.validate_common()?;
        let [mut lo, mut hi] = controls.heat_bracket_j;
        lo = lo.max(initial_flux.min(0.0));
        hi = hi.min(initial_flux.max(0.0));
        if let Some(source) = wall_source {
            let limits = source.heat_limits(contact, tg0)?;
            lo = lo.max(limits[0]);
            hi = hi.min(limits[1]);
        } else if let WallThermal::Lumped {
            heat_capacity_j_per_k,
            ..
        } = *wall
        {
            let equalized = heat_capacity_j_per_k * (initial.temperature - tg0);
            lo = lo.max(equalized.min(0.0));
            hi = hi.min(equalized.max(0.0));
        }
        if lo >= hi {
            return Err(Error::Invalid("mixture thermal heat bracket is empty"));
        }
        let mut previous: Option<(f64, f64)> = None;
        let mut equation = |q: f64, cancel: &mut dyn FnMut() -> bool| {
            let candidate = trial(parcel, *wall, wall_source, q, controls, &mut work, cancel)?;
            let tg = candidate.parcel.temperature();
            let tw = candidate.wall.temperature;
            let flux = coefficient * (tw - tg);
            let residual = q - flux;
            if !residual.is_finite() {
                return Err(Error::Invalid("mixture thermal trial flux"));
            }
            // A numerical secant accelerates this bounded scalar search. It is
            // neither an equilibrium heat capacity nor a physical tangent.
            let slope = previous.and_then(|(old_q, old_r)| {
                let slope = (residual - old_r) / (q - old_q);
                positive(slope).then_some(slope)
            });
            previous = Some((q, residual));
            // The scalar owner returns its last evaluated value. Retain that
            // full admitted trial so acceptance performs no second caloric solve.
            accepted = Some(candidate);
            Ok((residual, slope, q.abs().max(flux.abs())))
        };
        let scalar = WallControl {
            absolute_energy_tolerance_j: 0.5 * controls.thermal.absolute_energy_tolerance_j,
            relative_energy_tolerance: 0.5 * controls.thermal.relative_energy_tolerance,
            ..controls.thermal
        };
        scalar.validate()?;
        solve_scalar(lo, hi, 0.0, scalar, &mut equation, &mut cancelled)?
    };
    let accepted = accepted.ok_or(Error::Invalid("mixture thermal endpoint missing"))?;
    if accepted.heat.to_bits() != heat.to_bits() {
        return Err(Error::Invalid(
            "mixture thermal endpoint disagrees with scalar root",
        ));
    }
    let candidate = accepted.parcel;
    let equilibrium = accepted.equilibrium;
    let endpoint = accepted.wall;
    let actual = equilibrium.heat_change_j;
    let external = if matches!(*wall, WallThermal::Reservoir { .. }) {
        actual
    } else {
        0.0
    };
    let energy_residual = (equilibrium.internal_energy_change_j + endpoint.change)
        + equilibrium.boundary_work_j
        - external;
    let flux = coefficient * (endpoint.temperature - candidate.temperature());
    let thermal_residual = actual - flux;
    let tolerance = controls.thermal.tolerance(actual.abs().max(flux.abs()));
    if !energy_residual.is_finite()
        || !thermal_residual.is_finite()
        || energy_residual.abs() > tolerance
        || thermal_residual.abs() > tolerance
    {
        return Err(Error::Balance {
            energy_residual_j: energy_residual,
            thermal_residual_j: thermal_residual,
        });
    }
    let report = MixtureThermalStep {
        heat: WallCellReport {
            elapsed: duration,
            wall_impulse_ns: 0.0,
            mechanical_work_to_gas_j: if equilibrium.boundary_work_j == 0.0 {
                0.0
            } else {
                -equilibrium.boundary_work_j
            },
            drag_heating_j: 0.0,
            heat_to_gas_j: actual,
            wall_energy_change_j: endpoint.change,
            external_heat_to_gas_j: external,
            energy_residual_j: energy_residual,
            thermal_residual_j: thermal_residual,
            iterations,
        },
        equilibrium: equilibrium.source,
        initial_temperature_k: [tg0, initial.temperature],
        endpoint_temperature_k: [candidate.temperature(), endpoint.temperature],
        heat_evaluations: work.heat,
        pt_evaluations: work.pt,
        pressure_evaluations: work.pressure,
        eos_evaluation_upper_bound: work.eos,
    };
    poll(&mut cancelled)?;
    *parcel = candidate;
    *wall = endpoint.state;
    Ok(report)
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
