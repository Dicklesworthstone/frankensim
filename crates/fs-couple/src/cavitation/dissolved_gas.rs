//! Finite solvent heat closes the original two-film gas absorption history.
//!
//! One interval uses the original exact isothermal mass update at T0, then
//! solves the common temperature from total internal energy. This conservative
//! first-order split requires time refinement: a small energy residual does
//! not bound its kinetic error. Refreshed original H(T), D(T), gas pressure and
//! calorics feed the next interval. No additional dissolution heat is added.
//!
//! Both phases are uniform and share one temperature throughout the reduced
//! model. Volumes are fixed, the solvent is incompressible and nonvolatile,
//! and dissolved partial volume is negligible. Source charts must apply over
//! the declared pressure/composition interval; an authored gas-film coefficient
//! must also apply over the temperature change. Henry and caloric sources must
//! be physically compatible; this composition does not certify that empirical
//! relation or infer a heat of solution from an equilibrium fit. It supplies no finite
//! interfacial heat-transfer resistance, moving bubble or solvent evaporation.
//! Heat effects in absorption: ETH Zurich, nonisothermal absorption,
//! https://www.hyper-tvt.ethz.ch/abstrip-nonisothermal.html .

use crate::cavity_phs::piston::thermal::feedback;
use fs_material::cavitation::dissolved_gas::{
    DissolutionFilm, DissolutionSample, DissolvedGasControls, DissolvedGasError,
    DissolvedGasInventory, HenryDissolution,
    mixture::{DissolutionMixtureStep, MAX_DISSOLVING_SPECIES, advance_mixture},
};
use fs_material::phase::{EquilibriumPhaseState, PhaseStateError, ResolvedHeatCapacityCurve};
use fs_material::thermochemical_gas::energy::{GasEnergyError, parcel::ClosedGasParcel};
use fs_math::dd::Dd;
use fs_qty::Temperature;
use fs_thermochem::mixture::FrozenIdealGasMixtureErrorV1;

/// Original pure-solvent sensible storage in one fixed liquid volume.
/// Its chart must remain entirely liquid with constant density. The solutes'
/// separate original partial caloric charts supply their storage exactly once.
#[derive(Clone, Debug, PartialEq)]
pub struct IncompressibleSolvent {
    calorics: ResolvedHeatCapacityCurve,
    volume_m3: f64,
    mass_kg: f64,
}
impl IncompressibleSolvent {
    pub fn new(
        calorics: ResolvedHeatCapacityCurve,
        volume_m3: f64,
    ) -> Result<Self, ThermalDissolutionError> {
        let knots = calorics.curve().knots();
        let density = knots[0].bulk_density_kg_m3;
        let mass_kg = density * volume_m3;
        if !positive(volume_m3)
            || !positive(mass_kg)
            || knots
                .iter()
                .any(|k| k.liquid_mass_fraction != 1.0 || k.bulk_density_kg_m3 != density)
        {
            return Err(invalid(
                "positive fixed solvent volume and one incompressible liquid chart required",
            ));
        }
        Ok(Self {
            calorics,
            volume_m3,
            mass_kg,
        })
    }
    pub const fn calorics(&self) -> &ResolvedHeatCapacityCurve {
        &self.calorics
    }
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
}

/// Whole-interval source, physical-energy and bounded thermal-root controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalDissolutionControls {
    pub exchange: DissolvedGasControls,
    /// Absolute source-supported temperatures [K], containing the initial T.
    pub temperature_bracket_k: [f64; 2],
    /// Physical step restriction, not a temporal truncation-error estimate.
    pub maximum_temperature_increment_k: f64,
    /// Complete gas + dissolved species + solvent energy residual [J].
    pub energy_tolerance_j: f64,
    /// At most 128 safeguarded scalar iterations, plus bounded endpoint trials.
    pub maximum_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ThermalDissolutionError {
    Invalid(&'static str),
    Dissolution(DissolvedGasError),
    Caloric(PhaseStateError),
    Gas(GasEnergyError),
    Thermochemistry(FrozenIdealGasMixtureErrorV1),
    OutsideBracket {
        lower_residual_j: f64,
        upper_residual_j: f64,
    },
    NotConverged {
        iterations: usize,
        residual_j: f64,
    },
    Energy {
        residual_j: f64,
        tolerance_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for ThermalDissolutionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite thermal gas absorption: {self:?}")
    }
}
impl std::error::Error for ThermalDissolutionError {}
impl From<DissolvedGasError> for ThermalDissolutionError {
    fn from(e: DissolvedGasError) -> Self {
        if e == DissolvedGasError::Cancelled {
            Self::Cancelled
        } else {
            Self::Dissolution(e)
        }
    }
}
impl From<GasEnergyError> for ThermalDissolutionError {
    fn from(e: GasEnergyError) -> Self {
        if e == GasEnergyError::Cancelled {
            Self::Cancelled
        } else {
            Self::Gas(e)
        }
    }
}
impl From<PhaseStateError> for ThermalDissolutionError {
    fn from(e: PhaseStateError) -> Self {
        Self::Caloric(e)
    }
}
impl From<FrozenIdealGasMixtureErrorV1> for ThermalDissolutionError {
    fn from(e: FrozenIdealGasMixtureErrorV1) -> Self {
        Self::Thermochemistry(e)
    }
}

/// Accepted material/species/temperature history and separate physical accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalDissolutionStep {
    pub initial_temperature_k: f64,
    pub final_temperature_k: f64,
    pub solvent_initial: EquilibriumPhaseState,
    pub solvent_state: EquilibriumPhaseState,
    /// Original kinetic substep at the initial T. Its thermostats are a
    /// numerical decomposition; no thermostat is attached to the closed pair.
    pub isothermal_exchange: DissolutionMixtureStep,
    /// Actual final bulk/interface sources, in canonical species order.
    pub endpoint_sources: Vec<DissolutionSample>,
    pub gas_sensible_energy_change_j: f64,
    pub dissolved_sensible_energy_change_j: f64,
    pub solvent_energy_change_j: f64,
    pub external_heat_j: f64,
    /// Actual constituent energy increments minus the prescribed heat.
    pub energy_residual_j: f64,
    pub temperature_trials: usize,
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn invalid(s: &'static str) -> ThermalDissolutionError {
    ThermalDissolutionError::Invalid(s)
}
fn finite(x: f64) -> Result<f64, ThermalDissolutionError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(invalid("nonfinite thermal absorption arithmetic"))
    }
}
fn total(values: impl IntoIterator<Item = f64>) -> Result<f64, ThermalDissolutionError> {
    finite(
        values
            .into_iter()
            .fold(Dd::ZERO, |s, x| s + Dd::from_f64(x))
            .to_f64(),
    )
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ThermalDissolutionError> {
    if cancelled() {
        Err(ThermalDissolutionError::Cancelled)
    } else {
        Ok(())
    }
}

/// Advance one finite, dilute gas/solution pair with common-temperature heat.
///
/// Slices retain one shared liquid volume/clock and the existing mixture
/// owner's source/species rules. `external_heat_j = 0` is adiabatic for the
/// COMBINED pair. The closed gas alone may heat or cool as solute moves.
/// Solvent density must not be inferred from a solute partial caloric chart.
///
/// The original mass solve is frozen at T0. Its exact integrated internal-energy
/// change plus source-native sensible increments closes the final temperature:
/// `Delta U_transfer(T0) + Delta U_sensible(T0,T1) = Q_external`.
/// The full species/caloric/source state publishes only after the final source,
/// mixture-interface, energy and cancellation checks. Partial proposals never
/// advance caller clocks, sources or inventory. Time and temperature step limits
/// do not certify splitting accuracy; refine both for a predictive trajectory.
pub fn advance_thermal_mixture(
    laws: &mut [HenryDissolution],
    gas: &mut ClosedGasParcel,
    liquids: &mut [DissolvedGasInventory],
    films: &[DissolutionFilm],
    solvent: &IncompressibleSolvent,
    end_time_s: f64,
    external_heat_j: f64,
    controls: ThermalDissolutionControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ThermalDissolutionStep, ThermalDissolutionError> {
    poll(&mut cancelled)?;
    let t0 = gas.state().state().temperature;
    let [lo, hi] = controls.temperature_bracket_k;
    if laws.is_empty()
        || laws.len() > MAX_DISSOLVING_SPECIES
        || laws.len() != liquids.len()
        || laws.len() != films.len()
        || liquids
            .iter()
            .any(|l| l.liquid_volume_m3() != solvent.volume_m3)
        || !positive(lo)
        || !positive(hi)
        || lo >= hi
        || t0 < lo
        || t0 > hi
        || lo < controls.exchange.gas.gas.minimum_temperature_k
        || hi > controls.exchange.gas.gas.maximum_temperature_k
        || !positive(controls.maximum_temperature_increment_k)
        || !positive(controls.energy_tolerance_j)
        || !external_heat_j.is_finite()
        || !(1..=128).contains(&controls.maximum_iterations)
    {
        return Err(invalid(
            "matching finite solvent/species, supported thermal bracket and bounded controls required",
        ));
    }
    let solvent_initial = solvent.calorics.state_at_temperature(t0)?;
    // Admit the full offered caloric bracket without assigning empirical
    // authority to independent Henry or film coefficients.
    for t in [lo, hi] {
        solvent.calorics.state_at_temperature(t)?;
    }
    for law in laws.iter() {
        for t in [lo, hi] {
            law.calorics().partial_enthalpy.state_at_temperature(t)?;
        }
    }
    let mut next_gas = gas.clone();
    let mut next_liquids = liquids.to_vec();
    let exchange = advance_mixture(
        laws,
        &mut next_gas,
        &mut next_liquids,
        films,
        end_time_s,
        controls.exchange,
        &mut cancelled,
    )?;
    let mut order: Vec<_> = (0..laws.len()).collect();
    order.sort_by(|&a, &b| laws[a].species().cmp(laws[b].species()));
    // This existing source-coefficient bound establishes increasing sensible
    // U on the entire bracket; it does not replace the actual caloric integrals.
    let mut gas_controls = controls.exchange.gas.gas;
    gas_controls.minimum_temperature_k = lo;
    gas_controls.maximum_temperature_k = hi;
    next_gas
        .model()
        .maximum_isochoric_temperature_slope(gas_controls)?;
    let root_tolerance = 0.5 * (controls.energy_tolerance_j - exchange.energy_residual_j.abs());
    if !positive(root_tolerance) {
        return Err(ThermalDissolutionError::Energy {
            residual_j: exchange.energy_residual_j,
            tolerance_j: controls.energy_tolerance_j,
        });
    }
    let dissolved_masses: Vec<_> = order
        .iter()
        .map(|&i| {
            let species = laws[i]
                .catalog()
                .thermodynamics()
                .components()
                .iter()
                .find(|s| s.species() == laws[i].species())
                .expect("admitted Henry species");
            finite(next_liquids[i].dissolved_moles() * species.molar_mass().value())
        })
        .collect::<Result<_, _>>()?;
    let target = finite(external_heat_j - exchange.thermostat_heat_into_system_j)?;
    let evaluate = |delta: f64, trials: usize, cancel: &mut _| {
        let t = finite(t0 + delta)?;
        let gas_sensible = finite(
            next_gas.mass_kg()
                * next_gas
                    .model()
                    .thermodynamics()
                    .specific_internal_energy_change(Temperature::new(t0), Temperature::new(t))?
                    .value(),
        )?;
        let mut dissolved = Dd::ZERO;
        for (&i, &mass) in order.iter().zip(&dissolved_masses) {
            poll(cancel)?;
            dissolved = dissolved
                + Dd::from_f64(finite(
                    mass * laws[i]
                        .calorics()
                        .partial_enthalpy
                        .specific_enthalpy_change(t0, t)?,
                )?);
        }
        let dissolved = finite(dissolved.to_f64())?;
        let solvent_change =
            finite(solvent.mass_kg * solvent.calorics.specific_enthalpy_change(t0, t)?)?;
        let residual = total([gas_sensible, dissolved, solvent_change, -target])?;
        Ok((
            (t, gas_sensible, dissolved, solvent_change, trials),
            residual,
        ))
    };
    let (t1, gas_sensible, dissolved_sensible, solvent_change, trials) = feedback::solve_scalar(
        [lo - t0, hi - t0],
        root_tolerance,
        controls.maximum_iterations,
        evaluate,
        &mut cancelled,
    )
    .map_err(|e| match e {
        feedback::ScalarError::Evaluation(e) => e,
        feedback::ScalarError::Invalid(s) => invalid(s),
        feedback::ScalarError::Cancelled => ThermalDissolutionError::Cancelled,
        feedback::ScalarError::OutsideBracket { lower, upper } => {
            ThermalDissolutionError::OutsideBracket {
                lower_residual_j: lower,
                upper_residual_j: upper,
            }
        }
        feedback::ScalarError::NotConverged {
            iterations,
            residual,
        } => ThermalDissolutionError::NotConverged {
            iterations,
            residual_j: residual,
        },
    })?;
    if (t1 - t0).abs() > controls.maximum_temperature_increment_k {
        return Err(invalid(
            "thermal absorption exceeds its declared temperature step",
        ));
    }
    let (accepted_gas, actual_gas_sensible) = next_gas.at_temperature(t1)?;
    let mut next_laws = laws.to_vec();
    let mut samples = Vec::with_capacity(laws.len());
    let mut concentration = Dd::ZERO;
    let mut interface_excess = Dd::ZERO;
    let mut active_pressure = Dd::ZERO;
    let mut gas_partial_sum = Dd::ZERO;
    let mut gas_fraction_limit: Option<f64> = None;
    for &i in &order {
        poll(&mut cancelled)?;
        let (law, sample) =
            laws[i].at_conditions(&accepted_gas, next_liquids[i], films[i], controls.exchange)?;
        concentration = concentration
            + Dd::from_f64(
                sample
                    .concentration_mol_m3
                    .max(sample.interface_concentration_mol_m3),
            );
        gas_partial_sum = gas_partial_sum
            + Dd::from_f64(
                sample
                    .partial_pressure_pa
                    .max(sample.interface_partial_pressure_pa),
            );
        if let Some(source) = law.gas_film_source() {
            gas_fraction_limit = Some(
                gas_fraction_limit.map_or(source.maximum_mole_fraction, |maximum| {
                    maximum.min(source.maximum_mole_fraction)
                }),
            );
        }
        if films[i].area_m2 > 0.0
            && law.diffusivity_m2_s() > 0.0
            && law.gas_film_diffusivity_m2_s().is_none_or(|d| d > 0.0)
            && films[i]
                .gas_side_transfer
                .is_none_or(|k| k.value_si() > 0.0)
        {
            interface_excess = interface_excess
                + Dd::from_f64(sample.interface_partial_pressure_pa - sample.partial_pressure_pa);
            active_pressure = active_pressure + Dd::from_f64(sample.partial_pressure_pa);
        }
        next_laws[i] = law;
        samples.push(sample);
    }
    let concentration = finite(concentration.to_f64())?;
    if concentration > controls.exchange.maximum_concentration_mol_m3 {
        return Err(DissolvedGasError::Concentration {
            actual_mol_m3: concentration,
            maximum_mol_m3: controls.exchange.maximum_concentration_mol_m3,
        }
        .into());
    }
    let pressure = accepted_gas.state().state().pressure;
    if let Some(maximum) = gas_fraction_limit {
        let fraction = finite(gas_partial_sum.to_f64() / pressure)?;
        if fraction > maximum {
            return Err(DissolvedGasError::GasFilmConcentration {
                actual_mole_fraction: fraction,
                maximum_mole_fraction: maximum,
            }
            .into());
        }
    }
    if finite(interface_excess.to_f64())?
        > pressure - active_pressure.to_f64() + 64.0 * f64::EPSILON * pressure
    {
        return Err(invalid(
            "thermal endpoint gas-film interfaces exceed common pressure",
        ));
    }
    let residual = total([
        exchange.thermostat_heat_into_system_j,
        exchange.energy_residual_j,
        actual_gas_sensible,
        dissolved_sensible,
        solvent_change,
        -external_heat_j,
    ])?;
    if residual.abs() > controls.energy_tolerance_j
        || (actual_gas_sensible - gas_sensible).abs() > root_tolerance
    {
        return Err(ThermalDissolutionError::Energy {
            residual_j: residual,
            tolerance_j: controls.energy_tolerance_j,
        });
    }
    let report = ThermalDissolutionStep {
        initial_temperature_k: t0,
        final_temperature_k: t1,
        solvent_initial,
        solvent_state: solvent.calorics.state_at_temperature(t1)?,
        isothermal_exchange: exchange,
        endpoint_sources: samples,
        gas_sensible_energy_change_j: actual_gas_sensible,
        dissolved_sensible_energy_change_j: dissolved_sensible,
        solvent_energy_change_j: solvent_change,
        external_heat_j,
        energy_residual_j: residual,
        temperature_trials: trials,
    };
    poll(&mut cancelled)?;
    *gas = accepted_gas;
    liquids.copy_from_slice(&next_liquids);
    laws.clone_from_slice(&next_laws);
    Ok(report)
}
