//! Molecular relaxation absorption and dispersion of an arbitrary gas mixture
//! (IRP M0.1-W2): one law for every composition.
//!
//! Each relaxing internal reservoir `j` of a species `s` contributes to the
//! complex mass-specific isochoric heat capacity
//!
//! ```text
//! c_v(ω) = c_v∞ + Σ_j C_j / (1 − iωτ_j)
//! C_j    = x_s R_s Σ_m g_m E(θ_m/T),   E(u) = u² e^u / (e^u − 1)²   (vibrational modes)
//! C_j    = x_s R_s n/2                                              (n classical rotational DOF)
//! c_v∞   = (c_p − R_s) − Σ_j C_j
//! c²(ω)  = R_s T (c_v(ω) + R_s) / c_v(ω)
//! k(ω)   = ω / c(ω) + i α_SK(ω)
//! ```
//!
//! with mole fraction `x_s`, specific gas constant `R_s`, Planck–Einstein mode
//! heat capacity (characteristic temperature `θ_m`, degeneracy `g_m`) and the
//! evaluated [`GasState`]'s local-equilibrium `c_p`. The sign convention is
//! fs-duct's `e^{−iωt}` with travelling waves `e^{i(kx − ωt)}`, so `Im k > 0`
//! decays. Strengths are never linearized: absorption and dispersion come from
//! one complex expression. Rotational relaxation is a reservoir like any other,
//! so no separate bulk viscosity is added. `α_SK` is
//! [`GasState::stokes_kirchhoff_absorption`] (shear viscosity and conduction),
//! counted exactly once. With no reservoirs, `Im k` equals it bitwise.
//!
//! Relaxation times are data ([`RelaxationRateLaw`]):
//! - collision partners, `1/(pτ_j) = Σ_k x_k / (pτ_jk)(T)`, with the
//!   Landau–Teller form `ln(pτ_jk) = a + b T^{-1/3}` per partner;
//! - a partner-independent power law `pτ = pτ_ref (T/T_ref)^n`;
//! - ISO 9613-1's empirical oxygen/nitrogen relaxation frequencies, as one
//!   declared variant that is valid only for terrestrial air
//!   (`τ = 1/(2π f_r)`).
//!
//! Refusals rather than defaults:
//! - a species present in the gas with no relaxation entry (non-relaxing species
//!   such as Ar declare zero reservoirs as data);
//! - a collision partner present in the gas with no rate;
//! - a state outside a rate's temperature validity;
//! - the ISO variant applied to anything but terrestrial air;
//! - more than [`MAX_RELAXATION_PROCESSES`] reservoirs;
//! - declared relaxing heat capacity exceeding the state's internal heat
//!   capacity above translation.
//!
//! Determinism: `fs_math::det` transcendentals, `fs_math::c64` complex
//! arithmetic, and fixed-order sums (components sorted by species id, each
//! species' reservoirs in declared order).

use fs_math::c64::C64;
use fs_math::det;
use fs_thermochem::SpeciesId;

use crate::MaterialError;
use crate::gas::GasState;
use crate::iso9613::iso9613_relaxation_frequencies_hz;
use crate::thermochemical_gas::ResolvedThermochemicalGasState;

/// Bounded native declarations of the same species/reservoir laws.
pub mod input;
/// Exact positive time-domain compression realization of the same reservoirs.
pub mod bulk;

/// The ISO terrestrial-air admission policy, shared with
/// [`crate::iso9613::TerrestrialAir`].
pub use crate::iso9613::{ISO9613_DRY_AIR_RELATIVE_TOLERANCE, ISO9613_MAX_DRY_CO2};

/// Most relaxing reservoirs one data set may declare across all species.
pub const MAX_RELAXATION_PROCESSES: usize = 32;
/// Most vibrational modes one reservoir may lump together.
pub const MAX_MODES_PER_PROCESS: usize = 16;
/// Most collision partners one reservoir's rate law may name.
pub const MAX_COLLISION_PARTNERS: usize = 16;
/// Most species one data set or medium may carry.
pub const MAX_RELAXING_SPECIES: usize = 128;

const COMPOSITION_SUM_TOLERANCE: f64 = 1.0e-9;

fn parameters(what: impl Into<String>) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}

fn state_error(what: impl Into<String>) -> MaterialError {
    MaterialError::State { what: what.into() }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

/// One vibrational mode: characteristic temperature `θ = h c ν̃ / k_B` [K] and
/// degeneracy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VibrationalMode {
    /// Characteristic temperature [K].
    pub theta_k: f64,
    /// Degeneracy (number of modes sharing `θ`).
    pub degeneracy: u8,
}

/// The internal degrees of freedom that relax together with one time.
#[derive(Debug, Clone, PartialEq)]
pub enum InternalReservoir {
    /// Vibrational modes lumped into one relaxation (for example the
    /// Fermi-coupled modes of CO2, or all modes of SF6), as the source does.
    Vibrational(Vec<VibrationalMode>),
    /// Classical rotation with `degrees_of_freedom` (2 linear, 3 nonlinear),
    /// valid far above the rotational characteristic temperature.
    Rotational {
        /// Rotational degrees of freedom.
        degrees_of_freedom: u8,
    },
}

/// Landau–Teller relaxation of one reservoir by one collision partner:
/// `ln(pτ / (Pa s)) = ln_p_tau_pa_s + landau_teller_slope_k13 · T^{-1/3}`.
#[derive(Debug, Clone, PartialEq)]
pub struct PartnerRate {
    /// Collision partner species.
    pub partner: SpeciesId,
    /// Intercept `a` of `ln(pτ / (Pa s))`.
    pub ln_p_tau_pa_s: f64,
    /// Slope `b` [K^{1/3}].
    pub landau_teller_slope_k13: f64,
    /// Lowest valid temperature [K].
    pub t_min_k: f64,
    /// Highest valid temperature [K].
    pub t_max_k: f64,
}

/// The two ISO 9613-1 relaxing species.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsoRelaxingSpecies {
    /// O2 vibration, `f_rO`.
    Oxygen,
    /// N2 vibration, `f_rN`.
    Nitrogen,
}

/// Data law for a reservoir's relaxation time.
#[derive(Debug, Clone, PartialEq)]
pub enum RelaxationRateLaw {
    /// Linear mole-fraction mixture rule over named collision partners. Every
    /// species present in the gas must have a rate here.
    CollisionPartners(Vec<PartnerRate>),
    /// `pτ = p_tau_ref_pa_s (T/t_ref_k)^exponent`, declared independent of the
    /// collision partner.
    PartnerIndependent {
        /// `pτ` at `t_ref_k` [Pa s].
        p_tau_ref_pa_s: f64,
        /// Reference temperature [K].
        t_ref_k: f64,
        /// Temperature exponent.
        exponent: f64,
        /// Lowest valid temperature [K].
        t_min_k: f64,
        /// Highest valid temperature [K].
        t_max_k: f64,
    },
    /// ISO 9613-1 empirical relaxation frequency for terrestrial air only.
    Iso9613TerrestrialAir(IsoRelaxingSpecies),
}

/// One relaxing reservoir and its rate law.
#[derive(Debug, Clone, PartialEq)]
pub struct RelaxationProcess {
    /// Degrees of freedom that relax together.
    pub reservoir: InternalReservoir,
    /// Relaxation-time law.
    pub rate: RelaxationRateLaw,
}

/// All relaxing reservoirs of one species (empty for a species that has none,
/// such as a monatomic gas).
#[derive(Debug, Clone, PartialEq)]
pub struct SpeciesRelaxation {
    /// Species identity, shared with the thermochemical mixture.
    pub species: SpeciesId,
    /// Reservoirs in declared order.
    pub processes: Vec<RelaxationProcess>,
}

/// Validated relaxation data for a set of species.
#[derive(Debug, Clone, PartialEq)]
pub struct RelaxationAbsorption {
    species: Vec<SpeciesRelaxation>,
}

/// A gas state together with the composition the law needs.
#[derive(Debug, Clone, PartialEq)]
pub struct RelaxingMedium {
    state: GasState,
    components: Vec<(SpeciesId, f64)>,
}

/// One reservoir evaluated at the medium's state.
#[derive(Debug, Clone, PartialEq)]
pub struct EvaluatedProcess {
    /// Owning species.
    pub species: SpeciesId,
    /// Index into that species' declared reservoirs.
    pub index: usize,
    /// Relaxing heat capacity `C_j` [J/(kg K)].
    pub heat_capacity_j_kg_k: f64,
    /// Relaxation time `τ_j` [s].
    pub relaxation_time_s: f64,
}

/// Frequency-independent evaluation of the law at one state and composition.
#[derive(Debug, Clone, PartialEq)]
pub struct RelaxationState {
    gas: GasState,
    cv_frozen: f64,
    processes: Vec<EvaluatedProcess>,
}

/// The same molecular reservoirs viewed as a local acoustic constitutive law.
/// Phasors use `exp(-i omega t)`: heat capacities have nonnegative imaginary
/// parts and the passive bulk modulus has a nonpositive imaginary part.
/// These are molecular losses only. Shear viscosity, heat conduction and
/// wall boundary layers remain separate physical owners.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RelaxationCaloricResponse {
    /// Dynamic mass-specific isochoric heat capacity [J/(kg K)].
    pub cv_j_kg_k: C64,
    /// Dynamic mass-specific isobaric heat capacity, `cv + R_specific` [J/(kg K)].
    pub cp_j_kg_k: C64,
    /// Dynamic adiabatic bulk modulus, `p cp/cv` [Pa].
    pub bulk_modulus_pa: C64,
}

fn check_temperature_range(t_min_k: f64, t_max_k: f64, what: &str) -> Result<(), MaterialError> {
    if positive(t_min_k) && t_max_k.is_finite() && t_min_k < t_max_k {
        Ok(())
    } else {
        Err(parameters(format!(
            "{what} needs a finite temperature range 0 < t_min < t_max"
        )))
    }
}

fn validate_process(species: &SpeciesId, process: &RelaxationProcess) -> Result<(), MaterialError> {
    let name = species.as_str();
    match &process.reservoir {
        InternalReservoir::Vibrational(modes) => {
            if modes.is_empty() || modes.len() > MAX_MODES_PER_PROCESS {
                return Err(parameters(format!(
                    "{name} vibrational reservoir needs 1..={MAX_MODES_PER_PROCESS} modes"
                )));
            }
            if modes
                .iter()
                .any(|mode| !positive(mode.theta_k) || mode.degeneracy == 0)
            {
                return Err(parameters(format!(
                    "{name} vibrational modes need positive finite theta and degeneracy >= 1"
                )));
            }
        }
        InternalReservoir::Rotational { degrees_of_freedom } => {
            if !matches!(degrees_of_freedom, 2 | 3) {
                return Err(parameters(format!(
                    "{name} rotational reservoir needs 2 (linear) or 3 (nonlinear) degrees of freedom"
                )));
            }
        }
    }
    match &process.rate {
        RelaxationRateLaw::CollisionPartners(partners) => {
            if partners.is_empty() || partners.len() > MAX_COLLISION_PARTNERS {
                return Err(parameters(format!(
                    "{name} collision-partner law needs 1..={MAX_COLLISION_PARTNERS} partners"
                )));
            }
            for (i, rate) in partners.iter().enumerate() {
                if partners[..i]
                    .iter()
                    .any(|other| other.partner == rate.partner)
                {
                    return Err(parameters(format!(
                        "{name} names collision partner {} twice",
                        rate.partner.as_str()
                    )));
                }
                if !(rate.ln_p_tau_pa_s.is_finite() && rate.landau_teller_slope_k13.is_finite()) {
                    return Err(parameters(format!(
                        "{name}-{} Landau-Teller coefficients must be finite",
                        rate.partner.as_str()
                    )));
                }
                check_temperature_range(rate.t_min_k, rate.t_max_k, "collision-partner rate")?;
            }
        }
        RelaxationRateLaw::PartnerIndependent {
            p_tau_ref_pa_s,
            t_ref_k,
            exponent,
            t_min_k,
            t_max_k,
        } => {
            if !(positive(*p_tau_ref_pa_s) && positive(*t_ref_k) && exponent.is_finite()) {
                return Err(parameters(format!(
                    "{name} partner-independent law needs positive p_tau_ref, t_ref and a finite exponent"
                )));
            }
            check_temperature_range(*t_min_k, *t_max_k, "partner-independent rate")?;
        }
        RelaxationRateLaw::Iso9613TerrestrialAir(iso) => {
            let expected = match iso {
                IsoRelaxingSpecies::Oxygen => "O2",
                IsoRelaxingSpecies::Nitrogen => "N2",
            };
            if name != expected || !matches!(process.reservoir, InternalReservoir::Vibrational(_)) {
                return Err(parameters(format!(
                    "the ISO 9613 {iso:?} rate law belongs to the vibrational reservoir of {expected}, not {name}"
                )));
            }
        }
    }
    Ok(())
}

impl RelaxationAbsorption {
    /// Validate and canonicalize relaxation data (species sorted by id).
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] for duplicate species, more than
    /// [`MAX_RELAXING_SPECIES`] species, more than [`MAX_RELAXATION_PROCESSES`]
    /// reservoirs in total, or any malformed reservoir or rate law.
    pub fn new(mut species: Vec<SpeciesRelaxation>) -> Result<Self, MaterialError> {
        if species.len() > MAX_RELAXING_SPECIES {
            return Err(parameters(format!(
                "relaxation data exceeds MAX_RELAXING_SPECIES = {MAX_RELAXING_SPECIES}"
            )));
        }
        species.sort_by(|a, b| a.species.cmp(&b.species));
        if species
            .windows(2)
            .any(|pair| pair[0].species == pair[1].species)
        {
            return Err(parameters("relaxation data names a species twice"));
        }
        let total: usize = species.iter().map(|entry| entry.processes.len()).sum();
        if total > MAX_RELAXATION_PROCESSES {
            return Err(parameters(format!(
                "relaxation data declares {total} reservoirs, above MAX_RELAXATION_PROCESSES = {MAX_RELAXATION_PROCESSES}"
            )));
        }
        for entry in &species {
            for process in &entry.processes {
                validate_process(&entry.species, process)?;
            }
        }
        Ok(Self { species })
    }

    /// ISO 9613-1's terrestrial air as relaxation data. It comprises:
    /// - the O2 and N2 vibrations, with ISO's characteristic temperatures
    ///   (2239.1 K and 3352 K) and ISO's relaxation frequencies
    ///   ([`RelaxationRateLaw::Iso9613TerrestrialAir`]);
    /// - Ar, CO2 and H2O declared without reservoirs. That is ISO's own
    ///   assumption inside its CO2 ≤ 1e-3 admission; ISO carries no reservoir
    ///   of its own for either species.
    ///
    /// Rotation is not a reservoir here. Its bulk-viscosity loss, which ISO's
    /// classical constant includes, is left out.
    ///
    /// Evaluating it refuses any gas that is not terrestrial air, naming the
    /// species.
    ///
    /// # Errors
    /// None for the fixed data; the signature forwards [`Self::new`].
    pub fn iso9613_terrestrial_air() -> Result<Self, MaterialError> {
        let id = |name: &str| {
            SpeciesId::new(name).map_err(|error| parameters(format!("species id {name}: {error}")))
        };
        let vibration = |theta_k: f64, iso: IsoRelaxingSpecies| RelaxationProcess {
            reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                theta_k,
                degeneracy: 1,
            }]),
            rate: RelaxationRateLaw::Iso9613TerrestrialAir(iso),
        };
        Self::new(vec![
            SpeciesRelaxation {
                species: id("O2")?,
                processes: vec![vibration(2239.1, IsoRelaxingSpecies::Oxygen)],
            },
            SpeciesRelaxation {
                species: id("N2")?,
                processes: vec![vibration(3352.0, IsoRelaxingSpecies::Nitrogen)],
            },
            SpeciesRelaxation {
                species: id("Ar")?,
                processes: Vec::new(),
            },
            SpeciesRelaxation {
                species: id("CO2")?,
                processes: Vec::new(),
            },
            SpeciesRelaxation {
                species: id("H2O")?,
                processes: Vec::new(),
            },
        ])
    }

    /// The law for a dry-air-like gas off ISO 9613-1's fractions, evaluated
    /// with ISO's O2 and N2 relaxation frequencies at the medium's own T, p
    /// and water fraction HELD as a declared assumption. The gas may hold N2
    /// and O2 in any proportion (either may be absent), with Ar, H2O, and CO2
    /// at or below ISO's dry-basis ceiling. ISO's frequencies are fits for
    /// terrestrial air; the O2-O2, O2-N2 and water partner rates of any other
    /// proportion are not sourced. Strengths follow the medium's own mole
    /// fractions and heat capacities.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] naming any other species present, CO2
    /// above the ceiling, a gas without N2, O2 or Ar, a state outside ISO's
    /// frequency window, and each refusal of [`Self::evaluate`].
    pub fn evaluate_iso9613_rates_held(medium: &RelaxingMedium) -> Result<RelaxationState, MaterialError> {
        let (mut dry, mut co2, mut water) = (0.0, 0.0, 0.0);
        for (species, fraction) in medium.components() {
            match species.as_str() {
                "N2" | "O2" | "Ar" => dry += fraction,
                "CO2" => {
                    dry += fraction;
                    co2 = *fraction;
                }
                "H2O" => water = *fraction,
                _ if *fraction == 0.0 => {}
                other => {
                    return Err(parameters(format!(
                        "ISO 9613 rates cannot be held for a gas containing {other}"
                    )));
                }
            }
        }
        if dry <= 0.0 {
            return Err(parameters("ISO 9613 rates cannot be held for a gas without dry air"));
        }
        if co2 / dry > ISO9613_MAX_DRY_CO2 {
            return Err(parameters(format!(
                "ISO 9613 rates cannot be held for CO2 at dry-basis {} above {ISO9613_MAX_DRY_CO2}",
                co2 / dry
            )));
        }
        let gas = medium.state();
        let (fr_o, fr_n) =
            iso9613_relaxation_frequencies_hz(gas.temperature, gas.pressure, 100.0 * water)?;
        let id = |name: &str| {
            SpeciesId::new(name).map_err(|error| parameters(format!("species id {name}: {error}")))
        };
        let vibration = |theta_k: f64, fr_hz: f64| RelaxationProcess {
            reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                theta_k,
                degeneracy: 1,
            }]),
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: gas.pressure / (core::f64::consts::TAU * fr_hz),
                t_ref_k: gas.temperature,
                exponent: 0.0,
                t_min_k: gas.temperature - 1.0,
                t_max_k: gas.temperature + 1.0,
            },
        };
        Self::new(vec![
            SpeciesRelaxation { species: id("O2")?, processes: vec![vibration(2239.1, fr_o)] },
            SpeciesRelaxation { species: id("N2")?, processes: vec![vibration(3352.0, fr_n)] },
            SpeciesRelaxation { species: id("Ar")?, processes: Vec::new() },
            SpeciesRelaxation { species: id("CO2")?, processes: Vec::new() },
            SpeciesRelaxation { species: id("H2O")?, processes: Vec::new() },
        ])?
        .evaluate(medium)
    }

    /// Canonical species data.
    #[must_use]
    pub fn species(&self) -> &[SpeciesRelaxation] {
        &self.species
    }

    /// Evaluate every reservoir's heat capacity and relaxation time at the
    /// medium's temperature, pressure and composition.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] for a present species without data, a
    /// present collision partner without a rate, a temperature outside a rate's
    /// validity, the ISO variant on a gas that is not terrestrial air, or relaxing
    /// heat capacity above the state's internal heat capacity beyond
    /// translation. [`MaterialError::State`] for unrepresentable results.
    pub fn evaluate(&self, medium: &RelaxingMedium) -> Result<RelaxationState, MaterialError> {
        let gas = medium.state;
        let t = gas.temperature;
        let p = gas.pressure;
        let r_s = gas.specific_gas_constant;
        let cv_eq = gas.specific_heat_cp - r_s;
        if !positive(cv_eq) {
            return Err(state_error(
                "gas state has no positive isochoric heat capacity",
            ));
        }
        let mut processes = Vec::new();
        let mut relaxing = 0.0;
        for (species, fraction) in &medium.components {
            if *fraction == 0.0 {
                continue;
            }
            let data = self
                .species
                .binary_search_by(|entry| entry.species.cmp(species))
                .map(|index| &self.species[index])
                .map_err(|_| {
                    parameters(format!(
                        "{} is present (x = {fraction}) with no relaxation data; a non-relaxing species declares zero reservoirs",
                        species.as_str()
                    ))
                })?;
            for (index, process) in data.processes.iter().enumerate() {
                let per_mole = match &process.reservoir {
                    InternalReservoir::Vibrational(modes) => {
                        let mut sum = 0.0;
                        for mode in modes {
                            sum += f64::from(mode.degeneracy) * planck_einstein(mode.theta_k / t);
                        }
                        sum
                    }
                    InternalReservoir::Rotational { degrees_of_freedom } => {
                        0.5 * f64::from(*degrees_of_freedom)
                    }
                };
                let heat_capacity = fraction * r_s * per_mole;
                let tau = relaxation_time(species, &process.rate, medium, t, p)?;
                if !(heat_capacity.is_finite() && heat_capacity >= 0.0 && positive(tau)) {
                    return Err(state_error(format!(
                        "{} reservoir {index} heat capacity or relaxation time is unrepresentable",
                        species.as_str()
                    )));
                }
                relaxing += heat_capacity;
                processes.push(EvaluatedProcess {
                    species: species.clone(),
                    index,
                    heat_capacity_j_kg_k: heat_capacity,
                    relaxation_time_s: tau,
                });
            }
        }
        let cv_frozen = cv_eq - relaxing;
        if !(cv_frozen.is_finite() && cv_frozen >= 1.5 * r_s * (1.0 - 1.0e-12)) {
            return Err(parameters(format!(
                "declared relaxing heat capacity {relaxing} J/(kg K) leaves a frozen c_v {cv_frozen} below translation (1.5 R_s = {})",
                1.5 * r_s
            )));
        }
        Ok(RelaxationState {
            gas,
            cv_frozen,
            processes,
        })
    }
}

/// Planck–Einstein heat capacity of one mode per `R`: `u² e^{−u} / (1 − e^{−u})²`.
fn planck_einstein(u: f64) -> f64 {
    let em = det::exp(-u);
    let one_minus = det::expm1(-u);
    u * u * em / (one_minus * one_minus)
}

fn relaxation_time(
    species: &SpeciesId,
    rate: &RelaxationRateLaw,
    medium: &RelaxingMedium,
    t: f64,
    p: f64,
) -> Result<f64, MaterialError> {
    let name = species.as_str();
    let in_range = |t_min: f64, t_max: f64, what: &str| {
        if (t_min..=t_max).contains(&t) {
            Ok(())
        } else {
            Err(parameters(format!(
                "{name} {what} is valid for T in [{t_min}, {t_max}] K, not {t} K"
            )))
        }
    };
    match rate {
        RelaxationRateLaw::CollisionPartners(partners) => {
            let t_cbrt_inv = det::exp(-det::ln(t) / 3.0);
            let mut inverse_p_tau = 0.0;
            for (partner, fraction) in &medium.components {
                if *fraction == 0.0 {
                    continue;
                }
                let rate = partners.iter().find(|rate| &rate.partner == partner).ok_or_else(|| {
                    parameters(format!(
                        "{name} relaxation has no rate for collision partner {} (x = {fraction})",
                        partner.as_str()
                    ))
                })?;
                in_range(rate.t_min_k, rate.t_max_k, "collision-partner rate")?;
                let p_tau =
                    det::exp(rate.ln_p_tau_pa_s + rate.landau_teller_slope_k13 * t_cbrt_inv);
                inverse_p_tau += fraction / p_tau;
            }
            Ok(1.0 / (p * inverse_p_tau))
        }
        RelaxationRateLaw::PartnerIndependent {
            p_tau_ref_pa_s,
            t_ref_k,
            exponent,
            t_min_k,
            t_max_k,
        } => {
            in_range(*t_min_k, *t_max_k, "partner-independent rate")?;
            Ok(p_tau_ref_pa_s * det::exp(exponent * det::ln(t / t_ref_k)) / p)
        }
        RelaxationRateLaw::Iso9613TerrestrialAir(iso) => {
            let h = iso_terrestrial_air_water_percent(medium)?;
            let (fr_o, fr_n) = iso9613_relaxation_frequencies_hz(t, p, h)?;
            let fr = match iso {
                IsoRelaxingSpecies::Oxygen => fr_o,
                IsoRelaxingSpecies::Nitrogen => fr_n,
            };
            Ok(1.0 / (core::f64::consts::TAU * fr))
        }
    }
}

/// Molar water percent of a terrestrial-air medium, or the refusal naming the
/// species, by the same admission as [`crate::iso9613::TerrestrialAir`].
fn iso_terrestrial_air_water_percent(medium: &RelaxingMedium) -> Result<f64, MaterialError> {
    let (_, h2o) = crate::iso9613::admit_composition(&medium.components)
        .map_err(|refusal| parameters(format!("the ISO 9613 rate law: {refusal}")))?;
    Ok(100.0 * h2o)
}

impl RelaxingMedium {
    /// Pair a coherent gas state with its explicit mole-fraction composition.
    /// Components are sorted by species id. Zero fractions are retained and
    /// ignored by the law.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] for an incoherent state, an empty or
    /// oversized composition, duplicate species, fractions outside [0, 1], or a
    /// sum differing from one by more than 1e-9.
    pub fn new(
        state: GasState,
        mut components: Vec<(SpeciesId, f64)>,
    ) -> Result<Self, MaterialError> {
        state.validate_coherence()?;
        if components.is_empty() || components.len() > MAX_RELAXING_SPECIES {
            return Err(parameters(format!(
                "relaxing medium needs 1..={MAX_RELAXING_SPECIES} components"
            )));
        }
        components.sort_by(|a, b| a.0.cmp(&b.0));
        if components.windows(2).any(|pair| pair[0].0 == pair[1].0) {
            return Err(parameters("relaxing medium names a species twice"));
        }
        if components
            .iter()
            .any(|(_, x)| !(x.is_finite() && (0.0..=1.0).contains(x)))
        {
            return Err(parameters(
                "relaxing medium mole fractions must be finite and in [0, 1]",
            ));
        }
        let mut sum = 0.0;
        for (_, x) in &components {
            sum += x;
        }
        if (sum - 1.0).abs() > COMPOSITION_SUM_TOLERANCE {
            return Err(parameters(format!(
                "relaxing medium mole fractions sum to {sum}, not 1"
            )));
        }
        Ok(Self { state, components })
    }

    /// The thermochemical state's own species and mole fractions.
    ///
    /// # Errors
    /// As [`Self::new`].
    pub fn from_thermochemical(
        resolved: &ResolvedThermochemicalGasState,
    ) -> Result<Self, MaterialError> {
        let receipt = resolved.thermodynamics().receipt();
        let components = receipt
            .species()
            .iter()
            .cloned()
            .zip(
                receipt
                    .mole_fraction_bits()
                    .iter()
                    .map(|bits| f64::from_bits(*bits)),
            )
            .collect();
        Self::new(*resolved.state(), components)
    }

    /// The gas state.
    #[must_use]
    pub const fn state(&self) -> &GasState {
        &self.state
    }

    /// Mole fractions sorted by species id.
    #[must_use]
    pub fn components(&self) -> &[(SpeciesId, f64)] {
        &self.components
    }
}

/// `1/(1 − i x)` without overflow for large `x`.
fn debye(x: f64) -> (f64, f64) {
    if x <= 1.0 {
        let d = 1.0 + x * x;
        (1.0 / d, x / d)
    } else {
        let y = 1.0 / x;
        let d = 1.0 + y * y;
        (y * y / d, y / d)
    }
}

impl RelaxationState {
    /// The evaluated gas state.
    #[must_use]
    pub const fn gas(&self) -> &GasState {
        &self.gas
    }

    /// Evaluated reservoirs in canonical order.
    #[must_use]
    pub fn processes(&self) -> &[EvaluatedProcess] {
        &self.processes
    }

    /// Frozen (all reservoirs unrelaxed) isochoric heat capacity [J/(kg K)].
    #[must_use]
    pub const fn frozen_cv_j_kg_k(&self) -> f64 {
        self.cv_frozen
    }

    /// High-frequency sound speed with every reservoir frozen [m/s].
    #[must_use]
    pub fn frozen_sound_speed_m_s(&self) -> f64 {
        let r_s = self.gas.specific_gas_constant;
        det::sqrt(r_s * self.gas.temperature * (self.cv_frozen + r_s) / self.cv_frozen)
    }

    /// Evaluate the original reservoir sum for a duct, cavity or pressure port.
    ///
    /// `cv = cv_frozen + sum C_j/(1 - i omega tau_j)`, `cp = cv + R`,
    /// and `K = p cp/cv`. At zero frequency the reservoirs are equilibrated;
    /// at high frequency they freeze. The density remains the original
    /// equilibrium density: the constitutive response is not a second EOS.
    /// In particular, a propagation attenuation coefficient is not converted
    /// into fictitious compressional loss for a uniformly breathing cavity.
    /// Work against this modulus dissipates `-omega Im(K) |strain|^2 / 2`
    /// per unit volume for peak-amplitude phasors.
    ///
    /// Evaluation is allocation-free and bounded by the admitted reservoir
    /// count. It retains the same fixed summation order as wave propagation.
    ///
    /// # Errors
    /// Refuses negative/nonfinite frequency or an unrepresentable response.
    pub fn caloric_response(&self, omega: f64) -> Result<RelaxationCaloricResponse, MaterialError> {
        if !omega.is_finite() || omega < 0.0 {
            return Err(parameters(
                "caloric response requires finite nonnegative omega",
            ));
        }
        let r_s = self.gas.specific_gas_constant;
        let mut cv = C64::from_re(self.cv_frozen);
        for process in &self.processes {
            let (re, im) = debye(omega * process.relaxation_time_s);
            let c = process.heat_capacity_j_kg_k;
            cv = cv + C64::new(c * re, c * im);
        }
        let cp = cv + C64::from_re(r_s);
        // 1 + R/cv avoids cancellation of cv.im in (cv + R)/cv when
        // a reservoir's loss is tiny, retaining its passive sign.
        let bulk = (C64::ONE + cv.recip().scale(r_s)).scale(self.gas.pressure);
        if !(positive(cv.re)
            && cv.im.is_finite()
            && cv.im >= 0.0
            && positive(cp.re)
            && cp.im.is_finite()
            && positive(bulk.re)
            && bulk.im.is_finite()
            && bulk.im <= 0.0)
        {
            return Err(state_error("molecular caloric response is unrepresentable"));
        }
        Ok(RelaxationCaloricResponse {
            cv_j_kg_k: cv,
            cp_j_kg_k: cp,
            bulk_modulus_pa: bulk,
        })
    }

    /// Complex wavenumber of the relaxing reservoirs alone, `ω/c(ω)` [1/m],
    /// without the viscothermal term.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] for nonpositive or non-finite `omega`;
    /// [`MaterialError::State`] for an unrepresentable result.
    pub fn relaxing_wavenumber(&self, omega: f64) -> Result<C64, MaterialError> {
        if !positive(omega) {
            return Err(parameters(format!(
                "relaxation omega {omega} must be positive and finite"
            )));
        }
        let response = self.caloric_response(omega)?;
        let c2 = response.bulk_modulus_pa.scale(1.0 / self.gas.density);
        let k = c2.sqrt().recip().scale(omega);
        if !(k.re.is_finite() && k.re > 0.0 && k.im.is_finite()) {
            return Err(state_error("relaxing wavenumber is unrepresentable"));
        }
        Ok(k)
    }

    /// Complete complex wavenumber: relaxation plus Stokes–Kirchhoff, once.
    ///
    /// # Errors
    /// As [`Self::relaxing_wavenumber`].
    pub fn complex_wavenumber(&self, omega: f64) -> Result<C64, MaterialError> {
        let k = self.relaxing_wavenumber(omega)?;
        let im = k.im + self.gas.stokes_kirchhoff_absorption(omega);
        if !(im.is_finite() && im >= 0.0) {
            return Err(state_error(
                "relaxation absorption is not finite and nonnegative",
            ));
        }
        Ok(C64::new(k.re, im))
    }

    /// Causal complete wavenumber for propagating a history:
    /// `k = (w / c(w)) (1 - i w delta)^{-1/2}` with
    /// `delta = (4 mu / 3 + (gamma - 1) kappa / c_p) / (rho c^2)` from the
    /// equilibrium state. The viscothermal factor is the first-order Kirchhoff
    /// dispersion relation. Its imaginary part reproduces Stokes-Kirchhoff to
    /// relative `O(w delta)` (`delta` is about 1e-10 s in air), and unlike adding
    /// `i alpha_SK` it is analytic in the upper half plane. A history
    /// propagated with it therefore stays causal. The losses are still counted
    /// once.
    ///
    /// # Errors
    /// As [`Self::relaxing_wavenumber`].
    pub fn causal_wavenumber(&self, omega: f64) -> Result<C64, MaterialError> {
        let k = self.relaxing_wavenumber(omega)?;
        let gas = &self.gas;
        let delta = (4.0 / 3.0 * gas.dynamic_viscosity
            + (gas.gamma - 1.0) * gas.thermal_conductivity / gas.specific_heat_cp)
            / (gas.density * gas.sound_speed * gas.sound_speed);
        let k = k * C64::new(1.0, -omega * delta).sqrt().recip();
        if !(k.re.is_finite() && k.re > 0.0 && k.im.is_finite() && k.im >= 0.0) {
            return Err(state_error("causal wavenumber is unrepresentable"));
        }
        Ok(k)
    }

    /// Amplitude absorption coefficient `Im k` [Np/m].
    ///
    /// # Errors
    /// As [`Self::complex_wavenumber`].
    pub fn absorption_np_per_m(&self, omega: f64) -> Result<f64, MaterialError> {
        Ok(self.complex_wavenumber(omega)?.im)
    }

    /// Phase speed `ω / Re k` [m/s].
    ///
    /// # Errors
    /// As [`Self::relaxing_wavenumber`].
    pub fn phase_speed_m_s(&self, omega: f64) -> Result<f64, MaterialError> {
        Ok(omega / self.relaxing_wavenumber(omega)?.re)
    }
}
