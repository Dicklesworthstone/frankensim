//! Gas-related losses: cylinder motion and free-field observer absorption.
//!
//! A generic gas uses its own Stokes-Kirchhoff viscothermal loss. The
//! air-specific ISO 9613 law needs an explicit choice and a
//! [`TerrestrialAir`] admission of the gas's own composition. Temperature and
//! pressure alone do not establish that an arbitrary mixture is air.

/// Refinable causal gas inertia and viscous memory for circular transverse motion.
pub mod cylinder;

use fs_fft::{C64 as FftC64, Fft};
use fs_material::gas::GasState;
use fs_material::iso9613::{TerrestrialAir, iso9613_absorption};
use fs_material::relaxation::RelaxationState;
use fs_math::det;
use crate::acoustic_realize::AcousticRealizeError;

/// Constitutive choice for observer-path attenuation, never inferred from a
/// material name, gamma, molar mass, or an unrelated humidity coordinate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GasAbsorptionModel {
    /// Classical viscothermal absorption from this gas's rho, c, mu, k and cp.
    /// Molecular relaxation, bulk viscosity and chemical relaxation are not
    /// supplied by this model. No air-specific coefficients are substituted.
    StokesKirchhoff,
    /// ISO 9613-1 for a gas admitted as terrestrial air. The ISO term already
    /// includes classical loss; it must not be added to Stokes-Kirchhoff.
    /// Out-of-window requests refuse rather than switch laws, and so does a
    /// gas state the admission was not made for.
    Iso9613Air {
        /// Relative humidity as a fraction in [0, 1], not water mole fraction.
        relative_humidity: f64,
        /// The admission of the gas's own composition, bound to its state.
        air: TerrestrialAir,
    },
    /// The gas's own molecular relaxation: `fs_material::relaxation` evaluated
    /// from its composition and its species' relaxation data, together with
    /// Stokes-Kirchhoff, counted once. It is applied as a causal propagation
    /// ([`propagate_pressure_history_relaxing`]), so absorption and dispersion
    /// come from one wavenumber. The evaluated state is supplied beside this
    /// tag: `RealizeOptions::observer_relaxation`.
    MolecularRelaxation,
    /// A declared lossless observation (`observer_law none`): nothing is
    /// absorbed and the history is left untouched. It exists for comparisons
    /// against a lossless analytic oracle, and is never a default. Its receipt
    /// labels everything the gas would absorb as omitted.
    Lossless,
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

/// Mechanical resistance per unit length [N s/m²] of a slender circular
/// cylinder oscillating transversely in stationary gas.
///
/// Uses `R = 2 pi mu (1 + r sqrt(2 rho omega / mu))`, the small-amplitude
/// boundary-layer approximation in Desvages (2018), section 3.2.2, equations
/// 3.8 and 3.12: <https://hdl.handle.net/1842/31273>.
/// Here `mu` is dynamic viscosity [Pa s], not kinematic viscosity [m²/s].
/// For linear mass density `m_l`, modal damping is `zeta = R/(2 m_l omega)`.
///
/// This is an estimate for a long cylinder, including the model's constant
/// low-frequency continuation. It does not solve end effects, fluid added mass,
/// turbulence, confinement, rarefied gas, or large-amplitude motion. Applying it
/// to a nonlinear oscillator retains a linear, modal loss approximation.
///
/// # Errors
/// Refuses nonpositive/nonfinite radius [m], frequency [rad/s], incoherent gas
/// states, and unrepresentable derived resistance.
pub fn oscillating_cylinder_air_resistance_per_length(
    radius_m: f64,
    omega_rad_s: f64,
    gas: &GasState,
) -> Result<f64, AcousticRealizeError> {
    gas.validate_coherence().map_err(|_| invalid("cylinder resistance needs a coherent gas state"))?;
    if [radius_m, omega_rad_s, gas.dynamic_viscosity, gas.density]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(invalid("cylinder air resistance needs positive finite radius, angular frequency, viscosity and gas density"));
    }
    let resistance = core::f64::consts::TAU
        * gas.dynamic_viscosity
        * (1.0 + radius_m * det::sqrt(2.0 * gas.density * omega_rad_s / gas.dynamic_viscosity));
    if !resistance.is_finite() || resistance <= 0.0 {
        return Err(invalid("cylinder air resistance is unrepresentable"));
    }
    Ok(resistance)
}

/// Apply an explicitly chosen absorption law to the entire pressure history.
///
/// The transfer is real and even: this is the existing zero-phase map on the
/// padded DFT grid, not a causal propagation or dispersion solver. A zero range
/// and histories shorter than four samples are no-ops after input admission.
/// Source/material receipts remain with the caller's evaluated gas model.
///
/// # Errors
/// Invalid samples, time/range, incoherent gas, a gas the ISO admission was
/// not made for ([`AcousticRealizeError::NonTerrestrialAir`]), unavailable ISO
/// conditions, or unrepresentable absorption/FFT output refuse WITHOUT
/// changing the history.
pub fn absorb_pressure_history_with_model(
    pressure_pa: &mut [f64],
    dt: f64,
    range_m: f64,
    gas: &GasState,
    model: GasAbsorptionModel,
) -> Result<(), AcousticRealizeError> {
    gas.validate_coherence().map_err(|_| invalid("observer absorption needs a coherent gas state"))?;
    if model == GasAbsorptionModel::MolecularRelaxation {
        return Err(invalid(
            "molecular relaxation is applied by propagate_pressure_history_relaxing with its evaluated state",
        ));
    }
    if model == GasAbsorptionModel::Lossless {
        if !dt.is_finite() || dt <= 0.0 || !range_m.is_finite() || range_m < 0.0
            || pressure_pa.iter().any(|sample| !sample.is_finite())
        {
            return Err(invalid("absorption needs finite pressure, positive dt and nonnegative range"));
        }
        return Ok(());
    }
    if let GasAbsorptionModel::Iso9613Air { relative_humidity, air } = model {
        // Admit the declared air and its conditions even for an empty or
        // zero-range history. A failed air law is not permission to use a
        // different law.
        air.check(gas).map_err(AcousticRealizeError::NonTerrestrialAir)?;
        iso9613_absorption(&air, gas, relative_humidity, 1.0)
            .map_err(|_| invalid("explicit ISO 9613 air conditions are outside the model domain"))?;
    }
    transform_history(pressure_pa, dt, range_m, |omega| match model {
        GasAbsorptionModel::StokesKirchhoff => Ok(gas.stokes_kirchhoff_absorption(omega)),
        GasAbsorptionModel::Iso9613Air { relative_humidity, air } => {
            iso9613_absorption(&air, gas, relative_humidity, omega)
                .map_err(|_| invalid("explicit ISO 9613 absorption could not be evaluated"))
        }
        GasAbsorptionModel::MolecularRelaxation | GasAbsorptionModel::Lossless => {
            Err(invalid("molecular relaxation and lossless observation are not zero-phase absorptions"))
        }
    })
}

/// What an omitted-relaxation figure is. It is labelled, so that a convention
/// is not read as a bound.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OmissionBound {
    /// Standard air's worst ISO 9613-1 absorption at the band top, times the
    /// range [dB]: the convention of the orchestrator ruling of 2026-10-23 on
    /// pexfp.18, for a gas whose relaxation rates are unsourced. It is not a
    /// bound for a gas that relaxes more strongly than air, such as a CO2-rich
    /// mixture.
    StandardAirConvention {
        /// The figure [dB].
        db: f64,
    },
}

/// Molecular relaxation that an observer path did not apply, with a labelled
/// figure for its effect at the listener. A render may omit relaxation only
/// with this receipt, never silently (orchestrator ruling of 2026-10-23 on
/// pexfp.18).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RelaxationOmission {
    /// Listener distance the figure covers [m].
    pub range_m: f64,
    /// Highest frequency of the observed band [Hz].
    pub band_hz: f64,
    /// The omitted attenuation over the band, and what kind of figure it is.
    pub bound: OmissionBound,
}

/// What an observer path applied for gas absorption.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ObserverAbsorptionReceipt {
    /// The law applied to the pressure history.
    pub law: GasAbsorptionModel,
    /// Molecular relaxation that the law leaves out, with its labelled figure.
    /// `None` when the law includes it (ISO 9613-1 for terrestrial air).
    pub relaxation_omitted: Option<RelaxationOmission>,
}

/// The standard-air convention for omitted molecular relaxation: the worst
/// ISO 9613-1 absorption of standard air at `band_hz`, times `range_m`, in dB
/// ([`OmissionBound::StandardAirConvention`]).
///
/// "Worst" is the maximum over ISO's tabulated domain (-20..+50 C, 10..100 %
/// RH, 1 atm), sampled on a 1 K by 1 % grid. ISO absorption rises
/// monotonically with frequency, so the band's top covers the whole band.
/// ISO's figure includes classical loss, so the figure also covers the
/// viscothermal part.
///
/// This is the ruling's convention for mixtures whose relaxation rates are
/// unsourced. It is not a bound for such mixtures: a CO2-rich gas can relax
/// more strongly than air.
///
/// # Errors
/// A negative or non-finite range, or a band that is not positive and finite.
pub fn relaxation_omission(range_m: f64, band_hz: f64) -> Result<RelaxationOmission, AcousticRealizeError> {
    if !(range_m.is_finite() && range_m >= 0.0 && band_hz.is_finite() && band_hz > 0.0) {
        return Err(invalid("relaxation omission needs a finite nonnegative range and positive band"));
    }
    let omega = core::f64::consts::TAU * band_hz;
    let mut worst: f64 = 0.0;
    for kelvin in 0..=70 {
        for percent in 10..=100 {
            let db_per_km = fs_material::iso9613::iso9613_absorption_db_per_km(
                253.15 + f64::from(kelvin), 101_325.0, f64::from(percent) / 100.0, omega,
            ).map_err(|_| invalid("standard-air absorption could not be evaluated for the bound"))?;
            worst = worst.max(db_per_km);
        }
    }
    let db = worst * 1e-3 * range_m;
    if !db.is_finite() {
        return Err(invalid("relaxation omission figure is unrepresentable"));
    }
    Ok(RelaxationOmission { range_m, band_hz, bound: OmissionBound::StandardAirConvention { db } })
}

/// The absorption a render may leave unapplied at the top of its band [dB]:
/// about the level JND for a broadband change, here falling only on the top
/// of the band, where absorption is largest (at 8 kHz the worst over ISO
/// 9613-1's domain is 0.26 dB at 1 m). A render whose omitted figure exceeds
/// it refuses unless a law is applied or declared (rulings 6498, 6510).
pub const OBSERVER_OMISSION_TOLERANCE_DB: f64 = 1.0;

/// The highest frequency an omitted figure is taken at [Hz]: the edge of the
/// audible band. Absorption above it changes nothing a listener can hear.
pub const OBSERVER_BAND_TOP_CAP_HZ: f64 = 20_000.0;

const NEPER_TO_DB: f64 = 20.0 / core::f64::consts::LN_10;

/// The band top an omitted figure is taken at: the path's declared band, the
/// clock's Nyquist frequency and [`OBSERVER_BAND_TOP_CAP_HZ`], whichever is
/// lowest.
#[must_use]
pub fn observer_band_top_hz(declared_band_hz: Option<f64>, sample_rate_hz: u32) -> f64 {
    let nyquist = 0.5 * f64::from(sample_rate_hz);
    declared_band_hz.map_or(nyquist, |band| band.min(nyquist)).min(OBSERVER_BAND_TOP_CAP_HZ)
}

/// The whole absorption a listener gas would apply over a range, by the law
/// its composition admits, with that gas's classical (Stokes-Kirchhoff) part.
#[derive(Debug, Clone, PartialEq)]
pub enum ObserverAbsorptionFigure {
    /// ISO 9613-1, for a gas admitted as terrestrial air.
    Iso9613 {
        /// Relative humidity as a fraction.
        relative_humidity: f64,
        /// The figure [dB].
        db: f64,
    },
    /// The relaxation law with ISO 9613-1's O2 and N2 relaxation frequencies
    /// HELD, a declared assumption, for N2 and O2 off ISO's fractions
    /// (`RelaxationAbsorption::evaluate_iso9613_rates_held`).
    HeldIsoRates {
        /// O2 mole fraction.
        o2: f64,
        /// N2 mole fraction.
        n2: f64,
        /// H2O mole fraction.
        water: f64,
        /// The figure [dB].
        db: f64,
    },
    /// A gas of monatomic species only: it has no internal reservoir, so
    /// Stokes-Kirchhoff is its complete law.
    ClassicalComplete {
        /// The figure [dB].
        db: f64,
    },
    /// A polyatomic species whose relaxation rates are not sourced: there is
    /// no figure, only the classical part as a floor.
    Unbounded {
        /// The first such species, by id.
        species: String,
    },
    /// The listener gas has no record naming its composition (a resolved gas
    /// field cell): there is no figure, only the classical part as a floor.
    Unrecorded,
}

/// The absorption between source and listener at a path's range and band
/// top, for its receipt and its admission.
#[derive(Debug, Clone, PartialEq)]
pub struct ObserverOmission {
    /// Listener range [m].
    pub range_m: f64,
    /// [`observer_band_top_hz`] of the path [Hz].
    pub band_top_hz: f64,
    /// The gas's whole absorption over the range, by the law it admits.
    pub figure: ObserverAbsorptionFigure,
    /// The gas's Stokes-Kirchhoff absorption over the range [dB].
    pub classical_db: f64,
}

impl ObserverOmission {
    /// The receipt for a listener gas without a record naming its composition.
    #[must_use]
    pub fn unrecorded(gas: &GasState, range_m: f64, band_top_hz: f64) -> Self {
        let classical = gas.stokes_kirchhoff_absorption(core::f64::consts::TAU * band_top_hz);
        Self {
            range_m,
            band_top_hz,
            figure: ObserverAbsorptionFigure::Unrecorded,
            classical_db: classical * range_m * NEPER_TO_DB,
        }
    }

    /// What a path leaves out [dB]: the whole figure when it applies nothing,
    /// the part beyond Stokes-Kirchhoff when it applies Stokes-Kirchhoff.
    /// `None` when unbounded.
    #[must_use]
    pub fn omitted_db(&self, stokes_kirchhoff_applied: bool) -> Option<f64> {
        let whole = match self.figure {
            ObserverAbsorptionFigure::Iso9613 { db, .. }
            | ObserverAbsorptionFigure::HeldIsoRates { db, .. }
            | ObserverAbsorptionFigure::ClassicalComplete { db } => db,
            ObserverAbsorptionFigure::Unbounded { .. } | ObserverAbsorptionFigure::Unrecorded => return None,
        };
        Some(if stokes_kirchhoff_applied { (whole - self.classical_db).max(0.0) } else { whole })
    }

    /// Admit a path that applies no law, or Stokes-Kirchhoff only. A declared
    /// law (`observer_law ...`, including `none`) is an explicit acceptance of
    /// what it omits and is always admitted, with this receipt. Otherwise an
    /// unbounded figure, or an omitted figure above
    /// [`OBSERVER_OMISSION_TOLERANCE_DB`], refuses.
    ///
    /// # Errors
    /// [`AcousticRealizeError::ObserverOmission`] carrying this receipt.
    pub fn admit(&self, stokes_kirchhoff_applied: bool, declared: bool) -> Result<(), AcousticRealizeError> {
        match self.omitted_db(stokes_kirchhoff_applied) {
            _ if declared => Ok(()),
            Some(db) if db <= OBSERVER_OMISSION_TOLERANCE_DB => Ok(()),
            _ => Err(AcousticRealizeError::ObserverOmission(Box::new(self.clone()))),
        }
    }

    /// The receipt as one JSON object for a render sidecar.
    #[must_use]
    pub fn to_json(&self) -> String {
        let kind = match &self.figure {
            ObserverAbsorptionFigure::Iso9613 { relative_humidity, db } => format!(
                "\"kind\":\"iso9613\",\"relative_humidity\":{relative_humidity:e},\"whole_db\":{db:e}"
            ),
            ObserverAbsorptionFigure::HeldIsoRates { o2, n2, water, db } => format!(
                "\"kind\":\"held-iso-rates\",\"o2\":{o2:e},\"n2\":{n2:e},\"water\":{water:e},\"whole_db\":{db:e},\
                 \"held\":\"ISO 9613-1 O2 and N2 relaxation frequencies at this T, p and water fraction (declared assumption)\""
            ),
            ObserverAbsorptionFigure::ClassicalComplete { db } => {
                format!("\"kind\":\"classical-complete\",\"whole_db\":{db:e}")
            }
            ObserverAbsorptionFigure::Unbounded { species } => {
                format!("\"kind\":\"unbounded\",\"species\":\"{species}\",\"whole_db\":null")
            }
            ObserverAbsorptionFigure::Unrecorded => "\"kind\":\"unrecorded\",\"whole_db\":null".to_owned(),
        };
        format!(
            "{{\"range_m\":{:e},\"band_top_hz\":{:e},{kind},\"classical_db\":{:e},\"tolerance_db\":{OBSERVER_OMISSION_TOLERANCE_DB:e}}}",
            self.range_m, self.band_top_hz, self.classical_db
        )
    }
}

impl core::fmt::Display for ObserverOmission {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match &self.figure {
            ObserverAbsorptionFigure::Unbounded { species } => write!(
                f,
                "the listener gas contains {species}, whose relaxation is not sourced: the absorption \
                 over {} m up to {} Hz is unbounded (classical floor {:.3} dB)",
                self.range_m, self.band_top_hz, self.classical_db
            ),
            ObserverAbsorptionFigure::Unrecorded => write!(
                f,
                "the listener gas has no record naming its composition: the absorption over {} m \
                 up to {} Hz is unbounded (classical floor {:.3} dB)",
                self.range_m, self.band_top_hz, self.classical_db
            ),
            _ => write!(
                f,
                "the listener gas absorbs {:.3} dB over {} m at {} Hz ({:.3} dB beyond Stokes-Kirchhoff), \
                 above the {OBSERVER_OMISSION_TOLERANCE_DB} dB tolerance",
                self.omitted_db(false).unwrap_or(f64::NAN),
                self.range_m,
                self.band_top_hz,
                self.omitted_db(true).unwrap_or(f64::NAN)
            ),
        }
    }
}

/// The authored relative humidity of a legacy (T, p, RH) record, recovered
/// from its state: `x_w = RH e_sat(T) / p` inverted, then accepted only if
/// [`GasState::try_new_moist_air`] reproduces the state bit for bit.
#[must_use]
pub fn legacy_relative_humidity(state: &GasState) -> Option<f64> {
    let (t, p) = (state.temperature, state.pressure);
    let guess = if state.water_mole_fraction == 0.0 {
        0.0
    } else {
        state.water_mole_fraction * p / fs_material::gas::saturation_pressure_water_pa(t).ok()?
    };
    let reproduces = |rh: f64| {
        (0.0..=1.0).contains(&rh) && GasState::try_new_moist_air(t, p, rh).ok() == Some(*state)
    };
    let (mut up, mut down) = (guess, guess);
    for _ in 0..=16 {
        if reproduces(up) {
            return Some(up);
        }
        if reproduces(down) {
            return Some(down);
        }
        up = up.next_up();
        down = down.next_down();
    }
    None
}

const MONATOMIC: [&str; 6] = ["He", "Ne", "Ar", "Kr", "Xe", "Rn"];
const HELD_RATE_SPECIES: [&str; 5] = ["N2", "O2", "Ar", "H2O", "CO2"];

/// The absorption a listener gas applies over `range_m` up to `band_top_hz`,
/// by the law its record admits: ISO 9613-1 for a legacy record or a
/// composition admitted as terrestrial air; Stokes-Kirchhoff, complete, for
/// monatomic species only; the relaxation law with ISO's rates held for N2
/// and O2 off ISO's fractions (with Ar, H2O and CO2 at or below ISO's
/// ceiling); otherwise unbounded, naming the species. The composition is read
/// from the record, never inferred from the state.
///
/// # Errors
/// A legacy record whose humidity cannot be recovered, a state outside a
/// law's window, or an invalid range or band.
pub fn omitted_observer_absorption(
    record: &fs_material::gas_species::input::GasRecord,
    range_m: f64,
    band_top_hz: f64,
) -> Result<ObserverOmission, AcousticRealizeError> {
    use fs_material::gas_species::input::GasRecord;
    use fs_material::relaxation::{RelaxationAbsorption, RelaxingMedium};
    if !(range_m.is_finite() && range_m >= 0.0 && band_top_hz.is_finite() && band_top_hz > 0.0) {
        return Err(invalid("an observer figure needs a finite nonnegative range and a positive band top"));
    }
    let gas = record.state();
    let omega = core::f64::consts::TAU * band_top_hz;
    let over_range = |alpha_np_per_m: f64| alpha_np_per_m * range_m * NEPER_TO_DB;
    let ambient = |error: fs_material::MaterialError| AcousticRealizeError::Ambient(error.to_string());
    let classical_db = over_range(gas.stokes_kirchhoff_absorption(omega));
    let iso = |relative_humidity: f64| {
        fs_material::iso9613::iso9613_absorption_neper_per_m(
            gas.temperature, gas.pressure, relative_humidity, omega,
        )
        .map(|alpha| ObserverAbsorptionFigure::Iso9613 { relative_humidity, db: over_range(alpha) })
        .map_err(ambient)
    };
    let medium = match record {
        GasRecord::Legacy(state) => {
            let rh = legacy_relative_humidity(state)
                .ok_or(invalid("a legacy observer gas record must reproduce its own humidity"))?;
            return Ok(ObserverOmission { range_m, band_top_hz, figure: iso(rh)?, classical_db });
        }
        GasRecord::Mixture(state, composition) => {
            let components = composition
                .iter()
                .map(|(species, fraction)| {
                    fs_thermochem::SpeciesId::new(species.symbol()).map(|id| (id, *fraction))
                })
                .collect::<Result<Vec<_>, _>>()
                .map_err(|error| AcousticRealizeError::Ambient(error.to_string()))?;
            RelaxingMedium::new(*state, components).map_err(ambient)?
        }
        GasRecord::Thermochemical(resolved) => {
            RelaxingMedium::from_thermochemical(resolved).map_err(ambient)?
        }
    };
    let present: Vec<_> = medium.components().iter().filter(|(_, x)| *x > 0.0).collect();
    let figure = if TerrestrialAir::admit(&medium).is_ok() {
        iso(fs_material::iso9613::iso9613_relative_humidity(
            gas.temperature, gas.pressure, gas.water_mole_fraction,
        )
        .map_err(ambient)?)?
    } else if present.iter().all(|(id, _)| MONATOMIC.contains(&id.as_str())) {
        ObserverAbsorptionFigure::ClassicalComplete { db: classical_db }
    } else if let Some((id, _)) = present.iter().find(|(id, _)| !HELD_RATE_SPECIES.contains(&id.as_str())) {
        ObserverAbsorptionFigure::Unbounded { species: id.as_str().to_owned() }
    } else {
        let fraction = |name: &str| present.iter().find(|(id, _)| id.as_str() == name).map_or(0.0, |(_, x)| *x);
        let dry = 1.0 - fraction("H2O");
        if fraction("CO2") / dry > fs_material::iso9613::ISO9613_MAX_DRY_CO2 {
            ObserverAbsorptionFigure::Unbounded { species: "CO2".to_owned() }
        } else {
            let state = RelaxationAbsorption::evaluate_iso9613_rates_held(&medium).map_err(ambient)?;
            let alpha = state.causal_wavenumber(omega).map_err(ambient)?.im;
            ObserverAbsorptionFigure::HeldIsoRates {
                o2: fraction("O2"),
                n2: fraction("N2"),
                water: fraction("H2O"),
                db: over_range(alpha),
            }
        }
    };
    Ok(ObserverOmission { range_m, band_top_hz, figure, classical_db })
}

/// Propagate a pressure history through the evaluated gas's molecular relaxation.
///
/// Uses the same energy-bounded causal FIR as [`ObserverFilter::molecular_relaxation`],
/// on the exact clock `1 / dt`, and linear convolution with zero input before
/// the record. A late impulse cannot wrap around into its beginning. Output
/// after the supplied slice is discarded; short histories are filtered too.
///
/// The FIR includes its Nyquist alignment and precursor-padding latency.
/// Returns the receiver delay still to apply: `r / c_0` minus the FIR's DC
/// group delay, exactly as [`ObserverFilterReceipt::receiver_delay_s`]. This
/// replaces the former bare `r / c_inf` return, which did not account for
/// realization latency. A negative return remains visible: that clock and
/// range cannot place the causal FIR at the physical arrival without an advance.
///
/// No receiver delay is applied here. The retained samples are the same
/// source-time history as the streaming filter, up to FFT roundoff. A zero
/// range leaves the history untouched and returns zero.
///
/// # Errors
/// Invalid samples, time/range, an unrepresentable response or output, or a
/// response exceeding the observer FIR's energy/length limits refuses WITHOUT
/// changing the history.
pub fn propagate_pressure_history_relaxing(
    pressure_pa: &mut [f64],
    dt: f64,
    range_m: f64,
    relaxation: &RelaxationState,
) -> Result<f64, AcousticRealizeError> {
    if !dt.is_finite() || dt <= 0.0 || !range_m.is_finite() || range_m < 0.0
        || pressure_pa.iter().any(|sample| !sample.is_finite())
    {
        return Err(invalid("propagation needs finite pressure, positive dt and nonnegative range"));
    }
    if range_m == 0.0 { return Ok(0.0); }
    let response = molecular_response(relaxation, range_m, 1.0 / dt)?;
    let delay_s = response.receiver_delay_s(relaxation.gas(), range_m);
    if !delay_s.is_finite() {
        return Err(invalid("relaxation receiver delay is unrepresentable"));
    }
    convolve_pressure_history(pressure_pa, &response.taps)?;
    Ok(delay_s)
}

// Overlap-add with N >= block_len + taps - 1 realizes a linear convolution.
// Only the scratch size depends on the FIR length; the separate output keeps
// refusal atomic. Taps beyond this record cannot affect any retained sample.
fn convolve_pressure_history(
    pressure_pa: &mut [f64],
    taps: &[f64],
) -> Result<(), AcousticRealizeError> {
    if pressure_pa.is_empty() { return Ok(()); }
    let m = taps.len().min(pressure_pa.len());
    let mut output = vec![0.0; pressure_pa.len()];
    if m <= 32 {
        for (i, value) in output.iter_mut().enumerate() {
            for (j, tap) in taps[..m.min(i + 1)].iter().enumerate() {
                *value += tap * pressure_pa[i - j];
            }
        }
    } else {
        let n_fft = m.checked_next_power_of_two()
            .and_then(|n| n.checked_mul(2))
            .ok_or_else(|| invalid("propagation FFT length is unrepresentable"))?;
        let block_len = n_fft - m + 1;
        let fft = Fft::new(n_fft);
        let mut kernel = vec![FftC64::new(0.0, 0.0); n_fft];
        let mut buf = kernel.clone();
        let mut scratch = kernel.clone();
        for (bin, &tap) in kernel.iter_mut().zip(&taps[..m]) {
            bin.re = tap;
        }
        fft.forward(&mut kernel, &mut scratch);
        for start in (0..pressure_pa.len()).step_by(block_len) {
            let count = block_len.min(pressure_pa.len() - start);
            buf.fill(FftC64::new(0.0, 0.0));
            for (bin, &sample) in buf.iter_mut().zip(&pressure_pa[start..start + count]) {
                bin.re = sample;
            }
            fft.forward(&mut buf, &mut scratch);
            for (bin, h) in buf.iter_mut().zip(&kernel) {
                *bin = FftC64::new(
                    bin.re * h.re - bin.im * h.im,
                    bin.re * h.im + bin.im * h.re,
                );
            }
            fft.inverse(&mut buf, &mut scratch);
            let kept = (count + m - 1).min(output.len() - start);
            for (out, value) in output[start..start + kept].iter_mut().zip(&buf) {
                *out += value.re;
            }
        }
    }
    if output.iter().any(|value| !value.is_finite()) {
        return Err(invalid("propagated pressure is unrepresentable"));
    }
    pressure_pa.copy_from_slice(&output);
    Ok(())
}

/// The relaxation law's transfer over `range_m` at `omega`, referenced to the
/// frozen speed `c_inf`, for the DFT's `e^{+iwt}` synthesis: `(Re, Im)` of
/// `exp(-Im k r) exp(-i (Re k - w / c_inf) r)`.
fn relaxing_transfer(
    relaxation: &RelaxationState,
    omega: f64,
    range_m: f64,
    c_inf: f64,
) -> Result<(f64, f64), AcousticRealizeError> {
    let wave = relaxation.causal_wavenumber(omega)
        .map_err(|_| invalid("relaxation wavenumber could not be evaluated"))?;
    let magnitude = det::exp(-wave.im * range_m);
    let (sin, cos) = det::sin_cos(-(wave.re - omega / c_inf) * range_m);
    Ok((magnitude * cos, magnitude * sin))
}

/// Largest causal FIR an observer law may need.
pub const MAX_OBSERVER_TAPS: usize = 16_384;
/// Impulse-response energy an observer FIR may leave out, relative to its
/// total (the truncated tail plus the precursor before the declared latency),
/// at most. The certificate a receipt applies is the smaller of this and
/// [`OBSERVER_UNCERTAINTY_MARGIN`] times what ISO 9613-1's stated +-10 % on
/// alpha is in the same measure (orchestrator ruling of 2026-10-23, 6352).
pub const OBSERVER_TRUNCATION_ENERGY: f64 = 1e-9;
/// How far below the law's own alpha uncertainty the dropped energy must sit.
pub const OBSERVER_UNCERTAINTY_MARGIN: f64 = 1e-3;
/// `omega tau_v` at and above which the classical viscothermal law
/// (`alpha ~ omega^2`) is used outside its assumption `omega tau_v << 1`.
pub const CLASSICAL_VALIDITY_OMEGA_TAU: f64 = 0.01;
const MAX_OBSERVER_DESIGN_LEN: usize = 1 << 20;

/// What a streaming observer filter realizes, with its certificates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ObserverFilterReceipt {
    /// The law and any molecular relaxation it leaves out.
    pub absorption: ObserverAbsorptionReceipt,
    /// Listener distance the filter covers [m].
    pub range_m: f64,
    /// Clock of the filtered samples [Hz].
    pub sample_rate_hz: u32,
    /// FIR length.
    pub taps: usize,
    /// Impulse-response energy left out (precursor before the padding plus
    /// the truncated tail), relative to the total; at most
    /// [`Self::budget_energy`].
    pub truncated_energy: f64,
    /// The certificate this receipt applies:
    /// `min(OBSERVER_TRUNCATION_ENERGY, OBSERVER_UNCERTAINTY_MARGIN *
    /// alpha_uncertainty_energy)`.
    pub budget_energy: f64,
    /// What ISO 9613-1's stated +-10 % accuracy on alpha is in the same
    /// measure: `sum (|H| - |H|^1.1)^2 / sum |H|^2` over the design grid,
    /// with this law's own magnitude. The certificate is meaningful while
    /// `truncated_energy` sits orders of magnitude below it.
    pub alpha_uncertainty_energy: f64,
    /// Sub-sample alignment `d` in [0, 1) that makes the Nyquist bin real.
    pub alignment_samples: f64,
    /// Whole samples of padding `D`: the fewest for the precursor to meet
    /// half the budget.
    pub padding_samples: usize,
    /// DC group delay of the realized FIR [s], its first moment
    /// `sum n h[n] / sum h[n]`: the law's low-frequency excess over its front
    /// (`r (1/c_0 - 1/c_inf)` for a relaxing gas) plus [`Self::latency_s`].
    pub low_frequency_delay_s: f64,
    /// Declared delay that is not physics [s]. A causal transfer sampled only
    /// to Nyquist rings before `t = 0`, so a dispersive law's FIR is delayed by
    /// a sub-sample alignment (making its Nyquist bin real) plus whole samples
    /// of padding, just enough for the precursor to meet the certificate. It is
    /// zero for a minimum-phase law.
    pub latency_s: f64,
    /// Propagation delay a receiver must still apply [s] so that its delay
    /// plus the filter's low-frequency delay is `r / c_0` exactly once:
    /// `r / c_0 - low_frequency_delay_s`. Zero or negative when no receiver at
    /// this range can absorb the filter; a baffled receiver then refuses.
    pub receiver_delay_s: f64,
    /// `r / receiver_delay_s`, the speed a receiver's propagation delays use,
    /// when that delay is positive. It is `c_0` for a lossless law. A receiver
    /// with no propagation delay (the legacy proxy) leaves the filter's delay
    /// to its consumer.
    pub delay_speed_m_s: Option<f64>,
    /// Lowest frequency [Hz] below the clock's Nyquist at which the classical
    /// viscothermal term is used with `omega tau_v` at or above
    /// [`CLASSICAL_VALIDITY_OMEGA_TAU`], `tau_v = (4 mu / 3 + (gamma - 1)
    /// kappa / c_p) / (rho c^2)`. `None` when the assumption holds over the
    /// whole band or the law has no classical `omega^2` term.
    pub classical_validity_exceeded_from_hz: Option<f64>,
}

/// A causal FIR realization of an observer law over a fixed range, for
/// renderers that emit pressure sample by sample. The sum runs in a fixed
/// order, so output is bitwise independent of how samples are blocked.
///
/// How each law is made causal:
/// - Molecular relaxation and ISO 9613-1's relaxation terms carry their own
///   dispersion: the relaxing law's causal wavenumber, and the exact complex
///   Debye pair of each ISO term. Both are referenced to their front.
/// - The classical `omega^2` terms (Stokes-Kirchhoff, ISO's classical term)
///   are realized as the in-band minimum phase of their magnitude.
///
/// A minimum phase built from the in-band magnitude alone is NOT used where
/// dispersion exists. It cannot see absorption above Nyquist (O2 relaxes near
/// 35 kHz at RH 50 %). For air at 30 m against the direct causal transfer, its
/// worst in-band complex error was 2.0 / 1.07 / 0.11 on 24 / 48 / 192 kHz
/// Nyquist grids, and its low-frequency delay 20.2 / 27.4 / 33.8 us against
/// the physical 34.5 us.
#[derive(Debug, Clone)]
pub struct ObserverFilter {
    taps: Vec<f64>,
    history: Vec<f64>,
    head: usize,
    receipt: ObserverFilterReceipt,
}

// Response design is independent of the streaming receipt's integer clock.
// Batch propagation supplies its exact reciprocal dt, without quantizing it.
struct ObserverResponse {
    taps: Vec<f64>,
    truncated_energy: f64,
    budget_energy: f64,
    alpha_uncertainty_energy: f64,
    alignment_samples: f64,
    padding_samples: usize,
    low_frequency_delay_s: f64,
    latency_s: f64,
}

impl ObserverResponse {
    fn receiver_delay_s(&self, gas: &GasState, range_m: f64) -> f64 {
        range_m / gas.sound_speed - self.low_frequency_delay_s
    }
}

fn molecular_response(
    relaxation: &RelaxationState,
    range_m: f64,
    rate: f64,
) -> Result<ObserverResponse, AcousticRealizeError> {
    let c_inf = relaxation.frozen_sound_speed_m_s();
    design_observer_response(range_m, rate, |n| {
        let mut bins = vec![FftC64::new(1.0, 0.0); n / 2 + 1];
        for (k, bin) in bins.iter_mut().enumerate().skip(1) {
            let omega = core::f64::consts::TAU * k as f64 * rate / n as f64;
            let (re, im) = relaxing_transfer(relaxation, omega, range_m, c_inf)?;
            *bin = FftC64::new(re, im);
        }
        Ok(bins)
    })
}

// Build the positive-frequency transfer on growing DFT grids, align its
// Nyquist bin to be real, and keep the shortest FIR (precursor padding D
// plus causal part) whose left-out energy meets the certificate.
fn design_observer_response(
    range_m: f64,
    rate: f64,
    mut spectrum: impl FnMut(usize) -> Result<Vec<FftC64>, AcousticRealizeError>,
) -> Result<ObserverResponse, AcousticRealizeError> {
    if !(range_m.is_finite() && range_m >= 0.0 && rate.is_finite() && rate > 0.0) {
        return Err(invalid("an observer filter needs a finite nonnegative range and a positive finite rate"));
    }
    if range_m == 0.0 {
        return Ok(ObserverResponse {
            taps: vec![1.0],
            truncated_energy: 0.0,
            budget_energy: 0.0,
            alpha_uncertainty_energy: 0.0,
            alignment_samples: 0.0,
            padding_samples: 0,
            low_frequency_delay_s: 0.0,
            latency_s: 0.0,
        });
    }
    let mut n = 1usize << 14;
    loop {
        let half = spectrum(n)?;
        let (mut magnitude_energy, mut uncertainty) = (0.0, 0.0);
        for (k, h) in half.iter().enumerate() {
            let weight = if k == 0 || k == n / 2 { 1.0 } else { 2.0 };
            let magnitude = h.re.hypot(h.im);
            let perturbed = det::exp(1.1 * det::ln(magnitude));
            magnitude_energy += weight * magnitude * magnitude;
            uncertainty += weight * (magnitude - perturbed) * (magnitude - perturbed);
        }
        let alpha_uncertainty_energy = uncertainty / magnitude_energy;
        let budget = OBSERVER_TRUNCATION_ENERGY.min(OBSERVER_UNCERTAINTY_MARGIN * alpha_uncertainty_energy);
        // A sub-sample delay in [0, 1) that makes the Nyquist bin real, so
        // the periodic spectrum has no jump there. A minimum phase already
        // has a real Nyquist bin.
        let nyquist = half[n / 2];
        let align = if nyquist.im.abs() <= 1e-12 * nyquist.re.hypot(nyquist.im) {
            0.0
        } else {
            nyquist.im.atan2(nyquist.re).rem_euclid(core::f64::consts::PI) / core::f64::consts::PI
        };
        let mut buf = vec![FftC64::new(0.0, 0.0); n];
        buf[0] = FftC64::new(1.0, 0.0);
        for k in 1..=n / 2 {
            let (sin, cos) = det::sin_cos(-core::f64::consts::TAU * k as f64 * align / n as f64);
            let h = half[k];
            let im = if k == n / 2 { 0.0 } else { h.re * sin + h.im * cos };
            buf[k] = FftC64::new(h.re * cos - h.im * sin, im);
        }
        for k in 1..n / 2 {
            let c = buf[k];
            buf[n - k] = FftC64::new(c.re, -c.im);
        }
        let mut scratch = vec![FftC64::new(0.0, 0.0); n];
        Fft::new(n).inverse(&mut buf, &mut scratch);
        let response: Vec<f64> = buf.iter().map(|value| value.re).collect();
        if response.iter().any(|value| !value.is_finite()) {
            return Err(invalid("observer impulse response is unrepresentable"));
        }
        let energy = |range: core::ops::Range<usize>| -> f64 {
            response[range].iter().map(|v| v * v).sum()
        };
        let total = energy(0..n);
        if !(total > 0.0) {
            return Err(invalid("observer impulse response has no energy"));
        }
        // The grid is long enough once its middle half is quiet: the causal
        // tail has decayed and the precursor has not yet risen.
        if energy(n / 4..3 * n / 4) > budget * total / 4.0 {
            if n >= MAX_OBSERVER_DESIGN_LEN {
                return Err(invalid("observer response does not decay within the design grid"));
            }
            n *= 2;
            continue;
        }
        // Pad the fewest samples D for the precursor (negative time, the
        // grid's second half) to meet half the certificate.
        let mut padding = 0;
        let mut precursor = energy(n / 2..n);
        while precursor > budget * total / 2.0 && padding < n / 4 {
            padding += 1;
            precursor -= response[n - padding] * response[n - padding];
        }
        if precursor > budget * total / 2.0 {
            if n >= MAX_OBSERVER_DESIGN_LEN {
                return Err(invalid("observer precursor does not meet the certificate within the design grid"));
            }
            n *= 2;
            continue;
        }
        let mut tail = 0.0;
        let mut causal = n / 2;
        while causal > 1 {
            let next = tail + response[causal - 1] * response[causal - 1];
            if precursor + next > budget * total {
                break;
            }
            tail = next;
            causal -= 1;
        }
        let taps = padding + causal;
        if taps > MAX_OBSERVER_TAPS {
            return Err(invalid("observer response needs more taps than MAX_OBSERVER_TAPS"));
        }
        let mut fir = response[n - padding..].to_vec();
        fir.extend_from_slice(&response[..causal]);
        let latency_s = (align + padding as f64) / rate;
        let (moment, gain) = fir.iter().enumerate().fold((0.0, 0.0), |(m, g), (j, tap)| {
            (m + j as f64 * tap, g + tap)
        });
        let low_frequency_delay_s = moment / gain / rate;
        if !latency_s.is_finite() || !low_frequency_delay_s.is_finite() {
            return Err(invalid("observer response delay is unrepresentable"));
        }
        return Ok(ObserverResponse {
            taps: fir,
            truncated_energy: (precursor + tail) / total,
            budget_energy: budget,
            alpha_uncertainty_energy,
            alignment_samples: align,
            padding_samples: padding,
            low_frequency_delay_s,
            latency_s,
        });
    }
}

/// The excess transfer of one Debye relaxation with absorption
/// `strength f_r x^2 / (1 + x^2)` [Np/m], `x = f / f_r`, for the DFT's
/// `e^{+iwt}` synthesis. Its complex wavenumber is the exact Kramers-Kronig
/// pair `strength f_r x / (1 - i x)` (`e^{-iwt}`, pole in the lower half plane),
/// referenced to the front: the low-frequency slowness is `strength / 2 pi`.
fn debye_transfer(strength: f64, fr_hz: f64, f_hz: f64, range_m: f64) -> (f64, f64) {
    let x = f_hz / fr_hz;
    let denominator = 1.0 + x * x;
    let re_k = strength * fr_hz * x / denominator;
    let im_k = strength * fr_hz * x * x / denominator;
    let magnitude = det::exp(-im_k * range_m);
    let (sin, cos) = det::sin_cos(-re_k * range_m);
    (magnitude * cos, magnitude * sin)
}

/// The minimum-phase spectrum, bins `0..=n/2` for the DFT's `e^{+iwt}`
/// synthesis, of a real even log-magnitude given on those bins, by the real
/// cepstrum folded onto positive quefrency.
fn minimum_phase(log_magnitude: &[f64], n: usize) -> Vec<FftC64> {
    let fft = Fft::new(n);
    let mut buf: Vec<FftC64> = (0..n)
        .map(|k| FftC64::new(log_magnitude[if k <= n / 2 { k } else { n - k }], 0.0))
        .collect();
    let mut scratch = vec![FftC64::new(0.0, 0.0); n];
    fft.inverse(&mut buf, &mut scratch);
    for (m, value) in buf.iter_mut().enumerate() {
        let weight = if m == 0 || m == n / 2 { 1.0 } else if m < n / 2 { 2.0 } else { 0.0 };
        *value = FftC64::new(weight * value.re, 0.0);
    }
    fft.forward(&mut buf, &mut scratch);
    buf.iter()
        .take(n / 2 + 1)
        .map(|z| {
            let magnitude = det::exp(z.re);
            let (sin, cos) = det::sin_cos(z.im);
            FftC64::new(magnitude * cos, magnitude * sin)
        })
        .collect()
}

/// Viscothermal relaxation time `(4 mu / 3 + (gamma - 1) kappa / c_p) /
/// (rho c^2)` and the lowest frequency below Nyquist where `omega tau_v`
/// reaches [`CLASSICAL_VALIDITY_OMEGA_TAU`].
fn classical_validity(gas: &GasState, sample_rate_hz: u32) -> Option<f64> {
    let diffusivity = 4.0 / 3.0 * gas.dynamic_viscosity
        + (gas.gamma - 1.0) * gas.thermal_conductivity / gas.specific_heat_cp;
    let tau_v = diffusivity / (gas.density * gas.sound_speed * gas.sound_speed);
    let limit_hz = CLASSICAL_VALIDITY_OMEGA_TAU / (core::f64::consts::TAU * tau_v);
    (limit_hz < 0.5 * f64::from(sample_rate_hz)).then_some(limit_hz)
}

impl ObserverFilter {
    /// The gas's own molecular relaxation over `range_m`: the causal transfer
    /// of [`propagate_pressure_history_relaxing`], referenced to the frozen
    /// speed.
    ///
    /// # Errors
    /// A non-finite or negative range, a zero rate, an unevaluable wavenumber,
    /// or a response that needs more than [`MAX_OBSERVER_TAPS`] taps.
    pub fn molecular_relaxation(
        relaxation: &RelaxationState,
        range_m: f64,
        sample_rate_hz: u32,
    ) -> Result<Self, AcousticRealizeError> {
        let response = molecular_response(relaxation, range_m, f64::from(sample_rate_hz))?;
        let absorption = ObserverAbsorptionReceipt {
            law: GasAbsorptionModel::MolecularRelaxation,
            relaxation_omitted: None,
        };
        Ok(Self::from_response(
            response, absorption, relaxation.gas(), range_m, sample_rate_hz, None,
        ))
    }

    /// ISO 9613-1 for an admitted terrestrial-air `gas`: each O2 and N2 term
    /// by its exact complex Debye pair (magnitude ISO's exactly, dispersion its
    /// own), and the classical term as the in-band minimum phase.
    ///
    /// # Errors
    /// [`AcousticRealizeError::NonTerrestrialAir`] when `air` was not admitted
    /// for `gas`, ISO's window refusals, and the design refusals of
    /// [`Self::molecular_relaxation`].
    pub fn iso9613(
        air: &TerrestrialAir,
        gas: &GasState,
        relative_humidity: f64,
        range_m: f64,
        sample_rate_hz: u32,
    ) -> Result<Self, AcousticRealizeError> {
        air.check(gas).map_err(AcousticRealizeError::NonTerrestrialAir)?;
        let terms = fs_material::iso9613::iso9613_terms(gas.temperature, gas.pressure, relative_humidity)
            .map_err(|_| invalid("explicit ISO 9613 air conditions are outside the model domain"))?;
        let rate = f64::from(sample_rate_hz);
        let absorption = ObserverAbsorptionReceipt {
            law: GasAbsorptionModel::Iso9613Air { relative_humidity, air: *air },
            relaxation_omitted: None,
        };
        let validity = classical_validity(gas, sample_rate_hz);
        Self::design(absorption, gas, range_m, sample_rate_hz, validity, |n| {
            let log_classical: Vec<f64> = (0..=n / 2)
                .map(|k| {
                    let f = k as f64 * rate / n as f64;
                    -terms.classical_np_per_m_hz2 * f * f * range_m
                })
                .collect();
            let mut bins = minimum_phase(&log_classical, n);
            for (k, bin) in bins.iter_mut().enumerate().skip(1) {
                let f = k as f64 * rate / n as f64;
                for (strength, fr_hz) in terms.relaxation {
                    let (re, im) = debye_transfer(strength, fr_hz, f, range_m);
                    *bin = FftC64::new(bin.re * re - bin.im * im, bin.re * im + bin.im * re);
                }
            }
            Ok(bins)
        })
    }

    /// Stokes-Kirchhoff, the gas's classical viscothermal law, as the in-band
    /// minimum phase of its magnitude.
    ///
    /// # Errors
    /// An incoherent gas and the design refusals.
    pub fn stokes_kirchhoff(
        gas: &GasState,
        range_m: f64,
        sample_rate_hz: u32,
    ) -> Result<Self, AcousticRealizeError> {
        gas.validate_coherence().map_err(|_| invalid("observer absorption needs a coherent gas state"))?;
        let rate = f64::from(sample_rate_hz);
        let absorption = ObserverAbsorptionReceipt {
            law: GasAbsorptionModel::StokesKirchhoff,
            relaxation_omitted: Some(relaxation_omission(range_m, 0.5 * rate)?),
        };
        let validity = classical_validity(gas, sample_rate_hz);
        Self::design(absorption, gas, range_m, sample_rate_hz, validity, |n| {
            let log_magnitude: Vec<f64> = (0..=n / 2)
                .map(|k| -gas.stokes_kirchhoff_absorption(core::f64::consts::TAU * k as f64 * rate / n as f64) * range_m)
                .collect();
            Ok(minimum_phase(&log_magnitude, n))
        })
    }

    /// A declared lossless observation (`observer_law none`), for comparisons
    /// against a lossless analytic oracle. It is receipted like any other law:
    /// everything the gas would absorb is labelled as omitted.
    ///
    /// # Errors
    /// An invalid range or rate.
    pub fn lossless(gas: &GasState, range_m: f64, sample_rate_hz: u32) -> Result<Self, AcousticRealizeError> {
        if !(range_m.is_finite() && range_m >= 0.0) || sample_rate_hz == 0 {
            return Err(invalid("an observer filter needs a finite nonnegative range and a positive rate"));
        }
        let absorption = ObserverAbsorptionReceipt {
            law: GasAbsorptionModel::Lossless,
            relaxation_omitted: Some(relaxation_omission(range_m, 0.5 * f64::from(sample_rate_hz))?),
        };
        Ok(Self {
            taps: vec![1.0],
            history: vec![0.0],
            head: 0,
            receipt: ObserverFilterReceipt {
                absorption,
                range_m,
                sample_rate_hz,
                taps: 1,
                truncated_energy: 0.0,
                low_frequency_delay_s: 0.0,
                latency_s: 0.0,
                budget_energy: 0.0,
                alpha_uncertainty_energy: 0.0,
                alignment_samples: 0.0,
                padding_samples: 0,
                receiver_delay_s: range_m / gas.sound_speed,
                delay_speed_m_s: Some(gas.sound_speed),
                classical_validity_exceeded_from_hz: None,
            },
        })
    }

    /// The filter for `law`, with the relaxation state that molecular
    /// relaxation needs.
    ///
    /// # Errors
    /// Molecular relaxation without its state, and each law's own refusals.
    pub fn for_law(
        law: GasAbsorptionModel,
        gas: &GasState,
        relaxation: Option<&RelaxationState>,
        range_m: f64,
        sample_rate_hz: u32,
    ) -> Result<Self, AcousticRealizeError> {
        match law {
            GasAbsorptionModel::StokesKirchhoff => Self::stokes_kirchhoff(gas, range_m, sample_rate_hz),
            GasAbsorptionModel::Iso9613Air { relative_humidity, air } => {
                Self::iso9613(&air, gas, relative_humidity, range_m, sample_rate_hz)
            }
            GasAbsorptionModel::MolecularRelaxation => {
                let state = relaxation
                    .ok_or_else(|| invalid("molecular relaxation needs its evaluated relaxation state"))?;
                if state.gas() != gas {
                    return Err(invalid("the observer relaxation state was evaluated for a different gas state"));
                }
                Self::molecular_relaxation(state, range_m, sample_rate_hz)
            }
            GasAbsorptionModel::Lossless => Self::lossless(gas, range_m, sample_rate_hz),
        }
    }

    fn design(
        absorption: ObserverAbsorptionReceipt,
        gas: &GasState,
        range_m: f64,
        sample_rate_hz: u32,
        classical_validity_exceeded_from_hz: Option<f64>,
        spectrum: impl FnMut(usize) -> Result<Vec<FftC64>, AcousticRealizeError>,
    ) -> Result<Self, AcousticRealizeError> {
        let response = design_observer_response(range_m, f64::from(sample_rate_hz), spectrum)?;
        Ok(Self::from_response(
            response, absorption, gas, range_m, sample_rate_hz,
            classical_validity_exceeded_from_hz,
        ))
    }

    fn from_response(
        response: ObserverResponse,
        absorption: ObserverAbsorptionReceipt,
        gas: &GasState,
        range_m: f64,
        sample_rate_hz: u32,
        classical_validity_exceeded_from_hz: Option<f64>,
    ) -> Self {
        let taps = response.taps.len();
        let receiver_delay_s = response.receiver_delay_s(gas, range_m);
        let delay_speed_m_s = if range_m == 0.0 {
            Some(gas.sound_speed)
        } else {
            (receiver_delay_s > 0.0).then(|| range_m / receiver_delay_s)
        };
        Self {
            taps: response.taps,
            history: vec![0.0; taps],
            head: 0,
            receipt: ObserverFilterReceipt {
                absorption,
                range_m,
                sample_rate_hz,
                taps,
                truncated_energy: response.truncated_energy,
                budget_energy: response.budget_energy,
                alpha_uncertainty_energy: response.alpha_uncertainty_energy,
                alignment_samples: response.alignment_samples,
                padding_samples: response.padding_samples,
                low_frequency_delay_s: response.low_frequency_delay_s,
                latency_s: response.latency_s,
                receiver_delay_s,
                delay_speed_m_s,
                classical_validity_exceeded_from_hz,
            },
        }
    }

    /// The filter's receipt.
    #[must_use]
    pub const fn receipt(&self) -> &ObserverFilterReceipt {
        &self.receipt
    }

    /// The FIR, earliest tap first.
    #[must_use]
    pub fn taps(&self) -> &[f64] {
        &self.taps
    }

    /// Filter one sample on the receipt's clock.
    pub fn process(&mut self, sample: f64) -> f64 {
        let len = self.taps.len();
        self.history[self.head] = sample;
        let mut sum = 0.0;
        for (j, tap) in self.taps.iter().enumerate() {
            let index = if j <= self.head { self.head - j } else { self.head + len - j };
            sum += tap * self.history[index];
        }
        self.head = if self.head + 1 == len { 0 } else { self.head + 1 };
        sum
    }
}

// One transform implementation for explicit and legacy model selection.
fn transform_history(
    pressure_pa: &mut [f64],
    dt: f64,
    range_m: f64,
    mut absorption: impl FnMut(f64) -> Result<f64, AcousticRealizeError>,
) -> Result<(), AcousticRealizeError> {
    if !dt.is_finite() || dt <= 0.0 || !range_m.is_finite() || range_m < 0.0
        || pressure_pa.iter().any(|sample| !sample.is_finite())
    {
        return Err(invalid("absorption needs finite pressure, positive dt and nonnegative range"));
    }
    if pressure_pa.len() < 4 || range_m == 0.0 { return Ok(()); }
    let n = pressure_pa.len();
    // Cover the complete history, including samples beyond the old FFT cap.
    let n_fft = n.checked_next_power_of_two()
        .ok_or_else(|| invalid("absorption FFT length is unrepresentable"))?.max(8);
    let fft = Fft::new(n_fft);
    let mut buf: Vec<FftC64> = (0..n_fft)
        .map(|i| FftC64::new(if i < n { pressure_pa[i] } else { 0.0 }, 0.0))
        .collect();
    let mut scratch = vec![FftC64::new(0.0, 0.0); n_fft];
    fft.forward(&mut buf, &mut scratch);
    for (k, bin) in buf.iter_mut().enumerate().take(n_fft / 2 + 1) {
        if k == 0 { continue; }
        let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
        if !omega.is_finite() || omega <= 0.0 {
            return Err(invalid("absorption frequency is unrepresentable"));
        }
        let alpha = absorption(omega)?;
        if !alpha.is_finite() || alpha < 0.0 {
            return Err(invalid("absorption must be finite and nonnegative"));
        }
        let atten = det::exp(-alpha * range_m);
        *bin = FftC64::new(bin.re * atten, bin.im * atten);
    }
    for k in 1..n_fft / 2 {
        let c = buf[k];
        buf[n_fft - k] = FftC64::new(c.re, -c.im);
    }
    fft.inverse(&mut buf, &mut scratch);
    if buf[..n].iter().any(|value| !value.re.is_finite()) {
        return Err(invalid("absorbed pressure is unrepresentable"));
    }
    for (i, sample) in pressure_pa.iter_mut().enumerate() { *sample = buf[i].re; }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::gas::{GasSpec, GasState};

    /// Dry USSA air at `t` and 1 atm with its terrestrial-air admission, and
    /// the ISO law at 50 % RH.
    fn iso_air(t: f64) -> (GasState, GasAbsorptionModel) {
        let (gas, air) = TerrestrialAir::moist_air(t, 101_325.0, 0.0).expect("air");
        (gas, GasAbsorptionModel::Iso9613Air { relative_humidity: 0.50, air })
    }

    #[test]
    fn farther_path_kills_high_frequency_more() {
        let (gas, iso) = iso_air(288.15);
        let n = 256;
        let dt = 1.0 / 8_000.0;
        let mut near: Vec<f64> = (0..n)
            .map(|i| (core::f64::consts::TAU * 2_000.0 * f64::from(i) * dt).sin())
            .collect();
        let mut far = near.clone();
        absorb_pressure_history_with_model(&mut near, dt, 1.0, &gas, iso).unwrap();
        absorb_pressure_history_with_model(&mut far, dt, 200.0, &gas, iso).unwrap();
        let en: f64 = near.iter().map(|x| x * x).sum();
        let ef: f64 = far.iter().map(|x| x * x).sum();
        assert!(ef < en * 0.98, "200 m must absorb more 2 kHz energy than 1 m");
    }

    #[test]
    fn g0_history_beyond_legacy_fft_cap_is_transformed_in_full() {
        let (gas, iso) = iso_air(288.15);
        let mut pressure = vec![0.0; 8_193];
        pressure[8_192] = 1.0;
        absorb_pressure_history_with_model(&mut pressure, 1.0 / 48_000.0, 200.0, &gas, iso).unwrap();
        assert_eq!(pressure.len(), 8_193);
        assert!(pressure.iter().all(|sample| sample.is_finite()));
        assert!(pressure[8_192].abs() < 1.0);
    }

    #[test]
    fn generic_transfer_uses_the_actual_material_for_air_and_steam() {
        let dt = 1.0 / 8000.0;
        let omega = core::f64::consts::TAU * 1000.0;
        let input: Vec<_> = (0..256).map(|i| (omega * f64::from(i) * dt).cos()).collect();
        let mut outputs = Vec::new();
        for spec in [GasSpec::dry_air_ussa1976(), GasSpec::water_vapor_nist()] {
            let gas = GasState::try_new(&spec, 450.0, 100000.0).unwrap();
            let expected = det::exp(-gas.stokes_kirchhoff_absorption(omega) * 100.0);
            let mut signal = input.clone();
            absorb_pressure_history_with_model(&mut signal, dt, 100.0, &gas,
                GasAbsorptionModel::StokesKirchhoff).unwrap();
            for (&actual, &before) in signal.iter().zip(&input) {
                assert!((actual - expected * before).abs() < 1e-10);
            }
            outputs.push(signal);
        }
        assert_ne!(outputs[0], outputs[1]);
    }

    /// ISO applied to steam refuses with the named refusal, leaving the history
    /// unchanged: an air admission does not transfer to another gas.
    #[test]
    fn explicit_air_refuses_instead_of_silently_using_another_law() {
        let gas = GasState::try_new(&GasSpec::water_vapor_nist(), 450.0, 100000.0).unwrap();
        let (_, iso) = iso_air(293.15);
        let original = vec![1.0, -1.0, 0.5, -0.5, 1.0, 0.0, -1.0, 0.0];
        let mut signal = original.clone();
        assert_eq!(
            absorb_pressure_history_with_model(&mut signal, 1.0 / 8000.0, 100.0, &gas, iso),
            Err(AcousticRealizeError::NonTerrestrialAir(
                fs_material::iso9613::Iso9613Refusal::DifferentState
            ))
        );
        assert_eq!(signal, original);
        let mut forged = gas;
        forged.density *= 2.0;
        assert!(absorb_pressure_history_with_model(&mut signal, 1.0 / 8000.0, 100.0,
            &forged, GasAbsorptionModel::StokesKirchhoff).is_err());
        assert_eq!(signal, original);
    }

    /// The admitted ISO path multiplies a single tone by exactly the
    /// formula's `exp(-alpha r)`.
    #[test]
    fn admitted_iso_transfer_is_the_formula() {
        let dt = 1.0 / 8000.0;
        let omega = core::f64::consts::TAU * 1000.0;
        let input: Vec<_> = (0..256).map(|i| (omega * f64::from(i) * dt).cos()).collect();
        let (gas, iso) = iso_air(300.0);
        let alpha = fs_material::iso9613::iso9613_absorption_neper_per_m(300.0, 101_325.0, 0.5, omega)
            .unwrap();
        let expected = det::exp(-alpha * 100.0);
        let mut signal = input.clone();
        absorb_pressure_history_with_model(&mut signal, dt, 100.0, &gas, iso).unwrap();
        for (&actual, &before) in signal.iter().zip(&input) {
            assert!((actual - expected * before).abs() < 1e-10);
        }
    }

    use fs_material::relaxation::{
        InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
        RelaxationState, RelaxingMedium, SpeciesRelaxation, VibrationalMode,
    };
    use fs_thermochem::SpeciesId;

    /// One synthetic reservoir (theta 960 K, degeneracy 2, tau 1e-4 s), or none,
    /// in dry USSA air: a test medium for the transfer, not gas data.
    fn relaxing_air(reservoir: bool) -> RelaxationState {
        let p = 101_325.0;
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, p).unwrap();
        let id = SpeciesId::new("A").unwrap();
        let processes = if reservoir {
            vec![RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode { theta_k: 960.0, degeneracy: 2 }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 1e-4 * p,
                    t_ref_k: 293.15,
                    exponent: 0.0,
                    t_min_k: 100.0,
                    t_max_k: 2000.0,
                },
            }]
        } else {
            Vec::new()
        };
        let law = RelaxationAbsorption::new(vec![SpeciesRelaxation { species: id.clone(), processes }]).unwrap();
        law.evaluate(&RelaxingMedium::new(gas, vec![(id, 1.0)]).unwrap()).unwrap()
    }

    /// G1: the batch FIR carries the causal law's complex phase, including
    /// its declared realization latency, within the streaming FIR's 1e-3
    /// absolute transfer tolerance. Receiver delay restores r/c0 once.
    #[test]
    fn relaxation_transfer_is_the_law_and_causal() {
        let state = relaxing_air(true);
        let (n, dt, range, n0) = (1usize << 16, 1.0 / 48_000.0, 10.0, 1usize << 14);
        let filter = ObserverFilter::molecular_relaxation(&state, range, 48_000).unwrap();
        let receipt = *filter.receipt();
        let mut x = vec![0.0; n];
        x[n0] = 1.0;
        let delay = propagate_pressure_history_relaxing(&mut x, dt, range, &state).unwrap();
        assert_eq!(delay.to_bits(), receipt.receiver_delay_s.to_bits());
        assert!((delay + receipt.low_frequency_delay_s - range / state.gas().sound_speed).abs() < 1e-14);
        let total: f64 = x.iter().map(|v| v * v).sum();
        let early: f64 = x[..n0].iter().map(|v| v * v).sum();
        assert!(early <= 1e-20 * total, "pre-arrival energy ratio {:e}", early / total);
        let fft = Fft::new(n);
        let mut y: Vec<FftC64> = x.iter().map(|&v| FftC64::new(v, 0.0)).collect();
        let mut scratch = vec![FftC64::new(0.0, 0.0); n];
        fft.forward(&mut y, &mut scratch);
        let c_inf = state.frozen_sound_speed_m_s();
        let mut worst: f64 = 0.0;
        for k in (1..n / 2).step_by(997) {
            let omega = core::f64::consts::TAU * k as f64 / (n as f64 * dt);
            let wave = state.causal_wavenumber(omega).unwrap();
            let magnitude = (-wave.im * range).exp();
            let phase = -(wave.re - omega / c_inf) * range
                - omega * receipt.latency_s
                - core::f64::consts::TAU * (k * n0) as f64 / n as f64;
            let (er, ei) = (magnitude * phase.cos(), magnitude * phase.sin());
            worst = worst.max((y[k].re - er).hypot(y[k].im - ei));
        }
        assert!(worst <= 1e-3, "complex transfer error {worst:e}");
    }

    /// G3: finite records are zero-extended, never treated as periodic.
    #[test]
    fn relaxation_late_impulse_cannot_wrap_into_the_record_start() {
        let state = relaxing_air(true);
        let filter = ObserverFilter::molecular_relaxation(&state, 10.0, 48_000).unwrap();
        for n in [1, 2, 3, 17, 257, 2_048, 2_053] {
            let mut pressure = vec![0.0; n];
            pressure[n - 1] = 1.0;
            propagate_pressure_history_relaxing(&mut pressure, 1.0 / 48_000.0, 10.0, &state).unwrap();
            assert!(pressure[..n - 1].iter().all(|x| x.abs() < 1e-12), "n={n}");
            assert!((pressure[n - 1] - filter.taps()[0]).abs() < 1e-12, "n={n}");
        }
    }

    /// G3: the offline and streaming paths share one FIR, including across
    /// overlap-add block boundaries and for non-power-of-two record lengths.
    #[test]
    fn relaxation_batch_matches_streaming_and_preserves_an_exact_clock() {
        let state = relaxing_air(true);
        let filter = ObserverFilter::molecular_relaxation(&state, 10.0, 48_000).unwrap();
        let block = 2 * filter.taps().len().next_power_of_two() - filter.taps().len() + 1;
        for n in [0, 1, 3, 17, 2 * block + 3] {
            let mut streaming = filter.clone();
            let input: Vec<_> = (0..n).map(|i| (i % 29) as f64 / 14.0 - 1.0).collect();
            let expected: Vec<_> = input.iter().map(|&x| streaming.process(x)).collect();
            let mut batch = input;
            propagate_pressure_history_relaxing(&mut batch, 1.0 / 48_000.0, 10.0, &state).unwrap();
            assert!(batch.iter().zip(&expected).all(|(a, b)| (a - b).abs() < 1e-12), "n={n}");
        }

        // A fractional-Hz clock is supported without converting it to u32.
        let dt = 1.0 / 44_100.25;
        let response = molecular_response(&state, 10.0, 1.0 / dt).unwrap();
        let mut impulse = vec![0.0; response.taps.len() + 17];
        impulse[0] = 1.0;
        let delay = propagate_pressure_history_relaxing(&mut impulse, dt, 10.0, &state).unwrap();
        assert!(impulse.iter().zip(&response.taps).all(|(a, b)| (a - b).abs() < 1e-12));
        assert!(impulse[response.taps.len()..].iter().all(|x| x.abs() < 1e-12));
        assert!((delay + response.low_frequency_delay_s - 10.0 / state.gas().sound_speed).abs() < 1e-14);
    }

    /// G0: admission/design and both convolution paths commit only finite
    /// output; an overflowing intermediate cannot partly replace the input.
    #[test]
    fn relaxation_batch_refusal_is_atomic() {
        let state = relaxing_air(true);
        let original = [1.0, -2.0, 3.0];
        let mut pressure = original;
        assert!(propagate_pressure_history_relaxing(
            &mut pressure, f64::from_bits(1), 1.0, &state,
        ).is_err());
        assert_eq!(pressure, original);
        for n in [3, 65] {
            let mut pressure = vec![f64::MAX; n];
            let original = pressure.clone();
            assert!(convolve_pressure_history(&mut pressure, &vec![2.0; n]).is_err());
            assert_eq!(pressure, original);
        }
    }

    fn convention_db(omission: RelaxationOmission) -> f64 {
        let OmissionBound::StandardAirConvention { db } = omission.bound;
        db
    }

    /// The omission figure is labelled the standard-air convention: the worst
    /// standard-air ISO absorption at the band top times the range, linear in
    /// range, rising with the band, never below ISO at any grid state it
    /// covers, and refusing invalid inputs.
    #[test]
    fn relaxation_omission_is_the_standard_air_convention() {
        let one = convention_db(relaxation_omission(1.0, 24_000.0).unwrap());
        let two = convention_db(relaxation_omission(2.0, 24_000.0).unwrap());
        assert_eq!(two.to_bits(), (2.0 * one).to_bits());
        let lower = convention_db(relaxation_omission(1.0, 10_000.0).unwrap());
        assert!(lower < one);
        let omega = core::f64::consts::TAU * 24_000.0;
        for (t, rh) in [(253.15, 0.1), (293.15, 0.5), (303.15, 0.2), (323.15, 1.0)] {
            let iso = fs_material::iso9613::iso9613_absorption_db_per_km(t, 101_325.0, rh, omega).unwrap() * 1e-3;
            assert!(one >= iso, "convention {one} below ISO {iso} at {t} K, RH {rh}");
        }
        assert_eq!(convention_db(relaxation_omission(0.0, 24_000.0).unwrap()), 0.0);
        for (r, f) in [(-1.0, 1_000.0), (f64::NAN, 1_000.0), (1.0, 0.0), (1.0, f64::INFINITY)] {
            assert!(relaxation_omission(r, f).is_err());
        }
        println!(
            "{{\"bead\":\"M0.1-W2\",\"test\":\"relaxation_omission\",\"label\":\"standard-air convention\",\"band_hz\":24000,\"db_1m\":{one},\"db_2m\":{two},\"db_1m_10kHz\":{lower}}}"
        );
    }

    /// With no relaxing reservoir the causal transfer's loss is Stokes-Kirchhoff
    /// to relative O(w delta), and a non-dispersive gas adds only that small
    /// causal phase.
    #[test]
    fn classical_only_transfer_reduces_to_stokes_kirchhoff() {
        let state = relaxing_air(false);
        let gas = *state.gas();
        for f in [100.0, 1_000.0, 10_000.0, 20_000.0] {
            let omega = core::f64::consts::TAU * f;
            let k = state.causal_wavenumber(omega).unwrap();
            let rel = (k.im / gas.stokes_kirchhoff_absorption(omega) - 1.0).abs();
            assert!(rel <= 1e-6, "f={f}: causal loss off Stokes-Kirchhoff by {rel:e}");
        }
        let original = vec![1.0, -1.0, 0.5, -0.5, 1.0, 0.0, -1.0, 0.0];
        let mut signal = original.clone();
        assert!(propagate_pressure_history_relaxing(&mut signal, 0.0, 1.0, &state).is_err());
        assert!(propagate_pressure_history_relaxing(&mut signal, 1e-4, -1.0, &state).is_err());
        let mut poisoned = original.clone();
        poisoned[3] = f64::NAN;
        assert!(propagate_pressure_history_relaxing(&mut poisoned, 1e-4, 1.0, &state).is_err());
        assert_eq!(signal, original);
        propagate_pressure_history_relaxing(&mut signal, 1e-4, 0.0, &state).unwrap();
        assert_eq!(signal, original, "zero range is a no-op");
    }

    /// `sum_n taps[n] e^{-i w n dt}`: the FIR's own transfer at `f`.
    fn fir_transfer(taps: &[f64], rate: f64, f: f64) -> (f64, f64) {
        taps.iter().enumerate().fold((0.0, 0.0), |(re, im), (n, tap)| {
            let phase = -core::f64::consts::TAU * f * n as f64 / rate;
            (re + tap * phase.cos(), im + tap * phase.sin())
        })
    }

    /// The streaming relaxation filter is the causal transfer it was designed
    /// from, delayed by its declared latency: per frequency within 1e-3 (the
    /// scale of sqrt(dropped taps x 1e-9), declared before the run). Its
    /// low-frequency delay minus the latency is the dispersion
    /// `r (1/c_0 - 1/c_inf)`, and the receiver delay it names restores
    /// `r / c_0` exactly once.
    #[test]
    fn relaxation_filter_is_the_causal_transfer_and_names_its_delay_reference() {
        let state = relaxing_air(true);
        let gas = *state.gas();
        let (range, rate) = (10.0, 48_000u32);
        let filter = ObserverFilter::molecular_relaxation(&state, range, rate).unwrap();
        let receipt = *filter.receipt();
        assert!(receipt.truncated_energy <= OBSERVER_TRUNCATION_ENERGY);
        assert_eq!(receipt.absorption.law, GasAbsorptionModel::MolecularRelaxation);
        assert!((0.0..64.0 / f64::from(rate)).contains(&receipt.latency_s), "{}", receipt.latency_s);
        let c_inf = state.frozen_sound_speed_m_s();
        let mut worst: f64 = 0.0;
        for f in [100.0, 1_000.0, 5_000.0, 15_000.0, 23_000.0] {
            let w = core::f64::consts::TAU * f;
            let (re, im) = fir_transfer(filter.taps(), f64::from(rate), f);
            let (er, ei) = relaxing_transfer(&state, w, range, c_inf).unwrap();
            let (sin, cos) = (-w * receipt.latency_s).sin_cos();
            let (er, ei) = (er * cos - ei * sin, er * sin + ei * cos);
            worst = worst.max((re - er).hypot(im - ei));
        }
        assert!(worst <= 1e-3, "FIR off its transfer by {worst:e}");
        let dispersion = range * (1.0 / gas.sound_speed - 1.0 / c_inf);
        let physical = receipt.low_frequency_delay_s - receipt.latency_s;
        assert!(
            (physical / dispersion - 1.0).abs() <= 1e-3,
            "delay {physical} against r(1/c0 - 1/c_inf) = {dispersion}"
        );
        let total = range / receipt.delay_speed_m_s.unwrap() + receipt.low_frequency_delay_s;
        assert!((total * gas.sound_speed / range - 1.0).abs() <= 1e-12);
        println!(
            "{{\"bead\":\"pexfp.18\",\"test\":\"relaxation_filter\",\"taps\":{},\"truncated_energy\":{:e},\"latency_s\":{:e},\"dispersion_s\":{dispersion:e},\"delay_speed_m_s\":{},\"c0\":{},\"c_inf\":{c_inf},\"worst_transfer_error\":{worst:e}}}",
            receipt.taps, receipt.truncated_energy, receipt.latency_s,
            receipt.delay_speed_m_s.unwrap(), gas.sound_speed
        );
    }

    /// The ISO companion has ISO 9613-1's magnitude at every frequency (within
    /// 1e-3 absolute, declared before the run), and the receiver delay it names restores
    /// `r / c_0` once.
    #[test]
    fn iso_companion_filter_has_iso_magnitude() {
        let (gas, air) = TerrestrialAir::moist_air(293.15, 101_325.0, 0.5).unwrap();
        let (range, rate) = (100.0, 48_000u32);
        let filter = ObserverFilter::iso9613(&air, &gas, 0.5, range, rate).unwrap();
        let receipt = *filter.receipt();
        assert!(receipt.truncated_energy <= OBSERVER_TRUNCATION_ENERGY);
        assert_eq!(receipt.absorption.relaxation_omitted, None);
        assert_eq!(receipt.classical_validity_exceeded_from_hz, None);
        for f in [50.0, 500.0, 2_000.0, 8_000.0, 20_000.0] {
            let (re, im) = fir_transfer(filter.taps(), f64::from(rate), f);
            let alpha = fs_material::iso9613::iso9613_absorption_neper_per_m(
                293.15, 101_325.0, 0.5, core::f64::consts::TAU * f,
            )
            .unwrap();
            let expected = det::exp(-alpha * range);
            assert!(
                (re.hypot(im) - expected).abs() <= 1e-3,
                "f={f}: |H| {} against ISO {expected}",
                re.hypot(im)
            );
        }
        let total = range / receipt.delay_speed_m_s.unwrap() + receipt.low_frequency_delay_s;
        assert!((total * gas.sound_speed / range - 1.0).abs() <= 1e-12);
        assert!(receipt.low_frequency_delay_s > 0.0);
        // An admission does not transfer to another gas.
        let (other, _) = TerrestrialAir::moist_air(293.15, 101_325.0, 0.2).unwrap();
        assert!(matches!(
            ObserverFilter::iso9613(&air, &other, 0.2, range, rate),
            Err(AcousticRealizeError::NonTerrestrialAir(_))
        ));
        println!(
            "{{\"bead\":\"pexfp.18\",\"test\":\"iso_companion_filter\",\"range_m\":{range},\"taps\":{},\"truncated_energy\":{:e},\"low_frequency_delay_s\":{:e},\"latency_s\":{:e},\"delay_speed_m_s\":{}}}",
            receipt.taps, receipt.truncated_energy, receipt.low_frequency_delay_s, receipt.latency_s, receipt.delay_speed_m_s.unwrap()
        );
    }

    /// G1, ruling 6297, band declared before the run. The ISO companion's
    /// complex excess wavenumber (absorption AND dispersion) differs from the
    /// relaxing law's with ISO rates exactly by the landed strength + dilution
    /// decomposition, now applied to the complex quantity, to within 1e-3 of the
    /// companion: the same residual band as the landed absorption gate
    /// (fs-material `g1_exact_minus_iso_is_strength_plus_dilution`). Per process
    /// the two share the Debye form, so a wrong companion dispersion fails here.
    #[test]
    fn g1_iso_companion_is_the_relaxing_law_up_to_strength_and_dilution() {
        let law = fs_material::relaxation::RelaxationAbsorption::iso9613_terrestrial_air().unwrap();
        let (p, theta) = (101_325.0, [2239.1, 3352.0]);
        let iso_x = [0.209, 0.781];
        let mut worst: f64 = 0.0;
        for t in [273.15_f64, 293.15, 313.15] {
            for rh in [0.2, 0.5, 0.8] {
                let c_sat = -6.8346 * (273.16 / t).powf(1.261) + 4.6151;
                let x_w = rh * 10.0_f64.powf(c_sat) * 101_325.0 / p;
                let dry = 1.0 - x_w;
                let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), t, p).unwrap();
                let parts = [("N2", 0.781 * dry), ("O2", 0.209 * dry), ("Ar", 0.010 * dry), ("H2O", x_w)];
                let medium = RelaxingMedium::new(
                    gas,
                    parts.iter().map(|(name, x)| (SpeciesId::new(*name).unwrap(), *x)).collect(),
                )
                .unwrap();
                let state = law.evaluate(&medium).unwrap();
                let terms = fs_material::iso9613::iso9613_terms(t, p, rh).unwrap();
                let x = [0.209 * dry, 0.781 * dry];
                let r_s = gas.specific_gas_constant;
                let cv = r_s / (gas.gamma - 1.0);
                let c_j: Vec<f64> = theta
                    .iter()
                    .zip(x)
                    .map(|(&th, xj)| {
                        let u = th / t;
                        xj * r_s * u * u * u.exp() / (u.exp() - 1.0).powi(2)
                    })
                    .collect();
                let cvi = cv - c_j.iter().sum::<f64>();
                let cpi = cvi + r_s;
                let c_inf = (r_s * t * cpi / cvi).sqrt();
                let a = r_s / cpi;
                for i in 0..13 {
                    let f = 50.0 * 200.0_f64.powf(f64::from(i) / 12.0);
                    let w = core::f64::consts::TAU * f;
                    // Companion per process: strength f_r x / (1 - i x) = s f_r x (1 + i x) / (1 + x^2).
                    let iso_j: Vec<(f64, f64)> = terms
                        .relaxation
                        .iter()
                        .map(|(s, fr)| {
                            let xx = f / fr;
                            let d = 1.0 + xx * xx;
                            (s * fr * xx / d, s * fr * xx * xx / d)
                        })
                        .collect();
                    let iso = (iso_j[0].0 + iso_j[1].0, iso_j[0].1 + iso_j[1].1);
                    // u = sum_j C_j D_j / c_vinf with D = (1 + i w tau)/(1 + (w tau)^2).
                    let (mut u_re, mut u_im) = (0.0, 0.0);
                    for j in 0..2 {
                        let wt = f / terms.relaxation[j].1;
                        u_re += c_j[j] / cvi / (1.0 + wt * wt);
                        u_im += c_j[j] / cvi * wt / (1.0 + wt * wt);
                    }
                    let (u2_re, u2_im) = (u_re * u_re - u_im * u_im, 2.0 * u_re * u_im);
                    let q = 3.0 * a * a / 8.0 - a / 2.0;
                    let series = (w / c_inf * (a / 2.0 * u_re + q * u2_re), w / c_inf * (a / 2.0 * u_im + q * u2_im));
                    let undiluted = (0..2).fold((0.0, 0.0), |s, j| {
                        (s.0 + iso_j[j].0 * x[j] / iso_x[j], s.1 + iso_j[j].1 * x[j] / iso_x[j])
                    });
                    let dilution = (0..2).fold((0.0, 0.0), |s, j| {
                        let m = x[j] / iso_x[j] - 1.0;
                        (s.0 + iso_j[j].0 * m, s.1 + iso_j[j].1 * m)
                    });
                    let strength = (series.0 - undiluted.0, series.1 - undiluted.1);
                    let k = state.relaxing_wavenumber(w).unwrap();
                    let exact = (k.re - w / state.frozen_sound_speed_m_s(), k.im);
                    let residual = (
                        exact.0 - (iso.0 + strength.0 + dilution.0),
                        exact.1 - (iso.1 + strength.1 + dilution.1),
                    );
                    let relative = residual.0.hypot(residual.1) / iso.0.hypot(iso.1);
                    worst = worst.max(relative);
                    assert!(relative <= 1e-3, "T={t} RH={rh} f={f}: residual {relative:e}");
                }
            }
        }
        println!(
            "{{\"bead\":\"pexfp.18\",\"test\":\"g1_iso_companion_decomposition\",\"worst_residual\":{worst:e},\"band\":1e-3}}"
        );
    }

    /// Rulings 6324 and 6352: the certificate is derived, not asserted. Each
    /// receipt's budget is min(1e-9, 1e-3 x what ISO 9613-1's stated +-10 % on
    /// alpha is in the same response-energy measure), so the realized dropped
    /// energy sits at least three orders below the law's own uncertainty at
    /// every range. Prints d, D, the budget's two terms and the realized
    /// energy, and the smallest on-axis baffled range at 48 kHz whose own
    /// propagation delay absorbs the filter's (`r / c_0` samples against the
    /// filter's delay plus one): the floor the budget derives, stated.
    #[test]
    fn certificate_sits_three_orders_below_iso_alpha_uncertainty() {
        let (gas, air) = TerrestrialAir::moist_air(293.15, 101_325.0, 0.5).unwrap();
        let water = gas.water_mole_fraction;
        let dry = 1.0 - water;
        let parts = [("N2", 0.7808 * dry), ("O2", 0.2095 * dry), ("Ar", 0.0093 * dry), ("CO2", 0.0004 * dry), ("H2O", water)];
        let medium = RelaxingMedium::new(
            gas,
            parts.iter().map(|(name, x)| (SpeciesId::new(*name).unwrap(), *x)).collect(),
        )
        .unwrap();
        let state = fs_material::relaxation::RelaxationAbsorption::iso9613_terrestrial_air()
            .unwrap()
            .evaluate(&medium)
            .unwrap();
        let design = |law: &str, range: f64| match law {
            "iso9613" => ObserverFilter::iso9613(&air, &gas, 0.5, range, 48_000).unwrap(),
            _ => ObserverFilter::molecular_relaxation(&state, range, 48_000).unwrap(),
        };
        for range in [0.05, 0.5, 30.0] {
            for law in ["iso9613", "relaxation"] {
                let receipt = *design(law, range).receipt();
                let uncertainty_term = OBSERVER_UNCERTAINTY_MARGIN * receipt.alpha_uncertainty_energy;
                println!(
                    "{{\"bead\":\"pexfp.18\",\"test\":\"observer_certificate\",\"law\":\"{law}\",\"range_m\":{range},\"alignment_samples\":{},\"padding_samples\":{},\"taps\":{},\"budget_energy\":{:e},\"budget_terms\":[{:e},{uncertainty_term:e}],\"truncated_energy\":{:e},\"alpha_uncertainty_energy\":{:e}}}",
                    receipt.alignment_samples, receipt.padding_samples, receipt.taps, receipt.budget_energy,
                    OBSERVER_TRUNCATION_ENERGY, receipt.truncated_energy, receipt.alpha_uncertainty_energy
                );
                assert_eq!(receipt.budget_energy, OBSERVER_TRUNCATION_ENERGY.min(uncertainty_term));
                assert!(receipt.truncated_energy <= receipt.budget_energy);
                assert!(receipt.truncated_energy <= 1e-3 * receipt.alpha_uncertainty_energy);
            }
        }
        for law in ["iso9613", "relaxation"] {
            let smallest = (10..=100).map(|mm5| f64::from(mm5) * 0.005).find(|&range| {
                let receipt = *design(law, range).receipt();
                range / gas.sound_speed * 48_000.0 >= receipt.low_frequency_delay_s * 48_000.0 + 1.0
            });
            println!(
                "{{\"bead\":\"pexfp.18\",\"test\":\"observer_admitted_floor\",\"law\":\"{law}\",\"rate_hz\":48000,\"air\":\"293.15 K, RH 0.5\",\"smallest_admitted_on_axis_range_m\":{smallest:?},\"scan_m\":[0.05,0.5,0.005]}}"
            );
            assert!(smallest.is_some_and(|range| range > 0.05), "{law}: the 5 cm receiver must refuse and some range up to 50 cm must admit");
        }
    }

    /// Stokes-Kirchhoff is realized as its in-band minimum phase, carries the
    /// omission receipt, and prints where `omega tau_v << 1` fails: never in
    /// the audio band at 1 atm, below Nyquist at 100 Pa. `none` is the
    /// identity with the equilibrium delay.
    #[test]
    fn classical_and_lossless_filters_declare_validity_and_delay() {
        let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
        let classical = ObserverFilter::stokes_kirchhoff(&air, 30.0, 48_000).unwrap();
        let receipt = *classical.receipt();
        assert_eq!(receipt.absorption.law, GasAbsorptionModel::StokesKirchhoff);
        assert!(receipt.absorption.relaxation_omitted.is_some());
        assert_eq!(receipt.classical_validity_exceeded_from_hz, None);
        let (re, im) = fir_transfer(classical.taps(), 48_000.0, 10_000.0);
        let expected = det::exp(-air.stokes_kirchhoff_absorption(core::f64::consts::TAU * 10_000.0) * 30.0);
        assert!((re.hypot(im) - expected).abs() <= 1e-3);
        assert_eq!(receipt.latency_s, 0.0, "a minimum phase needs no latency");
        let thin = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 100.0).unwrap();
        let limit = ObserverFilter::stokes_kirchhoff(&thin, 0.1, 48_000)
            .unwrap()
            .receipt()
            .classical_validity_exceeded_from_hz;
        assert!(limit.is_some_and(|hz| hz < 24_000.0), "{limit:?}");
        let mut lossless = ObserverFilter::lossless(&air, 30.0, 48_000).unwrap();
        assert_eq!(lossless.taps(), &[1.0]);
        assert_eq!(lossless.receipt().delay_speed_m_s, Some(air.sound_speed));
        assert_eq!(lossless.receipt().absorption.law, GasAbsorptionModel::Lossless);
        assert!(lossless.receipt().absorption.relaxation_omitted.is_some());
        assert_eq!(lossless.process(0.25).to_bits(), 0.25_f64.to_bits());
        println!(
            "{{\"bead\":\"pexfp.18\",\"test\":\"classical_filter\",\"taps\":{},\"low_frequency_delay_s\":{:e},\"validity_limit_100Pa_hz\":{limit:?}}}",
            receipt.taps, receipt.low_frequency_delay_s
        );
    }

    fn record(text: &str) -> fs_material::gas_species::input::GasRecord {
        fs_material::gas_species::input::read_gas_record(&mut text.split_ascii_whitespace(), None, None)
            .expect("gas record")
    }

    /// Ruling 6510: the 1.0 dB tolerance at a 20 kHz band top spares a 1 m
    /// listener in standard air and refuses a 2 m one. The library's ISO
    /// figures are cross-checked against an independent transcription of ISO
    /// 9613-1's formula (Python, 2026-10-23, printed to three decimals).
    #[test]
    fn omission_tolerance_spares_one_metre_and_refuses_two_in_standard_air() {
        assert_eq!(observer_band_top_hz(None, 48_000), OBSERVER_BAND_TOP_CAP_HZ);
        assert_eq!(observer_band_top_hz(Some(4_800.0), 48_000), 4_800.0);
        assert_eq!(observer_band_top_hz(None, 32_000), 16_000.0);
        assert_eq!(observer_band_top_hz(Some(30_000.0), 96_000), OBSERVER_BAND_TOP_CAP_HZ);
        // dB per metre at 20 C, 1 atm: RH 20 %, 50 %, 80 %.
        let independent = [
            (8_000.0, [0.217, 0.105, 0.069]),
            (16_000.0, [0.435, 0.365, 0.252]),
            (20_000.0, [0.506, 0.524, 0.380]),
            (24_000.0, [0.565, 0.690, 0.527]),
        ];
        for (band_hz, row) in independent {
            for (rh, expected) in [0.2, 0.5, 0.8].into_iter().zip(row) {
                let gas = record(&format!("293.15 101325 {rh}"));
                let omission = omitted_observer_absorption(&gas, 1.0, band_hz).unwrap();
                let ObserverAbsorptionFigure::Iso9613 { relative_humidity, db } = omission.figure else {
                    panic!("a legacy record is ISO air: {omission:?}")
                };
                assert_eq!(relative_humidity, rh);
                let library = fs_material::iso9613::iso9613_absorption_db_per_km(
                    293.15, 101_325.0, rh, core::f64::consts::TAU * band_hz,
                ).unwrap() * 1e-3;
                assert!((db - library).abs() <= 1e-12 * library, "{db} vs {library}");
                assert!((db - expected).abs() <= 1e-3, "{band_hz} Hz RH {rh}: {db} vs independent {expected}");
            }
        }
        let standard = record("293.15 101325 0.5");
        let one = omitted_observer_absorption(&standard, 1.0, OBSERVER_BAND_TOP_CAP_HZ).unwrap();
        let two = omitted_observer_absorption(&standard, 2.0, OBSERVER_BAND_TOP_CAP_HZ).unwrap();
        println!(
            "{{\"bead\":\"pexfp.18\",\"test\":\"observer_omission_tolerance\",\"band_top_hz\":{OBSERVER_BAND_TOP_CAP_HZ},\"db_1m\":{},\"db_2m\":{},\"tolerance_db\":{OBSERVER_OMISSION_TOLERANCE_DB}}}",
            one.omitted_db(false).unwrap(), two.omitted_db(false).unwrap()
        );
        assert!(one.omitted_db(false).unwrap() < OBSERVER_OMISSION_TOLERANCE_DB);
        assert!(two.omitted_db(false).unwrap() > OBSERVER_OMISSION_TOLERANCE_DB);
        assert!(one.admit(false, false).is_ok());
        let refused = two.admit(false, false).unwrap_err();
        assert_eq!(refused, AcousticRealizeError::ObserverOmission(Box::new(two.clone())));
        assert!(refused.to_string().contains("observer_law"), "{refused}");
        assert!(two.admit(false, true).is_ok(), "a declared law accepts the omission");
    }

    /// Rulings 6510 and 6532: the figure follows the recorded composition.
    /// Terrestrial air is ISO; monatomic gas is Stokes-Kirchhoff, complete;
    /// N2 and O2 off ISO's fractions take ISO's rates held, which at ISO's own
    /// fractions IS ISO's rate law; heliox and a CO2-rich gas are unbounded,
    /// naming the species, and refuse unless a law is declared.
    #[test]
    fn omission_figure_follows_the_recorded_composition() {
        use fs_material::relaxation::{RelaxationAbsorption, RelaxingMedium};
        let band = OBSERVER_BAND_TOP_CAP_HZ;
        let terrestrial = record("mixture 293.15 101325 mole 3 N2 0.7809 O2 0.2095 Ar 0.0096");
        let figure = omitted_observer_absorption(&terrestrial, 1.0, band).unwrap();
        assert!(matches!(figure.figure, ObserverAbsorptionFigure::Iso9613 { .. }), "{figure:?}");
        // The held law is ISO's rate law at ISO's fractions.
        let fs_material::gas_species::input::GasRecord::Mixture(state, composition) = &terrestrial else {
            panic!("mixture record")
        };
        let medium = RelaxingMedium::new(*state, composition.iter().map(|(species, x)|
            (fs_thermochem::SpeciesId::new(species.symbol()).unwrap(), *x)).collect()).unwrap();
        let held = RelaxationAbsorption::evaluate_iso9613_rates_held(&medium).unwrap();
        let iso = RelaxationAbsorption::iso9613_terrestrial_air().unwrap().evaluate(&medium).unwrap();
        for omega in [1e3, 1e4, 1e5] {
            let (a, b) = (held.causal_wavenumber(omega).unwrap(), iso.causal_wavenumber(omega).unwrap());
            assert!((a.im - b.im).abs() <= 1e-12 * b.im && (a.re - b.re).abs() <= 1e-12 * b.re, "{a:?} vs {b:?}");
        }

        let enriched = omitted_observer_absorption(&record("mixture 293.15 100000 mole 2 N2 0.7 O2 0.3"), 1.0, band).unwrap();
        let ObserverAbsorptionFigure::HeldIsoRates { o2, n2, water, db } = enriched.figure else {
            panic!("N2/O2 0.7/0.3 holds ISO's rates: {enriched:?}")
        };
        assert_eq!((o2, n2, water), (0.3, 0.7, 0.0));
        assert!(db > enriched.classical_db, "relaxation adds to the classical part");
        assert!(enriched.admit(true, false).is_ok() && enriched.admit(false, false).is_ok());

        let helium = omitted_observer_absorption(&record("mixture 293.15 100000 mole 1 He 1"), 1.0, band).unwrap();
        assert_eq!(helium.figure, ObserverAbsorptionFigure::ClassicalComplete { db: helium.classical_db });
        assert_eq!(helium.omitted_db(true), Some(0.0), "Stokes-Kirchhoff is helium's whole law");

        for (text, species) in [
            ("mixture 293.15 100000 mole 2 He 0.8 O2 0.2", "He"),
            ("mixture 293.15 100000 mole 3 N2 0.4 O2 0.2 CO2 0.4", "CO2"),
        ] {
            let omission = omitted_observer_absorption(&record(text), 0.5, band).unwrap();
            assert_eq!(omission.figure, ObserverAbsorptionFigure::Unbounded { species: species.to_owned() });
            assert_eq!(omission.omitted_db(true), None);
            let refused = omission.admit(true, false).unwrap_err().to_string();
            assert!(refused.contains(species) && refused.contains("observer_law"), "{refused}");
            assert!(omission.admit(true, true).is_ok(), "a declared law accepts the omission");
            assert!(omission.to_json().contains(&format!("\"kind\":\"unbounded\",\"species\":\"{species}\"")));
            println!("{{\"bead\":\"pexfp.18\",\"test\":\"observer_omission_unbounded\",\"gas\":\"{text}\",\"receipt\":{}}}", omission.to_json());
        }

        let unrecorded = ObserverOmission::unrecorded(state, 0.5, band);
        assert_eq!(unrecorded.omitted_db(true), None);
        assert!(unrecorded.admit(true, false).unwrap_err().to_string().contains("no record naming its composition"));
        assert!(unrecorded.admit(true, true).is_ok() && unrecorded.to_json().contains("\"kind\":\"unrecorded\""));

        let foreign = RelaxingMedium::new(*state, vec![
            (fs_thermochem::SpeciesId::new("N2").unwrap(), 0.9),
            (fs_thermochem::SpeciesId::new("SF6").unwrap(), 0.1),
        ]).unwrap();
        let refused = RelaxationAbsorption::evaluate_iso9613_rates_held(&foreign).unwrap_err().to_string();
        assert!(refused.contains("SF6"), "{refused}");
    }
}
