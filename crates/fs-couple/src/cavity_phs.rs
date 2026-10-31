//! Time-domain plate × lumped Helmholtz cavity.
//!
//! The plate is an `fs-phs` modal bank. The cavity is the flow-driven
//! Helmholtz pHS. They join through a `transformer` (area already in
//! the plate's flow port). A bottle and a vented enclosure are the
//! same objects. There is no guitar type.

/// Finite source gas and the pressure modes of a material-derived elastic wall.
pub mod modal;

/// Nonlinear closed-gas cavity coupled to an inertial moving wall.
pub mod piston;

use crate::acoustic_realize::AcousticRealizeError;
use crate::air_path::{GasAbsorptionModel, propagate_pressure_history_with_model};
use crate::thin_plate::certified_radiators;
use fs_material::gas::GasState;
use fs_material::relaxation::RelaxationState;
use fs_phs::linear_exact::{LinearPropagator, LinearSplitBudget, split_quadratic, step_split};
use fs_phs::{PortHamiltonian, QuadraticStorage, helmholtz_resonator_flow, transformer};

/// Typed refusal of the plate–cavity clock.
#[derive(Debug, Clone, PartialEq)]
pub enum CavityPhsError {
    /// Non-physical description.
    Invalid {
        /// Which check failed.
        what: &'static str,
    },
    /// pHS admission or step refusal.
    Phs(String),
    /// Material, geometry or modal admission refused by the plate owner.
    Plate(AcousticRealizeError),
}

impl core::fmt::Display for CavityPhsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid { what } => write!(f, "FS-COUPLE-CAVITY: {what}"),
            Self::Phs(e) => write!(f, "FS-COUPLE-CAVITY-PHS: {e}"),
            Self::Plate(e) => write!(f, "FS-COUPLE-CAVITY-PLATE: {e}"),
        }
    }
}

impl std::error::Error for CavityPhsError {}

/// A modal plate facing a lumped Helmholtz volume.
#[derive(Debug, Clone)]
pub struct PlateCavitySpec {
    /// Plate modal frequencies [rad/s].
    pub omegas: Vec<f64>,
    /// Complete modal damping ratios, including any resolved material loss.
    /// The cavity never adds an inferred material contribution.
    pub zetas: Vec<f64>,
    /// Drive weights at the external force (mass-normalized).
    pub drive: Vec<f64>,
    /// Mass-normalized monopole areas [m²/√kg] (volume velocity = A · v).
    pub areas: Vec<f64>,
    /// Cavity volume [m³].
    pub volume_m3: f64,
    /// Neck radius [m].
    pub neck_radius_m: f64,
    /// Neck length [m] (end corrections are applied inside the pHS).
    pub neck_length_m: f64,
    /// Baffle condition where the neck opens into the cavity.
    pub inner_end: fs_scenario::NeckEnd,
    /// Baffle condition where the neck opens into the exterior.
    pub outer_end: fs_scenario::NeckEnd,
    /// Admitted ambient state used consistently by the cavity, radiation,
    /// observer and propagation path. Plate temperature is independent.
    pub gas: GasState,
    /// Relative humidity in `[0, 1]`, used by [`realize_plate_cavity_in_air`].
    /// The generic-gas path validates this legacy metadata but does not use it
    /// to infer composition or select an air-only absorption law.
    pub relative_humidity: f64,
}

impl PlateCavitySpec {
    /// Replace all modal mechanics and loss coefficients from the ordinary
    /// thin-plate owner, including material-bound thermal/viscous damping and
    /// prestress dilution. Both force and volume-flow ports use the same
    /// mass-normalized modal basis. Cavity geometry and ambient gas stay explicit.
    ///
    /// A source-bound caller passes `ResolvedPlateSpecimen::plate()` and retains
    /// that specimen's receipts. This numeric adapter adds no source authority,
    /// thermal evolution, or exterior plate-radiation back-reaction.
    ///
    /// # Errors
    /// Nonlinear plates require a nonlinear owner. Other plate admission errors
    /// propagate. Any refusal leaves all previously bound coefficients unchanged.
    pub fn bind_plate(&mut self, plate: fs_scenario::ThinPlate) -> Result<(), CavityPhsError> {
        if plate.geometric_nonlinearity {
            return Err(CavityPhsError::Invalid {
                what: "the linear plate-cavity modal bank cannot bind a nonlinear plate",
            });
        }
        let modes = certified_radiators(plate).map_err(CavityPhsError::Plate)?;
        self.omegas = modes.iter().map(|mode| mode.omega).collect();
        self.zetas = modes.iter().map(|mode| mode.zeta).collect();
        self.drive = modes
            .iter()
            .map(|mode| mode.drive_participation / mode.mass_kg.sqrt())
            .collect();
        self.areas = modes
            .iter()
            .map(|mode| mode.area_m2 / mode.mass_kg.sqrt())
            .collect();
        Ok(())
    }
}

/// Realize observer pressure of a driven plate on a Helmholtz volume.
///
/// `force_n[i]` is the external force at sample `i`. The far-field
/// observer is the baffled compact sum `ρ A ÿ / (2π r)`.
/// The supplied gas drives cavity storage, radiation and classical
/// causal Stokes-Kirchhoff propagation. Terrestrial-air ISO absorption requires
/// [`realize_plate_cavity_in_air`]; an evaluated molecular relaxation law can
/// be supplied through [`realize_plate_cavity_with_observer`].
///
/// # Errors
/// Empty/mismatched modes, incoherent gas, non-physical cavity, or a pHS/
/// absorption refusal. Generic steam is not subjected to humid-air fit limits.
pub fn realize_plate_cavity(
    spec: &PlateCavitySpec,
    force_n: &[f64],
    sample_rate_hz: u32,
    listener_m: f64,
) -> Result<Vec<f64>, CavityPhsError> {
    realize_plate_cavity_with_observer(
        spec,
        force_n,
        sample_rate_hz,
        listener_m,
        GasAbsorptionModel::StokesKirchhoff,
        None,
    )
}

/// Explicit terrestrial-air realization using the spec's relative humidity.
/// `air` admits the spec's gas state as terrestrial air from its own
/// composition. Gas identity is never inferred from coefficients. Cavity
/// mechanics still use the same supplied gas state.
///
/// # Errors
/// All [`realize_plate_cavity`] refusals, an admission not made for
/// `spec.gas` (`AcousticRealizeError::NonTerrestrialAir`), and out-of-domain
/// ISO conditions. An unavailable air law is never silently replaced with a
/// different law.
pub fn realize_plate_cavity_in_air(
    spec: &PlateCavitySpec,
    air: &fs_material::iso9613::TerrestrialAir,
    force_n: &[f64],
    sample_rate_hz: u32,
    listener_m: f64,
) -> Result<Vec<f64>, CavityPhsError> {
    realize_plate_cavity_with_observer(
        spec,
        force_n,
        sample_rate_hz,
        listener_m,
        GasAbsorptionModel::Iso9613Air {
            relative_humidity: spec.relative_humidity,
            air: *air,
        },
        None,
    )
}

/// Realize the same plate/cavity mechanics with an explicit exterior gas law.
///
/// Pressure uses the shared causal complex transfer and linear convolution:
/// molecular absorption and dispersion come from one evaluated wavenumber,
/// ISO air keeps its Debye phase, and classical loss uses its declared in-band
/// minimum-phase realization. `Lossless` preserves the unfiltered pressure.
/// A molecular state must belong to `spec.gas` and is accepted only with
/// `MolecularRelaxation`; no air composition or relaxation is inferred.
///
/// The returned history follows source time and includes FIR implementation
/// latency, as the other offline acoustic realizers do. The remaining receiver
/// delay is not applied here. The finite cavity retains its existing equilibrium
/// lumped gas model; the selected law governs exterior observer propagation.
///
/// # Errors
/// All [`realize_plate_cavity`] refusals, missing/mismatched molecular data,
/// a nonmatching terrestrial-air admission, and the shared response designer's
/// finite-band/FIR-length refusals. No partial pressure history is returned.
pub fn realize_plate_cavity_with_observer(
    spec: &PlateCavitySpec,
    force_n: &[f64],
    sample_rate_hz: u32,
    listener_m: f64,
    absorption: GasAbsorptionModel,
    relaxation: Option<&RelaxationState>,
) -> Result<Vec<f64>, CavityPhsError> {
    let n_m = spec.omegas.len();
    if n_m == 0 || spec.zetas.len() != n_m || spec.drive.len() != n_m || spec.areas.len() != n_m {
        return Err(CavityPhsError::Invalid {
            what: "plate modal vectors must be non-empty and aligned",
        });
    }
    if !(sample_rate_hz > 0 && listener_m > 0.0 && listener_m.is_finite()) {
        return Err(CavityPhsError::Invalid {
            what: "sample rate and listener distance must be positive",
        });
    }
    if spec.gas.validate_coherence().is_err() {
        return Err(CavityPhsError::Invalid {
            what: "gas state must remain inside its admitted physical domain",
        });
    }
    if !(spec.relative_humidity.is_finite() && (0.0..=1.0).contains(&spec.relative_humidity)) {
        return Err(CavityPhsError::Invalid {
            what: "relative humidity must be finite and inside [0, 1]",
        });
    }
    let plate = plate_force_and_flow(&spec.omegas, &spec.zetas, &spec.drive, &spec.areas)
        .map_err(|e| CavityPhsError::Phs(e.to_string()))?;
    // One neck owner (DEF-09): declared ends, outer radiation + Bessel loss.
    let geometry = fs_scenario::HelmholtzCavity {
        volume_m3: spec.volume_m3,
        neck_radius_m: spec.neck_radius_m,
        neck_length_m: spec.neck_length_m,
        inner_end: spec.inner_end,
        outer_end: spec.outer_end,
    };
    let resistance = crate::cavity_modal::resonance_neck_resistance(&geometry, &spec.gas)
        .map_err(CavityPhsError::Plate)?;
    let cavity = helmholtz_resonator_flow(
        spec.volume_m3,
        spec.neck_radius_m,
        spec.neck_length_m,
        spec.gas.density,
        spec.gas.sound_speed,
        resistance,
        crate::cavity_modal::neck_ends(&geometry),
    )
    .map_err(|e| CavityPhsError::Phs(e.to_string()))?;
    // Area is already in the plate's second port. The transformer is
    // the power-conserving F = p A, U = A v join — not a staggered
    // clock.
    let sys =
        transformer(plate, cavity, 1, 0, 1.0).map_err(|e| CavityPhsError::Phs(e.to_string()))?;
    // Prepare the complete coupled operator, not two independently exact
    // components connected by a trapezoidal reaction.
    let linear = split_quadratic(
        &sys,
        LinearSplitBudget {
            max_modes: sys.state_dim(),
            max_omega_dt: core::f64::consts::PI,
        },
    )
    .map_err(|e| CavityPhsError::Phs(e.to_string()))?;
    let dt = 1.0 / f64::from(sample_rate_hz);
    let mut x = vec![0.0; sys.state_dim()];
    let mut out = Vec::with_capacity(force_n.len());
    let two_pi = 2.0 * core::f64::consts::PI;
    let n_plate = 2 * n_m;
    for &f_ext in force_n {
        let rec = step_split(&sys, &linear, LinearPropagator::ExactZoh, &x, &[f_ext], dt)
            .map_err(|e| CavityPhsError::Phs(e.to_string()))?;
        let mut p_obs = 0.0;
        for k in 0..n_m {
            let acc = (rec.x[2 * k + 1] - x[2 * k + 1]) / dt;
            p_obs += spec.gas.density * spec.areas[k] * acc / (two_pi * listener_m);
        }
        debug_assert!(rec.x.len() >= n_plate);
        x = rec.x;
        if !p_obs.is_finite() {
            return Err(CavityPhsError::Invalid {
                what: "observer pressure left the finite set",
            });
        }
        out.push(p_obs);
    }
    propagate_pressure_history_with_model(
        &mut out, dt, listener_m, &spec.gas, absorption, relaxation,
    )
    .map_err(CavityPhsError::Plate)?;
    Ok(out)
}

fn plate_force_and_flow(
    omegas: &[f64],
    zetas: &[f64],
    drive: &[f64],
    areas: &[f64],
) -> Result<PortHamiltonian, fs_phs::PhsError> {
    let nm = omegas.len();
    let n = 2 * nm;
    let mut q = vec![0.0; n * n];
    let mut j = vec![0.0; n * n];
    let mut r = vec![0.0; n * n];
    let mut g = vec![0.0; n * 2];
    for i in 0..nm {
        if omegas[i] <= 0.0 || zetas[i] < 0.0 {
            return Err(fs_phs::PhsError::NotPsd {
                what: "modal parameters",
            });
        }
        let (qi, pi) = (2 * i, 2 * i + 1);
        q[qi * n + qi] = omegas[i] * omegas[i];
        q[pi * n + pi] = 1.0;
        j[qi * n + pi] = 1.0;
        j[pi * n + qi] = -1.0;
        r[pi * n + pi] = 2.0 * zetas[i] * omegas[i];
        g[pi * 2] = drive[i];
        g[pi * 2 + 1] = areas[i];
    }
    let storage = Box::new(QuadraticStorage::new(q, n)?);
    PortHamiltonian::new(n, 2, j, r, g, storage)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::gas::GasSpec;

    fn air_at(temperature_k: f64) -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature_k, 101_325.0)
            .expect("admitted dry air")
    }

    #[test]
    fn full_plate_cavity_transformer_has_unwarped_coupled_eigenfrequencies() {
        use fs_math::det;
        let gas = air_at(293.15);
        let plate = plate_force_and_flow(&[6_000.0], &[0.0], &[1.0], &[0.08]).unwrap();
        let ends = fs_phs::NeckEnds {
            inner: fs_phs::MouthFlange::Flanged,
            outer: fs_phs::MouthFlange::Flanged,
        };
        let cavity =
            helmholtz_resonator_flow(0.02, 0.02, 0.06, gas.density, gas.sound_speed, 0.0, ends)
                .unwrap();
        let sys = transformer(plate, cavity, 1, 0, 1.0).unwrap();
        let split = split_quadratic(
            &sys,
            LinearSplitBudget {
                max_modes: 4,
                max_omega_dt: core::f64::consts::PI,
            },
        )
        .unwrap();
        let frequencies = split.angular_frequencies_rad_s();
        let (low, high) = (frequencies[0], frequencies[3]);
        assert!(low > 0.0 && high > low);
        for z in [0.1, 1.0, 2.5] {
            let dt = z / high;
            let mut map = [[0.0; 4]; 4];
            for column in 0..4 {
                let mut initial = [0.0; 4];
                initial[column] = 1.0;
                let record = step_split(
                    &sys,
                    &split,
                    LinearPropagator::ExactZoh,
                    &initial,
                    &[0.0],
                    dt,
                )
                .unwrap();
                for row in 0..4 {
                    map[row][column] = record.x[row];
                }
            }
            // A lossless four-state map has poles exp(+-i theta_0),
            // exp(+-i theta_1). Recover both cosines from trace E and trace E².
            let trace: f64 = (0..4).map(|i| map[i][i]).sum();
            let square_trace: f64 = (0..4)
                .map(|i| (0..4).map(|j| map[i][j] * map[j][i]).sum::<f64>())
                .sum();
            let sum = 0.5 * trace;
            let sum_squares = 0.25 * (square_trace + 4.0);
            let difference = det::sqrt(2.0 * sum_squares - sum * sum);
            let measured = [
                det::acos(0.5 * (sum + difference)) / dt,
                det::acos(0.5 * (sum - difference)) / dt,
            ];
            for (observed, expected) in measured.into_iter().zip([low, high]) {
                let cents = 1200.0 * det::ln(observed / expected) / core::f64::consts::LN_2;
                assert!(
                    cents.abs() <= 0.01,
                    "coupled plate/cavity defect: {cents} cents at omega*dt={z}"
                );
            }
        }
    }

    fn panel(volume: f64) -> PlateCavitySpec {
        PlateCavitySpec {
            omegas: vec![6.0e2],
            zetas: vec![0.006],
            drive: vec![1.0],
            areas: vec![0.08],
            volume_m3: volume,
            neck_radius_m: 0.02,
            neck_length_m: 0.03,
            // A port through the panel's own baffle: flanged on both sides.
            inner_end: fs_scenario::NeckEnd::Flanged,
            outer_end: fs_scenario::NeckEnd::Flanged,
            gas: air_at(293.15),
            relative_humidity: 0.0,
        }
    }

    fn pulse(n: usize, peak: f64) -> Vec<f64> {
        (0..n)
            .map(|i| {
                let t = i as f64 / 8_000.0;
                if t < 0.002 {
                    peak * (core::f64::consts::PI * t / 0.002).sin()
                } else {
                    0.0
                }
            })
            .collect()
    }

    fn period(x: &[f64]) -> f64 {
        let mut prev = x[x.len() / 4];
        let mut times = Vec::new();
        for (i, &s) in x.iter().enumerate().skip(x.len() / 4 + 1) {
            if prev > 0.0 && s <= 0.0 {
                times.push(i as f64);
            }
            prev = s;
        }
        assert!(times.len() >= 3, "need crossings, got {}", times.len());
        (times[times.len() - 1] - times[0]) / (times.len() - 1) as f64
    }

    #[test]
    fn larger_cavity_lowers_the_coupled_frequency() {
        let f = pulse(2_000, 8.0);
        let small = realize_plate_cavity(&panel(0.004), &f, 8_000, 1.0).expect("small V");
        let large = realize_plate_cavity(&panel(0.016), &f, 8_000, 1.0).expect("large V");
        let ts = period(&small);
        let tl = period(&large);
        assert!(
            tl > ts * 1.02,
            "Helmholtz ω ~ 1/√V must lengthen the period ({tl:.2} vs {ts:.2})"
        );
    }

    #[test]
    fn declared_modal_loss_changes_the_waveform() {
        let f = pulse(800, 2.0);
        let mut with = panel(0.008);
        let mut bare = with.clone();
        bare.zetas[0] = 0.0;
        with.zetas[0] = 0.02;
        let a = realize_plate_cavity(&bare, &f, 8_000, 1.0).expect("bare");
        let b = realize_plate_cavity(&with, &f, 8_000, 1.0).expect("te");
        let err: f64 = a.iter().zip(&b).map(|(x, y)| (x - y).abs()).sum();
        assert!(err > 1.0e-10, "declared modal loss must move the waveform");
    }

    #[test]
    fn admitted_gas_state_sets_the_cavity_period_without_a_hidden_air_default() {
        let f = pulse(2_000, 8.0);
        let mut cold = panel(0.008);
        cold.gas = air_at(273.15);
        let mut warm = cold.clone();
        warm.gas = air_at(323.15);

        let cold_period =
            period(&realize_plate_cavity(&cold, &f, 8_000, 1.0).expect("cold cavity realizes"));
        let warm_period =
            period(&realize_plate_cavity(&warm, &f, 8_000, 1.0).expect("warm cavity realizes"));
        assert!(
            warm_period < cold_period * 0.99,
            "warmer admitted air must raise c and shorten the cavity period ({warm_period:.2} vs {cold_period:.2})"
        );
    }

    #[test]
    fn forged_gas_state_is_refused_before_cavity_compute() {
        let mut spec = panel(0.008);
        let admitted = spec.gas;
        let mut cases: Vec<(&str, GasState)> = Vec::new();
        for (name, mutate) in [
            (
                "temperature",
                (|gas: &mut GasState| gas.temperature = f64::NAN) as fn(&mut GasState),
            ),
            ("pressure", |gas: &mut GasState| gas.pressure = 0.0),
            ("density", |gas: &mut GasState| gas.density = f64::NAN),
            ("sound speed", |gas: &mut GasState| gas.sound_speed = 0.0),
            ("dynamic viscosity", |gas: &mut GasState| {
                gas.dynamic_viscosity = f64::NAN
            }),
            ("thermal conductivity", |gas: &mut GasState| {
                gas.thermal_conductivity = 0.0
            }),
            ("gamma", |gas: &mut GasState| gas.gamma = f64::NAN),
            ("gamma ceiling", |gas: &mut GasState| {
                gas.gamma = 5.0 / 3.0 + 1.0e-9
            }),
            ("specific gas constant", |gas: &mut GasState| {
                gas.specific_gas_constant = 0.0
            }),
            ("specific heat", |gas: &mut GasState| {
                gas.specific_heat_cp = f64::NAN
            }),
            ("Prandtl number", |gas: &mut GasState| gas.prandtl = 0.0),
            ("characteristic impedance", |gas: &mut GasState| {
                gas.characteristic_impedance = f64::NAN
            }),
            ("water fraction", |gas: &mut GasState| {
                gas.water_mole_fraction = -1.0
            }),
            ("incoherent moist temperature", |gas: &mut GasState| {
                gas.temperature = 200.0;
                gas.water_mole_fraction = 0.01;
            }),
            ("finite inconsistent pressure", |gas: &mut GasState| {
                gas.pressure *= 2.0
            }),
            ("finite inconsistent impedance", |gas: &mut GasState| {
                gas.characteristic_impedance *= 2.0
            }),
        ] {
            let mut forged = admitted;
            mutate(&mut forged);
            cases.push((name, forged));
        }

        for (name, forged) in cases {
            spec.gas = forged;
            let result = realize_plate_cavity(&spec, &[1.0], 8_000, 1.0);
            assert!(
                matches!(result, Err(CavityPhsError::Invalid { what }) if what == "gas state must remain inside its admitted physical domain"),
                "forged {name} must be refused, got {result:?}"
            );
        }

        spec.gas = admitted;
        spec.relative_humidity = f64::NAN;
        assert!(matches!(
            realize_plate_cavity(&spec, &[1.0], 8_000, 1.0),
            Err(CavityPhsError::Invalid { what })
                if what == "relative humidity must be finite and inside [0, 1]"
        ));

        let (moist, air) = fs_material::iso9613::TerrestrialAir::moist_air(293.15, 101_325.0, 0.5)
            .expect("admitted moist air");
        spec.gas = moist;
        spec.relative_humidity = 0.5;
        assert!(
            realize_plate_cavity(&spec, &[1.0], 8_000, 1.0).is_ok(),
            "the shared guard must preserve an admitted moist state"
        );
        assert!(realize_plate_cavity_in_air(&spec, &air, &[1.0], 8_000, 1.0).is_ok());
    }

    /// The actual cavity renderer must preserve phase and cannot produce
    /// pressure before a late driving force through periodic FFT wraparound.
    #[test]
    fn cavity_observers_share_streaming_dispersion_without_wraparound() {
        use crate::air_path::ObserverFilter;
        let (gas, air) =
            fs_material::iso9613::TerrestrialAir::moist_air(293.15, 101_325.0, 0.5).unwrap();
        let mut spec = panel(0.008);
        spec.gas = gas;
        spec.relative_humidity = 0.5;
        let mut input = vec![0.0; 513];
        input[400] = 1.0;
        let raw = realize_plate_cavity_with_observer(
            &spec,
            &input,
            8_000,
            30.0,
            GasAbsorptionModel::Lossless,
            None,
        )
        .unwrap();
        let peak = raw
            .iter()
            .fold(0.0_f64, |peak, value| peak.max(value.abs()));
        assert!(peak > 0.0);
        for law in [
            GasAbsorptionModel::StokesKirchhoff,
            GasAbsorptionModel::Iso9613Air {
                relative_humidity: 0.5,
                air,
            },
        ] {
            let actual = match law {
                GasAbsorptionModel::StokesKirchhoff => {
                    realize_plate_cavity(&spec, &input, 8_000, 30.0).unwrap()
                }
                _ => realize_plate_cavity_in_air(&spec, &air, &input, 8_000, 30.0).unwrap(),
            };
            let mut observer = ObserverFilter::for_law(law, &gas, None, 30.0, 8_000).unwrap();
            let expected: Vec<_> = raw.iter().map(|&value| observer.process(value)).collect();
            assert!(
                actual
                    .iter()
                    .zip(expected)
                    .all(|(a, b)| (a - b).abs() < 1e-12 * peak)
            );
            assert!(actual[..400].iter().all(|value| value.abs() < 1e-12 * peak));
            assert!(actual[400..].iter().any(|value| value.abs() > 1e-6 * peak));
        }
    }

    #[test]
    fn cavity_molecular_observer_uses_the_matching_evaluated_gas() {
        use crate::air_path::ObserverFilter;
        use fs_material::relaxation::{
            InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
            RelaxingMedium, SpeciesRelaxation, VibrationalMode,
        };
        let spec = panel(0.008);
        // Synthetic reservoir tests the physical consumer; it is not air data.
        let species = fs_thermochem::SpeciesId::new("cavity-test-gas").unwrap();
        let law = RelaxationAbsorption::new(vec![SpeciesRelaxation {
            species: species.clone(),
            processes: vec![RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 960.0,
                    degeneracy: 2,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 1e-4 * spec.gas.pressure,
                    t_ref_k: 293.15,
                    exponent: 0.0,
                    t_min_k: 100.0,
                    t_max_k: 2000.0,
                },
            }],
        }])
        .unwrap();
        let relaxation = law
            .evaluate(&RelaxingMedium::new(spec.gas, vec![(species, 1.0)]).unwrap())
            .unwrap();
        let mut input = vec![0.0; 513];
        input[400] = 1.0;
        let raw = realize_plate_cavity_with_observer(
            &spec,
            &input,
            8_000,
            30.0,
            GasAbsorptionModel::Lossless,
            None,
        )
        .unwrap();
        let actual = realize_plate_cavity_with_observer(
            &spec,
            &input,
            8_000,
            30.0,
            GasAbsorptionModel::MolecularRelaxation,
            Some(&relaxation),
        )
        .unwrap();
        let mut observer = ObserverFilter::molecular_relaxation(&relaxation, 30.0, 8_000).unwrap();
        let expected: Vec<_> = raw.iter().map(|&value| observer.process(value)).collect();
        let peak = raw
            .iter()
            .fold(0.0_f64, |peak, value| peak.max(value.abs()));
        assert!(peak > 0.0);
        assert!(
            actual
                .iter()
                .zip(expected)
                .all(|(a, b)| (a - b).abs() < 1e-12 * peak)
        );
        assert!(actual[..400].iter().all(|value| value.abs() < 1e-12 * peak));
        assert!(
            actual
                .iter()
                .zip(&raw)
                .any(|(a, b)| (a - b).abs() > 1e-4 * peak)
        );
        assert!(
            realize_plate_cavity_with_observer(
                &spec,
                &[],
                8_000,
                30.0,
                GasAbsorptionModel::MolecularRelaxation,
                None,
            )
            .is_err()
        );
        assert!(
            realize_plate_cavity_with_observer(
                &spec,
                &[],
                8_000,
                30.0,
                GasAbsorptionModel::StokesKirchhoff,
                Some(&relaxation),
            )
            .is_err()
        );
        let mut other = spec;
        other.gas = air_at(323.15);
        assert!(
            realize_plate_cavity_with_observer(
                &other,
                &[],
                8_000,
                30.0,
                GasAbsorptionModel::MolecularRelaxation,
                Some(&relaxation),
            )
            .is_err()
        );
    }

    #[test]
    fn steam_and_air_share_cavity_mechanics_without_air_only_limits() {
        use fs_qty::semantic::{Composition, CompositionBasis};
        let mut spec = panel(0.008);
        let input = pulse(64, 1.0);
        let air = realize_plate_cavity(&spec, &input, 8000, 1.0).unwrap();
        let composition = Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap();
        spec.gas = GasState::try_new_mixture(
            &[GasSpec::water_vapor_nist()],
            &composition,
            450.0,
            100000.0,
            Some(0),
        )
        .unwrap();
        assert_eq!(spec.gas.water_mole_fraction, 1.0);
        let steam = realize_plate_cavity(&spec, &input, 8000, 1.0).unwrap();
        assert_eq!(steam.len(), input.len());
        assert!(steam.iter().all(|sample| sample.is_finite()));
        assert!(steam.iter().any(|sample| *sample != 0.0));
        assert_ne!(air, steam);
        // ISO asked of steam refuses: an air admission does not transfer.
        let (_, admission) =
            fs_material::iso9613::TerrestrialAir::moist_air(293.15, 101_325.0, 0.0)
                .expect("air admission");
        assert!(matches!(
            realize_plate_cavity_in_air(&spec, &admission, &input, 8000, 1.0),
            Err(CavityPhsError::Plate(
                AcousticRealizeError::NonTerrestrialAir(_)
            ))
        ));
    }
}
