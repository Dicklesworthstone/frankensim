//! Gas-related losses: cylinder motion and free-field observer absorption.
//!
//! A generic gas uses its own Stokes-Kirchhoff viscothermal loss. The
//! air-specific ISO 9613 law requires an explicit choice; temperature and
//! pressure alone do not establish that an arbitrary mixture is air.

/// Refinable causal gas inertia and viscous memory for circular transverse motion.
pub mod cylinder;

use fs_fft::{C64 as FftC64, Fft};
use fs_material::gas::GasState;
use fs_material::iso9613::iso9613_absorption;
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
    /// Explicit caller assertion of terrestrial air and ISO applicability.
    /// The ISO term already includes classical loss; it must not be added to
    /// Stokes-Kirchhoff. Out-of-window requests refuse rather than switch laws.
    Iso9613Air {
        /// Relative humidity as a fraction in [0, 1], not water mole fraction.
        relative_humidity: f64,
    },
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
/// Invalid samples, time/range, incoherent gas, unavailable ISO conditions, or
/// unrepresentable absorption/FFT output refuse WITHOUT changing the history.
pub fn absorb_pressure_history_with_model(
    pressure_pa: &mut [f64],
    dt: f64,
    range_m: f64,
    gas: &GasState,
    model: GasAbsorptionModel,
) -> Result<(), AcousticRealizeError> {
    gas.validate_coherence().map_err(|_| invalid("observer absorption needs a coherent gas state"))?;
    if let GasAbsorptionModel::Iso9613Air { relative_humidity } = model {
        // Admit the declared air conditions even for an empty or zero-range
        // history. A failed air law is not permission to use a different law.
        iso9613_absorption(gas, relative_humidity, 1.0)
            .map_err(|_| invalid("explicit ISO 9613 air conditions are outside the model domain"))?;
    }
    transform_history(pressure_pa, dt, range_m, |omega| match model {
        GasAbsorptionModel::StokesKirchhoff => Ok(gas.stokes_kirchhoff_absorption(omega)),
        GasAbsorptionModel::Iso9613Air { relative_humidity } => {
            iso9613_absorption(gas, relative_humidity, omega)
                .map_err(|_| invalid("explicit ISO 9613 absorption could not be evaluated"))
        }
    })
}

/// Compatibility path for callers that explicitly supply terrestrial-air
/// humidity. Retains the legacy ISO-to-classical fallback, including NaN
/// humidity used by older generic-gas callers. New generic consumers should
/// use [`absorb_pressure_history_with_model`] and [`GasAbsorptionModel`] so
/// physical model selection is not encoded in an invalid numerical value.
///
/// `range_m <= 0` and histories shorter than four samples remain no-ops.
/// Invalid/unrepresentable transformations leave the history unchanged.
pub fn absorb_pressure_history(
    pressure_pa: &mut [f64],
    dt: f64,
    range_m: f64,
    gas: &GasState,
    relative_humidity: f64,
) {
    if pressure_pa.len() < 4 || !(dt > 0.0 && range_m > 0.0) {
        return;
    }
    let _ = transform_history(pressure_pa, dt, range_m, |omega| {
        Ok(iso9613_absorption(gas, relative_humidity, omega)
            .unwrap_or_else(|_| gas.stokes_kirchhoff_absorption(omega)))
    });
}

/// Molecular relaxation that an observer path did not apply, with an upper
/// bound on its effect at the listener. A render may omit relaxation only with
/// this receipt, never silently (orchestrator ruling of 2026-10-23 on pexfp.18).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RelaxationOmission {
    /// Listener distance the bound covers [m].
    pub range_m: f64,
    /// Highest frequency of the observed band [Hz].
    pub band_hz: f64,
    /// Upper bound on the omitted attenuation over the band [dB].
    pub bound_db: f64,
}

/// What an observer path applied for gas absorption.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ObserverAbsorptionReceipt {
    /// The law applied to the pressure history.
    pub law: GasAbsorptionModel,
    /// Molecular relaxation that the law leaves out, with its bound. `None`
    /// when the law includes it (ISO 9613-1 for terrestrial air).
    pub relaxation_omitted: Option<RelaxationOmission>,
}

/// Bound the attenuation that omitting molecular relaxation can leave out:
/// the worst ISO 9613-1 absorption of standard air at `band_hz`, times
/// `range_m`, in dB.
///
/// "Worst" is the maximum over ISO's tabulated domain (-20..+50 C, 10..100 %
/// RH, 1 atm), sampled on a 1 K by 1 % grid. ISO absorption rises
/// monotonically with frequency, so the band's top bounds the whole band.
/// ISO's figure includes classical loss, so the bound also covers the
/// viscothermal part.
///
/// This bound is the ruling's convention for mixtures whose relaxation rates
/// are unsourced. It is not a proof for such mixtures: a CO2-rich gas can
/// relax more strongly than air.
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
    let bound_db = worst * 1e-3 * range_m;
    if !bound_db.is_finite() {
        return Err(invalid("relaxation omission bound is unrepresentable"));
    }
    Ok(RelaxationOmission { range_m, band_hz, bound_db })
}

/// Propagate a pressure history over `range_m` through a gas whose molecular
/// relaxation has been evaluated, as a causal complex transfer.
///
/// Every positive-frequency DFT bin is multiplied by
/// `exp(-Im k r) exp(-i (Re k - w / c_inf) r)`, where `k` is
/// [`RelaxationState::causal_wavenumber`] and `c_inf` is the frozen sound
/// speed. This is fs-duct's `e^{-iwt}` wave `e^{i(kr - wt)}` written for the
/// DFT's `e^{+iwt}` synthesis. Unlike the zero-phase map above, absorption and
/// dispersion come from one wavenumber, counted once, and the frequency
/// response is that of a causal medium.
///
/// The bulk delay `r / c_inf` is not applied; it is returned so the caller can
/// place the arrival. The Nyquist bin of the even transform keeps only the real
/// part of its transfer, as a real signal requires. Sampling the response on
/// the DFT grid is exact for the bins; between them it inherits the usual
/// band-limiting of a response that is not itself band-limited. A zero range
/// and histories shorter than four samples are no-ops.
///
/// # Errors
/// Invalid samples, time/range, or an unrepresentable wavenumber or output
/// refuse WITHOUT changing the history.
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
    let c_inf = relaxation.frozen_sound_speed_m_s();
    let delay_s = range_m / c_inf;
    if !(delay_s.is_finite() && delay_s >= 0.0) {
        return Err(invalid("relaxation bulk delay is unrepresentable"));
    }
    if pressure_pa.len() < 4 || range_m == 0.0 { return Ok(delay_s); }
    let n = pressure_pa.len();
    let n_fft = n.checked_next_power_of_two()
        .ok_or_else(|| invalid("propagation FFT length is unrepresentable"))?.max(8);
    let fft = Fft::new(n_fft);
    let mut buf: Vec<FftC64> = (0..n_fft)
        .map(|i| FftC64::new(if i < n { pressure_pa[i] } else { 0.0 }, 0.0))
        .collect();
    let mut scratch = vec![FftC64::new(0.0, 0.0); n_fft];
    fft.forward(&mut buf, &mut scratch);
    for (k, bin) in buf.iter_mut().enumerate().take(n_fft / 2 + 1) {
        if k == 0 { continue; }
        let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
        let wave = relaxation.causal_wavenumber(omega)
            .map_err(|_| invalid("relaxation wavenumber could not be evaluated"))?;
        let magnitude = det::exp(-wave.im * range_m);
        let (sin, cos) = det::sin_cos(-(wave.re - omega / c_inf) * range_m);
        let (h_re, h_im) = if k == n_fft / 2 { (magnitude * cos, 0.0) } else { (magnitude * cos, magnitude * sin) };
        *bin = FftC64::new(bin.re * h_re - bin.im * h_im, bin.re * h_im + bin.im * h_re);
    }
    for k in 1..n_fft / 2 {
        let c = buf[k];
        buf[n_fft - k] = FftC64::new(c.re, -c.im);
    }
    fft.inverse(&mut buf, &mut scratch);
    if buf[..n].iter().any(|value| !value.re.is_finite()) {
        return Err(invalid("propagated pressure is unrepresentable"));
    }
    for (i, sample) in pressure_pa.iter_mut().enumerate() { *sample = buf[i].re; }
    Ok(delay_s)
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

    #[test]
    fn farther_path_kills_high_frequency_more() {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 288.15, 101_325.0).expect("air");
        let n = 256;
        let dt = 1.0 / 8_000.0;
        let mut near: Vec<f64> = (0..n)
            .map(|i| (core::f64::consts::TAU * 2_000.0 * f64::from(i) * dt).sin())
            .collect();
        let mut far = near.clone();
        absorb_pressure_history(&mut near, dt, 1.0, &gas, 0.50);
        absorb_pressure_history(&mut far, dt, 200.0, &gas, 0.50);
        let en: f64 = near.iter().map(|x| x * x).sum();
        let ef: f64 = far.iter().map(|x| x * x).sum();
        assert!(ef < en * 0.98, "200 m must absorb more 2 kHz energy than 1 m");
    }

    #[test]
    fn g0_history_beyond_legacy_fft_cap_is_transformed_in_full() {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 288.15, 101_325.0).expect("air");
        let mut pressure = vec![0.0; 8_193];
        pressure[8_192] = 1.0;
        absorb_pressure_history(&mut pressure, 1.0 / 48_000.0, 200.0, &gas, 0.50);
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

    #[test]
    fn explicit_air_refuses_instead_of_silently_using_another_law() {
        let gas = GasState::try_new(&GasSpec::water_vapor_nist(), 450.0, 100000.0).unwrap();
        let original = vec![1.0, -1.0, 0.5, -0.5, 1.0, 0.0, -1.0, 0.0];
        let mut signal = original.clone();
        assert!(absorb_pressure_history_with_model(&mut signal, 1.0 / 8000.0, 100.0,
            &gas, GasAbsorptionModel::Iso9613Air { relative_humidity: 0.5 }).is_err());
        assert_eq!(signal, original);
        let mut forged = gas;
        forged.density *= 2.0;
        assert!(absorb_pressure_history_with_model(&mut signal, 1.0 / 8000.0, 100.0,
            &forged, GasAbsorptionModel::StokesKirchhoff).is_err());
        assert_eq!(signal, original);
    }

    #[test]
    fn old_nan_selector_and_explicit_generic_share_the_same_operator() {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap();
        for (humidity, model) in [
            (f64::NAN, GasAbsorptionModel::StokesKirchhoff),
            (0.5, GasAbsorptionModel::Iso9613Air { relative_humidity: 0.5 }),
        ] {
            let mut legacy = vec![1.0, -1.0, 0.5, -0.5, 1.0, 0.0, -1.0, 0.0];
            let mut explicit = legacy.clone();
            absorb_pressure_history(&mut legacy, 1.0 / 8000.0, 100.0, &gas, humidity);
            absorb_pressure_history_with_model(&mut explicit, 1.0 / 8000.0, 100.0, &gas, model).unwrap();
            assert_eq!(legacy, explicit);
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

    /// The transfer is exactly the law's wavenumber in every bin, and the
    /// propagated impulse of a strongly dispersive reservoir arrives no earlier
    /// than the bulk delay: pre-arrival energy at most 1e-12 of the total.
    #[test]
    fn relaxation_transfer_is_the_law_per_bin_and_causal() {
        let state = relaxing_air(true);
        let (n, dt, range, n0) = (1usize << 16, 1.0 / 48_000.0, 10.0, 1usize << 14);
        let mut x = vec![0.0; n];
        x[n0] = 1.0;
        let delay = propagate_pressure_history_relaxing(&mut x, dt, range, &state).unwrap();
        assert_eq!(delay.to_bits(), (range / state.frozen_sound_speed_m_s()).to_bits());
        let total: f64 = x.iter().map(|v| v * v).sum();
        let early: f64 = x[..n0].iter().map(|v| v * v).sum();
        println!(
            "{{\"bead\":\"M0.1-W2\",\"test\":\"relaxation_transfer_causal\",\"range_m\":{range},\"pre_arrival_energy_ratio\":{:e},\"tolerance\":1e-12}}",
            early / total
        );
        assert!(early <= 1e-12 * total, "pre-arrival energy ratio {:e}", early / total);
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
            // Bins attenuated below 1e-6 hold only transform roundoff.
            if magnitude < 1e-6 {
                continue;
            }
            let phase = -(wave.re - omega / c_inf) * range - core::f64::consts::TAU * (k * n0) as f64 / n as f64;
            let (er, ei) = (magnitude * phase.cos(), magnitude * phase.sin());
            worst = worst.max(((y[k].re - er).hypot(y[k].im - ei)) / magnitude);
        }
        assert!(worst <= 1e-9, "per-bin transfer error {worst:e}");
    }

    /// The omission bound is the worst standard-air ISO absorption at the band
    /// top times the range: linear in range, rising with the band, never below
    /// ISO at any grid state it covers, and refusing invalid inputs.
    #[test]
    fn relaxation_omission_bounds_standard_air_absorption() {
        let one = relaxation_omission(1.0, 24_000.0).unwrap();
        let two = relaxation_omission(2.0, 24_000.0).unwrap();
        assert_eq!(two.bound_db.to_bits(), (2.0 * one.bound_db).to_bits());
        let lower = relaxation_omission(1.0, 10_000.0).unwrap();
        assert!(lower.bound_db < one.bound_db);
        let omega = core::f64::consts::TAU * 24_000.0;
        for (t, rh) in [(253.15, 0.1), (293.15, 0.5), (303.15, 0.2), (323.15, 1.0)] {
            let iso = fs_material::iso9613::iso9613_absorption_db_per_km(t, 101_325.0, rh, omega).unwrap() * 1e-3;
            assert!(one.bound_db >= iso, "bound {} below ISO {iso} at {t} K, RH {rh}", one.bound_db);
        }
        assert_eq!(relaxation_omission(0.0, 24_000.0).unwrap().bound_db, 0.0);
        for (r, f) in [(-1.0, 1_000.0), (f64::NAN, 1_000.0), (1.0, 0.0), (1.0, f64::INFINITY)] {
            assert!(relaxation_omission(r, f).is_err());
        }
        println!(
            "{{\"bead\":\"M0.1-W2\",\"test\":\"relaxation_omission\",\"band_hz\":24000,\"bound_db_1m\":{},\"bound_db_2m\":{},\"bound_db_1m_10kHz\":{}}}",
            one.bound_db, two.bound_db, lower.bound_db
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
}
