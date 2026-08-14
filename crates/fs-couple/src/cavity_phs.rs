//! Time-domain plate × lumped Helmholtz cavity.
//!
//! The plate is an `fs-phs` modal bank. The cavity is the flow-driven
//! Helmholtz pHS. They join through a `transformer` (area already in
//! the plate's flow port). A bottle and a vented enclosure are the
//! same objects. There is no guitar type.

use fs_material::visco::{ThermoelasticZener, loss_factor_to_zeta};
use fs_phs::{
    MouthFlange, PortHamiltonian, QuadraticStorage, compact_radiation_impedance,
    helmholtz_resonator_flow, step, transformer,
};

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
}

impl core::fmt::Display for CavityPhsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid { what } => write!(f, "FS-COUPLE-CAVITY: {what}"),
            Self::Phs(e) => write!(f, "FS-COUPLE-CAVITY-PHS: {e}"),
        }
    }
}

impl std::error::Error for CavityPhsError {}

/// A modal plate facing a lumped Helmholtz volume.
#[derive(Debug, Clone)]
pub struct PlateCavitySpec {
    /// Plate modal frequencies [rad/s].
    pub omegas: Vec<f64>,
    /// Authored viscous ratios (thermoelastic is added on top).
    pub zetas: Vec<f64>,
    /// Drive weights at the external force (mass-normalized).
    pub drive: Vec<f64>,
    /// Monopole areas [m²] (volume velocity = A · v).
    pub areas: Vec<f64>,
    /// Plate thickness [m]. Zero skips thermoelastic loss.
    pub thickness_m: f64,
    /// Plate density [kg/m³] (selects Al vs steel thermoelastic).
    pub plate_density_kg_m3: f64,
    /// Cavity volume [m³].
    pub volume_m3: f64,
    /// Neck radius [m].
    pub neck_radius_m: f64,
    /// Neck length [m] (end correction is applied inside the pHS).
    pub neck_length_m: f64,
    /// Gas density [kg/m³].
    pub density: f64,
    /// Sound speed [m/s].
    pub sound_speed: f64,
    /// Gas temperature [K] (thermoelastic T₀).
    pub temperature_k: f64,
    /// Relative humidity in `[0, 1]` for the observer path.
    pub relative_humidity: f64,
}

/// Realize observer pressure of a driven plate on a Helmholtz volume.
///
/// `force_n[i]` is the external force at sample `i`. The far-field
/// observer is the baffled compact sum `ρ A ÿ / (2π r)`.
///
/// # Errors
/// Empty/mismatched modes, non-physical cavity, or a pHS refusal.
pub fn realize_plate_cavity(
    spec: &PlateCavitySpec,
    force_n: &[f64],
    sample_rate_hz: u32,
    listener_m: f64,
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
    let mut zetas = spec.zetas.clone();
    if spec.thickness_m > 0.0 && spec.plate_density_kg_m3 > 0.0 {
        let te = if spec.plate_density_kg_m3 > 5_000.0 {
            ThermoelasticZener::structural_steel(spec.temperature_k)
        } else {
            ThermoelasticZener::aluminum(spec.temperature_k)
        }
        .map_err(|e| CavityPhsError::Phs(e.to_string()))?;
        for (z, &w) in zetas.iter_mut().zip(&spec.omegas) {
            *z += loss_factor_to_zeta(te.loss_factor(w, spec.thickness_m));
        }
    }
    let plate = plate_force_and_flow(&spec.omegas, &zetas, &spec.drive, &spec.areas)
        .map_err(|e| CavityPhsError::Phs(e.to_string()))?;
    let pi = core::f64::consts::PI;
    let neck_area = pi * spec.neck_radius_m * spec.neck_radius_m;
    let l_eff = spec.neck_length_m + 2.0 * (8.0 / (3.0 * pi)) * spec.neck_radius_m;
    let omega0 = spec.sound_speed * (neck_area / (spec.volume_m3 * l_eff)).sqrt();
    let r_rad = compact_radiation_impedance(
        spec.density,
        spec.sound_speed,
        spec.neck_radius_m,
        omega0,
        MouthFlange::Unflanged,
    )
    .map(|(r, _)| r)
    .unwrap_or(0.0);
    let cavity = helmholtz_resonator_flow(
        spec.volume_m3,
        spec.neck_radius_m,
        spec.neck_length_m,
        spec.density,
        spec.sound_speed,
        r_rad,
    )
    .map_err(|e| CavityPhsError::Phs(e.to_string()))?;
    // Area is already in the plate's second port. The transformer is
    // the power-conserving F = p A, U = A v join — not a staggered
    // clock.
    let sys =
        transformer(plate, cavity, 1, 0, 1.0).map_err(|e| CavityPhsError::Phs(e.to_string()))?;
    let dt = 1.0 / f64::from(sample_rate_hz);
    let mut x = vec![0.0; sys.state_dim()];
    let mut out = Vec::with_capacity(force_n.len());
    let two_pi = 2.0 * core::f64::consts::PI;
    let n_plate = 2 * n_m;
    for &f_ext in force_n {
        let rec = step(&sys, &x, &[f_ext], dt).map_err(|e| CavityPhsError::Phs(e.to_string()))?;
        let mut p_obs = 0.0;
        for k in 0..n_m {
            let acc = (rec.x[2 * k + 1] - x[2 * k + 1]) / dt;
            p_obs += spec.density * spec.areas[k] * acc / (two_pi * listener_m);
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
    if let Ok(gas) = fs_material::gas::GasState::try_new(
        &fs_material::gas::GasSpec::dry_air_ussa1976(),
        spec.temperature_k,
        101_325.0,
    ) {
        crate::air_path::absorb_pressure_history(
            &mut out,
            dt,
            listener_m,
            &gas,
            spec.relative_humidity,
        );
    }
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

    fn panel(volume: f64) -> PlateCavitySpec {
        PlateCavitySpec {
            omegas: vec![6.0e2],
            zetas: vec![0.006],
            drive: vec![1.0],
            areas: vec![0.08],
            thickness_m: 0.002,
            plate_density_kg_m3: 2_700.0,
            volume_m3: volume,
            neck_radius_m: 0.02,
            neck_length_m: 0.03,
            density: 1.2,
            sound_speed: 343.0,
            temperature_k: 293.0,
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
    fn thermoelastic_changes_the_waveform() {
        let f = pulse(800, 2.0);
        let mut with = panel(0.008);
        let mut bare = with.clone();
        bare.thickness_m = 0.0;
        with.thickness_m = 0.003;
        let a = realize_plate_cavity(&bare, &f, 8_000, 1.0).expect("bare");
        let b = realize_plate_cavity(&with, &f, 8_000, 1.0).expect("te");
        let err: f64 = a.iter().zip(&b).map(|(x, y)| (x - y).abs()).sum();
        assert!(err > 1.0e-10, "thermoelastic ζ must move the waveform");
    }
}
