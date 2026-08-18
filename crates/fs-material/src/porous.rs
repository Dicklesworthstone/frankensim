//! Johnson–Champoux–Allard(–Lafarge) porous-absorber facility (music
//! bead `frankensim-music-v8-root-3ez8g.5.4`).
//!
//! The LININGS claim — a material that is mostly air needs its own
//! acoustics — as an equivalent-fluid model: five JCA transport
//! parameters (open porosity, static airflow resistivity, tortuosity,
//! viscous and thermal characteristic lengths), plus the optional
//! Lafarge static thermal permeability `k0'` (JCAL). The card yields
//! the effective dynamic density and bulk modulus, characteristic
//! impedance and wavenumber, the SURFACE IMPEDANCE of a
//! finite-thickness layer on a rigid backing, and the
//! normal-incidence absorption coefficient.
//!
//! D23 consumption seam, decided here: the model lives in
//! fs-material next to `gas`/`visco` (a lining is a material law);
//! fs-duct's locally-reacting wall path and fs-phs's wall-pin family
//! consume `surface_impedance_rigid` as a frequency-dependent wall
//! admittance through thin adapters (displacement note in
//! CONTRACT.md — no new crate, no parallel implementation).
//!
//! Cards (licensing-first, no figure digitization anywhere):
//! - `melamine_uf_nguyen2024`: MEASURED characterization of the UF
//!   (Gray) melamine foam from Nguyen et al., Acta Acustica 8:54
//!   (2024), doi:10.1051/aacus/2024046, CC-BY-4.0 (HAL hal-04750927):
//!   Table 1 characterization row (k0 = 31.47e-10 m² so
//!   σ = η/k0, k0' = 75e-10 m², Λ = 141 µm, Λ' = 298 µm,
//!   α∞ = 1.06) with open porosity φ = 0.992 (pressure/mass method,
//!   Section 4).
//! - `basotect_foam02` / `pinta_fleece_foam02`: the inverse-identified
//!   JCAL rows from the FOAM 02 dataset (Zenodo 18242697, CC-BY-4.0;
//!   Acta Acustica 9(50) 2025, doi:10.1051/aacus/2025033) for the
//!   thickest/first-diameter configurations, whose fits are physical
//!   (thin-sample rows in that dataset pin parameters at fit bounds
//!   and are NOT card material — disclosed). Pinta Plano Polar is a
//!   polyester-fleece FELT-CLASS nonwoven, not wool piano felt; the
//!   wool-felt characterization is a recorded upgrade vein.
//!
//! Oracle: the committed FOAM 02 measured mean spectra
//! (`data/vv-corpus/acoustic/foam02-*-alpha.tsv`, machine-readable
//! CC-BY data, mean of 36 measurements each) — this implementation
//! must reproduce the measured absorption within an authored band
//! that covers the dataset's own identification residual (measured
//! max |Δα| 0.037 for Basotect, 0.023 for Pinta).
//!
//! NOT the hammer (compression is the hammer stack's claim), no
//! elastic-frame Biot waves (equivalent fluid only; the frame-wave
//! escalation is a named successor), no room acoustics.

use fs_math::c64::C64;
use fs_math::det;

use crate::MaterialError;

/// Ambient gas seen by the equivalent fluid (SI).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PorousAmbient {
    /// Density ρ0 [kg/m³].
    pub rho0: f64,
    /// Sound speed c0 [m/s].
    pub c0: f64,
    /// Dynamic viscosity η [Pa·s].
    pub eta: f64,
    /// Ratio of specific heats γ.
    pub gamma: f64,
    /// Ambient pressure P0 [Pa].
    pub p0: f64,
    /// Prandtl number.
    pub prandtl: f64,
}

impl PorousAmbient {
    /// Air at ~20 °C, matching `fs_couple::vibroacoustic`'s medium
    /// and the `GasState` primitive it is asserted against.
    #[must_use]
    pub const fn air() -> PorousAmbient {
        PorousAmbient {
            rho0: 1.204,
            c0: 343.0,
            eta: 1.81e-5,
            gamma: 1.4,
            p0: 101_325.0,
            prandtl: 0.71,
        }
    }
}

/// A JCA(L) porous-absorber card: typed parameters, named provenance,
/// validity band (matdb-shaped; migration must stay possible).
#[derive(Debug, Clone, PartialEq)]
pub struct PorousCard {
    /// Card name (registry spelling).
    pub name: &'static str,
    /// Open porosity φ ∈ (0, 1].
    pub porosity: f64,
    /// Static airflow resistivity σ [N·s/m⁴].
    pub flow_resistivity: f64,
    /// Tortuosity α∞ ≥ 1.
    pub tortuosity: f64,
    /// Viscous characteristic length Λ [m].
    pub viscous_length_m: f64,
    /// Thermal characteristic length Λ' [m].
    pub thermal_length_m: f64,
    /// Lafarge static thermal permeability k0' [m²]; `None` degrades
    /// to the Champoux–Allard thermal function.
    pub thermal_permeability_m2: Option<f64>,
    /// Where the numbers come from, verbatim-quotable.
    pub provenance: &'static str,
    /// Frequency band the card's source supports [Hz].
    pub valid_hz: (f64, f64),
}

impl PorousCard {
    /// Validate a card's parameters by name.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] with the offending field.
    pub fn validate(&self) -> Result<(), MaterialError> {
        let checks = [
            (
                self.porosity > 0.0 && self.porosity <= 1.0,
                "open porosity must lie in (0, 1]",
            ),
            (
                self.flow_resistivity.is_finite() && self.flow_resistivity > 0.0,
                "flow resistivity must be positive",
            ),
            (
                self.tortuosity.is_finite() && self.tortuosity >= 1.0,
                "tortuosity must be at least 1",
            ),
            (
                self.viscous_length_m.is_finite() && self.viscous_length_m > 0.0,
                "viscous characteristic length must be positive",
            ),
            (
                self.thermal_length_m.is_finite() && self.thermal_length_m > 0.0,
                "thermal characteristic length must be positive",
            ),
            (
                self.valid_hz.0 > 0.0 && self.valid_hz.1 > self.valid_hz.0,
                "validity band must be ordered and positive",
            ),
        ];
        for (ok, what) in checks {
            if !ok {
                return Err(MaterialError::Parameters {
                    what: what.to_string(),
                });
            }
        }
        if let Some(k0p) = self.thermal_permeability_m2 {
            if !(k0p.is_finite() && k0p > 0.0) {
                return Err(MaterialError::Parameters {
                    what: "thermal permeability must be positive when given".to_string(),
                });
            }
        }
        Ok(())
    }

    /// MEASURED melamine characterization: Nguyen et al. 2024 (Acta
    /// Acustica 8:54, CC-BY-4.0), UF (Gray) foam, Table 1
    /// characterization row + Section-4 porosity. σ = η/k0 with the
    /// card ambient's viscosity and the measured permeability
    /// k0 = 31.47e-10 m².
    #[must_use]
    pub fn melamine_uf_nguyen2024() -> PorousCard {
        PorousCard {
            name: "melamine-uf-nguyen2024",
            porosity: 0.992,
            flow_resistivity: 1.81e-5 / 31.47e-10,
            tortuosity: 1.06,
            viscous_length_m: 141.0e-6,
            thermal_length_m: 298.0e-6,
            thermal_permeability_m2: Some(75.0e-10),
            provenance: "Nguyen et al., Acta Acustica 8:54 (2024), doi:10.1051/aacus/2024046, \
                         CC-BY-4.0; Table 1 UF (Gray) characterization: k0 31.47±3.69e-10 m², \
                         k0' 75±12e-10 m², Λ 141±14 µm, Λ' 298±104 µm, α∞ 1.06±0.03; \
                         φ 0.992±0.005 (pressure/mass method)",
            valid_hz: (100.0, 4000.0),
        }
    }

    /// FOAM 02 Basotect (melamine) inverse-identified JCAL row,
    /// thickness 60 mm / diameter 98 mm — the physical fit (thin
    /// samples in the dataset pin parameters at bounds; refused).
    #[must_use]
    pub fn basotect_foam02() -> PorousCard {
        PorousCard {
            name: "basotect-foam02-60mm",
            porosity: 0.9987,
            flow_resistivity: 6332.0,
            tortuosity: 1.1188,
            viscous_length_m: 9.6526e-5,
            thermal_length_m: 1.5890e-3,
            thermal_permeability_m2: Some(5.7417e-9),
            provenance: "FOAM 02 dataset (Zenodo 18242697, CC-BY-4.0; Acta Acustica 9(50) 2025, \
                         doi:10.1051/aacus/2025033), JCAL_params_Basotect.csv row \
                         thickness=60,diameter=98",
            valid_hz: (150.0, 1600.0),
        }
    }

    /// FOAM 02 Pinta Plano Polar (polyester fleece, FELT-CLASS
    /// nonwoven — not wool piano felt; disclosed) inverse-identified
    /// JCAL row, thickness 50 mm / diameter 98 mm.
    #[must_use]
    pub fn pinta_fleece_foam02() -> PorousCard {
        PorousCard {
            name: "pinta-fleece-foam02-50mm",
            porosity: 0.9829,
            flow_resistivity: 11_003.0,
            tortuosity: 1.2897,
            viscous_length_m: 3.4931e-4,
            thermal_length_m: 1.2793e-4,
            thermal_permeability_m2: Some(1.0e-8),
            provenance: "FOAM 02 dataset (Zenodo 18242697, CC-BY-4.0; Acta Acustica 9(50) 2025, \
                         doi:10.1051/aacus/2025033), JCAL_params_Pinta.csv row \
                         thickness=50,diameter=98",
            valid_hz: (150.0, 1600.0),
        }
    }

    /// Effective dynamic density ρ_eff(ω) (Johnson et al.).
    ///
    /// # Errors
    /// Card validation; non-positive frequency.
    pub fn effective_density(
        &self,
        ambient: PorousAmbient,
        omega: f64,
    ) -> Result<C64, MaterialError> {
        self.validate()?;
        if !(omega.is_finite() && omega > 0.0) {
            return Err(MaterialError::State {
                what: "angular frequency must be positive".to_string(),
            });
        }
        let (phi, sigma, a_inf) = (self.porosity, self.flow_resistivity, self.tortuosity);
        let lam = self.viscous_length_m;
        let j = C64::new(0.0, 1.0);
        let g = (C64::new(1.0, 0.0)
            + j.scale(4.0 * a_inf * a_inf * ambient.eta * ambient.rho0 * omega)
                / C64::new(sigma * sigma * lam * lam * phi * phi, 0.0))
        .sqrt();
        let visc = C64::new(sigma * phi, 0.0) / (j.scale(omega * ambient.rho0 * a_inf)) * g;
        Ok((C64::new(1.0, 0.0) + visc).scale(a_inf * ambient.rho0 / phi))
    }

    /// Effective dynamic bulk modulus K_eff(ω): Lafarge when `k0'` is
    /// given, Champoux–Allard otherwise.
    ///
    /// # Errors
    /// Card validation; non-positive frequency.
    pub fn effective_bulk(&self, ambient: PorousAmbient, omega: f64) -> Result<C64, MaterialError> {
        self.validate()?;
        if !(omega.is_finite() && omega > 0.0) {
            return Err(MaterialError::State {
                what: "angular frequency must be positive".to_string(),
            });
        }
        let phi = self.porosity;
        let lam_t = self.thermal_length_m;
        let j = C64::new(0.0, 1.0);
        // Champoux–Allard's thermal function is the Lafarge form with
        // k0' = φ Λ'² / 8 — the degeneracy is explicit, not a fork.
        let k0p = self
            .thermal_permeability_m2
            .unwrap_or(phi * lam_t * lam_t / 8.0);
        let gp = (C64::new(1.0, 0.0)
            + j.scale(4.0 * k0p * k0p * ambient.prandtl * ambient.rho0 * omega)
                / C64::new(ambient.eta * lam_t * lam_t * phi * phi, 0.0))
        .sqrt();
        let denom = C64::new(1.0, 0.0)
            + C64::new(ambient.eta * phi, 0.0)
                / (j.scale(omega * ambient.rho0 * ambient.prandtl * k0p))
                * gp;
        let inner = C64::new(ambient.gamma, 0.0) - C64::new(ambient.gamma - 1.0, 0.0) / denom;
        Ok(C64::new(ambient.gamma * ambient.p0 / phi, 0.0) / inner)
    }

    /// Characteristic impedance Zc and complex wavenumber k.
    ///
    /// # Errors
    /// As [`PorousCard::effective_density`].
    pub fn characteristic(
        &self,
        ambient: PorousAmbient,
        omega: f64,
    ) -> Result<(C64, C64), MaterialError> {
        let rho = self.effective_density(ambient, omega)?;
        let k_mod = self.effective_bulk(ambient, omega)?;
        Ok(((rho * k_mod).sqrt(), (rho / k_mod).sqrt().scale(omega)))
    }

    /// Surface impedance of a thickness-`d` layer on a RIGID backing:
    /// `Zs = Zc · coth(j k d)` (e^{+jωt} convention).
    ///
    /// # Errors
    /// As [`PorousCard::effective_density`]; non-positive thickness.
    pub fn surface_impedance_rigid(
        &self,
        ambient: PorousAmbient,
        omega: f64,
        thickness_m: f64,
    ) -> Result<C64, MaterialError> {
        if !(thickness_m.is_finite() && thickness_m > 0.0) {
            return Err(MaterialError::Parameters {
                what: "layer thickness must be positive".to_string(),
            });
        }
        let (zc, k) = self.characteristic(ambient, omega)?;
        let gamma_d = C64::new(0.0, 1.0) * k.scale(thickness_m);
        Ok(zc * coth_c(gamma_d))
    }

    /// Normal-incidence absorption coefficient of the rigid-backed
    /// layer. Refuses frequencies outside the card's validity band.
    ///
    /// # Errors
    /// [`MaterialError::State`] outside `valid_hz`; parameter errors.
    pub fn absorption_normal(
        &self,
        ambient: PorousAmbient,
        f_hz: f64,
        thickness_m: f64,
    ) -> Result<f64, MaterialError> {
        if !(f_hz >= self.valid_hz.0 && f_hz <= self.valid_hz.1) {
            return Err(MaterialError::State {
                what: format!(
                    "frequency {f_hz} Hz outside the card's validity band \
                     [{}, {}] Hz",
                    self.valid_hz.0, self.valid_hz.1
                ),
            });
        }
        let omega = core::f64::consts::TAU * f_hz;
        let zs = self.surface_impedance_rigid(ambient, omega, thickness_m)?;
        let z0 = C64::new(ambient.rho0 * ambient.c0, 0.0);
        let r = (zs - z0) / (zs + z0);
        Ok(1.0 - r.norm_sq())
    }
}

/// Complex hyperbolic cotangent via `tanh(a+jb) =
/// (tanh a + j tan b) / (1 + j tanh a tan b)` (det-routed).
fn coth_c(z: C64) -> C64 {
    let (ta, tb) = (det::tanh(z.re), det::tan(z.im));
    let tanh = C64::new(ta, tb) / C64::new(1.0, ta * tb);
    tanh.recip()
}

#[cfg(test)]
mod porous_tests {
    use super::*;

    const AIR: PorousAmbient = PorousAmbient::air();

    fn read_alpha_tsv(name: &str) -> Vec<(f64, f64, f64)> {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(std::path::Path::parent)
            .expect("workspace root")
            .to_path_buf();
        let text = std::fs::read_to_string(root.join("data/vv-corpus/acoustic").join(name))
            .expect("committed FOAM 02 mean spectrum");
        text.lines()
            .filter(|l| !l.starts_with('#') && !l.starts_with("freq_hz"))
            .map(|l| {
                let mut it = l.split('\t');
                let f: f64 = it.next().expect("f").parse().expect("f");
                let m: f64 = it.next().expect("m").parse().expect("m");
                let s: f64 = it.next().expect("s").parse().expect("s");
                (f, m, s)
            })
            .collect()
    }

    fn spectrum_gate(card: &PorousCard, thickness_m: f64, tsv: &str, case: &str) {
        let rows = read_alpha_tsv(tsv);
        assert!(rows.len() > 700, "spectrum rows: {}", rows.len());
        let mut worst = 0.0f64;
        let mut sum = 0.0f64;
        for (f, m, _) in &rows {
            let a = card.absorption_normal(AIR, *f, thickness_m).expect("alpha");
            let d = (a - m).abs();
            worst = worst.max(d);
            sum += d;
        }
        let mean = sum / rows.len() as f64;
        // Authored from measurement (max |Δα| 0.037 / 0.023; mean
        // 0.015 / 0.008): the model must reproduce the MEASURED data
        // within the dataset's own identification-residual class.
        assert!(worst < 0.06, "worst |model - measured| {worst:.4}");
        assert!(mean < 0.03, "mean |model - measured| {mean:.4}");
        println!(
            "{{\"suite\":\"fs-material\",\"case\":\"{case}\",\"rows\":{},\
             \"worst_abs_dev\":{worst:.4},\"mean_abs_dev\":{mean:.4},\"verdict\":\"pass\"}}",
            rows.len()
        );
    }

    #[test]
    fn pa_001_basotect_reproduces_the_measured_spectrum() {
        spectrum_gate(
            &PorousCard::basotect_foam02(),
            0.060,
            "foam02-basotect-60mm-alpha.tsv",
            "pa-001-basotect",
        );
    }

    #[test]
    fn pa_002_pinta_reproduces_the_measured_spectrum() {
        spectrum_gate(
            &PorousCard::pinta_fleece_foam02(),
            0.050,
            "foam02-pinta-50mm-alpha.tsv",
            "pa-002-pinta",
        );
    }

    #[test]
    fn pa_003_analytic_degeneracies() {
        let card = PorousCard::melamine_uf_nguyen2024();
        // Thin-layer limit: absorption vanishes.
        let a_thin = card.absorption_normal(AIR, 500.0, 1.0e-6).expect("thin");
        assert!(a_thin < 1.0e-4, "thin-layer alpha {a_thin:.2e}");
        // Low-frequency stiffness-controlled limit of the rigid-backed
        // layer: Zs -> K_eff/(j omega d) — the layer is a pure
        // compliance. Compare against the closed form at kd << 1.
        let omega = core::f64::consts::TAU * 150.0;
        let d = 0.002;
        let zs = card.surface_impedance_rigid(AIR, omega, d).expect("zs");
        let k_eff = card.effective_bulk(AIR, omega).expect("k");
        let compliance = k_eff / (C64::new(0.0, 1.0).scale(omega * d));
        let rel = ((zs - compliance).abs()) / compliance.abs();
        assert!(rel < 0.01, "compliance degeneracy off by {rel:.4}");
        // Air-like card degenerates to the free medium: Zc -> rho0 c0.
        let air_like = PorousCard {
            name: "air-like",
            porosity: 1.0,
            flow_resistivity: 1.0e-3,
            tortuosity: 1.0,
            viscous_length_m: 1.0,
            thermal_length_m: 1.0,
            thermal_permeability_m2: None,
            provenance: "synthetic degeneracy",
            valid_hz: (1.0, 1.0e5),
        };
        let (zc, k) = air_like.characteristic(AIR, omega).expect("air");
        let z0 = AIR.rho0 * AIR.c0;
        // 2e-3 band: the ambient c0 = 343 is not exactly
        // sqrt(gamma P0 / rho0) = 343.25 (the constants are pinned to
        // the vibroacoustic medium, not derived) — the degeneracy
        // check absorbs that disclosed 0.07% mismatch.
        assert!((zc.re / z0 - 1.0).abs() < 2.0e-3 && zc.im.abs() / z0 < 2.0e-3);
        assert!((k.re / (omega / AIR.c0) - 1.0).abs() < 2.0e-3);
        println!(
            "{{\"suite\":\"fs-material\",\"case\":\"pa-003-degeneracies\",\
             \"thin_alpha\":{a_thin:.2e},\"compliance_rel\":{rel:.2e},\
             \"air_zc_re\":{:.2},\"verdict\":\"pass\"}}",
            zc.re
        );
    }

    #[test]
    fn pa_004_refusals_fire_by_name() {
        let mut card = PorousCard::melamine_uf_nguyen2024();
        card.porosity = 1.5;
        assert!(matches!(
            card.effective_density(AIR, 1000.0),
            Err(MaterialError::Parameters { .. })
        ));
        let card = PorousCard::basotect_foam02();
        // Outside the FOAM 02 band: refuse, never extrapolate quietly.
        assert!(matches!(
            card.absorption_normal(AIR, 20.0, 0.06),
            Err(MaterialError::State { .. })
        ));
        assert!(matches!(
            card.surface_impedance_rigid(AIR, 1000.0, -0.01),
            Err(MaterialError::Parameters { .. })
        ));
        println!("{{\"suite\":\"fs-material\",\"case\":\"pa-004-refusals\",\"verdict\":\"pass\"}}");
    }

    #[test]
    fn pa_005_lined_tube_composition() {
        // The composed demonstration: a plane-wave tube with the
        // melamine lining at one end. (a) A simulated TWO-MICROPHONE
        // measurement rig (the FOAM 02 method: two pressure samples in
        // the standing field, transfer-function decomposition) must
        // recover the direct-formula absorption — the composition
        // catches sign/convention drift the direct formula cannot.
        // (b) The lined tube's modal decay, predicted from |R|, is
        // FINITE and orders below the rigid tube's — the
        // decay-changes-as-published demonstration, logged.
        let card = PorousCard::melamine_uf_nguyen2024();
        let d = 0.030;
        let (length, mic1, mic2) = (0.5f64, 0.4f64, 0.35f64);
        let mut worst = 0.0f64;
        for f_hz in [300.0f64, 700.0, 1200.0] {
            let omega = core::f64::consts::TAU * f_hz;
            let k0 = omega / AIR.c0;
            let zs = card.surface_impedance_rigid(AIR, omega, d).expect("zs");
            let z0 = C64::new(AIR.rho0 * AIR.c0, 0.0);
            let r = (zs - z0) / (zs + z0);
            // Standing field p(x) = e^{-jkx} + R e^{+jkx} with x
            // measured from the sample face toward the source.
            let p_at = |x: f64| -> C64 {
                let (c, s) = (det::cos(k0 * x), det::sin(k0 * x));
                C64::new(c, -s) + r * C64::new(c, s)
            };
            let h12 = p_at(mic1) / p_at(mic2);
            // Two-microphone decomposition (ISO 10534-2 algebra).
            let s_ = mic1 - mic2;
            let (cs, ss) = (det::cos(k0 * s_), det::sin(k0 * s_));
            let e_minus = C64::new(cs, -ss);
            let e_plus = C64::new(cs, ss);
            let r_rig = (h12 - e_minus) / (e_plus - h12)
                * C64::new(det::cos(2.0 * k0 * mic2), -det::sin(2.0 * k0 * mic2));
            let alpha_rig = 1.0 - r_rig.norm_sq();
            let alpha_direct = 1.0 - r.norm_sq();
            let dev = (alpha_rig - alpha_direct).abs();
            worst = worst.max(dev);
            assert!(
                dev < 1.0e-10,
                "two-microphone rig vs direct at {f_hz} Hz: {dev:.3e}"
            );
        }
        // Modal decay of the fundamental with one lined end: energy
        // decays by |R|^2 per round trip, T = 2L/c.
        let f_mode = AIR.c0 / (2.0 * length);
        let omega = core::f64::consts::TAU * f_mode;
        let zs = card.surface_impedance_rigid(AIR, omega, d).expect("zs");
        let z0 = C64::new(AIR.rho0 * AIR.c0, 0.0);
        let r_mag = ((zs - z0) / (zs + z0)).abs();
        let decay_per_s = -AIR.c0 / (2.0 * length) * det::ln(r_mag * r_mag);
        let t60 = 6.9077 / decay_per_s;
        assert!(
            t60.is_finite() && t60 > 0.0 && t60 < 5.0,
            "lined-tube T60 {t60:.3} s"
        );
        println!(
            "{{\"suite\":\"fs-material\",\"case\":\"pa-005-lined-tube\",\
             \"rig_vs_direct_worst\":{worst:.3e},\"fundamental_hz\":{f_mode:.0},\
             \"lined_r_mag\":{r_mag:.4},\"lined_t60_s\":{t60:.3},\
             \"rigid_t60\":\"infinite (|R| = 1)\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn pa_006_cards_validate_and_are_deterministic() {
        for card in [
            PorousCard::melamine_uf_nguyen2024(),
            PorousCard::basotect_foam02(),
            PorousCard::pinta_fleece_foam02(),
        ] {
            card.validate().expect("card validates");
            let a1 = card.absorption_normal(AIR, 500.0, 0.03).expect("a");
            let a2 = card.absorption_normal(AIR, 500.0, 0.03).expect("a");
            assert!(a1.to_bits() == a2.to_bits(), "bitwise determinism");
            assert!(a1 > 0.0 && a1 < 1.0);
        }
        println!("{{\"suite\":\"fs-material\",\"case\":\"pa-006-cards\",\"verdict\":\"pass\"}}");
    }
}
