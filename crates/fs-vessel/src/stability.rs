//! Stages 1–2: the vessel profile and the Orr–Sommerfeld stability
//! objective. The vessel of revolution is a Chebyshev radius profile
//! r(z) on [0, 1] with a scalar lip-channel width; pouring drives a
//! quasi-steady thin film along the lip whose LOCAL Reynolds number
//! varies along the pour path with the film thickness. Each station
//! feeds a plane-Poiseuille Orr–Sommerfeld eigenproblem (fs-cheb —
//! the machinery validated against Orszag's Re_c = 5772.22), and the
//! objective is the MIN-MAX certified modal growth over stations and
//! modes: negative everywhere = every mode decays = the differentiable
//! laminarity proxy of Appendix C. The level-set lip (topology change)
//! is the recorded successor; the scalar channel is the smoke knob.

use fs_cheb::Cheb1;
use fs_cheb::orr_sommerfeld::growth_rates;

/// A vessel of revolution: radius profile + lip channel.
pub struct VesselProfile {
    /// Radius r(z), z ∈ [0, 1] (base to lip).
    pub radius: Cheb1,
    /// Lip channel width factor (thins the film as it narrows).
    pub lip_width: f64,
}

impl VesselProfile {
    /// A smooth default carafe: wide base, narrowing neck, flared lip.
    ///
    /// # Panics
    /// Never at fixture scale (`Cheb1::build` resolves the smooth
    /// profile well below its degree cap).
    #[must_use]
    pub fn carafe(lip_width: f64) -> VesselProfile {
        let radius = Cheb1::build(
            &|z: f64| {
                // det::sin: platform trig is a build-mode hazard (xo2k).
                0.05f64.mul_add(
                    fs_math::det::sin(2.5 * std::f64::consts::PI * z),
                    1.0 - 0.45 * z * z,
                )
            },
            0.0,
            1.0,
            64,
        );
        VesselProfile { radius, lip_width }
    }

    /// Film Reynolds numbers along the pour path (quasi-steady smoke
    /// proxy, DOCUMENTED as such): two competing branches make the
    /// design problem real — the THICKNESS branch (wide lip → thick
    /// film → high Re, thinned by viscosity) and the VELOCITY branch
    /// (narrow lip → the fixed pour rate accelerates through the
    /// channel → high Re, viscosity-blind). U-shaped in the lip width,
    /// so the growth objective has an interior optimum that MOVES with
    /// the fluid — the robustification stage has something to do.
    #[must_use]
    pub fn film_reynolds(&self, pour_rate: f64, viscosity: f64, stations: usize) -> Vec<f64> {
        let dr = self.radius.differentiate();
        (0..stations)
            .map(|k| {
                #[allow(clippy::cast_precision_loss)]
                let z = 0.7 + 0.3 * (k as f64 + 0.5) / stations as f64; // the lip run
                let slope = dr.eval(z).abs().max(0.05);
                let film = self.lip_width / slope.sqrt();
                // 1000× scale puts the design range astride the
                // Orr–Sommerfeld transition (Re ≈ 3200–6000 nominal).
                // Measured rejections: a clamped single branch pinned
                // every design at Re = 50 (constant objective); a
                // monotone unclamped branch pinned both optimizers at
                // the grid floor (no robustification story).
                (1000.0 * pour_rate * (film / viscosity + 2.5 / film)).max(50.0)
            })
            .collect()
    }
}

/// The min-max modal growth objective: the LARGEST real growth rate
/// over the pour-path stations and the first `modes` Orr–Sommerfeld
/// modes. Negative = certified-decaying everywhere (laminar proxy).
///
/// # Panics
/// If the eigensolver fails to converge (typed failure surfaced at
/// fixture scale).
#[must_use]
pub fn growth_objective(
    profile: &VesselProfile,
    pour_rate: f64,
    viscosity: f64,
    stations: usize,
    modes: usize,
) -> f64 {
    let alpha = 1.020_56; // the classic most-dangerous wavenumber
    let n = 32; // collocation size: Re_c reproduced to 4 digits at 48
    let mut worst = f64::NEG_INFINITY;
    for re in profile.film_reynolds(pour_rate, viscosity, stations) {
        let sig = growth_rates(re, alpha, n, modes).expect("OS eigensolve converges");
        for s in sig {
            worst = worst.max(s.re);
        }
    }
    worst
}
