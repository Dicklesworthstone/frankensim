//! Geometry-derived compact baffled-piston load on the existing wave-port owner.
//!
//! Rayleigh gives Z/Zc = (ka)^2/2 - i alpha ka + O((ka)^3), alpha=8/(3*pi),
//! under the workspace exp(-i omega t) convention. Matching these TWO leading
//! terms to the positive-real R*s/(s+p) gives p=2*alpha*c/a and R=2*alpha^2*Zc.
//! This is the existing fs-vfit relaxation impedance, not a new time integrator,
//! fitted reflection, constant damping or an outlet end-length adjustment.
//!
//! Both the continuous model and its actual bilinear response are checked against
//! fs-phs's full Rayleigh piston on a fixed frequency grid before publication.
//! The checks are finite-sample estimates, not a continuum error enclosure. The
//! compact model is intentionally unavailable outside ka<=0.5; no full-band,
//! unbaffled-pipe or measured-instrument accuracy is implied. Its stored near-field
//! energy is separate from the positive loss representing escaped acoustic energy.
//! It is not an M2.5-W3b realization and is excluded from W3b acceptance; that
//! is fs-vfit's passive multiport model, audited in `tests/closed_loop_radiation.rs`.
//!
//! An open end may instead be declared [`RadiationFlange::Unflanged`] (an explicit,
//! per-aperture opt-in; the default is the baffled piston above). The same two-term
//! match is then made to the planar-mode radiation of a thin-wall UNFLANGED circular
//! pipe: Z/Zc = (ka)^2/4 - i eta ka + O((ka)^3) with eta = 0.6133, i.e.
//! |R| = 1 - beta (ka)^2 with beta = 1/2 and end correction L/a = eta (Silva,
//! Guillemain, Kergomard, Mallaroni & Norris, J. Sound Vib. 322 (2009) 255-263,
//! Eqs. (10) and (11), whose static values are Levine & Schwinger's). Writing the
//! general two-term match Z/Zc = -i eta ka + (beta/2)(ka)^2 gives p = 2 eta c/(beta a)
//! and R = 2 eta^2 Zc/beta; the baffled piston is eta = alpha, beta = 1. For the
//! unflanged pipe p = 4 eta c/a and R = 4 eta^2 Zc. The same 32-point 5 % gates
//! apply, with fs-duct's exact Levine-Schwinger evaluation as the oracle. This is
//! the law the tone-hole end correction t_r of Garcia Mayen (2021, Eq. 14) assumes
//! ("a tube without flange"); no hole-in-a-cylinder or intermediate-flange formula
//! is implemented because none was available from an opened source.

/// A shared exterior gas for compact radiation reaction and microphone propagation.
pub mod environment;

use crate::acoustic_realize::AcousticRealizeError;
use super::network::NetworkNode;
use fs_exec::CancelGate;
use fs_vfit::impedance::SeriesImpedanceSpec;
use fs_vfit::relaxation::{RelaxationImpedanceSpec, RelaxationTerm};

/// Exterior law of one compact open end. `Baffled` is the original piston in an
/// infinite rigid baffle; `Unflanged` is a thin-wall unflanged circular pipe end
/// (see the module notes). Selecting `Unflanged` changes the REACTION only: the
/// exterior pressure observation of the same opening stays the baffled half-space
/// monopole, so reaction and observation then use different exterior models.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RadiationFlange {
    /// Piston in an infinite rigid baffle: alpha = 8/(3 pi), beta = 1.
    #[default]
    Baffled,
    /// Unflanged pipe, planar mode: eta = 0.6133, beta = 1/2.
    Unflanged,
}

/// Static end correction of an unflanged circular pipe, in radii (Silva et al. 2009,
/// Eq. 11; the same historical four-digit value fs-duct's `UnflangedOpen` termination
/// uses). fs-duct's exact Levine-Schwinger evaluation gives 0.6127; the 0.1 % difference
/// is far inside the 5 % admission gates below.
const UNFLANGED_END_CORRECTION: f64 = 0.6133;

/// Admitted physical source, numerical clock and checked response of the load.
/// Fields are read-only so geometry cannot drift from the generated coefficients.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaffledRadiationLoad {
    radius_m: f64,
    density_kg_m3: f64,
    sound_speed_m_s: f64,
    time_step_s: f64,
    maximum_frequency_hz: f64,
    max_complex_relative_error: f64,
    max_resistance_relative_error: f64,
    load: RelaxationImpedanceSpec,
    flange: RadiationFlange,
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

impl BaffledRadiationLoad {
    /// The original baffled-piston load: exactly [`Self::with_flange`] with
    /// [`RadiationFlange::Baffled`], which keeps every byte of its coefficients.
    ///
    /// # Errors
    /// As [`Self::with_flange`].
    pub fn new(
        radius_m: f64,
        density_kg_m3: f64,
        sound_speed_m_s: f64,
        time_step_s: f64,
        maximum_frequency_hz: f64,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        Self::with_flange(radius_m, density_kg_m3, sound_speed_m_s, time_step_s,
            maximum_frequency_hz, RadiationFlange::Baffled, gate)
    }

    /// Derive the sole relaxation branch from the ACTUAL mouth radius and fluid.
    /// All 32 positive frequencies spanning (0, maximum_frequency_hz] must meet
    /// fixed 5% relative complex AND radiation-resistance checks for the analogue
    /// and discrete responses. The analytical zero-frequency limits are exact.
    /// The `Baffled` oracle is the shared closed-form Bessel–Struve piston; it has
    /// no spatial quadrature error. The analogue-model test requires at most
    /// 2.2% complex and 4.2% resistance discrepancy on ka ∈ (0, 0.5]
    /// (a sampled estimate pending native qualification;
    /// `compact_load_two_term_exact_and_error_recorded` records the grid).
    /// The `Unflanged` oracle is fs-duct's exact Levine-Schwinger evaluation
    /// (256 Gauss-Legendre nodes); its recorded discrepancy is in
    /// `unflanged_load_two_term_exact_and_oracle_error_recorded`.
    ///
    /// `dt*fmax <= 0.1` limits temporal warping and `ka_max <= 0.5` declares the
    /// low-frequency model domain. They do not prove that an arbitrary transient
    /// has no higher-frequency content. No coefficients are retuned after checks.
    ///
    /// # Errors
    /// Invalid/unrepresentable geometry or medium, cancellation, an unresolved
    /// clock, oracle failure, or either checked discrepancy above 5%.
    #[allow(clippy::too_many_arguments)]
    pub fn with_flange(
        radius_m: f64,
        density_kg_m3: f64,
        sound_speed_m_s: f64,
        time_step_s: f64,
        maximum_frequency_hz: f64,
        flange: RadiationFlange,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        let unflanged = flange == RadiationFlange::Unflanged;
        if ![radius_m,density_kg_m3,sound_speed_m_s,time_step_s,maximum_frequency_hz]
            .iter().all(|v| v.is_finite() && *v > 0.0)
        {
            return Err(invalid(if unflanged {
                "unflanged radiation requires positive finite geometry, fluid, clock and band"
            } else {
                "baffled radiation requires positive finite geometry, fluid, clock and band"
            }));
        }
        let ka=core::f64::consts::TAU*maximum_frequency_hz/sound_speed_m_s*radius_m;
        if !ka.is_finite() || ka>0.5 || time_step_s*maximum_frequency_hz>0.1 {
            return Err(invalid(if unflanged {
                "compact unflanged radiation requires ka <= 0.5 and dt*fmax <= 0.1"
            } else {
                "compact baffled radiation requires ka <= 0.5 and dt*fmax <= 0.1"
            }));
        }
        let area=core::f64::consts::PI*radius_m*radius_m;
        let z=density_kg_m3*sound_speed_m_s/area;
        let term=match flange {
            RadiationFlange::Baffled => {
                let alpha=8.0/(3.0*core::f64::consts::PI);
                RelaxationTerm {resistance_pa_s_m3:2.0*alpha*alpha*z,
                    rate_per_s:2.0*alpha*sound_speed_m_s/radius_m}
            }
            RadiationFlange::Unflanged => {
                let eta=UNFLANGED_END_CORRECTION;
                RelaxationTerm {resistance_pa_s_m3:4.0*eta*eta*z,
                    rate_per_s:4.0*eta*sound_speed_m_s/radius_m}
            }
        };
        let load=RelaxationImpedanceSpec::new(SeriesImpedanceSpec {
            resistance_pa_s_m3:0.0,inertance_pa_s2_m3:0.0,compliance_m3_pa:None,
        }, &[term]).map_err(|e|AcousticRealizeError::Nonlinear(e.to_string()))?;
        // The existing owner also checks representability of its stored energy
        // and discrete coefficients before any physical system is constructed.
        fs_vfit::relaxation::RelaxationImpedance::new(load,z,time_step_s)
            .map_err(|e|AcousticRealizeError::Nonlinear(e.to_string()))?;
        let mut result=Self {radius_m,density_kg_m3,sound_speed_m_s,time_step_s,
            maximum_frequency_hz,max_complex_relative_error:0.0,max_resistance_relative_error:0.0,load,flange};
        for i in 1..=32 {
            if gate.is_requested() {return Err(invalid("baffled radiation construction cancelled"));}
            let frequency=maximum_frequency_hz*f64::from(i)/32.0;
            let omega=core::f64::consts::TAU*frequency;
            let (r,x)=match flange {
                RadiationFlange::Baffled => fs_phs::baffled_piston_impedance(
                    density_kg_m3,sound_speed_m_s,radius_m,omega)
                    .map_err(|e|AcousticRealizeError::Nonlinear(e.to_string()))?,
                RadiationFlange::Unflanged => {
                    let ka_i=omega/sound_speed_m_s*radius_m;
                    let normalized=fs_duct::levine_schwinger::levine_schwinger(ka_i,
                        fs_duct::levine_schwinger::LsQuadrature {nodes:256})
                        .map_err(|e|AcousticRealizeError::Nonlinear(e.to_string()))?
                        .normalized_impedance();
                    (normalized.re*z,normalized.im*z)
                }
            };
            if !r.is_finite() || r<=0.0 || !x.is_finite() {
                return Err(invalid(if unflanged {
                    "Levine-Schwinger radiation oracle returned unresolved resistance or reactance"
                } else {
                    "Rayleigh radiation oracle returned unresolved resistance or reactance"
                }));
            }
            for discrete in [false,true] {
                let (rr,xx)=result.impedance_at(frequency,discrete)?;
                let complex=(rr-r).hypot(xx-x)/r.hypot(x);
                let resistance=(rr-r).abs()/r;
                if !complex.is_finite() || !resistance.is_finite() || complex>0.05 || resistance>0.05 {
                    return Err(invalid(if unflanged {
                        "compact unflanged load exceeds the fixed 5% complex or resistance check"
                    } else {
                        "compact baffled load exceeds the fixed 5% complex or resistance check"
                    }));
                }
                result.max_complex_relative_error=result.max_complex_relative_error.max(complex);
                result.max_resistance_relative_error=result.max_resistance_relative_error.max(resistance);
            }
        }
        if gate.is_requested() {return Err(invalid("baffled radiation construction cancelled"));}
        Ok(result)
    }

    /// Existing passive terminal model; port storage/loss remain in fs-vfit.
    #[must_use]
    pub const fn termination(&self) -> NetworkNode {NetworkNode::Relaxation {load:self.load}}
    /// Complete coefficients in acoustic impedance units, not a signal gain.
    #[must_use]
    pub const fn load(&self) -> RelaxationImpedanceSpec {self.load}
    /// Exterior law of this open end (`Baffled` unless a card selected `Unflanged`).
    #[must_use]
    pub const fn flange(&self) -> RadiationFlange {self.flange}
    /// Actual uniform mouth radius [m].
    #[must_use]
    pub const fn radius_m(&self) -> f64 {self.radius_m}
    /// Homogeneous gas density [kg/m^3].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {self.density_kg_m3}
    /// Homogeneous sound speed [m/s].
    #[must_use]
    pub const fn sound_speed_m_s(&self) -> f64 {self.sound_speed_m_s}
    /// Mechanical step used by the checked bilinear response [s].
    #[must_use]
    pub const fn time_step_s(&self) -> f64 {self.time_step_s}
    /// Declared low-frequency use band, not an output low-pass filter.
    #[must_use]
    pub const fn maximum_frequency_hz(&self) -> f64 {self.maximum_frequency_hz}
    /// Largest checked complex discrepancy, including the discrete-time response.
    #[must_use]
    pub const fn max_complex_relative_error(&self) -> f64 {self.max_complex_relative_error}
    /// Largest checked resistance discrepancy, independently of stronger reactance.
    #[must_use]
    pub const fn max_resistance_relative_error(&self) -> f64 {self.max_resistance_relative_error}

    /// Analogue or actually realized bilinear (R,X), under exp(-i omega t).
    /// Querying beyond the admitted source band refuses rather than extrapolating
    /// a compact approximation as an accurate full-band radiation model.
    ///
    /// # Errors
    /// Nonfinite/negative/out-of-band frequency or unrepresentable response.
    pub fn impedance_at(&self, frequency_hz:f64, discrete:bool)->Result<(f64,f64),AcousticRealizeError> {
        if !frequency_hz.is_finite() || frequency_hz<0.0 || frequency_hz>self.maximum_frequency_hz {
            return Err(invalid("baffled load impedance query is outside its admitted frequency band"));
        }
        if frequency_hz==0.0 {return Ok((0.0,0.0));}
        let mut omega=core::f64::consts::TAU*frequency_hz;
        if discrete {
            let angle=omega*(0.5*self.time_step_s);
            omega=2.0*(fs_math::det::sin(angle)/fs_math::det::cos(angle))/self.time_step_s;
        }
        let term=self.load.terms()[0];
        let ratio=omega/term.rate_per_s;
        let x=term.resistance_pa_s_m3*ratio/(1.0+ratio*ratio);
        let r=x*ratio;
        if ![r,x].iter().all(|v|v.is_finite()) {return Err(invalid("baffled impedance response overflowed"));}
        Ok((r,-x))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// G3: coefficients derived from the two exact low-ka terms, then a
    /// fixed 512-point full-domain discrepancy measurement, without fitting.
    #[test]
    fn compact_load_two_term_exact_and_error_recorded() {
        let start = std::time::Instant::now();
        let c = core::f64::consts::TAU;
        let a = BaffledRadiationLoad::new(1.0, 1.0, c, 1e-4, 0.5, &CancelGate::new()).unwrap();
        let alpha = 8.0 / (3.0 * core::f64::consts::PI);
        let z0 = c / core::f64::consts::PI;
        let term = a.load.terms()[0];
        let inertance = term.resistance_pa_s_m3 / term.rate_per_s;
        let resistance_coefficient = term.resistance_pa_s_m3 / (term.rate_per_s * term.rate_per_s);
        assert!((inertance / (z0 * alpha / c) - 1.0).abs() < 1e-14);
        assert!((resistance_coefficient / (z0 * 0.5 / (c * c)) - 1.0).abs() < 1e-14);
        let (mut worst_complex, mut worst_r) = (0.0_f64, 0.0_f64);
        for i in 1..=512 {
            let ka = 0.5 * f64::from(i) / 512.0;
            let (r, x) = fs_phs::baffled_piston_impedance(1.0, c, 1.0, ka * c).unwrap();
            let (rr, xx) = a.impedance_at(ka, false).unwrap();
            worst_complex = worst_complex.max((rr - r).hypot(xx - x) / r.hypot(x));
            worst_r = worst_r.max((rr / r - 1.0).abs());
        }
        assert!(worst_complex < 0.022 && worst_r < 0.042);
        println!(
            "{{\"case\":\"compact_load_two_term_exact_and_error_recorded\",\"inputs\":\"analytic piston; ka in (0,0.5], 512 frequencies\",\"observed\":{worst_complex:e},\"resistance_error\":{worst_r:e},\"oracle\":\"Bessel-Struve\",\"tol\":0.022,\"resistance_tol\":0.042,\"verdict\":\"pass\",\"elapsed_ms\":{}}}",
            start.elapsed().as_millis()
        );
    }

    /// G0: exact existing domain admission; cap accepted and cap+1 refused.
    #[test]
    fn compact_baffled_load_domain_boundaries() {
        let make = |dt, f| {
            BaffledRadiationLoad::new(1.0, 1.0, core::f64::consts::TAU, dt, f, &CancelGate::new())
        };
        assert!(make(0.2, 0.5).is_ok()); // ka=0.5, dt*fmax=0.1.
        assert!(make(0.2, 0.5 * (1.0 + 1e-9)).is_err());
        assert!(make(0.2 * (1.0 + 1e-9), 0.5).is_err());
    }

    /// Unflanged load, test 1 (MB 92): the two low-ka terms are exact (end correction
    /// 0.6133 a, resistance (ka)^2/4 of Silva et al. 2009 Eqs. 10-11) and the checked
    /// discrepancy against fs-duct's exact Levine-Schwinger evaluation is recorded.
    /// Predicted before this code existed, from three of the exact LS rows printed in
    /// fs-duct (mpmath): 1.33 % complex and 0.88 % resistance at ka = 0.5, 0.10 % at ka 0.01
    /// (0.6133 against the exact 0.6127). The 512-point grid measured on ts2 (2026-11-04)
    /// gave 1.33 % complex (as predicted, at ka 0.5) and 1.74 % resistance, which the
    /// three-row prediction missed because the resistance error peaks at an interior ka.
    /// The bounds below (2 %) are regression bounds over that measured grid, not a physics
    /// claim; the engine's admission gate stays the 5 % check in `with_flange`.
    #[test]
    fn unflanged_load_two_term_exact_and_oracle_error_recorded() {
        use fs_duct::levine_schwinger::{LsQuadrature, levine_schwinger};
        let start = std::time::Instant::now();
        let c = core::f64::consts::TAU;
        let a = BaffledRadiationLoad::with_flange(1.0, 1.0, c, 1e-4, 0.5,
            RadiationFlange::Unflanged, &CancelGate::new()).unwrap();
        assert_eq!(a.flange(), RadiationFlange::Unflanged);
        let z0 = c / core::f64::consts::PI;
        let term = a.load.terms()[0];
        let inertance = term.resistance_pa_s_m3 / term.rate_per_s;
        let resistance_coefficient = term.resistance_pa_s_m3 / (term.rate_per_s * term.rate_per_s);
        assert!((inertance / (z0 * UNFLANGED_END_CORRECTION / c) - 1.0).abs() < 1e-14);
        assert!((resistance_coefficient / (z0 * 0.25 / (c * c)) - 1.0).abs() < 1e-14);
        let (mut worst_complex, mut worst_r, mut at_r) = (0.0_f64, 0.0_f64, 0.0_f64);
        for i in 1..=512 {
            let ka = 0.5 * f64::from(i) / 512.0;
            let exact = levine_schwinger(ka, LsQuadrature { nodes: 256 }).unwrap().normalized_impedance();
            let (r, x) = (exact.re * z0, exact.im * z0);
            let (rr, xx) = a.impedance_at(ka, false).unwrap();
            worst_complex = worst_complex.max((rr - r).hypot(xx - x) / r.hypot(x));
            let resistance = (rr / r - 1.0).abs();
            if resistance > worst_r { (worst_r, at_r) = (resistance, ka); }
        }
        assert!(worst_complex < 0.02 && worst_r < 0.02,
            "complex {worst_complex}, resistance {worst_r} at ka {at_r}");
        println!(
            "{{\"case\":\"unflanged_load_two_term_exact_and_oracle_error_recorded\",\"inputs\":\"fs-duct Levine-Schwinger, 256 nodes; ka in (0,0.5], 512 frequencies\",\"observed\":{worst_complex:e},\"resistance_error\":{worst_r:e},\"resistance_error_at_ka\":{at_r:e},\"oracle\":\"Levine-Schwinger\",\"tol\":0.02,\"resistance_tol\":0.02,\"verdict\":\"pass\",\"elapsed_ms\":{}}}",
            start.elapsed().as_millis()
        );
    }

    /// Unflanged load: the same domain admission as the baffled load (ka <= 0.5 and dt*fmax <= 0.1
    /// are refused one part in 1e9 beyond the cap). At the corner ka = 0.5 AND dt*fmax = 0.1 together
    /// the bilinear warping puts the unflanged resistance check within a fraction of a percent of its
    /// 5 % gate (Python port: 4.97 % against Silva's Pade reference), so admission THERE is not
    /// asserted either way; the accepted case uses dt*fmax = 0.01, the recorder's order.
    #[test]
    fn compact_unflanged_load_domain_boundaries() {
        let make = |dt, f| BaffledRadiationLoad::with_flange(1.0, 1.0, core::f64::consts::TAU, dt, f,
            RadiationFlange::Unflanged, &CancelGate::new());
        assert!(make(0.02, 0.5).is_ok()); // ka=0.5, dt*fmax=0.01.
        assert!(make(0.2, 0.5 * (1.0 + 1e-9)).is_err());
        assert!(make(0.2 * (1.0 + 1e-9), 0.5).is_err());
        assert!(BaffledRadiationLoad::with_flange(0.0, 1.0, 1.0, 1e-4, 0.5,
            RadiationFlange::Unflanged, &CancelGate::new()).is_err());
    }

    /// Absent selection keeps the original bytes; the planted negative (an unflanged
    /// load must differ from the baffled one, and by the half resistance, not a retune).
    #[test]
    fn baffled_flange_is_the_original_load_and_unflanged_differs() {
        let gate = CancelGate::new();
        // A recorder tone hole: radius 3 mm, breath-warm gas, 96 kHz mechanics, ka_max 0.16.
        let (a, rho, c, dt, fmax) = (0.003, 1.17, 345.36, 1.0 / 96_000.0, 3_000.0);
        let original = BaffledRadiationLoad::new(a, rho, c, dt, fmax, &gate).unwrap();
        let named = BaffledRadiationLoad::with_flange(a, rho, c, dt, fmax,
            RadiationFlange::Baffled, &gate).unwrap();
        assert_eq!(original, named);
        assert_eq!(original.flange(), RadiationFlange::Baffled);
        assert_eq!(RadiationFlange::default(), RadiationFlange::Baffled);
        let unflanged = BaffledRadiationLoad::with_flange(a, rho, c, dt, fmax,
            RadiationFlange::Unflanged, &gate).unwrap();
        assert_ne!(original, unflanged);
        assert_ne!(original.termination(), unflanged.termination());
        let coefficient = |load: &BaffledRadiationLoad| {
            let t = load.load.terms()[0];
            t.resistance_pa_s_m3 / (t.rate_per_s * t.rate_per_s)
        };
        // Low-ka radiation resistance: unflanged (ka)^2/4, baffled (ka)^2/2.
        assert!((coefficient(&unflanged) / coefficient(&original) - 0.5).abs() < 1e-12);
        // And a lower pole rate / longer inertance on the baffled side: 0.6133 a against 0.8488 a.
        let inertance = |load: &BaffledRadiationLoad| {
            let t = load.load.terms()[0];
            t.resistance_pa_s_m3 / t.rate_per_s
        };
        let ratio = inertance(&unflanged) / inertance(&original);
        assert!((ratio - UNFLANGED_END_CORRECTION / (8.0 / (3.0 * core::f64::consts::PI))).abs() < 1e-12);
    }

    /// Unflanged load, test 3 (MB 92): at low ka the option on the UNMODIFIED chimney equals the
    /// shortened-chimney emulation (baffled load, chimney shortened by (alpha - 0.6133) b) in
    /// reactance, hence in pitch. Predicted before the code (Python port, same formulas): reactance
    /// of the chimney input impedance differs by 6e-6 / 2.5e-5 / 1.6e-4 / 7.5e-4 at ka_hole
    /// 0.01 / 0.02 / 0.05 / 0.1, i.e. 0.003-0.010 cents at the passive first mode. The RESISTANCE
    /// differs by a factor 2 (the baffled resistance is twice the unflanged one), which is the
    /// planted negative: the emulation is not the same loss.
    #[test]
    fn unflanged_hole_matches_shortened_chimney_in_reactance_at_low_ka() {
        let gate = CancelGate::new();
        let (a, rho, c) = (0.002675, 1.17, 345.36);
        let (dt, fmax) = (1.0e-5, 0.2 * c / (core::f64::consts::TAU * a)); // ka_max 0.2, dt*fmax 0.04
        let z0 = rho * c / (core::f64::consts::PI * a * a);
        let option = BaffledRadiationLoad::with_flange(a, rho, c, dt, fmax,
            RadiationFlange::Unflanged, &gate).unwrap();
        let emulation = BaffledRadiationLoad::new(a, rho, c, dt, fmax, &gate).unwrap();
        let alpha = 8.0 / (3.0 * core::f64::consts::PI);
        let chimney = 0.0049 + 0.0070; // chimney plus inner length of the card's first hole [m]
        // exp(-i omega t): Z_in = (Z_L - i tan kl) / (1 - i Z_L tan kl), all divided by Zc.
        let z_in = |z_l: (f64, f64), tan_kl: f64| {
            let (n_re, n_im) = (z_l.0, z_l.1 - tan_kl);
            let (d_re, d_im) = (1.0 + z_l.1 * tan_kl, -z_l.0 * tan_kl);
            let d = d_re * d_re + d_im * d_im;
            ((n_re * d_re + n_im * d_im) / d, (n_im * d_re - n_re * d_im) / d)
        };
        for ka in [0.01, 0.02, 0.05, 0.1] {
            let frequency = ka * c / (core::f64::consts::TAU * a);
            let k = core::f64::consts::TAU * frequency / c;
            let (r1, x1) = option.impedance_at(frequency, false).unwrap();
            let (r2, x2) = emulation.impedance_at(frequency, false).unwrap();
            let a_in = z_in((r1 / z0, x1 / z0), (k * chimney).tan());
            let e_in = z_in((r2 / z0, x2 / z0), (k * (chimney - (alpha - UNFLANGED_END_CORRECTION) * a)).tan());
            let reactance = (a_in.1 - e_in.1).abs() / a_in.1.abs();
            let resistance = (a_in.0 - e_in.0).abs() / a_in.0.abs();
            assert!(reactance < 2.0e-3, "ka {ka}: reactance differs by {reactance}");
            assert!(resistance > 0.9, "ka {ka}: resistance differs by only {resistance}");
        }
    }
}
