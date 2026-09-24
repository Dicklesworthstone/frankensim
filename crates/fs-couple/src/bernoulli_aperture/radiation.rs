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

use crate::acoustic_realize::AcousticRealizeError;
use super::network::NetworkNode;
use fs_exec::CancelGate;
use fs_vfit::impedance::SeriesImpedanceSpec;
use fs_vfit::relaxation::{RelaxationImpedanceSpec, RelaxationTerm};

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
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

impl BaffledRadiationLoad {
    /// Derive the sole relaxation branch from the ACTUAL mouth radius and fluid.
    /// All 32 positive frequencies spanning (0, maximum_frequency_hz] must meet
    /// fixed 5% relative complex AND radiation-resistance checks for the analogue
    /// and discrete responses. The analytical zero-frequency limits are exact.
    /// The disk oracle uses its existing eight-ring positive-area quadrature;
    /// that oracle's spatial error is not removed from the discrepancy report.
    ///
    /// `dt*fmax <= 0.1` limits temporal warping and `ka_max <= 0.5` declares the
    /// low-frequency model domain. They do not prove that an arbitrary transient
    /// has no higher-frequency content. No coefficients are retuned after checks.
    ///
    /// # Errors
    /// Invalid/unrepresentable geometry or medium, cancellation, an unresolved
    /// clock, oracle failure, or either checked discrepancy above 5%.
    pub fn new(radius_m: f64, density_kg_m3: f64, sound_speed_m_s: f64,
        time_step_s: f64, maximum_frequency_hz: f64, gate: &CancelGate)
        -> Result<Self, AcousticRealizeError>
    {
        if ![radius_m,density_kg_m3,sound_speed_m_s,time_step_s,maximum_frequency_hz]
            .iter().all(|v| v.is_finite() && *v > 0.0)
        { return Err(invalid("baffled radiation requires positive finite geometry, fluid, clock and band")); }
        let ka=core::f64::consts::TAU*maximum_frequency_hz/sound_speed_m_s*radius_m;
        if !ka.is_finite() || ka>0.5 || time_step_s*maximum_frequency_hz>0.1 {
            return Err(invalid("compact baffled radiation requires ka <= 0.5 and dt*fmax <= 0.1"));
        }
        let area=core::f64::consts::PI*radius_m*radius_m;
        let z=density_kg_m3*sound_speed_m_s/area;
        let alpha=8.0/(3.0*core::f64::consts::PI);
        let term=RelaxationTerm {resistance_pa_s_m3:2.0*alpha*alpha*z,
            rate_per_s:2.0*alpha*sound_speed_m_s/radius_m};
        let load=RelaxationImpedanceSpec::new(SeriesImpedanceSpec {
            resistance_pa_s_m3:0.0,inertance_pa_s2_m3:0.0,compliance_m3_pa:None,
        }, &[term]).map_err(|e|AcousticRealizeError::Nonlinear(e.to_string()))?;
        // The existing owner also checks representability of its stored energy
        // and discrete coefficients before any physical system is constructed.
        fs_vfit::relaxation::RelaxationImpedance::new(load,z,time_step_s)
            .map_err(|e|AcousticRealizeError::Nonlinear(e.to_string()))?;
        let mut result=Self {radius_m,density_kg_m3,sound_speed_m_s,time_step_s,
            maximum_frequency_hz,max_complex_relative_error:0.0,max_resistance_relative_error:0.0,load};
        for i in 1..=32 {
            if gate.is_requested() {return Err(invalid("baffled radiation construction cancelled"));}
            let frequency=maximum_frequency_hz*f64::from(i)/32.0;
            let omega=core::f64::consts::TAU*frequency;
            let (r,x)=fs_phs::baffled_piston_impedance(density_kg_m3,sound_speed_m_s,radius_m,omega,8)
                .map_err(|e|AcousticRealizeError::Nonlinear(e.to_string()))?;
            if !r.is_finite() || r<=0.0 || !x.is_finite() {
                return Err(invalid("Rayleigh radiation oracle returned unresolved resistance or reactance"));
            }
            for discrete in [false,true] {
                let (rr,xx)=result.impedance_at(frequency,discrete)?;
                let complex=(rr-r).hypot(xx-x)/r.hypot(x);
                let resistance=(rr-r).abs()/r;
                if !complex.is_finite() || !resistance.is_finite() || complex>0.05 || resistance>0.05 {
                    return Err(invalid("compact baffled load exceeds the fixed 5% complex or resistance check"));
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
