//! Small-signal Helmholtz cavity termination for the existing tube network.
//!
//! L = rho l_eff / S and C = V / (rho c^2), with explicitly supplied R.
//! The effective neck length includes the caller's chosen end correction; no
//! hidden correction, radiation resistance or material-name inference occurs.
//! See UNSW, "Helmholtz Resonance", and Scientific Reports 10 (2020),
//! doi:10.1038/s41598-020-67608-z, equations (1)-(2).
//!
//! This is the lumped, adiabatic, linear cavity/plug-neck limit. The caller must
//! justify uniform cavity pressure and a neck short compared with wavelength;
//! volume alone cannot establish that frequency range for arbitrary shapes.
//! Distributed viscothermal losses, radiation, moving walls, nonlinear jet loss
//! and experimental identification are NOT inferred from these coefficients.

use super::network::NetworkNode;
use crate::acoustic_realize::AcousticRealizeError;
pub use fs_vfit::impedance::SeriesImpedanceSpec;

/// Explicit physical geometry and measured/declared small-signal resistance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HelmholtzLoadSpec {
    /// Enclosed cavity volume [m^3].
    pub volume_m3: f64,
    /// Circular neck radius [m].
    pub neck_radius_m: f64,
    /// Effective acoustic neck length [m], INCLUDING all intended end corrections.
    /// The same neck must not also be included as a propagating TubeSection.
    pub effective_neck_length_m: f64,
    /// Constant series resistance [Pa s/m^3], nonnegative and explicitly supplied.
    /// This is not an automatically predicted broadband viscothermal resistance.
    pub resistance_pa_s_m3: f64,
}

impl HelmholtzLoadSpec {
    /// Derive the acoustic series impedance from the declared fluid and geometry.
    ///
    /// # Errors
    /// Nonpositive/nonfinite geometry/fluid, negative resistance or overflow.
    pub fn impedance(&self, density_kg_m3: f64, sound_speed_m_s: f64)
        -> Result<SeriesImpedanceSpec, AcousticRealizeError>
    {
        let invalid = || AcousticRealizeError::InvalidDescription {
            what: "Helmholtz load requires positive finite geometry/fluid and nonnegative finite resistance",
        };
        if ![self.volume_m3, self.neck_radius_m, self.effective_neck_length_m,
            density_kg_m3, sound_speed_m_s].iter().all(|x| x.is_finite() && *x > 0.0)
            || !self.resistance_pa_s_m3.is_finite() || self.resistance_pa_s_m3 < 0.0
        {
            return Err(invalid());
        }
        let area = core::f64::consts::PI * self.neck_radius_m * self.neck_radius_m;
        let inertance = density_kg_m3 * self.effective_neck_length_m / area;
        let compliance = ((self.volume_m3 / density_kg_m3) / sound_speed_m_s) / sound_speed_m_s;
        if ![area, inertance, compliance].iter().all(|x| x.is_finite() && *x > 0.0) {
            return Err(invalid());
        }
        Ok(SeriesImpedanceSpec {
            resistance_pa_s_m3: self.resistance_pa_s_m3,
            inertance_pa_s2_m3: inertance,
            compliance_m3_pa: Some(compliance),
        })
    }

    /// Build the reactive terminal used directly by TubeNetworkSpec::nodes.
    /// The network owns its persistent pressure/flow state and energy budget.
    ///
    /// # Errors
    /// Same physical admission as `impedance`.
    pub fn termination(&self, density_kg_m3: f64, sound_speed_m_s: f64)
        -> Result<NetworkNode, AcousticRealizeError>
    {
        Ok(NetworkNode::Impedance { load: self.impedance(density_kg_m3, sound_speed_m_s)? })
    }

    /// Undamped continuous-model resonance [Hz], not the loaded network's peak.
    /// The runtime additionally has explicit bilinear frequency warping.
    ///
    /// # Errors
    /// Invalid physical parameters or an unrepresentable resonance.
    pub fn resonance_hz(&self, density_kg_m3: f64, sound_speed_m_s: f64)
        -> Result<f64, AcousticRealizeError>
    {
        let load = self.impedance(density_kg_m3, sound_speed_m_s)?;
        let product = load.inertance_pa_s2_m3 * load.compliance_m3_pa.expect("cavity compliance");
        let hz = 1.0 / (2.0 * core::f64::consts::PI * fs_math::det::sqrt(product));
        if !hz.is_finite() || hz <= 0.0 {
            return Err(AcousticRealizeError::InvalidDescription { what: "Helmholtz resonance is not representable" });
        }
        Ok(hz)
    }
}
