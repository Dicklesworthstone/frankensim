//! Explicit conversion of planar plastic strain ratios to the native Hill metric.
//!
//! This fixes the rolling-direction yield normalization, not an independent
//! least-squares fit to both directional stresses and plastic strain ratios.

use super::{Hill48Coefficients, HillError, Modes, positive};

impl Hill48Coefficients {
    /// Calibrate planar anisotropy from Lankford plastic width/thickness ratios
    /// `[r_0, r_45, r_90]`, measured at 0, 45 and 90 degrees to rolling.
    /// The reference yield stress MUST be the rolling-direction uniaxial yield
    /// stress: this convention fixes `G + H = 1`.
    ///
    /// The remaining out-of-plane shear yield ratios `[R_23, R_31]` are REQUIRED
    /// independent data, relative to `reference_yield / sqrt(3)`. Planar r-values
    /// cannot identify them. No transverse shear isotropy is silently assumed.
    /// All five ratios must be finite and strictly positive.
    ///
    /// r-values are ratios of PLASTIC strain increments, not total strain,
    /// engineering elongation, fracture strain or a scalar malleability rating.
    /// The quadratic associated-flow model then predicts other directional
    /// yield stresses; independently measured yields may not fit the same
    /// surface. Coefficients and orientation remain fixed during a solve.
    ///
    /// Reference: Abaqus, Hill anisotropic yield, the planar/general strain-ratio
    /// conversion: `r_0=H/G`, `r_90=H/F`, `r_45=N/(F+G)-1/2`.
    /// <https://docs.software.vt.edu/abaqusv2025/English/SIMACAEMATRefMap/simamat-c-anisoyield.htm>
    pub fn from_planar_r_values(
        planar_r: [f64; 3],
        transverse_shear_yield_ratios: [f64; 2],
    ) -> Result<Self, HillError> {
        if planar_r
            .iter()
            .chain(transverse_shear_yield_ratios.iter())
            .any(|value| !positive(*value))
        {
            return Err(HillError::Parameters(
                "planar plastic strain ratios and both transverse shear yield ratios must be finite and positive",
            ));
        }
        let [r_0, r_45, r_90] = planar_r;
        // r/(1+r) avoids losing H when r is small; do not use H=1-G.
        let g = 1.0 / (1.0 + r_0);
        let h = r_0 / (1.0 + r_0);
        let f = h / r_90;
        let value = Self {
            f,
            g,
            h,
            l: 1.5 * (1.0 / transverse_shear_yield_ratios[0]).powi(2),
            m: 1.5 * (1.0 / transverse_shear_yield_ratios[1]).powi(2),
            n: (r_45 + 0.5) * (f + g),
        };
        // Convexity alone does not guarantee the calibration is retained: F
        // can underflow while G/H still admit a metric, silently losing r90.
        // Likewise, positive r45 must not round to zero width contraction.
        if [value.f, value.g, value.h, value.l, value.m, value.n]
            .iter()
            .any(|coefficient| !positive(*coefficient))
            || value.n <= 0.5 * (value.f + value.g)
        {
            return Err(HillError::Parameters(
                "planar calibration is not representable",
            ));
        }
        Modes::new(value)?;
        Ok(value)
    }
}
