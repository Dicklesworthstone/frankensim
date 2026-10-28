//! Unit-consistent, allocation-free SGS closures independent of grid storage.

use fs_math::det;

/// An explicitly selected closure and its declared, nonnegative coefficient.
/// Zero is the exact disabled control for refinement and same-trajectory tests.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SgsModel {
    Smagorinsky { cs: f64 },
    Wale { cw: f64 },
}

impl SgsModel {
    pub fn validate(self) -> Result<(), SgsError> {
        let coefficient = match self {
            Self::Smagorinsky { cs } => cs,
            Self::Wale { cw } => cw,
        };
        if coefficient.is_finite() && coefficient >= 0.0 {
            Ok(())
        } else {
            Err(SgsError::InvalidCoefficient)
        }
    }
}

/// Closures accelerate exploration; they cannot produce acceptance renders.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SgsUse {
    Exploration,
    AcceptanceRender,
}

/// Lattice-unit controls. `filter_width_lu` is explicit (usually one cell).
/// `maximum_eddy_to_molecular` is a refusal boundary, never a viscosity clamp.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SgsConfig {
    pub model: SgsModel,
    pub filter_width_lu: f64,
    pub maximum_eddy_to_molecular: f64,
    pub usage: SgsUse,
}

impl SgsConfig {
    pub fn validate(self) -> Result<(), SgsError> {
        if self.usage == SgsUse::AcceptanceRender {
            return Err(SgsError::AcceptanceRenderForbidden);
        }
        self.model.validate()?;
        if !(self.filter_width_lu.is_finite() && self.filter_width_lu > 0.0) {
            return Err(SgsError::InvalidFilterWidth);
        }
        if !(self.maximum_eddy_to_molecular.is_finite() && self.maximum_eddy_to_molecular >= 0.0) {
            return Err(SgsError::InvalidViscosityRatio);
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SgsError {
    InvalidCoefficient,
    InvalidFilterWidth,
    InvalidViscosityRatio,
    NonFiniteGradient,
    NonFiniteViscosity,
    InvalidMolecularViscosity {
        cell: usize,
    },
    AcceptanceRenderForbidden,
    ViscosityRatioExceeded {
        cell: usize,
        ratio: f64,
        maximum: f64,
    },
}

impl core::fmt::Display for SgsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-SGS3: {self:?}")
    }
}
impl std::error::Error for SgsError {}

/// Evaluate `nu_t` from `gradient[i][j] = d u_i / d x_j`.
/// Length and gradient units must agree: metres and inverse seconds return
/// m^2/s; lattice lengths and gradients return lattice kinematic viscosity.
/// Fixed-size, allocation-free arithmetic with normalized invariants avoids
/// overflow/underflow in WALE's fourth- and tenth-degree intermediate powers.
pub fn eddy_viscosity(
    model: SgsModel,
    filter_width: f64,
    gradient: [[f64; 3]; 3],
) -> Result<f64, SgsError> {
    model.validate()?;
    if !(filter_width.is_finite() && filter_width > 0.0) {
        return Err(SgsError::InvalidFilterWidth);
    }
    let mut scale: f64 = 0.0;
    for row in gradient {
        for value in row {
            if !value.is_finite() {
                return Err(SgsError::NonFiniteGradient);
            }
            scale = scale.max(value.abs());
        }
    }
    let coefficient = match model {
        SgsModel::Smagorinsky { cs } => cs,
        SgsModel::Wale { cw } => cw,
    };
    if scale == 0.0 || coefficient == 0.0 {
        return Ok(0.0);
    }
    let g = gradient.map(|row| row.map(|value| value / scale));
    let mut strain2 = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            let s = 0.5 * (g[i][j] + g[j][i]);
            strain2 += s * s;
        }
    }
    let rate = match model {
        SgsModel::Smagorinsky { .. } => det::sqrt(2.0 * strain2),
        SgsModel::Wale { .. } => {
            let square: [[f64; 3]; 3] = core::array::from_fn(|i| {
                core::array::from_fn(|j| (0..3).map(|k| g[i][k] * g[k][j]).sum())
            });
            let trace_third = (square[0][0] + square[1][1] + square[2][2]) / 3.0;
            let mut deviator2 = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    let d = 0.5 * (square[i][j] + square[j][i])
                        - if i == j { trace_third } else { 0.0 };
                    deviator2 += d * d;
                }
            }
            if deviator2 == 0.0 {
                return Ok(0.0);
            }
            let root = det::sqrt(deviator2);
            deviator2 * root
                / (strain2 * strain2 * det::sqrt(strain2) + deviator2 * det::sqrt(root))
        }
    };
    let length = coefficient * filter_width;
    let viscosity = length * (length * (scale * rate));
    if viscosity.is_finite() {
        Ok(viscosity)
    } else {
        Err(SgsError::NonFiniteViscosity)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn closures_obey_shear_wall_and_rotation_laws() {
        let smag = SgsModel::Smagorinsky { cs: 0.17 };
        let wale = SgsModel::Wale { cw: 0.5 };
        let shear = [[0.0, 0.3, 0.0], [0.0; 3], [0.0; 3]];
        let expected = (0.17_f64 * 2.0).powi(2) * 0.3;
        assert!((eddy_viscosity(smag, 2.0, shear).unwrap() - expected).abs() < 1e-15);
        // G^2, not G G^T: WALE vanishes in simple laminar shear.
        assert_eq!(eddy_viscosity(wale, 2.0, shear).unwrap(), 0.0);
        let rotation = [[0.0, -0.2, 0.0], [0.2, 0.0, 0.0], [0.0; 3]];
        assert_eq!(eddy_viscosity(smag, 1.0, rotation).unwrap(), 0.0);
        for model in [smag, wale] {
            assert_eq!(eddy_viscosity(model, 1.0, [[0.0; 3]; 3]).unwrap(), 0.0);
        }
        // Incompressible smooth 3-D no-slip field: WALE's near-wall cubic
        // asymptote must emerge without a wall distance or damping knob.
        let wall = |y: f64| {
            [
                [0.2 * y, 1.0, 0.0],
                [0.0, -0.5 * y, 0.0],
                [0.0, 0.4, 0.3 * y],
            ]
        };
        let viscosity = [0.01, 0.005, 0.0025].map(|y| eddy_viscosity(wale, 1.0, wall(y)).unwrap());
        for pair in viscosity.windows(2) {
            let order = (pair[0] / pair[1]).ln() / 2.0_f64.ln();
            assert!((order - 3.0).abs() < 0.01, "near-wall order {order}");
        }
    }

    #[test]
    fn invariants_are_frame_invariant_and_scale_safe() {
        let g = [[0.3, 0.7, -0.2], [-0.1, -0.4, 0.6], [0.8, -0.5, 0.1]];
        let r = [[0.6, -0.8, 0.0], [0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];
        let rotated = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                (0..3)
                    .flat_map(|a| (0..3).map(move |b| r[i][a] * g[a][b] * r[j][b]))
                    .sum()
            })
        });
        for model in [
            SgsModel::Smagorinsky { cs: 0.17 },
            SgsModel::Wale { cw: 0.5 },
        ] {
            let reference = eddy_viscosity(model, 1.0, g).unwrap();
            assert!((eddy_viscosity(model, 1.0, rotated).unwrap() / reference - 1.0).abs() < 1e-14);
            assert!((eddy_viscosity(model, 0.5, g).unwrap() / reference - 0.25).abs() < 1e-15);
            for scale in [1e-150, 1e150] {
                let scaled = g.map(|row| row.map(|v| v * scale));
                let ratio = eddy_viscosity(model, 1.0, scaled).unwrap() / (scale * reference);
                assert!((ratio - 1.0).abs() < 1e-14, "scale={scale} ratio={ratio}");
            }
        }
    }

    #[test]
    fn invalid_inputs_and_acceptance_render_refuse() {
        let base = SgsConfig {
            model: SgsModel::Wale { cw: 0.5 },
            filter_width_lu: 1.0,
            maximum_eddy_to_molecular: 1.0,
            usage: SgsUse::Exploration,
        };
        assert_eq!(base.validate(), Ok(()));
        assert_eq!(
            SgsConfig {
                usage: SgsUse::AcceptanceRender,
                ..base
            }
            .validate(),
            Err(SgsError::AcceptanceRenderForbidden)
        );
        for coefficient in [-1.0, f64::NAN, f64::INFINITY] {
            for model in [
                SgsModel::Smagorinsky { cs: coefficient },
                SgsModel::Wale { cw: coefficient },
            ] {
                assert_eq!(
                    eddy_viscosity(model, 1.0, [[0.0; 3]; 3]),
                    Err(SgsError::InvalidCoefficient)
                );
            }
        }
        for width in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert_eq!(
                eddy_viscosity(base.model, width, [[0.0; 3]; 3]),
                Err(SgsError::InvalidFilterWidth)
            );
        }
        assert_eq!(
            eddy_viscosity(base.model, 1.0, [[f64::NAN; 3]; 3]),
            Err(SgsError::NonFiniteGradient)
        );
    }
}
