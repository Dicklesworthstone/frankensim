//! Source-integrated molecular Fourier flux in the original NSF energy equation.
//!
//! A face holds its actual pressure and other source coordinates fixed. Its
//! Kirchhoff primitive is evaluated by the original WMS quadrature; this is
//! an estimated transport discretization, not a new caloric or energy owner.

use super::{NsfError, NsfGasMedium, NsfLaw, gas_error};
use core::cell::Cell;
use fs_material::thermochemical_gas::conductivity::{
    GasConductivityError, GasConductivityQuadrature, PreparedThermochemicalConductivity,
};

/// Optional original-source integration of molecular conductivity in NSF.
///
/// The default remains the historical midpoint coefficient. This option
/// integrates the retained pure-species transport sources and original WMS
/// mixture at the actual face pressure. Quadrature tolerances bound an
/// embedded numerical estimate in W/(m K), not an enclosure or joule defect.
/// The source's fixed-composition and material-support restrictions remain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfSourceConductivityControl {
    pub quadrature: GasConductivityQuadrature,
    /// Actual WMS point evaluations shared by every integral in one face
    /// stencil or one complete film/characteristic temperature root.
    /// Independent face budgets keep admission independent of worker count.
    pub maximum_evaluations_per_face: usize,
}

impl NsfSourceConductivityControl {
    pub(super) fn validate(self) -> Result<(), NsfError> {
        if !self.quadrature.validate() || self.maximum_evaluations_per_face == 0 {
            return Err(NsfError::Input("invalid NSF source conductivity controls"));
        }
        Ok(())
    }
}

enum FaceSource<'a> {
    Midpoint {
        gas: &'a NsfGasMedium,
        pressure_pa: f64,
    },
    // A genuinely isothermal face needs only the original point law. No
    // fictitious epsilon interval or extrapolation is introduced.
    Point {
        theta: f64,
        conductivity_lu: f64,
    },
    Integral {
        source: PreparedThermochemicalConductivity,
        controls: NsfSourceConductivityControl,
        evaluations: Cell<usize>,
    },
}

/// A short-lived face owner. Source selection happens once, including across
/// every candidate surface temperature of a nonlinear boundary root.
pub(in crate::nsf3) struct NsfFaceConductivity<'a> {
    law: &'a NsfLaw,
    source: FaceSource<'a>,
}

impl NsfGasMedium {
    pub(in crate::nsf3) fn prepare_face_conductivity<'a>(
        &'a self,
        law: &'a NsfLaw,
        pressure_pa: f64,
        theta_range: [f64; 2],
    ) -> Result<NsfFaceConductivity<'a>, NsfError> {
        let NsfGasMedium::Thermochemical(gas) = self else {
            return Ok(NsfFaceConductivity {
                law,
                source: FaceSource::Midpoint {
                    gas: self,
                    pressure_pa,
                },
            });
        };
        let Some(controls) = gas.source_conductivity else {
            return Ok(NsfFaceConductivity {
                law,
                source: FaceSource::Midpoint {
                    gas: self,
                    pressure_pa,
                },
            });
        };
        if !(pressure_pa.is_finite() && pressure_pa > 0.0)
            || theta_range.iter().any(|t| !t.is_finite() || *t <= 0.0)
        {
            return Err(NsfError::Input("invalid NSF conductivity face state"));
        }
        let lower = theta_range[0].min(theta_range[1]);
        let upper = theta_range[0].max(theta_range[1]);
        if lower == upper {
            let conductivity_lu = self.face_conductivity(law, lower, pressure_pa)?;
            return Ok(NsfFaceConductivity {
                law,
                source: FaceSource::Point {
                    theta: lower,
                    conductivity_lu,
                },
            });
        }
        let context = gas
            .context
            .clone()
            .with("pressure", pressure_pa)
            .map_err(gas_error)?;
        let source = PreparedThermochemicalConductivity::new(
            &gas.model,
            &context,
            [law.kelvin_of_theta(lower), law.kelvin_of_theta(upper)],
            || Ok(()),
        )
        .map_err(gas_error)?;
        Ok(NsfFaceConductivity {
            law,
            source: FaceSource::Integral {
                source,
                controls,
                evaluations: Cell::new(0),
            },
        })
    }
}

impl NsfFaceConductivity<'_> {
    /// The legacy query receives the caller's original midpoint expression
    /// verbatim, retaining its arithmetic even in a nonlinear film root.
    pub(in crate::nsf3) fn mean_lattice(
        &self,
        from: f64,
        to: f64,
        midpoint: f64,
    ) -> Result<f64, NsfError> {
        match &self.source {
            FaceSource::Midpoint { gas, pressure_pa } => {
                gas.face_conductivity(self.law, midpoint, *pressure_pa)
            }
            FaceSource::Point {
                theta,
                conductivity_lu,
            } => {
                if from != *theta || to != *theta {
                    return Err(NsfError::Input("NSF point conductivity interval changed"));
                }
                Ok(*conductivity_lu)
            }
            FaceSource::Integral {
                source,
                controls,
                evaluations,
            } => {
                let mean = source
                    .mean_conductivity(
                        self.law.kelvin_of_theta(from),
                        self.law.kelvin_of_theta(to),
                        controls.quadrature,
                        || {
                            let used = evaluations.get();
                            if used >= controls.maximum_evaluations_per_face {
                                return Err(GasConductivityError::WorkLimit);
                            }
                            evaluations.set(used + 1);
                            Ok(())
                        },
                    )
                    .map_err(gas_error)?;
                let value = mean.mean_w_m_k * self.law.transport_scale() / self.law.gas_constant;
                if !(value.is_finite() && value > 0.0) {
                    return Err(NsfError::Input("unrepresentable NSF mean conductivity"));
                }
                Ok(value)
            }
        }
    }

    /// Positive-axis conductivity times temperature gradient. The source
    /// option applies the native wide stencil to K(T; p_face), where K'=k.
    /// Multiplying a secant conductivity by the wide T stencil would not
    /// differentiate that same nonlinear primitive.
    pub(in crate::nsf3) fn gradient_lattice(
        &self,
        near: [f64; 2],
        wide: Option<[f64; 2]>,
        near_mean: f64,
        turbulent: f64,
    ) -> Result<f64, NsfError> {
        let difference = near[1] - near[0];
        let gradient = match wide {
            None => (near_mean + turbulent) * difference,
            Some(far) => {
                let far_difference = far[1] - far[0];
                let dt4 = (27.0 * difference - far_difference) / 24.0;
                if matches!(&self.source, FaceSource::Midpoint { .. }) {
                    // Preserve the original coefficient sum and product.
                    (near_mean + turbulent) * dt4
                } else {
                    let far_mean = self.mean_lattice(far[0], far[1], 0.5 * (far[0] + far[1]))?;
                    let molecular =
                        (27.0 * (near_mean * difference) - far_mean * far_difference) / 24.0;
                    molecular + turbulent * dt4
                }
            }
        };
        if !gradient.is_finite() {
            return Err(NsfError::Input("unrepresentable NSF conductivity gradient"));
        }
        Ok(gradient)
    }
}

#[cfg(test)]
mod tests;
