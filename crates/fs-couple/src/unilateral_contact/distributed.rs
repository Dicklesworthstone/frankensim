//! Preserve the existing unit-slit arithmetic while admitting distributed scalar
//! contact through fs-dcontact's pointwise secant and local unloading clamps.
use super::{DContactError, Obstacle, SlitContactStep};
use fs_dcontact::{AffineContactResponse, AffineContactStep};

/// Explicit bound on quadrature work in one scalar aperture trial.
pub(crate) const MAX_APERTURE_CONTACT_POINTS: usize = 4096;

pub(crate) enum ApertureContactStep<'a> {
    Unit(SlitContactStep),
    Distributed(AffineContactStep<'a>),
}
impl<'a> ApertureContactStep<'a> {
    pub(crate) fn new(law: &'a Obstacle, before: f64) -> Result<Self, DContactError> {
        let distributed = AffineContactStep::new(law, before, MAX_APERTURE_CONTACT_POINTS)?;
        if law.n_points() == 1 && law.collocation() == [-1.0] {
            Ok(Self::Unit(SlitContactStep::new(law,before)?))
        } else {
            Ok(Self::Distributed(distributed))
        }
    }
    pub(crate) fn response(&self, after: f64, velocity: f64) -> Result<AffineContactResponse, DContactError> {
        match self {
            Self::Distributed(law) => law.response(after,velocity),
            Self::Unit(law) => {
                let (elastic, damping) = law.coefficients(after)?;
                let force = (elastic-damping*velocity).max(0.0);
                let power = (elastic-force)*velocity;
                if !velocity.is_finite() || !force.is_finite() || !power.is_finite() {
                    return Err(DContactError::Parameter { what: "unit slit response overflowed" });
                }
                Ok(AffineContactResponse { elastic_force: elastic, force, dissipated_power: power })
            }
        }
    }
}
