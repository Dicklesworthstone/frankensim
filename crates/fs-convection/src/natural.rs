//! Geometry and signed buoyancy select a natural-convection card, never a sum
//! of independently fitted heat-transfer coefficients.

use super::{CorrelationError, CorrelationId, CorrelationInputs, NusseltEvaluation, evaluate};

/// Which side of a horizontal plate is exposed, relative to gravity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HorizontalPlateFacing {
    /// Outward normal points opposite gravity.
    Up,
    /// Outward normal points along gravity.
    Down,
}

/// Geometry for isolated, isothermal natural convection in a quiescent fluid.
/// The caller must supply the geometry's characteristic length consistently
/// in both Rayleigh number and `h = Nu k/L`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NaturalConvectionGeometry {
    /// Churchill-Chu vertical plate; characteristic length is height.
    VerticalPlate,
    /// Churchill-Chu long horizontal cylinder; characteristic length is diameter.
    HorizontalCylinder,
    /// Horizontal plate; characteristic length is full plate area/perimeter.
    /// The face and signed buoyancy determine whether the boundary layer is
    /// stable or unstable. Heat still flows from the hotter body to the colder.
    HorizontalPlate(HorizontalPlateFacing),
    /// Churchill isolated isothermal sphere; characteristic length is diameter.
    /// Includes the steady unbounded-fluid conductive limit Nu=2 at Ra=0.
    Sphere,
}

impl NaturalConvectionGeometry {
    /// Evaluate the appropriate existing source card using
    /// `signed_rayleigh = g beta (T_wall - T_ambient) L^3 / (nu alpha)`.
    /// Gravity is a magnitude; signed expansion is retained. For a horizontal
    /// face it is the sign of density departure, not merely hot/cold wording,
    /// that selects stable versus unstable buoyancy.
    ///
    /// The unstable plate uses T9.7 through Ra=8e6 and T9.8 above it. Those
    /// empirical branches have a roughly 4.46% jump, which is retained rather
    /// than hidden by an unsourced blend. Smooth derivatives across this
    /// switch, extrapolation below the cards, confinement, and mixed forced/
    /// natural convection are not supplied by this selector.
    ///
    /// # Errors
    /// Refuses nonfinite groups and out-of-card ranges through the ordinary
    /// evaluator. Zero driving is admitted only for the sphere's physical
    /// steady conduction limit Nu=2. Other geometries leave exact thermal
    /// equilibrium to the coupling owner's zero-heat treatment.
    pub fn evaluate(
        self,
        signed_rayleigh: f64,
        prandtl: f64,
    ) -> Result<NusseltEvaluation, CorrelationError> {
        let magnitude = signed_rayleigh.abs();
        let id = match self {
            Self::VerticalPlate => CorrelationId::ChurchillChuVerticalPlate,
            Self::HorizontalCylinder => CorrelationId::ChurchillChuHorizontalCylinder,
            Self::Sphere => CorrelationId::ChurchillSphere,
            Self::HorizontalPlate(facing) => {
                let unstable = match facing {
                    HorizontalPlateFacing::Up => signed_rayleigh > 0.0,
                    HorizontalPlateFacing::Down => signed_rayleigh < 0.0,
                };
                if !unstable {
                    CorrelationId::HorizontalPlateStable
                } else if magnitude <= 8.0e6 {
                    CorrelationId::HorizontalPlateUnstableLaminar
                } else {
                    CorrelationId::HorizontalPlateUnstableTurbulent
                }
            }
        };
        evaluate(id, CorrelationInputs::natural(magnitude, prandtl))
    }
}
