//! Exact elastic Hertz work and mean force along an affine indentation path.
//!
//! The compiled coefficient retains the original finite-patch geometry and
//! material receipts. Positive polynomial evaluation avoids subtracting nearly
//! equal endpoint energies, including at first contact and during unloading.
use super::SourcePointContact;
use crate::interface_binding::BoundNormalContactLaw;
use crate::normal_patch::NormalPatchError;

/// Allocation-free normal response along an admitted elastic indentation path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceElasticPathAverage {
    /// Time-averaged compressive force for affine indentation [N].
    pub mean_force_n: f64,
    /// Derivative of the mean force with respect to the final indentation [N/m].
    /// The initial indentation is held fixed; at equal endpoints this is half
    /// the instantaneous Hertz tangent, not the full instantaneous tangent.
    pub endpoint_derivative_n_per_m: f64,
    /// Reversible normal work, positive on loading and negative on unloading [J].
    pub work_j: f64,
}

impl SourcePointContact {
    /// Integrate the original elastic Hertz force over a straight indentation
    /// path without quadrature, allocation, or endpoint-energy cancellation.
    ///
    /// Both endpoints and their constant rate must lie inside the source-owned
    /// operating box. A Hunt--Crossley source is refused, even if its current
    /// dissipation coefficient is zero: this API represents the elastic model.
    /// The result is suitable for a discrete-gradient normal momentum solve;
    /// the caller still owns geometry, native mobility, and time integration.
    pub fn elastic_path_average(
        &self,
        initial_indentation_m: f64,
        final_indentation_m: f64,
        dt_s: f64,
    ) -> Result<SourceElasticPathAverage, NormalPatchError> {
        if self.model().law() != BoundNormalContactLaw::ElasticHertz {
            return Err(NormalPatchError::InvalidInput {
                field: "elastic source path requires ElasticHertz",
            });
        }
        self.validate_kinematics(initial_indentation_m, 0.0, dt_s)?;
        let increment = final_indentation_m - initial_indentation_m;
        let rate = increment / dt_s;
        self.validate_kinematics(initial_indentation_m, rate, dt_s)?;
        self.validate_kinematics(final_indentation_m, rate, dt_s)?;
        self.path_average(initial_indentation_m, final_indentation_m)
    }

    /// Mean force of the original reversible Hertz storage, without assigning
    /// a physical duration or rate to the path. This applies to both the elastic
    /// Hertz and Hunt--Crossley owners: their conservative potential is the
    /// same, while the latter's irreversible force remains a separate port.
    ///
    /// Only indentation is admitted here. A mechanical consumer must still
    /// check its actual rate and duration through `validate_kinematics`; this
    /// storage integral cannot establish dynamic applicability or supply loss.
    pub fn conservative_path_average(
        &self,
        initial_indentation_m: f64,
        final_indentation_m: f64,
    ) -> Result<SourceElasticPathAverage, NormalPatchError> {
        for indentation_m in [initial_indentation_m, final_indentation_m] {
            if !indentation_m.is_finite() || indentation_m < 0.0 {
                return Err(NormalPatchError::InvalidInput {
                    field: "source point indentation",
                });
            }
            let limit = self.envelope().maximum_indentation_m;
            if indentation_m > limit {
                return Err(NormalPatchError::OutsideApplicability {
                    ratio: "source point indentation",
                    value: indentation_m,
                    limit,
                });
            }
        }
        self.path_average(initial_indentation_m, final_indentation_m)
    }

    fn path_average(
        &self,
        initial_indentation_m: f64,
        final_indentation_m: f64,
    ) -> Result<SourceElasticPathAverage, NormalPatchError> {
        let increment = final_indentation_m - initial_indentation_m;
        let depth = initial_indentation_m.max(final_indentation_m);
        if depth == 0.0 {
            return Ok(SourceElasticPathAverage {
                mean_force_n: 0.0,
                endpoint_derivative_n_per_m: 0.0,
                work_j: 0.0,
            });
        }
        let scale = depth.sqrt();
        let a = initial_indentation_m.sqrt() / scale;
        let b = final_indentation_m.sqrt() / scale;
        let a2 = a * a;
        let b2 = b * b;
        let sum = a + b;
        let mean_shape = 0.4 * (b2 * b2 + b2 * b * a + b2 * a2 + b * a2 * a + a2 * a2) / sum;
        let derivative_shape =
            0.2 * (3.0 * b2 * b + 6.0 * a * b2 + 4.0 * a2 * b + 2.0 * a2 * a) / (sum * sum);
        let stiffness_scale = self.stiffness_n_per_m_3_2() * scale;
        let mean_force_n = (stiffness_scale * mean_shape) * depth;
        let endpoint_derivative_n_per_m = stiffness_scale * derivative_shape;
        let work_j = mean_force_n * increment;
        if !mean_force_n.is_finite()
            || !endpoint_derivative_n_per_m.is_finite()
            || !work_j.is_finite()
        {
            return Err(NormalPatchError::Overflow {
                field: "elastic source path average",
            });
        }
        Ok(SourceElasticPathAverage {
            mean_force_n,
            endpoint_derivative_n_per_m,
            work_j,
        })
    }
}
