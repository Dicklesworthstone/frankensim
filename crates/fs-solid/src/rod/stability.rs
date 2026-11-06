//! Conservative stability of the native clamped rod at a dead-force equilibrium.
//!
//! A small force residual locates a stationary configuration; it does not make
//! that configuration a minimum. The original analytical tangent includes the
//! geometric stiffness of compression. Its symmetric second variation must be
//! positive on all free coordinates before this screen can report stability.
//! For equilibrium symmetry and the external-moment restriction, see section
//! 3.2.1.1 of <https://mediatum.ub.tum.de/doc/1306287/1306287.pdf>.

use super::{Rod, RodEnergy, SolidError, TipLoad, UniformEnergy};
use fs_la::factor::cholesky;

impl Rod {
    /// Screen a converged uniform rod for local stability under a dead tip force.
    ///
    /// Node zero remains clamped in position and frame, as in [`Self::solve_static`].
    /// The actual free force/couple residual must satisfy `equilibrium_tolerance`
    /// in that solver's norm convention. A nonzero tip moment refuses: the
    /// follower-couple and world-couple laws do not supply the conservative
    /// potential required by this test.
    ///
    /// The original analytical body/right Jacobian is symmetric at equilibrium
    /// with zero applied moments. Its symmetric part is the second variation
    /// in the local exponential coordinates; symmetrization here does not alter
    /// the Newton owner. A positive diagonal congruence removes the different
    /// translation/rotation scales, then Cholesky tests strict positivity after
    /// a `64 * ndof * EPSILON` dimensionless margin. `false` includes negative,
    /// neutral and numerically unresolved tangents. This is a floating-point
    /// screen of this discretized equilibrium, not a continuum buckling load,
    /// dynamic stability proof, path continuation, or material-validity claim.
    ///
    /// Work is bounded to 1..=64 segments (384 free DOFs). Like the legacy raw
    /// rod statics entry, this bounded query has no cancellation context.
    /// It never mutates the rod.
    ///
    /// # Errors
    /// Refuses inadmissible geometry/section parameters, nonfinite loads,
    /// nonzero tip moments, an unconverged state, or an unrepresentable tangent.
    pub fn is_stable_under_dead_load(
        &self,
        load: &TipLoad,
        equilibrium_tolerance: f64,
    ) -> Result<bool, SolidError> {
        let invalid = |what: &str| SolidError::InvalidInput { what: what.into() };
        if !(2..=65).contains(&self.positions.len())
            || self.frames.len() != self.positions.len()
            || !self.l0.is_finite()
            || self.l0 <= 0.0
            || self
                .positions
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
            || [
                self.section.ea,
                self.section.ga,
                self.section.gj,
                self.section.ei,
            ]
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
            || !equilibrium_tolerance.is_finite()
            || equilibrium_tolerance <= 0.0
            || load.force.iter().any(|value| !value.is_finite())
            || load.moment.iter().any(|value| *value != 0.0)
        {
            return Err(invalid(
                "rod stability requires 1..=64 finite positive-stiffness segments, finite dead force, zero tip moment and positive tolerance",
            ));
        }
        let dimension = self.ndof();
        let free: Vec<_> = (0..dimension).collect();
        let exact = UniformEnergy
            .linearize(self, &free, true, &|| Ok(()))?
            .ok_or_else(|| SolidError::InternalInvariant {
                what: "native rod omitted its analytical stability tangent".into(),
            })?;
        let mut residual = 0.0;
        for (node, gradient) in exact.gradient.iter().enumerate().skip(1) {
            for (component, &value) in gradient.iter().enumerate() {
                let applied = if node + 1 == self.positions.len() && component < 3 {
                    load.force[component]
                } else {
                    0.0
                };
                residual = fs_math::det::hypot(residual, value - applied);
            }
        }
        if !residual.is_finite() || residual > equilibrium_tolerance {
            return Err(invalid(
                "rod stability requires a converged dead-force equilibrium",
            ));
        }
        let mut tangent = exact.tangent.ok_or_else(|| SolidError::InternalInvariant {
            what: "native rod omitted its requested stability tangent".into(),
        })?;
        let mut scales = Vec::with_capacity(dimension);
        for row in 0..dimension {
            let diagonal = tangent[row * dimension + row];
            if diagonal <= 0.0 {
                return Ok(false);
            }
            scales.push(fs_math::det::sqrt(diagonal));
        }
        // For a scalar potential on R^3 x SO(3), the body/right connection
        // correction is skew; only the symmetric Jacobian contributes to the
        // directional second variation. Dead forces have no load Hessian.
        for row in 0..dimension {
            for column in 0..=row {
                let a = tangent[row * dimension + column] / scales[row] / scales[column];
                let b = tangent[column * dimension + row] / scales[column] / scales[row];
                let value = 0.5 * a + 0.5 * b;
                if !value.is_finite() {
                    return Err(invalid("rod stability tangent cannot be equilibrated"));
                }
                tangent[row * dimension + column] = value;
                tangent[column * dimension + row] = value;
            }
        }
        // The dimension is bounded above by 384, exactly representable in f64.
        #[allow(clippy::cast_precision_loss)]
        let margin = 64.0 * dimension as f64 * f64::EPSILON;
        for row in 0..dimension {
            tangent[row * dimension + row] -= margin;
        }
        Ok(cholesky(&tangent, dimension).is_ok())
    }
}
