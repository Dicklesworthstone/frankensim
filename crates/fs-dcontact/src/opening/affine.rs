//! Pointwise contact work for a distributed obstacle on one generalized coordinate.
//! Each point keeps its original potential and local normal velocity. The
//! non-attractive unloading clamp is applied BEFORE modal force accumulation.
use super::{DContactError, Obstacle, OpeningContactStep};

/// Reaction in the caller's scalar coordinate (not necessarily a local gap).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AffineContactResponse {
    /// Conservative generalized force; its displacement work is minus the
    /// change of the original summed contact potential, up to arithmetic error.
    pub elastic_force: f64,
    /// Generalized force including pointwise Hunt--Crossley loss and unloading.
    pub force: f64,
    /// Nonnegative dissipative power in the same coordinate and time level.
    pub dissipated_power: f64,
}

/// Immutable trial adapter for p_i = b_i q - c_i. There is ONE mechanical
/// coordinate, but any admitted number of independent contact quadrature points.
/// Coefficients, weights, gaps and provenance remain in the original obstacle.
/// No trial allocation, stored-state mutation or new contact law is introduced.
pub struct AffineContactStep<'a> {
    law: &'a Obstacle,
    before: f64,
}
impl<'a> AffineContactStep<'a> {
    /// Validate the complete scalar collocation and caller's point-work budget.
    /// Zero and mixed-sign rows are legal: physical orientation is not discarded.
    ///
    /// # Errors
    /// Malformed raw parts, nonfinite/invalid coefficients, empty or over-budget
    /// quadrature, or a nonfinite initial scalar coordinate.
    pub fn new(law: &'a Obstacle, before: f64, max_points: usize) -> Result<Self, DContactError> {
        let n = law.n_points();
        if n == 0 || n > max_points || law.collocation().len() != n
            || law.gaps().len() != n || law.weights().len() != n {
            return Err(DContactError::Shape { what: "affine contact needs bounded nonempty scalar collocation" });
        }
        if !before.is_finite() || law.collocation().iter().chain(law.gaps()).chain(law.weights())
            .any(|x| !x.is_finite()) || law.weights().iter().any(|x| *x < 0.0)
            || !law.stiffness().is_finite() || law.stiffness() < 0.0
            || !law.alpha().is_finite() || law.alpha() < 1.0
            || !law.internal_loss().is_finite() || law.internal_loss() < 0.0 {
            return Err(bad());
        }
        Ok(Self { law, before })
    }

    /// Evaluate each existing stable opening-potential secant in its LOCAL
    /// opening -b_i q, then pull the local force back with -b_i. At generalized
    /// velocity v, local opening velocity is -b_i v. Dissipation is summed
    /// pointwise; a clamp on the total force would be wrong for unequal rows.
    /// The caller must use v = (after-before)/dt for a discrete energy balance.
    ///
    /// # Errors
    /// Nonfinite coordinates, velocity, secants or accumulated reactions.
    pub fn response(&self, after: f64, velocity: f64) -> Result<AffineContactResponse, DContactError> {
        if !after.is_finite() || !velocity.is_finite() { return Err(bad()); }
        let mut result = AffineContactResponse { elastic_force: 0.0, force: 0.0, dissipated_power: 0.0 };
        for (i, &b) in self.law.collocation().iter().enumerate() {
            let local_before = -b * self.before;
            let local_after = -b * after;
            let local_velocity = -b * velocity;
            if ![local_before, local_after, local_velocity].iter().all(|x| x.is_finite()) { return Err(bad()); }
            // Fields belong to this same owner. No temporary Obstacle, duplicate
            // power law or subtractive energy secant is constructed here.
            let point = OpeningContactStep {
                before: local_before, gap: self.law.gaps()[i],
                stiffness: self.law.stiffness(), weight: self.law.weights()[i],
                alpha: self.law.alpha(), loss: self.law.internal_loss(),
            };
            let (elastic, damping) = point.coefficients(local_after)?;
            let viscous_force = damping * local_velocity;
            if !viscous_force.is_finite() { return Err(bad()); }
            // Retain the loss increment before adding the elastic reaction:
            // subtracting two nearly equal total forces would erase weak
            // physical dissipation. Cap unloading locally at zero attraction.
            let loss_force = viscous_force.min(elastic);
            let normal_force = elastic - loss_force;
            if !normal_force.is_finite() { return Err(bad()); }
            result.elastic_force -= b * elastic;
            result.force -= b * normal_force;
            result.dissipated_power += loss_force * local_velocity;
        }
        if ![result.elastic_force, result.force, result.dissipated_power].iter().all(|x| x.is_finite())
            || result.dissipated_power < 0.0 { return Err(bad()); }
        Ok(result)
    }
}
fn bad() -> DContactError {
    DContactError::Parameter { what: "affine contact requires finite physical coefficients, states and reactions" }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn law(rows: Vec<f64>, gaps: Vec<f64>, weights: Vec<f64>, chi: f64) -> Obstacle {
        Obstacle::new(rows, gaps.len(), 1, gaps, weights, 3e6, 2.0, "affine contact test".into())
            .unwrap().with_internal_loss(chi).unwrap()
    }
    fn energy(law: &Obstacle, q: f64) -> f64 {
        law.collocation().iter().zip(law.gaps()).zip(law.weights()).map(|((&b, &c), &w)|
            w * law.stiffness() * (b*q-c).max(0.0).powi(3)/3.0).sum()
    }
    #[test]
    fn mixed_orientations_and_contact_crossings_balance_original_potential_work() {
        for chi in [0.0, 0.2, 5.0] {
            let ob = law(vec![-2.0, -0.3, 0.7, 0.0], vec![0.0001, -0.0001, 0.0002, -0.0003], vec![0.2,0.4,0.3,0.1], chi);
            for a in [-0.003, -0.0001, 0.0, 0.002] {
                for b in [-0.003, -0.0001, 0.0, 0.002] {
                    let dt = 0.002;
                    let response = AffineContactStep::new(&ob, a, 4).unwrap().response(b, (b-a)/dt).unwrap();
                    let change = energy(&ob,b)-energy(&ob,a);
                    let residual = change + response.force*(b-a) + response.dissipated_power*dt;
                    let scale = energy(&ob,a)+energy(&ob,b)+(response.force*(b-a)).abs()+response.dissipated_power*dt;
                    assert!(residual.abs() <= 1e-12*scale.max(1e-30), "{a} {b}: {residual:e}");
                    assert!(response.dissipated_power >= 0.0);
                }
            }
        }
    }
    #[test]
    fn local_unloading_clamps_precede_generalized_force_accumulation() {
        let ob = law(vec![-2.0,-0.5],vec![0.0,0.0],vec![1.0,1.0],1.0);
        let response = AffineContactStep::new(&ob,-0.01,2).unwrap().response(-0.01,0.75).unwrap();
        // The first point unloads beyond its cap, the second still reacts.
        let expected = 0.5 * 3e6 * 0.005_f64.powi(2) * (1.0-0.5*0.75);
        assert!((response.force-expected).abs() < 1e-12*expected);
        assert!(response.dissipated_power > 0.0);
    }
    #[test]
    fn quadrature_splitting_keeps_the_same_physical_law() {
        let whole = law(vec![-0.7],vec![0.0002],vec![0.8],0.5);
        let split = law(vec![-0.7;4],vec![0.0002;4],vec![0.2;4],0.5);
        let a = AffineContactStep::new(&whole,-0.002,4).unwrap().response(-0.001,0.1).unwrap();
        let b = AffineContactStep::new(&split,-0.002,4).unwrap().response(-0.001,0.1).unwrap();
        assert_eq!(a,b);
    }
    #[test]
    fn weak_contact_loss_survives_rounding_into_the_elastic_reaction() {
        let ob = law(vec![-1.0], vec![0.0], vec![1.0], 1e-20);
        for velocity in [-0.1, 0.1] {
            let response = AffineContactStep::new(&ob, -0.001, 1).unwrap()
                .response(-0.001, velocity).unwrap();
            // G1: the instantaneous HC loss is chi*K*p^alpha*v^2.
            // It is well represented even though adding the viscous increment
            // to the much larger conservative force rounds it away.
            let expected = 1e-20 * 3e6 * 0.001_f64.powi(2) * velocity * velocity;
            assert_eq!(response.force, response.elastic_force);
            assert!(response.dissipated_power > 0.0);
            assert!((response.dissipated_power - expected).abs() < 1e-13 * expected);
        }
    }
    #[test]
    fn capped_unloading_and_loading_accumulate_before_mixed_sign_projection() {
        let ob = law(vec![-2.0, 0.5], vec![-0.004, -0.004], vec![0.2, 0.8], 1.0);
        let response = AffineContactStep::new(&ob, 0.0, 2).unwrap()
            .response(0.0, 0.75).unwrap();
        // First point unloads at 1.5 m/s and caps to zero reaction;
        // second point loads at -0.375 m/s. Each must dissipate positively
        // even though their generalized conservative forces cancel here.
        let first_elastic = 0.2 * 3e6 * 0.004_f64.powi(2);
        let second_elastic = 0.8 * 3e6 * 0.004_f64.powi(2);
        let expected = first_elastic * 1.5 + second_elastic * 0.375_f64.powi(2);
        assert!((response.dissipated_power - expected).abs() < 1e-13 * expected);
        assert_eq!(response.elastic_force, 0.0);
    }
    #[test]
    fn raw_shape_budget_and_nonfinite_trials_refuse_before_indexing() {
        let ob = law(vec![-1.0,-2.0],vec![0.0;2],vec![1.0;2],0.2);
        assert!(AffineContactStep::new(&ob,0.0,1).is_err());
        let raw = Obstacle::from_raw_parts(vec![-1.0],usize::MAX,vec![],vec![],1.0,2.0,"raw".into());
        assert!(AffineContactStep::new(&raw,0.0,usize::MAX).is_err());
        let step = AffineContactStep::new(&ob,0.0,2).unwrap();
        for bad in [f64::NAN,f64::INFINITY] {
            assert!(step.response(bad,0.0).is_err()); assert!(step.response(0.0,bad).is_err());
        }
        assert_eq!(step.response(0.0,0.0).unwrap().force,0.0);
    }
}
