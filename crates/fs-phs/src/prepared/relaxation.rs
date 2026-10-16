//! Exact scalar condensation of the existing Maxwell midpoint equation.
//! This eliminates only a linear internal coordinate, not a time integrator.
use crate::{PhsError, RelaxationBranch};
use fs_math::det;

impl RelaxationBranch {
    fn scalar_coupling(&self) -> Result<f64, PhsError> {
        if self.projection.len() != 1 || !self.stiffness.is_finite() || self.stiffness <= 0.0
            || !self.relaxation_time_s.is_finite() || self.relaxation_time_s <= 0.0
            || !(1.0 / self.relaxation_time_s).is_finite()
        { return Err(PhsError::RelaxationParameters); }
        let c = det::sqrt(self.stiffness) * self.projection[0];
        if !c.is_finite() || c == 0.0 { return Err(PhsError::RelaxationParameters); }
        Ok(c)
    }

    /// Energy-normalized memory z for a fully relaxed scalar coordinate x.
    /// The arm's energy is `(sqrt(k) * projection[0] * x - z)^2 / 2`.
    /// Equilibrium stiffness belongs to the original system, not this arm.
    ///
    /// # Errors
    /// Non-scalar/nonfinite branch or an unrepresentable physical state.
    pub fn scalar_relaxed_memory(&self, x: f64) -> Result<f64, PhsError> {
        let z = self.scalar_coupling()? * x;
        if !x.is_finite() || !z.is_finite() { return Err(PhsError::RelaxationParameters); }
        Ok(z)
    }

    /// Stored arm energy at scalar coordinate x and energy-normalized memory z.
    /// This is the same quadratic storage as `with_relaxation_branches`.
    ///
    /// # Errors
    /// Invalid branch/state or nonfinite energy.
    pub fn scalar_stored_energy(&self, x: f64, z: f64) -> Result<f64, PhsError> {
        let extension = self.scalar_relaxed_memory(x)? - z;
        let energy = 0.5 * extension * extension;
        if !z.is_finite() || !energy.is_finite() { return Err(PhsError::RelaxationParameters); }
        Ok(energy)
    }

    /// Eliminate z from this owner's implicit midpoint Maxwell equation.
    /// Returns `(z_next, restoring_force, stored_energy_next, dissipated_energy)`.
    /// The force is conjugate to x and opposes increasing elastic extension;
    /// a mechanical balance therefore subtracts this force. It is evaluated at
    /// the SAME midpoint as the hosting mechanics, not lagged or split afterward.
    /// `force * (x_next-x_old) = H_next-H_old + dissipated` to rounding.
    ///
    /// No state is mutated. A host may reject a nonlinear trial or another
    /// participant before publishing z_next. This is midpoint, NOT exponential
    /// relaxation: dt >> tau is passive but can alternate memory signs. Accuracy
    /// and material-use bandwidth must be admitted by the physical consumer.
    ///
    /// # Errors
    /// Non-scalar branch, invalid state, nonpositive dt, or arithmetic overflow.
    pub fn scalar_midpoint_step(&self, x_old: f64, x_next: f64, z_old: f64, dt: f64)
        -> Result<(f64, f64, f64, f64), PhsError>
    {
        if ![x_old, x_next, z_old, dt].iter().all(|x| x.is_finite()) || dt <= 0.0 {
            return Err(PhsError::RelaxationParameters);
        }
        let c = self.scalar_coupling()?;
        let tau = self.relaxation_time_s;
        // a = tau/(tau+dt/2), b = (dt/2)/(tau+dt/2). Scale both first so
        // neither very stiff nor very slow arms overflow the denominator.
        let half_dt = 0.5 * dt;
        let scale = tau.max(half_dt);
        let (t, h) = (tau / scale, half_dt / scale);
        let (a, b) = (t / (t+h), h / (t+h));
        let mismatch = c * f64::midpoint(x_old, x_next) - z_old;
        let midpoint_effort = a * mismatch;
        let z_next = z_old + 2.0 * b * mismatch;
        let force = c * midpoint_effort;
        let energy = self.scalar_stored_energy(x_next, z_next)?;
        let dissipated = 2.0 * midpoint_effort * (b * mismatch);
        if ![z_next, force, energy, dissipated].iter().all(|x| x.is_finite())
            || dissipated < 0.0
        { return Err(PhsError::RelaxationParameters); }
        Ok((z_next, force, energy, dissipated))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{PortHamiltonian, QuadraticStorage, step};
    fn branch(k: f64, tau: f64, b: f64) -> RelaxationBranch {
        RelaxationBranch { projection: vec![b], stiffness:k, relaxation_time_s:tau }
    }
    fn near(a: f64, b: f64, scale: f64) {
        assert!((a-b).abs() <= 3e-10*scale.max(1e-12), "{a:e} != {b:e}");
    }
    #[test]
    fn condensed_coordinate_matches_the_full_owner_state_force_and_work() {
        for orientation in [-0.7, 1.0, 2.3] {
            let arm=branch(7.0,0.03,orientation);
            // x'=u is an ideal prescribed-strain port; its output is the exact
            // work-conjugate material force. The original owner appends z.
            let system=PortHamiltonian::new(1,1,vec![0.0],vec![0.0],vec![1.0],
                Box::new(QuadraticStorage::new(vec![0.0],1).unwrap())).unwrap()
                .with_relaxation_branches(vec![arm.clone()]).unwrap();
            let (mut x,mut z)=(0.2,-0.1);
            for n in 0..40 {
                let dt=0.002;
                let next=if n<20 { x+0.001 } else { x-0.002 };
                let full=step(&system,&[x,z],&[(next-x)/dt],dt).unwrap();
                let (zn,force,energy,loss)=arm.scalar_midpoint_step(x,next,z,dt).unwrap();
                near(zn,full.x[1],zn.abs()+full.x[1].abs());
                near(force,full.y[0],force.abs()+full.y[0].abs());
                near(energy,system.hamiltonian(&full.x),energy.abs());
                near(loss,full.dissipated,energy.abs());
                near(force*(next-x),full.supplied,energy.abs());
                x=next; z=zn;
            }
        }
    }
    #[test]
    fn holds_release_and_reversals_obey_the_unmodified_quadratic_storage() {
        let arm=branch(3.0,0.2,1.0);
        for ratio in [1e-9,1e-3,1.0,1e3,1e9] {
            let (x0,x1,z0)=(0.7,-0.4,0.2);
            let old=arm.scalar_stored_energy(x0,z0).unwrap();
            let (_,force,new,loss)=arm.scalar_midpoint_step(x0,x1,z0,ratio*0.2).unwrap();
            near(new-old+loss,force*(x1-x0),old+new+loss+force.abs());
            assert!(loss>=0.0);
        }
        let (mut z,mut previous)=(0.0,arm.scalar_stored_energy(0.3,0.0).unwrap());
        for _ in 0..100 {
            let (next,_,energy,loss)=arm.scalar_midpoint_step(0.3,0.3,z,0.02).unwrap();
            assert!(energy<previous); near(previous-energy,loss,previous);
            z=next;previous=energy;
        }
        let relaxed=arm.scalar_relaxed_memory(0.3).unwrap();
        assert_eq!(arm.scalar_midpoint_step(0.3,0.3,relaxed,0.02).unwrap(),(relaxed,0.0,0.0,0.0));
    }
    #[test]
    fn signed_coordinate_scaling_changes_force_but_not_energy_or_memory() {
        let a=branch(10.0,0.003,0.7); let b=branch(10.0,0.003,-0.35);
        let x=a.scalar_midpoint_step(0.1,-0.2,0.17,0.001).unwrap();
        let y=b.scalar_midpoint_step(-0.2,0.4,0.17,0.001).unwrap();
        assert_eq!(x.0,y.0); assert_eq!(x.2,y.2); assert_eq!(x.3,y.3);
        assert_eq!(x.1,-2.0*y.1);
    }
    #[test]
    fn malformed_scalar_branches_and_nonfinite_trials_are_refused() {
        let good=branch(1.0,1.0,1.0);
        for bad in [f64::NAN,f64::INFINITY] {
            assert!(good.scalar_midpoint_step(bad,0.0,0.0,1.0).is_err());
            assert!(good.scalar_midpoint_step(0.0,bad,0.0,1.0).is_err());
            assert!(good.scalar_midpoint_step(0.0,0.0,bad,1.0).is_err());
        }
        for dt in [0.0,-1.0,f64::INFINITY] { assert!(good.scalar_midpoint_step(0.0,1.0,0.0,dt).is_err()); }
        for arm in [branch(0.0,1.0,1.0),branch(1.0,0.0,1.0),branch(1.0,1.0,0.0),
            branch(1.0,1.0,f64::NAN),RelaxationBranch {projection:vec![1.0,0.0],..good.clone()}] {
            assert!(arm.scalar_midpoint_step(0.0,1.0,0.0,0.1).is_err());
        }
    }
}
