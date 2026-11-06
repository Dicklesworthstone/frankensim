//! Inviscid terms of the material vorticity equation.
//!
//! Curling rho Du/Dt = -grad(p) gives
//! Domega/Dt = (omega.grad)u - omega div(u) + grad(rho) x grad(p)/rho^2.
//! Viscous stress and nonconservative body-force curls are separate terms.
//! See MIT 12.800 (2004), lecture 13, and NISTIR 6783, section 2.3.

use super::{FlowKinematics3, KinematicsError};

/// Resolved local inviscid contributions, in one coherent unit system.
/// Native lattice queries use lattice time; `to_physical` converts once to SI.
///
/// These terms describe a material derivative. They do not include the
/// Eulerian transport term -u.grad(omega), viscous diffusion, variable-viscosity
/// stress curl, body-force curl, or a turbulence model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VorticityBudget3 {
    /// Full velocity gradient, including dilation.
    pub flow: FlowKinematics3,
    /// (omega.grad)u, including tilting.
    pub stretching: [f64; 3],
    /// -omega div(u); compressible expansion attenuates vorticity.
    pub dilatation: [f64; 3],
    /// grad(rho) x grad(p) / rho^2, with thermodynamic absolute pressure.
    pub baroclinic: [f64; 3],
}

impl FlowKinematics3 {
    /// Evaluate the inviscid terms using density and scalar gradients from
    /// the SAME state and coordinate/unit system as this velocity gradient.
    /// Pressure is thermodynamic pressure, not a kinematic pressure p/rho.
    ///
    /// # Errors
    /// Rejects nonpositive/nonfinite density, nonfinite gradients, or any
    /// unrepresentable normalized gradient, rate or enstrophy production.
    pub fn inviscid_vorticity_budget(
        self,
        density: f64,
        density_gradient: [f64; 3],
        pressure_gradient: [f64; 3],
    ) -> Result<VorticityBudget3, KinematicsError> {
        if !(density.is_finite() && density > 0.0)
            || density_gradient.iter().chain(&pressure_gradient).any(|v| !v.is_finite())
        {
            return Err(KinematicsError::Input(
                "vorticity budget needs positive density and finite scalar gradients",
            ));
        }
        // Normalize BEFORE the cross product: rho^2 need not be representable
        // for the physical baroclinic rate to be finite.
        let a = density_gradient.map(|v| v / density);
        let b = pressure_gradient.map(|v| v / density);
        if a.iter().chain(&b).any(|v| !v.is_finite()) {
            return Err(KinematicsError::NonFiniteGradient);
        }
        let budget = VorticityBudget3 {
            flow: self,
            stretching: self.vortex_stretching(),
            dilatation: self.vorticity().map(|w| -w * self.divergence()),
            baroclinic: [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ],
        };
        budget.validate()?;
        Ok(budget)
    }
}

impl VorticityBudget3 {
    /// Sum of the three represented inviscid material-rate terms.
    /// Units are s^-2 after physical conversion.
    #[must_use]
    pub fn inviscid_rate(self) -> [f64; 3] {
        core::array::from_fn(|a| {
            self.stretching[a] + self.dilatation[a] + self.baroclinic[a]
        })
    }

    /// Inviscid local rate of |omega|^2/2; s^-3 after physical conversion.
    /// This is not a volume-integrated or viscous enstrophy balance.
    #[must_use]
    pub fn inviscid_enstrophy_rate(self) -> f64 {
        let w = self.flow.vorticity();
        let rate = self.inviscid_rate();
        w[0] * rate[0] + w[1] * rate[1] + w[2] * rate[2]
    }

    /// Convert a native lattice budget once using that state's timestep.
    /// Curl scales as 1/dt; every represented source term scales as 1/dt^2.
    /// The lattice spacing and reference density cancel, including in the
    /// baroclinic term. No density or pressure reference is substituted.
    pub fn to_physical(self, dt_s: f64) -> Result<Self, KinematicsError> {
        let budget = Self {
            flow: self.flow.to_physical(dt_s)?,
            stretching: self.stretching.map(|v| (v / dt_s) / dt_s),
            dilatation: self.dilatation.map(|v| (v / dt_s) / dt_s),
            baroclinic: self.baroclinic.map(|v| (v / dt_s) / dt_s),
        };
        budget.validate()?;
        Ok(budget)
    }

    fn validate(self) -> Result<(), KinematicsError> {
        if self.stretching.iter()
            .chain(&self.dilatation)
            .chain(&self.baroclinic)
            .chain(&self.inviscid_rate())
            .any(|v| !v.is_finite())
            || !self.inviscid_enstrophy_rate().is_finite()
        {
            return Err(KinematicsError::NonFiniteGradient);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: f64, b: f64) {
        assert!((a - b).abs() <= 2e-14 * (1.0 + b.abs()), "{a} != {b}");
    }

    #[test]
    fn expansion_rotation_and_baroclinicity_have_distinct_signs() {
        // u = (a*x-b*y, b*x+a*y, a*z): omega_z=2b,
        // stretching_z=2ab, dilation_z=-6ab, net_z=-4ab.
        let flow = FlowKinematics3::from_gradient([
            [0.2, -0.3, 0.0],
            [0.3, 0.2, 0.0],
            [0.0, 0.0, 0.2],
        ]).unwrap();
        let q = flow.inviscid_vorticity_budget(2.0, [0.5, 0.0, 0.0], [0.0, 4.0, 0.0]).unwrap();
        near(q.stretching[2], 0.12);
        near(q.dilatation[2], -0.36);
        near(q.baroclinic[2], 0.5);
        near(q.inviscid_rate()[2], 0.26);
        near(q.inviscid_enstrophy_rate(), 0.156);
        let physical = q.to_physical(0.25).unwrap();
        for axis in 0..3 {
            near(physical.flow.vorticity()[axis], 4.0 * flow.vorticity()[axis]);
            near(physical.inviscid_rate()[axis], 16.0 * q.inviscid_rate()[axis]);
        }
        near(physical.inviscid_enstrophy_rate(), 64.0 * q.inviscid_enstrophy_rate());
    }

    #[test]
    fn barotropic_and_density_unit_changes_preserve_the_physics() {
        let rest = FlowKinematics3::from_gradient([[0.0; 3]; 3]).unwrap();
        let q = rest.inviscid_vorticity_budget(2.0, [1.0, -2.0, 3.0], [4.0, -8.0, 12.0]).unwrap();
        assert_eq!(q.inviscid_rate(), [0.0; 3]);
        for scale in [1e-150, 1.0, 1e150] {
            let q = rest.inviscid_vorticity_budget(
                2.0 * scale, [scale, 0.0, 0.0], [0.0, 3.0 * scale, 0.0],
            ).unwrap();
            near(q.baroclinic[2], 0.75);
        }
    }

    #[test]
    fn invalid_density_gradients_and_scaling_refuse() {
        let rest = FlowKinematics3::from_gradient([[0.0; 3]; 3]).unwrap();
        for rho in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(rest.inviscid_vorticity_budget(rho, [0.0; 3], [0.0; 3]).is_err());
        }
        assert!(rest.inviscid_vorticity_budget(1.0, [f64::NAN; 3], [0.0; 3]).is_err());
        assert!(rest.inviscid_vorticity_budget(1.0, [0.0; 3], [f64::INFINITY; 3]).is_err());
        let q = rest.inviscid_vorticity_budget(1.0, [0.0; 3], [0.0; 3]).unwrap();
        for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(q.to_physical(dt).is_err());
        }
    }
}
