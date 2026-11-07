//! Native accepted-state compressible vorticity sources.
//!
//! This observes the thermal momentum owner; it does not add another force.
//! Thermal-pressure forcing already generates baroclinic vorticity through
//! the conservative momentum equation.

use super::{Lattice3, NsfGrid3, half_force_velocity, raw_moments};
use crate::kinematics::{FlowKinematics3, KinematicsError, VorticityBudget3};

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Inviscid material-vorticity terms from the accepted thermal flow, in
    /// lattice units. The physical velocity includes the current half-force
    /// correction; pressure is the accepted thermodynamic p = rho*theta.
    /// The original gas source already supplied theta during state acceptance.
    ///
    /// A centered stencil uses both fluid neighbours along every axis,
    /// wrapping periodic faces exactly as the native owner does. A solid or
    /// nonperiodic boundary interrupting this stencil refuses: neither a zero
    /// ghost nor an invented thermodynamic wall condition is substituted.
    /// The existing `kinematics` field remains available at those boundaries.
    /// One-cell periodic axes are homogeneous; two-cell axes refuse because
    /// their centered neighbors alias and cannot resolve a derivative.
    ///
    /// Seven fixed-size population reductions, no allocation or state change.
    /// For SI use `sample.to_physical(self.law().dt_s)`; all source terms
    /// scale as 1/dt^2. Viscous stress curl, imposed body-force curl and source
    /// mass/momentum injection are outside this inviscid budget.
    ///
    /// # Errors
    /// Out-of-grid/solid cells, an incomplete fluid stencil, invalid accepted
    /// data or unrepresentable derivatives/rates.
    pub fn vorticity_budget(
        &self,
        cell: [usize; 3],
    ) -> Result<VorticityBudget3, KinematicsError> {
        if cell.iter().zip(self.spec.dims).any(|(&i, n)| i >= n) {
            return Err(KinematicsError::Input("vorticity cell outside the grid"));
        }
        let center = self.index(cell);
        if self.solid[center] {
            return Err(KinematicsError::Input("solid cell has no fluid vorticity"));
        }
        let (density, _, _) = self.vorticity_sample(center)?;
        let mut gradient = [[0.0; 3]; 3];
        let mut density_gradient = [0.0; 3];
        let mut pressure_gradient = [0.0; 3];
        for axis in 0..3 {
            if self.spec.dims[axis] == 2 {
                return Err(KinematicsError::Input("two-cell vorticity derivative is unresolved"));
            }
            let low = self.fluid_at(cell, axis, -1).ok_or(KinematicsError::Input(
                "vorticity sources require two fluid neighbours per axis",
            ))?;
            let high = self.fluid_at(cell, axis, 1).ok_or(KinematicsError::Input(
                "vorticity sources require two fluid neighbours per axis",
            ))?;
            let (rho_low, velocity_low, pressure_low) = self.vorticity_sample(low)?;
            let (rho_high, velocity_high, pressure_high) = self.vorticity_sample(high)?;
            density_gradient[axis] = 0.5 * (rho_high - rho_low);
            pressure_gradient[axis] = 0.5 * (pressure_high - pressure_low);
            for component in 0..3 {
                gradient[component][axis] =
                    0.5 * (velocity_high[component] - velocity_low[component]);
            }
        }
        FlowKinematics3::from_gradient(gradient)?.inviscid_vorticity_budget(
            density,
            density_gradient,
            pressure_gradient,
        )
    }

    /// Vorticity production curl(a) of the held spatial body acceleration,
    /// in physical s^-2. Uses the original accepted SI field and lattice dx.
    ///
    /// This is the acceleration curl, not curl(rho*a)/rho: a constant
    /// acceleration creates no vorticity merely because density varies.
    /// Uniform gravity has zero curl and is already separate. Manufactured
    /// mass/momentum sources and viscous stress are not part of this value.
    ///
    /// Uses the same centered fluid-neighbor topology as the inviscid query.
    /// A one-cell periodic axis is a homogeneous extrusion. A two-cell axis,
    /// solid neighbor or missing nonperiodic neighbor refuses rather than
    /// reporting an unresolved derivative as zero.
    pub fn body_acceleration_curl_si(
        &self,
        cell: [usize; 3],
    ) -> Result<[f64; 3], KinematicsError> {
        if cell.iter().zip(self.spec.dims).any(|(&i, n)| i >= n) {
            return Err(KinematicsError::Input("vorticity cell outside the grid"));
        }
        if self.solid[self.index(cell)] {
            return Err(KinematicsError::Input("solid cell has no fluid vorticity"));
        }
        let mut gradient = [[0.0; 3]; 3];
        for axis in 0..3 {
            if self.spec.dims[axis] == 2 {
                return Err(KinematicsError::Input("two-cell body-acceleration derivative is unresolved"));
            }
            let low = self.fluid_at(cell, axis, -1).ok_or(KinematicsError::Input(
                "body-acceleration curl requires two fluid neighbours per axis",
            ))?;
            let high = self.fluid_at(cell, axis, 1).ok_or(KinematicsError::Input(
                "body-acceleration curl requires two fluid neighbours per axis",
            ))?;
            let a_low = self.body_acceleration.get(low)
                .map_or([0.0; 3], |a| a.acceleration_m_s2);
            let a_high = self.body_acceleration.get(high)
                .map_or([0.0; 3], |a| a.acceleration_m_s2);
            for component in 0..3 {
                gradient[component][axis] =
                    (0.5 * (a_high[component] - a_low[component])) / self.law.dx_m;
            }
        }
        let curl = [
            gradient[2][1] - gradient[1][2],
            gradient[0][2] - gradient[2][0],
            gradient[1][0] - gradient[0][1],
        ];
        if gradient.iter().flatten().chain(&curl).any(|v| !v.is_finite()) {
            return Err(KinematicsError::NonFiniteGradient);
        }
        Ok(curl)
    }

    fn vorticity_sample(&self, i: usize) -> Result<(f64, [f64; 3], f64), KinematicsError> {
        let (rho, momentum) = raw_moments::<L>(&self.f[i]);
        let velocity = half_force_velocity(rho, momentum, self.force[i]);
        let pressure = rho * self.theta[i];
        if !(rho.is_finite() && rho > 0.0 && pressure.is_finite() && pressure > 0.0)
            || velocity.iter().any(|v| !v.is_finite())
        {
            return Err(KinematicsError::InvalidState { cell: self.at(i) });
        }
        Ok((rho, velocity, pressure))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q27::D3q27;
    use crate::kinematics::KinematicsWorkspace3;
    use crate::nsf3::{NsfSpec, THETA0, WallThermal};
    use fs_exec::CancelGate;
    use fs_material::gas::GasSpec;

    fn spec(dims: [usize; 3], periodic: [bool; 3]) -> NsfSpec {
        NsfSpec {
            dims,
            periodic,
            wall: WallThermal::Adiabatic,
            gas: GasSpec::dry_air_ussa1976(),
            temperature_k: 293.15,
            pressure_pa: 101_325.0,
            dx_m: 5e-5,
            temperature_window_k: [273.15, 333.15],
            maximum_mach: 0.2,
            maximum_temperature_ratio: 0.2,
            bulk_to_shear_ratio: 0.73,
            bulk_ratio_source: "Cramer 2012, Phys. Fluids 24, 066102",
            maximum_cells: 1024,
            workers: 1,
            transport_points: 9,
        }
    }

    fn near(a: f64, b: f64, tolerance: f64) {
        assert!((a - b).abs() <= tolerance, "{a} != {b}");
    }

    #[test]
    fn accepted_pressure_density_misalignment_has_the_analytic_source() {
        // rho = 1 + a*x, p = theta0*(1+b*y); both gradients are
        // affine, so the centered stencil is exact at the interior probe.
        let grid = NsfGrid3::<D3q27>::new(
            spec([5; 3], [false; 3]),
            |_| false,
            |x| {
                let rho = 1.0 + 0.01 * (x[0] as f64 + 0.5);
                let pressure_factor = 1.0 + 0.01 * (x[1] as f64 + 0.5);
                (rho, [0.0; 3], 293.15 * pressure_factor / rho)
            },
        ).unwrap();
        let before = grid.state_hash();
        let q = grid.vorticity_budget([2; 3]).unwrap();
        let expected = 0.01 * 0.01 * THETA0 / 1.025_f64.powi(2);
        near(q.baroclinic[2], expected, 2e-17);
        near(q.baroclinic[0], 0.0, 2e-17);
        near(q.baroclinic[1], 0.0, 2e-17);
        for w in q.flow.vorticity() {
            near(w, 0.0, 2e-16);
        }
        let physical = q.to_physical(grid.law().dt_s).unwrap();
        near(physical.baroclinic[2] / ((expected / grid.law().dt_s) / grid.law().dt_s), 1.0, 2e-12);
        assert_eq!(before, grid.state_hash());
        assert_eq!(grid.steps(), 0);
        assert!(grid.vorticity_budget([0, 2, 2]).is_err());
        assert!(grid.vorticity_budget([5, 2, 2]).is_err());
    }

    #[test]
    fn periodic_sources_share_the_native_force_centered_curl() {
        let dims = [4, 8, 3];
        let grid = NsfGrid3::<D3q27>::new(
            spec(dims, [true; 3]),
            |_| false,
            |x| {
                let y = core::f64::consts::TAU * (x[1] as f64 + 0.5) / 8.0;
                (1.0, [0.02 * y.sin(), 0.0, 0.0], 293.15)
            },
        ).unwrap();
        let mut workspace = KinematicsWorkspace3::new(
            dims, KinematicsWorkspace3::required_bytes(dims).unwrap(),
        ).unwrap();
        let field = grid.kinematics(&CancelGate::new_clock_free(), &mut workspace).unwrap();
        for y in 0..8 {
            let cell = [1, y, 1];
            let q = grid.vorticity_budget(cell).unwrap();
            assert_eq!(q.flow, field[grid.index(cell)].unwrap());
            let angle = core::f64::consts::TAU * (y as f64 + 0.5) / 8.0;
            let expected = -0.02 * (core::f64::consts::TAU / 8.0).sin() * angle.cos();
            near(q.flow.vorticity()[2], expected, 2e-16);
            for rate in q.inviscid_rate() {
                near(rate, 0.0, 2e-17);
            }
        }
    }

    #[test]
    fn periodic_two_cell_axes_refuse_the_aliased_source_stencil() {
        for axis in 0..3 {
            let mut dims = [3; 3];
            dims[axis] = 2;
            let grid = NsfGrid3::<D3q27>::new(
                spec(dims, [true; 3]),
                |_| false,
                |x| (1.0, [0.01 * x[axis] as f64, 0.0, 0.0], 293.15),
            ).unwrap();
            let before = grid.state_hash();
            assert!(grid.vorticity_budget([0; 3]).is_err());
            assert_eq!(before, grid.state_hash());
        }
        let homogeneous = NsfGrid3::<D3q27>::new(
            spec([1, 3, 3], [true; 3]),
            |_| false,
            |_| (1.0, [0.0; 3], 293.15),
        ).unwrap();
        assert_eq!(homogeneous.vorticity_budget([0; 3]).unwrap().flow.vorticity(), [0.0; 3]);
    }

    #[test]
    fn a_solid_breaks_the_thermodynamic_stencil_without_inventing_wall_data() {
        let grid = NsfGrid3::<D3q27>::new(
            spec([5; 3], [true; 3]),
            |x| x == [3, 2, 2],
            |_| (1.0, [0.0; 3], 293.15),
        ).unwrap();
        let before = grid.state_hash();
        assert!(grid.vorticity_budget([2; 3]).is_err());
        assert!(grid.vorticity_budget([3, 2, 2]).is_err());
        assert!(grid.vorticity_budget([0; 3]).is_ok());
        assert_eq!(before, grid.state_hash());
    }
}
