//! Uniform physical gravity on the compressible momentum and energy owners.
//!
//! Momentum receives rho*g; total energy receives rho*g.u. There is no
//! Boussinesq density subtraction: pressure and density remain the original
//! compressible fields. The same endpoint quadrature supplies the energy
//! update and its independent work account. Gravity does not heat the gas.
//! Guo population rebinding follows https://doi.org/10.1103/PhysRevE.65.046308.

use super::{NsfError, NsfGrid3, NsfLaw, check_state, raw_moments, rebase_force};
use crate::lattice3::Lattice3;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Gravity {
    pub(super) acceleration_m_s2: [f64; 3],
    pub(super) acceleration_lattice: [f64; 3],
}

impl Gravity {
    fn resolve(law: &NsfLaw, acceleration_m_s2: [f64; 3]) -> Result<Option<Self>, NsfError> {
        if acceleration_m_s2.iter().any(|g| !g.is_finite()) {
            return Err(NsfError::Input("gravity must have finite SI components"));
        }
        if acceleration_m_s2 == [0.0; 3] {
            return Ok(None);
        }
        let acceleration_lattice =
            acceleration_m_s2.map(|g| g * law.dt_s * (law.dt_s / law.dx_m));
        if acceleration_m_s2.iter().zip(acceleration_lattice).any(|(&g, a)| {
            !a.is_finite() || (g != 0.0 && a == 0.0)
        }) {
            return Err(NsfError::Input("gravity is not representable in the original lattice units"));
        }
        Ok(Some(Self { acceleration_m_s2, acceleration_lattice }))
    }

    pub(super) fn work(self, rho0: f64, u0: [f64; 3], rho1: f64, u1: [f64; 3]) -> f64 {
        let power = |rho: f64, u: [f64; 3]| {
            rho * (self.acceleration_lattice[0] * u[0]
                + self.acceleration_lattice[1] * u[1]
                + self.acceleration_lattice[2] * u[2])
        };
        0.5 * (power(rho0, u0) + power(rho1, u1))
    }
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Attach constant signed gravitational acceleration [m/s²] in grid axes.
    ///
    /// Uses the SAME compressible density as inertia and the EOS at every
    /// force epoch, including thermal and characteristic correctors. The
    /// energy equation receives the matching trapezoid of rho*g.u, separately
    /// reported by `EnergyLedger::gravity_work_j`. Do not duplicate these
    /// loads in a manufactured momentum or energy source.
    ///
    /// The existing Mach, material and thermal bounds still apply. This is
    /// not a hydrostatic reconstruction or an exact discrete conservation
    /// law for fluid plus gravitational potential energy. A constant force
    /// on a periodic axis describes an externally driven periodic domain;
    /// no single-valued periodic gravitational potential is implied.
    ///
    /// # Errors
    /// Refuses nonfinite/unrepresentable gravity or a force chart whose
    /// physical state and finite flux history cannot be retained.
    pub fn with_gravity_si(mut self, acceleration_m_s2: [f64; 3]) -> Result<Self, NsfError> {
        self.set_gravity_si(acceleration_m_s2)?;
        Ok(self)
    }

    /// Accepted signed gravitational acceleration [m/s²]; zero when disabled.
    #[must_use]
    pub const fn gravity_m_s2(&self) -> [f64; 3] {
        match self.gravity {
            Some(gravity) => gravity.acceleration_m_s2,
            None => [0.0; 3],
        }
    }

    /// Change held gravity without an instantaneous impulse or heat transfer.
    ///
    /// Rebinds the full Guo half-source at fixed physical velocity and stress.
    /// Existing thermal/SGS/source histories are retained; the conservative
    /// flux history is refreshed from that same physical state. All
    /// admission precedes publication. Passing zero to a gravity-free grid
    /// is an exact no-op. The next step uses the new held acceleration.
    ///
    /// # Errors
    /// The same admission failures as [`Self::with_gravity_si`]; a refusal
    /// retains all previously accepted fields, controls and accounts.
    pub fn set_gravity_si(&mut self, acceleration_m_s2: [f64; 3]) -> Result<(), NsfError> {
        let gravity = Gravity::resolve(&self.law, acceleration_m_s2)?;
        if gravity == self.gravity {
            return Ok(());
        }
        let old = self.gravity.map_or([0.0; 3], |g| g.acceleration_lattice);
        let new = gravity.map_or([0.0; 3], |g| g.acceleration_lattice);
        let (rho, momentum) = self.macro_fields(&self.f);
        let velocity = self.velocity_field(&rho, &momentum, &self.force);
        let rest = L::E.iter().position(|direction| *direction == [0; 3])
            .ok_or(NsfError::Input("gravity rebinding requires a rest population"))?;
        let mut force = self.force.clone();
        let mut f = self.f.clone();
        for i in 0..self.n {
            if self.solid[i] {
                continue;
            }
            for a in 0..3 {
                force[i][a] += rho[i] * (new[a] - old[a]);
            }
            let difference = core::array::from_fn(|a| self.force[i][a] - force[i][a]);
            if !force[i].iter().chain(&difference).all(|v| v.is_finite()) {
                return Err(NsfError::NonFinite { cell: i });
            }
            rebase_force::<L>(&mut f[i], velocity[i], difference);
            if f[i].iter().any(|v| !v.is_finite()) {
                return Err(NsfError::NonFinite { cell: i });
            }
            // The analytic source has zero mass. Keep its summation roundoff
            // in the rest population, which carries neither momentum nor stress.
            let mass_correction = rho[i] - f[i].iter().sum::<f64>();
            f[i][rest] += mass_correction;
            let (r, m) = raw_moments::<L>(&f[i]);
            let u = super::half_force_velocity(r, m, force[i]);
            let gamma = self.material.get(i).and_then(Option::as_ref)
                .map_or(self.law.gamma, |material| material.state.gamma);
            check_state(&self.spec, gamma, i, r, u, self.theta[i])?;
            if (r - rho[i]).abs() > 128.0 * f64::EPSILON * rho[i]
                || u.iter().zip(velocity[i]).any(|(&actual, expected)| {
                    (actual - expected).abs()
                        > 128.0 * f64::EPSILON * expected.abs().max(1.0)
                })
            {
                return Err(NsfError::Input("gravity force chart cannot retain the physical state"));
            }
        }
        let stress = self.stress_field(
            &f, &force, &rho, &self.theta, &self.material, &self.sgs_transport,
        );
        let (div, wall, open) = self.flux_divergence(
            &self.e, &rho, &velocity, &self.theta, &stress,
            &self.sgs_transport, self.open_faces.as_ref(),
        )?;
        if !wall.is_finite() || div.iter().any(|value| !value.is_finite()) {
            return Err(NsfError::Input("gravity force-chart flux history overflowed"));
        }
        self.gravity = gravity;
        self.f = f;
        self.force = force;
        self.div = div;
        self.wall_out = wall;
        self.open_out = open;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q27::D3q27;
    use crate::nsf3::{NsfSpec, WallThermal};
    use fs_material::gas::GasSpec;

    #[test]
    fn rebind_retains_evolved_force_corrected_stress() {
        let mut flow = NsfGrid3::<D3q27>::new(NsfSpec {
            dims: [4, 3, 2], periodic: [true; 3],
            gas: GasSpec::dry_air_ussa1976(), temperature_k: 300.0,
            pressure_pa: 100_000.0, dx_m: 0.001,
            temperature_window_k: [250.0, 400.0], maximum_mach: 0.2,
            maximum_temperature_ratio: 0.3, bulk_to_shear_ratio: 0.73,
            bulk_ratio_source: "Cramer 2012 N2 ratio; force chart regression",
            maximum_cells: 24, wall: WallThermal::Adiabatic,
            workers: 1, transport_points: 17,
        }, |_| false, |x| (1.0, [0.01, 0.002 * x[0] as f64, -0.003], 300.0))
            .unwrap();
        for _ in 0..3 { flow.step().unwrap(); }
        let (rho, _) = flow.macro_fields(&flow.f);
        let before = flow.stress_field(
            &flow.f, &flow.force, &rho, &flow.theta, &flow.material, &flow.sgs_transport,
        );
        assert!(before.iter().flatten().any(|stress| stress.abs() > 1e-12));
        for gravity in [[1.0, -9.80665, 2.0], [-3.0, 5.0, -7.0], [0.0; 3]] {
            flow.set_gravity_si(gravity).unwrap();
            let (rho, _) = flow.macro_fields(&flow.f);
            let after = flow.stress_field(
                &flow.f, &flow.force, &rho, &flow.theta, &flow.material, &flow.sgs_transport,
            );
            for (old, new) in before.iter().flatten().zip(after.iter().flatten()) {
                assert!((old - new).abs() < 5e-14, "{old:e} != {new:e}");
            }
        }
    }
}
