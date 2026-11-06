//! Held spatial acceleration on the original compressible force/energy owner.
//! The field is prescribed in SI at cell centers and held during native steps.
//! Its momentum supply is rho*a; its mechanical power is rho*a.u.

use super::{
    NsfError, NsfGrid3, NsfLaw, check_state, half_force_velocity, raw_moments, rebase_force,
};
use crate::lattice3::Lattice3;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct BodyAcceleration {
    pub(super) acceleration_m_s2: [f64; 3],
    pub(super) acceleration_lattice: [f64; 3],
}

impl BodyAcceleration {
    fn resolve(law: &NsfLaw, acceleration_m_s2: [f64; 3]) -> Result<Self, NsfError> {
        if acceleration_m_s2.iter().any(|a| !a.is_finite()) {
            return Err(NsfError::Input(
                "body acceleration must have finite SI components",
            ));
        }
        let acceleration_lattice = acceleration_m_s2.map(|a| a * law.dt_s * (law.dt_s / law.dx_m));
        if acceleration_m_s2
            .iter()
            .zip(acceleration_lattice)
            .any(|(&a, value)| !value.is_finite() || (a != 0.0 && value == 0.0))
        {
            return Err(NsfError::Input(
                "body acceleration is not representable in lattice units",
            ));
        }
        Ok(Self {
            acceleration_m_s2,
            acceleration_lattice,
        })
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

fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), NsfError> {
    if cancelled() {
        Err(NsfError::Input("body acceleration update cancelled"))
    } else {
        Ok(())
    }
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Attach a held, spatially varying acceleration [m/s²] in grid axes.
    ///
    /// Supply one vector per native cell in x-fastest order. Solid entries
    /// are ignored. An empty slice or an all-zero fluid field disables it.
    /// The existing cell cap bounds its 48-byte-per-cell numerical storage;
    /// no field is stored for the disabled route. All fluid SI components
    /// must convert to representable lattice accelerations.
    ///
    /// The SAME accepted compressible density multiplies the field at every
    /// force epoch. The native energy corrector adds the matching endpoint
    /// trapezoid of rho*a.u, reported separately as body_force_work_j. Gravity
    /// and manufactured sources retain their own supplies; do not duplicate
    /// this field in those controls. There is no separate vorticity forcing.
    ///
    /// This is an externally prescribed, held field, not an implicit
    /// velocity-dependent force or a rotating-frame integrator. Call the
    /// setter between steps to change it. Spatial curl naturally generates
    /// vorticity through the original momentum update.
    pub fn with_body_acceleration_si(
        mut self,
        acceleration_m_s2: &[[f64; 3]],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, NsfError> {
        self.set_body_acceleration_si(acceleration_m_s2, cancelled)?;
        Ok(self)
    }

    /// Change the held acceleration without an instantaneous impulse or heat.
    ///
    /// Rebinds the full Guo half-source at fixed physical velocity and stress,
    /// then refreshes the original conservative flux history. A refusal or
    /// cancellation leaves all fields, source controls and accounts unchanged.
    /// Polls before staging, every 64 input/rebinding cells, around the native
    /// history reconstruction and immediately before publication.
    ///
    /// # Errors
    /// Wrong field length, allocation failure, cancellation, invalid SI or
    /// lattice values, or inability to preserve the finite accepted state.
    pub fn set_body_acceleration_si(
        &mut self,
        acceleration_m_s2: &[[f64; 3]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), NsfError> {
        checkpoint(&mut cancelled)?;
        if !acceleration_m_s2.is_empty() && acceleration_m_s2.len() != self.n {
            return Err(NsfError::Input(
                "body acceleration field length differs from the grid",
            ));
        }
        let mut staged = Vec::new();
        if !acceleration_m_s2.is_empty() {
            staged
                .try_reserve_exact(self.n)
                .map_err(|_| NsfError::Input("body acceleration field allocation failed"))?;
            let mut active = false;
            for (i, &acceleration) in acceleration_m_s2.iter().enumerate() {
                if i % 64 == 0 {
                    checkpoint(&mut cancelled)?;
                }
                let value = BodyAcceleration::resolve(
                    &self.law,
                    if self.solid[i] {
                        [0.0; 3]
                    } else {
                        acceleration
                    },
                )?;
                active |= value.acceleration_lattice != [0.0; 3];
                staged.push(value);
            }
            if !active {
                staged = Vec::new();
            }
        }
        if staged == self.body_acceleration {
            checkpoint(&mut cancelled)?;
            return Ok(());
        }
        self.rebind_held_acceleration(
            |grid, i| {
                let old = grid
                    .body_acceleration
                    .get(i)
                    .map_or([0.0; 3], |a| a.acceleration_lattice);
                let new = staged.get(i).map_or([0.0; 3], |a| a.acceleration_lattice);
                core::array::from_fn(|axis| new[axis] - old[axis])
            },
            &mut cancelled,
        )?;
        self.body_acceleration = staged;
        Ok(())
    }

    /// The original held SI value at a fluid cell; zero when disabled.
    pub fn body_acceleration_m_s2(&self, cell: [usize; 3]) -> Result<[f64; 3], NsfError> {
        if cell.iter().zip(self.spec.dims).any(|(&i, n)| i >= n) {
            return Err(NsfError::Input("body acceleration cell outside the grid"));
        }
        let i = self.index(cell);
        if self.solid[i] {
            return Err(NsfError::Input("solid cell has no body acceleration"));
        }
        Ok(self
            .body_acceleration
            .get(i)
            .map_or([0.0; 3], |a| a.acceleration_m_s2))
    }

    // Shared by uniform gravity and a held spatial body drive. Only the
    // force chart/history is published here; the caller then installs its
    // already admitted control, with no remaining fallible work.
    pub(super) fn rebind_held_acceleration(
        &mut self,
        delta: impl Fn(&Self, usize) -> [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), NsfError> {
        checkpoint(cancelled)?;
        let (rho, momentum) = self.macro_fields(&self.f);
        let velocity = self.velocity_field(&rho, &momentum, &self.force);
        checkpoint(cancelled)?;
        let rest = L::E
            .iter()
            .position(|direction| *direction == [0; 3])
            .ok_or(NsfError::Input(
                "acceleration rebinding requires a rest population",
            ))?;
        let mut force = self.force.clone();
        let mut f = self.f.clone();
        for i in 0..self.n {
            if i % 64 == 0 {
                checkpoint(cancelled)?;
            }
            if self.solid[i] {
                continue;
            }
            let acceleration = delta(self, i);
            for a in 0..3 {
                force[i][a] += rho[i] * acceleration[a];
            }
            let difference = core::array::from_fn(|a| self.force[i][a] - force[i][a]);
            if !force[i].iter().chain(&difference).all(|v| v.is_finite()) {
                return Err(NsfError::NonFinite { cell: i });
            }
            rebase_force::<L>(&mut f[i], velocity[i], difference);
            if f[i].iter().any(|v| !v.is_finite()) {
                return Err(NsfError::NonFinite { cell: i });
            }
            // The half-source has zero mass analytically. The rest population
            // retains summation roundoff without changing momentum or stress.
            let mass_correction = rho[i] - f[i].iter().sum::<f64>();
            f[i][rest] += mass_correction;
            let (r, m) = raw_moments::<L>(&f[i]);
            let u = half_force_velocity(r, m, force[i]);
            let gamma = self
                .material
                .get(i)
                .and_then(Option::as_ref)
                .map_or(self.law.gamma, |material| material.state.gamma);
            check_state(&self.spec, gamma, i, r, u, self.theta[i])?;
            if (r - rho[i]).abs() > 128.0 * f64::EPSILON * rho[i]
                || u.iter().zip(velocity[i]).any(|(&actual, expected)| {
                    (actual - expected).abs() > 128.0 * f64::EPSILON * expected.abs().max(1.0)
                })
            {
                return Err(NsfError::Input(
                    "acceleration chart cannot retain the physical state",
                ));
            }
        }
        checkpoint(cancelled)?;
        let stress = self.stress_field(
            &f,
            &force,
            &rho,
            &self.theta,
            &self.material,
            &self.sgs_transport,
        );
        let (div, wall, open) = self.flux_divergence(
            &self.e,
            &rho,
            &velocity,
            &self.theta,
            &stress,
            &self.sgs_transport,
            self.open_faces.as_ref(),
        )?;
        if !wall.is_finite() || div.iter().any(|value| !value.is_finite()) {
            return Err(NsfError::Input(
                "acceleration force-chart flux history overflowed",
            ));
        }
        checkpoint(cancelled)?;
        self.f = f;
        self.force = force;
        self.div = div;
        self.wall_out = wall;
        self.open_out = open;
        Ok(())
    }
}

#[cfg(test)]
mod tests;
