//! Prescribed spatial acceleration on the original compressible force/energy owner.
//! SI cell-center fields may be held or supplied at the next native time level.
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

// Borrows the accepted OLD field and, for an unsteady step, the independently
// admitted NEW field. Collision still uses the accepted self.force; only
// endpoint correctors consume field(). No callback or solver copy is involved.
pub(super) struct BodyForceStep<'a> {
    old: &'a [BodyAcceleration],
    new: Option<&'a [BodyAcceleration]>,
}

impl<'a> BodyForceStep<'a> {
    pub(super) fn new(
        old: &'a [BodyAcceleration],
        new: Option<&'a [BodyAcceleration]>,
    ) -> Self {
        Self { old, new }
    }

    pub(super) fn field(&self) -> &'a [BodyAcceleration] {
        self.new.unwrap_or(self.old)
    }

    pub(super) fn active(&self) -> bool {
        !self.old.is_empty() || !self.field().is_empty()
    }

    pub(super) fn work(&self, cell: usize, rho0: f64, u0: [f64; 3], rho1: f64, u1: [f64; 3]) -> f64 {
        let old = self.old.get(cell);
        let Some(new) = self.new else {
            // Keep every operation of the held route, including signed zero.
            return old.map_or(0.0, |a| a.work(rho0, u0, rho1, u1));
        };
        let a0 = old.map_or([0.0; 3], |a| a.acceleration_lattice);
        let a1 = new.get(cell).map_or([0.0; 3], |a| a.acceleration_lattice);
        // Power is acceleration dotted with physical momentum, rho*u.
        // Center the two factors separately: a_bar . (m0 + m1)/2.
        // For constant rho, the native Guo impulse is rho*(a0+a1)/2,
        // so this is exactly its kinetic work, including force reversal.
        // An endpoint-power trapezoid instead adds Delta(a).Delta(m)/4
        // as spurious heat. This centered quadrature is second order and
        // continuously reduces to the original held-field work.
        let acceleration: [f64; 3] = core::array::from_fn(|a| 0.5 * a0[a] + 0.5 * a1[a]);
        let power = |rho: f64, u: [f64; 3]| {
            rho * (acceleration[0] * u[0] + acceleration[1] * u[1] + acceleration[2] * u[2])
        };
        0.5 * (power(rho0, u0) + power(rho1, u1))
    }
}

pub(super) fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), NsfError> {
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
    /// This is an externally prescribed field, not an implicit
    /// velocity-dependent force or a rotating-frame integrator. Use the
    /// setter to rebind it at fixed physical state; use
    /// [`Self::step_with_body_acceleration_si`] to advance a time-varying
    /// drive with distinct old/new force epochs. Spatial curl naturally
    /// generates vorticity through the original momentum update.
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
        let staged = self.resolve_body_acceleration(acceleration_m_s2, &mut cancelled)?;
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

    /// Advance one native step using the prescribed acceleration at its endpoint.
    ///
    /// The accepted field is the OLD acceleration and `acceleration_m_s2`
    /// is the NEW one, one vector per cell in x-fastest order. Empty/all-zero
    /// fluid fields mean zero acceleration; solid entries are ignored.
    /// The native collision retains the OLD force, while every endpoint
    /// thermal/SGS/contact corrector uses NEW `rho*a`. After acceptance the
    /// NEW field remains bound for subsequent held steps. Identical endpoint
    /// fields execute the original held arithmetic exactly.
    ///
    /// Total energy receives `(a_old+a_new).(rho_old*u_old+rho_new*u_new)/4`
    /// in lattice units, using physical half-force velocities. This
    /// second-order centered mechanical work reduces to the held quadrature
    /// and equals the kinetic change for uniform constant-density forcing,
    /// without artificial heating under a changing or reversing drive.
    /// It is reported once in `EnergyLedger::body_force_work_j`.
    ///
    /// The input is admitted once; there are no acceleration callbacks,
    /// extra source samples, grid copies or separate vorticity integrators.
    /// Gravity and manufactured sources retain their existing supplies.
    /// Time samples must resolve the prescribed drive; discontinuities,
    /// implicit velocity dependence and rotating-frame dynamics require
    /// their own integration treatment.
    ///
    /// # Errors
    /// Field-length, SI conversion, allocation, cancellation and all native
    /// step refusals leave the accepted state and acceleration unchanged.
    /// Cancellation is polled through field admission, native step stages,
    /// corrector cells and immediately before publication. As for `step`,
    /// an external source callback's own side effects cannot be rolled back.
    pub fn step_with_body_acceleration_si(
        &mut self,
        acceleration_m_s2: &[[f64; 3]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<super::EnergyLedger, NsfError> {
        let next = self.resolve_body_acceleration(acceleration_m_s2, &mut cancelled)?;
        let next = if next == self.body_acceleration {
            None
        } else {
            Some(next)
        };
        self.step_body_acceleration(next, &mut cancelled)
    }

    fn resolve_body_acceleration(
        &self,
        acceleration_m_s2: &[[f64; 3]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<BodyAcceleration>, NsfError> {
        checkpoint(cancelled)?;
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
                    checkpoint(cancelled)?;
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
        Ok(staged)
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
