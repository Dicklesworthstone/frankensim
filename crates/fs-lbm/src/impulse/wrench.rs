//! Local force-and-couple deposition on actual fluid inertia.
//!
//! A rigid velocity increment is the least mass-weighted change that realizes
//! a prescribed linear and angular impulse on a fixed cloud of cell centers.
//! The population owner retains responsibility for its Guo frame, speed/stress
//! admission, signed work, and transactional publication.

use super::{ImpulseControls, ImpulseError, ImpulseReport, poll};
use crate::angular_momentum::angular_correction;
use fs_exec::{CancelGate, reduce::Compensated};

/// An instantaneous wrench in lattice units, applied TO the fluid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WrenchImpulse3 {
    /// Point about which angular impulse is specified, in lattice cell lengths.
    pub origin: [f64; 3],
    /// Total linear impulse; this is not an impulse per cell or a force.
    pub linear: [f64; 3],
    /// Total angular impulse about origin, including any free couple.
    pub angular: [f64; 3],
}

/// Native impulse limits and a separate physical torque-readback allowance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WrenchImpulseControls3 {
    pub impulse: ImpulseControls,
    /// L1 deposition defect plus local/global angular readback defects.
    pub maximum_angular_momentum_defect: f64,
}

/// The actual accepted cell impulses, torque and work. No body reaction is
/// applied implicitly; its owner must receive the opposite wrench.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct WrenchImpulseReport3 {
    pub requested: WrenchImpulse3,
    /// Existing population-level linear impulse, mass and signed kinetic work.
    pub impulse: ImpulseReport<3>,
    pub stencil_mass: f64,
    pub stencil_center_of_mass: [f64; 3],
    /// Torque of the deposited cell impulses about the requested origin.
    pub applied_angular_impulse: [f64; 3],
    /// Independently read physical angular momentum change at those positions.
    pub angular_momentum_change: [f64; 3],
    /// Deposition defect plus the native local/global linear readback defect.
    pub absolute_linear_momentum_defect: f64,
    /// Deposition defect plus local/global angular readback defects.
    pub absolute_angular_momentum_defect: f64,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct WrenchCell3 {
    pub index: usize,
    pub position: [f64; 3],
    pub mass: f64,
}

/// Fixed-size flow map; the caller owns the compact original-cell stencil.
pub(crate) struct WrenchPlan3 {
    requested: WrenchImpulse3,
    mass: f64,
    anchor: [f64; 3],
    center_relative: [f64; 3],
    translation: [f64; 3],
    rotation: [f64; 3],
}

impl WrenchPlan3 {
    pub(crate) fn prepare(
        cells: &[WrenchCell3],
        requested: WrenchImpulse3,
        controls: WrenchImpulseControls3,
        gate: &CancelGate,
    ) -> Result<Self, ImpulseError> {
        poll(gate)?;
        controls.impulse.admit(cells.len())?;
        if cells.is_empty() || cells.len() > controls.impulse.max_cells {
            return Err(ImpulseError::Invalid(
                "wrench needs a nonempty bounded cell stencil",
            ));
        }
        if !controls.maximum_angular_momentum_defect.is_finite()
            || controls.maximum_angular_momentum_defect < 0.0
            || !requested
                .origin
                .into_iter()
                .chain(requested.linear)
                .chain(requested.angular)
                .all(f64::is_finite)
        {
            return Err(ImpulseError::Invalid(
                "invalid wrench or angular balance budget",
            ));
        }
        let anchor = cells[0].position;
        let mut mass = Compensated::zero();
        let mut first = [Compensated::zero(); 3];
        for (i, cell) in cells.iter().enumerate() {
            if i % 64 == 0 {
                poll(gate)?;
            }
            if !cell.mass.is_finite()
                || cell.mass <= 0.0
                || !cell.position.into_iter().all(f64::is_finite)
                || (i != 0 && cells[i - 1].index >= cell.index)
            {
                return Err(ImpulseError::Invalid(
                    "wrench cells need positive actual inertia and increasing unique indices",
                ));
            }
            mass = mass.accumulate(cell.mass);
            for a in 0..3 {
                first[a] = first[a].accumulate(cell.mass * (cell.position[a] - anchor[a]));
            }
        }
        let mass = mass.value();
        let center_relative = first.map(|s| s.value() / mass);
        let translation = requested.linear.map(|j| j / mass);
        if !mass.is_finite()
            || mass <= 0.0
            || !center_relative
                .into_iter()
                .chain(translation)
                .all(f64::is_finite)
        {
            return Err(ImpulseError::Invalid(
                "unrepresentable wrench stencil inertia",
            ));
        }
        let mut inertia = [Compensated::zero(); 9];
        for (i, cell) in cells.iter().enumerate() {
            if i % 64 == 0 {
                poll(gate)?;
            }
            let r: [f64; 3] =
                core::array::from_fn(|a| (cell.position[a] - anchor[a]) - center_relative[a]);
            for a in 0..3 {
                for b in a..3 {
                    let value = if a == b {
                        let u = r[(a + 1) % 3];
                        let v = r[(a + 2) % 3];
                        cell.mass * u.mul_add(u, v * v)
                    } else {
                        -(cell.mass * r[a]) * r[b]
                    };
                    inertia[3 * a + b] = inertia[3 * a + b].accumulate(value);
                    if a != b {
                        inertia[3 * b + a] = inertia[3 * b + a].accumulate(value);
                    }
                }
            }
        }
        let arm = core::array::from_fn(|a| (anchor[a] - requested.origin[a]) + center_relative[a]);
        let orbital = cross(arm, requested.linear);
        let angular = core::array::from_fn(|a| requested.angular[a] - orbital[a]);
        let rotation = angular_correction(
            &inertia.map(Compensated::value),
            angular,
            0.25 * controls.maximum_angular_momentum_defect,
        )
        .ok_or(ImpulseError::Invalid(
            "selected fluid cells cannot resolve the requested angular impulse",
        ))?;
        poll(gate)?;
        Ok(Self {
            requested,
            mass,
            anchor,
            center_relative,
            translation,
            rotation,
        })
    }

    pub(crate) fn impulse(&self, cell: WrenchCell3) -> [f64; 3] {
        let r =
            core::array::from_fn(|a| (cell.position[a] - self.anchor[a]) - self.center_relative[a]);
        let rotation = cross(self.rotation, r);
        core::array::from_fn(|a| cell.mass * (self.translation[a] + rotation[a]))
    }

    /// Observe original and staged physical momenta before either owner changes.
    /// The observer supplies the SAME force-corrected readback as its native
    /// impulse report, including any floating-point source-frame rounding.
    pub(crate) fn finish(
        &self,
        cells: &[WrenchCell3],
        report: ImpulseReport<3>,
        controls: WrenchImpulseControls3,
        gate: &CancelGate,
        mut momentum_change: impl FnMut(usize) -> [f64; 3],
    ) -> Result<WrenchImpulseReport3, ImpulseError> {
        let mut applied = [Compensated::zero(); 3];
        let mut change = [Compensated::zero(); 3];
        let mut local_defect = Compensated::zero();
        for (i, &cell) in cells.iter().enumerate() {
            if i % 64 == 0 {
                poll(gate)?;
            }
            let arm = core::array::from_fn(|a| cell.position[a] - self.requested.origin[a]);
            let prescribed = cross(arm, self.impulse(cell));
            let observed = cross(arm, momentum_change(cell.index));
            for a in 0..3 {
                applied[a] = applied[a].accumulate(prescribed[a]);
                change[a] = change[a].accumulate(observed[a]);
                local_defect = local_defect.accumulate((observed[a] - prescribed[a]).abs());
            }
        }
        let applied = applied.map(Compensated::value);
        let change = change.map(Compensated::value);
        let linear_defect = report.absolute_momentum_defect
            + l1_difference(report.applied_impulse, self.requested.linear);
        let angular_defect = local_defect.value()
            + l1_difference(change, applied)
            + l1_difference(applied, self.requested.angular);
        let center = core::array::from_fn(|a| self.anchor[a] + self.center_relative[a]);
        if !center
            .into_iter()
            .chain(applied)
            .chain(change)
            .chain([linear_defect, angular_defect])
            .all(f64::is_finite)
            || linear_defect > controls.impulse.maximum_momentum_defect
            || angular_defect > controls.maximum_angular_momentum_defect
        {
            return Err(ImpulseError::Balance);
        }
        poll(gate)?;
        Ok(WrenchImpulseReport3 {
            requested: self.requested,
            impulse: report,
            stencil_mass: self.mass,
            stencil_center_of_mass: center,
            applied_angular_impulse: applied,
            angular_momentum_change: change,
            absolute_linear_momentum_defect: linear_defect,
            absolute_angular_momentum_defect: angular_defect,
        })
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -a[2] * b[1]),
        a[2].mul_add(b[0], -a[0] * b[2]),
        a[0].mul_add(b[1], -a[1] * b[0]),
    ]
}

fn l1_difference(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0]).abs() + (a[1] - b[1]).abs() + (a[2] - b[2]).abs()
}

#[cfg(test)]
mod tests;
