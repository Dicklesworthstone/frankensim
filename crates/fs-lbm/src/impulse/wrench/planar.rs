//! Planar force-and-couple deposition on the represented liquid inertia.

use super::super::{ImpulseControls, ImpulseError, ImpulseReport, poll};
use crate::{E, Q};
use fs_exec::{CancelGate, reduce::Compensated};

/// An instantaneous planar wrench in lattice units, applied TO the fluid.
/// Linear and angular impulses are per unit depth of the two-dimensional owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WrenchImpulse2 {
    /// Point about which angular impulse is specified, in lattice cell lengths.
    pub origin: [f64; 2],
    /// Total linear impulse, not an impulse per cell or a force.
    pub linear: [f64; 2],
    /// Total out-of-plane angular impulse; positive is counterclockwise.
    pub angular: f64,
}

/// Native impulse limits and an independent angular-momentum allowance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WrenchImpulseControls2 {
    pub impulse: ImpulseControls,
    /// Deposition defect plus local/global angular readback defects.
    pub maximum_angular_momentum_defect: f64,
}

/// Accepted linear/angular loading and signed kinetic work. A coupled body's
/// opposite reaction remains the responsibility of its own momentum owner.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct WrenchImpulseReport2 {
    pub requested: WrenchImpulse2,
    pub impulse: ImpulseReport<2>,
    pub stencil_mass: f64,
    pub stencil_center_of_mass: [f64; 2],
    /// Moment of the prescribed cell impulses about the requested origin.
    pub applied_angular_impulse: f64,
    /// Independent physical momentum readback at the original cell centers.
    pub angular_momentum_change: f64,
    /// Deposition defect plus the native local/global linear readback defect.
    pub absolute_linear_momentum_defect: f64,
    pub absolute_angular_momentum_defect: f64,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct WrenchCell2 {
    pub index: usize,
    pub position: [f64; 2],
    pub mass: f64,
}

/// Scalar planar inertia solve; only the caller's compact cell stencil is held.
pub(crate) struct WrenchPlan2 {
    requested: WrenchImpulse2,
    mass: f64,
    anchor: [f64; 2],
    center_relative: [f64; 2],
    translation: [f64; 2],
    rotation: f64,
}

impl WrenchPlan2 {
    pub(crate) fn prepare(
        cells: &[WrenchCell2],
        requested: WrenchImpulse2,
        controls: WrenchImpulseControls2,
        gate: &CancelGate,
    ) -> Result<Self, ImpulseError> {
        poll(gate)?;
        controls.impulse.admit(cells.len())?;
        if !controls.maximum_angular_momentum_defect.is_finite()
            || controls.maximum_angular_momentum_defect < 0.0
            || !requested
                .origin
                .into_iter()
                .chain(requested.linear)
                .chain([requested.angular])
                .all(f64::is_finite)
        {
            return Err(ImpulseError::Invalid(
                "invalid planar wrench or angular balance budget",
            ));
        }
        let anchor = cells[0].position;
        let mut mass = Compensated::zero();
        let mut first = [Compensated::zero(); 2];
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
            for a in 0..2 {
                first[a] = first[a].accumulate(cell.mass * (cell.position[a] - anchor[a]));
            }
        }
        let mass = mass.value();
        let center_relative = first.map(|sum| sum.value() / mass);
        let translation = requested.linear.map(|impulse| impulse / mass);
        if !mass.is_finite()
            || mass <= 0.0
            || !center_relative
                .into_iter()
                .chain(translation)
                .all(f64::is_finite)
        {
            return Err(ImpulseError::Invalid(
                "unrepresentable planar stencil inertia",
            ));
        }
        let mut inertia = Compensated::zero();
        for (i, cell) in cells.iter().enumerate() {
            if i % 64 == 0 {
                poll(gate)?;
            }
            let r: [f64; 2] =
                core::array::from_fn(|a| (cell.position[a] - anchor[a]) - center_relative[a]);
            inertia = inertia.accumulate(cell.mass * r[0].mul_add(r[0], r[1] * r[1]));
        }
        let inertia = inertia.value();
        let arm = core::array::from_fn(|a| (anchor[a] - requested.origin[a]) + center_relative[a]);
        let angular = requested.angular - cross(arm, requested.linear);
        let rotation =
            if inertia == 0.0 && angular.abs() <= 0.25 * controls.maximum_angular_momentum_defect {
                0.0
            } else {
                angular / inertia
            };
        if !inertia.is_finite() || inertia < 0.0 || !rotation.is_finite() {
            return Err(ImpulseError::Invalid(
                "selected fluid cells cannot resolve the requested planar angular impulse",
            ));
        }
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

    pub(crate) fn impulse(&self, cell: WrenchCell2) -> [f64; 2] {
        let r: [f64; 2] =
            core::array::from_fn(|a| (cell.position[a] - self.anchor[a]) - self.center_relative[a]);
        [
            cell.mass * (-self.rotation).mul_add(r[1], self.translation[0]),
            cell.mass * self.rotation.mul_add(r[0], self.translation[1]),
        ]
    }

    /// Read the staged physical momenta before the population owner's swap.
    pub(crate) fn finish(
        &self,
        cells: &[WrenchCell2],
        report: ImpulseReport<2>,
        controls: WrenchImpulseControls2,
        gate: &CancelGate,
        mut momentum_change: impl FnMut(usize) -> [f64; 2],
    ) -> Result<WrenchImpulseReport2, ImpulseError> {
        let mut applied = Compensated::zero();
        let mut change = Compensated::zero();
        let mut local_defect = Compensated::zero();
        for (i, &cell) in cells.iter().enumerate() {
            if i % 64 == 0 {
                poll(gate)?;
            }
            let arm = core::array::from_fn(|a| cell.position[a] - self.requested.origin[a]);
            let prescribed = cross(arm, self.impulse(cell));
            let observed = cross(arm, momentum_change(cell.index));
            applied = applied.accumulate(prescribed);
            change = change.accumulate(observed);
            local_defect = local_defect.accumulate((observed - prescribed).abs());
        }
        let applied = applied.value();
        let change = change.value();
        let linear_defect = report.absolute_momentum_defect
            + (report.applied_impulse[0] - self.requested.linear[0]).abs()
            + (report.applied_impulse[1] - self.requested.linear[1]).abs();
        let angular_defect = local_defect.value()
            + (change - applied).abs()
            + (applied - self.requested.angular).abs();
        let center = core::array::from_fn(|a| self.anchor[a] + self.center_relative[a]);
        if !center
            .into_iter()
            .chain([applied, change, linear_defect, angular_defect])
            .all(f64::is_finite)
            || linear_defect > controls.impulse.maximum_momentum_defect
            || angular_defect > controls.maximum_angular_momentum_defect
        {
            return Err(ImpulseError::Balance);
        }
        poll(gate)?;
        Ok(WrenchImpulseReport2 {
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

fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0].mul_add(b[1], -a[1] * b[0])
}

/// Same Guo-frame readback as the shared impulse kernel. Liquid owners supply
/// original held mass, whereas a fully occupied grid owns population momentum.
pub(crate) fn planar_momentum_change(
    before: &[f64; Q],
    after: &[f64; Q],
    acceleration: [f64; 2],
    external_force: [f64; 2],
    held_mass: Option<f64>,
) -> [f64; 2] {
    let observe = |populations: &[f64; Q]| {
        let rho = populations.iter().sum::<f64>();
        let directions = E.map(|(x, y)| [f64::from(x), f64::from(y)]);
        let momentum = core::array::from_fn::<_, 2, _>(|a| {
            populations
                .iter()
                .zip(directions)
                .map(|(f, e)| f * e[a])
                .sum::<f64>()
                + 0.5 * acceleration[a].mul_add(rho, external_force[a])
        });
        (rho, momentum)
    };
    let (rho_before, before) = observe(before);
    let (rho_after, after) = observe(after);
    core::array::from_fn(|a| match held_mass {
        Some(mass) => mass * (after[a] / rho_after - before[a] / rho_before),
        None => after[a] - before[a],
    })
}
