//! Capacity-bounded signed mass and momenta on native D2Q9 conversions.

use super::*;
use fs_exec::CancelGate;
use std::collections::BTreeMap;

mod angular;

/// Local numerical conversion policy for conservative momentum stepping.
/// Legacy native step and pressure-only stepping retain their original
/// conversion behavior; this selects the explicitly conservative owner.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ConversionRedistribution2 {
    /// Preserve signed mass and physical linear momentum on adjacent routes.
    #[default]
    ConservativeLocal,
    /// Also preserve out-of-plane angular momentum independently in each
    /// route-connected region. Periodic domains and unrepresentable spin refuse.
    ConservativeLocalAngular,
}

/// Conversion balance in lattice units, excluding collision, streaming,
/// pressure traction and imposed forces. These are numerical remap terms.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ConversionMomentumReport2 {
    /// Sum of absolute mass transfers between adjacent interface cells.
    pub redistributed_mass: f64,
    /// Absolute source mass remainders, plus component closure in angular mode.
    pub absolute_mass_defect: f64,
    /// Absolute local physical momentum discrepancies, plus component closure
    /// after the angular correction when that mode is selected.
    pub absolute_momentum_defect: [f64; 2],
    /// Absolute out-of-plane angular defects summed over actual route-connected
    /// regions, about each region's first row-major cell centre. Linear
    /// residuals are separate; translating the angular residual to another
    /// origin includes their moment arm. None for the default linear policy;
    /// Some(0.0) for angular mode without transfers.
    pub absolute_angular_momentum_defect: Option<f64>,
    /// Signed resolved kinetic change during conversion, never physical heat.
    pub kinetic_energy_change: f64,
    /// Nonzero adjacent-cell transfers. Periodic aliases count only once.
    pub routes: usize,
}

/// A conservative proposal failed; the accepted surface remains unchanged.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FreeSurfaceMomentumError2 {
    /// The existing cancellation gate was requested before publication.
    Cancelled,
    /// Invalid state, insufficient local capacity, or unresolved momentum.
    Invalid(&'static str),
}
impl core::fmt::Display for FreeSurfaceMomentumError2 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "D2Q9 conservative free surface: {self:?}")
    }
}
impl std::error::Error for FreeSurfaceMomentumError2 {}

pub(super) fn poll(gate: &CancelGate) -> Result<(), FreeSurfaceMomentumError2> {
    if gate.is_requested() {
        Err(FreeSurfaceMomentumError2::Cancelled)
    } else {
        Ok(())
    }
}
fn invalid(reason: &'static str) -> FreeSurfaceMomentumError2 {
    FreeSurfaceMomentumError2::Invalid(reason)
}
fn admit_velocity(u: [f64; 2]) -> Result<(), FreeSurfaceMomentumError2> {
    if u[0].hypot(u[1]) / CS2.sqrt() < crate::MACH_LIMIT {
        Ok(())
    } else {
        Err(invalid(
            "conversion exceeds the finite low-Mach lattice domain; refine",
        ))
    }
}

pub(super) fn admit_initialized_cell(
    populations: &[f64; Q],
    force: [f64; 2],
) -> Result<(), FreeSurfaceMomentumError2> {
    let rho = populations.iter().sum::<f64>();
    if !rho.is_finite()
        || rho <= 0.0
        || populations.iter().any(|f| !f.is_finite() || *f <= 0.0)
        || force.iter().any(|f| !f.is_finite())
    {
        return Err(invalid(
            "new interface populations or held force are unrepresentable",
        ));
    }
    let u = core::array::from_fn(|d| {
        let raw = populations
            .iter()
            .zip(E)
            .map(|(f, (x, y))| f * f64::from([x, y][d]))
            .sum::<f64>();
        (raw + 0.5 * force[d]) / rho
    });
    admit_velocity(u)
}

impl FreeSurface {
    /// Policy used by momentum stepping, including its gas-pressure variants.
    #[must_use]
    pub const fn conversion_redistribution(&self) -> ConversionRedistribution2 {
        self.redistribution
    }

    /// Select local conversion closure without changing accepted populations.
    /// The angular policy requires nonperiodic coordinates: no moment about a
    /// torus is inferred. The same condition is checked again at every step.
    /// Native legacy step and pressure-only methods are unaffected; select
    /// conservative momentum stepping to use this policy.
    pub fn set_conversion_redistribution(
        &mut self,
        policy: ConversionRedistribution2,
    ) -> Result<(), FreeSurfaceMomentumError2> {
        if policy == ConversionRedistribution2::ConservativeLocalAngular
            && (self.grid.periodic_x || self.grid.periodic_y)
        {
            return Err(invalid("angular conversion requires a nonperiodic domain"));
        }
        self.redistribution = policy;
        Ok(())
    }

    /// Advance the original collision/streaming owner, then carry each signed
    /// conversion mass with its donor's physical velocity. Only adjacent
    /// interface cells with actual mass/capacity participate; disconnected
    /// droplets cannot exchange mass or momentum through a global pool.
    ///
    /// One staged surface makes cancellation, missing capacity and population
    /// refusal atomic. The legacy `step` retains its historical replay behavior.
    /// The receipt covers conversion alone: it is not a whole-step conservation
    /// certificate or physical dissipation. The selected angular policy closes
    /// both momenta locally, but cell-centre relocation does not preserve first
    /// mass moments or define a Galilean-invariant geometric remap.
    pub fn step_with_momentum(
        &mut self,
        gate: &CancelGate,
    ) -> Result<ConversionMomentumReport2, FreeSurfaceMomentumError2> {
        let (candidate, report) = self.momentum_candidate(gate)?;
        poll(gate)?;
        *self = candidate;
        Ok(report)
    }

    /// Advance with local signed mass, linear and cell-centre angular momentum
    /// preserved during conversion. Only cells connected by actual mass routes
    /// share the minimum mass-weighted velocity correction; detached fronts
    /// cannot supply one another's corrective impulse or torque.
    ///
    /// Positive final masses supply the centroid and planar inertia. An
    /// unresolved spin, inadmissible speed or nonpositive corrected population
    /// refuses the whole candidate. Density and physical viscous stress retain
    /// the same Guo-preserving retarget as linear conversion.
    ///
    /// This call forces angular closure for one step without changing the
    /// stored conversion policy. Pressure/momentum stepping uses that stored
    /// policy, so select ConservativeLocalAngular to compose pressure and spin.
    /// Periodic domains refuse; relocated first mass moments do not define a
    /// Galilean-equivariant remap. The receipt excludes full-step torque/heat.
    pub fn step_with_angular_momentum(
        &mut self,
        gate: &CancelGate,
    ) -> Result<ConversionMomentumReport2, FreeSurfaceMomentumError2> {
        let (candidate, report) = self.angular_momentum_candidate(gate)?;
        poll(gate)?;
        *self = candidate;
        Ok(report)
    }

    // Material units and contact-line histories must admit the SAME native
    // proposal before publication. Return its one private copy rather than
    // nesting a second transactional clone in each physical adapter.
    pub(crate) fn momentum_candidate(
        &self,
        gate: &CancelGate,
    ) -> Result<(Self, ConversionMomentumReport2), FreeSurfaceMomentumError2> {
        self.momentum_candidate_with_stream(gate, |candidate| candidate.stream_step(Some(gate)))
    }

    pub(crate) fn angular_momentum_candidate(
        &self,
        gate: &CancelGate,
    ) -> Result<(Self, ConversionMomentumReport2), FreeSurfaceMomentumError2> {
        poll(gate)?;
        if self.grid.periodic_x || self.grid.periodic_y {
            return Err(invalid(
                "angular conversion requires nonperiodic lattice coordinates",
            ));
        }
        // Force the one-call policy only in the shared owner's single private
        // clone. Its admission, pressure-free stream, conversion and readback
        // are unchanged; accepting this entry does not reconfigure later steps.
        let (mut candidate, report) = self.momentum_candidate_with_stream(gate, |candidate| {
            candidate.redistribution = ConversionRedistribution2::ConservativeLocalAngular;
            candidate.stream_step(Some(gate))
        })?;
        candidate.redistribution = self.redistribution;
        Ok((candidate, report))
    }

    // A prescribed gas load changes only the original streaming boundary.
    // Both entry points must retain this same admission/conversion sequence
    // and one private copy; the caller's error type preserves pressure errors.
    pub(super) fn momentum_candidate_with_stream<Error>(
        &self,
        gate: &CancelGate,
        stream: impl FnOnce(&mut Self) -> Result<(), Error>,
    ) -> Result<(Self, ConversionMomentumReport2), Error>
    where
        Error: From<FreeSurfaceMomentumError2>,
    {
        self.admit_momentum_step(gate)?;
        let mut candidate = self.clone();
        poll(gate)?;
        stream(&mut candidate)?;
        candidate.admit_momentum_step(gate)?;
        let report =
            if candidate.redistribution == ConversionRedistribution2::ConservativeLocalAngular {
                candidate.convert_angular_momentum(gate)?
            } else {
                candidate.convert_momentum(gate)?
            };
        candidate.admit_momentum_step(gate)?;
        poll(gate)?;
        Ok((candidate, report))
    }

    fn admit_momentum_step(&self, gate: &CancelGate) -> Result<(), FreeSurfaceMomentumError2> {
        poll(gate)?;
        if self.redistribution == ConversionRedistribution2::ConservativeLocalAngular
            && (self.grid.periodic_x || self.grid.periodic_y)
        {
            return Err(invalid("angular conversion requires a nonperiodic domain"));
        }
        let n = self
            .grid
            .nx
            .checked_mul(self.grid.ny)
            .filter(|&n| n > 0)
            .ok_or_else(|| invalid("invalid grid extent"))?;
        if self.grid.f.len() != n
            || self.grid.flags.len() != n
            || self.grid.tau.len() != n
            || self.grid.fext.len() != n
            || self.mass.len() != n
            || self.fill_raw.len() != n
            || self.fill_smooth.len() != n
            || !self.grid.g.iter().all(|g| g.is_finite())
            || self.carry != 0.0
            || !crate::admit_explicit_capillary_coefficient(self.sigma)
        {
            return Err(invalid(
                "complete bounded grid with zero unlocated mass carry required",
            ));
        }
        if let crate::core2::CollisionModel2::Trt { magic } = self.grid.collision {
            if !magic.is_finite() || magic <= 0.0 {
                return Err(invalid("finite positive TRT magic parameter required"));
            }
        }
        let wall = self
            .planar_contact
            .map(|c| c.wall)
            .or_else(|| self.planar_contact_field.as_ref().map(|c| c.wall));
        if let Some(wall) = wall {
            if self.contact != ContactModel::Neutral || self.admit_planar_wall(wall).is_err() {
                return Err(invalid("planar contact wall or neutral fallback changed"));
            }
        }
        for i in 0..n {
            if i % 256 == 0 {
                poll(gate)?;
            }
            // Dry neighbours may become wet during the conversion cascade.
            // Their held material force and relaxation must be usable before
            // another new interface samples them as an initialization donor.
            if self.grid.flags[i] != Cell::Wall
                && (!self.grid.tau[i].is_finite()
                    || self.grid.tau[i] <= 0.5
                    || self.grid.fext[i].iter().any(|f| !f.is_finite()))
            {
                return Err(invalid(
                    "finite force and relaxation required in potential wet cells",
                ));
            }
            if !matches!(self.grid.flags[i], Cell::Fluid | Cell::Interface) {
                continue;
            }
            if !self.mass[i].is_finite()
                || !self.grid.tau[i].is_finite()
                || self.grid.tau[i] <= 0.5
                || self.grid.f[i].iter().any(|f| !f.is_finite() || *f <= 0.0)
                || self.grid.fext[i].iter().any(|f| !f.is_finite())
            {
                return Err(invalid(
                    "finite positive wet-cell populations and relaxation required",
                ));
            }
            let rho = self.grid.f[i].iter().sum::<f64>();
            if !rho.is_finite() || rho <= 0.0 {
                return Err(invalid("finite positive wet-cell density required"));
            }
            if !crate::admit_explicit_capillary_density(self.sigma, rho) {
                return Err(invalid(
                    "surface tension exceeds conservative explicit capillary-wave domain at current density",
                ));
            }
            admit_velocity(self.grid.moments(i).u)?;
        }
        Ok(())
    }

    fn convert_angular_momentum(
        &mut self,
        gate: &CancelGate,
    ) -> Result<ConversionMomentumReport2, FreeSurfaceMomentumError2> {
        let policy = self.redistribution;
        self.set_conversion_redistribution(ConversionRedistribution2::ConservativeLocalAngular)?;
        let result = self.convert_momentum(gate);
        self.redistribution = policy;
        result
    }

    fn convert_momentum(
        &mut self,
        gate: &CancelGate,
    ) -> Result<ConversionMomentumReport2, FreeSurfaceMomentumError2> {
        poll(gate)?;
        let mut angular = (self.redistribution
            == ConversionRedistribution2::ConservativeLocalAngular)
            .then(angular::AngularConversion2::default);
        let excess = self.conversion_cascade(true, |surface, cell| {
            if let Some(angular) = &mut angular {
                // Capture before class/mass changes, including a source which
                // survives the cascade and later receives or repays mass.
                angular.capture(cell, surface.mass[cell], surface.grid.moments(cell).u);
            }
        })?;
        let mut recipients = BTreeMap::<usize, Recipient>::new();
        let mut report = ConversionMomentumReport2::default();
        for (source, amount) in excess {
            poll(gate)?;
            if amount == 0.0 {
                continue;
            }
            // The cascade changes flags/mass, never the converted donor's
            // populations or its held force. Even a removed negative-mass
            // cell therefore retains the physical velocity of its signed debt.
            let donor = self.grid.moments(source);
            admit_velocity(donor.u)?;
            let mut neighbors = [(0, 0.0, 0.0); Q - 1];
            let mut count = 0;
            let mut capacity = 0.0;
            let mut mass_scale = donor.rho + amount.abs();
            for q in 1..Q {
                let Some(nb) =
                    neighbor(&self.grid, source % self.grid.nx, source / self.grid.nx, q)
                else {
                    continue;
                };
                if self.grid.flags[nb] != Cell::Interface
                    || neighbors[..count].iter().any(|&(i, _, _)| i == nb)
                {
                    continue;
                }
                let rho = self.grid.f[nb].iter().sum::<f64>();
                let mass = self.mass[nb];
                if !rho.is_finite() || rho <= 0.0 || !mass.is_finite() || mass < 0.0 || mass > rho {
                    return Err(invalid("bounded adjacent inventories required; refine"));
                }
                let available = if amount > 0.0 { rho - mass } else { mass };
                if available > 0.0 {
                    neighbors[count] = (nb, rho, available);
                    count += 1;
                    capacity += available;
                    mass_scale += rho;
                }
            }
            let mass_roundoff = 128.0 * f64::EPSILON * mass_scale;
            if count == 0
                || !capacity.is_finite()
                || !mass_roundoff.is_finite()
                || amount.abs() > capacity + mass_roundoff
            {
                return Err(invalid(
                    "insufficient adjacent mass/capacity; refine or remap",
                ));
            }
            report.kinetic_energy_change -= 0.5 * amount * norm2(donor.u);
            let mut remaining = amount;
            let mut remaining_capacity = capacity;
            for (j, &(nb, rho, available)) in neighbors[..count].iter().enumerate() {
                let magnitude = if amount > 0.0 {
                    remaining.max(0.0)
                } else {
                    (-remaining).max(0.0)
                };
                let wanted = if j + 1 == count {
                    magnitude
                } else {
                    magnitude * (available / remaining_capacity.max(available))
                };
                let selected = wanted.min(available);
                let old_mass = self.mass[nb];
                let next_mass = if selected == available {
                    if amount > 0.0 { rho } else { 0.0 }
                } else {
                    old_mass + selected.copysign(amount)
                };
                let actual = next_mass - old_mass;
                if !next_mass.is_finite()
                    || next_mass < 0.0
                    || next_mass > rho
                    || actual * amount < 0.0
                {
                    return Err(invalid("conversion mass is unrepresentable"));
                }
                if actual != 0.0 {
                    let entry = recipients.entry(nb).or_insert_with(|| Recipient {
                        mass: old_mass,
                        velocity: self.grid.moments(nb).u,
                        relative_momentum: [0.0; 2],
                        carried_momentum: [0.0; 2],
                        transferred_mass: 0.0,
                    });
                    if let Some(angular) = &mut angular {
                        angular.link(source, nb, old_mass, entry.velocity, actual);
                    }
                    for d in 0..2 {
                        entry.relative_momentum[d] += actual * (donor.u[d] - entry.velocity[d]);
                        entry.carried_momentum[d] += actual * donor.u[d];
                    }
                    entry.transferred_mass += actual.abs();
                    self.mass[nb] = next_mass;
                    report.redistributed_mass += actual.abs();
                    report.routes += 1;
                }
                remaining -= actual;
                remaining_capacity -= available;
            }
            if !remaining.is_finite() || remaining.abs() > mass_roundoff {
                return Err(invalid("local conversion mass defect exceeds roundoff"));
            }
            report.absolute_mass_defect += remaining.abs();
            for d in 0..2 {
                report.absolute_momentum_defect[d] += (remaining * donor.u[d]).abs();
            }
        }
        // Round each recipient's populations only once after all simultaneous
        // transfers. Relative velocities retain an exactly co-moving state.
        for (i, recipient) in recipients {
            poll(gate)?;
            let mass = self.mass[i];
            if mass > 0.0 {
                let next = core::array::from_fn(|d| {
                    recipient.velocity[d] + recipient.relative_momentum[d] / mass
                });
                admit_velocity(next)?;
                let rho = self.grid.f[i].iter().sum::<f64>();
                let force =
                    core::array::from_fn(|d| self.grid.g[d].mul_add(rho, self.grid.fext[i][d]));
                self.grid.f[i] =
                    retarget_velocity(self.grid.f[i], recipient.velocity, next, force)?;
            }
            let actual = self.grid.moments(i).u;
            let dm = mass - recipient.mass;
            let du: [f64; 2] = core::array::from_fn(|d| actual[d] - recipient.velocity[d]);
            let roundoff =
                256.0 * f64::EPSILON * (mass + recipient.mass + recipient.transferred_mass);
            for d in 0..2 {
                let defect = (dm * recipient.velocity[d] + mass * du[d]
                    - recipient.carried_momentum[d])
                    .abs();
                if !defect.is_finite() || defect > roundoff {
                    return Err(invalid(
                        "recipient momentum cannot fit positive bounded cells; refine",
                    ));
                }
                report.absolute_momentum_defect[d] += defect;
            }
            report.kinetic_energy_change += 0.5 * dm * norm2(recipient.velocity)
                + mass
                    * (0..2)
                        .map(|d| (recipient.velocity[d] + 0.5 * du[d]) * du[d])
                        .sum::<f64>();
        }
        if let Some(angular) = angular {
            angular.finish(self, &mut report, gate)?;
        }
        if [
            report.redistributed_mass,
            report.absolute_mass_defect,
            report.kinetic_energy_change,
        ]
        .iter()
        .chain(&report.absolute_momentum_defect)
        .chain(report.absolute_angular_momentum_defect.iter())
        .any(|x| !x.is_finite())
        {
            return Err(invalid("conversion balance overflowed"));
        }
        Ok(report)
    }
}

struct Recipient {
    mass: f64,
    velocity: [f64; 2],
    relative_momentum: [f64; 2],
    carried_momentum: [f64; 2],
    transferred_mass: f64,
}
fn norm2(u: [f64; 2]) -> f64 {
    u[0] * u[0] + u[1] * u[1]
}

// Retain density and physical nonequilibrium stress, using the original Guo
// source-frame owner. Equilibrium replacement would erase viscous stress.
fn retarget_velocity(
    mut f: [f64; Q],
    previous: [f64; 2],
    next: [f64; 2],
    force: [f64; 2],
) -> Result<[f64; Q], FreeSurfaceMomentumError2> {
    if previous == next {
        return Ok(f);
    }
    let rho = f.iter().sum::<f64>();
    crate::thermal::forcing::reframe2(&mut f, previous, force, [0.0; 2]);
    let old_eq = equilibrium(rho, previous[0], previous[1]);
    let new_eq = equilibrium(rho, next[0], next[1]);
    for q in 0..Q {
        f[q] += new_eq[q] - old_eq[q];
    }
    crate::thermal::forcing::reframe2(&mut f, next, [0.0; 2], force);
    // Match the native ordered density reduction exactly: its bits also decide
    // whether an interface is full. Rest corrections carry no momentum/stress.
    f[0] += rho - f.iter().sum::<f64>();
    let rest = f[0];
    let below = f.iter().sum::<f64>() < rho;
    for _ in 0..8 {
        if f.iter().sum::<f64>() == rho {
            break;
        }
        f[0] = if below {
            f[0].next_up()
        } else {
            f[0].next_down()
        };
    }
    if f.iter().sum::<f64>() != rho {
        f[0] = rest;
        let last = rho - f[..Q - 1].iter().sum::<f64>();
        if (last - f[Q - 1]).abs() > 64.0 * f64::EPSILON * rho {
            return Err(invalid("population density closure exceeds roundoff"));
        }
        f[Q - 1] = last;
    }
    if f.iter().sum::<f64>() != rho || f.iter().any(|x| !x.is_finite() || *x <= 0.0) {
        return Err(invalid(
            "momentum retarget needs representable positive populations",
        ));
    }
    Ok(f)
}

#[cfg(test)]
mod tests;
