//! Cell momentum on the already accepted material mass transfers.

use super::{
    AcceptedMassTransport, CancelGate, ChtError, InletLocation, MaterialAdvectionStep, Route, Sum,
    bicgstab_ilu0, finite, finite_positive, poll, poll_every, refuse, residual,
};

/// Original extensive cell momentum and an explicitly authored impulse split.
#[derive(Clone, Copy, Debug)]
pub struct MaterialMomentumInput<'a> {
    /// Original momentum [kg m/s] per original grid cell; zero in solids.
    /// Keep returned momenta between steps instead of reconstructing M*u.
    pub momentum_kg_m_s: &'a [[f64; 3]],
    /// Prescribed impulse [kg m/s] applied to the original inventory BEFORE
    /// donor advection. This is not inferred from the face pressure correction.
    pub impulse_kg_m_s: &'a [[f64; 3]],
    /// Velocity [m/s] of locally injected material, indexed by grid cell.
    /// Every actually injected route needs a value; extraction uses its donor.
    pub source_velocity_m_s: &'a [Option<[f64; 3]>],
    /// Incoming/backflow velocities on x-, x+, y-, y+, z-, z+ boundaries.
    pub boundary_velocity_m_s: [Option<[f64; 3]>; 6],
}

/// Native solve and independent physical-unit readback budgets.
#[derive(Clone, Copy, Debug)]
pub struct MaterialMomentumControls {
    pub max_iterations: usize,
    pub relative_tolerance: f64,
    /// Sum of absolute local/global momentum defects and impulse
    /// representation defects over all three components [kg m/s].
    pub maximum_momentum_defect_kg_m_s: f64,
    /// Absolute kinetic identity defect, including retained mass roundoff [J].
    pub maximum_energy_defect_j: f64,
}

/// Actual inventory, external transfers, and numerical kinetic-energy loss.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialMomentumReport {
    pub momentum_before_kg_m_s: [f64; 3],
    pub momentum_after_kg_m_s: [f64; 3],
    /// Actual representable change from the prescribed pre-advection impulse.
    pub represented_impulse_kg_m_s: [f64; 3],
    pub carried_momentum_kg_m_s: [f64; 3],
    pub absolute_momentum_defect_kg_m_s: f64,
    pub kinetic_energy_before_j: f64,
    pub kinetic_energy_after_j: f64,
    /// Actual inventory increment, evaluated around component velocity datums.
    pub kinetic_energy_change_j: f64,
    /// Work of the represented pre-advection impulse on the original masses.
    pub impulse_work_j: f64,
    /// Signed boundary/injection/extraction kinetic-energy transfer.
    pub carried_kinetic_energy_j: f64,
    /// Nonnegative backward-Euler and donor-mixing loss. It remains numerical
    /// and never changes the accepted material enthalpy.
    pub numerical_transport_loss_j: f64,
    /// Signed contribution of retained mass-roundoff closure to the kinetic
    /// identity. This is a numerical defect, NOT dissipation or physical work.
    pub mass_roundoff_energy_j: f64,
    /// delta_K - impulse_work - kinetic_carry + numerical_transport_loss.
    pub energy_defect_j: f64,
    pub absolute_energy_defect_j: f64,
    pub iterations: [usize; 3],
}

/// Proposed cell momentum, without changing the material or staggered flow.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialMomentumStep {
    pub momentum_kg_m_s: Vec<[f64; 3]>,
    /// Derived from the actual returned momentum and the accepted cell masses.
    pub velocity_m_s: Vec<[f64; 3]>,
    pub report: MaterialMomentumReport,
}

#[derive(Clone, Copy, Default)]
struct Account {
    change: Sum,
    work: Sum,
    carry: Sum,
    loss: Sum,
    roundoff: Sum,
    mass_change: Sum,
    mass_carry: Sum,
    momentum_change: [Sum; 3],
    impulse: [Sum; 3],
    momentum_carry: [Sum; 3],
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|k| a[k] - b[k])
}

fn inlet_velocity(route: &Route, input: MaterialMomentumInput<'_>) -> Result<[f64; 3], ChtError> {
    match route.inlet_location {
        Some(InletLocation::Boundary(face)) => input.boundary_velocity_m_s[face],
        Some(InletLocation::Source(cell)) => input.source_velocity_m_s[cell],
        None => None,
    }
    .ok_or_else(|| refuse("incoming mass route lacks an authored velocity"))
}

fn admit(
    transport: &AcceptedMassTransport,
    input: MaterialMomentumInput<'_>,
    controls: MaterialMomentumControls,
    gate: &CancelGate,
) -> Result<(), ChtError> {
    let n = transport.row_of.len();
    if input.momentum_kg_m_s.len() != n
        || input.impulse_kg_m_s.len() != n
        || input.source_velocity_m_s.len() != n
        || controls.max_iterations == 0
    {
        return Err(refuse("invalid momentum inventory shape or solve budget"));
    }
    finite_positive("flow.momentum_tolerance", controls.relative_tolerance)?;
    for bound in [
        controls.maximum_momentum_defect_kg_m_s,
        controls.maximum_energy_defect_j,
    ] {
        if !bound.is_finite() || bound < 0.0 {
            return Err(refuse("invalid momentum or kinetic-energy balance budget"));
        }
    }
    for value in input.boundary_velocity_m_s.into_iter().flatten().flatten() {
        finite("flow.incoming_velocity", value)?;
    }
    for c in 0..n {
        poll_every(gate, c)?;
        for value in input.momentum_kg_m_s[c]
            .into_iter()
            .chain(input.impulse_kg_m_s[c])
            .chain(input.source_velocity_m_s[c].into_iter().flatten())
        {
            finite("flow.momentum_input", value)?;
        }
        if transport.row_of[c] == usize::MAX
            && (input.momentum_kg_m_s[c] != [0.0; 3]
                || input.impulse_kg_m_s[c] != [0.0; 3]
                || input.source_velocity_m_s[c].is_some())
        {
            return Err(refuse(
                "solid grid cells cannot carry material momentum or impulse",
            ));
        }
    }
    for (index, route) in transport.routes.iter().enumerate() {
        poll_every(gate, index)?;
        if route.from.is_none() {
            inlet_velocity(route, input)?;
        }
    }
    Ok(())
}

impl MaterialAdvectionStep {
    /// Transport actual cell momentum on this candidate's EXACT accepted mass
    /// transfers, including candidates returned by phase-volume evolution.
    /// No density solve, flux reconstruction, or material source is added.
    ///
    /// The existing mass*h operator solves three momentum components, centered
    /// in each connected component. Represented prescribed impulses act first.
    /// Paired donor transfers then form the authoritative endpoint momentum;
    /// its actual P/M is independently checked in the original donor equations.
    /// Zero transfers and zero impulses preserve the original momentum exactly.
    ///
    /// The independently checked kinetic identity includes backward-Euler
    /// temporal loss, upwind mixing loss, external kinetic carry and represented
    /// impulse work. Retained mass roundoff has its own defect account. None of
    /// these numerical losses is deposited in material enthalpy.
    ///
    /// This is a conservative cell advection/impulse substep, not a complete
    /// Navier-Stokes update. The native pressure correction and its work remain
    /// owned by the staggered face projection; they are NOT applied a second
    /// time to these cell momenta. The enclosing owner supplies the force split,
    /// admissible phase mobility and any cell/face momentum reconstruction.
    /// Changed public endpoint masses, missing incoming velocities, invalid
    /// budgets, failed independent readback or cancellation return no candidate.
    pub fn advect_momentum(
        &self,
        input: MaterialMomentumInput<'_>,
        controls: MaterialMomentumControls,
        gate: &CancelGate,
    ) -> Result<MaterialMomentumStep, ChtError> {
        poll(gate)?;
        let t = self.accepted_mass_transport(gate)?;
        admit(t, input, controls, gate)?;
        let n = t.cells.len();
        let mut p_star = vec![[0.0; 3]; n];
        let mut u_before = vec![[0.0; 3]; n];
        let mut u_star = vec![[0.0; 3]; n];
        let mut references = vec![None; n];
        let mut momentum_defect = Sum::default();
        for (row, &c) in t.cells.iter().enumerate() {
            poll_every(gate, row)?;
            for k in 0..3 {
                let p = input.momentum_kg_m_s[c][k];
                p_star[row][k] = p + input.impulse_kg_m_s[c][k];
                finite("flow.impulsed_momentum", p_star[row][k])?;
                momentum_defect.add((p_star[row][k] - p - input.impulse_kg_m_s[c][k]).abs());
                u_before[row][k] = p / t.initial_mass_kg[row];
                u_star[row][k] = p_star[row][k] / t.initial_mass_kg[row];
                finite("flow.original_velocity", u_before[row][k])?;
                finite("flow.impulsed_velocity", u_star[row][k])?;
            }
            references[t.parents[row]].get_or_insert(u_star[row]);
        }
        let reference: Vec<_> = t
            .parents
            .iter()
            .map(|&r| references[r].expect("component inventory"))
            .collect();
        let mut donor_velocity = vec![[0.0; 3]; n];
        let mut right_hand_sides: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; n]);
        let mut iterations = [0; 3];
        for k in 0..3 {
            let rhs = &mut right_hand_sides[k];
            let mut velocity = vec![0.0; n];
            for row in 0..n {
                poll_every(gate, row)?;
                velocity[row] = u_star[row][k] - reference[row][k];
                rhs[row] = (-t.initial_mass_kg[row]).mul_add(reference[row][k], p_star[row][k])
                    + t.closure[row] * reference[row][k];
            }
            for (index, route) in t.routes.iter().enumerate() {
                poll_every(gate, index)?;
                if route.from.is_none() {
                    let to = route.to.expect("incoming receiver");
                    rhs[to] += route.mass * (inlet_velocity(route, input)?[k] - reference[to][k]);
                }
            }
            for row in 0..n {
                poll_every(gate, row)?;
                rhs[row] *= 1.0 / t.diagonal[row];
                finite("flow.momentum_rhs", rhs[row])?;
                finite("flow.momentum_guess", velocity[row])?;
            }
            iterations[k] = bicgstab_ilu0(
                "material momentum transport",
                &t.operator,
                rhs,
                &mut velocity,
                controls.relative_tolerance,
                controls.max_iterations,
                gate,
            )?
            .iterations;
            for row in 0..n {
                poll_every(gate, row)?;
                donor_velocity[row][k] = velocity[row] + reference[row][k];
                finite("flow.momentum_donor_velocity", donor_velocity[row][k])?;
            }
        }
        let mut deltas = vec![[Sum::default(); 3]; n];
        for (index, route) in t.routes.iter().enumerate() {
            poll_every(gate, index)?;
            let velocity = match route.from {
                Some(from) => donor_velocity[from],
                None => inlet_velocity(route, input)?,
            };
            for k in 0..3 {
                let transfer = route.mass * velocity[k];
                finite("flow.transferred_momentum", transfer)?;
                if let Some(from) = route.from {
                    deltas[from][k].add(-transfer);
                }
                if let Some(to) = route.to {
                    deltas[to][k].add(transfer);
                }
            }
        }
        let mut momentum = vec![[0.0; 3]; t.row_of.len()];
        let mut velocity = vec![[0.0; 3]; t.row_of.len()];
        for (row, &c) in t.cells.iter().enumerate() {
            poll_every(gate, row)?;
            for k in 0..3 {
                momentum[c][k] = p_star[row][k] + deltas[row][k].value;
                velocity[c][k] = momentum[c][k] / t.endpoint_mass_kg[row];
                finite("flow.accepted_momentum", momentum[c][k])?;
                finite("flow.accepted_velocity", velocity[c][k])?;
                momentum_defect
                    .add((deltas[row][k].value - (momentum[c][k] - p_star[row][k])).abs());
            }
        }
        for k in 0..3 {
            let actual: Vec<_> = t
                .cells
                .iter()
                .enumerate()
                .map(|(r, &c)| velocity[c][k] - reference[r][k])
                .collect();
            momentum_defect.add(residual(
                &t.operator,
                &right_hand_sides[k],
                &actual,
                &t.diagonal,
                gate,
            )?);
        }
        let mut accounts = vec![Account::default(); n];
        let mut before = [Sum::default(); 3];
        let mut after = [Sum::default(); 3];
        let mut kinetic_before = Sum::default();
        let mut kinetic_after = Sum::default();
        let mut absolute_roundoff = Sum::default();
        for (row, &c) in t.cells.iter().enumerate() {
            poll_every(gate, row)?;
            let account = &mut accounts[t.parents[row]];
            let old_relative = sub(u_before[row], reference[row]);
            let new_relative = sub(velocity[c], reference[row]);
            let star_relative = sub(u_star[row], reference[row]);
            let impulse = sub(p_star[row], input.momentum_kg_m_s[c]);
            account
                .change
                .add(0.5 * t.endpoint_mass_kg[row] * dot(new_relative, new_relative));
            account
                .change
                .add(-0.5 * t.initial_mass_kg[row] * dot(old_relative, old_relative));
            account.work.add(dot(
                impulse,
                std::array::from_fn(|k| 0.5 * (old_relative[k] + star_relative[k])),
            ));
            let temporal = sub(velocity[c], u_star[row]);
            account
                .loss
                .add(0.5 * t.initial_mass_kg[row] * dot(temporal, temporal));
            let roundoff = 0.5 * t.closure[row] * dot(velocity[c], velocity[c]);
            account.roundoff.add(roundoff);
            absolute_roundoff.add(roundoff.abs());
            account
                .mass_change
                .add(t.endpoint_mass_kg[row] - t.initial_mass_kg[row]);
            kinetic_before.add(0.5 * dot(input.momentum_kg_m_s[c], u_before[row]));
            kinetic_after.add(0.5 * dot(momentum[c], velocity[c]));
            for k in 0..3 {
                before[k].add(input.momentum_kg_m_s[c][k]);
                after[k].add(momentum[c][k]);
                account.momentum_change[k].add(momentum[c][k] - input.momentum_kg_m_s[c][k]);
                account.impulse[k].add(impulse[k]);
            }
        }
        for (index, route) in t.routes.iter().enumerate() {
            poll_every(gate, index)?;
            let row = route.from.or(route.to).expect("route endpoint");
            let account = &mut accounts[t.parents[row]];
            let donor = match route.from {
                Some(from) => donor_velocity[from],
                None => inlet_velocity(route, input)?,
            };
            if let Some(to) = route.to {
                let difference = sub(donor, velocity[t.cells[to]]);
                account
                    .loss
                    .add(0.5 * route.mass * dot(difference, difference));
            }
            if route.from.is_none() || route.to.is_none() {
                let mass = if route.from.is_none() {
                    route.mass
                } else {
                    -route.mass
                };
                let relative = sub(donor, reference[row]);
                account.mass_carry.add(mass);
                account.carry.add(0.5 * mass * dot(relative, relative));
                for k in 0..3 {
                    account.momentum_carry[k].add(mass * donor[k]);
                }
            }
        }
        let mut impulse = [Sum::default(); 3];
        let mut carried_momentum = [Sum::default(); 3];
        let mut change = Sum::default();
        let mut work = Sum::default();
        let mut carry = Sum::default();
        let mut loss = Sum::default();
        let mut roundoff = Sum::default();
        let mut energy_defect = Sum::default();
        let mut absolute_energy_defect = absolute_roundoff;
        for (r, account) in accounts.iter().enumerate() {
            poll_every(gate, r)?;
            let Some(datum) = references[r] else {
                continue;
            };
            let mut component_defect =
                account.change.value - account.work.value - account.carry.value
                    + account.loss.value;
            let half_speed2 = 0.5 * dot(datum, datum);
            let mut physical_change =
                account.change.value - half_speed2 * account.mass_change.value;
            let mut physical_carry = account.carry.value - half_speed2 * account.mass_carry.value;
            let mut physical_work = account.work.value;
            for k in 0..3 {
                let defect = account.momentum_change[k].value
                    - account.impulse[k].value
                    - account.momentum_carry[k].value;
                momentum_defect.add(defect.abs());
                component_defect = datum[k].mul_add(defect, component_defect);
                physical_change =
                    datum[k].mul_add(account.momentum_change[k].value, physical_change);
                physical_carry = datum[k].mul_add(account.momentum_carry[k].value, physical_carry);
                physical_work = datum[k].mul_add(account.impulse[k].value, physical_work);
                impulse[k].add(account.impulse[k].value);
                carried_momentum[k].add(account.momentum_carry[k].value);
            }
            component_defect -=
                half_speed2 * (account.mass_change.value - account.mass_carry.value);
            energy_defect.add(component_defect);
            absolute_energy_defect.add((component_defect - account.roundoff.value).abs());
            change.add(physical_change);
            work.add(physical_work);
            carry.add(physical_carry);
            loss.add(account.loss.value);
            roundoff.add(account.roundoff.value);
        }
        for value in [
            kinetic_before.value,
            kinetic_after.value,
            change.value,
            work.value,
            carry.value,
            loss.value,
            roundoff.value,
            energy_defect.value,
            absolute_energy_defect.value,
            momentum_defect.value,
        ] {
            finite("flow.momentum_account", value)?;
        }
        for values in [before, after, impulse, carried_momentum] {
            for value in values {
                finite("flow.momentum_account", value.value)?;
            }
        }
        if momentum_defect.value > controls.maximum_momentum_defect_kg_m_s
            || absolute_energy_defect.value > controls.maximum_energy_defect_j
        {
            return Err(refuse(
                "actual momentum or kinetic-energy transport exceeds its balance budget",
            ));
        }
        poll(gate)?;
        Ok(MaterialMomentumStep {
            momentum_kg_m_s: momentum,
            velocity_m_s: velocity,
            report: MaterialMomentumReport {
                momentum_before_kg_m_s: before.map(|s| s.value),
                momentum_after_kg_m_s: after.map(|s| s.value),
                represented_impulse_kg_m_s: impulse.map(|s| s.value),
                carried_momentum_kg_m_s: carried_momentum.map(|s| s.value),
                absolute_momentum_defect_kg_m_s: momentum_defect.value,
                kinetic_energy_before_j: kinetic_before.value,
                kinetic_energy_after_j: kinetic_after.value,
                kinetic_energy_change_j: change.value,
                impulse_work_j: work.value,
                carried_kinetic_energy_j: carry.value,
                numerical_transport_loss_j: loss.value,
                mass_roundoff_energy_j: roundoff.value,
                energy_defect_j: energy_defect.value,
                absolute_energy_defect_j: absolute_energy_defect.value,
                iterations,
            },
        })
    }
}

#[cfg(test)]
mod tests;
