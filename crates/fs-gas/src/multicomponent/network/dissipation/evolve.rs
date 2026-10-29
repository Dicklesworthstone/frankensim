use super::{
    DissipationControl, DissipationError, DissipationReport, MAX_COMPONENTS, MixtureNetwork1d,
    NetworkDissipationOperator, NetworkDissipationReport, NetworkTotals, PortExchange, Time,
    Totals, TransportLaw, poll, positive,
};
use crate::multicomponent::dissipation::evolve as molecular;
use crate::multicomponent::evolve::totals;
use crate::multicomponent::{MAX_SPECIES, sum};

#[derive(Clone, Copy)]
enum Input {
    Current,
    Stage,
}

pub(super) fn advance<L: TransportLaw>(
    owner: &mut NetworkDissipationOperator<L>,
    network: &mut MixtureNetwork1d,
    duration: Time,
    controls: DissipationControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<NetworkDissipationReport, DissipationError> {
    poll(cancelled)?;
    controls.validate()?;
    if !duration.value().is_finite() || duration.value() < 0.0 {
        return Err(DissipationError::Invalid(
            "nonnegative finite molecular duration required",
        ));
    }
    if owner.gas != *network.gas()
        || owner.law.source() != owner.gas.source()
        || owner.connections != network.connections
        || owner.endpoint_connections != network.endpoint_connections
        || owner.branches.len() != network.branches.len()
        || owner.plenum_volumes.len() != network.plenums.len()
    {
        return Err(DissipationError::Invalid(
            "molecular operator and network ownership differ",
        ));
    }
    let n = owner.gas.species_count();
    let m = n + 2;
    for (index, (branch, shape)) in network.branches.iter().zip(&owner.branches).enumerate() {
        poll(cancelled)?;
        if branch.gas != owner.gas
            || branch.mesh != shape.mesh
            || branch.boundaries != shape.boundaries
        {
            return Err(DissipationError::Invalid(
                "molecular branch geometry or boundary changed",
            ));
        }
        let destination = &mut owner.work.current
            [shape.cell_start * m..(shape.cell_start + branch.cell_count()) * m];
        for (source, target) in branch.states.chunks(1024).zip(destination.chunks_mut(1024)) {
            poll(cancelled)?;
            target.copy_from_slice(source);
        }
        owner.work.branch_before[index] =
            totals(&branch.states, shape.mesh.volumes(), n, cancelled)?;
    }
    for (index, plenum) in network.plenums.iter().enumerate() {
        poll(cancelled)?;
        if plenum.volume_m3.to_bits() != owner.plenum_volumes[index].to_bits() {
            return Err(DissipationError::Invalid(
                "molecular plenum geometry changed",
            ));
        }
        let start = (owner.plenum_cell_start + index) * m;
        owner.work.current[start..start + n].copy_from_slice(&plenum.state.partial_densities);
        owner.work.current[start + n] = 0.0;
        owner.work.current[start + n + 1] = plenum.state.energy;
    }
    let before = network_totals(owner, &owner.work.current, cancelled)?;
    owner.work.endpoint_changes.fill(0.0);
    owner.work.endpoint_corrections.fill(0.0);
    let mut elapsed = 0.0;
    let mut accepted_steps = 0;
    let mut rejected_steps = 0;
    let mut maximum_parabolic = 0.0_f64;
    let mut maximum_gradient = 0.0_f64;
    let mut maximum_mass_frame = 0.0_f64;
    while elapsed < duration.value() {
        poll(cancelled)?;
        let remaining = duration.value() - elapsed;
        let mut step = remaining.min(controls.maximum_step_s);
        let mut accepted = false;
        let mut last = DissipationError::Invalid("no admissible molecular network step");
        for _ in 0..=controls.maximum_halvings {
            poll(cancelled)?;
            if accepted_steps + rejected_steps >= controls.maximum_steps {
                return Err(DissipationError::WorkLimit {
                    elapsed_time_s: elapsed,
                    accepted_steps,
                    rejected_steps,
                });
            }
            let a = spatial(owner, Input::Current, controls, cancelled)?;
            if a.rate > 0.0 {
                step = step.min(controls.parabolic_cfl / a.rate);
            }
            if !positive(step)
                || elapsed + step == elapsed
                || (step < controls.minimum_step_s && remaining > controls.minimum_step_s)
            {
                return Err(DissipationError::StepFailure {
                    elapsed_time_s: elapsed,
                    step_s: step,
                    last: Box::new(last),
                });
            }
            let attempt = (|| {
                forward(owner, Input::Current, step, cancelled)?;
                let b = spatial(owner, Input::Stage, controls, cancelled)?;
                let number = step * a.rate.max(b.rate);
                if !number.is_finite()
                    || number > controls.parabolic_cfl * (1.0 + 8.0 * f64::EPSILON)
                {
                    return Err(DissipationError::ParabolicNumber {
                        observed: number,
                        maximum: controls.parabolic_cfl,
                    });
                }
                forward(owner, Input::Stage, step, cancelled)?;
                for (index, q) in owner.work.candidate.chunks_exact_mut(m).enumerate() {
                    poll(cancelled)?;
                    for k in 0..m {
                        q[k] = 0.5 * owner.work.current[index * m + k] + 0.5 * q[k];
                    }
                    owner.work.primitive[index] = owner.gas.decode(
                        q,
                        Some(owner.work.primitive[index].temperature),
                        cancelled,
                    )?;
                }
                Ok::<_, DissipationError>(b)
            })();
            match attempt {
                Ok(b) => {
                    for (index, shape) in owner.branches.iter().enumerate() {
                        poll(cancelled)?;
                        let area = shape.mesh.face_areas()[0];
                        for end in 0..2 {
                            let face = shape.face_start
                                + if end == 0 { 0 } else { shape.mesh.cell_count() };
                            let sign = if end == 0 { 1.0 } else { -1.0 };
                            for k in 0..m {
                                let target = (2 * index + end) * m + k;
                                let change = 0.5
                                    * step
                                    * area
                                    * sign
                                    * (owner.work.endpoints_a[target]
                                        + owner.work.faces[face * m + k]);
                                compensated_add(
                                    &mut owner.work.endpoint_changes[target],
                                    &mut owner.work.endpoint_corrections[target],
                                    change,
                                )?;
                            }
                        }
                    }
                    maximum_parabolic = maximum_parabolic.max(step * a.rate.max(b.rate));
                    maximum_gradient = maximum_gradient
                        .max(a.gradient_residual)
                        .max(b.gradient_residual);
                    maximum_mass_frame = maximum_mass_frame
                        .max(a.mass_frame_residual)
                        .max(b.mass_frame_residual);
                    core::mem::swap(&mut owner.work.current, &mut owner.work.candidate);
                    elapsed = if step == remaining {
                        duration.value()
                    } else {
                        elapsed + step
                    };
                    accepted_steps += 1;
                    accepted = true;
                    break;
                }
                Err(DissipationError::Cancelled) => return Err(DissipationError::Cancelled),
                Err(error) => {
                    last = error;
                    rejected_steps += 1;
                    step *= 0.5;
                }
            }
        }
        if !accepted {
            return Err(DissipationError::StepFailure {
                elapsed_time_s: elapsed,
                step_s: step,
                last: Box::new(last),
            });
        }
    }
    let after = network_totals(owner, &owner.work.current, cancelled)?;
    let mut external = super::super::Accumulator::new(n);
    for (index, shape) in owner.branches.iter().enumerate() {
        poll(cancelled)?;
        let branch_before = owner.work.branch_before[index];
        let branch_after = totals(
            &owner.work.current
                [shape.cell_start * m..(shape.cell_start + shape.mesh.cell_count()) * m],
            shape.mesh.volumes(),
            n,
            cancelled,
        )?;
        let mut boundary = [0.0; MAX_COMPONENTS];
        for end in 0..2 {
            let change =
                &owner.work.endpoint_changes[(2 * index + end) * m..(2 * index + end + 1) * m];
            for k in 0..m {
                boundary[k] += change[k];
            }
            if owner.endpoint_connections[index][end].is_none() {
                external.add(&change[..n], change[n + 1])?;
            }
        }
        let mut residual = [0.0; MAX_COMPONENTS];
        for k in 0..m {
            residual[k] = branch_after.get(k) - branch_before.get(k) - boundary[k];
        }
        let absolute_energy = sum(core::iter::once(residual[n + 1])
            .chain((0..n).map(|k| residual[k] * owner.gas.energy_references()[k])));
        if !boundary[..m]
            .iter()
            .chain(&residual[..m])
            .all(|x| x.is_finite())
            || !absolute_energy.is_finite()
        {
            return Err(DissipationError::Invalid(
                "non-finite molecular branch balance",
            ));
        }
        owner.work.branch_report_candidate[index] = DissipationReport {
            elapsed: duration,
            accepted_steps,
            rejected_steps,
            maximum_parabolic_number: maximum_parabolic,
            maximum_gradient_residual_per_m: maximum_gradient,
            maximum_mass_frame_residual_kg_m2_s: maximum_mass_frame,
            before: branch_before,
            after: branch_after,
            boundary_change: Totals::from_values(&boundary[..m], n),
            conservation_residual: Totals::from_values(&residual[..m], n),
            absolute_energy_residual_j: absolute_energy,
        };
    }
    for (index, connection) in owner.connections.iter().enumerate() {
        poll(cancelled)?;
        let start = (2 * connection.branch + connection.end.index()) * m;
        let change = &owner.work.endpoint_changes[start..start + m];
        let mut species = [0.0; MAX_SPECIES];
        for k in 0..n {
            species[k] = -change[k];
        }
        owner.work.port_exchange_candidate[index] = PortExchange {
            species_masses_to_plenum: species,
            species_count: n,
            momentum_to_plenum: -change[n],
            energy_to_plenum: -change[n + 1],
        };
    }
    let boundary_change = external.totals();
    let mut residual = NetworkTotals {
        species_masses: [0.0; MAX_SPECIES],
        species_count: n,
        energy: after.energy - before.energy - boundary_change.energy,
    };
    for k in 0..n {
        residual.species_masses[k] =
            after.species_masses[k] - before.species_masses[k] - boundary_change.species_masses[k];
    }
    let absolute_energy_residual_j = residual.absolute_energy(&owner.gas)?;
    if !residual.species_masses[..n].iter().all(|x| x.is_finite()) {
        return Err(DissipationError::Invalid(
            "non-finite molecular network balance",
        ));
    }
    let report = NetworkDissipationReport {
        elapsed: duration,
        accepted_steps,
        rejected_steps,
        maximum_parabolic_number: maximum_parabolic,
        maximum_gradient_residual_per_m: maximum_gradient,
        maximum_mass_frame_residual_kg_m2_s: maximum_mass_frame,
        before,
        after,
        boundary_change,
        conservation_residual: residual,
        absolute_energy_residual_j,
    };
    poll(cancelled)?;
    // All fallible work and cancellation checks precede one complete publication.
    for (branch, shape) in network.branches.iter_mut().zip(&owner.branches) {
        let start = shape.cell_start * m;
        let end = start + branch.cell_count() * m;
        branch
            .states
            .copy_from_slice(&owner.work.current[start..end]);
    }
    for (index, plenum) in network.plenums.iter_mut().enumerate() {
        let start = (owner.plenum_cell_start + index) * m;
        plenum
            .state
            .partial_densities
            .copy_from_slice(&owner.work.current[start..start + n]);
        plenum.state.momentum = 0.0;
        plenum.state.energy = owner.work.current[start + n + 1];
        if duration.value() > 0.0 {
            network.work.plenum_primitive[index] =
                owner.work.primitive[owner.plenum_cell_start + index];
        }
    }
    if duration.value() > 0.0 {
        network.last_report = None;
    }
    Ok(report)
}

fn spatial<L: TransportLaw>(
    owner: &mut NetworkDissipationOperator<L>,
    input: Input,
    controls: DissipationControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<molecular::Stage, DissipationError> {
    let NetworkDissipationOperator {
        gas,
        branches,
        connections,
        plenum_volumes,
        plenum_lengths,
        plenum_cell_start,
        law,
        maxwell_stefan,
        mass_ratio,
        work,
        ..
    } = owner;
    let n = gas.species_count();
    let m = n + 2;
    let states = match input {
        Input::Current => &work.current,
        Input::Stage => &work.stage,
    };
    let mut aggregate = molecular::Stage {
        left: [0.0; MAX_COMPONENTS],
        right: [0.0; MAX_COMPONENTS],
        rate: 0.0,
        gradient_residual: 0.0,
        mass_frame_residual: 0.0,
    };
    for shape in branches.iter() {
        poll(cancelled)?;
        let count = shape.mesh.cell_count();
        let cs = shape.cell_start;
        let fs = shape.face_start;
        let stage = molecular::spatial(
            gas,
            &shape.mesh,
            shape.molecular_boundaries,
            &shape.widths,
            *mass_ratio,
            law,
            maxwell_stefan,
            &states[cs * m..(cs + count) * m],
            &mut work.primitive[cs..cs + count],
            &mut work.heat_capacity_v[cs..cs + count],
            &mut work.faces[fs * m..(fs + count + 1) * m],
            &mut work.face_rates[fs..fs + count + 1],
            &mut work.binary_diffusivities,
            &mut work.diffusion,
            controls.diffusion,
            cancelled,
        )?;
        aggregate.gradient_residual = aggregate.gradient_residual.max(stage.gradient_residual);
        aggregate.mass_frame_residual =
            aggregate.mass_frame_residual.max(stage.mass_frame_residual);
    }
    work.plenum_rhs.fill(0.0);
    work.plenum_corrections.fill(0.0);
    work.plenum_rates.fill(0.0);
    for index in 0..plenum_volumes.len() {
        poll(cancelled)?;
        let cell = *plenum_cell_start + index;
        let q = &states[cell * m..(cell + 1) * m];
        let p = gas.decode(q, Some(work.primitive[cell].temperature), cancelled)?;
        work.primitive[cell] = p;
        work.heat_capacity_v[cell] = molecular::heat_capacity_density(gas, q, p)?;
    }
    for (index, connection) in connections.iter().enumerate() {
        poll(cancelled)?;
        let shape = &branches[connection.branch];
        let left_end = connection.end.index() == 0;
        let local = if left_end {
            0
        } else {
            shape.mesh.cell_count() - 1
        };
        let branch_cell = shape.cell_start + local;
        let plenum_cell = *plenum_cell_start + connection.plenum;
        let face = shape.face_start + if left_end { 0 } else { shape.mesh.cell_count() };
        let (left, right, wl, wr) = if left_end {
            (
                plenum_cell,
                branch_cell,
                2.0 * plenum_lengths[index].value(),
                shape.widths[local],
            )
        } else {
            (
                branch_cell,
                plenum_cell,
                shape.widths[local],
                2.0 * plenum_lengths[index].value(),
            )
        };
        if !positive(wl) || !positive(wr) {
            return Err(DissipationError::Invalid(
                "molecular port distance overflow",
            ));
        }
        let rate = molecular::interior_face(
            gas,
            law,
            maxwell_stefan,
            &mut work.diffusion,
            &states[left * m..(left + 1) * m],
            &states[right * m..(right + 1) * m],
            work.primitive[left],
            work.primitive[right],
            work.heat_capacity_v[left],
            work.heat_capacity_v[right],
            wl,
            wr,
            *mass_ratio,
            &mut work.binary_diffusivities,
            &mut work.faces[face * m..(face + 1) * m],
            controls.diffusion,
            &mut aggregate,
            cancelled,
        )?;
        work.face_rates[face] = rate;
        let area = shape.mesh.face_areas()[0];
        let sign = if left_end { -1.0 } else { 1.0 };
        for k in 0..m {
            if k == n {
                continue;
            }
            let target = connection.plenum * m + k;
            compensated_add(
                &mut work.plenum_rhs[target],
                &mut work.plenum_corrections[target],
                sign * area * work.faces[face * m + k],
            )?;
        }
        work.plenum_rates[connection.plenum] += (area / plenum_volumes[connection.plenum]) * rate;
    }
    for (index, shape) in branches.iter().enumerate() {
        poll(cancelled)?;
        let area = shape.mesh.face_areas()[0];
        for i in 0..shape.mesh.cell_count() {
            let rate = (area / shape.mesh.volumes()[i])
                * (work.face_rates[shape.face_start + i]
                    + work.face_rates[shape.face_start + i + 1]);
            if !rate.is_finite() || rate < 0.0 {
                return Err(DissipationError::Invalid(
                    "molecular branch parabolic scale overflow",
                ));
            }
            aggregate.rate = aggregate.rate.max(rate);
        }
        if matches!(input, Input::Current) {
            for end in 0..2 {
                let face = shape.face_start + if end == 0 { 0 } else { shape.mesh.cell_count() };
                work.endpoints_a[(2 * index + end) * m..(2 * index + end + 1) * m]
                    .copy_from_slice(&work.faces[face * m..(face + 1) * m]);
            }
        }
    }
    for &rate in &work.plenum_rates {
        if !rate.is_finite() || rate < 0.0 {
            return Err(DissipationError::Invalid(
                "compact molecular parabolic scale overflow",
            ));
        }
        aggregate.rate = aggregate.rate.max(rate);
    }
    Ok(aggregate)
}

fn forward<L: TransportLaw>(
    owner: &mut NetworkDissipationOperator<L>,
    input: Input,
    dt: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DissipationError> {
    let NetworkDissipationOperator {
        gas,
        branches,
        plenum_volumes,
        plenum_cell_start,
        work,
        ..
    } = owner;
    let n = gas.species_count();
    let m = n + 2;
    let (states, out) = match input {
        Input::Current => (&work.current, &mut work.stage),
        Input::Stage => (&work.stage, &mut work.candidate),
    };
    for shape in branches.iter() {
        let count = shape.mesh.cell_count();
        let cs = shape.cell_start;
        let fs = shape.face_start;
        molecular::forward(
            gas,
            &shape.mesh,
            &states[cs * m..(cs + count) * m],
            &work.faces[fs * m..(fs + count + 1) * m],
            dt,
            &mut out[cs * m..(cs + count) * m],
            cancelled,
        )?;
    }
    for (index, &volume) in plenum_volumes.iter().enumerate() {
        poll(cancelled)?;
        let cell = *plenum_cell_start + index;
        let destination = &mut out[cell * m..(cell + 1) * m];
        for k in 0..m {
            destination[k] = states[cell * m + k] + (dt / volume) * work.plenum_rhs[index * m + k];
        }
        destination[n] = 0.0;
        if !gas.admissible(destination) {
            return Err(DissipationError::Invalid(
                "molecular plenum stage left species/caloric domain",
            ));
        }
    }
    Ok(())
}

fn network_totals<L: TransportLaw>(
    owner: &NetworkDissipationOperator<L>,
    states: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<NetworkTotals, DissipationError> {
    let n = owner.gas.species_count();
    let m = n + 2;
    let mut total = super::super::Accumulator::new(n);
    for shape in &owner.branches {
        let value = totals(
            &states[shape.cell_start * m..(shape.cell_start + shape.mesh.cell_count()) * m],
            shape.mesh.volumes(),
            n,
            cancelled,
        )?;
        total.add(&value.species_masses[..n], value.energy)?;
    }
    for (index, &volume) in owner.plenum_volumes.iter().enumerate() {
        poll(cancelled)?;
        let start = (owner.plenum_cell_start + index) * m;
        let mut species = [0.0; MAX_SPECIES];
        for k in 0..n {
            species[k] = volume * states[start + k];
        }
        total.add(&species[..n], volume * states[start + n + 1])?;
    }
    Ok(total.totals())
}

fn compensated_add(
    sum: &mut f64,
    correction: &mut f64,
    value: f64,
) -> Result<(), DissipationError> {
    let y = value - *correction;
    let next = *sum + y;
    if !value.is_finite() || !next.is_finite() {
        return Err(DissipationError::Invalid("molecular exchange overflow"));
    }
    *correction = (next - *sum) - y;
    *sum = next;
    Ok(())
}
