use super::super::evolve::totals;
use super::super::{
    Boundaries, Boundary, MAX_SPECIES, Mesh1d, MixtureEuler1d, MixtureGas, Primitive, Totals, sum,
};
use super::reservoir::{ReservoirEnds, ReservoirTrace};
use super::{
    DissipationControl, DissipationError, DissipationOperator, DissipationReport, MAX_COMPONENTS,
    TransportLaw, poll, positive,
};
use fs_qty::{Temperature, Time};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::transport::maxwell_stefan::{
    MaxwellStefanSolveV1, MaxwellStefanV1, MaxwellStefanWorkspaceV1,
};

#[derive(Clone, Copy)]
pub(in crate::multicomponent) struct Stage {
    pub(in crate::multicomponent) left: [f64; MAX_COMPONENTS],
    pub(in crate::multicomponent) right: [f64; MAX_COMPONENTS],
    pub(in crate::multicomponent) rate: f64,
    pub(in crate::multicomponent) gradient_residual: f64,
    pub(in crate::multicomponent) mass_frame_residual: f64,
}

pub(super) fn advance<L: TransportLaw>(
    owner: &mut DissipationOperator<L>,
    flow: &mut MixtureEuler1d,
    duration: Time,
    controls: DissipationControl,
    reservoirs: ReservoirEnds<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<DissipationReport, DissipationError> {
    advance_inner(owner, flow, duration, controls, reservoirs, None, cancelled)
}

pub(super) fn advance_moving<L: TransportLaw>(
    owner: &mut DissipationOperator<L>,
    flow: &mut MixtureEuler1d,
    duration: Time,
    controls: DissipationControl,
    wall_velocity: [f64; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<DissipationReport, DissipationError> {
    if wall_velocity.iter().any(|v| !v.is_finite())
        || flow.boundaries
            != (Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            })
    {
        return Err(DissipationError::Invalid(
            "molecular moving trial requires two finite impermeable wall speeds",
        ));
    }
    advance_inner(
        owner,
        flow,
        duration,
        controls,
        [None; 2],
        Some(wall_velocity),
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn advance_inner<L: TransportLaw>(
    owner: &mut DissipationOperator<L>,
    flow: &mut MixtureEuler1d,
    duration: Time,
    controls: DissipationControl,
    reservoirs: ReservoirEnds<'_>,
    wall_velocity: Option<[f64; 2]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<DissipationReport, DissipationError> {
    poll(cancelled)?;
    controls.validate()?;
    if !duration.value().is_finite() || duration.value() < 0.0 {
        return Err(DissipationError::Invalid(
            "nonnegative finite molecular duration required",
        ));
    }
    let geometry = if wall_velocity.is_some() {
        owner.motion_mesh.as_ref().ok_or(DissipationError::Invalid(
            "molecular motion was not admitted",
        ))?
    } else {
        &owner.mesh
    };
    if flow.gas != owner.gas
        || flow.mesh != *geometry
        || flow.boundaries != owner.boundaries
        || owner.law.source() != owner.gas.source()
    {
        return Err(DissipationError::Invalid(
            "molecular operator and flow ownership differ",
        ));
    }
    let n = owner.gas.species_count();
    let m = n + 2;
    let before = totals(&flow.states, flow.mesh.volumes(), n, cancelled)?;
    let mut boundary = [0.0; MAX_COMPONENTS];
    let mut boundary_correction = [0.0; MAX_COMPONENTS];
    let mut ends = [[0.0; MAX_COMPONENTS]; 2];
    let mut end_corrections = [[0.0; MAX_COMPONENTS]; 2];
    let mut accepted_steps = 0;
    let mut rejected_steps = 0;
    let mut elapsed = 0.0;
    let mut maximum_parabolic = 0.0_f64;
    let mut maximum_gradient_residual = 0.0_f64;
    let mut maximum_mass_frame_residual = 0.0_f64;
    owner.work.current.copy_from_slice(&flow.states);
    let DissipationOperator {
        gas,
        mesh,
        motion_mesh,
        boundaries,
        law,
        maxwell_stefan,
        widths,
        mass_ratio,
        work,
    } = owner;
    let mesh: &Mesh1d = if wall_velocity.is_some() {
        motion_mesh.as_ref().expect("admitted moving geometry")
    } else {
        mesh
    };
    let area = mesh.face_areas()[0];
    // Widths are scratch, not accepted geometry: even an aborted moving trial
    // must not poison the next stationary update's gradients or parabolic bound.
    for (i, width) in widths.iter_mut().enumerate() {
        poll(cancelled)?;
        *width = mesh.uniform_cell_width().map_or_else(
            || mesh.positions()[i + 1] - mesh.positions()[i],
            |w| w.value(),
        );
        let expected = area * *width;
        if !positive(expected)
            || (mesh.volumes()[i] - expected).abs()
                > 64.0 * f64::EPSILON * expected.max(mesh.volumes()[i])
        {
            return Err(DissipationError::Invalid(
                "molecular trial volume differs from area times physical width",
            ));
        }
    }
    while elapsed < duration.value() {
        poll(cancelled)?;
        let remaining = duration.value() - elapsed;
        let mut step = remaining.min(controls.maximum_step_s);
        let mut accepted = false;
        let mut last = DissipationError::Invalid("no admissible molecular step");
        for _ in 0..=controls.maximum_halvings {
            poll(cancelled)?;
            if accepted_steps + rejected_steps >= controls.maximum_steps {
                return Err(DissipationError::WorkLimit {
                    elapsed_time_s: elapsed,
                    accepted_steps,
                    rejected_steps,
                });
            }
            let a = spatial_with_wall_motion(
                gas,
                mesh,
                *boundaries,
                reservoirs,
                widths,
                *mass_ratio,
                law,
                maxwell_stefan,
                &work.current,
                &mut work.primitive,
                &mut work.heat_capacity_v,
                &mut work.faces,
                &mut work.face_rates,
                &mut work.binary_diffusivities,
                &mut work.diffusion,
                controls.diffusion,
                wall_velocity.unwrap_or([0.0; 2]),
                cancelled,
            )?;
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
                forward(
                    gas,
                    mesh,
                    &work.current,
                    &work.faces,
                    step,
                    &mut work.stage,
                    cancelled,
                )?;
                let b = spatial_with_wall_motion(
                    gas,
                    mesh,
                    *boundaries,
                    reservoirs,
                    widths,
                    *mass_ratio,
                    law,
                    maxwell_stefan,
                    &work.stage,
                    &mut work.primitive,
                    &mut work.heat_capacity_v,
                    &mut work.faces,
                    &mut work.face_rates,
                    &mut work.binary_diffusivities,
                    &mut work.diffusion,
                    controls.diffusion,
                    wall_velocity.unwrap_or([0.0; 2]),
                    cancelled,
                )?;
                let number = step * a.rate.max(b.rate);
                if !number.is_finite()
                    || number > controls.parabolic_cfl * (1.0 + 8.0 * f64::EPSILON)
                {
                    return Err(DissipationError::ParabolicNumber {
                        observed: number,
                        maximum: controls.parabolic_cfl,
                    });
                }
                forward(
                    gas,
                    mesh,
                    &work.stage,
                    &work.faces,
                    step,
                    &mut work.candidate,
                    cancelled,
                )?;
                for (i, q) in work.candidate.chunks_exact_mut(m).enumerate() {
                    poll(cancelled)?;
                    for k in 0..m {
                        q[k] = 0.5 * work.current[i * m + k] + 0.5 * q[k];
                    }
                    gas.decode(q, Some(work.primitive[i].temperature), cancelled)?;
                }
                Ok::<_, DissipationError>(b)
            })();
            match attempt {
                Ok(b) => {
                    for k in 0..m {
                        let change = 0.5
                            * step
                            * area
                            * ((a.left[k] + b.left[k]) - (a.right[k] + b.right[k]));
                        let y = change - boundary_correction[k];
                        let next = boundary[k] + y;
                        boundary_correction[k] = (next - boundary[k]) - y;
                        boundary[k] = next;
                        for (end, transfer) in [
                            0.5 * step * area * (a.left[k] + b.left[k]),
                            -0.5 * step * area * (a.right[k] + b.right[k]),
                        ]
                        .into_iter()
                        .enumerate()
                        {
                            let y = transfer - end_corrections[end][k];
                            let next = ends[end][k] + y;
                            end_corrections[end][k] = (next - ends[end][k]) - y;
                            ends[end][k] = next;
                        }
                    }
                    maximum_parabolic = maximum_parabolic.max(step * a.rate.max(b.rate));
                    maximum_gradient_residual = maximum_gradient_residual
                        .max(a.gradient_residual)
                        .max(b.gradient_residual);
                    maximum_mass_frame_residual = maximum_mass_frame_residual
                        .max(a.mass_frame_residual)
                        .max(b.mass_frame_residual);
                    core::mem::swap(&mut work.current, &mut work.candidate);
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
    let after = totals(&work.current, mesh.volumes(), n, cancelled)?;
    let mut residual = [0.0; MAX_COMPONENTS];
    for k in 0..m {
        residual[k] = after.get(k) - before.get(k) - boundary[k];
    }
    let absolute_energy_residual = sum(core::iter::once(residual[n + 1])
        .chain((0..n).map(|k| residual[k] * gas.energy_references()[k])));
    if !boundary[..m]
        .iter()
        .chain(&residual[..m])
        .chain(&ends[0][..m])
        .chain(&ends[1][..m])
        .all(|x| x.is_finite())
        || !absolute_energy_residual.is_finite()
    {
        return Err(DissipationError::Invalid(
            "non-finite molecular extensive balance",
        ));
    }
    poll(cancelled)?;
    // No scratch clock is advanced: the composing time integrator owns time.
    flow.states.copy_from_slice(&work.current);
    Ok(DissipationReport {
        elapsed: duration,
        accepted_steps,
        rejected_steps,
        maximum_parabolic_number: maximum_parabolic,
        maximum_gradient_residual_per_m: maximum_gradient_residual,
        maximum_mass_frame_residual_kg_m2_s: maximum_mass_frame_residual,
        before,
        after,
        boundary_change: Totals::from_values(&boundary[..m], n),
        left_boundary_change: Totals::from_values(&ends[0][..m], n),
        right_boundary_change: Totals::from_values(&ends[1][..m], n),
        conservation_residual: Totals::from_values(&residual[..m], n),
        absolute_energy_residual_j: absolute_energy_residual,
    })
}

pub(in crate::multicomponent) fn forward(
    gas: &MixtureGas,
    mesh: &Mesh1d,
    states: &[f64],
    faces: &[f64],
    dt: f64,
    out: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DissipationError> {
    let n = gas.species_count();
    let m = n + 2;
    let area = mesh.face_areas()[0];
    for (i, q) in states.chunks_exact(m).enumerate() {
        poll(cancelled)?;
        let scale = dt * (area / mesh.volumes()[i]);
        let destination = &mut out[i * m..(i + 1) * m];
        for k in 0..m {
            destination[k] = q[k] + scale * (faces[i * m + k] - faces[(i + 1) * m + k]);
        }
        if !gas.admissible(destination) {
            return Err(DissipationError::Invalid(
                "molecular stage left the nonnegative species/caloric domain",
            ));
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(in crate::multicomponent) fn spatial<L: TransportLaw>(
    gas: &MixtureGas,
    mesh: &Mesh1d,
    boundaries: Boundaries,
    reservoirs: ReservoirEnds<'_>,
    widths: &[f64],
    mass_ratio: f64,
    law: &mut L,
    diffusion: &MaxwellStefanV1,
    states: &[f64],
    primitive: &mut [Primitive],
    cv: &mut [f64],
    faces: &mut [f64],
    face_rates: &mut [f64],
    binary: &mut [f64],
    diffusion_workspace: &mut MaxwellStefanWorkspaceV1,
    controls: MaxwellStefanSolveV1,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Stage, DissipationError> {
    spatial_with_wall_motion(
        gas,
        mesh,
        boundaries,
        reservoirs,
        widths,
        mass_ratio,
        law,
        diffusion,
        states,
        primitive,
        cv,
        faces,
        face_rates,
        binary,
        diffusion_workspace,
        controls,
        [0.0; 2],
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn spatial_with_wall_motion<L: TransportLaw>(
    gas: &MixtureGas,
    mesh: &Mesh1d,
    boundaries: Boundaries,
    reservoirs: ReservoirEnds<'_>,
    widths: &[f64],
    mass_ratio: f64,
    law: &mut L,
    diffusion: &MaxwellStefanV1,
    states: &[f64],
    primitive: &mut [Primitive],
    cv: &mut [f64],
    faces: &mut [f64],
    face_rates: &mut [f64],
    binary: &mut [f64],
    diffusion_workspace: &mut MaxwellStefanWorkspaceV1,
    controls: MaxwellStefanSolveV1,
    wall_velocity: [f64; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Stage, DissipationError> {
    let n = gas.species_count();
    let m = n + 2;
    let count = mesh.cell_count();
    for (i, q) in states.chunks_exact(m).enumerate() {
        poll(cancelled)?;
        let p = gas.decode(q, Some(primitive[i].temperature), cancelled)?;
        primitive[i] = p;
        cv[i] = heat_capacity_density(gas, q, p)?;
    }
    faces.fill(0.0);
    face_rates.fill(0.0);
    let mut stage = Stage {
        left: [0.0; MAX_COMPONENTS],
        right: [0.0; MAX_COMPONENTS],
        rate: 0.0,
        gradient_residual: 0.0,
        mass_frame_residual: 0.0,
    };
    for face in 1..count {
        poll(cancelled)?;
        let rate = interior_face(
            gas,
            law,
            diffusion,
            diffusion_workspace,
            &states[(face - 1) * m..face * m],
            &states[face * m..(face + 1) * m],
            primitive[face - 1],
            primitive[face],
            cv[face - 1],
            cv[face],
            widths[face - 1],
            widths[face],
            mass_ratio,
            binary,
            &mut faces[face * m..(face + 1) * m],
            controls,
            &mut stage,
            cancelled,
        )?;
        face_rates[face] = rate;
    }
    match boundaries {
        Boundaries::Periodic => {
            let rate = interior_face(
                gas,
                law,
                diffusion,
                diffusion_workspace,
                &states[(count - 1) * m..count * m],
                &states[..m],
                primitive[count - 1],
                primitive[0],
                cv[count - 1],
                cv[0],
                widths[count - 1],
                widths[0],
                mass_ratio,
                binary,
                &mut faces[..m],
                controls,
                &mut stage,
                cancelled,
            )?;
            for k in 0..m {
                faces[count * m + k] = faces[k];
            }
            face_rates[0] = rate;
            face_rates[count] = rate;
        }
        Boundaries::Ends { left, right } => {
            if matches!(left, Boundary::Reservoir(_)) || matches!(right, Boundary::Reservoir(_)) {
                return Err(DissipationError::Invalid(
                    "molecular reservoir requires species-bearing state and energy references",
                ));
            }
            if let Some(reservoir) = reservoirs[0] {
                face_rates[0] = reservoir_face(
                    gas,
                    law,
                    diffusion,
                    diffusion_workspace,
                    reservoir,
                    &states[..m],
                    primitive[0],
                    cv[0],
                    0.5 * widths[0],
                    true,
                    mass_ratio,
                    binary,
                    &mut faces[..m],
                    controls,
                    &mut stage,
                    cancelled,
                )?;
            } else if left == Boundary::Reflecting {
                face_rates[0] = wall_face(
                    gas,
                    law,
                    &states[..m],
                    primitive[0],
                    0.5 * widths[0],
                    true,
                    wall_velocity[0],
                    binary,
                    &mut faces[..m],
                    cancelled,
                )?;
            }
            if let Some(reservoir) = reservoirs[1] {
                face_rates[count] = reservoir_face(
                    gas,
                    law,
                    diffusion,
                    diffusion_workspace,
                    reservoir,
                    &states[(count - 1) * m..count * m],
                    primitive[count - 1],
                    cv[count - 1],
                    0.5 * widths[count - 1],
                    false,
                    mass_ratio,
                    binary,
                    &mut faces[count * m..(count + 1) * m],
                    controls,
                    &mut stage,
                    cancelled,
                )?;
            } else if right == Boundary::Reflecting {
                face_rates[count] = wall_face(
                    gas,
                    law,
                    &states[(count - 1) * m..count * m],
                    primitive[count - 1],
                    0.5 * widths[count - 1],
                    false,
                    wall_velocity[1],
                    binary,
                    &mut faces[count * m..(count + 1) * m],
                    cancelled,
                )?;
            }
        }
    }
    let area = mesh.face_areas()[0];
    for i in 0..count {
        let rate = (area / mesh.volumes()[i]) * (face_rates[i] + face_rates[i + 1]);
        if !rate.is_finite() || rate < 0.0 {
            return Err(DissipationError::Invalid("parabolic scale overflow"));
        }
        stage.rate = stage.rate.max(rate);
    }
    stage.left[..m].copy_from_slice(&faces[..m]);
    stage.right[..m].copy_from_slice(&faces[count * m..(count + 1) * m]);
    Ok(stage)
}

#[allow(clippy::too_many_arguments)]
fn reservoir_face<L: TransportLaw>(
    gas: &MixtureGas,
    law: &mut L,
    diffusion: &MaxwellStefanV1,
    workspace: &mut MaxwellStefanWorkspaceV1,
    reservoir: ReservoirTrace<'_>,
    state: &[f64],
    primitive: Primitive,
    cv_density: f64,
    distance: f64,
    boundary_on_left: bool,
    mass_ratio: f64,
    binary: &mut [f64],
    out: &mut [f64],
    controls: MaxwellStefanSolveV1,
    stage: &mut Stage,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, DissipationError> {
    let bath = super::face::Trace {
        densities: reservoir.cell.partial_densities,
        primitive: reservoir.primitive,
        transverse_velocity: 0.0,
        cv_density: reservoir.heat_capacity_v,
    };
    let cell = super::face::Trace {
        densities: &state[..gas.species_count()],
        primitive,
        transverse_velocity: 0.0,
        cv_density,
    };
    let (left, right) = if boundary_on_left {
        (bath, cell)
    } else {
        (cell, bath)
    };
    let evaluated = super::face::reservoir(
        gas,
        law,
        diffusion,
        workspace,
        left,
        right,
        distance,
        boundary_on_left,
        mass_ratio,
        binary,
        out,
        controls,
        cancelled,
    )?;
    stage.gradient_residual = stage
        .gradient_residual
        .max(evaluated.maximum_gradient_residual_per_m);
    stage.mass_frame_residual = stage
        .mass_frame_residual
        .max(evaluated.mass_frame_residual_kg_m2_s.abs());
    let rate = evaluated.diffusivity_m2_s / distance;
    if !rate.is_finite() || rate < 0.0 {
        return Err(DissipationError::Invalid(
            "molecular reservoir scale overflow",
        ));
    }
    Ok(rate)
}

#[allow(clippy::too_many_arguments)]
pub(in crate::multicomponent) fn interior_face<L: TransportLaw>(
    gas: &MixtureGas,
    law: &mut L,
    diffusion: &MaxwellStefanV1,
    workspace: &mut MaxwellStefanWorkspaceV1,
    left: &[f64],
    right: &[f64],
    pl: Primitive,
    pr: Primitive,
    cvl: f64,
    cvr: f64,
    wl: f64,
    wr: f64,
    mass_ratio: f64,
    binary: &mut [f64],
    out: &mut [f64],
    controls: MaxwellStefanSolveV1,
    stage: &mut Stage,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, DissipationError> {
    let n = gas.species_count();
    let evaluated = super::face::interior::<false, L>(
        gas,
        law,
        diffusion,
        workspace,
        super::face::Trace {
            densities: &left[..n],
            primitive: pl,
            transverse_velocity: 0.0,
            cv_density: cvl,
        },
        super::face::Trace {
            densities: &right[..n],
            primitive: pr,
            transverse_velocity: 0.0,
            cv_density: cvr,
        },
        wl,
        wr,
        [0.0; 2],
        mass_ratio,
        binary,
        out,
        controls,
        cancelled,
    )?;
    stage.gradient_residual = stage
        .gradient_residual
        .max(evaluated.maximum_gradient_residual_per_m);
    stage.mass_frame_residual = stage
        .mass_frame_residual
        .max(evaluated.mass_frame_residual_kg_m2_s.abs());
    let rate = evaluated.diffusivity_m2_s / (0.5 * wl + 0.5 * wr);
    if !rate.is_finite() || rate < 0.0 {
        return Err(DissipationError::Invalid("molecular face scale overflow"));
    }
    Ok(rate)
}

/// Original-source volumetric heat capacity shared by branches and compact gas.
pub(in crate::multicomponent) fn heat_capacity_density(
    gas: &MixtureGas,
    state: &[f64],
    primitive: Primitive,
) -> Result<f64, DissipationError> {
    let n = gas.species_count();
    let mut capacity = [0.0; MAX_SPECIES];
    for k in 0..n {
        if state[k] == 0.0 {
            continue;
        }
        let model = &gas.source().models()[k];
        capacity[k] = (state[k] / model.molar_mass().value())
            * (model
                .evaluate_properties(Temperature::new(primitive.temperature))
                .map_err(super::FlowError::from)?
                .cp()
                .value()
                - R);
    }
    let cv = sum(capacity[..n].iter().copied());
    if !positive(cv) {
        return Err(DissipationError::Invalid("unrepresentable cell Cv density"));
    }
    Ok(cv)
}

#[allow(clippy::too_many_arguments)]
fn wall_face<L: TransportLaw>(
    gas: &MixtureGas,
    law: &mut L,
    state: &[f64],
    p: Primitive,
    distance: f64,
    left: bool,
    wall_velocity: f64,
    binary: &mut [f64],
    out: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, DissipationError> {
    let n = gas.species_count();
    let gradient = if left {
        (p.velocity - wall_velocity) / distance
    } else {
        (wall_velocity - p.velocity) / distance
    };
    let evaluated = super::face::wall::<false, L>(
        gas,
        law,
        super::face::Trace {
            densities: &state[..n],
            primitive: p,
            transverse_velocity: 0.0,
            cv_density: 0.0,
        },
        gradient,
        0.0,
        binary,
        out,
        cancelled,
    )?;
    // The laboratory total-energy flux carries the SAME stress evaluated
    // above times the moving wall's speed. No additional viscous heat source.
    out[n + 1] = if wall_velocity == 0.0 {
        0.0
    } else {
        wall_velocity * out[n]
    };
    let rate = evaluated.diffusivity_m2_s / distance;
    if !rate.is_finite() || !out[n + 1].is_finite() {
        return Err(DissipationError::Invalid("wall traction overflow"));
    }
    Ok(rate)
}

pub(crate) fn fractions(
    gas: &MixtureGas,
    densities: &[f64],
    out: &mut [f64],
) -> Result<(), DissipationError> {
    for (k, model) in gas.source().models().iter().enumerate() {
        out[k] = densities[k] / model.molar_mass().value();
    }
    let concentration = sum(out.iter().copied());
    if !positive(concentration) {
        return Err(DissipationError::Invalid("molecular face concentration"));
    }
    for x in out.iter_mut() {
        *x /= concentration;
    }
    let pivot = (0..out.len())
        .max_by(|&a, &b| out[a].total_cmp(&out[b]))
        .expect("nonempty source bank");
    close_simplex(out, pivot)
}
pub(super) fn close_simplex(x: &mut [f64], pivot: usize) -> Result<(), DissipationError> {
    let remainder = 1.0 - sum((0..x.len()).filter(|&k| k != pivot).map(|k| x[k]));
    if !positive(remainder)
        || remainder > 1.0
        || (remainder - x[pivot]).abs() > 64.0 * f64::EPSILON * x.len() as f64
    {
        return Err(DissipationError::Invalid(
            "derived molecular face fractions do not close within roundoff",
        ));
    }
    x[pivot] = remainder;
    Ok(())
}
