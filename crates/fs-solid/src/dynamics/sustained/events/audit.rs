//! Original-operator, continuous phase and complete-coordinate work gates.

use super::*;

pub(super) fn rigid_response<E, N>(
    before: &State<'_>,
    damping: Option<&[f64]>,
    weights: [&[f64]; 3],
    response: &TraceResponse,
    cx: &Cx<'_>,
) -> Result<(), SolidSustainedError<E, N>> {
    let n = before.displacement.len();
    for (column, weight) in response.columns.iter().zip(weights) {
        for i in 0..n {
            poll(cx)?;
            let row = &before.mass[i * n..(i + 1) * n];
            let product = dot(row, column);
            let magnitude = row.iter().zip(column).fold(Dd::ZERO, |sum, (m, b)| {
                sum + (Dd::from_f64(*m) * Dd::from_f64(*b)).abs()
            });
            gate(
                "original contact mass response",
                product - Dd::from_f64(weight[i]),
                rounding(&[magnitude, Dd::from_f64(weight[i])]),
            )?;
            let stiffness = dot(&before.stiffness[i * n..(i + 1) * n], column);
            let viscosity = damping.map_or(Dd::ZERO, |c| dot(&c[i * n..(i + 1) * n], column));
            finite(stiffness)?;
            finite(viscosity)?;
            // An endpoint coincidence cannot certify an interior trace. This
            // exact nullspace requirement is deliberately stronger than the
            // native preparation's approximate operator residual allowance.
            if stiffness.hi != 0.0
                || stiffness.lo != 0.0
                || viscosity.hi != 0.0
                || viscosity.lo != 0.0
            {
                return Err("contact mass response must lie in the original stiffness and viscosity nullspaces".into());
            }
        }
    }
    Ok(())
}

pub(super) fn free_trace<E, N>(
    before: &State<'_>,
    q: &[f64],
    v: &[f64],
    c: &SolidSustainedContact<'_>,
    initial: [f64; 3],
    acceleration: [f64; 3],
    controls: SolidSustainedControls,
) -> Result<(), SolidSustainedError<E, N>> {
    let h = Dd::from_f64(before.dt_s);
    for a in 0..3 {
        let relative = q
            .iter()
            .zip(before.displacement)
            .zip(c.velocity_weights[a])
            .fold(
                -h * Dd::from_f64(c.prescribed_surface_velocity_mps[a]),
                |sum, ((q, old), w)| {
                    sum + Dd::from_f64(*w) * (Dd::from_f64(*q) - Dd::from_f64(*old))
                },
            );
        let expected =
            (Dd::from_f64(initial[a]) + Dd::from_f64(0.5) * Dd::from_f64(acceleration[a]) * h) * h;
        gate(
            "native rigid trace displacement",
            relative - expected,
            Dd::from_f64(controls.displacement_tolerance_m) + rounding(&[relative, expected]),
        )?;
        let actual =
            dot(c.velocity_weights[a], v) - Dd::from_f64(c.prescribed_surface_velocity_mps[a]);
        let expected = Dd::from_f64(initial[a]) + Dd::from_f64(acceleration[a]) * h;
        gate(
            "native rigid trace velocity",
            actual - expected,
            Dd::from_f64(controls.velocity_tolerance_mps) + rounding(&[actual, expected]),
        )?;
    }
    Ok(())
}

pub(super) fn phase<E, N>(
    p: &SolidCoulombPhase,
    axis: [f64; 2],
    id: u64,
    controls: SolidSustainedControls,
) -> Result<(), SolidSustainedError<E, N>> {
    let force = p.force_n;
    let cap = Dd::from_f64(p.friction_coefficient) * Dd::from_f64(force[0]);
    let capacity = finite(cap)?;
    let ft = force[1].hypot(force[2]);
    let travel = p.relative_displacement_m;
    let distance = travel[1].hypot(travel[2]);
    let v0 = p.relative_velocity_before_mps;
    let v1 = p.relative_velocity_after_mps;
    if force[0] <= 0.0
        || !ft.is_finite()
        || !distance.is_finite()
        || ft > capacity + controls.force_tolerance_n
    {
        return Err(mode(
            id,
            "compressive force inside the actual Coulomb cone required",
        ));
    }
    let displacement_allowance = controls.displacement_tolerance_m;
    normal_path(
        id,
        force[0],
        [p.normal_gap_before_m, p.normal_gap_after_m],
        [v0[0], v1[0]],
        p.duration_s,
        controls,
    )?;
    tangent_path(
        id,
        p.mode,
        [[v0[1], v0[2]], [v1[1], v1[2]]],
        [travel[1], travel[2]],
        p.duration,
        cap,
        controls,
    )?;
    let complementarity = Dd::from_f64(controls.complementarity_tolerance_j);
    for velocity in [v0, v1] {
        phase::collinear::<E, N>(id, axis, velocity, controls.velocity_tolerance_mps)?;
        let speed = velocity[1].hypot(velocity[2]);
        if !speed.is_finite() {
            return Err("unrepresentable event tangent speed".into());
        }
        if p.mode == CoulombPhaseMode::Sticking {
            if speed > controls.velocity_tolerance_mps || distance > displacement_allowance {
                return Err(mode(id, "static segment requires vanishing tangent motion"));
            }
        } else {
            let power = dot(&force[1..], &velocity[1..]);
            let boundary_power = cap * Dd::from_f64(speed);
            gate(
                "continuous kinetic sliding direction",
                power + boundary_power,
                cap * Dd::from_f64(controls.velocity_tolerance_mps)
                    + Dd::from_f64(controls.force_tolerance_n) * Dd::from_f64(speed)
                    + rounding(&[power, boundary_power]),
            )?;
        }
    }
    if p.mode == CoulombPhaseMode::Sliding {
        // The force is constant and velocity affine; opposing both endpoints
        // excludes an unlocated interior reversal in this collinear segment.
        gate(
            "kinetic force magnitude",
            Dd::from_f64(ft) - cap,
            Dd::from_f64(controls.force_tolerance_n) + rounding(&[cap]),
        )?;
    }
    let tangent_work = dot(&force[1..], &travel[1..]);
    let boundary_work = cap * Dd::from_f64(distance);
    gate(
        "event tangent displacement duality",
        tangent_work + boundary_work,
        complementarity + rounding(&[tangent_work, boundary_work]),
    )?;
    if finite(tangent_work)? > controls.maximum_work_error_j {
        return Err(mode(
            id,
            "relative friction creates energy during an event segment",
        ));
    }
    Ok(())
}

fn normal_path<E, N>(
    id: u64,
    force: f64,
    gap: [f64; 2],
    velocity: [f64; 2],
    duration: f64,
    controls: SolidSustainedControls,
) -> Result<(), SolidSustainedError<E, N>> {
    if velocity
        .iter()
        .any(|v| v.abs() > controls.velocity_tolerance_mps)
    {
        return Err(mode(
            id,
            "loaded contact must dwell throughout every event segment",
        ));
    }
    let mut gaps = [gap[0], gap[1], gap[0]];
    let change = velocity[1] - velocity[0];
    if change != 0.0 {
        let fraction = -velocity[0] / change;
        if fraction > 0.0 && fraction < 1.0 {
            // Check the stationary point for either sign of acceleration:
            // loaded contact must neither penetrate nor open and return.
            gaps[2] = finite(
                Dd::from_f64(gap[0])
                    + Dd::from_f64(0.5)
                        * Dd::from_f64(velocity[0])
                        * Dd::from_f64(duration)
                        * Dd::from_f64(fraction),
            )?;
        }
    }
    for gap in gaps {
        if gap.abs() > controls.displacement_tolerance_m {
            return Err(mode(
                id,
                "interior or endpoint normal gap violates loaded dwell",
            ));
        }
        let work = Dd::from_f64(force) * Dd::from_f64(gap);
        gate(
            "event normal force-gap complementarity",
            work,
            Dd::from_f64(controls.complementarity_tolerance_j) + rounding(&[work]),
        )?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn tangent_path<E, N>(
    id: u64,
    phase_mode: CoulombPhaseMode,
    velocity: [[f64; 2]; 2],
    travel: [f64; 2],
    duration: Dd,
    capacity: Dd,
    controls: SolidSustainedControls,
) -> Result<(), SolidSustainedError<E, N>> {
    let mut extent = travel.map(f64::abs);
    for a in 0..2 {
        let change = velocity[1][a] - velocity[0][a];
        if change == 0.0 {
            continue;
        }
        let fraction = -velocity[0][a] / change;
        if fraction <= 0.0 || fraction >= 1.0 {
            continue;
        }
        let fraction = Dd::from_f64(fraction);
        let mut vertex = [0.0; 2];
        for b in 0..2 {
            vertex[b] = finite(
                duration
                    * fraction
                    * (Dd::from_f64(velocity[0][b])
                        + Dd::from_f64(0.5)
                            * fraction
                            * (Dd::from_f64(velocity[1][b]) - Dd::from_f64(velocity[0][b]))),
            )?;
        }
        extent[a] = extent[a].max(vertex[a].abs());
        if phase_mode == CoulombPhaseMode::Sliding {
            // A tiny stop residual is not reset. It may only create unresolved
            // out-and-back travel, never a finite reversal hidden by net slip.
            let arc = Dd::from_f64(vertex[0].hypot(vertex[1]))
                + Dd::from_f64((travel[0] - vertex[0]).hypot(travel[1] - vertex[1]));
            let net = Dd::from_f64(travel[0].hypot(travel[1]));
            let excess = arc - net;
            if finite(excess)? > controls.displacement_tolerance_m {
                return Err(mode(
                    id,
                    "resolved interior kinetic reversal requires another event",
                ));
            }
            gate(
                "unlocated kinetic reversal work",
                capacity * excess,
                Dd::from_f64(controls.complementarity_tolerance_j)
                    + rounding(&[capacity * arc, capacity * net]),
            )?;
        }
    }
    // Component extrema give a conservative whole-path norm bound, including
    // residual transverse motion within the authored-axis velocity allowance.
    if phase_mode == CoulombPhaseMode::Sticking {
        let excursion = extent[0].hypot(extent[1]);
        if excursion > controls.displacement_tolerance_m {
            return Err(mode(id, "static segment contains resolved interior slip"));
        }
        let work = capacity * Dd::from_f64(excursion);
        gate(
            "static interior slip work",
            work,
            Dd::from_f64(controls.complementarity_tolerance_j) + rounding(&[work]),
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn controls() -> SolidSustainedControls {
        SolidSustainedControls {
            solve: Default::default(),
            maximum_workspace_bytes: 1_000_000,
            maximum_response_multiply_terms: 1_000_000,
            displacement_tolerance_m: 1e-10,
            velocity_tolerance_mps: 1e-9,
            force_tolerance_n: 1e-8,
            complementarity_tolerance_j: 1e-10,
            maximum_work_error_j: 1e-9,
        }
    }

    #[test]
    fn loaded_dwell_checks_interior_opening_as_well_as_penetration() {
        let controls = controls();
        for sign in [-1.0, 1.0] {
            let velocity = [sign * 8e-10, -sign * 8e-10];
            assert!(matches!(
                normal_path::<(), ()>(7, 20.0, [0.0; 2], velocity, 4.0, controls),
                Err(SolidSustainedError::ContactMode { .. })
            ));
            // Permissive geometry still cannot hide loaded separation from the
            // independent force-gap complementarity allowance.
            assert!(matches!(
                normal_path::<(), ()>(
                    7,
                    20.0,
                    [0.0; 2],
                    velocity,
                    4.0,
                    SolidSustainedControls {
                        displacement_tolerance_m: 1e-8,
                        ..controls
                    }
                ),
                Err(SolidSustainedError::Balance { .. })
            ));
        }
    }

    #[test]
    fn static_capture_cannot_hide_resolved_out_and_back_travel() {
        assert!(matches!(
            tangent_path::<(), ()>(
                7,
                CoulombPhaseMode::Sticking,
                [[8e-10, 0.0], [-8e-10, 0.0]],
                [0.0; 2],
                Dd::from_f64(4.0),
                Dd::from_f64(12.0),
                controls()
            ),
            Err(SolidSustainedError::ContactMode { .. })
        ));
        assert!(matches!(
            tangent_path::<(), ()>(
                7,
                CoulombPhaseMode::Sticking,
                [[8e-10, 0.0], [-8e-10, 0.0]],
                [0.0; 2],
                Dd::from_f64(4.0),
                Dd::from_f64(12.0),
                SolidSustainedControls {
                    displacement_tolerance_m: 1e-8,
                    ..controls()
                }
            ),
            Err(SolidSustainedError::Balance { .. })
        ));
    }

    #[test]
    fn kinetic_stop_residual_cannot_hide_resolved_reversal_loss() {
        let velocity = [[4e-10, 0.0], [-8e-10, 0.0]];
        assert!(matches!(
            tangent_path::<(), ()>(
                7,
                CoulombPhaseMode::Sliding,
                velocity,
                [-8e-10, 0.0],
                Dd::from_f64(4.0),
                Dd::from_f64(6.0),
                controls()
            ),
            Err(SolidSustainedError::ContactMode { .. })
        ));
        assert!(matches!(
            tangent_path::<(), ()>(
                7,
                CoulombPhaseMode::Sliding,
                velocity,
                [-8e-10, 0.0],
                Dd::from_f64(4.0),
                Dd::from_f64(6.0),
                SolidSustainedControls {
                    displacement_tolerance_m: 1e-8,
                    ..controls()
                }
            ),
            Err(SolidSustainedError::Balance { .. })
        ));
    }
}

pub(super) struct Measured {
    pub momentum: Vec<f64>,
    pub distance: f64,
    pub external: f64,
    pub surface: f64,
    pub normal: f64,
    pub friction: f64,
    pub supplied: f64,
    pub energy_error: f64,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn final_path<E, N>(
    before: &State<'_>,
    q: &[f64],
    v: &[f64],
    baseline_velocity: &[f64],
    du: &[f64],
    c: &SolidSustainedContact<'_>,
    phases: &[SolidCoulombPhase],
    impulse: [Dd; 3],
    external: &[f64],
    viscosity: f64,
    controls: SolidSustainedControls,
    cx: &Cx<'_>,
) -> Result<Measured, SolidSustainedError<E, N>> {
    let n = before.displacement.len();
    let h = Dd::from_f64(before.dt_s);
    let (mut surface, mut normal, mut friction, mut distance) =
        (Dd::ZERO, Dd::ZERO, Dd::ZERO, Dd::ZERO);
    let mut relative = [Dd::ZERO; 3];
    let mut elapsed = Dd::ZERO;
    for p in phases {
        poll(cx)?;
        gate(
            "contiguous event segments",
            Dd::from_f64(p.start_offset_s) - elapsed,
            rounding(&[elapsed]),
        )?;
        let dt = p.duration;
        elapsed = elapsed + dt;
        surface = surface + dt * dot(&p.force_n, &c.prescribed_surface_velocity_mps);
        normal = normal + Dd::from_f64(p.force_n[0]) * Dd::from_f64(p.relative_displacement_m[0]);
        friction = friction + dot(&p.force_n[1..], &p.relative_displacement_m[1..]);
        if p.mode == CoulombPhaseMode::Sliding {
            distance = distance
                + Dd::from_f64(p.relative_displacement_m[1].hypot(p.relative_displacement_m[2]));
        }
        for (a, value) in relative.iter_mut().enumerate() {
            *value = *value + Dd::from_f64(p.relative_displacement_m[a]);
        }
    }
    gate(
        "complete event duration",
        elapsed - h,
        rounding(&[elapsed, h]),
    )?;
    let last = phases.last().ok_or("missing actual Coulomb phase")?;
    for a in 0..3 {
        let actual =
            dot(c.velocity_weights[a], du) - h * Dd::from_f64(c.prescribed_surface_velocity_mps[a]);
        gate(
            "complete native event displacement",
            actual - relative[a],
            Dd::from_f64(controls.displacement_tolerance_m) + rounding(&[actual, relative[a]]),
        )?;
        let actual =
            dot(c.velocity_weights[a], v) - Dd::from_f64(c.prescribed_surface_velocity_mps[a]);
        let expected = Dd::from_f64(last.relative_velocity_after_mps[a]);
        gate(
            "complete native event velocity",
            actual - expected,
            Dd::from_f64(controls.velocity_tolerance_mps) + rounding(&[actual, expected]),
        )?;
    }
    let mut momentum = Vec::with_capacity(n);
    let mut change = Dd::ZERO;
    for i in 0..n {
        poll(cx)?;
        let (mut dp, mut contact_dp, mut elastic, mut magnitude) =
            (Dd::ZERO, Dd::ZERO, Dd::ZERO, Dd::ZERO);
        for j in 0..n {
            let mass = Dd::from_f64(before.mass[i * n + j]);
            dp = dp + mass * (Dd::from_f64(v[j]) - Dd::from_f64(before.velocity[j]));
            contact_dp =
                contact_dp + mass * (Dd::from_f64(v[j]) - Dd::from_f64(baseline_velocity[j]));
            magnitude = magnitude
                + mass.abs()
                    * (Dd::from_f64(v[j]).abs() + Dd::from_f64(baseline_velocity[j]).abs());
            elastic = elastic
                + Dd::from_f64(before.stiffness[i * n + j])
                    * (Dd::from_f64(q[j]) - Dd::from_f64(before.displacement[j]));
        }
        let expected = (0..3).fold(Dd::ZERO, |sum, a| {
            sum + Dd::from_f64(c.velocity_weights[a][i]) * impulse[a]
        });
        gate(
            "original mass and event impulse",
            contact_dp - expected,
            Dd::from_f64(controls.force_tolerance_n) * h + rounding(&[magnitude, expected]),
        )?;
        momentum.push(finite(dp)?);
        change = change
            + Dd::from_f64(0.5) * (Dd::from_f64(v[i]) + Dd::from_f64(before.velocity[i])) * dp
            + Dd::from_f64(0.5)
                * (Dd::from_f64(q[i]) + Dd::from_f64(before.displacement[i]))
                * elastic;
    }
    let external = dot(external, du);
    let supplied = external + surface + normal + friction;
    let energy_error = change + Dd::from_f64(viscosity) - supplied;
    gate(
        "original event energy and all force work",
        energy_error,
        Dd::from_f64(controls.maximum_work_error_j)
            + rounding(&[
                change,
                Dd::from_f64(viscosity),
                external,
                surface,
                normal,
                friction,
            ]),
    )?;
    Ok(Measured {
        momentum,
        distance: finite(distance)?,
        external: finite(external)?,
        surface: finite(surface)?,
        normal: finite(normal)?,
        friction: finite(friction)?,
        supplied: finite(supplied)?,
        energy_error: finite(energy_error)?,
    })
}
