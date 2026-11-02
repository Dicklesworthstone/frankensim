//! One shared multiple-scattering engine. K selects the original real or
//! complex outgoing kernel. Pressure-coordinate scaling keeps real rho*c_eq;
//! the physical derivative remains q=i*rho*omega*v.
//! Native real-k loops and arithmetic order are retained unchanged.
use super::*;

/// Original body metadata shared by real and molecular transmission owners.
pub(in crate::helmholtz::transmission) trait BoundaryBlock {
    fn surface(&self) -> &SpherePanels;
    fn first_panel(&self) -> usize;
    fn first_unknown(&self) -> usize;
    fn contrast(&self) -> &[C64];
}

/// Geometry and peak storage admission before any quadratic allocation.
pub(in crate::helmholtz::transmission) fn admit<'a>(
    surfaces: impl Iterator<Item = &'a SpherePanels>,
    body_count: usize,
    descriptor_bytes: usize,
    minimum_body_clearance_m: f64,
    options: Options,
    gate: &CancelGate,
) -> Result<(usize, usize), TransmissionError> {
    let mut count = 0usize;
    let mut local_entries = 0usize;
    let mut boxes: Vec<([f64; 3], [f64; 3])> = Vec::with_capacity(body_count);
    for surface in surfaces {
        poll(gate)?;
        let n = surface.areas().len();
        count = count
            .checked_add(n)
            .ok_or_else(|| bad("multiple-fluid panel overflow"))?;
        if count > options.maximum_panels {
            return Err(HelmholtzError::WorkCap { panels: count }.into());
        }
        surface
            .validate()
            .map_err(|_| bad("invalid multiple-fluid surface"))?;
        let triangles = surface
            .triangles()
            .ok_or_else(|| bad("multiple-fluid bodies require retained triangles"))?;
        if near_field::components(triangles)?.len() != 1 {
            return Err(bad("each fluid body requires one closed outward component").into());
        }
        let mut lower = [f64::INFINITY; 3];
        let mut upper = [f64::NEG_INFINITY; 3];
        for vertex in triangles.iter().flatten() {
            for axis in 0..3 {
                lower[axis] = lower[axis].min(vertex[axis]);
                upper[axis] = upper[axis].max(vertex[axis]);
            }
        }
        for &(lo, hi) in &boxes {
            let mut gap2 = 0.0;
            for axis in 0..3 {
                let gap = (lower[axis] - hi[axis])
                    .max(lo[axis] - upper[axis])
                    .max(0.0);
                gap2 += gap * gap;
            }
            let scale = lower
                .iter()
                .chain(&upper)
                .chain(&lo)
                .chain(&hi)
                .fold(0.0_f64, |a, &b| a.max(b.abs()));
            if gap2.sqrt() < minimum_body_clearance_m.max(128.0 * f64::EPSILON * scale) {
                return Err(bad("fluid body bounding boxes overlap or violate separation").into());
            }
        }
        boxes.push((lower, upper));
        local_entries += 4 * n * n;
    }
    drop(boxes);
    let m = 2 * count;
    let workers = DENSE_WORKERS.with(core::cell::Cell::get);
    let scratch = lu_complex_blocked_workspace_bytes(m, workers)
        .ok_or_else(|| bad("multiple-fluid LU workspace byte count overflow"))?;
    let complex_bytes = core::mem::size_of::<C64>();
    let transfer_entries = if body_count == 1 { 0 } else { m * m };
    // Factorization clones A; retained C_i and E are live throughout.
    // The solve peak includes incident, solution, exciting, output traces.
    let fixed = (local_entries + transfer_entries) * complex_bytes
        + m * core::mem::size_of::<usize>()
        + count * 16 * core::mem::size_of::<f64>()
        + body_count * descriptor_bytes;
    let factor_peak = fixed
        .checked_add((2 * m * m + m) * complex_bytes)
        .and_then(|x| x.checked_add(scratch))
        .ok_or_else(|| bad("multiple-fluid factor byte count overflow"))?;
    let solve_peak = fixed
        .checked_add((m * m + 8 * m) * complex_bytes)
        .and_then(|x| {
            x.checked_add(
                body_count * core::mem::size_of::<(RadiationPowerDiagnostics, (f64, f64))>(),
            )
        })
        .ok_or_else(|| bad("multiple-fluid solve byte count overflow"))?;
    let required_bytes = factor_peak.max(solve_peak);
    if required_bytes > options.maximum_bytes {
        return Err(HelmholtzError::CacheBudget {
            required_bytes,
            maximum_bytes: options.maximum_bytes,
        }
        .into());
    }

    Ok((count, workers))
}

/// One cached coupled factorization, preserving every mutual reflection.
pub(in crate::helmholtz::transmission) fn prepare<B: BoundaryBlock, K: WaveNumber>(
    bodies: &[B],
    wave: K,
    k_reference: f64,
    m: usize,
    workers: usize,
    gate: &CancelGate,
) -> Result<(Vec<C64>, Option<LuComplex>, f64), HelmholtzError> {
    let transparent = bodies.iter().all(|body| body.contrast().is_empty());
    let (transfer, factor, condition_lower_bound) = if transparent {
        (Vec::new(), None, 1.0)
    } else {
        let transfer = if bodies.len() == 1 {
            Vec::new()
        } else {
            assemble_transfer(bodies, wave, k_reference, m, gate)?
        };
        let mut matrix = zero_matrix(m * m, "multiple-fluid matrix")?;
        for body in bodies {
            let mi = 2 * body.surface().areas().len();
            for row in 0..mi {
                poll(gate)?;
                let global_row = body.first_unknown() + row;
                matrix[global_row * m + global_row] = C64::ONE;
                if body.contrast().is_empty() {
                    continue;
                }
                for col in 0..mi {
                    matrix[global_row * m + body.first_unknown() + col] = matrix
                        [global_row * m + body.first_unknown() + col]
                        + body.contrast()[row * mi + col];
                }
                for other in bodies {
                    // An identical-fluid body's exact equation is s_j=0.
                    // Eliminate those zero columns algebraically; its
                    // interior still receives every active body's field.
                    if other.first_unknown() == body.first_unknown() || other.contrast().is_empty()
                    {
                        continue;
                    }
                    let mj = 2 * other.surface().areas().len();
                    for col in 0..mj {
                        if col % 64 == 0 {
                            poll(gate)?;
                        }
                        let global_col = other.first_unknown() + col;
                        let mut value = C64::ZERO;
                        for inner in 0..mi {
                            value = value
                                + body.contrast()[row * mi + inner]
                                    * transfer[(body.first_unknown() + inner) * m + global_col];
                        }
                        matrix[global_row * m + global_col] = value;
                    }
                }
            }
        }
        if matrix.iter().any(|&z| !finite(z)) {
            return Err(bad("nonfinite multiple-fluid coupled operator"));
        }
        let factor = lu_complex_blocked_with_cancel(&matrix, m, workers, || gate.is_requested())
            .map_err(lu_error)?;
        poll(gate)?;
        let condition = condition_lower_bound(&matrix, m, &factor);
        if !condition.is_finite() {
            return Err(bad("nonfinite multiple-fluid condition diagnostic"));
        }
        (transfer, Some(factor), condition)
    };

    Ok((transfer, factor, condition_lower_bound))
}

/// Local outgoing coefficients are distinct from the complete physical field.
pub(in crate::helmholtz::transmission) struct Traces {
    pub(in crate::helmholtz::transmission) pressure: Vec<C64>,
    pub(in crate::helmholtz::transmission) velocity: Vec<C64>,
    pub(in crate::helmholtz::transmission) scattered_pressure: Vec<C64>,
    pub(in crate::helmholtz::transmission) scattered_velocity: Vec<C64>,
    pub(in crate::helmholtz::transmission) scattered_work: RadiationPowerDiagnostics,
    pub(in crate::helmholtz::transmission) scattered_work_interval: (f64, f64),
    pub(in crate::helmholtz::transmission) total_pressure: Vec<C64>,
    pub(in crate::helmholtz::transmission) total_velocity: Vec<C64>,
    pub(in crate::helmholtz::transmission) work: Vec<(RadiationPowerDiagnostics, (f64, f64))>,
}

#[allow(clippy::too_many_arguments)]
pub(in crate::helmholtz::transmission) fn solve<B: BoundaryBlock>(
    bodies: &[B],
    surface: &SpherePanels,
    exterior: Medium,
    pressure_impedance: f64,
    transfer: &[C64],
    factor: Option<&LuComplex>,
    incident_pressure: &[C64],
    incident_velocity: &[C64],
    gate: &CancelGate,
) -> Result<Traces, HelmholtzError> {
    let n = surface.areas().len();
    let m = 2 * n;
    let z = pressure_impedance;
    let mut incoming = vec![C64::ZERO; m];
    for body in bodies {
        let ni = body.surface().areas().len();
        for j in 0..ni {
            incoming[body.first_unknown() + j] = incident_pressure[body.first_panel() + j];
            incoming[body.first_unknown() + ni + j] =
                incident_velocity[body.first_panel() + j] * C64::new(0.0, z);
        }
    }
    if incoming.iter().any(|&value| !finite(value)) {
        return Err(bad(
            "nonfinite pressure-scaled multiple-fluid incident field",
        ));
    }
    let mut outgoing = vec![C64::ZERO; m];
    if let Some(factor) = factor {
        for body in bodies {
            if body.contrast().is_empty() {
                continue;
            }
            let mi = 2 * body.surface().areas().len();
            for row in 0..mi {
                poll(gate)?;
                for col in 0..mi {
                    outgoing[body.first_unknown() + row] = outgoing[body.first_unknown() + row]
                        - body.contrast()[row * mi + col] * incoming[body.first_unknown() + col];
                }
            }
        }
        if outgoing.iter().any(|&v| !finite(v)) {
            return Err(bad("nonfinite multiple-fluid contrast drive"));
        }
        factor
            .solve_with_cancel(&mut outgoing, || gate.is_requested())
            .map_err(lu_error)?;
    }
    let mut total_scattered = outgoing.clone();
    if !transfer.is_empty() {
        for row in 0..m {
            poll(gate)?;
            for (col, &value) in outgoing.iter().enumerate() {
                total_scattered[row] = total_scattered[row] + transfer[row * m + col] * value;
            }
        }
    }
    let mut pressure = vec![C64::ZERO; n];
    let mut velocity = vec![C64::ZERO; n];
    let mut scattered_pressure = vec![C64::ZERO; n];
    let mut scattered_velocity = vec![C64::ZERO; n];
    let mut total_pressure = vec![C64::ZERO; n];
    let mut total_velocity = vec![C64::ZERO; n];
    let mut work = Vec::with_capacity(bodies.len());
    for body in bodies {
        poll(gate)?;
        let ni = body.surface().areas().len();
        for j in 0..ni {
            let (p, q) = (body.first_unknown() + j, body.first_unknown() + ni + j);
            let target = body.first_panel() + j;
            pressure[target] = outgoing[p];
            velocity[target] = if body.contrast().is_empty() {
                C64::ZERO
            } else {
                C64::new(outgoing[q].im / z, -outgoing[q].re / z)
            };
            scattered_pressure[target] = total_scattered[p];
            scattered_velocity[target] = if factor.is_none() {
                C64::ZERO
            } else {
                C64::new(total_scattered[q].im / z, -total_scattered[q].re / z)
            };
            if factor.is_none() {
                total_pressure[target] = incident_pressure[target];
                total_velocity[target] = incident_velocity[target];
            } else {
                total_pressure[target] = total_scattered[p] + incoming[p];
                let total_q = total_scattered[q] + incoming[q];
                total_velocity[target] = C64::new(total_q.im / z, -total_q.re / z);
            }
        }
        let range = body.first_panel()..body.first_panel() + ni;
        let (_, interval, diagnostics) = evaluate_radiated_power(
            &total_pressure[range.clone()],
            &total_velocity[range],
            body.surface().areas(),
            exterior,
        );
        if !finite(diagnostics.surface_power) || !diagnostics.apparent_power.is_finite() {
            return Err(bad("nonfinite multiple-fluid interface work"));
        }
        work.push((diagnostics, interval));
    }
    if pressure
        .iter()
        .chain(&velocity)
        .chain(&total_pressure)
        .chain(&total_velocity)
        .chain(&scattered_pressure)
        .chain(&scattered_velocity)
        .any(|&v| !finite(v))
    {
        return Err(bad("nonfinite multiple-fluid solution traces"));
    }

    let (_, scattered_work_interval, scattered_work) = evaluate_radiated_power(
        &scattered_pressure,
        &scattered_velocity,
        surface.areas(),
        exterior,
    );
    if !finite(scattered_work.surface_power) || !scattered_work.apparent_power.is_finite() {
        return Err(bad("nonfinite complete multiple-fluid scattered work"));
    }

    poll(gate)?;
    Ok(Traces {
        pressure,
        velocity,
        scattered_pressure,
        scattered_velocity,
        scattered_work,
        scattered_work_interval,
        total_pressure,
        total_velocity,
        work,
    })
}

fn zero_matrix(entries: usize, what: &'static str) -> Result<Vec<C64>, HelmholtzError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(entries)
        .map_err(|_| HelmholtzError::AllocationFailed { what })?;
    values.resize(entries, C64::ZERO);
    Ok(values)
}

fn assemble_transfer<B: BoundaryBlock, K: WaveNumber>(
    bodies: &[B],
    k: K,
    k_reference: f64,
    m: usize,
    gate: &CancelGate,
) -> Result<Vec<C64>, HelmholtzError> {
    let mut transfer = zero_matrix(m * m, "multiple-fluid exterior transfer")?;
    for target in bodies {
        let a = target.surface();
        let ni = a.areas().len();
        for source in bodies {
            if source.first_unknown() == target.first_unknown() || source.contrast().is_empty() {
                continue;
            }
            let b = source.surface();
            let nj = b.areas().len();
            for i in 0..ni {
                poll(gate)?;
                for j in 0..nj {
                    let (s, d, dp, n) = kernels(
                        k,
                        a.centroids()[i],
                        a.normals()[i],
                        b.centroids()[j],
                        b.normals()[j],
                    );
                    let (r, c) = (target.first_unknown() + i, source.first_unknown() + j);
                    let area = b.areas()[j];
                    transfer[r * m + c] = d.scale(area);
                    transfer[r * m + c + nj] = s.scale(-k_reference * area);
                    transfer[(r + ni) * m + c] = n.scale(area / k_reference);
                    transfer[(r + ni) * m + c + nj] = dp.scale(-area);
                }
            }
        }
    }
    if transfer.iter().any(|&v| !finite(v)) {
        return Err(bad("nonfinite multiple-fluid exterior transfer"));
    }
    Ok(transfer)
}
