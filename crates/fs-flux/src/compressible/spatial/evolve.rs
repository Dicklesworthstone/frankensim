use super::*;

#[allow(clippy::too_many_arguments)]
pub(super) fn attempt(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: &reservoir::Reservoirs,
    states: &[State],
    work: &mut Workspace,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    force_controls: Option<forcing::ForceControl3d>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StageStats, Error> {
    let first = rk_stage(
        gas,
        grid,
        boundaries,
        reservoirs.accepted(),
        states,
        work,
        dt,
        cfl,
        false,
        reconstruction,
        force_controls,
        cancelled,
    )?;
    let second = rk_stage(
        gas,
        grid,
        boundaries,
        reservoirs.accepted(),
        states,
        work,
        dt,
        cfl,
        true,
        reconstruction,
        force_controls,
        cancelled,
    )?;
    Ok(combine(first, second))
}

/// The start profile has already been staged for the initial time estimate.
/// Reuse that same scratch for the endpoint only after its first residual has
/// been consumed; accepted profiles are never swapped into speculative work.
#[allow(clippy::too_many_arguments)]
pub(super) fn attempt_with_reservoir_profiles(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: &mut reservoir::Reservoirs,
    states: &[State],
    work: &mut Workspace,
    source: &mut reservoir::TimedSource<'_>,
    time: f64,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StageStats, Error> {
    let first_rate = reservoirs
        .staged()
        .maximum_rate(gas, grid, states, cancelled)?;
    check_cfl::<Error>(first_rate * dt, cfl)?;
    let first = rk_stage(
        gas,
        grid,
        boundaries,
        reservoirs.staged(),
        states,
        work,
        dt,
        cfl,
        false,
        reconstruction,
        None,
        cancelled,
    )?;
    source.stage(&mut reservoirs.profiles, gas, grid, time + dt, cancelled)?;
    // The actual endpoint trace can be much faster or hotter than the start.
    // Its full unsplit exterior bound and the resulting native face-wave bound
    // must both admit this same dt, including the final SSP residual.
    let second_rate = reservoirs
        .staged()
        .maximum_rate(gas, grid, &work.stage, cancelled)?;
    check_cfl::<Error>(second_rate * dt, cfl)?;
    let second = rk_stage(
        gas,
        grid,
        boundaries,
        reservoirs.staged(),
        states,
        work,
        dt,
        cfl,
        true,
        reconstruction,
        None,
        cancelled,
    )?;
    let mut stats = combine(first, second);
    stats.rate = stats.rate.max(first_rate).max(second_rate);
    Ok(stats)
}

fn combine(first: StageStats, second: StageStats) -> StageStats {
    StageStats {
        rate: first.rate.max(second.rate),
        fallback: first.fallback + second.fallback,
        limited: first.limited + second.limited,
    }
}

#[allow(clippy::too_many_arguments)]
fn rk_stage(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: reservoir::View<'_>,
    states: &[State],
    work: &mut Workspace,
    dt: f64,
    cfl: f64,
    final_stage: bool,
    reconstruction: Reconstruction,
    force_controls: Option<forcing::ForceControl3d>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StageStats, Error> {
    let input = if final_stage { &work.stage[..] } else { states };
    let stats = spatial(
        gas,
        grid,
        boundaries,
        reservoirs,
        input,
        &mut work.primitive,
        &mut work.slopes,
        &mut work.flux,
        reconstruction,
        force_controls,
        cancelled,
    )?;
    check_cfl::<Error>(stats.rate * dt, cfl)?;
    if final_stage {
        update(
            gas,
            grid,
            states,
            &work.stage,
            &work.flux,
            &mut work.candidate,
            dt,
            true,
            force_controls,
            cancelled,
        )?;
    } else {
        update(
            gas,
            grid,
            states,
            states,
            &work.flux,
            &mut work.stage,
            dt,
            false,
            force_controls,
            cancelled,
        )?;
    }
    capture_boundaries(
        grid,
        &work.flux,
        0.5 * dt,
        &mut work.boundary_trial,
        final_stage,
        cancelled,
    )?;
    Ok(stats)
}

#[allow(clippy::too_many_arguments)]
fn spatial(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: reservoir::View<'_>,
    states: &[State],
    primitive: &mut [Primitive],
    slopes: &mut [[[f64; 5]; 3]],
    fluxes: &mut [Vec<SpatialFaceFlux>; 3],
    reconstruction: Reconstruction,
    force_controls: Option<forcing::ForceControl3d>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StageStats, Error> {
    let mut stats = StageStats::default();
    for (p, q) in primitive.iter_mut().zip(states) {
        *p = q.primitive(gas, &mut *cancelled)?;
        if let Some(controls) = force_controls {
            controls.admit_velocity::<Error>(p.velocity, p.sound_speed)?;
        }
    }
    stats.limited = cartesian::slopes::<Error>(
        grid,
        boundaries,
        5,
        [1, 2, 3],
        reconstruction,
        |k, d| primitive[k].array()[d],
        |_, _, _, _| Ok(0),
        slopes.as_flattened_mut().as_flattened_mut(),
        |k, _, sign, slope| {
            let center = primitive[k].array();
            let p: [f64; 5] = std::array::from_fn(|d| center[d] + sign * slope[d]);
            let temperature = p[4] / (gas.gas_constant() * p[0]);
            p.into_iter().all(f64::is_finite)
                && positive(p[0])
                && positive(p[4])
                && temperature >= gas.domain().minimum_temperature
                && temperature <= gas.domain().maximum_temperature
        },
        cancelled,
    )?;
    // Every axis uses a contiguous, x-fast face lattice. Internal neighbours
    // consume the exact same stored flux; periodic endpoints share its bits.
    cartesian::faces::<Error, _>(grid, boundaries, cancelled, |face, cancelled| {
        let axis = face.axis;
        let flux = &mut fluxes[axis];
        let [t1, t2] = tangents(axis);
        if let cartesian::FaceKind::PeriodicCopy { first } = face.kind {
            flux[face.index] = flux[first];
            return Ok(());
        }
        let reconstruct = |k: usize, sign: f64| -> SpatialPrimitive {
            let center = primitive[k].array();
            let p: [f64; 5] = std::array::from_fn(|d| center[d] + sign * slopes[k][axis][d]);
            SpatialPrimitive {
                density: p[0],
                normal_velocity: p[axis + 1],
                transverse_velocity: [p[t1 + 1], p[t2 + 1]],
                pressure: p[4],
            }
        };
        let first_order = reconstruction == Reconstruction::FirstOrder;
        let f = match face.kind {
            cartesian::FaceKind::Interior { left, right } => {
                if first_order {
                    gas.hllc_spatial_flux(
                        states[left].normal(axis),
                        states[right].normal(axis),
                        &mut *cancelled,
                    )?
                } else {
                    gas.hllc_spatial_flux_from_primitive(
                        reconstruct(left, 0.5),
                        reconstruct(right, -0.5),
                        &mut *cancelled,
                    )?
                }
            }
            cartesian::FaceKind::Boundary {
                cell: k,
                sign,
                condition,
            } => {
                if let Some(outlet) = reservoirs.outlet(axis, sign) {
                    if first_order {
                        gas.pressure_outlet_spatial_flux(
                            states[k].normal(axis),
                            sign,
                            outlet.pressure_pa,
                            &mut *cancelled,
                        )?
                    } else {
                        gas.pressure_outlet_spatial_flux_from_primitive(
                            reconstruct(k, 0.5 * sign),
                            sign,
                            outlet.pressure_pa,
                            &mut *cancelled,
                        )?
                    }
                } else if let Some(exterior) = reservoirs.get_at(grid, axis, sign, k) {
                    if first_order {
                        let inside = states[k].normal(axis);
                        let outside = exterior.state.normal(axis);
                        let (left, right) = if sign < 0.0 {
                            (outside, inside)
                        } else {
                            (inside, outside)
                        };
                        gas.hllc_spatial_flux(left, right, &mut *cancelled)?
                    } else {
                        let inside = reconstruct(k, 0.5 * sign);
                        let outside = exterior.normal(axis);
                        let (left, right) = if sign < 0.0 {
                            (outside, inside)
                        } else {
                            (inside, outside)
                        };
                        gas.hllc_spatial_flux_from_primitive(left, right, &mut *cancelled)?
                    }
                } else {
                    match (first_order, condition) {
                        (true, BoundaryCondition::Reflecting) => gas.stationary_wall_spatial_flux(
                            states[k].normal(axis),
                            sign,
                            &mut *cancelled,
                        )?,
                        (true, BoundaryCondition::Transmissive) => {
                            let q = states[k].normal(axis);
                            gas.hllc_spatial_flux(q, q, &mut *cancelled)?
                        }
                        (false, BoundaryCondition::Reflecting) => gas
                            .stationary_wall_spatial_flux_from_primitive(
                                reconstruct(k, 0.5 * sign),
                                sign,
                                &mut *cancelled,
                            )?,
                        (false, BoundaryCondition::Transmissive) => {
                            let p = reconstruct(k, 0.5 * sign);
                            gas.hllc_spatial_flux_from_primitive(p, p, &mut *cancelled)?
                        }
                        (_, BoundaryCondition::PressureOutlet { .. }) => {
                            unreachable!("admitted pressure outlet has an immutable source binding")
                        }
                    }
                }
            }
            cartesian::FaceKind::PeriodicCopy { .. } => unreachable!(),
        };
        stats.fallback += usize::from(f.used_hll_fallback);
        flux[face.index] = f;
        Ok(())
    })?;
    stats.rate = cartesian::rate::<Error>(
        grid,
        |axis, face| fluxes[axis][face].signal_speed,
        cancelled,
    )?;
    Ok(stats)
}

pub(super) fn neighbour(
    p: &[Primitive],
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    c: [usize; 3],
    axis: usize,
    upper: bool,
) -> [f64; 5] {
    let (k, reflected) = cartesian::neighbour(grid, boundaries, grid.index(c), axis, upper);
    let mut q = p[k].array();
    if reflected {
        q[axis + 1] = -q[axis + 1];
    }
    q
}

pub(super) fn flux_array(f: SpatialFaceFlux, axis: usize) -> [f64; 5] {
    let [a, b] = tangents(axis);
    let mut q = [0.0; 5];
    q[0] = f.mass;
    q[axis + 1] = f.normal_momentum;
    q[a + 1] = f.transverse_momentum[0];
    q[b + 1] = f.transverse_momentum[1];
    q[4] = f.energy;
    q
}

#[allow(clippy::too_many_arguments)]
pub(super) fn update(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    original: &[State],
    stage: &[State],
    fluxes: &[Vec<SpatialFaceFlux>; 3],
    output: &mut [State],
    dt: f64,
    final_stage: bool,
    force_controls: Option<forcing::ForceControl3d>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    cartesian::update::<5, Error, _>(
        grid,
        5,
        dt,
        final_stage,
        |k, d| original[k].array()[d],
        |k, d| stage[k].array()[d],
        |axis, face, d| flux_array(fluxes[axis][face], axis)[d],
        cancelled,
        |k, row, cancelled| {
            let q = State::from_array(row.try_into().expect("five conserved components"));
            let primitive = q.primitive(gas, &mut *cancelled)?;
            if let Some(controls) = force_controls {
                controls.admit_velocity::<Error>(primitive.velocity, primitive.sound_speed)?;
            }
            output[k] = q;
            Ok(())
        },
    )
}

pub(super) fn capture_boundaries(
    grid: &CartesianGrid,
    fluxes: &[Vec<SpatialFaceFlux>; 3],
    duration: f64,
    output: &mut [BoundaryImpulse],
    add: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    cartesian::capture::<Error>(
        grid,
        5,
        duration,
        add,
        |axis, face, d| flux_array(fluxes[axis][face], axis)[d],
        |boundary, d, delta, add| {
            let mut q = output[boundary].array();
            q[d] = (if add { q[d] } else { 0.0 }) + delta;
            output[boundary] = Totals::from_array(q);
        },
        cancelled,
    )
}
