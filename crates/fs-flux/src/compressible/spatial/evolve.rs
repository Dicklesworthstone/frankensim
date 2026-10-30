use super::*;

#[allow(clippy::too_many_arguments)]
pub(super) fn attempt(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    states: &[State],
    work: &mut Workspace,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StageStats, Error> {
    let first = spatial(
        gas,
        grid,
        boundaries,
        states,
        &mut work.primitive,
        &mut work.slopes,
        &mut work.flux,
        reconstruction,
        cancelled,
    )?;
    check_cfl::<Error>(first.rate * dt, cfl)?;
    update(
        gas,
        grid,
        states,
        states,
        &work.flux,
        &mut work.stage,
        dt,
        false,
        cancelled,
    )?;
    capture_boundaries(
        grid,
        &work.flux,
        0.5 * dt,
        &mut work.boundary_trial,
        false,
        cancelled,
    )?;
    let second = spatial(
        gas,
        grid,
        boundaries,
        &work.stage,
        &mut work.primitive,
        &mut work.slopes,
        &mut work.flux,
        reconstruction,
        cancelled,
    )?;
    check_cfl::<Error>(second.rate * dt, cfl)?;
    update(
        gas,
        grid,
        states,
        &work.stage,
        &work.flux,
        &mut work.candidate,
        dt,
        true,
        cancelled,
    )?;
    capture_boundaries(
        grid,
        &work.flux,
        0.5 * dt,
        &mut work.boundary_trial,
        true,
        cancelled,
    )?;
    Ok(StageStats {
        rate: first.rate.max(second.rate),
        fallback: first.fallback + second.fallback,
        limited: first.limited + second.limited,
    })
}

#[allow(clippy::too_many_arguments)]
fn spatial(
    gas: &FrozenGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    states: &[State],
    primitive: &mut [Primitive],
    slopes: &mut [[[f64; 5]; 3]],
    fluxes: &mut [Vec<SpatialFaceFlux>; 3],
    reconstruction: Reconstruction,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StageStats, Error> {
    let mut stats = StageStats::default();
    for (p, q) in primitive.iter_mut().zip(states) {
        *p = q.primitive(gas, &mut *cancelled)?;
    }
    stats.limited = cartesian::slopes::<Error>(
        grid,
        boundaries,
        5,
        [1, 2, 3],
        reconstruction,
        |k, d| primitive[k].array()[d],
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
            } => match (first_order, condition) {
                (true, BoundaryCondition::Reflecting) => {
                    gas.stationary_wall_spatial_flux(states[k].normal(axis), sign, &mut *cancelled)?
                }
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
            },
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

fn flux_array(f: SpatialFaceFlux, axis: usize) -> [f64; 5] {
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
            q.primitive(gas, &mut *cancelled)?;
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
