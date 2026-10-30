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
    for k in 0..states.len() {
        if k % 64 == 0 {
            poll(cancelled)?;
        }
        slopes[k] = [[0.0; 5]; 3];
        if reconstruction == Reconstruction::FirstOrder {
            continue;
        }
        let coordinates = grid.coordinates(k);
        let center = primitive[k].array();
        for axis in 0..3 {
            let lower = neighbour(primitive, grid, boundaries, coordinates, axis, false);
            let upper = neighbour(primitive, grid, boundaries, coordinates, axis, true);
            slopes[k][axis] =
                std::array::from_fn(|d| mc(center[d] - lower[d], upper[d] - center[d]));
        }
        // One cell factor limits all six face states without clipping an average.
        // The final flux constructor still checks full conservative precision.
        let mut theta = 1.0;
        for _ in 0..32 {
            let valid = (0..3).all(|axis| {
                [-0.5, 0.5].into_iter().all(|sign| {
                    let p: [f64; 5] =
                        std::array::from_fn(|d| center[d] + sign * theta * slopes[k][axis][d]);
                    let temperature = p[4] / (gas.gas_constant() * p[0]);
                    p.into_iter().all(f64::is_finite)
                        && positive(p[0])
                        && positive(p[4])
                        && temperature >= gas.domain().minimum_temperature
                        && temperature <= gas.domain().maximum_temperature
                })
            });
            if valid {
                break;
            }
            theta *= 0.5;
        }
        if theta < 1.0 {
            stats.limited += 1;
            if theta <= 2.0_f64.powi(-32) {
                theta = 0.0;
            }
            for slope in &mut slopes[k] {
                for value in slope {
                    *value *= theta;
                }
            }
        }
    }
    // Every axis uses a contiguous, x-fast face lattice. Internal neighbours
    // consume the exact same stored flux; periodic endpoints share its bits.
    for axis in 0..3 {
        let flux = &mut fluxes[axis];
        let count = grid.shape[axis];
        let pair = boundaries.axis(axis);
        let [t1, t2] = tangents(axis);
        let mut face_shape = grid.shape;
        face_shape[axis] += 1;
        for index in 0..flux.len() {
            if index % 64 == 0 {
                poll(cancelled)?;
            }
            let mut c = [
                index % face_shape[0],
                (index / face_shape[0]) % face_shape[1],
                index / (face_shape[0] * face_shape[1]),
            ];
            let face = c[axis];
            if face == count && pair == AxisBoundaries::Periodic {
                c[axis] = 0;
                flux[index] = flux[grid.cell_face(axis, c)];
                continue;
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
            let f = if face > 0 && face < count || pair == AxisBoundaries::Periodic {
                let mut left = c;
                let mut right = c;
                if face == 0 {
                    left[axis] = count - 1;
                    right[axis] = 0;
                } else {
                    left[axis] -= 1;
                }
                let left = grid.index(left);
                let right = grid.index(right);
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
            } else {
                let AxisBoundaries::Ends { lower, upper } = pair else {
                    unreachable!()
                };
                let at_lower = face == 0;
                let boundary = if at_lower { lower } else { upper };
                let sign = if at_lower { -1.0 } else { 1.0 };
                c[axis] = if at_lower { 0 } else { count - 1 };
                let k = grid.index(c);
                match (first_order, boundary) {
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
                }
            };
            stats.fallback += usize::from(f.used_hll_fallback);
            flux[index] = f;
        }
    }
    let strides = grid.face_strides();
    for k in 0..states.len() {
        if k % 64 == 0 {
            poll(cancelled)?;
        }
        let c = grid.coordinates(k);
        let mut rate = 0.0;
        for axis in 0..3 {
            let low = grid.cell_face(axis, c);
            rate += fluxes[axis][low]
                .signal_speed
                .max(fluxes[axis][low + strides[axis]].signal_speed)
                / grid.spacing[axis];
        }
        if !rate.is_finite() {
            return Err(Error::Invalid("unrepresentable three-axis face CFL"));
        }
        stats.rate = stats.rate.max(rate);
    }
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
    let position = c[axis];
    let count = grid.shape[axis];
    let at_edge = if upper {
        position + 1 == count
    } else {
        position == 0
    };
    let mut next = c;
    if !at_edge {
        next[axis] = if upper { position + 1 } else { position - 1 };
        p[grid.index(next)].array()
    } else if boundaries.axis(axis) == AxisBoundaries::Periodic {
        next[axis] = if upper { 0 } else { count - 1 };
        p[grid.index(next)].array()
    } else {
        let mut q = p[grid.index(c)].array();
        let AxisBoundaries::Ends { lower, upper: high } = boundaries.axis(axis) else {
            unreachable!()
        };
        if (if upper { high } else { lower }) == BoundaryCondition::Reflecting {
            q[axis + 1] = -q[axis + 1];
        }
        q
    }
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
    let strides = grid.face_strides();
    let scale = grid.spacing.map(|h| dt / h);
    for k in 0..output.len() {
        if k % 64 == 0 {
            poll(cancelled)?;
        }
        let c = grid.coordinates(k);
        let mut delta = [0.0; 5];
        for axis in 0..3 {
            let low = grid.cell_face(axis, c);
            let lower = flux_array(fluxes[axis][low], axis);
            let upper = flux_array(fluxes[axis][low + strides[axis]], axis);
            for d in 0..5 {
                delta[d] += scale[axis] * (lower[d] - upper[d]);
            }
        }
        let old = original[k].array();
        let current = stage[k].array();
        output[k] = State::from_array(std::array::from_fn(|d| {
            if final_stage {
                0.5 * old[d] + 0.5 * current[d] + 0.5 * delta[d]
            } else {
                current[d] + delta[d]
            }
        }));
        output[k].primitive(gas, &mut *cancelled)?;
    }
    Ok(())
}

pub(super) fn capture_boundaries(
    grid: &CartesianGrid,
    fluxes: &[Vec<SpatialFaceFlux>; 3],
    duration: f64,
    output: &mut [BoundaryImpulse],
    add: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let face = if side.sign() < 0.0 {
            0
        } else {
            grid.shape[axis]
        };
        let scale = side.sign() * duration * grid.face_area(side);
        for (line, target) in output[grid.boundary_range(side)].iter_mut().enumerate() {
            if line % 64 == 0 {
                poll(cancelled)?;
            }
            let flux = flux_array(fluxes[axis][grid.boundary_face(axis, face, line)], axis);
            let previous = if add { target.array() } else { [0.0; 5] };
            *target = Totals::from_array(std::array::from_fn(|d| previous[d] + scale * flux[d]));
        }
    }
    Ok(())
}
