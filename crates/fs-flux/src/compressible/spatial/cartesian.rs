//! Shared three-dimensional topology and conservative spatial arithmetic.
//! Concrete source/state views monomorphize; this module owns no gas physics.

pub(super) use super::super::cartesian::{Face, FaceKind};
use super::super::{
    Reconstruction,
    cartesian::{Failure, mc, poll},
};
use super::{AxisBoundaries, Boundaries3d, BoundaryCondition, BoundarySide, CartesianGrid};

/// Original cell and whether a stationary-wall ghost reflects normal velocity.
pub(super) fn neighbour(
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    cell: usize,
    axis: usize,
    upper: bool,
) -> (usize, bool) {
    let mut c = grid.coordinates(cell);
    let position = c[axis];
    let count = grid.shape[axis];
    let edge = if upper {
        position + 1 == count
    } else {
        position == 0
    };
    if !edge {
        c[axis] = if upper { position + 1 } else { position - 1 };
        return (grid.index(c), false);
    }
    match boundaries.axis(axis) {
        AxisBoundaries::Periodic => {
            c[axis] = if upper { 0 } else { count - 1 };
            (grid.index(c), false)
        }
        AxisBoundaries::Ends { lower, upper: high } => (
            cell,
            (if upper { high } else { lower }) == BoundaryCondition::Reflecting,
        ),
    }
}

/// One source-admitted contraction for all components on all six traces.
#[allow(clippy::too_many_arguments)]
pub(super) fn slopes<E: Failure>(
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    width: usize,
    velocity_components: [usize; 3],
    reconstruction: Reconstruction,
    value: impl Fn(usize, usize) -> f64,
    output: &mut [f64],
    mut admits: impl FnMut(usize, usize, f64, &[f64]) -> bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<usize, E> {
    let mut limited = 0;
    for cell in 0..grid.len() {
        poll::<E>(cancelled)?;
        let start = 3 * cell * width;
        output[start..start + 3 * width].fill(0.0);
        if reconstruction == Reconstruction::FirstOrder {
            continue;
        }
        for axis in 0..3 {
            let (lo, reflect_lo) = neighbour(grid, boundaries, cell, axis, false);
            let (hi, reflect_hi) = neighbour(grid, boundaries, cell, axis, true);
            let offset = start + axis * width;
            for d in 0..width {
                let lower = value(lo, d)
                    * if reflect_lo && d == velocity_components[axis] {
                        -1.0
                    } else {
                        1.0
                    };
                let upper = value(hi, d)
                    * if reflect_hi && d == velocity_components[axis] {
                        -1.0
                    } else {
                        1.0
                    };
                output[offset + d] = mc(value(cell, d) - lower, upper - value(cell, d));
            }
        }
        let mut theta = 1.0;
        for _ in 0..32 {
            poll::<E>(cancelled)?;
            if (0..3).all(|axis| {
                let offset = start + axis * width;
                [-0.5, 0.5]
                    .into_iter()
                    .all(|sign| admits(cell, axis, sign * theta, &output[offset..offset + width]))
            }) {
                break;
            }
            theta *= 0.5;
        }
        if theta < 1.0 {
            limited += 1;
            if theta <= 2.0_f64.powi(-32) {
                theta = 0.0;
            }
            for s in &mut output[start..start + 3 * width] {
                *s *= theta;
            }
        }
    }
    Ok(limited)
}

/// X-fast face lattices; periodic upper endpoints copy their lower partner.
pub(super) fn faces<E: Failure, C: FnMut() -> bool>(
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    cancelled: &mut C,
    mut visit: impl FnMut(Face, &mut C) -> Result<(), E>,
) -> Result<(), E> {
    for axis in 0..3 {
        let count = grid.shape[axis];
        let pair = boundaries.axis(axis);
        let mut shape = grid.shape;
        shape[axis] += 1;
        for index in 0..grid.face_len(axis) {
            poll::<E>(cancelled)?;
            let mut c = [
                index % shape[0],
                (index / shape[0]) % shape[1],
                index / (shape[0] * shape[1]),
            ];
            let face = c[axis];
            let kind = if face == count && pair == AxisBoundaries::Periodic {
                c[axis] = 0;
                FaceKind::PeriodicCopy {
                    first: grid.cell_face(axis, c),
                }
            } else if face > 0 && face < count || pair == AxisBoundaries::Periodic {
                let mut left = c;
                let mut right = c;
                if face == 0 {
                    left[axis] = count - 1;
                    right[axis] = 0;
                } else {
                    left[axis] -= 1;
                }
                FaceKind::Interior {
                    left: grid.index(left),
                    right: grid.index(right),
                }
            } else {
                let AxisBoundaries::Ends { lower, upper } = pair else {
                    unreachable!()
                };
                let at_lower = face == 0;
                c[axis] = if at_lower { 0 } else { count - 1 };
                FaceKind::Boundary {
                    cell: grid.index(c),
                    sign: if at_lower { -1.0 } else { 1.0 },
                    condition: if at_lower { lower } else { upper },
                }
            };
            visit(Face { axis, index, kind }, cancelled)?;
        }
    }
    Ok(())
}

pub(super) fn rate<E: Failure>(
    grid: &CartesianGrid,
    speed: impl Fn(usize, usize) -> f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, E> {
    let strides = grid.face_strides();
    let mut maximum: f64 = 0.0;
    for k in 0..grid.len() {
        poll::<E>(cancelled)?;
        let c = grid.coordinates(k);
        let mut rate = 0.0;
        for axis in 0..3 {
            let low = grid.cell_face(axis, c);
            rate += speed(axis, low).max(speed(axis, low + strides[axis])) / grid.spacing[axis];
        }
        if !rate.is_finite() || rate < 0.0 {
            return Err(E::invalid("volumetric face CFL overflow"));
        }
        maximum = maximum.max(rate);
    }
    Ok(maximum)
}

/// The same conservative divergence and SSP arithmetic for every state layout.
#[allow(clippy::too_many_arguments)]
pub(super) fn update<const MAX: usize, E: Failure, C: FnMut() -> bool>(
    grid: &CartesianGrid,
    width: usize,
    dt: f64,
    final_stage: bool,
    original: impl Fn(usize, usize) -> f64,
    stage: impl Fn(usize, usize) -> f64,
    flux: impl Fn(usize, usize, usize) -> f64,
    cancelled: &mut C,
    mut write: impl FnMut(usize, &[f64], &mut C) -> Result<(), E>,
) -> Result<(), E> {
    let strides = grid.face_strides();
    let scale = grid.spacing.map(|h| dt / h);
    let mut row = [0.0; MAX];
    for k in 0..grid.len() {
        poll::<E>(cancelled)?;
        let c = grid.coordinates(k);
        row[..width].fill(0.0);
        for axis in 0..3 {
            let low = grid.cell_face(axis, c);
            for (d, delta) in row[..width].iter_mut().enumerate() {
                *delta += scale[axis] * (flux(axis, low, d) - flux(axis, low + strides[axis], d));
            }
        }
        for (d, out) in row[..width].iter_mut().enumerate() {
            *out = if final_stage {
                0.5 * original(k, d) + 0.5 * stage(k, d) + 0.5 * *out
            } else {
                stage(k, d) + *out
            };
        }
        write(k, &row[..width], cancelled)?;
    }
    Ok(())
}

pub(super) fn capture<E: Failure>(
    grid: &CartesianGrid,
    width: usize,
    duration: f64,
    add: bool,
    flux: impl Fn(usize, usize, usize) -> f64,
    mut accumulate: impl FnMut(usize, usize, f64, bool),
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), E> {
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let face = if side.sign() < 0.0 {
            0
        } else {
            grid.shape[axis]
        };
        let scale = side.sign() * duration * grid.face_area(side);
        for (line, boundary) in grid.boundary_range(side).enumerate() {
            poll::<E>(cancelled)?;
            let index = grid.boundary_face(axis, face, line);
            for d in 0..width {
                accumulate(boundary, d, scale * flux(axis, index, d), add);
            }
        }
    }
    Ok(())
}
