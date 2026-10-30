//! Shared Cartesian discretization; thermodynamics and Riemann physics stay
//! with the source owners. Concrete frozen/species views are monomorphized.

use super::{
    AxisBoundaries, Boundaries2d, BoundaryCondition, CartesianGrid, Reconstruction, SIDES,
    StepControl, positive,
};

pub(super) trait Failure: Sized {
    fn invalid(reason: &'static str) -> Self;
    fn cancelled() -> Self;
    fn cfl(observed: f64, maximum: f64) -> Self;
    fn retries(attempts: usize, last: Self) -> Self;
    fn is_cancelled(&self) -> bool;
    fn is_cfl(&self) -> bool;
}
impl Failure for super::Error {
    fn invalid(reason: &'static str) -> Self {
        Self::Invalid(reason)
    }
    fn cancelled() -> Self {
        Self::Cancelled
    }
    fn cfl(observed: f64, maximum: f64) -> Self {
        Self::CflExceeded { observed, maximum }
    }
    fn retries(attempts: usize, last: Self) -> Self {
        Self::StepRetries {
            attempts,
            last: Box::new(last),
        }
    }
    fn is_cancelled(&self) -> bool {
        matches!(self, Self::Cancelled)
    }
    fn is_cfl(&self) -> bool {
        matches!(self, Self::CflExceeded { .. })
    }
}
impl Failure for fs_gas::multicomponent::Error {
    fn invalid(reason: &'static str) -> Self {
        Self::Invalid(reason)
    }
    fn cancelled() -> Self {
        Self::Cancelled
    }
    fn cfl(observed: f64, maximum: f64) -> Self {
        Self::CflExceeded { observed, maximum }
    }
    fn retries(attempts: usize, last: Self) -> Self {
        Self::StepRetries {
            attempts,
            last: Box::new(last),
        }
    }
    fn is_cancelled(&self) -> bool {
        matches!(self, Self::Cancelled)
    }
    fn is_cfl(&self) -> bool {
        matches!(self, Self::CflExceeded { .. })
    }
}
pub(super) fn poll<E: Failure>(cancelled: &mut impl FnMut() -> bool) -> Result<(), E> {
    if cancelled() {
        Err(E::cancelled())
    } else {
        Ok(())
    }
}
pub(super) fn check_cfl<E: Failure>(observed: f64, maximum: f64) -> Result<(), E> {
    if !observed.is_finite() || observed > maximum * (1.0 + 16.0 * f64::EPSILON) {
        Err(E::cfl(observed, maximum))
    } else {
        Ok(())
    }
}

#[derive(Clone, Copy, Default)]
pub(super) struct StageStats {
    pub rate: f64,
    pub fallback: usize,
    pub limited: usize,
}
pub(super) struct AcceptedStep {
    pub dt: f64,
    pub halvings: usize,
    pub reconstruction: Reconstruction,
    pub stats: StageStats,
}

/// Identical CFL/reconstruction retry policy for both conservative layouts.
#[allow(clippy::too_many_arguments)]
pub(super) fn advance<E: Failure, C: FnMut() -> bool>(
    maximum_duration: f64,
    time: f64,
    initial_rate: f64,
    controls: StepControl,
    cancelled: &mut C,
    mut attempt: impl FnMut(f64, Reconstruction, &mut C) -> Result<StageStats, E>,
) -> Result<AcceptedStep, E> {
    poll::<E>(cancelled)?;
    if !positive(maximum_duration)
        || !positive(controls.cfl)
        || controls.cfl > 0.45
        || controls.maximum_halvings > 32
    {
        return Err(E::invalid("Cartesian duration, CFL or retry controls"));
    }
    let mut dt = maximum_duration.min(controls.cfl / initial_rate);
    let mut reconstruction = controls.reconstruction;
    for halvings in 0..=controls.maximum_halvings {
        poll::<E>(cancelled)?;
        if !positive(dt) || !time.is_finite() || time + dt <= time || !(time + dt).is_finite() {
            return Err(E::invalid("unrepresentable Cartesian time increment"));
        }
        match attempt(dt, reconstruction, cancelled) {
            Ok(stats) => {
                return Ok(AcceptedStep {
                    dt,
                    halvings,
                    reconstruction,
                    stats,
                });
            }
            Err(error) if error.is_cancelled() => return Err(error),
            Err(error) => {
                if halvings == controls.maximum_halvings {
                    return Err(E::retries(halvings + 1, error));
                }
                if !error.is_cfl() {
                    reconstruction = Reconstruction::FirstOrder;
                }
                dt *= 0.5;
            }
        }
    }
    unreachable!("bounded retries return on the last attempt")
}

/// Stage false writes the forward state; stage true writes its SSP combination.
pub(super) fn rk2<E, C: FnMut() -> bool>(
    cancelled: &mut C,
    mut stage: impl FnMut(bool, &mut C) -> Result<StageStats, E>,
) -> Result<StageStats, E> {
    let first = stage(false, cancelled)?;
    let second = stage(true, cancelled)?;
    Ok(StageStats {
        rate: first.rate.max(second.rate),
        fallback: first.fallback + second.fallback,
        limited: first.limited + second.limited,
    })
}

/// Cell index and whether the normal velocity must be reflected in a ghost.
pub(super) fn neighbour(
    grid: &CartesianGrid,
    boundaries: Boundaries2d,
    cell: usize,
    axis: usize,
    upper: bool,
) -> (usize, bool) {
    let [nx, _] = grid.shape;
    let position = if axis == 0 { cell % nx } else { cell / nx };
    let count = grid.shape[axis];
    let edge = if upper {
        position + 1 == count
    } else {
        position == 0
    };
    if !edge {
        let stride = if axis == 0 { 1 } else { nx };
        return (if upper { cell + stride } else { cell - stride }, false);
    }
    match boundaries.axis(axis) {
        AxisBoundaries::Periodic => {
            let stride = if axis == 0 { 1 } else { nx };
            (
                if upper {
                    cell - (count - 1) * stride
                } else {
                    cell + (count - 1) * stride
                },
                false,
            )
        }
        AxisBoundaries::Ends { lower, upper: high } => (
            cell,
            (if upper { high } else { lower }) == BoundaryCondition::Reflecting,
        ),
    }
}

/// MC reconstruction with one source-admitted theta for every component/axis.
/// A source may limit a coupled leading block before scalar velocity/pressure.
#[allow(clippy::too_many_arguments)]
pub(super) fn slopes<E: Failure>(
    grid: &CartesianGrid,
    boundaries: Boundaries2d,
    width: usize,
    velocity_components: [usize; 2],
    reconstruction: Reconstruction,
    value: impl Fn(usize, usize) -> f64,
    mut coupled: impl FnMut(usize, usize, usize, &mut [f64]) -> Result<usize, E>,
    output: &mut [f64],
    mut admits: impl FnMut(usize, usize, f64, &[f64]) -> bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<usize, E> {
    output.fill(0.0);
    if reconstruction == Reconstruction::FirstOrder {
        return Ok(0);
    }
    let mut limited = 0;
    for cell in 0..grid.len() {
        poll::<E>(cancelled)?;
        for axis in 0..2 {
            let (lo, reflect_lo) = neighbour(grid, boundaries, cell, axis, false);
            let (hi, reflect_hi) = neighbour(grid, boundaries, cell, axis, true);
            let offset = (2 * cell + axis) * width;
            let leading = coupled(cell, lo, hi, &mut output[offset..offset + width])?;
            for d in leading..width {
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
            if (0..2).all(|axis| {
                let offset = (2 * cell + axis) * width;
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
            for s in &mut output[2 * cell * width..2 * (cell + 1) * width] {
                *s *= theta;
            }
        }
    }
    Ok(limited)
}
pub(super) fn mc(left: f64, right: f64) -> f64 {
    if left > 0.0 && right > 0.0 {
        (0.5 * left + 0.5 * right).min(2.0 * left).min(2.0 * right)
    } else if left < 0.0 && right < 0.0 {
        (0.5 * left + 0.5 * right).max(2.0 * left).max(2.0 * right)
    } else {
        0.0
    }
}

#[derive(Clone, Copy)]
pub(super) enum FaceKind {
    Interior {
        left: usize,
        right: usize,
    },
    Boundary {
        cell: usize,
        sign: f64,
        condition: BoundaryCondition,
    },
    PeriodicCopy {
        first: usize,
    },
}
#[derive(Clone, Copy)]
pub(super) struct Face {
    pub axis: usize,
    pub index: usize,
    pub kind: FaceKind,
}
pub(super) fn face_index(grid: &CartesianGrid, axis: usize, face: usize, line: usize) -> usize {
    if axis == 0 {
        line * (grid.shape[0] + 1) + face
    } else {
        face * grid.shape[0] + line
    }
}
pub(super) fn faces<E: Failure, C: FnMut() -> bool>(
    grid: &CartesianGrid,
    boundaries: Boundaries2d,
    cancelled: &mut C,
    mut visit: impl FnMut(Face, &mut C) -> Result<(), E>,
) -> Result<(), E> {
    for axis in 0..2 {
        let count = grid.shape[axis];
        let pair = boundaries.axis(axis);
        for line in 0..grid.shape[1 - axis] {
            let cell = |position| {
                if axis == 0 {
                    line * grid.shape[0] + position
                } else {
                    position * grid.shape[0] + line
                }
            };
            for f in 0..=count {
                poll::<E>(cancelled)?;
                let kind = if f == count && pair == AxisBoundaries::Periodic {
                    FaceKind::PeriodicCopy {
                        first: face_index(grid, axis, 0, line),
                    }
                } else if f > 0 && f < count {
                    FaceKind::Interior {
                        left: cell(f - 1),
                        right: cell(f),
                    }
                } else if pair == AxisBoundaries::Periodic {
                    FaceKind::Interior {
                        left: cell(count - 1),
                        right: cell(0),
                    }
                } else {
                    let AxisBoundaries::Ends { lower, upper } = pair else {
                        unreachable!()
                    };
                    FaceKind::Boundary {
                        cell: cell(if f == 0 { 0 } else { count - 1 }),
                        sign: if f == 0 { -1.0 } else { 1.0 },
                        condition: if f == 0 { lower } else { upper },
                    }
                };
                visit(
                    Face {
                        axis,
                        index: face_index(grid, axis, f, line),
                        kind,
                    },
                    cancelled,
                )?;
            }
        }
    }
    Ok(())
}
pub(super) fn rate<E: Failure>(
    grid: &CartesianGrid,
    speed: impl Fn(usize, usize) -> f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, E> {
    let [nx, ny] = grid.shape;
    let mut maximum: f64 = 0.0;
    for j in 0..ny {
        for i in 0..nx {
            poll::<E>(cancelled)?;
            let x = j * (nx + 1) + i;
            let y = j * nx + i;
            let rate = speed(0, x).max(speed(0, x + 1)) / grid.spacing[0]
                + speed(1, y).max(speed(1, y + nx)) / grid.spacing[1];
            if !rate.is_finite() || rate < 0.0 {
                return Err(E::invalid("Cartesian face rate overflow"));
            }
            maximum = maximum.max(rate);
        }
    }
    Ok(maximum)
}

/// The same conservative divergence and SSP arithmetic order for all components.
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
    let [nx, ny] = grid.shape;
    let mut row = [0.0; MAX];
    for j in 0..ny {
        for i in 0..nx {
            poll::<E>(cancelled)?;
            let k = j * nx + i;
            let x = j * (nx + 1) + i;
            for (d, out) in row[..width].iter_mut().enumerate() {
                let delta = dt / grid.spacing[0] * (flux(0, x, d) - flux(0, x + 1, d))
                    + dt / grid.spacing[1] * (flux(1, k, d) - flux(1, k + nx, d));
                *out = if final_stage {
                    0.5 * original(k, d) + 0.5 * stage(k, d) + 0.5 * delta
                } else {
                    stage(k, d) + delta
                };
            }
            write(k, &row[..width], cancelled)?;
        }
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
    for side in SIDES {
        let axis = side.axis();
        let f = if side.sign() < 0.0 {
            0
        } else {
            grid.shape[axis]
        };
        let scale = side.sign() * duration * grid.face_area(side);
        for (line, boundary) in grid.boundary_range(side).enumerate() {
            poll::<E>(cancelled)?;
            let index = face_index(grid, axis, f, line);
            for d in 0..width {
                accumulate(boundary, d, scale * flux(axis, index, d), add);
            }
        }
    }
    Ok(())
}
