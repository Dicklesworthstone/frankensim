//! Vorticity and the inviscid material-rate budget of accepted Euler cells.
//!
//! For smooth flow, with G_ij = du_i/dx_j and omega = curl(u),
//! D(omega)/Dt = G*omega - omega*div(u) + grad(rho) x grad(p)/rho^2.
//! These are three distinct contributions: stretching, compressible dilation,
//! and baroclinicity. Conservative pressure fluxes already generate the last
//! term; this module observes them and never applies another momentum source.
//!
//! Gradients use physical axis widths and second-order centered differences.
//! Periodic axes wrap; a single periodic cell denotes homogeneous extrusion.
//! At nonperiodic ends a quadratic through three interior cell centers supplies
//! a second-order one-sided derivative. This is an interior reconstruction,
//! not a wall/reservoir ghost prescription or an estimate of boundary traction.
//! Two-cell axes, and single-cell nonperiodic axes, cannot resolve this stencil
//! and are refused. No absent neighbor is silently replaced by zero.
//!
//! The budget is Estimated for smooth, resolved cell-average fields. It is not
//! an exact discrete circulation identity or a shock derivative. Viscous torque,
//! nonconservative body-force curl and unresolved turbulent terms are excluded;
//! the material rate is not the Eulerian partial time derivative.

use super::{
    AxisBoundaries, Boundaries3d, CartesianGrid, Error, Euler3d, Primitive, allocated, poll,
};
use crate::compressible::cartesian::Failure;

/// Local physical curl and separately observable inviscid production terms.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct VorticityBudget3d {
    /// curl(u) [1/s], with right-handed x/y/z components.
    pub vorticity: [f64; 3],
    /// div(u) [1/s]; expansion is positive.
    pub divergence: f64,
    /// (omega dot grad)u = G*omega [1/s^2].
    pub stretching: [f64; 3],
    /// -omega*div(u) [1/s^2], including compressible expansion/compression.
    pub dilation: [f64; 3],
    /// grad(rho) cross grad(p) / rho^2 [1/s^2].
    pub baroclinic: [f64; 3],
}

impl VorticityBudget3d {
    /// Sum of the three resolved inviscid contributions to D(omega)/Dt [1/s^2].
    ///
    /// Advection is already on the left-hand side of this material derivative.
    /// Molecular, applied-force and SGS contributions are not included.
    pub fn inviscid_material_rate(self) -> [f64; 3] {
        std::array::from_fn(|a| self.stretching[a] + self.dilation[a] + self.baroclinic[a])
    }
}

/// Reusable, explicitly budgeted scratch for source recovery and native curl.
///
/// Output is private and returned only after a complete successful evaluation.
/// A failed evaluation may change scratch, but cannot expose partial results
/// or change the flow's cells, time, accepted loads or step report.
pub struct VorticityWorkspace3d {
    pub(super) grid: CartesianGrid,
    pub(super) primitive: Vec<Primitive>,
    pub(super) output: Vec<VorticityBudget3d>,
}

impl VorticityWorkspace3d {
    /// Allocated numerical capacity in bytes, excluding fixed metadata and
    /// allocator overhead. Evaluation allocates no additional field buffers.
    pub fn workspace_bytes(&self) -> usize {
        self.primitive.capacity() * core::mem::size_of::<Primitive>()
            + self.output.capacity() * core::mem::size_of::<VorticityBudget3d>()
    }
}

impl Euler3d {
    /// Required numerical storage for one primitive and one budget per cell.
    pub fn required_vorticity_workspace_bytes(&self) -> Result<usize, Error> {
        self.states
            .len()
            .checked_mul(
                core::mem::size_of::<Primitive>() + core::mem::size_of::<VorticityBudget3d>(),
            )
            .ok_or(Error::Invalid("vorticity workspace size overflow"))
    }

    /// Allocate native vorticity scratch within an explicit numerical-byte cap.
    ///
    /// Each resolved axis needs at least three cells; a singleton periodic axis
    /// is admitted as homogeneous extrusion. Nonperiodic derivatives use only
    /// the interior, independently of exterior reservoir or wall prescriptions.
    pub fn vorticity_workspace(&self, maximum_bytes: usize) -> Result<VorticityWorkspace3d, Error> {
        admit_axes::<Error>(&self.grid, self.boundaries)?;
        if self.required_vorticity_workspace_bytes()? > maximum_bytes {
            return Err(Error::Invalid(
                "vorticity workspace exceeds its byte budget",
            ));
        }
        let work = VorticityWorkspace3d {
            grid: self.grid,
            primitive: allocated(self.states.len(), Primitive::default())?,
            output: allocated(self.states.len(), VorticityBudget3d::default())?,
        };
        if work.workspace_bytes() > maximum_bytes {
            return Err(Error::Invalid(
                "vorticity workspace capacity exceeds its byte budget",
            ));
        }
        Ok(work)
    }

    /// Evaluate the accepted field in x-fast cell order, without evolving it.
    ///
    /// Pressure and velocity are recovered once per cell from this owner's
    /// original gas EOS and accepted conservative state, never from stale RK
    /// scratch or a substituted constant-gamma law. Cancellation is polled
    /// during source recovery, at most every 64 gradient cells, and immediately
    /// before publication. Repeated calls reuse the same numerical capacity.
    pub fn vorticity<'a>(
        &self,
        work: &'a mut VorticityWorkspace3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&'a [VorticityBudget3d], Error> {
        poll(&mut cancelled)?;
        admit_axes::<Error>(&self.grid, self.boundaries)?;
        if work.grid != self.grid {
            return Err(Error::Invalid(
                "vorticity workspace and physical grid differ",
            ));
        }
        for (cell, (primitive, state)) in work.primitive.iter_mut().zip(&self.states).enumerate() {
            if cell % 64 == 0 {
                poll(&mut cancelled)?;
            }
            *primitive = state.primitive(&self.gas, &mut cancelled)?;
        }
        for (cell, output) in work.output.iter_mut().enumerate() {
            if cell % 64 == 0 {
                poll(&mut cancelled)?;
            }
            *output = budget::<Error>(&self.grid, self.boundaries, &work.primitive, cell)?;
        }
        poll(&mut cancelled)?;
        Ok(&work.output)
    }
}

pub(super) fn admit_axes<E: Failure>(
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
) -> Result<(), E> {
    for axis in 0..3 {
        let count = grid.shape[axis];
        if count < 3 && !(count == 1 && boundaries.axis(axis) == AxisBoundaries::Periodic) {
            return Err(E::invalid(
                "vorticity needs three cells per resolved axis or one periodic cell",
            ));
        }
    }
    Ok(())
}

/// Shared physical derivative for cell-centered primitive or acceleration fields.
pub(super) fn derivative<const N: usize>(
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    cell: usize,
    axis: usize,
    value: impl Fn(usize) -> [f64; N],
) -> [f64; N] {
    let c = grid.coordinates(cell);
    let count = grid.shape[axis];
    if count == 1 {
        return [0.0; N];
    }
    let at = |position| {
        let mut other = c;
        other[axis] = position;
        value(grid.index(other))
    };
    let h = grid.spacing[axis];
    let position = c[axis];
    if boundaries.axis(axis) == AxisBoundaries::Periodic || (position > 0 && position + 1 < count) {
        let lower = at(if position == 0 {
            count - 1
        } else {
            position - 1
        });
        let upper = at(if position + 1 == count {
            0
        } else {
            position + 1
        });
        return std::array::from_fn(|d| (0.5 * (upper[d] - lower[d])) / h);
    }
    let center = value(cell);
    let (near, far, sign) = if position == 0 {
        (at(1), at(2), 1.0)
    } else {
        (at(count - 2), at(count - 3), -1.0)
    };
    // Differences remove uniform offsets before the one-sided combination.
    std::array::from_fn(|d| {
        sign * (1.5 * ((near[d] - center[d]) / h) - 0.5 * ((far[d] - near[d]) / h))
    })
}

pub(super) fn budget<E: Failure>(
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    primitive: &[Primitive],
    cell: usize,
) -> Result<VorticityBudget3d, E> {
    let mut gradient = [[0.0; 3]; 5];
    for axis in 0..3 {
        let column = derivative(grid, boundaries, cell, axis, |other| {
            primitive[other].array()
        });
        for (component, row) in gradient.iter_mut().enumerate() {
            row[axis] = column[component];
        }
    }
    if !gradient.iter().flatten().all(|value| value.is_finite()) {
        return Err(E::invalid("unrepresentable physical vorticity gradient"));
    }
    let omega = [
        gradient[3][1] - gradient[2][2],
        gradient[1][2] - gradient[3][0],
        gradient[2][0] - gradient[1][1],
    ];
    let divergence = gradient[1][0] + gradient[2][1] + gradient[3][2];
    let stretching = std::array::from_fn(|i| (0..3).map(|j| gradient[i + 1][j] * omega[j]).sum());
    // Normalize both gradients before the cross product; rho^2 need never
    // overflow or underflow merely because the caller's density scale changes.
    let density = primitive[cell].density;
    let r = gradient[0].map(|value| value / density);
    let p = gradient[4].map(|value| value / density);
    let result = VorticityBudget3d {
        vorticity: omega,
        divergence,
        stretching,
        dilation: omega.map(|value| -value * divergence),
        baroclinic: [
            r[1] * p[2] - r[2] * p[1],
            r[2] * p[0] - r[0] * p[2],
            r[0] * p[1] - r[1] * p[0],
        ],
    };
    if !result
        .vorticity
        .into_iter()
        .chain([result.divergence])
        .chain(result.stretching)
        .chain(result.dilation)
        .chain(result.baroclinic)
        .chain(result.inviscid_material_rate())
        .all(f64::is_finite)
    {
        return Err(E::invalid("unrepresentable physical vorticity budget"));
    }
    Ok(result)
}

#[cfg(test)]
mod tests;
