//! Moment-preserving redistribution of planar vortex circulation.
//!
//! The tensor-product M4' kernel deposits each original circulation onto at
//! most 4x4 nodes of an explicit Cartesian grid. Its discrete zeroth, first
//! and second moments reproduce those of each source particle, retaining
//! circulation and linear/angular vortex-impulse moments without a global
//! correction. See Siemaszko et al., Int. J. Numer. Methods Fluids (2025),
//! equations 34–35, <https://doi.org/10.1002/fld.70002>.
//!
//! The interpolation has signed outer lobes: positivity, maximum principles
//! and enstrophy contraction are not claimed. The velocity core is unchanged.
//! Full nonzero stencil support must fit the authored grid; walls, wrapping,
//! clipped tails and particle-strength thresholds are not inferred.

use crate::VortexParticle;

/// Explicit axis-aligned Cartesian quadrature nodes in metres.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CartesianGrid2 {
    origin: [f64; 2],
    spacing: [f64; 2],
    shape: [usize; 2],
    nodes: usize,
    area: f64,
}

impl CartesianGrid2 {
    /// Admit a grid with at least four nodes per axis, finite distinct node
    /// positions, and a representable positive quadrature area `dx*dy`.
    /// The origin is node [0,0], and x varies fastest in the output.
    pub fn new(
        origin_m: [f64; 2],
        spacing_m: [f64; 2],
        shape: [usize; 2],
    ) -> Result<Self, RemeshError> {
        let nodes = shape[0]
            .checked_mul(shape[1])
            .ok_or(RemeshError::Invalid("grid size overflow"))?;
        let area = spacing_m[0] * spacing_m[1];
        if !area.is_finite() || area <= 0.0 {
            return Err(RemeshError::Invalid(
                "positive representable grid area required",
            ));
        }
        for axis in 0..2 {
            if !origin_m[axis].is_finite()
                || !spacing_m[axis].is_finite()
                || spacing_m[axis] <= 0.0
                || shape[axis] < 4
                || shape[axis] as u128 > (1_u128 << 52)
            {
                return Err(RemeshError::Invalid(
                    "finite grid origin, positive spacing and at least four nodes required",
                ));
            }
            let span = spacing_m[axis] * (shape[axis] - 1) as f64;
            let end = origin_m[axis] + span;
            let magnitude = origin_m[axis].abs().max(end.abs()).max(span);
            if !span.is_finite()
                || !end.is_finite()
                || spacing_m[axis] <= 2.0 * f64::EPSILON * magnitude
            {
                return Err(RemeshError::Invalid(
                    "grid nodes are not representable at this origin and spacing",
                ));
            }
        }
        Ok(Self {
            origin: origin_m,
            spacing: spacing_m,
            shape,
            nodes,
            area,
        })
    }

    /// Original first grid node [m].
    pub const fn origin_m(self) -> [f64; 2] {
        self.origin
    }
    /// Independent x/y node spacing [m].
    pub const fn spacing_m(self) -> [f64; 2] {
        self.spacing
    }
    /// Number of x/y nodes.
    pub const fn shape(self) -> [usize; 2] {
        self.shape
    }
    /// Total output node count, including zero-strength diffusion support.
    pub const fn node_count(self) -> usize {
        self.nodes
    }
    /// Uniform quadrature area [m²], independent of any velocity core.
    pub const fn quadrature_area_m2(self) -> f64 {
        self.area
    }

    fn coordinate(self, axis: usize, index: usize) -> f64 {
        self.origin[axis] + self.spacing[axis] * index as f64
    }

    fn stencil(
        self,
        position: f64,
        axis: usize,
        particle: usize,
    ) -> Result<([usize; 4], [f64; 4]), RemeshError> {
        let mut coordinate = (position - self.origin[axis]) / self.spacing[axis];
        let end = (self.shape[axis] - 1) as f64;
        if coordinate.is_finite() {
            let nearest = coordinate.round();
            // Recognize the exact stored grid coordinate, without a tolerance
            // or snapping a distinct point. This keeps grid-to-itself identity
            // even when a decimal spacing does not invert exactly in f64.
            if nearest >= 0.0
                && nearest <= end
                && self.coordinate(axis, nearest as usize) == position
            {
                coordinate = nearest;
            }
        }
        if !coordinate.is_finite() || coordinate < 0.0 || coordinate > end {
            return Err(RemeshError::SupportOutsideGrid { particle, axis });
        }
        let base = coordinate.floor() as usize;
        let t = coordinate - base as f64;
        if t == 0.0 && self.coordinate(axis, base) != position {
            // A distinct displacement can disappear in the normalization
            // (notably when position/spacing underflows). It must not be
            // published as an exact node or bypass the full-support check.
            return Err(RemeshError::Invalid(
                "particle offset is not representable in grid coordinates",
            ));
        }
        let u = 1.0 - t;
        let weights = [
            -0.5 * t * u * u,
            u * (1.0 + t - 1.5 * t * t),
            t * (1.0 + u - 1.5 * u * u),
            -0.5 * u * t * t,
        ];
        let candidates = [
            base.checked_sub(1),
            Some(base),
            base.checked_add(1),
            base.checked_add(2),
        ];
        let mut indices = [0; 4];
        for k in 0..4 {
            match candidates[k].filter(|&index| index < self.shape[axis]) {
                Some(index) => indices[k] = index,
                None if weights[k] == 0.0 => {}
                None => return Err(RemeshError::SupportOutsideGrid { particle, axis }),
            }
        }
        Ok((indices, weights))
    }
}

/// Local resource bounds admitted before either output buffer is allocated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RemeshBudget {
    /// Maximum borrowed input particles.
    pub max_input_particles: usize,
    /// Maximum output grid nodes, including zero-strength support.
    pub max_grid_particles: usize,
    /// Maximum `17*N_input + 2*N_grid` logical work units.
    pub max_work_units: u64,
    /// Maximum logical particle and area payload bytes, excluding borrowed
    /// input, allocator metadata and capacity rounding.
    pub max_live_bytes: usize,
}

/// A complete remeshed cloud, directly accepted by the existing PSE entries.
#[derive(Debug, Clone, PartialEq)]
pub struct RemeshedCloud2 {
    /// Canonical row-major grid, including zero-circulation nodes. No strength
    /// threshold may remove these nodes without changing diffusion support.
    pub particles: Vec<VortexParticle>,
    /// One `dx*dy` quadrature area per output particle, suitable for
    /// [`crate::diffusion::simulate_viscous_with_cx`].
    pub areas_m2: Vec<f64>,
    /// Authored target grid.
    pub grid: CartesianGrid2,
    /// Exact source-to-stencil visits, including zero-weight visits: `16*N`.
    pub stencil_deposits: u64,
    /// Completed input/stencil/grid visits: `17*N_input + 2*N_grid`.
    pub work_units: u64,
    /// Logical peak owned payload. Compensated-sum scratch becomes `areas_m2`
    /// at publication, so no third grid-sized buffer is allocated.
    pub peak_live_bytes: usize,
}

/// A refusal publishes no output cloud and never modifies the borrowed input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RemeshError {
    /// Invalid geometry or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Declared local resource cap was exceeded.
    Limit {
        resource: &'static str,
        required: u128,
        maximum: u128,
    },
    /// A nonzero interpolation lobe would escape the authored grid.
    SupportOutsideGrid { particle: usize, axis: usize },
    /// Fallible reservation of an admitted output buffer failed.
    Allocation,
    /// Caller cancelled at an input/grid/stencil checkpoint or publication.
    Cancelled,
}
impl core::fmt::Display for RemeshError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar vortex remeshing refused: {self:?}")
    }
}
impl core::error::Error for RemeshError {}

/// Redistribute circulation with the signed M4' tensor kernel.
///
/// Restores an explicitly selected regular quadrature without changing the
/// velocity core or inferring material viscosity. The result preserves
/// `sum Gamma`, `sum Gamma*x`, `sum Gamma*y` and all quadratic position
/// moments to floating-point roundoff. These are point-vortex moments, not
/// momentum inside an arbitrary finite fluid volume. Subsequent PSE retains
/// its own finite-cloud quadrature and domain-padding requirements.
///
/// Uses O(16N+G) work and two output-sized buffers. Compensated circulation
/// deposition reduces signed-cancellation loss; its scratch is reused for
/// output areas. Zero-strength inputs are validated but do not require stencil
/// support, since they transfer nothing. Nonzero tails are never clipped.
/// Cancellation is polled at entry, at most every 256 logical visits, and
/// immediately before publication. A callback can poll the caller's original
/// execution scope. Input is immutable and no partial output escapes.
pub fn remesh(
    particles: &[VortexParticle],
    grid: CartesianGrid2,
    budget: RemeshBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RemeshedCloud2, RemeshError> {
    poll(&mut cancelled)?;
    limit(
        "input particles",
        particles.len() as u128,
        budget.max_input_particles as u128,
    )?;
    limit(
        "grid particles",
        grid.nodes as u128,
        budget.max_grid_particles as u128,
    )?;
    let deposits = (particles.len() as u128) * 16;
    let work = deposits + particles.len() as u128 + 2 * grid.nodes as u128;
    limit("work units", work, budget.max_work_units as u128)?;
    let payload = grid.nodes as u128
        * (core::mem::size_of::<VortexParticle>() + core::mem::size_of::<f64>()) as u128;
    limit("live bytes", payload, budget.max_live_bytes as u128)?;
    let mut output = Vec::new();
    let mut scratch = Vec::new();
    output
        .try_reserve_exact(grid.nodes)
        .map_err(|_| RemeshError::Allocation)?;
    scratch
        .try_reserve_exact(grid.nodes)
        .map_err(|_| RemeshError::Allocation)?;
    let mut completed = 0_u64;
    for y in 0..grid.shape[1] {
        for x in 0..grid.shape[0] {
            tick(&mut completed, &mut cancelled)?;
            output.push(VortexParticle::new(
                [grid.coordinate(0, x), grid.coordinate(1, y)],
                0.0,
            ));
            scratch.push(0.0);
        }
    }
    for (particle, source) in particles.iter().enumerate() {
        tick(&mut completed, &mut cancelled)?;
        if !source.circulation.is_finite() || source.pos.iter().any(|p| !p.is_finite()) {
            return Err(RemeshError::Invalid(
                "finite particle coordinates and circulation required",
            ));
        }
        let (xs, ys) = if source.circulation == 0.0 {
            (([0; 4], [0.0; 4]), ([0; 4], [0.0; 4]))
        } else {
            (
                grid.stencil(source.pos[0], 0, particle)?,
                grid.stencil(source.pos[1], 1, particle)?,
            )
        };
        for j in 0..4 {
            for i in 0..4 {
                tick(&mut completed, &mut cancelled)?;
                // Multiply circulation first: two tiny weights may underflow
                // as a product even when their circulation deposit is finite.
                let deposit = (source.circulation * xs.1[i]) * ys.1[j];
                if deposit != 0.0 {
                    let index = ys.0[j] * grid.shape[0] + xs.0[i];
                    compensated_add(&mut output[index].circulation, &mut scratch[index], deposit)?;
                }
            }
        }
    }
    for (particle, area) in output.iter_mut().zip(&mut scratch) {
        tick(&mut completed, &mut cancelled)?;
        particle.circulation += *area;
        if !particle.circulation.is_finite() {
            return Err(RemeshError::Invalid("unrepresentable remeshed circulation"));
        }
        *area = grid.area;
    }
    debug_assert_eq!(u128::from(completed), work);
    poll(&mut cancelled)?;
    Ok(RemeshedCloud2 {
        particles: output,
        areas_m2: scratch,
        grid,
        stencil_deposits: deposits as u64,
        work_units: completed,
        peak_live_bytes: payload as usize,
    })
}

fn compensated_add(sum: &mut f64, correction: &mut f64, value: f64) -> Result<(), RemeshError> {
    let next = *sum + value;
    let error = if sum.abs() >= value.abs() {
        (*sum - next) + value
    } else {
        (value - next) + *sum
    };
    *correction += error;
    if !next.is_finite() || !correction.is_finite() {
        return Err(RemeshError::Invalid(
            "unrepresentable circulation accumulation",
        ));
    }
    *sum = next;
    Ok(())
}
fn limit(resource: &'static str, required: u128, maximum: u128) -> Result<(), RemeshError> {
    if required > maximum {
        Err(RemeshError::Limit {
            resource,
            required,
            maximum,
        })
    } else {
        Ok(())
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RemeshError> {
    if cancelled() {
        Err(RemeshError::Cancelled)
    } else {
        Ok(())
    }
}
fn tick(completed: &mut u64, cancelled: &mut impl FnMut() -> bool) -> Result<(), RemeshError> {
    if *completed % 256 == 0 {
        poll(cancelled)?;
    }
    *completed += 1;
    Ok(())
}
