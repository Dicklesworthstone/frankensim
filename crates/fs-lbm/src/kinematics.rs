//! Resolved vorticity from the accepted velocity of each native flow owner.
//!
//! These are derivatives of the resolved velocity, not an added vorticity
//! evolution equation or a turbulence closure. Native adapters reuse their
//! existing full-gradient stencils, including force corrections and boundary
//! geometry. No wall value is silently replaced by a zero-velocity solid node.

use fs_exec::CancelGate;

/// A full velocity gradient `gradient[velocity_component][spatial_axis]`.
/// Native grids return lattice units. [`Self::to_physical`] converts a lattice
/// gradient to SI using the actual timestep of the isotropic lattice:
/// `d(u dx/dt)/d(x dx) = grad(u)/dt`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlowKinematics3 {
    gradient: [[f64; 3]; 3],
}

impl FlowKinematics3 {
    /// Admit a finite full gradient and finite derived invariants. Symmetric
    /// stress alone is insufficient: it discards precisely the curl information.
    pub fn from_gradient(gradient: [[f64; 3]; 3]) -> Result<Self, KinematicsError> {
        let value = Self { gradient };
        if gradient.iter().flatten().any(|x| !x.is_finite())
            || !value.divergence().is_finite()
            || !value.enstrophy_density().is_finite()
            || !value.q_criterion().is_finite()
            || value.vortex_stretching().iter().any(|x| !x.is_finite())
        {
            return Err(KinematicsError::NonFiniteGradient);
        }
        Ok(value)
    }

    #[must_use]
    pub fn gradient(self) -> [[f64; 3]; 3] {
        self.gradient
    }

    /// `curl(u)`, twice the local angular velocity of a rigidly rotating fluid.
    #[must_use]
    pub fn vorticity(self) -> [f64; 3] {
        let g = self.gradient;
        [g[2][1] - g[1][2], g[0][2] - g[2][0], g[1][0] - g[0][1]]
    }

    #[must_use]
    pub fn divergence(self) -> f64 {
        self.gradient[0][0] + self.gradient[1][1] + self.gradient[2][2]
    }

    /// Local `|curl(u)|^2/2`, without density or cell-volume weighting.
    #[must_use]
    pub fn enstrophy_density(self) -> f64 {
        0.5 * self.vorticity().iter().map(|x| x * x).sum::<f64>()
    }

    /// Rotation excess `Q = (Omega:Omega - S:S)/2`, using the FULL symmetric
    /// strain, including dilation. Positive Q indicates rotation exceeding
    /// strain locally; it is not by itself a coherent-vortex certificate.
    #[must_use]
    pub fn q_criterion(self) -> f64 {
        let g = self.gradient;
        let mut trace_square = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                trace_square += g[i][j] * g[j][i];
            }
        }
        -0.5 * trace_square
    }

    /// `(omega dot grad) u`, the stretching/tilting term in the vorticity
    /// equation. Compressible dilation, baroclinic and viscous terms are NOT
    /// included, so this is not the full material derivative of vorticity.
    #[must_use]
    pub fn vortex_stretching(self) -> [f64; 3] {
        let w = self.vorticity();
        self.gradient
            .map(|row| row[0] * w[0] + row[1] * w[1] + row[2] * w[2])
    }

    /// Convert a lattice result to SI: curl/divergence are s^-1;
    /// enstrophy, Q and stretching are s^-2. Call exactly once, with the
    /// timestep belonging to the sampled state; `dx` cancels from the gradient.
    pub fn to_physical(self, dt_s: f64) -> Result<Self, KinematicsError> {
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(KinematicsError::Input("timestep must be finite and positive"));
        }
        Self::from_gradient(self.gradient.map(|row| row.map(|x| x / dt_s)))
    }
}

/// Refusals never change the flow owner. Scratch is rebuilt on retry and is
/// exposed only through a successful owner's `kinematics` return value.
#[derive(Clone, Debug, PartialEq)]
pub enum KinematicsError {
    Input(&'static str),
    Budget { requested: usize, maximum: usize },
    Allocation,
    InvalidState { cell: [usize; 3] },
    NonFiniteGradient,
    Cancelled,
}

impl core::fmt::Display for KinematicsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-KINEMATICS3: {self:?}")
    }
}
impl std::error::Error for KinematicsError {}

/// Reused velocity cache and native-gradient output, with an explicit byte
/// budget. Every successful field is x-fastest, `(z*ny+y)*nx+x`; solid cells
/// are `None`, rather than fictitious zero-vorticity fluid. No SGS is enabled.
pub struct KinematicsWorkspace3 {
    pub(crate) dimensions: [usize; 3],
    pub(crate) velocity: Vec<[f64; 3]>,
    pub(crate) field: Vec<Option<FlowKinematics3>>,
}

impl KinematicsWorkspace3 {
    /// Numerical buffer bytes (including Option layout), excluding allocator
    /// bookkeeping. Checked before either buffer allocation.
    pub fn required_bytes(dimensions: [usize; 3]) -> Result<usize, KinematicsError> {
        if dimensions.contains(&0) {
            return Err(KinematicsError::Input("zero grid dimension"));
        }
        dimensions
            .into_iter()
            .try_fold(1usize, |n, d| n.checked_mul(d))
            .and_then(|n| {
                n.checked_mul(
                    core::mem::size_of::<[f64; 3]>()
                        + core::mem::size_of::<Option<FlowKinematics3>>(),
                )
            })
            .ok_or(KinematicsError::Input("kinematics workspace size overflow"))
    }

    pub fn new(dimensions: [usize; 3], maximum_bytes: usize) -> Result<Self, KinematicsError> {
        let requested = Self::required_bytes(dimensions)?;
        if requested > maximum_bytes {
            return Err(KinematicsError::Budget {
                requested,
                maximum: maximum_bytes,
            });
        }
        let cells = dimensions.into_iter().product();
        let mut velocity = Vec::new();
        let mut field = Vec::new();
        velocity
            .try_reserve_exact(cells)
            .map_err(|_| KinematicsError::Allocation)?;
        field
            .try_reserve_exact(cells)
            .map_err(|_| KinematicsError::Allocation)?;
        velocity.resize(cells, [0.0; 3]);
        field.resize(cells, None);
        Ok(Self {
            dimensions,
            velocity,
            field,
        })
    }

    pub(crate) fn admit(&self, dimensions: [usize; 3], gate: &CancelGate) -> Result<(), KinematicsError> {
        checkpoint(gate)?;
        if self.dimensions != dimensions {
            return Err(KinematicsError::Input(
                "kinematics workspace dimensions differ",
            ));
        }
        Ok(())
    }
}

pub(crate) fn checkpoint(gate: &CancelGate) -> Result<(), KinematicsError> {
    if gate.is_requested() {
        Err(KinematicsError::Cancelled)
    } else {
        Ok(())
    }
}
