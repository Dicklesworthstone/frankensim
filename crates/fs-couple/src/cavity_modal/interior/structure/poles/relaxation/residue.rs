//! Excitation and design slopes of one simple full-field memory-loaded pole.
//!
//! The shared Maxwell owner resolves BOTH original right and transpose modes.
//! The reciprocal structure/fluid velocity blocks make the pencil generally
//! nonsymmetric: the left mode must never be replaced by the right mode or its
//! conjugate. Subsequent load and parameter batches reuse that same pair.

use super::{
    C64, Cx, MaxwellModeError, PreparedStructuralInteriorRelaxationModes, StructuralInterior,
    StructuralInteriorRelaxationError, StructuralInteriorRelaxationPole, invalid, poll, response,
};
use fs_modal::relaxation_modes::{MaxwellModeResidue, MaxwellPoleDirection};

pub use super::super::parameter::StructuralInteriorPoleDirection;
pub use fs_modal::relaxation_modes::MaxwellBranchDirection;

mod fields;
pub use fields::StructuralInteriorRelaxationForceResidue;

/// Additional physical assembly, observation and batch ceilings. Native
/// paired mode refinement retains the original preparation's PoleOptions.
/// Each method admits its complete batch before allocating or solving; these
/// ceilings include its retained physical baseline but exclude the original
/// source/native pair, caller arrays, previous results and allocator overhead.
#[derive(Clone, Copy, Debug)]
pub struct RelaxationResidueOptions {
    /// Both original derivative cancellation and equilibrated alignment must
    /// reach this floor in (0,1). It is not an eigenvalue enclosure.
    pub minimum_derivative_coupling: f64,
    pub max_loads: usize,
    pub max_directions: usize,
    pub max_physical_bytes: usize,
    /// Conservative arithmetic score, not elapsed time or an instruction count.
    pub max_physical_work: u64,
}

impl Default for RelaxationResidueOptions {
    fn default() -> Self {
        Self {
            minimum_derivative_coupling: 1e-10,
            max_loads: 32,
            max_directions: 32,
            max_physical_bytes: 64 * 1024 * 1024,
            max_physical_work: 1_000_000_000,
        }
    }
}

/// One complex Laplace amplitude in the same original input coordinates as
/// the forced cavity source. Positive acoustic flow enters the cavity.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorRelaxationLoad<'d> {
    /// Generalized force conjugate to every original structural displacement.
    pub structural_force: &'d [C64],
    /// Generalized inward volume flow, one entry per pressure test function.
    /// These are assembled trace loads, not pointwise surface velocities.
    pub inward_flow_m3_s: &'d [C64],
}

/// Complete infinitesimal physical direction per caller-declared parameter.
/// The relaxed spatial matrices and explicit branch laws remain separate.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorRelaxationDirection<'d> {
    pub spatial: StructuralInteriorPoleDirection<'d>,
    /// One entry per original branch, including explicit zeros. Each strain
    /// direction covers all original [u,psi] coordinates. Time changes use
    /// d ln(tau)=d(tau)/tau, rather than a dimensional relaxation-time delta.
    pub branches: &'d [MaxwellBranchDirection<'d>],
}

/// Local pole slope; no mode-shape or forced-residue derivative is implied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StructuralInteriorRelaxationPoleDerivative {
    /// ds/dtheta in the original fixed coordinate and time units.
    pub laplace_pole_derivative_per_s: C64,
    /// -Re(ds/dtheta), retaining the signed change in amplitude decay.
    pub amplitude_decay_rate_derivative_per_s: f64,
    /// d|Im(s)|/dtheta for an oscillatory pole. A real pole returns None:
    /// absolute oscillation frequency has no chosen derivative there.
    pub oscillation_frequency_derivative_rad_s: Option<f64>,
}

/// Source-bound original right/transpose pair and its independently observed
/// physical free mode. Loads and parameters do not repeat an eigensolve,
/// form a dense inverse, or perturb the source model.
pub struct PreparedStructuralInteriorRelaxationResidue<'p, 'm, 'a, 'f, 'c> {
    model: &'m StructuralInterior<'a, 'f, 'c>,
    native: MaxwellModeResidue<'p, 'm>,
    mode: StructuralInteriorRelaxationPole,
    options: RelaxationResidueOptions,
}

impl<'m, 'a, 'f, 'c> PreparedStructuralInteriorRelaxationModes<'m, 'a, 'f, 'c> {
    /// Resolve a single simple pole and its excitation normalization once.
    ///
    /// The transpose seeds target T(s)^T, without conjugation. Both solves
    /// use the original source; the transpose residual is independently
    /// rechecked at the right pole. Material poles s=-1/tau, pure material
    /// dark modes, unresolved/repeated poles and inadequate coupling refuse.
    ///
    /// Physical preparation reserves `320*n+128*b+1024` bytes and work score
    /// `256*n*n+512*n+512*b`, n=ns+na, b=branch count. Native PoleOptions
    /// independently bound the shared paired solves and their retained pair.
    #[allow(clippy::too_many_arguments)]
    pub fn refine_residue<'p>(
        &'p self,
        laplace_seed_per_s: C64,
        structural_displacement_seed: &[C64],
        pressure_potential_seed_pa_s: &[C64],
        transpose_structural_seed: &[C64],
        transpose_potential_seed: &[C64],
        options: RelaxationResidueOptions,
        cx: &Cx<'_>,
    ) -> Result<
        PreparedStructuralInteriorRelaxationResidue<'p, 'm, 'a, 'f, 'c>,
        StructuralInteriorRelaxationError,
    > {
        poll(cx)?;
        let ns = self.model.structural_dofs;
        let na = self.model.pressure_dofs();
        let n = (ns + na) as u128;
        let b = self.branches.len() as u128;
        if options.max_loads == 0 || options.max_directions == 0 {
            return Err(invalid("positive load and parameter capacities required"));
        }
        if [structural_displacement_seed, transpose_structural_seed]
            .iter()
            .any(|v| v.len() != ns)
            || [pressure_potential_seed_pa_s, transpose_potential_seed]
                .iter()
                .any(|v| v.len() != na)
        {
            return Err(invalid(
                "right and transpose seeds must cover every original coupled coordinate",
            ));
        }
        options.admit(320 * n + 128 * b + 1024, 256 * n * n + 512 * n + 512 * b)?;
        let mut right = Vec::with_capacity(ns + na);
        right.extend_from_slice(structural_displacement_seed);
        right.extend_from_slice(pressure_potential_seed_pa_s);
        let mut left = Vec::with_capacity(ns + na);
        left.extend_from_slice(transpose_structural_seed);
        left.extend_from_slice(transpose_potential_seed);
        let native = self.modes.refine_residue(
            laplace_seed_per_s,
            &right,
            &left,
            options.minimum_derivative_coupling,
            || cx.checkpoint().is_err(),
        )?;
        // Observe the genuine paired owner's right mode. There is no second
        // right refinement and no synthetic mode or residual receipt.
        let mode = self.observe(native.right_mode().clone(), cx)?;
        poll(cx)?;
        Ok(PreparedStructuralInteriorRelaxationResidue {
            model: self.model,
            native,
            mode,
            options,
        })
    }
}

impl PreparedStructuralInteriorRelaxationResidue<'_, '_, '_, '_, '_> {
    #[must_use]
    pub fn mode(&self) -> &StructuralInteriorRelaxationPole {
        &self.mode
    }

    /// Complete transpose displacement/potential vector in original [u,psi]
    /// order. This is the bilinear left mode, not a second physical free mode.
    #[must_use]
    pub fn transpose_displacement_potential(&self) -> &[C64] {
        self.native.transpose_displacement()
    }

    /// Original-unit ell^T T'(s0)q used by all loads and parameter slopes.
    #[must_use]
    pub fn derivative_normalization(&self) -> C64 {
        self.native.derivative_normalization()
    }

    #[must_use]
    pub fn derivative_alignment(&self) -> f64 {
        self.native.derivative_alignment()
    }

    #[must_use]
    pub fn source_cancellation_ratio(&self) -> f64 {
        self.native.source_cancellation_ratio()
    }

    /// Original transposed mechanical and branch backward residuals, both
    /// evaluated at the published RIGHT pole.
    #[must_use]
    pub fn transpose_source_residuals(&self) -> (f64, f64) {
        self.native.transpose_source_residuals()
    }

    /// Laurent coefficients for a complete nonempty load batch. Each result
    /// multiplies 1/(s-s0); analytic background, initial history, other poles
    /// and a complete forced response are outside this local operation.
    ///
    /// For k loads, additional physical payload is bounded by
    /// `256*n+128*b+1024+k*(80*n+32*b+512)` bytes and work score by
    /// `k*(64*ns*na+512*n+128*b+512)`. Native force application is linear in
    /// n+b; only the existing original interface observation traverses ns*na.
    /// Refusal or cancellation returns no partial batch or changed baseline.
    pub fn force_residues(
        &self,
        loads: &[StructuralInteriorRelaxationLoad<'_>],
        cx: &Cx<'_>,
    ) -> Result<Vec<StructuralInteriorRelaxationForceResidue>, StructuralInteriorRelaxationError>
    {
        poll(cx)?;
        if loads.is_empty() || loads.len() > self.options.max_loads {
            return Err(invalid(
                "nonempty load batch within the declared cap required",
            ));
        }
        let ns = self.model.structural_dofs;
        let na = self.model.pressure_dofs();
        let n = (ns + na) as u128;
        let b = self.mode.memory.len() as u128;
        let k = loads.len() as u128;
        self.options.admit(
            256 * n + 128 * b + 1024 + k * (80 * n + 32 * b + 512),
            k * (64 * ns as u128 * na as u128 + 512 * n + 128 * b + 512),
        )?;
        let mut force = Vec::with_capacity(ns + na);
        let mut results = Vec::with_capacity(loads.len());
        for load in loads {
            poll(cx)?;
            if load.structural_force.len() != ns || load.inward_flow_m3_s.len() != na {
                return Err(invalid(
                    "load must cover every original force and inward flow",
                ));
            }
            force.clear();
            force.extend_from_slice(load.structural_force);
            force.extend_from_slice(load.inward_flow_m3_s);
            let native = self
                .native
                .force_residue(&force, || cx.checkpoint().is_err())?;
            results.push(fields::observe(
                &self.model.coupling,
                ns,
                na,
                native,
                &mut || cx.checkpoint().is_err(),
            )?);
        }
        poll(cx)?;
        Ok(results)
    }

    /// Differentiate the same pole through complete original M/C/K, both
    /// reciprocal dB blocks, and every branch stiffness, log-time and strain.
    /// The native bilinear source contraction uses factored double-double
    /// arithmetic. No eigenvector, equilibration or Newton iteration is
    /// differentiated. Signed directions do not replace positive source laws.
    ///
    /// One direction is assembled at a time, so payload is
    /// `24*n*n+256*n+128*b+1024+64*k` bytes for k returned slopes; work score
    /// is `k*(1024*n*n*(b+1)+128*n+256*b+512)`, including native contractions.
    pub fn pole_derivatives(
        &self,
        directions: &[StructuralInteriorRelaxationDirection<'_>],
        cx: &Cx<'_>,
    ) -> Result<Vec<StructuralInteriorRelaxationPoleDerivative>, StructuralInteriorRelaxationError>
    {
        poll(cx)?;
        if directions.is_empty() || directions.len() > self.options.max_directions {
            return Err(invalid(
                "nonempty parameter batch within the declared cap required",
            ));
        }
        let ns = self.model.structural_dofs;
        let na = self.model.pressure_dofs();
        let n = (ns + na) as u128;
        let b = self.mode.memory.len() as u128;
        let k = directions.len() as u128;
        self.options.admit(
            24 * n * n + 256 * n + 128 * b + 1024 + 64 * k,
            k * (1024 * n * n * (b + 1) + 128 * n + 256 * b + 512),
        )?;
        let mut results = Vec::with_capacity(directions.len());
        for direction in directions {
            let blocks = direction
                .spatial
                .assemble(ns, na, &mut || cx.checkpoint().is_err())?;
            let slope = self.native.pole_derivative(
                &MaxwellPoleDirection {
                    mass: &blocks[0],
                    damping: &blocks[1],
                    relaxed_stiffness: &blocks[2],
                    branches: direction.branches,
                },
                || cx.checkpoint().is_err(),
            )?;
            let im = self.mode.laplace_pole_per_s.im;
            results.push(StructuralInteriorRelaxationPoleDerivative {
                laplace_pole_derivative_per_s: slope,
                amplitude_decay_rate_derivative_per_s: -slope.re,
                oscillation_frequency_derivative_rad_s: (im != 0.0)
                    .then_some(im.signum() * slope.im),
            });
        }
        poll(cx)?;
        Ok(results)
    }
}

impl RelaxationResidueOptions {
    fn admit(self, bytes: u128, work: u128) -> Result<(), StructuralInteriorRelaxationError> {
        for (quantity, requested, cap) in [
            (
                "coupled relaxation residue bytes",
                bytes,
                self.max_physical_bytes as u128,
            ),
            (
                "coupled relaxation residue work",
                work,
                self.max_physical_work as u128,
            ),
        ] {
            if requested > cap {
                return Err(MaxwellModeError::Budget {
                    quantity,
                    requested,
                    cap,
                }
                .into());
            }
        }
        Ok(())
    }
}
