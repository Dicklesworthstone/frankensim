//! Reduced flat-plate flexible-base response for the Euler-disc E2E ladder.
//!
//! The rung projects the production `fs-solid` flat plate operators onto one
//! deterministic, load-shaped mode and advances that scalar state under a
//! moving *nodal* normal load. It is deliberately not a resolved finite-patch
//! contact solve, curved shell solve, or multi-patch coupling path. This
//! synchronous bounded rung has no `Cx`/cancellation surface and therefore
//! does not conform to FrankenSim's cancellable hot-kernel invariant; a future
//! cancellable API is required before it can make that project-level claim.

use core::fmt;

use fs_solid::{OperatorDiagnostics, ShellAssembly, ShellError, ShellPlate, ShellSupport};

/// Largest retained integration length for this synchronous campaign rung.
///
/// The largest single synchronous trajectory is 400 steps. Three-level
/// refinement admits only coarse trajectories of at most 100 steps, followed
/// by 200- and 400-step trajectories over the same horizon. This is
/// deliberately not a general long-running integration API because it has no
/// cancellation point.
pub const MAX_BASE_RESPONSE_STEPS: u32 = 400;

/// Largest accepted dimensionless residual of the scaled supported static solve.
///
/// With translation scale `L`, the reported residual is
/// `||D (K u - f)||_2 / ||D f||_2`, where `D` has `L` on translational DOFs and
/// one on rotational DOFs. Both norms are energies, so the quotient is
/// dimensionless.
pub const MAX_REDUCED_SOLVE_SCALED_RESIDUAL: f64 = 1.0e-8;

/// Explicit geometry scope of the reduced response request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseGeometryScope {
    /// One flat, consistently oriented plate assembled by `fs-solid`.
    FlatSinglePatch,
    /// Refused: curved shells require the IGA shell path.
    CurvedShell,
    /// Refused: multi-patch/mortar coupling is not part of this rung.
    MultiPatch,
}

/// Explicit contact scope of the applied force.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactLoadScope {
    /// The normal force moves linearly between two plate nodes.
    NodalNormalLoad,
    /// Refused: a resolved finite patch requires the contact/coupling path.
    ResolvedFinitePatch,
}

/// Level and three-point support declaration consumed by the plate assembly.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LevelSupportInput {
    /// Three pointwise pinned supports; must match `plate.support` exactly.
    pub support: ShellSupport,
    /// Declared upward level normal in the common Cartesian frame.
    pub level_normal: [f64; 3],
    /// Largest admissible tilt from world +Z, in radians.
    pub maximum_tilt_rad: f64,
}

/// A compressive normal force moving linearly from one node to another.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MovingContactLoad {
    /// Start node in `ShellPlate::nodes`.
    pub start_node: usize,
    /// End node in `ShellPlate::nodes`.
    pub end_node: usize,
    /// Compressive load magnitude [N], applied into the base opposite its
    /// declared upward level normal. Zero is useful for free decay.
    pub normal_force_n: f64,
}

/// Complete bounded reduced-base request.
#[derive(Debug, Clone, PartialEq)]
pub struct BaseResponseInput {
    /// Production flat plate/operator source.
    pub plate: ShellPlate,
    /// Exact support and level declaration.
    pub level_support: LevelSupportInput,
    /// Declared geometry scope.
    pub geometry_scope: BaseGeometryScope,
    /// Declared contact scope.
    pub contact_scope: ContactLoadScope,
    /// Moving nodal normal contact load.
    pub load: MovingContactLoad,
    /// Initial generalized displacement [m].
    pub initial_modal_displacement_m: f64,
    /// Initial generalized velocity [m/s].
    pub initial_modal_velocity_m_per_s: f64,
    /// Fixed timestep [s].
    pub timestep_s: f64,
    /// Number of retained integration steps.
    pub steps: u32,
}

/// Stable refusal from the reduced flexible-base rung.
#[derive(Debug, Clone, PartialEq)]
pub enum BaseResponseError {
    /// A required scalar or index was malformed.
    InvalidInput { field: &'static str },
    /// The request tried to claim a higher-fidelity geometry/contact path.
    UnsupportedScope { scope: &'static str },
    /// The level declaration and support card do not identify the same base.
    SupportMismatch,
    /// Production `fs-solid` assembly refused the plate before response work.
    PlateAssembly { detail: String },
    /// The supported plate cannot produce a positive one-mode reduction.
    ModalReduction { detail: &'static str },
    /// The retained trajectory would exceed the declared bounded work budget.
    StepBudgetExceeded,
    /// The timestep is outside the conservative modal-resolution envelope.
    ///
    /// The scalar implicit-midpoint update itself is linearly stable for this
    /// damped linear system; this refusal bounds fixture-scale resolution, not
    /// unconditional stability.
    TimestepOutsideResolution {
        nondimensional_step: f64,
        limit: f64,
    },
    /// The reduced supported static solve did not meet its scaled residual bound.
    ReducedSolveResidual { scaled_residual: f64, limit: f64 },
    /// A finer member of the retained refinement ladder did not improve a
    /// terminal component or total-energy difference.
    RefinementNotImproved { component: &'static str },
    /// A derived quantity became non-finite.
    NonFiniteDerived { field: &'static str },
}

impl fmt::Display for BaseResponseError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for BaseResponseError {}

/// One retained state, energy, and support-reaction row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaseResponseSample {
    /// Elapsed time [s].
    pub time_s: f64,
    /// Deterministic load progress from zero to one.
    pub load_progress: f64,
    /// Generalized moving-load force [N].
    pub modal_force_n: f64,
    /// Modal displacement [m].
    pub modal_displacement_m: f64,
    /// Modal velocity [m/s].
    pub modal_velocity_m_per_s: f64,
    /// Kinetic energy [J].
    pub modal_kinetic_energy_j: f64,
    /// Elastic strain energy [J].
    pub elastic_energy_j: f64,
    /// Cumulative Rayleigh damping work [J].
    pub damping_work_j: f64,
    /// Cumulative external moving-load work [J].
    pub external_work_j: f64,
    /// Euclidean norm of reactions at the three constrained support nodes [N].
    pub support_reaction_norm_n: f64,
}

/// Fixed operator and level diagnostics retained with the complete run.
#[derive(Debug, Clone, PartialEq)]
pub struct BaseResponseDiagnostics {
    /// `fs-solid` raw algebraic/nullity/symmetry report for the full plate.
    pub operator: OperatorDiagnostics,
    /// Actual angle between declared level normal and world +Z [rad].
    pub level_tilt_rad: f64,
    /// Generalized modal mass [kg].
    pub modal_mass_kg: f64,
    /// Projected stiffness [N/m].
    pub modal_stiffness_n_per_m: f64,
    /// Projected Rayleigh damping [N s/m].
    pub modal_damping_n_s_per_m: f64,
    /// Translation amplitude used to make the modal shape dimensionless [m].
    ///
    /// Translational shape entries are dimensionless and rotational entries are
    /// in 1/m, so multiplying the shape by the generalized displacement gives
    /// physical translations [m] and small rotations [rad].
    pub modal_shape_translation_scale_m: f64,
    /// Dimensionless damping ratio `C / (2 sqrt(K M))` for the retained mode.
    pub modal_damping_ratio: f64,
    /// Largest admitted dimensionless timestep for the retained damping ratio.
    pub nondimensional_timestep_limit: f64,
    /// Dimensionless residual of the scaled supported static load-shape solve.
    pub reduced_solve_scaled_residual: f64,
    /// Declared dimensionless acceptance limit for `reduced_solve_scaled_residual`.
    pub reduced_solve_scaled_residual_limit: f64,
    /// Full-DOF supported static unit-load shape before modal normalization.
    ///
    /// Translation entries are [m] and rotation entries are [rad]. It is
    /// retained solely to allow independent reconstruction of the scaled
    /// supported-solve residual; it is not a dynamic state or a continuum
    /// shape claim.
    pub supported_static_shape: Vec<f64>,
    /// Undamped angular frequency [rad/s].
    pub modal_frequency_rad_s: f64,
}

/// Complete deterministic run. The residual is
/// `E_final - E_initial - external_work + damping_work`.
#[derive(Debug, Clone, PartialEq)]
pub struct BaseResponseRun {
    /// Retained initial and per-step states.
    pub samples: Vec<BaseResponseSample>,
    /// Energy closure residual [J].
    pub energy_closure_residual_j: f64,
    /// Positive energy scale used to normalize the closure residual [J].
    pub energy_closure_scale_j: f64,
    /// Dimensionless `abs(energy_closure_residual_j) / energy_closure_scale_j`.
    pub normalized_energy_closure_residual: f64,
    /// Damping applicability and conditioning/support information.
    pub diagnostics: BaseResponseDiagnostics,
}

impl BaseResponseRun {
    /// Final retained sample, if the run originated from this integrator.
    #[must_use]
    pub fn final_sample(&self) -> Option<BaseResponseSample> {
        self.samples.last().copied()
    }
}

/// Three-level evidence for one deterministic timestep-refinement ladder.
#[derive(Debug, Clone, PartialEq)]
pub struct BaseResponseRefinement {
    /// Requested timestep run.
    pub coarse: BaseResponseRun,
    /// Half-timestep, doubled-step run over the same horizon.
    pub medium: BaseResponseRun,
    /// Quarter-timestep, four-times-step run over the same horizon.
    pub fine: BaseResponseRun,
    /// Absolute difference in terminal modal displacement [m].
    pub terminal_displacement_difference_m: f64,
    /// Absolute difference in terminal elastic energy [J].
    pub terminal_elastic_energy_difference_j: f64,
    /// Dimensionless difference in terminal total modal energy.
    pub terminal_normalized_energy_difference: f64,
    /// The medium-to-fine terminal displacement difference did not increase
    /// relative to the coarse-to-medium difference.
    pub displacement_refinement_improved: bool,
    /// The medium-to-fine terminal total-energy difference did not increase
    /// relative to the coarse-to-medium difference.
    pub energy_refinement_improved: bool,
}

/// Assemble and integrate the one-mode production flexible-base rung.
pub fn run_reduced_base_response(
    input: &BaseResponseInput,
) -> Result<BaseResponseRun, BaseResponseError> {
    validate_input(input)?;
    let assembly = input
        .plate
        .assemble()
        .map_err(|error| BaseResponseError::PlateAssembly {
            detail: shell_error_detail(error),
        })?;
    let level_tilt_rad = validate_level_support(input)?;
    let mode = load_shaped_mode(input, &assembly)?;
    let modal_frequency_rad_s = (mode.stiffness / mode.mass).sqrt();
    let modal_damping_ratio = mode.damping / (2.0 * (mode.stiffness * mode.mass).sqrt());
    finite(modal_frequency_rad_s, "modal_frequency")?;
    finite(modal_damping_ratio, "modal_damping_ratio")?;
    let nondimensional_timestep_limit = 0.2 / (1.0 + modal_damping_ratio);
    let diagnostics = BaseResponseDiagnostics {
        operator: assembly.diagnostics.clone(),
        level_tilt_rad,
        modal_mass_kg: mode.mass,
        modal_stiffness_n_per_m: mode.stiffness,
        modal_damping_n_s_per_m: mode.damping,
        modal_shape_translation_scale_m: mode.translation_scale_m,
        modal_damping_ratio,
        nondimensional_timestep_limit,
        reduced_solve_scaled_residual: mode.scaled_solve_residual,
        reduced_solve_scaled_residual_limit: MAX_REDUCED_SOLVE_SCALED_RESIDUAL,
        supported_static_shape: mode.full_static_shape.clone(),
        modal_frequency_rad_s,
    };
    let nondimensional_step = input.timestep_s * diagnostics.modal_frequency_rad_s;
    if !(nondimensional_step.is_finite()
        && nondimensional_step <= diagnostics.nondimensional_timestep_limit)
    {
        return Err(BaseResponseError::TimestepOutsideResolution {
            nondimensional_step,
            limit: diagnostics.nondimensional_timestep_limit,
        });
    }

    let mut displacement = input.initial_modal_displacement_m;
    let mut velocity = input.initial_modal_velocity_m_per_s;
    let mut acceleration = modal_acceleration(&mode, input, 0.0, displacement, velocity);
    finite(acceleration, "initial_acceleration")?;
    let mut damping_work = 0.0;
    let mut external_work = 0.0;
    let mut samples = Vec::with_capacity(input.steps as usize + 1);
    samples.push(sample(
        input,
        &assembly,
        &mode,
        0,
        displacement,
        velocity,
        acceleration,
        damping_work,
        external_work,
    )?);

    for step in 1..=input.steps {
        // Implicit midpoint uses one midpoint state for acceleration, damping
        // work, external work, and the state update. Its linear solve is
        // scalar because this rung retains exactly one supported mode.
        let timestep = input.timestep_s;
        let midpoint_force = mode.force_at(input, (step as f64 - 0.5) / input.steps as f64);
        let midpoint_velocity = (midpoint_force + 2.0 * mode.mass * velocity / timestep
            - mode.stiffness * displacement)
            / (2.0 * mode.mass / timestep + mode.damping + 0.5 * mode.stiffness * timestep);
        let next_displacement = displacement + timestep * midpoint_velocity;
        let next_velocity = 2.0 * midpoint_velocity - velocity;
        let progress = step as f64 / input.steps as f64;
        let next_acceleration =
            modal_acceleration(&mode, input, progress, next_displacement, next_velocity);
        finite(next_displacement, "modal_displacement")?;
        finite(next_velocity, "modal_velocity")?;
        finite(next_acceleration, "modal_acceleration")?;
        damping_work += mode.damping * midpoint_velocity * midpoint_velocity * timestep;
        external_work += midpoint_force * (next_displacement - displacement);
        finite(damping_work, "damping_work")?;
        finite(external_work, "external_work")?;
        displacement = next_displacement;
        velocity = next_velocity;
        acceleration = next_acceleration;
        samples.push(sample(
            input,
            &assembly,
            &mode,
            step,
            displacement,
            velocity,
            acceleration,
            damping_work,
            external_work,
        )?);
    }
    let initial_energy = samples[0].modal_kinetic_energy_j + samples[0].elastic_energy_j;
    let final_energy = samples
        .last()
        .map(|value| value.modal_kinetic_energy_j + value.elastic_energy_j)
        .ok_or(BaseResponseError::NonFiniteDerived { field: "samples" })?;
    let energy_closure_residual_j = final_energy - initial_energy - external_work + damping_work;
    finite(energy_closure_residual_j, "energy_closure_residual")?;
    let energy_closure_scale_j = initial_energy
        .abs()
        .max(final_energy.abs())
        .max(external_work.abs())
        .max(damping_work.abs())
        .max(f64::MIN_POSITIVE);
    let normalized_energy_closure_residual =
        energy_closure_residual_j.abs() / energy_closure_scale_j;
    finite(
        normalized_energy_closure_residual,
        "normalized_energy_closure_residual",
    )?;
    Ok(BaseResponseRun {
        samples,
        energy_closure_residual_j,
        energy_closure_scale_j,
        normalized_energy_closure_residual,
        diagnostics,
    })
}

/// Run a deterministic three-level timestep-halving ladder over one horizon.
///
/// This retains component and energy improvement checks but makes no
/// convergence-order claim.
pub fn refine_reduced_base_response(
    input: &BaseResponseInput,
) -> Result<BaseResponseRefinement, BaseResponseError> {
    if input.steps > MAX_BASE_RESPONSE_STEPS / 4 {
        return Err(BaseResponseError::StepBudgetExceeded);
    }
    let coarse = run_reduced_base_response(input)?;
    let mut medium_input = input.clone();
    medium_input.timestep_s *= 0.5;
    medium_input.steps *= 2;
    let medium = run_reduced_base_response(&medium_input)?;
    let mut fine_input = medium_input;
    fine_input.timestep_s *= 0.5;
    fine_input.steps *= 2;
    let fine = run_reduced_base_response(&fine_input)?;
    let coarse_terminal = coarse
        .final_sample()
        .ok_or(BaseResponseError::NonFiniteDerived {
            field: "coarse samples",
        })?;
    let medium_terminal = medium
        .final_sample()
        .ok_or(BaseResponseError::NonFiniteDerived {
            field: "medium samples",
        })?;
    let fine_terminal = fine
        .final_sample()
        .ok_or(BaseResponseError::NonFiniteDerived {
            field: "fine samples",
        })?;
    let coarse_terminal_energy =
        coarse_terminal.modal_kinetic_energy_j + coarse_terminal.elastic_energy_j;
    let medium_terminal_energy =
        medium_terminal.modal_kinetic_energy_j + medium_terminal.elastic_energy_j;
    let fine_terminal_energy =
        fine_terminal.modal_kinetic_energy_j + fine_terminal.elastic_energy_j;
    let coarse_medium_displacement =
        (coarse_terminal.modal_displacement_m - medium_terminal.modal_displacement_m).abs();
    let medium_fine_displacement =
        (medium_terminal.modal_displacement_m - fine_terminal.modal_displacement_m).abs();
    let coarse_medium_energy = (coarse_terminal_energy - medium_terminal_energy).abs();
    let medium_fine_energy = (medium_terminal_energy - fine_terminal_energy).abs();
    let displacement_scale = coarse_terminal
        .modal_displacement_m
        .abs()
        .max(medium_terminal.modal_displacement_m.abs())
        .max(fine_terminal.modal_displacement_m.abs())
        .max(f64::MIN_POSITIVE);
    let energy_scale = coarse_terminal_energy
        .abs()
        .max(medium_terminal_energy.abs())
        .max(fine_terminal_energy.abs())
        .max(f64::MIN_POSITIVE);
    let displacement_refinement_improved = refinement_improved(
        coarse_medium_displacement,
        medium_fine_displacement,
        displacement_scale,
    );
    if !displacement_refinement_improved {
        return Err(BaseResponseError::RefinementNotImproved {
            component: "terminal displacement",
        });
    }
    let energy_refinement_improved =
        refinement_improved(coarse_medium_energy, medium_fine_energy, energy_scale);
    if !energy_refinement_improved {
        return Err(BaseResponseError::RefinementNotImproved {
            component: "terminal total energy",
        });
    }
    let terminal_normalized_energy_difference = medium_fine_energy / energy_scale;
    finite(
        terminal_normalized_energy_difference,
        "terminal_normalized_energy_difference",
    )?;
    Ok(BaseResponseRefinement {
        terminal_displacement_difference_m: medium_fine_displacement,
        terminal_elastic_energy_difference_j: (medium_terminal.elastic_energy_j
            - fine_terminal.elastic_energy_j)
            .abs(),
        terminal_normalized_energy_difference,
        displacement_refinement_improved,
        energy_refinement_improved,
        coarse,
        medium,
        fine,
    })
}

fn refinement_improved(coarse_difference: f64, fine_difference: f64, scale: f64) -> bool {
    fine_difference.is_finite()
        && coarse_difference.is_finite()
        && fine_difference <= coarse_difference + 64.0 * f64::EPSILON * scale
}

#[derive(Debug, Clone)]
struct ReducedMode {
    full_shape: Vec<f64>,
    full_static_shape: Vec<f64>,
    translation_scale_m: f64,
    scaled_solve_residual: f64,
    mass: f64,
    stiffness: f64,
    damping: f64,
}

impl ReducedMode {
    fn force_at(&self, input: &BaseResponseInput, progress: f64) -> f64 {
        dot(&self.full_shape, &full_load(input, progress))
    }
}

fn modal_acceleration(
    mode: &ReducedMode,
    input: &BaseResponseInput,
    progress: f64,
    displacement: f64,
    velocity: f64,
) -> f64 {
    (mode.force_at(input, progress) - mode.damping * velocity - mode.stiffness * displacement)
        / mode.mass
}

fn validate_input(input: &BaseResponseInput) -> Result<(), BaseResponseError> {
    match input.geometry_scope {
        BaseGeometryScope::FlatSinglePatch => {}
        BaseGeometryScope::CurvedShell => {
            return Err(BaseResponseError::UnsupportedScope {
                scope: "curved shell",
            });
        }
        BaseGeometryScope::MultiPatch => {
            return Err(BaseResponseError::UnsupportedScope {
                scope: "multi-patch shell",
            });
        }
    }
    if input.contact_scope != ContactLoadScope::NodalNormalLoad {
        return Err(BaseResponseError::UnsupportedScope {
            scope: "resolved finite-patch contact",
        });
    }
    for (value, field) in [
        (input.load.normal_force_n, "normal_force_n"),
        (
            input.initial_modal_displacement_m,
            "initial_modal_displacement_m",
        ),
        (
            input.initial_modal_velocity_m_per_s,
            "initial_modal_velocity_m_per_s",
        ),
        (input.timestep_s, "timestep_s"),
        (input.level_support.maximum_tilt_rad, "maximum_tilt_rad"),
    ] {
        if !value.is_finite()
            || (field == "timestep_s" && value <= 0.0)
            || (field == "normal_force_n" && value < 0.0)
            || (field == "maximum_tilt_rad" && value < 0.0)
        {
            return Err(BaseResponseError::InvalidInput { field });
        }
    }
    if input.steps == 0 || input.steps > MAX_BASE_RESPONSE_STEPS {
        return Err(BaseResponseError::StepBudgetExceeded);
    }
    if input.load.start_node >= input.plate.nodes.len()
        || input.load.end_node >= input.plate.nodes.len()
    {
        return Err(BaseResponseError::InvalidInput { field: "load node" });
    }
    Ok(())
}

fn validate_level_support(input: &BaseResponseInput) -> Result<f64, BaseResponseError> {
    if input.plate.support != Some(input.level_support.support) {
        return Err(BaseResponseError::SupportMismatch);
    }
    let normal = input.level_support.level_normal;
    let length = norm(normal);
    if !length.is_finite() || (length - 1.0).abs() > 1.0e-10 {
        return Err(BaseResponseError::InvalidInput {
            field: "level_normal",
        });
    }
    let support_alignment = dot(&normal, &input.level_support.support.normal);
    if support_alignment < 1.0 - 1.0e-10 {
        return Err(BaseResponseError::SupportMismatch);
    }
    let tilt = normal[2].clamp(-1.0, 1.0).acos();
    if tilt > input.level_support.maximum_tilt_rad {
        return Err(BaseResponseError::InvalidInput {
            field: "level tilt",
        });
    }
    Ok(tilt)
}

fn load_shaped_mode(
    input: &BaseResponseInput,
    assembly: &ShellAssembly,
) -> Result<ReducedMode, BaseResponseError> {
    // Shape the mode with a unit normal load so a zero-load free trajectory
    // remains a valid energy-consistency case.
    let shape_load = unit_load(input, 0.5);
    let reduced_load: Vec<f64> = assembly
        .free_dofs
        .iter()
        .map(|&index| shape_load[index])
        .collect();
    if norm_slice(&reduced_load) <= 1.0e-14 {
        return Err(BaseResponseError::ModalReduction {
            detail: "load is entirely constrained",
        });
    }
    let mixed_unit_scales = mixed_unit_scales(input, &assembly.free_dofs)?;
    let (reduced_static_shape, scaled_solve_residual) = solve_mixed_unit_system(
        assembly.stiffness.values(),
        assembly.stiffness.dimension(),
        &reduced_load,
        &mixed_unit_scales,
    )
    .ok_or(BaseResponseError::ModalReduction {
        detail: "reduced stiffness is singular",
    })?;
    if scaled_solve_residual > MAX_REDUCED_SOLVE_SCALED_RESIDUAL {
        return Err(BaseResponseError::ReducedSolveResidual {
            scaled_residual: scaled_solve_residual,
            limit: MAX_REDUCED_SOLVE_SCALED_RESIDUAL,
        });
    }
    let translation_scale_m = reduced_static_shape
        .iter()
        .zip(&assembly.free_dofs)
        .filter(|(_, dof)| **dof % 6 < 3)
        .map(|(value, _)| value.abs())
        .fold(0.0_f64, f64::max);
    if !(translation_scale_m.is_finite() && translation_scale_m > 0.0) {
        return Err(BaseResponseError::ModalReduction {
            detail: "mode has no translational displacement scale",
        });
    }
    // `K^-1 f` has translations in metres and rotations in radians.  Dividing
    // by a translational metre scale leaves a dimensionless displacement shape
    // (and rotational entries in 1/m), preserving q[m], M[kg], K[N/m], C[Ns/m]
    // and phi^T F[N].  Mass normalization would destroy those units.
    let reduced_shape: Vec<f64> = reduced_static_shape
        .iter()
        .copied()
        .map(|value| value / translation_scale_m)
        .collect();
    let mut full_static_shape = vec![0.0; assembly.full_mass.dimension()];
    let mut full_shape = vec![0.0; assembly.full_mass.dimension()];
    for ((&index, &static_value), &value) in assembly
        .free_dofs
        .iter()
        .zip(&reduced_static_shape)
        .zip(&reduced_shape)
    {
        full_static_shape[index] = static_value;
        full_shape[index] = value;
    }
    let mass = dot(&reduced_shape, &assembly.mass.apply(&reduced_shape));
    let stiffness = dot(&reduced_shape, &assembly.stiffness.apply(&reduced_shape));
    let damping = assembly.damping.as_ref().map_or(0.0, |matrix| {
        dot(&reduced_shape, &matrix.apply(&reduced_shape))
    });
    if !(mass.is_finite()
        && mass > 0.0
        && stiffness.is_finite()
        && stiffness > 0.0
        && damping.is_finite()
        && damping >= 0.0)
    {
        return Err(BaseResponseError::ModalReduction {
            detail: "non-positive projected operator",
        });
    }
    Ok(ReducedMode {
        full_shape,
        full_static_shape,
        translation_scale_m,
        scaled_solve_residual,
        mass,
        stiffness,
        damping,
    })
}

#[allow(clippy::too_many_arguments)]
fn sample(
    input: &BaseResponseInput,
    assembly: &ShellAssembly,
    mode: &ReducedMode,
    step: u32,
    displacement: f64,
    velocity: f64,
    acceleration: f64,
    damping_work: f64,
    external_work: f64,
) -> Result<BaseResponseSample, BaseResponseError> {
    let progress = step as f64 / input.steps as f64;
    let force = mode.force_at(input, progress);
    let full_displacement: Vec<f64> = mode
        .full_shape
        .iter()
        .map(|value| value * displacement)
        .collect();
    let full_velocity: Vec<f64> = mode
        .full_shape
        .iter()
        .map(|value| value * velocity)
        .collect();
    let full_acceleration: Vec<f64> = mode
        .full_shape
        .iter()
        .map(|value| value * acceleration)
        .collect();
    let mut reaction = assembly.full_mass.apply(&full_acceleration);
    let stiffness = assembly.full_stiffness.apply(&full_displacement);
    for (value, addend) in reaction.iter_mut().zip(stiffness) {
        *value += addend;
    }
    if let Some(damping) = &assembly.full_damping {
        for (value, addend) in reaction.iter_mut().zip(damping.apply(&full_velocity)) {
            *value += addend;
        }
    }
    for (value, applied) in reaction.iter_mut().zip(full_load(input, progress)) {
        *value -= applied;
    }
    let support_reaction_norm_n = input
        .level_support
        .support
        .node_indices
        .iter()
        .flat_map(|&node| reaction[(node * 6)..(node * 6 + 3)].iter())
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    let sample = BaseResponseSample {
        time_s: step as f64 * input.timestep_s,
        load_progress: progress,
        modal_force_n: force,
        modal_displacement_m: displacement,
        modal_velocity_m_per_s: velocity,
        modal_kinetic_energy_j: 0.5 * mode.mass * velocity * velocity,
        elastic_energy_j: 0.5 * mode.stiffness * displacement * displacement,
        damping_work_j: damping_work,
        external_work_j: external_work,
        support_reaction_norm_n,
    };
    for (value, field) in [
        (sample.modal_kinetic_energy_j, "kinetic energy"),
        (sample.elastic_energy_j, "elastic energy"),
        (sample.support_reaction_norm_n, "support reaction"),
    ] {
        finite(value, field)?;
    }
    Ok(sample)
}

fn full_load(input: &BaseResponseInput, progress: f64) -> Vec<f64> {
    full_load_with_magnitude(input, progress, input.load.normal_force_n)
}

fn unit_load(input: &BaseResponseInput, progress: f64) -> Vec<f64> {
    full_load_with_magnitude(input, progress, 1.0)
}

fn full_load_with_magnitude(
    input: &BaseResponseInput,
    progress: f64,
    normal_force_n: f64,
) -> Vec<f64> {
    let mut force = vec![0.0; input.plate.nodes.len() * 6];
    let start_weight = 1.0 - progress;
    let end_weight = progress;
    for (node, weight) in [
        (input.load.start_node, start_weight),
        (input.load.end_node, end_weight),
    ] {
        for component in 0..3 {
            force[node * 6 + component] +=
                -normal_force_n * weight * input.level_support.level_normal[component];
        }
    }
    force
}

fn mixed_unit_scales(
    input: &BaseResponseInput,
    free_dofs: &[usize],
) -> Result<Vec<f64>, BaseResponseError> {
    let mut translation_scale_m = 0.0_f64;
    for (index, node) in input.plate.nodes.iter().enumerate() {
        for other in input.plate.nodes.iter().skip(index + 1) {
            let separation = [
                node.position_m[0] - other.position_m[0],
                node.position_m[1] - other.position_m[1],
                node.position_m[2] - other.position_m[2],
            ];
            translation_scale_m = translation_scale_m.max(norm(separation));
        }
    }
    if !(translation_scale_m.is_finite() && translation_scale_m > 0.0) {
        return Err(BaseResponseError::ModalReduction {
            detail: "plate has no finite translational scaling length",
        });
    }
    Ok(free_dofs
        .iter()
        .map(|dof| {
            if *dof % 6 < 3 {
                translation_scale_m
            } else {
                1.0
            }
        })
        .collect())
}

fn solve_mixed_unit_system(
    values: &[f64],
    dimension: usize,
    rhs: &[f64],
    scales: &[f64],
) -> Option<(Vec<f64>, f64)> {
    if values.len() != dimension * dimension || rhs.len() != dimension || scales.len() != dimension
    {
        return None;
    }
    let scaled_rhs: Vec<f64> = rhs
        .iter()
        .zip(scales)
        .map(|(value, scale)| value * scale)
        .collect();
    let scaled_values: Vec<f64> = values
        .iter()
        .enumerate()
        .map(|(index, value)| value * scales[index / dimension] * scales[index % dimension])
        .collect();
    let scaled_shape = solve_dense(&scaled_values, dimension, &scaled_rhs)?;
    let shape: Vec<f64> = scaled_shape
        .iter()
        .zip(scales)
        .map(|(value, scale)| value * scale)
        .collect();
    let residual_norm_sq: f64 = (0..dimension)
        .map(|row| {
            let residual = values[(row * dimension)..((row + 1) * dimension)]
                .iter()
                .zip(&shape)
                .map(|(value, shape_value)| value * shape_value)
                .sum::<f64>()
                - rhs[row];
            let scaled_residual = scales[row] * residual;
            scaled_residual * scaled_residual
        })
        .sum();
    let force_norm_sq: f64 = scaled_rhs.iter().map(|value| value * value).sum();
    if !(force_norm_sq.is_finite() && force_norm_sq > 0.0) {
        return None;
    }
    let scaled_residual = residual_norm_sq.sqrt() / force_norm_sq.sqrt();
    scaled_residual
        .is_finite()
        .then_some((shape, scaled_residual))
}

fn solve_dense(values: &[f64], dimension: usize, rhs: &[f64]) -> Option<Vec<f64>> {
    if dimension == 0 || values.len() != dimension * dimension || rhs.len() != dimension {
        return None;
    }
    let mut matrix = values.to_vec();
    let mut vector = rhs.to_vec();
    let matrix_scale = matrix
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    if !(matrix_scale.is_finite() && matrix_scale > 0.0) {
        return None;
    }
    let pivot_tolerance = f64::EPSILON * (dimension as f64).max(1.0) * matrix_scale;
    for column in 0..dimension {
        let pivot = (column..dimension).max_by(|&left, &right| {
            matrix[left * dimension + column]
                .abs()
                .total_cmp(&matrix[right * dimension + column].abs())
        })?;
        if matrix[pivot * dimension + column].abs() <= pivot_tolerance {
            return None;
        }
        for index in column..dimension {
            matrix.swap(column * dimension + index, pivot * dimension + index);
        }
        vector.swap(column, pivot);
        let pivot_value = matrix[column * dimension + column];
        for row in (column + 1)..dimension {
            let factor = matrix[row * dimension + column] / pivot_value;
            matrix[row * dimension + column] = 0.0;
            for index in (column + 1)..dimension {
                matrix[row * dimension + index] -= factor * matrix[column * dimension + index];
            }
            vector[row] -= factor * vector[column];
        }
    }
    let mut solution = vec![0.0; dimension];
    for row in (0..dimension).rev() {
        let residual: f64 = matrix[(row * dimension + row + 1)..((row + 1) * dimension)]
            .iter()
            .zip(&solution[(row + 1)..])
            .map(|(a, x)| a * x)
            .sum();
        solution[row] = (vector[row] - residual) / matrix[row * dimension + row];
    }
    solution
        .iter()
        .all(|value| value.is_finite())
        .then_some(solution)
}

fn shell_error_detail(error: ShellError) -> String {
    error.to_string()
}
fn finite(value: f64, field: &'static str) -> Result<(), BaseResponseError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(BaseResponseError::NonFiniteDerived { field })
    }
}
fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}
fn norm(value: [f64; 3]) -> f64 {
    dot(&value, &value).sqrt()
}
fn norm_slice(value: &[f64]) -> f64 {
    dot(value, value).sqrt()
}
