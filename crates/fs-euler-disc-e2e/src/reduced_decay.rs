//! Small-angle, late-stage Euler-disc reduced decay.
//!
//! This numerical reference uses the declared relations
//! `Omega^2 sin(theta) = 4 g / R` and `E approximately equals 3 m g R theta / 2`.
//! It stops at a caller-declared positive validity cutoff before `theta = 0`.
//! It is neither a full rigid-body/contact solve nor a fit to a video or Mould
//! outcome, thin-gap solution, or resolved CFD. The dry channel calls
//! `fs-tribo` directly; the separately named Bildsten boundary-layer channel
//! is energy-only and supplies no force/wrench.

use core::fmt;

use fs_tribo::{
    ConstantContourForce, InputAuthority, InterfaceMedium, InterfaceSystemRef, ResistanceInput,
    ResistanceLaw,
};

/// Standard gravitational acceleration [m/s^2].
pub const STANDARD_GRAVITY_M_PER_S2: f64 = 9.806_65;
/// Retention bound for deterministic one-shot integration.
pub const MAX_REDUCED_DECAY_STEPS: u32 = 200_000;
/// Declared small-angle applicability ceiling [rad].
///
/// This is an input refusal boundary for the reduced reference, not a
/// certification that every state below it is physically resolved.
pub const MAX_SMALL_ANGLE_THETA_RAD: f64 = 0.2;

/// Refusal from the bounded late-stage reference model.
#[derive(Debug, Clone, PartialEq)]
pub enum ReducedDecayError {
    /// A scalar input was non-finite or outside its declared domain.
    InvalidInput { field: &'static str },
    /// A required source identity was blank.
    MissingIdentity { field: &'static str },
    /// The initial state is not above the positive validity cutoff.
    InitialStateOutsideValidity {
        theta_rad: f64,
        cutoff_theta_rad: f64,
    },
    /// No active dissipation channel was declared.
    NoActiveChannel,
    /// The real dry-law boundary refused the supplied interface/input.
    DryLawRefusal { detail: String },
    /// A checked arithmetic result was not finite.
    NonFiniteDerived { field: &'static str },
    /// Refinement would exceed the retained step bound.
    RefinementStepBudgetOverflow,
    /// A caller supplied a malformed run with no retained samples.
    MissingSample,
}

impl fmt::Display for ReducedDecayError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for ReducedDecayError {}

fn finite_positive(value: f64, field: &'static str) -> Result<(), ReducedDecayError> {
    if !(value.is_finite() && value > 0.0) {
        return Err(ReducedDecayError::InvalidInput { field });
    }
    Ok(())
}

fn finite_nonnegative(value: f64, field: &'static str) -> Result<(), ReducedDecayError> {
    if !(value.is_finite() && value >= 0.0) {
        return Err(ReducedDecayError::InvalidInput { field });
    }
    Ok(())
}

fn nonblank(value: &str, field: &'static str) -> Result<(), ReducedDecayError> {
    if value.trim().is_empty() {
        return Err(ReducedDecayError::MissingIdentity { field });
    }
    Ok(())
}

/// Caller-declared dry contour channel.  Its force is supplied to the generic
/// `fs-tribo::ConstantContourForce` law; it is not calibrated here.
#[derive(Debug, Clone, PartialEq)]
pub struct DryContourChannel {
    /// Ordered dry interface and history identity.
    pub interface: InterfaceSystemRef,
    /// Caller-declared normal force [N] retained by `fs-tribo` input checking.
    pub normal_force_n: f64,
    /// Caller-declared contour resistance magnitude [N].
    pub contour_force_n: f64,
}

impl DryContourChannel {
    fn validate(&self) -> Result<(), ReducedDecayError> {
        finite_nonnegative(self.normal_force_n, "dry_contour.normal_force_n")?;
        finite_nonnegative(self.contour_force_n, "dry_contour.contour_force_n")
    }
}

/// Bildsten-style rotating-disc boundary-layer energy closure.
///
/// With `nu = mu/rho`, its power is
/// `C * mu * R^4 * Omega^(5/2) / sqrt(nu)`.  Under the declared late-stage
/// relation, this is proportional to `theta^(-5/4)`, yielding the reference
/// inclination exponent `4/9` while this closure remains applicable.
#[derive(Debug, Clone, PartialEq)]
pub struct BildstenBoundaryLayerChannel {
    /// Named correlation/source identity; no authority is upgraded here.
    pub source_id: String,
    /// Gas density [kg/m^3].
    pub density_kg_per_m3: f64,
    /// Dynamic viscosity [Pa s].
    pub dynamic_viscosity_pa_s: f64,
    /// Caller-declared non-negative dimensionless prefactor, not a fitted
    /// Euler-outcome coefficient.
    pub dimensionless_prefactor: f64,
}

impl BildstenBoundaryLayerChannel {
    fn validate(&self) -> Result<(), ReducedDecayError> {
        nonblank(&self.source_id, "bildsten.source_id")?;
        finite_positive(self.density_kg_per_m3, "bildsten.density_kg_per_m3")?;
        finite_positive(
            self.dynamic_viscosity_pa_s,
            "bildsten.dynamic_viscosity_pa_s",
        )?;
        finite_nonnegative(
            self.dimensionless_prefactor,
            "bildsten.dimensionless_prefactor",
        )
    }

    fn power_w(&self, radius_m: f64, omega_rad_s: f64) -> Result<f64, ReducedDecayError> {
        let kinematic_viscosity = self.dynamic_viscosity_pa_s / self.density_kg_per_m3;
        let power = self.dimensionless_prefactor
            * self.dynamic_viscosity_pa_s
            * radius_m.powi(4)
            * omega_rad_s.powf(2.5)
            / kinematic_viscosity.sqrt();
        if !(power.is_finite() && power >= 0.0) {
            return Err(ReducedDecayError::NonFiniteDerived {
                field: "bildsten.power_w",
            });
        }
        Ok(power)
    }
}

/// Bounded input for the late-stage reference integration.
#[derive(Debug, Clone, PartialEq)]
pub struct ReducedDecayInput {
    /// Disc mass [kg].
    pub mass_kg: f64,
    /// Disc radius [m].
    pub radius_m: f64,
    /// Declared gravitational acceleration [m/s^2].
    pub gravity_m_per_s2: f64,
    /// Initial inclination [rad], strictly above the cutoff and at most the
    /// declared small-angle ceiling.
    pub initial_theta_rad: f64,
    /// Positive terminal validity cutoff [rad].
    pub validity_cutoff_theta_rad: f64,
    /// Fixed requested integration step [s].
    pub timestep_s: f64,
    /// Hard maximum retained integration steps.
    pub maximum_steps: u32,
    /// Optional real dry contour channel.
    pub dry_contour: Option<DryContourChannel>,
    /// Optional energy-only boundary-layer channel.
    pub bildsten_boundary_layer: Option<BildstenBoundaryLayerChannel>,
}

impl ReducedDecayInput {
    /// A named numerical-reference setup with both channels active.  It is not
    /// a material card, physical prediction, or outcome fit.
    pub fn nominal_reference() -> Result<Self, ReducedDecayError> {
        Ok(Self {
            mass_kg: 0.12,
            radius_m: 0.038,
            gravity_m_per_s2: STANDARD_GRAVITY_M_PER_S2,
            initial_theta_rad: 0.05,
            validity_cutoff_theta_rad: 0.001,
            timestep_s: 1.0e-5,
            maximum_steps: 100_000,
            dry_contour: Some(DryContourChannel {
                interface: InterfaceSystemRef::new(
                    "reduced-decay/disc->support",
                    "reduced-decay/reference-history",
                    "caller-declared/numerical-reference",
                    InputAuthority::CallerDeclared,
                    InterfaceMedium::Dry,
                )
                .map_err(|error| ReducedDecayError::DryLawRefusal {
                    detail: error.to_string(),
                })?,
                normal_force_n: 0.12 * STANDARD_GRAVITY_M_PER_S2,
                contour_force_n: 0.002,
            }),
            bildsten_boundary_layer: Some(BildstenBoundaryLayerChannel {
                source_id: "caller-declared/bildsten-energy-only-v1".to_owned(),
                density_kg_per_m3: 1.2,
                dynamic_viscosity_pa_s: 1.8e-5,
                dimensionless_prefactor: 1.0,
            }),
        })
    }

    fn validate(&self) -> Result<(), ReducedDecayError> {
        finite_positive(self.mass_kg, "mass_kg")?;
        finite_positive(self.radius_m, "radius_m")?;
        finite_positive(self.gravity_m_per_s2, "gravity_m_per_s2")?;
        finite_positive(self.initial_theta_rad, "initial_theta_rad")?;
        finite_positive(self.validity_cutoff_theta_rad, "validity_cutoff_theta_rad")?;
        finite_positive(self.timestep_s, "timestep_s")?;
        if self.maximum_steps == 0 || self.maximum_steps > MAX_REDUCED_DECAY_STEPS {
            return Err(ReducedDecayError::InvalidInput {
                field: "maximum_steps",
            });
        }
        if self.initial_theta_rad <= self.validity_cutoff_theta_rad {
            return Err(ReducedDecayError::InitialStateOutsideValidity {
                theta_rad: self.initial_theta_rad,
                cutoff_theta_rad: self.validity_cutoff_theta_rad,
            });
        }
        if self.initial_theta_rad > MAX_SMALL_ANGLE_THETA_RAD {
            return Err(ReducedDecayError::InvalidInput {
                field: "initial_theta_rad_small_angle",
            });
        }
        if let Some(channel) = &self.dry_contour {
            channel.validate()?;
        }
        if let Some(channel) = &self.bildsten_boundary_layer {
            channel.validate()?;
        }
        if self.dry_contour.is_none() && self.bildsten_boundary_layer.is_none() {
            return Err(ReducedDecayError::NoActiveChannel);
        }
        Ok(())
    }

    fn energy_slope_j_per_rad(&self) -> f64 {
        1.5 * self.mass_kg * self.gravity_m_per_s2 * self.radius_m
    }

    /// Late-stage precession rate from `Omega^2 sin(theta) = 4g/R` [rad/s].
    pub fn omega_rad_s(&self, theta_rad: f64) -> Result<f64, ReducedDecayError> {
        finite_positive(theta_rad, "theta_rad")?;
        let value = (4.0 * self.gravity_m_per_s2 / (self.radius_m * theta_rad.sin())).sqrt();
        if !(value.is_finite() && value > 0.0) {
            return Err(ReducedDecayError::NonFiniteDerived {
                field: "omega_rad_s",
            });
        }
        Ok(value)
    }

    /// Approximate late-stage energy `3 m g R theta / 2` [J].
    pub fn energy_j(&self, theta_rad: f64) -> Result<f64, ReducedDecayError> {
        finite_positive(theta_rad, "theta_rad")?;
        let value = self.energy_slope_j_per_rad() * theta_rad;
        if !(value.is_finite() && value > 0.0) {
            return Err(ReducedDecayError::NonFiniteDerived { field: "energy_j" });
        }
        Ok(value)
    }
}

/// Separately retained powers at one state.  `bildsten_boundary_layer_w` is an
/// energy-only closure, and no `fs-flux` exterior-wrench channel is inserted
/// without that distinct generic receipt/dependency being cleanly registered.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChannelPowers {
    /// Dry contour dissipation [W].
    pub dry_contour_w: f64,
    /// Bildsten boundary-layer energy-only dissipation [W].
    pub bildsten_boundary_layer_w: f64,
}

impl ChannelPowers {
    /// Sum [W].
    #[must_use]
    pub const fn total_w(self) -> f64 {
        self.dry_contour_w + self.bildsten_boundary_layer_w
    }
}

/// Per-channel accumulated work [J].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChannelWork {
    /// Dry contour work [J].
    pub dry_contour_j: f64,
    /// Bildsten boundary-layer energy-only work [J].
    pub bildsten_boundary_layer_j: f64,
}

impl ChannelWork {
    /// Sum [J].
    #[must_use]
    pub const fn total_j(self) -> f64 {
        self.dry_contour_j + self.bildsten_boundary_layer_j
    }
}

/// One deterministic retained point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReducedDecaySample {
    /// Elapsed integration time [s].
    pub time_s: f64,
    /// Inclination [rad].
    pub theta_rad: f64,
    /// Derived precession rate [rad/s].
    pub omega_rad_s: f64,
    /// Reduced energy [J].
    pub energy_j: f64,
    /// Active channel powers [W].
    pub powers: ChannelPowers,
    /// Cumulative channel work [J].
    pub work: ChannelWork,
}

/// Why the integration stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReducedDecayTerminal {
    /// The positive validity cutoff was reached; no theta-zero claim is made.
    ValidityCutoff,
    /// The declared bounded work budget was exhausted before the cutoff.
    StepBudgetExhausted,
}

/// Complete numerical-reference output.
#[derive(Debug, Clone, PartialEq)]
pub struct ReducedDecayRun {
    /// Retained deterministic state samples, including initial/final states.
    pub samples: Vec<ReducedDecaySample>,
    /// Structured terminal condition.
    pub terminal: ReducedDecayTerminal,
    /// `initial_energy - final_energy - channel_work` [J].
    pub energy_closure_residual_j: f64,
}

impl ReducedDecayRun {
    /// Final retained sample, refusing externally constructed malformed runs.
    pub fn final_sample(&self) -> Result<&ReducedDecaySample, ReducedDecayError> {
        self.samples.last().ok_or(ReducedDecayError::MissingSample)
    }
}

/// Coarse/fine deterministic comparison.  This is numerical refinement
/// evidence only, not physical validation.
#[derive(Debug, Clone, PartialEq)]
pub struct RefinementEvidence {
    /// Requested-step run.
    pub coarse: ReducedDecayRun,
    /// Half-step run with doubled bounded step budget.
    pub fine: ReducedDecayRun,
    /// Final elapsed-time difference [s].
    pub terminal_time_difference_s: f64,
    /// Final cumulative-work difference [J].
    pub total_work_difference_j: f64,
}

/// Produces the stable, single-record numerical-reference runner output.
///
/// The record carries computed channel accounting and refinement differences;
/// it does not imply experimental validation or a resolved flow field.
#[must_use]
pub fn structured_runner_output(
    run: &ReducedDecayRun,
    refinement: &RefinementEvidence,
) -> Result<String, ReducedDecayError> {
    let final_sample = run.final_sample()?;
    Ok(format!(
        "schema=reduced-decay-v1 terminal={:?} time_s={:.12e} theta_rad={:.12e} energy_j={:.12e} dry_work_j={:.12e} bildsten_work_j={:.12e} closure_residual_j={:.12e}\nrefinement_terminal_time_difference_s={:.12e} refinement_total_work_difference_j={:.12e}",
        run.terminal,
        final_sample.time_s,
        final_sample.theta_rad,
        final_sample.energy_j,
        final_sample.work.dry_contour_j,
        final_sample.work.bildsten_boundary_layer_j,
        run.energy_closure_residual_j,
        refinement.terminal_time_difference_s,
        refinement.total_work_difference_j,
    ))
}

/// Runs the fixed-step reference with exact per-step energy decrement.
pub fn run_reduced_decay(input: &ReducedDecayInput) -> Result<ReducedDecayRun, ReducedDecayError> {
    input.validate()?;
    let initial_energy_j = input.energy_j(input.initial_theta_rad)?;
    let mut theta_rad = input.initial_theta_rad;
    let mut time_s = 0.0;
    let mut work = ChannelWork {
        dry_contour_j: 0.0,
        bildsten_boundary_layer_j: 0.0,
    };
    let capacity = usize::try_from(input.maximum_steps)
        .map_err(|_| ReducedDecayError::InvalidInput {
            field: "maximum_steps",
        })?
        .saturating_add(1);
    let mut samples = Vec::with_capacity(capacity);
    let initial_powers = powers_at(input, theta_rad)?;
    samples.push(sample(input, time_s, theta_rad, initial_powers, work)?);

    for _ in 0..input.maximum_steps {
        let powers = powers_at(input, theta_rad)?;
        let total_power_w = powers.total_w();
        if !(total_power_w.is_finite() && total_power_w > 0.0) {
            return Err(ReducedDecayError::InvalidInput {
                field: "active_channel_power_w",
            });
        }
        let time_to_cutoff_s = (theta_rad - input.validity_cutoff_theta_rad)
            * input.energy_slope_j_per_rad()
            / total_power_w;
        let dt_s = input.timestep_s.min(time_to_cutoff_s);
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(ReducedDecayError::NonFiniteDerived {
                field: "integration_step_s",
            });
        }
        work.dry_contour_j += powers.dry_contour_w * dt_s;
        work.bildsten_boundary_layer_j += powers.bildsten_boundary_layer_w * dt_s;
        time_s += dt_s;
        theta_rad -= total_power_w * dt_s / input.energy_slope_j_per_rad();
        if time_to_cutoff_s <= input.timestep_s {
            theta_rad = input.validity_cutoff_theta_rad;
        }
        let next_powers = powers_at(input, theta_rad)?;
        samples.push(sample(input, time_s, theta_rad, next_powers, work)?);
        if theta_rad <= input.validity_cutoff_theta_rad {
            let final_energy_j = input.energy_j(theta_rad)?;
            return Ok(ReducedDecayRun {
                samples,
                terminal: ReducedDecayTerminal::ValidityCutoff,
                energy_closure_residual_j: initial_energy_j - final_energy_j - work.total_j(),
            });
        }
    }
    let final_energy_j = input.energy_j(theta_rad)?;
    Ok(ReducedDecayRun {
        samples,
        terminal: ReducedDecayTerminal::StepBudgetExhausted,
        energy_closure_residual_j: initial_energy_j - final_energy_j - work.total_j(),
    })
}

/// Runs the requested and half-step numerical references.
pub fn refinement_evidence(
    input: &ReducedDecayInput,
) -> Result<RefinementEvidence, ReducedDecayError> {
    input.validate()?;
    let fine_steps = input
        .maximum_steps
        .checked_mul(2)
        .filter(|steps| *steps <= MAX_REDUCED_DECAY_STEPS)
        .ok_or(ReducedDecayError::RefinementStepBudgetOverflow)?;
    let coarse = run_reduced_decay(input)?;
    let mut fine_input = input.clone();
    fine_input.timestep_s *= 0.5;
    fine_input.maximum_steps = fine_steps;
    let fine = run_reduced_decay(&fine_input)?;
    Ok(RefinementEvidence {
        terminal_time_difference_s: (coarse.final_sample()?.time_s - fine.final_sample()?.time_s)
            .abs(),
        total_work_difference_j: (coarse.final_sample()?.work.total_j()
            - fine.final_sample()?.work.total_j())
        .abs(),
        coarse,
        fine,
    })
}

fn powers_at(
    input: &ReducedDecayInput,
    theta_rad: f64,
) -> Result<ChannelPowers, ReducedDecayError> {
    let omega_rad_s = input.omega_rad_s(theta_rad)?;
    let dry_contour_w = if let Some(channel) = &input.dry_contour {
        let contour_speed_m_per_s = input.radius_m * theta_rad.cos() * omega_rad_s;
        let response = ConstantContourForce {
            force_n: channel.contour_force_n,
        }
        .evaluate(
            &channel.interface,
            ResistanceInput {
                normal_force_n: channel.normal_force_n,
                angular_speed_rad_s: omega_rad_s,
                contour_speed_mps: contour_speed_m_per_s,
            },
        )
        .map_err(|error| ReducedDecayError::DryLawRefusal {
            detail: error.to_string(),
        })?;
        response.dissipated_power_w()
    } else {
        0.0
    };
    let bildsten_boundary_layer_w = if let Some(channel) = &input.bildsten_boundary_layer {
        channel.power_w(input.radius_m, omega_rad_s)?
    } else {
        0.0
    };
    Ok(ChannelPowers {
        dry_contour_w,
        bildsten_boundary_layer_w,
    })
}

fn sample(
    input: &ReducedDecayInput,
    time_s: f64,
    theta_rad: f64,
    powers: ChannelPowers,
    work: ChannelWork,
) -> Result<ReducedDecaySample, ReducedDecayError> {
    Ok(ReducedDecaySample {
        time_s,
        theta_rad,
        omega_rad_s: input.omega_rad_s(theta_rad)?,
        energy_j: input.energy_j(theta_rad)?,
        powers,
        work,
    })
}
