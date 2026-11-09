//! Original flow-dependent pressure sources condensed into native radial flow.
//!
//! Q = F(dp + H(Q)) gives dQ/ddp = G/(1-G H'), where G is the ORIGINAL
//! radial tangent. Nonincreasing H retains the positive network Schur matrix.
//! No mean-flow inertia, fitted viscosity or separate momentum integrator enters.

use super::*;

/// The original source's pressure rise and its local slope in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlowHeadResponse {
    pub pressure: Pressure,
    /// dH/dQ [Pa s/m3]. Must be finite and nonpositive.
    pub pressure_slope_pa_s_m3: f64,
}

/// A pure, bounded original equipment law. Implementations own their source
/// identity and applicability; the numerical caller never fits or extends it.
/// The declared signed flow interval is inclusive. H must be continuous and
/// nonincreasing throughout it. At a continuous corner, return an explicitly
/// chosen one-sided slope; convergence is checked independently of that slope.
pub trait FlowHeadLaw: core::fmt::Debug + Send + Sync {
    fn flow_interval_m3_s(&self) -> [f64; 2];
    fn evaluate(
        &self,
        flow: VolumetricFlowRate,
    ) -> Result<FlowHeadResponse, Box<dyn std::error::Error + Send + Sync>>;
    /// Density-dependent equipment conversion can require the same fixed fluid
    /// density as the native pipe. Pressure-only source data need not infer one.
    fn reference_density_kg_m3(&self) -> Option<f64> {
        None
    }
}

/// Shared ownership of the exact original source, without cloning source data.
#[derive(Clone, Debug)]
pub struct FlowHeadSource(Arc<dyn FlowHeadLaw>);
impl FlowHeadSource {
    pub fn new<L: FlowHeadLaw + 'static>(source: Arc<L>) -> Self {
        Self(source)
    }
    pub fn original(&self) -> &dyn FlowHeadLaw {
        self.0.as_ref()
    }
}
impl PartialEq for FlowHeadSource {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

/// Separate constitutive pressure acceptance and bounded source work.
#[derive(Clone, Copy, Debug)]
pub struct FlowHeadControls {
    pub absolute_pressure_tolerance: Pressure,
    pub relative_pressure_tolerance: f64,
    /// Admitted radial candidates per branch/pressure trial, not network sweeps.
    pub maximum_iterations: usize,
    /// Every original head query across the complete native graph step.
    pub maximum_source_evaluations: u64,
}
impl Default for FlowHeadControls {
    fn default() -> Self {
        Self {
            absolute_pressure_tolerance: Pressure::new(1e-8),
            relative_pressure_tolerance: 1e-11,
            maximum_iterations: 64,
            maximum_source_evaluations: 131_072,
        }
    }
}
impl FlowHeadControls {
    pub(super) fn validate(self) -> Result<(), GraphError> {
        if !self.absolute_pressure_tolerance.value().is_finite()
            || self.absolute_pressure_tolerance.value() <= 0.0
            || !self.relative_pressure_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_pressure_tolerance)
            || !(1..=1024).contains(&self.maximum_iterations)
            || self.maximum_source_evaluations == 0
        {
            return Err(GraphError::InvalidInput(
                "invalid original head solve controls",
            ));
        }
        Ok(())
    }
}

#[derive(Debug)]
pub enum FlowHeadError {
    Invalid(&'static str),
    Source(Box<dyn std::error::Error + Send + Sync>),
    OutsideHistory {
        initial_m3_s: f64,
        final_m3_s: f64,
        support: [f64; 2],
    },
    NotConverged {
        pressure_residual_pa: Option<f64>,
    },
    SourceBudget,
}
impl core::fmt::Display for FlowHeadError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original flow-dependent head: {self:?}")
    }
}
impl std::error::Error for FlowHeadError {}

/// Source and its actual midpoint evaluation at the accepted radial state.
#[derive(Clone, Debug, PartialEq)]
pub struct AcceptedFlowHead {
    pub source: FlowHeadSource,
    pub flow_interval_m3_s: [f64; 2],
    pub midpoint_flow: VolumetricFlowRate,
    pub response: FlowHeadResponse,
    /// Actual radial drive minus nodal/held head minus H(Q_mid), in Pa.
    pub pressure_residual: Pressure,
    pub effective_conductance_m3_s_pa: f64,
    pub work_j: f64,
    /// Constitutive mismatch work is reported, never deposited as viscous heat.
    pub residual_work_j: f64,
}

fn refusal(branch: usize, source: FlowHeadError) -> TransientCarreauError {
    TransientCarreauError::Head { branch, source }
}

fn query(
    source: &FlowHeadSource,
    flow: f64,
    branch: usize,
    controls: FlowHeadControls,
    evaluations: &mut u64,
    cx: &Cx<'_>,
) -> Result<FlowHeadResponse, TransientCarreauError> {
    checkpoint(cx)?;
    if *evaluations == controls.maximum_source_evaluations {
        return Err(refusal(branch, FlowHeadError::SourceBudget));
    }
    *evaluations += 1;
    let response = source
        .original()
        .evaluate(VolumetricFlowRate::new(flow))
        .map_err(|error| refusal(branch, FlowHeadError::Source(error)))?;
    checkpoint(cx)?;
    if !response.pressure.value().is_finite()
        || !response.pressure_slope_pa_s_m3.is_finite()
        || response.pressure_slope_pa_s_m3 > 0.0
    {
        return Err(refusal(
            branch,
            FlowHeadError::Invalid("finite nonincreasing original pressure-head response required"),
        ));
    }
    Ok(response)
}

/// Condense one original equipment relation into the existing pressure graph.
/// Only source values at admitted flows are queried. Support endpoints bound
/// pressure mathematically; their potentially extreme radial flows are NOT
/// required as trials. Rejected native proposals backtrack to the last admitted
/// radial candidate. Every radial attempt and source query consumes its budget.
#[allow(clippy::too_many_arguments)]
pub(super) fn solve_branch(
    source: &FlowHeadSource,
    pipe: &RadialCarreauPipe,
    initial: &RadialCarreauState,
    output: &mut RadialCarreauState,
    workspace: &mut RadialCarreauWorkspace,
    base_pressure: f64,
    duration_s: f64,
    branch: usize,
    controls: TransientCarreauControls,
    branch_solves: &mut u64,
    head_evaluations: &mut u64,
    cx: &Cx<'_>,
) -> Result<(RadialCarreauStep, AcceptedFlowHead), TransientCarreauError> {
    let support @ [q_min, q_max] = source.original().flow_interval_m3_s();
    if !q_min.is_finite()
        || !q_max.is_finite()
        || q_min >= q_max
        || source
            .original()
            .reference_density_kg_m3()
            .is_some_and(|density| density != pipe.density_kg_m3())
    {
        return Err(refusal(
            branch,
            FlowHeadError::Invalid(
                "finite increasing signed flow support and matching original density required",
            ),
        ));
    }
    let at_min = query(source, q_min, branch, controls.head, head_evaluations, cx)?;
    let at_max = query(source, q_max, branch, controls.head, head_evaluations, cx)?;
    if at_min.pressure.value() < at_max.pressure.value() {
        return Err(refusal(
            branch,
            FlowHeadError::Invalid("original head increases across its support"),
        ));
    }
    let mut lower = finite(
        base_pressure + at_max.pressure.value(),
        "head pressure bracket",
    )?;
    let mut upper = finite(
        base_pressure + at_min.pressure.value(),
        "head pressure bracket",
    )?;
    // Minimize forcing for the first radial trial; no head value is clamped.
    let mut pressure = 0.0_f64.clamp(lower, upper);
    let mut previous = None;
    let mut last_residual = None;
    for _ in 0..controls.head.maximum_iterations {
        checkpoint(cx)?;
        let mut backtracks = 0;
        let step = loop {
            checkpoint(cx)?;
            if *branch_solves == controls.maximum_branch_solves {
                return Err(TransientCarreauError::BranchSolveBudget);
            }
            *branch_solves += 1;
            match pipe.step_from(
                initial,
                output,
                workspace,
                pressure,
                duration_s,
                controls.radial,
                || cx.checkpoint().is_err(),
            ) {
                Ok(step) => break step,
                Err(RadialCarreauError::Interrupted) => {
                    return Err(TransientCarreauError::Radial {
                        branch,
                        source: RadialCarreauError::Interrupted,
                    });
                }
                Err(source) => {
                    let Some(admitted) = previous else {
                        return Err(TransientCarreauError::Radial { branch, source });
                    };
                    let reduced = f64::midpoint(admitted, pressure);
                    if backtracks == controls.maximum_backtracks
                        || reduced == admitted
                        || reduced == pressure
                    {
                        return Err(TransientCarreauError::Radial { branch, source });
                    }
                    pressure = reduced;
                    backtracks += 1;
                }
            }
        };
        let q = step.midpoint_volume_flow_m3_s;
        if step.initial_volume_flow_m3_s < q_min || step.initial_volume_flow_m3_s > q_max {
            return Err(refusal(
                branch,
                FlowHeadError::OutsideHistory {
                    initial_m3_s: step.initial_volume_flow_m3_s,
                    final_m3_s: step.final_volume_flow_m3_s,
                    support,
                },
            ));
        }
        let mut proposed = None;
        if q < q_min {
            // Monotonicity of original radial Q(p) locates any admissible
            // solution above this pressure, without evaluating H outside support.
            lower = lower.max(pressure);
        } else if q > q_max {
            upper = upper.min(pressure);
        } else {
            let response = query(source, q, branch, controls.head, head_evaluations, cx)?;
            let mut residual = Sum::default();
            residual.add(pressure)?;
            residual.add(-base_pressure)?;
            residual.add(-response.pressure.value())?;
            let residual = residual.value()?;
            last_residual = Some(residual);
            let g = step.midpoint_pressure_conductance_m3_s_pa;
            let denominator = finite(
                1.0 - g * response.pressure_slope_pa_s_m3,
                "original head Schur denominator",
            )?;
            let effective = finite(g / denominator, "original head Schur conductance")?;
            if denominator < 1.0 || effective <= 0.0 {
                return Err(refusal(
                    branch,
                    FlowHeadError::Invalid("positive original head Schur response required"),
                ));
            }
            let tolerance = finite(
                controls.head.absolute_pressure_tolerance.value()
                    + controls.head.relative_pressure_tolerance
                        * (pressure - base_pressure)
                            .abs()
                            .max(response.pressure.value().abs()),
                "original head pressure tolerance",
            )?;
            if residual.abs() <= tolerance {
                let accepted = AcceptedFlowHead {
                    source: source.clone(),
                    flow_interval_m3_s: support,
                    midpoint_flow: VolumetricFlowRate::new(q),
                    response,
                    pressure_residual: Pressure::new(residual),
                    effective_conductance_m3_s_pa: effective,
                    work_j: finite(
                        (response.pressure.value() * q) * duration_s,
                        "original head work",
                    )?,
                    residual_work_j: finite(
                        (residual * q) * duration_s,
                        "head constitutive residual work",
                    )?,
                };
                checkpoint(cx)?;
                return Ok((step, accepted));
            }
            if residual < 0.0 {
                lower = lower.max(pressure);
            } else {
                upper = upper.min(pressure);
            }
            let newton = pressure - residual / denominator;
            if newton.is_finite() && newton > lower && newton < upper {
                proposed = Some(newton);
            }
        }
        let next = proposed.unwrap_or_else(|| f64::midpoint(lower, upper));
        if next == pressure || !(lower < next && next < upper) {
            return Err(refusal(
                branch,
                FlowHeadError::NotConverged {
                    pressure_residual_pa: last_residual,
                },
            ));
        }
        previous = Some(pressure);
        pressure = next;
    }
    Err(refusal(
        branch,
        FlowHeadError::NotConverged {
            pressure_residual_pa: last_residual,
        },
    ))
}
