//! Wall-dependent convection on the original reciprocal body/gas transaction.

use super::*;

/// Independent wall-temperature solve inside each proposed gas/heat state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallConvectionSolve {
    /// Wall temperature minus its accepted initial temperature [K]. Must
    /// bracket zero; source-domain and thermal-owner refusals still apply.
    pub temperature_increment_bracket_k: [f64; 2],
    /// Maximum difference between queried and solved wall temperature [K].
    pub temperature_tolerance_k: f64,
    /// Maximum safeguarded iterations per gas proposal, in 1..=128.
    pub maximum_iterations: usize,
}

/// Both material evaluations at the converged physical wall endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct WallConvectionSample<S> {
    /// Original source evaluated at `temperature_k`; supplies the actual h.
    pub source: S,
    /// Independent reevaluation at the thermal owner's actual endpoint.
    pub endpoint_source: S,
    /// Wall temperature supplied to the coefficient law [K].
    pub temperature_k: f64,
    /// Queried minus actual thermal endpoint [K].
    pub temperature_residual_k: f64,
    /// dt A (h_used - h_at_endpoint) (T_wall,end - T_gas,end) [J].
    /// Must meet the enclosing exchange's original source tolerance.
    pub boundary_residual_j: f64,
    /// Wall proposals for the final gas/heat trial, at most iterations + 2.
    pub wall_trials: usize,
    /// All wall proposals, including rejected outer gas/heat trials.
    pub total_wall_trials: usize,
}

/// Physical source/thermal failure or bounded wall-feedback failure.
#[derive(Debug)]
pub enum WallConvectionError<E> {
    Source(E),
    Thermal(Error),
    Invalid(&'static str),
    /// Residuals are signed maxima of normalized temperature and work errors.
    OutsideBracket {
        lower_scaled_residual: f64,
        upper_scaled_residual: f64,
    },
    NotConverged {
        iterations: usize,
        scaled_residual: f64,
    },
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for WallConvectionError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "wall convection feedback: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for WallConvectionError<E> {}

impl FiniteGasThermalEnvironment<'_> {
    /// Resolve convection from the SAME gas and wall endpoint, including laws
    /// needing film temperature or the actual wall-to-gas temperature drop.
    ///
    /// Every trial reuses `LumpedThermalEnvironment::advance` from the original
    /// accepted body, so source caloric/radiation/transport and Biot admission
    /// remain authoritative. The existing outer solve retains reciprocal Q,
    /// source-selected gas EOS/composition, energy checks and atomic publication.
    /// The callback may reevaluate the original gas model at a declared film
    /// temperature; changing the bulk gas inventory to that temperature would
    /// instead violate its energy balance.
    ///
    /// Acceptance requires BOTH queried wall temperature and actual boundary
    /// work to converge. A small temperature residual cannot hide a steep h(T)
    /// error. The two source evaluations are retained, without source promotion.
    /// Callback effects are not rolled back; it must be a repeatable evaluator.
    /// Arbitrary/discontinuous laws have no convergence or uniqueness guarantee.
    /// This is still the original isothermal-body/well-mixed-gas physical rung,
    /// not a resolved boundary-layer or spatial-temperature model.
    pub fn advance_with_wall_convection<S, E>(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
        wall: WallConvectionSolve,
        mut resolve: impl FnMut(
            &ResolvedThermochemicalGasState,
            f64,
        ) -> Result<ConvectionCoefficient<S>, E>,
    ) -> Result<
        UniformEnthalpyStep<FiniteGasThermalReport<WallConvectionSample<S>>>,
        ConvectionError<WallConvectionError<E>>,
    > {
        let mut total_wall_trials = 0;
        let result = self.advance_with_convection_feedback(cx, input, |gas, _| {
            feedback::solve_scalar(
                wall.temperature_increment_bracket_k,
                1.0,
                wall.maximum_iterations,
                |increment, wall_trials, _| {
                    if !wall.temperature_tolerance_k.is_finite()
                        || wall.temperature_tolerance_k <= 0.0
                    {
                        return Err(WallConvectionError::Invalid(
                            "positive finite wall temperature tolerance required",
                        ));
                    }
                    total_wall_trials += 1;
                    let temperature_k = input.initial.temperature_k() + increment;
                    if !temperature_k.is_finite() || temperature_k <= 0.0 {
                        return Err(WallConvectionError::Invalid(
                            "wall trial temperature must be positive and finite",
                        ));
                    }
                    let coefficient =
                        resolve(gas, temperature_k).map_err(WallConvectionError::Source)?;
                    let mut environment = self.body.clone();
                    environment.temperature_k = gas.state().temperature;
                    environment.convection_w_per_m2_k = coefficient.coefficient_w_per_m2_k;
                    let body = environment.advance(cx, input).map_err(|error| {
                        if cx.checkpoint().is_err() {
                            WallConvectionError::Cancelled
                        } else {
                            WallConvectionError::Thermal(Error::Interface(error))
                        }
                    })?;
                    let endpoint_temperature_k = body.state.temperature_k();
                    let endpoint = resolve(gas, endpoint_temperature_k)
                        .map_err(WallConvectionError::Source)?;
                    if !endpoint.coefficient_w_per_m2_k.is_finite()
                        || endpoint.coefficient_w_per_m2_k < 0.0
                    {
                        return Err(WallConvectionError::Invalid(
                            "endpoint convection coefficient must be finite and nonnegative",
                        ));
                    }
                    let temperature_residual_k = temperature_k - endpoint_temperature_k;
                    let boundary_residual_j = input.duration_s
                        * input.surface_area_m2
                        * (coefficient.coefficient_w_per_m2_k - endpoint.coefficient_w_per_m2_k)
                        * (endpoint_temperature_k - gas.state().temperature);
                    let scaled = (temperature_residual_k.abs() / wall.temperature_tolerance_k)
                        .max(boundary_residual_j.abs() / self.exchange.source_tolerance_j);
                    if !boundary_residual_j.is_finite()
                        || !scaled.is_finite()
                        || (temperature_residual_k == 0.0 && boundary_residual_j != 0.0)
                    {
                        return Err(WallConvectionError::Invalid(
                            "nonfinite or nonrepeatable wall convection residual",
                        ));
                    }
                    Ok((
                        ConvectionCoefficient {
                            coefficient_w_per_m2_k: coefficient.coefficient_w_per_m2_k,
                            source: WallConvectionSample {
                                source: coefficient.source,
                                endpoint_source: endpoint.source,
                                temperature_k,
                                temperature_residual_k,
                                boundary_residual_j,
                                wall_trials,
                                total_wall_trials,
                            },
                        },
                        scaled.copysign(temperature_residual_k),
                    ))
                },
                &mut || cx.checkpoint().is_err(),
            )
            .map_err(|error| match error {
                feedback::ScalarError::Evaluation(error) => error,
                feedback::ScalarError::Invalid(message) => WallConvectionError::Invalid(message),
                feedback::ScalarError::Cancelled => WallConvectionError::Cancelled,
                feedback::ScalarError::OutsideBracket { lower, upper } => {
                    WallConvectionError::OutsideBracket {
                        lower_scaled_residual: lower,
                        upper_scaled_residual: upper,
                    }
                }
                feedback::ScalarError::NotConverged {
                    iterations,
                    residual,
                } => WallConvectionError::NotConverged {
                    iterations,
                    scaled_residual: residual,
                },
            })
        });
        match result {
            Err(ConvectionError::Convection(WallConvectionError::Cancelled)) => {
                Err(ConvectionError::Thermal(Error::Cancelled))
            }
            other => other,
        }
    }
}
