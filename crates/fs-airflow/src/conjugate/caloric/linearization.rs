//! Original-source coolant derivatives in the native ordered segment march.
//!
//! One forward recurrence and its reverse apply inlet, wall, mass-flow and
//! conductance responses without a dense Jacobian, material re-query or new solve.

use super::{AirMarch, AirflowError, CaloricAirPath, Cx, finite_coupling};

/// Directional changes of the three physical outputs on one original segment.
/// The same layout accepts objective seeds for the transpose operation.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CaloricSegmentTangent {
    pub outlet_temperature_k: f64,
    pub heat_rate_w: f64,
    pub reference_temperature_k: f64,
}

/// Directional changes of the final outlet and complete path heat transfer.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CaloricPathTangent {
    pub outlet_temperature_k: f64,
    pub total_heat_rate_w: f64,
}

/// Objective gradients with respect to the two shared stream inputs.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CaloricFlowGradient {
    /// Objective units per kelvin of inlet temperature.
    pub inlet_temperature: f64,
    /// Objective units per kg/s of prescribed mass flow.
    pub mass_flow: f64,
}

#[derive(Clone, Copy, Debug)]
struct SegmentDerivative {
    inlet_gain: f64,
    inlet_loss: f64,
    wall_enthalpy: f64,
    mean_inlet_enthalpy: f64,
    mean_wall: f64,
    outlet_cp: f64,
    exposure: f64,
    remaining_drive: f64,
    mean_log_exposure: f64,
}

/// Prepared derivatives of one original source-Cp inlet/wall operating point.
///
/// Pressure, composition and the original material source remain fixed.
/// Original material Cp remains variable along the stream. The borrowed path
/// and wall temperatures retain the immutable source and operating point.
/// Repeated actions use O(segments) work and caller-owned output, allocate
/// nothing on success, and never evaluate the material or invert calorics.
/// `apply` and `apply_transpose` hold flow and conductance fixed. Their
/// `with_flow` counterparts also differentiate prescribed mass flow and each
/// segment's total conductance h*A. Hydraulic, geometric and correlation owners
/// must supply their own control derivatives; none are inferred here.
#[derive(Debug)]
pub struct CaloricPathLinearization<'p> {
    path: &'p CaloricAirPath,
    wall_temperatures: &'p [f64],
    inlet_cp: f64,
    derivatives: Vec<SegmentDerivative>,
}

impl CaloricAirPath {
    /// March once and prepare inlet/wall and flow/conductance derivatives.
    ///
    /// The original exchange retains stable weak-transfer heat and Robin-mean
    /// derivatives, as well as the exponential equilibrium limit. A reservoir
    /// may lie outside material support, but all actual fluid endpoints must
    /// be strictly inside it for this two-sided linearization. Redundant Cp
    /// knots are admitted: continuous Cp has a continuous first derivative of
    /// the caloric inverse even when dCp/dT changes there.
    pub fn linearize<'p>(
        &'p self,
        cx: &Cx<'_>,
        wall_temperatures_k: &'p [f64],
    ) -> Result<(AirMarch, CaloricPathLinearization<'p>), AirflowError> {
        poll(cx, 0)?;
        let span = self.caloric_law().temperature_span_k();
        interior(self.inlet_temperature_k(), span)?;
        let inlet_cp = self
            .caloric_law()
            .specific_heat_cp_j_kg_k(self.inlet_temperature_k())
            .map_err(AirflowError::CaloricSource)?;
        let mut derivatives = Vec::with_capacity(self.segments().len());
        let (march, _) = self.march_with_response(wall_temperatures_k, |exchange, exposure| {
            poll(cx, derivatives.len())?;
            interior(exchange.temperature_k, span)?;
            let outlet_cp = self
                .caloric_law()
                .specific_heat_cp_j_kg_k(exchange.temperature_k)
                .map_err(AirflowError::CaloricSource)?;
            let mean_inlet_enthalpy = finite_coupling(
                "source-caloric mean inlet derivative",
                exchange.inlet_enthalpy_loss / exposure,
            )?;
            derivatives.push(SegmentDerivative {
                inlet_gain: exchange.inlet_enthalpy_gain,
                inlet_loss: exchange.inlet_enthalpy_loss,
                wall_enthalpy: exchange.reservoir_enthalpy_derivative_j_kg_k,
                mean_inlet_enthalpy,
                mean_wall: exchange.mean_reservoir_derivative,
                outlet_cp,
                exposure,
                remaining_drive: exchange.remaining_temperature_difference_k,
                mean_log_exposure: exchange.mean_log_exposure_derivative_k,
            });
            Ok(())
        })?;
        poll(cx, derivatives.len())?;
        Ok((
            march,
            CaloricPathLinearization {
                path: self,
                wall_temperatures: wall_temperatures_k,
                inlet_cp,
                derivatives,
            },
        ))
    }
}

impl<'p> CaloricPathLinearization<'p> {
    #[must_use]
    pub const fn path(&self) -> &'p CaloricAirPath {
        self.path
    }

    #[must_use]
    pub const fn wall_temperatures_k(&self) -> &'p [f64] {
        self.wall_temperatures
    }

    /// Apply inlet and ordered wall-temperature directions to original outputs.
    ///
    /// Output is scratch storage and may be partly written on refusal or
    /// cancellation; no physical state is mutated. The returned total uses a
    /// compensated sum of local heat changes, retaining weak exchange that
    /// subtracting two final enthalpy derivatives would lose.
    pub fn apply(
        &self,
        cx: &Cx<'_>,
        inlet_temperature_k: f64,
        wall_temperatures_k: &[f64],
        output: &mut [CaloricSegmentTangent],
    ) -> Result<CaloricPathTangent, AirflowError> {
        self.apply_directions(cx, inlet_temperature_k, wall_temperatures_k, None, output)
    }

    /// Differentiate the same march with prescribed flow and h*A controls.
    ///
    /// Directions have physical units: kelvin, kg/s, and W/K. Conductance
    /// directions follow the original segment order. The inlet/wall source Cp
    /// response and the direct change in transported heat m_dot*delta_h both
    /// participate. No finite differences or additional caloric solves occur.
    /// The returned Robin mean includes the conductance change, so differentiating
    /// G*(T_wall-T_mean) gives the same heat tangent as the stream enthalpy.
    pub fn apply_with_flow(
        &self,
        cx: &Cx<'_>,
        inlet_temperature_k: f64,
        mass_flow_kg_per_s: f64,
        wall_temperatures_k: &[f64],
        conductances_w_per_k: &[f64],
        output: &mut [CaloricSegmentTangent],
    ) -> Result<CaloricPathTangent, AirflowError> {
        self.apply_directions(
            cx,
            inlet_temperature_k,
            wall_temperatures_k,
            Some((mass_flow_kg_per_s, conductances_w_per_k)),
            output,
        )
    }

    fn apply_directions(
        &self,
        cx: &Cx<'_>,
        inlet_temperature_k: f64,
        wall_temperatures_k: &[f64],
        controls: Option<(f64, &[f64])>,
        output: &mut [CaloricSegmentTangent],
    ) -> Result<CaloricPathTangent, AirflowError> {
        poll(cx, 0)?;
        shape(wall_temperatures_k.len(), self.derivatives.len())?;
        shape(output.len(), self.derivatives.len())?;
        let mass_flow = self.path.mass_flow_kg_per_s();
        let relative_mass = if let Some((mass, conductances)) = controls {
            shape(conductances.len(), self.derivatives.len())?;
            finite_coupling(
                "source-caloric relative mass-flow direction",
                mass / mass_flow,
            )?
        } else {
            0.0
        };
        let mut enthalpy = finite_coupling(
            "source-caloric inlet direction",
            self.inlet_cp * inlet_temperature_k,
        )?;
        let (mut heat, mut correction) = (0.0, 0.0);
        for (index, ((d, &wall), out)) in self
            .derivatives
            .iter()
            .zip(wall_temperatures_k)
            .zip(output.iter_mut())
            .enumerate()
        {
            if index.is_multiple_of(256) {
                poll(cx, index)?;
            }
            finite_coupling("source-caloric wall direction", wall)?;
            let (enthalpy_control, heat_control, mean_control) =
                if let Some((_, directions)) = controls {
                    let segment = &self.path.segments()[index];
                    let conductance = segment.area_m2() * segment.htc_w_per_m2_k();
                    let relative_conductance = finite_coupling(
                        "source-caloric relative conductance direction",
                        directions[index] / conductance,
                    )?;
                    let relative_exposure = relative_conductance - relative_mass;
                    (
                        d.remaining_drive * (d.exposure * relative_exposure),
                        // dQ/dln(m_dot) = G*(T_out-T_mean). This positive
                        // relative moment retains weak flow sensitivity without
                        // subtracting two nearly equal enthalpy transfers.
                        conductance
                            * (d.remaining_drive * relative_conductance
                                + d.mean_log_exposure * relative_mass),
                        d.mean_log_exposure * relative_exposure,
                    )
                } else {
                    (0.0, 0.0, 0.0)
                };
            let next = finite_coupling(
                "source-caloric endpoint tangent",
                d.inlet_gain * enthalpy + d.wall_enthalpy * wall + enthalpy_control,
            )?;
            let local_heat = finite_coupling(
                "source-caloric heat tangent",
                mass_flow * (-d.inlet_loss * enthalpy + d.wall_enthalpy * wall) + heat_control,
            )?;
            *out = CaloricSegmentTangent {
                outlet_temperature_k: finite_coupling(
                    "source-caloric temperature tangent",
                    next / d.outlet_cp,
                )?,
                heat_rate_w: local_heat,
                reference_temperature_k: finite_coupling(
                    "source-caloric Robin tangent",
                    d.mean_inlet_enthalpy * enthalpy + d.mean_wall * wall + mean_control,
                )?,
            };
            let adjusted = local_heat - correction;
            let sum = finite_coupling("source-caloric total heat tangent", heat + adjusted)?;
            correction = (sum - heat) - adjusted;
            heat = sum;
            enthalpy = next;
        }
        poll(cx, self.derivatives.len())?;
        Ok(CaloricPathTangent {
            outlet_temperature_k: output
                .last()
                .expect("admitted nonempty path")
                .outlet_temperature_k,
            total_heat_rate_w: heat,
        })
    }

    /// Apply the transpose to segment-output and terminal-output objective seeds.
    ///
    /// Writes one gradient per original wall and returns the inlet-temperature
    /// gradient. Seeds carry the caller's objective units per corresponding
    /// output unit. Caller-owned output is scratch on error, as for `apply`.
    pub fn apply_transpose(
        &self,
        cx: &Cx<'_>,
        segments: &[CaloricSegmentTangent],
        terminal: CaloricPathTangent,
        wall_gradient: &mut [f64],
    ) -> Result<f64, AirflowError> {
        self.apply_pullback(cx, segments, terminal, wall_gradient, None)
            .map(|gradient| gradient.inlet_temperature)
    }

    /// Pull an objective back to inlet temperature, mass flow, walls and h*A.
    ///
    /// Each conductance gradient has objective units per W/K and follows the
    /// original segment order. Work is linear in segment count, with no dense
    /// Jacobian, source queries or allocations. Output buffers are scratch on
    /// refusal, exactly as for the original inlet/wall transpose.
    pub fn apply_transpose_with_flow(
        &self,
        cx: &Cx<'_>,
        segments: &[CaloricSegmentTangent],
        terminal: CaloricPathTangent,
        wall_gradient: &mut [f64],
        conductance_gradient: &mut [f64],
    ) -> Result<CaloricFlowGradient, AirflowError> {
        self.apply_pullback(
            cx,
            segments,
            terminal,
            wall_gradient,
            Some(conductance_gradient),
        )
    }

    fn apply_pullback(
        &self,
        cx: &Cx<'_>,
        segments: &[CaloricSegmentTangent],
        terminal: CaloricPathTangent,
        wall_gradient: &mut [f64],
        mut conductance_gradient: Option<&mut [f64]>,
    ) -> Result<CaloricFlowGradient, AirflowError> {
        poll(cx, 0)?;
        shape(segments.len(), self.derivatives.len())?;
        shape(wall_gradient.len(), self.derivatives.len())?;
        if let Some(gradient) = &conductance_gradient {
            shape(gradient.len(), self.derivatives.len())?;
        }
        let outlet_cp = self
            .derivatives
            .last()
            .expect("admitted nonempty path")
            .outlet_cp;
        let mut adjoint = finite_coupling(
            "source-caloric outlet seed",
            terminal.outlet_temperature_k / outlet_cp,
        )?;
        finite_coupling("source-caloric total heat seed", terminal.total_heat_rate_w)?;
        let mass_flow = self.path.mass_flow_kg_per_s();
        let (mut mass_gradient, mut mass_correction): (f64, f64) = (0.0, 0.0);
        for index in (0..self.derivatives.len()).rev() {
            if index.is_multiple_of(256) {
                poll(cx, index)?;
            }
            let (d, seed) = (self.derivatives[index], segments[index]);
            let endpoint = finite_coupling(
                "source-caloric endpoint seed",
                adjoint + seed.outlet_temperature_k / d.outlet_cp,
            )?;
            let heat = finite_coupling(
                "source-caloric heat seed",
                mass_flow * (seed.heat_rate_w + terminal.total_heat_rate_w),
            )?;
            finite_coupling("source-caloric Robin seed", seed.reference_temperature_k)?;
            wall_gradient[index] = finite_coupling(
                "source-caloric wall pullback",
                d.wall_enthalpy * (endpoint + heat) + d.mean_wall * seed.reference_temperature_k,
            )?;
            if let Some(gradient) = &mut conductance_gradient {
                let segment = &self.path.segments()[index];
                let conductance = segment.area_m2() * segment.htc_w_per_m2_k();
                let endpoint_control = d.exposure * d.remaining_drive;
                let mean_seed = d.mean_log_exposure * seed.reference_temperature_k;
                gradient[index] = finite_coupling(
                    "source-caloric conductance pullback",
                    (endpoint_control * (endpoint + heat) + mean_seed) / conductance,
                )?;
                let contribution = finite_coupling(
                    "source-caloric mass-flow pullback",
                    (-endpoint_control * endpoint + d.exposure * d.mean_log_exposure * heat
                        - mean_seed)
                        / mass_flow,
                )?;
                let sum = finite_coupling(
                    "source-caloric mass-flow gradient sum",
                    mass_gradient + contribution,
                )?;
                mass_correction += if mass_gradient.abs() >= contribution.abs() {
                    (mass_gradient - sum) + contribution
                } else {
                    (contribution - sum) + mass_gradient
                };
                mass_gradient = sum;
            }
            adjoint = finite_coupling(
                "source-caloric inlet pullback",
                d.inlet_gain * endpoint - d.inlet_loss * heat
                    + d.mean_inlet_enthalpy * seed.reference_temperature_k,
            )?;
        }
        poll(cx, self.derivatives.len())?;
        Ok(CaloricFlowGradient {
            inlet_temperature: finite_coupling(
                "source-caloric inlet temperature gradient",
                self.inlet_cp * adjoint,
            )?,
            mass_flow: finite_coupling(
                "source-caloric mass-flow gradient",
                mass_gradient + mass_correction,
            )?,
        })
    }
}

fn shape(found: usize, expected: usize) -> Result<(), AirflowError> {
    if found != expected {
        return Err(AirflowError::SolidResponseArity { expected, found });
    }
    Ok(())
}

fn interior(temperature: f64, span: [f64; 2]) -> Result<(), AirflowError> {
    if temperature <= span[0] || temperature >= span[1] {
        return Err(AirflowError::InvalidConjugateInput {
            field: "source-caloric derivative needs interior fluid temperatures",
            value_bits: temperature.to_bits(),
        });
    }
    Ok(())
}

fn poll(cx: &Cx<'_>, iteration: usize) -> Result<(), AirflowError> {
    cx.checkpoint().map_err(|_| AirflowError::Cancelled {
        iteration,
        references_k: Vec::new(),
    })
}
