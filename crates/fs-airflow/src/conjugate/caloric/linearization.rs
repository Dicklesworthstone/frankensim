//! Original-source coolant derivatives in the native ordered segment march.
//!
//! One forward recurrence and its reverse apply the complete inlet/wall
//! response without a dense segment Jacobian, material re-query or new solve.

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

#[derive(Clone, Copy, Debug)]
struct SegmentDerivative {
    inlet_gain: f64,
    inlet_loss: f64,
    wall_enthalpy: f64,
    mean_inlet_enthalpy: f64,
    mean_wall: f64,
    outlet_cp: f64,
}

/// Prepared derivatives of one original source-Cp inlet/wall operating point.
///
/// Mass flow, segment geometry/conductance, pressure and composition are fixed.
/// Original material Cp remains variable along the stream. The borrowed path
/// and wall temperatures retain the immutable source and operating point.
/// Repeated actions use O(segments) work and caller-owned output, allocate
/// nothing on success, and never evaluate the material or invert calorics.
/// Solid response, hydraulic feedback and mass-flow controls are separate.
#[derive(Debug)]
pub struct CaloricPathLinearization<'p> {
    path: &'p CaloricAirPath,
    wall_temperatures: &'p [f64],
    inlet_cp: f64,
    derivatives: Vec<SegmentDerivative>,
}

impl CaloricAirPath {
    /// March once and prepare the complete inlet/wall derivative of that march.
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
        poll(cx, 0)?;
        shape(wall_temperatures_k.len(), self.derivatives.len())?;
        shape(output.len(), self.derivatives.len())?;
        let mut enthalpy = finite_coupling(
            "source-caloric inlet direction",
            self.inlet_cp * inlet_temperature_k,
        )?;
        let mass_flow = self.path.mass_flow_kg_per_s();
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
            let next = finite_coupling(
                "source-caloric endpoint tangent",
                d.inlet_gain * enthalpy + d.wall_enthalpy * wall,
            )?;
            let local_heat = finite_coupling(
                "source-caloric heat tangent",
                mass_flow * (-d.inlet_loss * enthalpy + d.wall_enthalpy * wall),
            )?;
            *out = CaloricSegmentTangent {
                outlet_temperature_k: finite_coupling(
                    "source-caloric temperature tangent",
                    next / d.outlet_cp,
                )?,
                heat_rate_w: local_heat,
                reference_temperature_k: finite_coupling(
                    "source-caloric Robin tangent",
                    d.mean_inlet_enthalpy * enthalpy + d.mean_wall * wall,
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
        poll(cx, 0)?;
        shape(segments.len(), self.derivatives.len())?;
        shape(wall_gradient.len(), self.derivatives.len())?;
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
            adjoint = finite_coupling(
                "source-caloric inlet pullback",
                d.inlet_gain * endpoint - d.inlet_loss * heat
                    + d.mean_inlet_enthalpy * seed.reference_temperature_k,
            )?;
        }
        poll(cx, self.derivatives.len())?;
        finite_coupling(
            "source-caloric inlet temperature gradient",
            self.inlet_cp * adjoint,
        )
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
