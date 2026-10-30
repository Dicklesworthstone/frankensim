//! Source-bound Ohmic transport on the public constitutive graph.
//!
//! The native conductor owns resistance and Joule power. This adapter adds only
//! named graph ports and the chain rule along the original, continuously
//! supported material-temperature slice. It neither fits a temperature
//! coefficient nor interpolates reciprocal conductivity knots as resistivity.
//! Thermal storage and physical time remain with the caller's thermal owner.

use super::response::ScalarConductorResponse;
use super::{ConductorError, ElectricalPropertyBasis, ResolvedConductor, ResolvedConductorState};
use crate::graph::{
    Differentiability, EnergyBehavior, GraphError, LawNode, NodeDeclaration, NodeOutput, NodeRole,
    Port, TimeParity,
};
use crate::state_point::{
    MaterialStatePointError, ScalarAdmissibility, is_fixed_external_coordinate,
    pinned_temperature_line,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, MaterialCard};
use fs_qty::{Dims, Temperature};

const LAW: &str = "source-bound-uniform-ohmic";
const CURRENT: Dims = Dims([0, 0, 0, 0, 1, 0]);
const VOLTAGE: Dims = Dims([2, 1, -3, 0, -1, 0]);
const POWER: Dims = Dims([2, 1, -3, 0, 0, 0]);

/// One immutable material, geometry and bounded thermal response.
///
/// Inputs are `[temperature_change, current]`: kelvin relative to the original
/// source temperature and signed amperes. Outputs are `[voltage, heat]`: signed
/// volts and nonnegative watts. Zero input is the admitted, unpowered original
/// specimen, so graph admission never probes an unphysical absolute zero.
///
/// Both source resistivity and conductivity use the same native uniform
/// conductor; every query retains its exact key, claim and complete bundle.
/// Changes of material or fixed external conditions construct a new node.
#[derive(Clone, Debug)]
pub struct MaterialOhmicNode {
    response: ScalarConductorResponse,
    temperature_axis: String,
    reference_temperature_k: f64,
    temperature_range_k: [f64; 2],
    source_line: Vec<(f64, f64)>,
    length_m: f64,
    area_m2: f64,
    declaration: NodeDeclaration,
}

impl MaterialOhmicNode {
    /// Bind a complete original conductor to a continuously supported thermal
    /// range and explicit fixed geometry. All other source coordinates stay at
    /// their original values. The common material source-line owner admits
    /// scalar/curve laws and fixed-context slices of sampled surfaces, keeping
    /// every cell breakpoint and refusing holes or discrete-only spans.
    ///
    /// This is memoryless scalar DC transport. It does not supply contacts,
    /// thermal expansion, skin effect, field-dependent conductivity, phase
    /// evolution or measured authority beyond the original material card.
    /// At an interior source knot the tangent uses the right interval; the
    /// final endpoint uses the left interval. This is declared piecewise smooth.
    ///
    /// # Errors
    /// Refuses foreign/missing sources, unowned response coordinates, an invalid
    /// temperature interval, unsupported source paths, nonphysical coefficients,
    /// or unrepresentable geometry/derivatives. No source selection is changed.
    pub fn try_new(
        card: &MaterialCard,
        original: ResolvedConductorState,
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
        length_m: f64,
        area_m2: f64,
    ) -> Result<Self, ConductorError> {
        let response = ScalarConductorResponse::new(card, original, &[])?;
        if response.original().liquid_mass_fraction_axis().is_some() {
            return Err(invalid(
                "Ohmic graph node does not own phase-fraction evolution",
            ));
        }
        let material = response.original().material();
        let context = material.query();
        let reference_temperature_k = *context
            .axes()
            .get(temperature_axis)
            .ok_or_else(|| invalid("Ohmic graph source temperature is absent"))?;
        let quantity = context.axis_quantities().get(temperature_axis).copied();
        if quantity.is_some_and(|q| q.dims() != Temperature::DIMS)
            || !is_fixed_external_coordinate(
                temperature_axis,
                quantity,
                Some(reference_temperature_k),
            )
        {
            return Err(invalid(
                "Ohmic graph needs an original absolute-temperature coordinate",
            ));
        }
        let [low, high] = temperature_range_k;
        if !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || low >= high
            || !(low..=high).contains(&reference_temperature_k)
        {
            return Err(invalid(
                "Ohmic thermal range must contain the original positive temperature",
            ));
        }
        response
            .original()
            .clone()
            .into_uniform_conductor(length_m, area_m2)?;
        let lower = context
            .clone()
            .with(temperature_axis, low)
            .map_err(|_| invalid("Ohmic lower source coordinate is invalid"))?;
        let upper = context
            .with(temperature_axis, high)
            .map_err(|_| invalid("Ohmic upper source coordinate is invalid"))?;
        let mut source_line = None;
        for property in material.properties() {
            if !property.has_only_fixed_external_conditions(Some(temperature_axis), false) {
                return Err(invalid(
                    "Ohmic graph requires one thermal owner and fixed external source conditions",
                ));
            }
            let selected = property.answer().receipt.selected;
            let claim = card
                .claims()
                .claim(selected)
                .ok_or_else(|| invalid("Ohmic original source claim is absent"))?;
            let line = pinned_temperature_line(
                card,
                &claim.key,
                temperature_axis,
                &lower,
                &upper,
                selected,
            )?;
            for &(_, value) in &line {
                let requirement = property.requirement().admissibility();
                let admitted = value.is_finite()
                    && match requirement {
                        ScalarAdmissibility::Finite => true,
                        ScalarAdmissibility::StrictlyPositive => value > 0.0,
                        ScalarAdmissibility::NonNegative => value >= 0.0,
                        ScalarAdmissibility::OpenInterval { lower, upper } => {
                            lower < value && value < upper
                        }
                    };
                if !admitted {
                    return Err(MaterialStatePointError::OutsideConsumerDomain {
                        property: claim.key.name().to_owned(),
                        value,
                        admissibility: requirement,
                    }
                    .into());
                }
            }
            if &claim.key == response.original().source_key() {
                source_line = Some(line);
            }
        }
        let source_line = source_line.ok_or_else(|| invalid("Ohmic source line is absent"))?;
        if source_line.len() < 2
            || source_line.first().map(|p| p.0) != Some(low)
            || source_line.last().map(|p| p.0) != Some(high)
        {
            return Err(invalid(
                "Ohmic source line does not span the complete thermal range",
            ));
        }
        for pair in source_line.windows(2) {
            let width = pair[1].0 - pair[0].0;
            let slope = (pair[1].1 - pair[0].1) / width;
            if !width.is_finite() || width <= 0.0 || !slope.is_finite() {
                return Err(invalid("Ohmic source derivative is unrepresentable"));
            }
        }
        // The original requirement may be more permissive than a conductor's
        // positive, representable response. Admit native transport at every
        // extremum of the piecewise-affine source, including reciprocal and
        // geometry checks; admitting only the original temperature is weaker.
        for &(temperature, _) in &source_line {
            let query = material
                .query()
                .with(temperature_axis, temperature)
                .map_err(|_| invalid("Ohmic source knot coordinate is invalid"))?;
            response
                .at_query_point(&query)?
                .into_uniform_conductor(length_m, area_m2)?;
        }
        let port = |name: &str, dims, parity| Port {
            name: name.into(),
            dims,
            parity,
        };
        let declaration = NodeDeclaration {
            law: LawId(LAW.into()),
            law_version: 1,
            role: NodeRole::BulkTransport,
            inputs: vec![
                port("temperature_change", Temperature::DIMS, TimeParity::Even),
                port("current", CURRENT, TimeParity::Odd),
            ],
            outputs: vec![
                port("voltage", VOLTAGE, TimeParity::Even),
                port("heat", POWER, TimeParity::Even),
            ],
            state_slots: Vec::new(),
            state_schema_version: 1,
            calibration: ValidityDomain::unconstrained().with(
                "temperature_change",
                low - reference_temperature_k,
                high - reference_temperature_k,
            ),
            differentiability: Differentiability::PiecewiseSmooth,
            energy: EnergyBehavior::NonNegativeDissipation,
            tangent_claimed: true,
        };
        Ok(Self {
            response,
            temperature_axis: temperature_axis.into(),
            reference_temperature_k,
            temperature_range_k,
            source_line,
            length_m,
            area_m2,
            declaration,
        })
    }

    /// Original material, selected electrical source key and complete bundle.
    pub fn source(&self) -> &ScalarConductorResponse {
        &self.response
    }
    /// Original absolute temperature defining zero `temperature_change` [K].
    pub const fn reference_temperature_k(&self) -> f64 {
        self.reference_temperature_k
    }
    /// Admitted continuous interval in absolute kelvin.
    pub const fn temperature_range_k(&self) -> [f64; 2] {
        self.temperature_range_k
    }
    /// Fixed physical conductor length [m].
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }
    /// Fixed physical conductor area [m²].
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }

    /// Resolve a native conductor and all original receipts at an absolute
    /// temperature. This does not mutate the original source or commit history.
    pub fn resolve_at_temperature_k(
        &self,
        temperature_k: f64,
    ) -> Result<ResolvedConductor, ConductorError> {
        if !temperature_k.is_finite()
            || !(self.temperature_range_k[0]..=self.temperature_range_k[1]).contains(&temperature_k)
        {
            return Err(invalid(
                "Ohmic graph temperature outside its admitted source interval",
            ));
        }
        let query = self
            .response
            .original()
            .material()
            .query()
            .with(&self.temperature_axis, temperature_k)
            .map_err(|_| invalid("Ohmic graph temperature coordinate is invalid"))?;
        self.response
            .at_query_point(&query)?
            .into_uniform_conductor(self.length_m, self.area_m2)
    }

    fn response_at(
        &self,
        state: &[f64],
        inputs: &[f64],
    ) -> Result<([f64; 2], [f64; 4]), GraphError> {
        for (which, expected, found) in [("state", 0, state.len()), ("inputs", 2, inputs.len())] {
            if expected != found {
                return Err(GraphError::ArityMismatch {
                    node: LAW.into(),
                    which,
                    expected,
                    found,
                });
            }
        }
        if inputs.iter().any(|v| !v.is_finite()) {
            return Err(graph_error("Ohmic graph inputs must be finite"));
        }
        let temperature = self.reference_temperature_k + inputs[0];
        let conductor = self.resolve_at_temperature_k(temperature).map_err(|_| {
            graph_error("Ohmic material temperature or source bundle outside original support")
        })?;
        let current = inputs[1];
        let resistance = conductor.resistance_ohm();
        let voltage = current * resistance;
        let heat = conductor
            .joule_power_w(current)
            .map_err(|_| graph_error("Ohmic material Joule power is unrepresentable"))?;
        let segment = self
            .source_line
            .partition_point(|&(t, _)| t <= temperature)
            .saturating_sub(1)
            .min(self.source_line.len() - 2);
        let [left, right] = [self.source_line[segment], self.source_line[segment + 1]];
        let source_slope = (right.1 - left.1) / (right.0 - left.0);
        let relative = source_slope / conductor.source_property().value_si();
        let resistance_slope = resistance
            * relative
            * match conductor.source_basis() {
                ElectricalPropertyBasis::Resistivity => 1.0,
                ElectricalPropertyBasis::Conductivity => -1.0,
            };
        let voltage_slope = current * resistance_slope;
        let tangent = [
            voltage_slope,
            resistance,
            current * voltage_slope,
            2.0 * voltage,
        ];
        if !voltage.is_finite()
            || (current != 0.0 && voltage == 0.0)
            || !relative.is_finite()
            || !resistance_slope.is_finite()
            || (source_slope != 0.0 && (relative == 0.0 || resistance_slope == 0.0))
            || tangent.iter().any(|v| !v.is_finite())
        {
            return Err(graph_error(
                "Ohmic graph voltage or consistent tangent is unrepresentable",
            ));
        }
        Ok(([voltage, heat], tangent))
    }
}

impl LawNode for MaterialOhmicNode {
    fn declaration(&self) -> &NodeDeclaration {
        &self.declaration
    }

    fn evaluate(&self, state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        let (values, _) = self.response_at(state, inputs)?;
        Ok(NodeOutput {
            outputs: values.to_vec(),
            next_state: Vec::new(),
            dissipation_rate: Some(values[1]),
        })
    }

    fn tangent(&self, state: &[f64], inputs: &[f64]) -> Option<Vec<f64>> {
        self.response_at(state, inputs)
            .ok()
            .map(|(_, tangent)| tangent.to_vec())
    }
}

fn invalid(quantity: &'static str) -> ConductorError {
    ConductorError::InvalidInput { quantity }
}
fn graph_error(obligation: &'static str) -> GraphError {
    GraphError::CardMismatch {
        law: LAW.into(),
        obligation,
    }
}
