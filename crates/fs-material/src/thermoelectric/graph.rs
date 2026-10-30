//! Original Kelvin transport on the public constitutive graph.
//!
//! The native three-temperature Seebeck span owns every value and source query.
//! This adapter supplies named ports and its exact piecewise tangent. In
//! particular, the bulk Thomson exchange is retained alongside both terminal
//! Peltier terms; signed heat transfers are not reported as Joule dissipation.

use super::{ABSOLUTE_SEEBECK_PROPERTY, ResolvedSeebeckSpan, ThermoelectricError};
use crate::graph::{
    Differentiability, EnergyBehavior, GraphError, LawNode, NodeDeclaration, NodeOutput, NodeRole,
    Port, TimeParity,
};
use crate::state_point::{MaterialPropertySelection, pinned_temperature_line};
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, MaterialCard, QueryPoint};
use fs_qty::{Dims, Temperature};

const LAW: &str = "source-bound-scalar-kelvin-transfer";
const CURRENT: Dims = Dims([0, 0, 0, 0, 1, 0]);
const VOLTAGE: Dims = Dims([2, 1, -3, 0, -1, 0]);
const POWER: Dims = Dims([2, 1, -3, 0, 0, 0]);

/// Immutable original thermopower source with three actual thermal inputs.
///
/// Inputs, in order, are `temperature_change_a`, `temperature_change_bulk`,
/// `temperature_change_b` [K relative to the original context temperature], and
/// `current` [A, positive a to b]. Outputs are `voltage` [V drop a to b],
/// `heat_a`, `heat_bulk`, `heat_b` [W into their original caloric stores], and
/// `electrical_supply` [W supplied to the electrical owner]. Analytically,
/// sum(heat) + electrical_supply = 0 and sum(heat / actual temperature) = 0.
///
/// Fourier conduction, electrical resistance/Joule heat and caloric storage
/// remain with their existing owners. This node has no state or local storage.
#[derive(Clone, Debug)]
pub struct MaterialSeebeckNode {
    card: MaterialCard,
    context: QueryPoint,
    original: ResolvedSeebeckSpan,
    selection: MaterialPropertySelection,
    temperature_axis: String,
    reference_temperature_k: f64,
    temperature_range_k: [f64; 2],
    source_line: Vec<(f64, f64)>,
    declaration: NodeDeclaration,
}

impl MaterialSeebeckNode {
    /// Pin the original absolute thermopower claim and continuously supported
    /// temperature interval. All nonthermal source coordinates remain fixed.
    /// Scalar, piecewise-linear and supported fixed-context surface slices use
    /// the same native source owner. Relative thermopower is not substituted.
    ///
    /// At an isothermal edge the mean's two temperature derivatives are half
    /// the local source slope. At a source knot this chooses the right interval,
    /// except that the final admitted endpoint uses the left interval. Unequal
    /// temperatures integrate every crossed source interval exactly.
    ///
    /// # Errors
    /// Refuses unsupported spans, source/domain/quantity mismatch, unowned
    /// response coordinates, nonpositive temperatures or unrepresentable values
    /// and derivatives. No geometry, phase kinetics, magnetic coupling or
    /// thermal/electrical storage is inferred.
    pub fn try_new(
        card: &MaterialCard,
        context: &QueryPoint,
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
        selection: MaterialPropertySelection,
    ) -> Result<Self, ThermoelectricError> {
        let reference_temperature_k = *context
            .axes()
            .get(temperature_axis)
            .ok_or_else(|| invalid("Kelvin graph original temperature is absent"))?;
        let [low, high] = temperature_range_k;
        if !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || low >= high
            || !(high - low).is_finite()
            || !(low..=high).contains(&reference_temperature_k)
        {
            return Err(invalid(
                "Kelvin graph needs a positive interval containing the original temperature",
            ));
        }
        let original = ResolvedSeebeckSpan::resolve(
            card,
            context,
            temperature_axis,
            [reference_temperature_k; 2],
            selection,
        )?;
        let property = original
            .bulk_material()
            .property(ABSOLUTE_SEEBECK_PROPERTY)
            .ok_or_else(|| invalid("Kelvin graph original thermopower is absent"))?;
        if !property.has_only_fixed_external_conditions(Some(temperature_axis), false) {
            return Err(invalid(
                "Kelvin graph requires one temperature owner and fixed external source conditions",
            ));
        }
        let selected = property.answer().receipt.selected;
        let selection = MaterialPropertySelection::PinnedByProperty(vec![(
            ABSOLUTE_SEEBECK_PROPERTY.to_owned(),
            selected,
        )]);
        let claim = card
            .claims()
            .claim(selected)
            .ok_or_else(|| invalid("Kelvin graph original claim is absent"))?;
        let lower = context
            .clone()
            .with(temperature_axis, low)
            .map_err(|_| invalid("Kelvin graph lower temperature is invalid"))?;
        let upper = context
            .clone()
            .with(temperature_axis, high)
            .map_err(|_| invalid("Kelvin graph upper temperature is invalid"))?;
        let source_line =
            pinned_temperature_line(card, &claim.key, temperature_axis, &lower, &upper, selected)?;
        if source_line.len() < 2
            || source_line.first().map(|p| p.0) != Some(low)
            || source_line.last().map(|p| p.0) != Some(high)
        {
            return Err(invalid(
                "Kelvin graph source does not span its thermal interval",
            ));
        }
        for pair in source_line.windows(2) {
            let width = pair[1].0 - pair[0].0;
            let change = pair[1].1 - pair[0].1;
            let slope = change / width;
            if !width.is_finite()
                || width <= 0.0
                || pair.iter().any(|p| !p.1.is_finite())
                || !slope.is_finite()
                || change != 0.0 && slope == 0.0
            {
                return Err(invalid("Kelvin graph source derivative is unrepresentable"));
            }
        }
        ResolvedSeebeckSpan::resolve(
            card,
            context,
            temperature_axis,
            temperature_range_k,
            selection.clone(),
        )?;
        let port = |name: &str, dims, parity| Port {
            name: name.into(),
            dims,
            parity,
        };
        let mut calibration = ValidityDomain::unconstrained();
        for name in [
            "temperature_change_a",
            "temperature_change_bulk",
            "temperature_change_b",
        ] {
            calibration = calibration.with(
                name,
                low - reference_temperature_k,
                high - reference_temperature_k,
            );
        }
        let declaration = NodeDeclaration {
            law: LawId(LAW.into()),
            law_version: 1,
            role: NodeRole::ReversibleCoupling,
            inputs: vec![
                port("temperature_change_a", Temperature::DIMS, TimeParity::Even),
                port(
                    "temperature_change_bulk",
                    Temperature::DIMS,
                    TimeParity::Even,
                ),
                port("temperature_change_b", Temperature::DIMS, TimeParity::Even),
                port("current", CURRENT, TimeParity::Odd),
            ],
            outputs: vec![
                port("voltage", VOLTAGE, TimeParity::Even),
                port("heat_a", POWER, TimeParity::Odd),
                port("heat_bulk", POWER, TimeParity::Odd),
                port("heat_b", POWER, TimeParity::Odd),
                port("electrical_supply", POWER, TimeParity::Odd),
            ],
            state_slots: Vec::new(),
            state_schema_version: 1,
            calibration,
            differentiability: Differentiability::PiecewiseSmooth,
            energy: EnergyBehavior::ReversibleTransfer,
            tangent_claimed: true,
        };
        let node = Self {
            card: card.clone(),
            context: context.clone(),
            original,
            selection,
            temperature_axis: temperature_axis.into(),
            reference_temperature_k,
            temperature_range_k,
            source_line,
            declaration,
        };
        node.response_at(&[], &[0.0; 4])
            .map_err(|_| invalid("Kelvin graph original response is unrepresentable"))?;
        Ok(node)
    }

    /// Original selected material and complete source receipts at the reference.
    pub const fn source(&self) -> &ResolvedSeebeckSpan {
        &self.original
    }
    /// Common positive reference for all three temperature-change inputs [K].
    pub const fn reference_temperature_k(&self) -> f64 {
        self.reference_temperature_k
    }
    /// Original continuously supported absolute-temperature interval [K].
    pub const fn temperature_range_k(&self) -> [f64; 2] {
        self.temperature_range_k
    }

    /// Resolve the original native span at `[terminal a, bulk, terminal b]` [K].
    /// Every actual source receipt remains available without advancing history.
    pub fn resolve_at_temperatures_k(
        &self,
        temperature_k: [f64; 3],
    ) -> Result<ResolvedSeebeckSpan, ThermoelectricError> {
        if temperature_k.iter().any(|t| {
            !t.is_finite()
                || !(self.temperature_range_k[0]..=self.temperature_range_k[1]).contains(t)
        }) {
            return Err(invalid(
                "Kelvin graph temperature outside its admitted source interval",
            ));
        }
        let context = self
            .context
            .clone()
            .with(&self.temperature_axis, temperature_k[1])
            .map_err(|_| invalid("Kelvin graph bulk temperature coordinate is invalid"))?;
        ResolvedSeebeckSpan::resolve(
            &self.card,
            &context,
            &self.temperature_axis,
            [temperature_k[0], temperature_k[2]],
            self.selection.clone(),
        )
    }

    // Differentiate the integral mean as integral_0^1 S(x+u(y-x)) du.
    // Integrating the weights (1-u) and u avoids subtracting nearly equal
    // thermopowers and dividing by a tiny temperature difference.
    fn edge_derivatives(&self, x: f64, y: f64) -> [f64; 2] {
        if x == y {
            let i = self
                .source_line
                .partition_point(|&(t, _)| t <= x)
                .saturating_sub(1)
                .min(self.source_line.len() - 2);
            let [a, b] = [self.source_line[i], self.source_line[i + 1]];
            return [0.5 * ((b.1 - a.1) / (b.0 - a.0)); 2];
        }
        let mut result = [0.0; 2];
        let (low, high) = (x.min(y), x.max(y));
        let first = self
            .source_line
            .partition_point(|&(t, _)| t <= low)
            .saturating_sub(1)
            .min(self.source_line.len() - 2);
        for pair in self.source_line[first..].windows(2) {
            if pair[0].0 >= high {
                break;
            }
            let left = pair[0].0.max(low);
            let right = pair[1].0.min(high);
            let u0 = (left - x) / (y - x);
            let u1 = (right - x) / (y - x);
            let width = (u1 - u0).abs();
            let middle = 0.5 * u0 + 0.5 * u1;
            let slope = (pair[1].1 - pair[0].1) / (pair[1].0 - pair[0].0);
            result[0] += slope * width * (1.0 - middle);
            result[1] += slope * width * middle;
        }
        result
    }

    fn response_at(
        &self,
        state: &[f64],
        inputs: &[f64],
    ) -> Result<([f64; 5], [f64; 20]), GraphError> {
        for (which, expected, found) in [("state", 0, state.len()), ("inputs", 4, inputs.len())] {
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
            return Err(graph_error("Kelvin graph inputs must be finite"));
        }
        let temperatures = [0, 1, 2].map(|i| self.reference_temperature_k + inputs[i]);
        let span = self.resolve_at_temperatures_k(temperatures).map_err(|_| {
            graph_error("Kelvin material temperatures or sources outside original support")
        })?;
        let current = inputs[3];
        let transfer = span
            .transfer(current)
            .map_err(|_| graph_error("Kelvin transfer is unrepresentable"))?;
        let [a, bulk, b] = temperatures;
        let [left, right] = span.edge_seebeck_v_k();
        let [la, lm] = self.edge_derivatives(a, bulk);
        let [rm, rb] = self.edge_derivatives(bulk, b);
        let sa = span.material()[0]
            .property(ABSOLUTE_SEEBECK_PROPERTY)
            .ok_or_else(|| graph_error("Kelvin terminal a source disappeared"))?
            .value_si();
        let sb = span.material()[1]
            .property(ABSOLUTE_SEEBECK_PROPERTY)
            .ok_or_else(|| graph_error("Kelvin terminal b source disappeared"))?
            .value_si();
        let tangent = [
            -sa,
            0.0,
            sb,
            0.0,
            -current * (left + a * la),
            -current * a * lm,
            0.0,
            -a * left,
            current * bulk * la,
            current * (left - right + bulk * (lm - rm)),
            -current * bulk * rb,
            bulk * (left - right),
            0.0,
            current * b * rm,
            current * (right + b * rb),
            b * right,
            current * sa,
            0.0,
            -current * sb,
            -transfer.voltage_drop_v,
        ];
        if tangent.iter().any(|v| !v.is_finite()) {
            return Err(graph_error(
                "Kelvin graph consistent tangent is unrepresentable",
            ));
        }
        Ok((
            [
                transfer.voltage_drop_v,
                transfer.terminal_heat_w[0],
                transfer.bulk_heat_w,
                transfer.terminal_heat_w[1],
                transfer.electrical_supply_w,
            ],
            tangent,
        ))
    }
}

impl LawNode for MaterialSeebeckNode {
    fn declaration(&self) -> &NodeDeclaration {
        &self.declaration
    }
    fn evaluate(&self, state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        let (values, _) = self.response_at(state, inputs)?;
        Ok(NodeOutput {
            outputs: values.to_vec(),
            next_state: Vec::new(),
            dissipation_rate: Some(0.0),
        })
    }
    fn tangent(&self, state: &[f64], inputs: &[f64]) -> Option<Vec<f64>> {
        self.response_at(state, inputs)
            .ok()
            .map(|(_, tangent)| tangent.to_vec())
    }
}

fn invalid(reason: &'static str) -> ThermoelectricError {
    ThermoelectricError::Invalid(reason)
}
fn graph_error(obligation: &'static str) -> GraphError {
    GraphError::CardMismatch {
        law: LAW.into(),
        obligation,
    }
}
