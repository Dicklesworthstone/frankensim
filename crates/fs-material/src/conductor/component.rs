//! A material card's Ohmic conductor behind the composition contract's
//! resistor ports.
//!
//! Two electrical terminals and a heat port. Every proposal resolves the
//! card's resistance `R(T) = rho(T) L / A` again at the temperature the heat
//! port receives ([`MaterialOhmicNode`]), and all of the Joule heat leaves
//! through that port. A connected thermal owner therefore moves the
//! resistance with no pairing-specific loop. A temperature outside the node's
//! admitted interval is the node's refusal, never an extrapolation.
//!
//! This is memoryless DC conduction: no inductance, capacitance, skin effect,
//! contact resistance or heat storage of its own. The geometry is the node's.
//!
//! Refusals are the node's own [`ConductorError`], unchanged. A port value,
//! timestep or slot count the resistor cannot take, and a cancelled caller,
//! are `ConductorError::InvalidInput` naming what was refused: the crate
//! adds no error type of its own here (the enumeration ratchet counts them).

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::component::{
    Class, Component, EndOfStep, PortRole, PortSpec, PortTangent, Proposal, Receipt,
};
use fs_ports::{PortKind, StableId};

use super::graph::MaterialOhmicNode;
use super::{ConductorError, OhmicDrive};

/// Which of the two terminal quantities the resistor receives.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Terminals {
    /// Currents in at `a` and `b`, voltages out. The current through the
    /// conductor is `i_a - i_b`: the form of a resistor in a current loop.
    Current,
    /// Voltages in at `a` and `b`, currents out. The voltage across the
    /// conductor is `v_a - v_b`: the form of a branch between two nodes.
    Voltage,
}

/// What the resistor itself refuses.
const fn refused(quantity: &'static str) -> ConductorError {
    ConductorError::InvalidInput { quantity }
}

fn poll(cx: &Cx<'_>) -> Result<(), ConductorError> {
    cx.checkpoint()
        .map_err(|_| refused("Ohmic resistor cancelled"))
}

fn values(values: &[f64]) -> Result<[f64; 3], ConductorError> {
    match values {
        [a, b, c] if values.iter().all(|value| value.is_finite()) => Ok([*a, *b, *c]),
        _ => Err(refused("Ohmic resistor needs three finite port values")),
    }
}

/// The conductor's answer at one set of port values.
struct Response {
    /// The output at `a`: a voltage for current terminals, a current for
    /// voltage terminals. The output at `b` is its negative.
    output: f64,
    /// The current through the conductor [A].
    current_a: f64,
    /// The Joule power `I^2 R` [W], the node's own.
    power_w: f64,
    resistance_ohm: f64,
    /// `dR/dT` [ohm/K] at the heat port's temperature.
    slope_ohm_k: f64,
}

/// A source-owned resistor whose resistance follows its heat port.
///
/// Inputs, in port order, are the terminal currents or voltages at `a` and
/// `b` ([`Terminals`]) and the conductor's temperature at `heat` [K]. Outputs
/// are the complementary terminal quantities and the negative Joule power
/// [W]. The state is `()`: all heat belongs to the connected thermal owner.
#[derive(Clone, Debug)]
pub struct OhmicResistor {
    node: MaterialOhmicNode,
    terminals: Terminals,
    ports: [PortSpec; 3],
}

impl OhmicResistor {
    /// Put `node` behind electrical ports of the form `terminals` and a heat
    /// port.
    ///
    /// # Errors
    /// A node bound to a liquid mass fraction: the heat port carries a
    /// temperature only.
    pub fn new(node: MaterialOhmicNode, terminals: Terminals) -> Result<Self, ConductorError> {
        if node.liquid_mass_fraction_axis().is_some() {
            return Err(refused(
                "Ohmic resistor node bound to a liquid mass fraction (the heat port carries a temperature only)",
            ));
        }
        let role = match terminals {
            Terminals::Current => PortRole::FlowIn,
            Terminals::Voltage => PortRole::EffortIn,
        };
        let electrical = |name| {
            PortSpec::power(
                StableId::new(name).expect("static port name"),
                PortKind::ElectricalVoltageCurrent,
                role,
                NonZeroUsize::MIN,
            )
            .expect("scalar electrical power port")
        };
        Ok(Self {
            node,
            terminals,
            ports: [
                electrical("a"),
                electrical("b"),
                PortSpec::heat(
                    StableId::new("heat").expect("static port name"),
                    PortRole::EffortIn,
                )
                .expect("thermal outlet"),
            ],
        })
    }

    /// The unchanged node: card, source, temperature interval and geometry.
    #[must_use]
    pub const fn node(&self) -> &MaterialOhmicNode {
        &self.node
    }

    /// The terminal form.
    #[must_use]
    pub const fn terminals(&self) -> Terminals {
        self.terminals
    }

    fn response(&self, inputs: &[f64]) -> Result<Response, ConductorError> {
        let [a, b, kelvin] = values(inputs)?;
        let (conductor, slope_ohm_k) = self.node.resolve_and_slope(kelvin)?;
        let resistance_ohm = conductor.resistance_ohm();
        let difference = a - b;
        let (current_a, output) = match self.terminals {
            // The node's own voltage, `I R`.
            Terminals::Current => (difference, difference * resistance_ohm),
            Terminals::Voltage => {
                let current = conductor.current_a_for_drive(OhmicDrive::Voltage(difference))?;
                (current, current)
            }
        };
        let power_w = conductor.joule_power_w(current_a)?;
        if !output.is_finite() {
            return Err(refused("Ohmic resistor terminal value is unrepresentable"));
        }
        Ok(Response {
            output,
            current_a,
            power_w,
            resistance_ohm,
            slope_ohm_k,
        })
    }

    /// `d(outputs)/d(inputs)`, rows `a`, `b`, `heat` and columns the same.
    fn jacobian(&self, inputs: &[f64]) -> Result<[[f64; 3]; 3], ConductorError> {
        let response = self.response(inputs)?;
        let (r, slope, current) = (
            response.resistance_ohm,
            response.slope_ohm_k,
            response.current_a,
        );
        let matrix = match self.terminals {
            // V = R I with I = i_a - i_b; P = I^2 R.
            Terminals::Current => {
                let voltage = response.output;
                [
                    [r, -r, current * slope],
                    [-r, r, -current * slope],
                    [-2.0 * voltage, 2.0 * voltage, -current * current * slope],
                ]
            }
            // I = (v_a - v_b) / R; P = (v_a - v_b)^2 / R.
            Terminals::Voltage => {
                let conductance = 1.0 / r;
                let thermal = -current * slope / r;
                [
                    [conductance, -conductance, thermal],
                    [-conductance, conductance, -thermal],
                    [-2.0 * current, 2.0 * current, current * current * slope],
                ]
            }
        };
        if matrix.iter().flatten().any(|value| !value.is_finite()) {
            return Err(refused("Ohmic resistor tangent is unrepresentable"));
        }
        Ok(matrix)
    }
}

impl Component for OhmicResistor {
    type State = ();
    type Error = ConductorError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-material::conductor::component::OhmicResistor",
        state_family: "stateless-ohmic-resistor",
        schema_grammar: "unit state on the original Ohmic node, its card, interval and geometry",
        codec_grammar: "zero u64 words",
        codec_version: 1,
    });

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, _state: &()) -> Result<f64, ConductorError> {
        Ok(0.0)
    }

    fn encode_state(&self, _state: &(), _words: &mut Vec<u64>) -> bool {
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<()> {
        words.is_empty().then_some(())
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        // The card, the selected claim and the query it was resolved at; the
        // interval, the reference temperature and the geometry. The
        // temperature of each step is not bound: every proposal resolves the
        // card again at the heat port's.
        let source = self.node.source();
        let hash = |bytes: &[u8; 32]| {
            bytes
                .as_chunks::<8>()
                .0
                .iter()
                .map(|chunk| u64::from_le_bytes(*chunk))
                .collect::<Vec<_>>()
        };
        words.extend(hash(source.card().content_hash().as_bytes()));
        let receipt = &source.original().source_property().answer().receipt;
        words.extend(hash(receipt.content_hash().as_bytes()));
        let axis = self.node.temperature_axis().as_bytes();
        words.push(axis.len() as u64);
        words.extend(axis.iter().map(|&byte| u64::from(byte)));
        words.push(self.node.reference_temperature_k().to_bits());
        words.extend(self.node.temperature_range_k().map(f64::to_bits));
        words.extend([
            self.node.length_m().to_bits(),
            self.node.area_m2().to_bits(),
            u64::from(self.node.free_expansion().is_some()),
            match self.terminals {
                Terminals::Current => 0,
                Terminals::Voltage => 1,
            },
        ]);
    }

    fn end_of_step(&self, _port: usize) -> EndOfStep {
        // An algebraic, held-input response: no flow history from which an
        // endpoint relation could be reconstructed.
        EndOfStep::Untold
    }

    fn propose(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ConductorError> {
        poll(cx)?;
        if outputs.len() != 3 {
            return Err(refused("Ohmic resistor needs three output slots"));
        }
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(refused(
                "Ohmic resistor timestep must be positive and finite",
            ));
        }
        let response = self.response(inputs)?;
        let kelvin = inputs[2];
        let dissipated_j = dt_s * response.power_w;
        let energy_tolerance_j = 64.0
            * f64::EPSILON
            * (dt_s * (response.output * inputs[0]).abs()
                + dt_s * (response.output * inputs[1]).abs()
                + dissipated_j);
        let entropy = dissipated_j / kelvin;
        if ![dissipated_j, energy_tolerance_j, entropy]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(refused("Ohmic resistor energy account is unrepresentable"));
        }
        poll(cx)?;
        *new = ();
        outputs.copy_from_slice(&[response.output, -response.output, -response.power_w]);
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: entropy,
            energy_tolerance_j,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

impl PortTangent for OhmicResistor {
    fn jvp(
        &self,
        _old: &(),
        _dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ConductorError> {
        poll(cx)?;
        if d_outputs.len() != 3 {
            return Err(refused("Ohmic resistor needs three tangent output slots"));
        }
        let direction = values(d_inputs)?;
        let matrix = self.jacobian(inputs)?;
        let mut next = [0.0; 3];
        for (value, row) in next.iter_mut().zip(matrix) {
            *value = row.iter().zip(direction).map(|(j, d)| j * d).sum();
        }
        d_outputs.copy_from_slice(&next);
        Ok(())
    }

    fn vjp(
        &self,
        _old: &(),
        _dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ConductorError> {
        poll(cx)?;
        if bar_inputs.len() != 3 {
            return Err(refused("Ohmic resistor needs three tangent input slots"));
        }
        let direction = values(bar_outputs)?;
        let matrix = self.jacobian(inputs)?;
        let mut next = [0.0; 3];
        for (column, value) in next.iter_mut().enumerate() {
            *value = matrix
                .iter()
                .zip(direction)
                .map(|(row, d)| row[column] * d)
                .sum();
        }
        bar_inputs.copy_from_slice(&next);
        Ok(())
    }
}
