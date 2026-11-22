//! Execute an admitted lossless N-port law in the ordinary component graph.

use core::{fmt, num::NonZeroUsize};

use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_exec::{Cancelled, Cx};
use fs_ports::component::{
    Class, Component, ContractViolation, PortRole, PortSpec, PortTangent, Proposal, Receipt,
};
use fs_ports::{JunctionLaw, JunctionSchema, PortValueShape};

use super::roundoff;

/// Why an admitted junction cannot run as a scalar or vector Dirac component.
#[derive(Debug, Clone, PartialEq)]
pub enum DiracError {
    /// The supplied schema describes common effort or common flow instead.
    NotDirac,
    /// A tensor or field needs an explicit coordinate or metric adapter.
    UnsupportedShape {
        /// Zero-based endpoint in the retained schema order.
        port: usize,
        /// The unsupported endpoint shape.
        shape: PortValueShape,
    },
    /// The complete port-major input/output extent cannot fit in `usize`.
    WidthOverflow {
        /// Number of endpoints.
        ports: usize,
        /// Coordinates at each endpoint.
        width: usize,
    },
    /// A port, step, vector extent or arithmetic result violates the contract.
    Contract(ContractViolation),
    /// The execution context requested cancellation.
    Cancelled,
}

impl fmt::Display for DiracError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotDirac => write!(f, "the junction schema does not describe a Dirac law"),
            Self::UnsupportedShape { port, shape } => write!(
                f,
                "Dirac endpoint {port} has unsupported shape {shape:?}; expected scalar or vector"
            ),
            Self::WidthOverflow { ports, width } => write!(
                f,
                "Dirac junction with {ports} endpoints of width {width} overflows its layout"
            ),
            Self::Contract(error) => write!(f, "{error}"),
            Self::Cancelled => write!(f, "cancelled"),
        }
    }
}

impl std::error::Error for DiracError {}

impl From<ContractViolation> for DiracError {
    fn from(error: ContractViolation) -> Self {
        Self::Contract(error)
    }
}

impl From<Cancelled> for DiracError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

/// A stateless, lossless relation with any admitted number of scalar or
/// compatible Euclidean vector endpoints, including singular and zero skew
/// matrices. Each vector coordinate obeys `flows = J * efforts`: the
/// port-major operator is `J tensor I`, with no mixing between coordinates.
///
/// Ports retain the schema's identifiers and order and receive effort.
/// The schema endpoints describe the attached owners: a flow outward from
/// an owner is inward to this junction. A [`crate::SystemBuilder::bond`]
/// therefore gives that owner minus the corresponding junction output.
/// The law adds neither storage, dissipation nor an integration rule.
///
/// All endpoints use the common SI coordinates admitted by the schema.
/// Its basis, frame and clock are retained, not converted. Ordinary runtime
/// power ports do not carry those labels: the caller must attach owners in
/// these same coordinates, using explicit power charts where needed. This
/// law does not flatten tensors or implement a field metric or cross-domain chart.
#[derive(Debug, Clone)]
pub struct DiracJunction {
    schema: JunctionSchema,
    ports: Vec<PortSpec>,
    width: usize,
    coordinates: usize,
}

impl DiracJunction {
    /// Execute an already admitted scalar or vector [`JunctionLaw::Dirac`].
    ///
    /// Structural admission (finite exact-skew matrix, distinct identities
    /// and compatible endpoint coordinates) belongs to `JunctionSchema`.
    /// Its endpoints have the same shape and pairing; vector endpoints use
    /// the Euclidean dot product in their declared common coordinates.
    ///
    /// # Errors
    /// Refuses another law, a tensor or field endpoint, or an overflowing
    /// total extent. No matrix is repaired, symmetrized, inverted or
    /// independently re-admitted here.
    pub fn new(schema: JunctionSchema) -> Result<Self, DiracError> {
        if !matches!(schema.law(), JunctionLaw::Dirac { .. }) {
            return Err(DiracError::NotDirac);
        }
        // An admitted schema has at least two endpoints, with identical
        // shapes. Check the complete extent before creating runtime ports.
        let width = match schema.ports()[0].shape() {
            PortValueShape::Scalar => NonZeroUsize::MIN,
            PortValueShape::Vector(width) => width,
            shape => return Err(DiracError::UnsupportedShape { port: 0, shape }),
        };
        let coordinates = schema.ports().len().checked_mul(width.get()).ok_or(
            DiracError::WidthOverflow {
                ports: schema.ports().len(),
                width: width.get(),
            },
        )?;
        let mut ports = Vec::with_capacity(schema.ports().len());
        for port in schema.ports() {
            ports.push(PortSpec::power(
                port.id().clone(),
                port.kind(),
                PortRole::EffortIn,
                width,
            )?);
        }
        Ok(Self {
            schema,
            ports,
            width: width.get(),
            coordinates,
        })
    }

    /// The original admitted relation, endpoint identities and coordinates.
    #[must_use]
    pub fn schema(&self) -> &JunctionSchema {
        &self.schema
    }

    fn matrix(&self) -> &[f64] {
        let JunctionLaw::Dirac { incidence, .. } = self.schema.law() else {
            unreachable!("construction requires a Dirac law");
        };
        incidence
    }

    fn extent(&self, vector: &'static str, actual: usize) -> Result<(), DiracError> {
        let expected = self.coordinates;
        if actual != expected {
            return Err(ContractViolation::VectorLength {
                vector,
                expected,
                actual,
            }
            .into());
        }
        Ok(())
    }

    fn validate(&self, dt_s: f64, inputs: &[f64], cx: &Cx<'_>) -> Result<(), DiracError> {
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(ContractViolation::Step { dt_s }.into());
        }
        self.extent("Dirac inputs", inputs.len())?;
        finite(inputs, "Dirac inputs", cx)
    }

    fn apply(
        &self,
        inputs: &[f64],
        outputs: &mut [f64],
        transpose: bool,
        cx: &Cx<'_>,
    ) -> Result<(), DiracError> {
        self.extent("Dirac operand", inputs.len())?;
        self.extent("Dirac result", outputs.len())?;
        finite(inputs, "Dirac operand", cx)?;
        for (row, output) in self
            .matrix()
            .chunks_exact(self.ports.len())
            .zip(outputs.chunks_exact_mut(self.width))
        {
            for (coordinate, value_out) in output.iter_mut().enumerate() {
                // One dot is bounded by the schema's endpoint ceiling,
                // regardless of the vector width.
                cx.checkpoint()?;
                let value: f64 = row
                    .iter()
                    .zip(inputs.chunks_exact(self.width))
                    .map(|(j, effort)| j * effort[coordinate])
                    .sum();
                finite(&[value], "Dirac result", cx)?;
                // Exact skew admission gives J^T = -J, including singular J.
                *value_out = if transpose { -value } else { value };
            }
        }
        Ok(())
    }
}

fn finite(values: &[f64], quantity: &'static str, cx: &Cx<'_>) -> Result<(), DiracError> {
    // Validate before writing output, while bounding cancellation latency
    // even when a vector endpoint has many coordinates.
    for chunk in values.chunks(1024) {
        cx.checkpoint()?;
        if chunk.iter().any(|value| !value.is_finite()) {
            return Err(ContractViolation::NonFinite { quantity }.into());
        }
    }
    Ok(())
}

impl Component for DiracJunction {
    type State = ();
    type Error = DiracError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::elements::DiracJunction",
        state_family: "no-state",
        schema_grammar: "(): an admitted lossless Dirac relation stores nothing",
        codec_grammar: "no words",
        codec_version: 1,
    });

    fn encode_state(&self, (): &(), _words: &mut Vec<u64>) -> bool {
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<()> {
        words.is_empty().then_some(())
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        // Runtime ports omit the schema's frame, basis, clock and audit
        // roles. Bind these too, using the system identity's same Debug
        // convention, and retain every incidence coefficient as exact bits.
        let metadata = format!("{:?} {:?}", self.schema.id(), self.schema.ports());
        words.push(metadata.len() as u64);
        for chunk in metadata.as_bytes().chunks(8) {
            let mut bytes = [0; 8];
            bytes[..chunk.len()].copy_from_slice(chunk);
            words.push(u64::from_le_bytes(bytes));
        }
        words.push(self.ports.len() as u64);
        words.extend(self.matrix().iter().map(|value| value.to_bits()));
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, _state: &()) -> Result<f64, DiracError> {
        Ok(0.0)
    }

    fn propose(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        _new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, DiracError> {
        cx.checkpoint()?;
        self.validate(dt_s, inputs, cx)?;
        self.apply(inputs, outputs, false, cx)?;
        let mut scale_j = 0.0;
        for (i, row) in self.matrix().chunks_exact(self.ports.len()).enumerate() {
            for coordinate in 0..self.width {
                cx.checkpoint()?;
                let own_effort = inputs[i * self.width + coordinate];
                for (j, effort) in row.iter().zip(inputs.chunks_exact(self.width)) {
                    scale_j += (own_effort * (j * effort[coordinate])).abs() * dt_s;
                }
            }
        }
        // Each dot has n products; power sums span all n * width
        // coordinates. Bound rounding by term magnitudes, never by the
        // cancelling net power. Width one retains the scalar arithmetic.
        let tolerance_j = roundoff((self.coordinates as f64 + 2.0) * scale_j);
        finite(&[tolerance_j], "Dirac energy allowance", cx)?;
        cx.checkpoint()?;
        Ok(Proposal::Step(Receipt::lossless(0.0, tolerance_j)))
    }
}

impl PortTangent for DiracJunction {
    fn jvp(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), DiracError> {
        cx.checkpoint()?;
        self.validate(dt_s, inputs, cx)?;
        self.apply(d_inputs, d_outputs, false, cx)
    }

    fn vjp(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), DiracError> {
        cx.checkpoint()?;
        self.validate(dt_s, inputs, cx)?;
        self.apply(bar_outputs, bar_inputs, true, cx)
    }
}
