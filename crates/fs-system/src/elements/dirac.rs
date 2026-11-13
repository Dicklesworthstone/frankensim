//! Execute an admitted lossless N-port law in the ordinary component graph.

use core::{fmt, num::NonZeroUsize};

use fs_exec::{Cancelled, Cx};
use fs_ports::component::{
    Class, Component, ContractViolation, PortRole, PortSpec, PortTangent, Proposal, Receipt,
};
use fs_ports::{JunctionLaw, JunctionSchema, PortValueShape};

use super::roundoff;

/// Why an admitted junction cannot run as a scalar Dirac component.
#[derive(Debug, Clone, PartialEq)]
pub enum DiracError {
    /// The supplied schema describes common effort or common flow instead.
    NotDirac,
    /// This endpoint is not scalar; no field or coordinate map is inferred.
    NonScalar {
        /// Zero-based endpoint in the retained schema order.
        port: usize,
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
            Self::NonScalar { port } => write!(f, "Dirac endpoint {port} is not scalar"),
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

/// A stateless, lossless relation `flows = J * efforts` with any admitted
/// number of scalar endpoints, including singular and zero skew matrices.
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
/// scalar law does not implement a field metric or cross-domain chart.
#[derive(Debug, Clone)]
pub struct DiracJunction {
    schema: JunctionSchema,
    ports: Vec<PortSpec>,
}

impl DiracJunction {
    /// Execute an already admitted scalar [`JunctionLaw::Dirac`].
    ///
    /// Structural admission (finite exact-skew matrix, distinct identities
    /// and compatible endpoint coordinates) belongs to `JunctionSchema`.
    ///
    /// # Errors
    /// Refuses another law or a non-scalar endpoint. No matrix is repaired,
    /// symmetrized, inverted or independently re-admitted here.
    pub fn new(schema: JunctionSchema) -> Result<Self, DiracError> {
        if !matches!(schema.law(), JunctionLaw::Dirac { .. }) {
            return Err(DiracError::NotDirac);
        }
        let mut ports = Vec::with_capacity(schema.ports().len());
        for (index, port) in schema.ports().iter().enumerate() {
            if port.shape() != PortValueShape::Scalar {
                return Err(DiracError::NonScalar { port: index });
            }
            ports.push(PortSpec::power(
                port.id().clone(),
                port.kind(),
                PortRole::EffortIn,
                NonZeroUsize::MIN,
            )?);
        }
        Ok(Self { schema, ports })
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
        let expected = self.ports.len();
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

    fn validate(&self, dt_s: f64, inputs: &[f64]) -> Result<(), DiracError> {
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(ContractViolation::Step { dt_s }.into());
        }
        self.extent("Dirac inputs", inputs.len())?;
        finite(inputs, "Dirac inputs")
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
        finite(inputs, "Dirac operand")?;
        for (row, output) in self.matrix().chunks_exact(self.ports.len()).zip(outputs) {
            cx.checkpoint()?;
            let value: f64 = row.iter().zip(inputs).map(|(j, effort)| j * effort).sum();
            finite(&[value], "Dirac result")?;
            // Exact skew admission gives J^T = -J, including singular J.
            *output = if transpose { -value } else { value };
        }
        Ok(())
    }
}

fn finite(values: &[f64], quantity: &'static str) -> Result<(), DiracError> {
    if values.iter().any(|value| !value.is_finite()) {
        return Err(ContractViolation::NonFinite { quantity }.into());
    }
    Ok(())
}

impl Component for DiracJunction {
    type State = ();
    type Error = DiracError;

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
        self.validate(dt_s, inputs)?;
        self.apply(inputs, outputs, false, cx)?;
        let mut scale_j = 0.0;
        for (i, row) in self.matrix().chunks_exact(self.ports.len()).enumerate() {
            cx.checkpoint()?;
            for (j, effort) in row.iter().zip(inputs) {
                scale_j += (inputs[i] * (j * effort)).abs() * dt_s;
            }
        }
        // Each dot has n products; power sums add n more terms. Bound
        // rounding by term magnitudes, never by the cancelling net power.
        let tolerance_j = roundoff((self.ports.len() as f64 + 2.0) * scale_j);
        finite(&[tolerance_j], "Dirac energy allowance")?;
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
        self.validate(dt_s, inputs)?;
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
        self.validate(dt_s, inputs)?;
        self.apply(bar_outputs, bar_inputs, true, cx)
    }
}
