//! Typed boundary conditions with dimensional analysis per physics.
//! Every (physics, kind) pair declares what dimensions its value must
//! carry — a Dirichlet velocity is not a Dirichlet temperature — and
//! inlets carrying flux declare their compatibility regime so net-flux
//! consistency can be verified at admission, not discovered as a diverged
//! solve.

use crate::ScenarioError;
use crate::scenario::Violation;
use crate::signal::{ChebProfile, TimeSignal};
use fs_qty::{Dims, QtyAny};

/// Which physics a boundary condition speaks to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Physics {
    /// Incompressible flow (velocity/pressure).
    IncompressibleFlow,
    /// Heat conduction/convection.
    Thermal,
    /// Solid elasticity.
    Elasticity,
}

/// The boundary-condition kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BcKind {
    /// Prescribed primary field value.
    Dirichlet,
    /// Prescribed flux.
    Neumann,
    /// Mixed (transfer-coefficient) condition.
    Robin,
    /// Mass-flow inlet (kg/s across the region).
    MassFlowInlet,
    /// Pressure outlet.
    PressureOutlet,
    /// No-slip wall (no value).
    WallNoSlip,
    /// Free-slip wall (no value).
    WallSlip,
    /// Prescribed traction (Pa).
    Traction,
}

/// The declared compatibility regime of a flux-carrying condition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compat {
    /// Net mass flux over the boundary must balance (or a pressure outlet
    /// must exist to absorb it).
    Incompressible,
}

/// A boundary value: uniform, time-varying, or a spatial profile ("the
/// inlet profile is THIS FUNCTION").
#[derive(Debug, Clone, PartialEq)]
pub enum BcValue {
    /// One dimensioned value.
    Uniform(QtyAny),
    /// A time history.
    Signal(TimeSignal),
    /// A chebfun spatial profile over the region's parameter.
    Profile(ChebProfile),
}

impl BcValue {
    /// The dimensions this value carries.
    #[must_use]
    pub fn dims(&self) -> Dims {
        match self {
            BcValue::Uniform(q) => q.dims,
            BcValue::Signal(s) => s.dims(),
            BcValue::Profile(p) => p.dims,
        }
    }
}

/// One typed boundary condition attached to a named region patch.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundaryCondition {
    /// The region/patch name (resolved against fs-geom regions upstream).
    pub region: String,
    /// Which physics this condition constrains.
    pub physics: Physics,
    /// The condition kind.
    pub kind: BcKind,
    /// The value (None only for kinds that forbid one).
    pub value: Option<BcValue>,
    /// Declared compatibility regime (required for flux-carrying inlets).
    pub compatibility: Option<Compat>,
    /// The frame this condition's vector data is expressed in (0 = world).
    pub frame: u32,
}

/// What a (physics, kind) pair demands of its value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Expectation {
    /// A value with exactly these SI exponents.
    Value(Dims),
    /// No value at all (geometric/no-slip kinds).
    NoValue,
    /// The pair is not part of this physics' vocabulary.
    Unsupported,
}

/// SI exponents (m, kg, s, K, A, mol) for the common quantities.
mod dims {
    use fs_qty::Dims;
    pub const VELOCITY: Dims = Dims([1, 0, -1, 0, 0, 0]);
    pub const MASS_FLOW: Dims = Dims([0, 1, -1, 0, 0, 0]);
    pub const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
    pub const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
    pub const HEAT_FLUX: Dims = Dims([0, 1, -3, 0, 0, 0]);
    pub const HTC: Dims = Dims([0, 1, -3, -1, 0, 0]);
    pub const DISPLACEMENT: Dims = Dims([1, 0, 0, 0, 0, 0]);
}

/// The dimensional contract of every supported (physics, kind) pair.
#[must_use]
pub fn expectation(physics: Physics, kind: BcKind) -> Expectation {
    use Expectation::{NoValue, Unsupported, Value};
    match (physics, kind) {
        (Physics::IncompressibleFlow, BcKind::Dirichlet) => Value(dims::VELOCITY),
        (Physics::IncompressibleFlow, BcKind::MassFlowInlet) => Value(dims::MASS_FLOW),
        (Physics::IncompressibleFlow, BcKind::PressureOutlet | BcKind::Traction) => {
            Value(dims::PRESSURE)
        }
        (Physics::IncompressibleFlow, BcKind::WallNoSlip | BcKind::WallSlip) => NoValue,
        (Physics::Thermal, BcKind::Dirichlet) => Value(dims::TEMPERATURE),
        (Physics::Thermal, BcKind::Neumann) => Value(dims::HEAT_FLUX),
        (Physics::Thermal, BcKind::Robin) => Value(dims::HTC),
        (Physics::Elasticity, BcKind::Dirichlet) => Value(dims::DISPLACEMENT),
        (Physics::Elasticity, BcKind::Traction) => Value(dims::PRESSURE),
        _ => Unsupported,
    }
}

impl BoundaryCondition {
    /// Validate this condition against the dimensional contract.
    pub fn check(&self, out: &mut Vec<Violation>) {
        let ctx = format!(
            "bc on {:?} ({:?}/{:?})",
            self.region, self.physics, self.kind
        );
        if let Some(value) = &self.value {
            match value {
                BcValue::Uniform(quantity) => {
                    if !quantity.value.is_finite() {
                        out.push(Violation {
                            code: "bc-value-nonfinite",
                            what: format!("{ctx}: uniform value {} is non-finite", quantity.value),
                            fix: "replace the boundary value with a finite value".to_string(),
                        });
                    }
                }
                BcValue::Signal(signal) => signal.check(&ctx, out),
                BcValue::Profile(profile) => profile.check(&ctx, out),
            }
        }
        match expectation(self.physics, self.kind) {
            Expectation::Unsupported => out.push(Violation {
                code: "bc-kind-unsupported",
                what: format!("{ctx}: this kind is not in this physics' vocabulary"),
                fix: "use a supported kind for the physics (see fs_scenario::bc::expectation)"
                    .to_string(),
            }),
            Expectation::NoValue => {
                if self.value.is_some() {
                    out.push(Violation {
                        code: "bc-value-forbidden",
                        what: format!("{ctx}: geometric condition carries a value"),
                        fix: "remove the value; no-slip/slip walls take none".to_string(),
                    });
                }
            }
            Expectation::Value(expected) => match &self.value {
                None => out.push(Violation {
                    code: "bc-value-missing",
                    what: format!("{ctx}: no value supplied"),
                    fix: format!(
                        "supply a value with SI exponents {:?} (m, kg, s, K, A, mol)",
                        expected.0
                    ),
                }),
                Some(v) => {
                    if v.dims() != expected {
                        out.push(Violation {
                            code: "bc-dims",
                            what: format!(
                                "{ctx}: value has dimensions {:?}, contract demands {:?}",
                                v.dims().0,
                                expected.0
                            ),
                            fix: "express the value in the quantity the contract demands"
                                .to_string(),
                        });
                    }
                }
            },
        }
        if self.kind == BcKind::MassFlowInlet && self.compatibility.is_none() {
            out.push(Violation {
                code: "bc-compat-missing",
                what: format!("{ctx}: flux-carrying inlet declares no compatibility regime"),
                fix: "declare `compatibility: incompressible` so net-flux \
                      consistency can be verified at admission"
                    .to_string(),
            });
        }
        if self.kind != BcKind::MassFlowInlet && self.compatibility.is_some() {
            out.push(Violation {
                code: "bc-compat-forbidden",
                what: format!("{ctx}: only mass-flow inlets may declare a compatibility regime"),
                fix: "remove compatibility from this boundary condition".to_string(),
            });
        }
        if self.kind == BcKind::MassFlowInlet && matches!(&self.value, Some(BcValue::Profile(_))) {
            out.push(Violation {
                code: "bc-mass-flow-profile",
                what: format!(
                    "{ctx}: a spatial profile cannot yet certify a total mass-flow contribution"
                ),
                fix: "supply a uniform or time-signal total in kg/s; velocity-profile surface integration belongs at the geometry-bound solver boundary"
                    .to_string(),
            });
        }
    }

    /// Signed mass-flow contribution of this condition at time `t`
    /// (inlets positive, pressure outlets are flux-free), used by the
    /// net-flux compatibility check.
    ///
    /// # Errors
    /// Returns [`ScenarioError`] when a declared total-flow value has the wrong
    /// dimensions, is non-finite, cannot be evaluated, or is a spatial profile
    /// for which this layer has no geometry-backed surface integral.
    pub fn mass_flow_at(&self, t: f64) -> Result<Option<f64>, ScenarioError> {
        if self.kind != BcKind::MassFlowInlet {
            return Ok(None);
        }
        if self.physics != Physics::IncompressibleFlow {
            return Err(ScenarioError::Evaluate {
                what: format!(
                    "mass-flow inlet on {:?} is attached to unsupported physics {:?}",
                    self.region, self.physics
                ),
            });
        }
        if !t.is_finite() {
            return Err(ScenarioError::Evaluate {
                what: format!("mass-flow evaluation time {t} is non-finite"),
            });
        }
        match &self.value {
            Some(BcValue::Uniform(quantity)) => {
                if quantity.dims != dims::MASS_FLOW {
                    return Err(ScenarioError::Dimensions {
                        context: format!("mass-flow inlet on {:?}", self.region),
                        expected: dims::MASS_FLOW.0,
                        got: quantity.dims.0,
                    });
                }
                if !quantity.value.is_finite() {
                    return Err(ScenarioError::Evaluate {
                        what: format!(
                            "mass-flow inlet on {:?} has non-finite value {}",
                            self.region, quantity.value
                        ),
                    });
                }
                Ok(Some(quantity.value))
            }
            Some(BcValue::Signal(signal)) => {
                let quantity = signal.eval(t)?;
                if quantity.dims != dims::MASS_FLOW {
                    return Err(ScenarioError::Dimensions {
                        context: format!("mass-flow inlet signal on {:?}", self.region),
                        expected: dims::MASS_FLOW.0,
                        got: quantity.dims.0,
                    });
                }
                Ok(Some(quantity.value))
            }
            Some(BcValue::Profile(_)) => Err(ScenarioError::Evaluate {
                what: format!(
                    "mass-flow inlet on {:?} uses a spatial profile without a geometry-backed surface integral",
                    self.region
                ),
            }),
            None => Err(ScenarioError::Evaluate {
                what: format!("mass-flow inlet on {:?} has no declared value", self.region),
            }),
        }
    }
}
