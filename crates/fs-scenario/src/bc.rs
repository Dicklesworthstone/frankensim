//! Typed boundary conditions with dimensional analysis per physics.
//! Every (physics, kind) pair declares what dimensions its value must
//! carry — a Dirichlet velocity is not a Dirichlet temperature — and
//! inlets carrying flux declare their compatibility regime so net-flux
//! consistency can be verified at admission, not discovered as a diverged
//! solve.

use crate::ScenarioError;
use crate::payload::{Payload, PayloadKind};
use crate::scenario::Violation;
use crate::signal::{ChebProfile, TimeSignal};
use fs_qty::{Dims, QtyAny};
use std::fmt;

/// Which physics a boundary condition speaks to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Physics {
    /// Incompressible flow (velocity/pressure).
    IncompressibleFlow,
    /// Heat conduction/convection.
    Thermal,
    /// Solid elasticity.
    Elasticity,
    /// Magnetostatic or low-frequency magnetic fields.
    Magnetics,
    /// Electrostatic or conduction-current fields.
    Electrics,
    /// Multispecies gas exchange and characteristic boundaries.
    GasExchange,
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
    /// Magnetic vector potential (Wb/m), expressed as a vector payload.
    MagneticVectorPotential,
    /// Boundary-normal magnetic flux density (T).
    NormalMagneticFluxDensity,
    /// Electric potential (V).
    ElectricPotential,
    /// Boundary-normal current density (A/m²).
    NormalCurrentDensity,
    /// Species amount flux (mol/(m² s)).
    SpeciesAmountFlux,
    /// Species mass flux (kg/(m² s)).
    SpeciesMassFlux,
    /// Incoming gas characteristic state.
    GasCharacteristicInlet,
    /// Outgoing gas characteristic state.
    GasCharacteristicOutlet,
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
    /// A versioned payload carrying explicit kind, basis, frame, orientation,
    /// and continuity/reset semantics.
    Typed(Payload),
}

impl BcValue {
    /// Homogeneous dimensions carried by this value, or `None` for a
    /// heterogeneous characteristic payload.
    #[must_use]
    pub fn homogeneous_dims(&self) -> Option<Dims> {
        match self {
            BcValue::Uniform(q) => Some(q.dims),
            BcValue::Signal(s) => Some(s.dims()),
            BcValue::Profile(p) => Some(p.dims),
            BcValue::Typed(payload) => payload.homogeneous_dims(),
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

#[derive(Clone, Copy)]
struct BoundaryDiagnosticContext<'a> {
    region: &'a str,
    physics: Physics,
    kind: BcKind,
}

impl fmt::Display for BoundaryDiagnosticContext<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "bc on {:?} ({:?}/{:?})",
            self.region, self.physics, self.kind
        )
    }
}

/// What a (physics, kind) pair demands of its value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Expectation {
    /// A value with exactly these SI exponents.
    Value(Dims),
    /// A versioned payload with an exact structural kind and, when present,
    /// exact homogeneous dimensions.
    Typed {
        /// Required payload carrier.
        kind: PayloadKind,
        /// Required dimensions; `None` selects a heterogeneous carrier.
        dims: Option<Dims>,
    },
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
    pub const MAGNETIC_VECTOR_POTENTIAL: Dims = Dims([1, 1, -2, 0, -1, 0]);
    pub const MAGNETIC_FLUX_DENSITY: Dims = Dims([0, 1, -2, 0, -1, 0]);
    pub const ELECTRIC_POTENTIAL: Dims = Dims([2, 1, -3, 0, -1, 0]);
    pub const CURRENT_DENSITY: Dims = Dims([-2, 0, 0, 0, 1, 0]);
    pub const SPECIES_AMOUNT_FLUX: Dims = Dims([-2, 0, -1, 0, 0, 1]);
    pub const SPECIES_MASS_FLUX: Dims = Dims([-2, 1, -1, 0, 0, 0]);
}

/// The dimensional contract of every supported (physics, kind) pair.
#[must_use]
pub fn expectation(physics: Physics, kind: BcKind) -> Expectation {
    use Expectation::{NoValue, Typed, Unsupported, Value};
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
        (Physics::Magnetics, BcKind::MagneticVectorPotential) => Typed {
            kind: PayloadKind::Vector,
            dims: Some(dims::MAGNETIC_VECTOR_POTENTIAL),
        },
        (Physics::Magnetics, BcKind::NormalMagneticFluxDensity) => Typed {
            kind: PayloadKind::Scalar,
            dims: Some(dims::MAGNETIC_FLUX_DENSITY),
        },
        (Physics::Electrics, BcKind::ElectricPotential) => Typed {
            kind: PayloadKind::Scalar,
            dims: Some(dims::ELECTRIC_POTENTIAL),
        },
        (Physics::Electrics, BcKind::NormalCurrentDensity) => Typed {
            kind: PayloadKind::Scalar,
            dims: Some(dims::CURRENT_DENSITY),
        },
        (Physics::GasExchange, BcKind::SpeciesAmountFlux) => Typed {
            kind: PayloadKind::SpeciesBundle,
            dims: Some(dims::SPECIES_AMOUNT_FLUX),
        },
        (Physics::GasExchange, BcKind::SpeciesMassFlux) => Typed {
            kind: PayloadKind::SpeciesBundle,
            dims: Some(dims::SPECIES_MASS_FLUX),
        },
        (
            Physics::GasExchange,
            BcKind::GasCharacteristicInlet | BcKind::GasCharacteristicOutlet,
        ) => Typed {
            kind: PayloadKind::CharacteristicState,
            dims: None,
        },
        _ => Unsupported,
    }
}

impl BoundaryCondition {
    /// Validate this condition against the dimensional contract.
    pub fn check(&self, out: &mut Vec<Violation>) {
        let mut checkpoint = |_: &'static str| Ok::<(), core::convert::Infallible>(());
        match self.check_with_checkpoint(out, &mut checkpoint) {
            Ok(()) => {}
            Err(never) => match never {},
        }
    }

    pub(crate) fn check_with_checkpoint<E>(
        &self,
        out: &mut Vec<Violation>,
        checkpoint: &mut impl FnMut(&'static str) -> Result<(), E>,
    ) -> Result<(), E> {
        let ctx = BoundaryDiagnosticContext {
            region: self.region.as_str(),
            physics: self.physics,
            kind: self.kind,
        };
        if self.region.is_empty() {
            out.push(Violation {
                code: "bc-region-empty",
                what: format!("{ctx}: region identity is empty"),
                fix: "bind the condition to a nonempty exact UTF-8 region name".to_string(),
            });
        }
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
                BcValue::Signal(signal) => {
                    signal.check_with_checkpoint(&ctx, out, checkpoint)?;
                }
                BcValue::Profile(profile) => {
                    profile.check_with_checkpoint(&ctx, out, checkpoint)?;
                }
                BcValue::Typed(_) => {}
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
                Some(BcValue::Typed(payload)) => out.push(Violation {
                    code: if self.kind == BcKind::MassFlowInlet {
                        "bc-mass-flow-typed"
                    } else {
                        "bc-legacy-value-required"
                    },
                    what: format!(
                        "{ctx}: typed {:?} payload cannot stand in for this legacy scalar/signal/profile value",
                        payload.kind()
                    ),
                    fix: format!(
                        "supply a legacy value with SI exponents {:?}; typed payloads are accepted only by typed expectation rows",
                        expected.0
                    ),
                }),
                Some(v) => {
                    if v.homogeneous_dims() != Some(expected) {
                        out.push(Violation {
                            code: "bc-dims",
                            what: format!(
                                "{ctx}: value has dimensions {:?}, contract demands {:?}",
                                v.homogeneous_dims().map(|dims| dims.0),
                                expected.0
                            ),
                            fix: "express the value in the quantity the contract demands"
                                .to_string(),
                        });
                    }
                }
            },
            Expectation::Typed {
                kind,
                dims: expected_dims,
            } => match &self.value {
                None => out.push(Violation {
                    code: "bc-value-missing",
                    what: format!("{ctx}: no typed payload supplied"),
                    fix: format!(
                        "supply a {kind:?} payload with homogeneous dimensions {expected_dims:?}"
                    ),
                }),
                Some(BcValue::Typed(payload)) => {
                    if payload.kind() != kind {
                        out.push(Violation {
                            code: "bc-payload-kind",
                            what: format!(
                                "{ctx}: payload kind {:?} does not match required {kind:?}",
                                payload.kind()
                            ),
                            fix: format!("supply a {kind:?} payload for this boundary kind"),
                        });
                    }
                    if payload.homogeneous_dims() != expected_dims {
                        out.push(Violation {
                            code: "bc-payload-dims",
                            what: format!(
                                "{ctx}: payload dimensions {:?} do not match required {:?}",
                                payload.homogeneous_dims().map(|dims| dims.0),
                                expected_dims.map(|dims| dims.0)
                            ),
                            fix: "construct the payload with the exact six-base quantity contract required by the expectation row"
                                .to_string(),
                        });
                    }
                    if payload.meta().frame().0 != self.frame {
                        out.push(Violation {
                            code: "bc-payload-frame",
                            what: format!(
                                "{ctx}: payload frame {} does not match boundary frame {}",
                                payload.meta().frame().0,
                                self.frame
                            ),
                            fix: "construct the payload and boundary condition with the same declared frame id"
                                .to_string(),
                        });
                    }
                }
                Some(_) => out.push(Violation {
                    code: "bc-typed-payload-required",
                    what: format!(
                        "{ctx}: this boundary row requires a versioned {kind:?} payload"
                    ),
                    fix: "replace the legacy scalar/signal/profile value with BcValue::Typed"
                        .to_string(),
                }),
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
        Ok(())
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
        self.mass_flow_at_impl(t, false)
    }

    /// Evaluate after whole-scenario validation has already scanned this
    /// boundary condition's dynamic signal payload.
    pub(crate) fn mass_flow_at_prevalidated(&self, t: f64) -> Result<Option<f64>, ScenarioError> {
        self.mass_flow_at_impl(t, true)
    }

    /// Prevalidated mass-flow evaluation with bounded Chebyshev cancellation.
    /// The outer result reports checkpoint refusal while the inner result
    /// retains the established evaluation errors.
    pub(crate) fn mass_flow_at_prevalidated_with_checkpoint<E>(
        &self,
        t: f64,
        checkpoint: &mut impl FnMut() -> Result<(), E>,
    ) -> Result<Result<Option<f64>, ScenarioError>, E> {
        let signal = match &self.value {
            Some(BcValue::Signal(signal @ TimeSignal::Chebfun(_))) => signal,
            _ => return Ok(self.mass_flow_at_prevalidated(t)),
        };
        if self.kind != BcKind::MassFlowInlet {
            return Ok(Ok(None));
        }
        if self.physics != Physics::IncompressibleFlow {
            return Ok(Err(ScenarioError::Evaluate {
                what: format!(
                    "mass-flow inlet on {:?} is attached to unsupported physics {:?}",
                    self.region, self.physics
                ),
            }));
        }
        if !t.is_finite() {
            return Ok(Err(ScenarioError::Evaluate {
                what: format!("mass-flow evaluation time {t} is non-finite"),
            }));
        }
        let quantity = match signal.eval_prevalidated_with_checkpoint(t, checkpoint)? {
            Ok(quantity) => quantity,
            Err(error) => return Ok(Err(error)),
        };
        if quantity.dims != dims::MASS_FLOW {
            return Ok(Err(ScenarioError::Dimensions {
                context: format!("mass-flow inlet signal on {:?}", self.region),
                expected: dims::MASS_FLOW.0,
                got: quantity.dims.0,
            }));
        }
        Ok(Ok(Some(quantity.value)))
    }

    fn mass_flow_at_impl(
        &self,
        t: f64,
        signal_prevalidated: bool,
    ) -> Result<Option<f64>, ScenarioError> {
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
                let quantity = if signal_prevalidated {
                    signal.eval_prevalidated(t)?
                } else {
                    signal.eval(t)?
                };
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
            Some(BcValue::Typed(payload)) => Err(ScenarioError::Evaluate {
                what: format!(
                    "mass-flow inlet on {:?} uses typed {:?} payload data, not an evaluable declared total in kg/s",
                    self.region,
                    payload.kind()
                ),
            }),
            None => Err(ScenarioError::Evaluate {
                what: format!("mass-flow inlet on {:?} has no declared value", self.region),
            }),
        }
    }
}

#[cfg(test)]
mod validation_internal_tests {
    use super::{BcKind, BoundaryCondition, BoundaryDiagnosticContext, Compat, Physics, dims};
    use crate::bc::BcValue;
    use fs_qty::QtyAny;

    #[test]
    fn diagnostic_context_is_borrowed_and_output_stable() {
        let region = String::from("inlet");
        let context = BoundaryDiagnosticContext {
            region: region.as_str(),
            physics: Physics::IncompressibleFlow,
            kind: BcKind::MassFlowInlet,
        };
        assert_eq!(
            format!("{context}"),
            "bc on \"inlet\" (IncompressibleFlow/MassFlowInlet)"
        );

        let condition = BoundaryCondition {
            region,
            physics: Physics::IncompressibleFlow,
            kind: BcKind::MassFlowInlet,
            value: Some(BcValue::Uniform(QtyAny::new(1.0, dims::MASS_FLOW))),
            compatibility: Some(Compat::Incompressible),
            frame: 0,
        };
        let mut findings = Vec::new();
        condition.check(&mut findings);
        assert!(findings.is_empty());
    }
}
