//! The `Scenario` root value: frames + BCs + load cases + combinations +
//! ensembles + environment + contact laws, with whole-scenario validity
//! checking. Violations are STRUCTURED FIXES (code, what, fix) — the
//! agent-facing refusal format — never free-text panics. Environmental
//! fields are explicit constructor arguments: nothing is defaulted
//! silently.

use crate::bc::{BcKind, BoundaryCondition, Compat, Physics};
use crate::ensemble::StochasticEnsemble;
use crate::frame::FrameTree;
use fs_qty::{Dims, QtyAny};

const ACCEL_DIMS: Dims = Dims([1, 0, -2, 0, 0]);
const TEMP_DIMS: Dims = Dims([0, 0, 0, 1, 0]);
const PRESSURE_DIMS: Dims = Dims([-1, 1, -2, 0, 0]);
/// Net-flux tolerance relative to the gross flux magnitude.
const FLUX_REL_TOL: f64 = 1e-9;

/// One validity finding with its structured fix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Violation {
    /// Stable machine-readable code.
    pub code: &'static str,
    /// What is wrong, with context.
    pub what: String,
    /// How to fix it.
    pub fix: String,
}

/// Explicit environmental fields — required, never silently defaulted.
#[derive(Debug, Clone, PartialEq)]
pub struct Environment {
    /// Gravity vector (m/s², world frame).
    pub gravity: [QtyAny; 3],
    /// Ambient temperature (K).
    pub ambient_temperature: QtyAny,
    /// Ambient pressure (Pa).
    pub ambient_pressure: QtyAny,
}

impl Environment {
    /// Standard Earth laboratory environment (explicitly chosen, not a
    /// hidden default: the call site names it).
    #[must_use]
    pub fn earth_lab() -> Self {
        Environment {
            gravity: [
                QtyAny::new(0.0, ACCEL_DIMS),
                QtyAny::new(0.0, ACCEL_DIMS),
                QtyAny::new(-9.806_65, ACCEL_DIMS),
            ],
            ambient_temperature: QtyAny::new(293.15, TEMP_DIMS),
            ambient_pressure: QtyAny::new(101_325.0, PRESSURE_DIMS),
        }
    }

    fn check(&self, out: &mut Vec<Violation>) {
        for (i, g) in self.gravity.iter().enumerate() {
            if g.dims != ACCEL_DIMS {
                out.push(Violation {
                    code: "env-gravity-dims",
                    what: format!("gravity component {i} has dimensions {:?}", g.dims.0),
                    fix: "express gravity in m/s² (SI exponents [1,0,-2,0,0])".to_string(),
                });
            }
        }
        if self.ambient_temperature.dims != TEMP_DIMS {
            out.push(Violation {
                code: "env-temperature-dims",
                what: format!(
                    "ambient temperature has dimensions {:?}",
                    self.ambient_temperature.dims.0
                ),
                fix: "express ambient temperature in kelvin".to_string(),
            });
        }
        if self.ambient_pressure.dims != PRESSURE_DIMS {
            out.push(Violation {
                code: "env-pressure-dims",
                what: format!(
                    "ambient pressure has dimensions {:?}",
                    self.ambient_pressure.dims.0
                ),
                fix: "express ambient pressure in pascals".to_string(),
            });
        }
    }
}

/// A named set of loads applied together.
#[derive(Debug, Clone, PartialEq)]
pub struct LoadCase {
    /// Case name (referenced by combinations).
    pub name: String,
    /// The case's boundary conditions / loads.
    pub bcs: Vec<BoundaryCondition>,
}

/// A factored combination of load cases (code-style, e.g. `1.2D + 1.6L`).
#[derive(Debug, Clone, PartialEq)]
pub struct Combination {
    /// Combination name.
    pub name: String,
    /// `(case name, factor)` terms.
    pub terms: Vec<(String, f64)>,
}

/// A contact/friction pairing between two named regions.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactLaw {
    /// First region.
    pub region_a: String,
    /// Second region.
    pub region_b: String,
    /// The friction model.
    pub model: ContactModel,
}

/// Supported contact models.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ContactModel {
    /// No tangential resistance.
    Frictionless,
    /// Coulomb friction with static/kinetic coefficients.
    Coulomb {
        /// Static coefficient.
        mu_s: f64,
        /// Kinetic coefficient (≤ static).
        mu_k: f64,
    },
    /// Fully tied (no relative motion).
    Tied,
}

/// The scenario root value.
#[derive(Debug, Clone, PartialEq)]
pub struct Scenario {
    /// Scenario name (IR identity).
    pub name: String,
    /// Study seed (one of the Five Explicits).
    pub seed: u64,
    /// Reference frames (world implicit).
    pub frames: FrameTree,
    /// Always-active boundary conditions.
    pub base_bcs: Vec<BoundaryCondition>,
    /// Named load cases.
    pub cases: Vec<LoadCase>,
    /// Factored combinations over the cases.
    pub combinations: Vec<Combination>,
    /// Stochastic ensembles.
    pub ensembles: Vec<StochasticEnsemble>,
    /// Contact laws.
    pub contacts: Vec<ContactLaw>,
    /// The explicit environment.
    pub environment: Environment,
}

impl Scenario {
    /// A scenario with no loads yet. The environment is a REQUIRED
    /// argument — there is deliberately no default.
    #[must_use]
    pub fn new(name: &str, seed: u64, environment: Environment) -> Self {
        Scenario {
            name: name.to_string(),
            seed,
            frames: FrameTree::new(),
            base_bcs: Vec::new(),
            cases: Vec::new(),
            combinations: Vec::new(),
            ensembles: Vec::new(),
            contacts: Vec::new(),
            environment,
        }
    }

    /// Validate the whole scenario; empty result means admissible.
    #[must_use]
    pub fn validate(&self) -> Vec<Violation> {
        let mut out = Vec::new();
        self.environment.check(&mut out);
        self.frames.check(&mut out);
        for bc in &self.base_bcs {
            bc.check(&mut out);
            self.check_bc_frame(bc, &mut out);
        }
        for (i, case) in self.cases.iter().enumerate() {
            if self.cases[..i].iter().any(|c| c.name == case.name) {
                out.push(Violation {
                    code: "case-name-duplicate",
                    what: format!("load case {:?} defined twice", case.name),
                    fix: "give every load case a unique name".to_string(),
                });
            }
            for bc in &case.bcs {
                bc.check(&mut out);
                self.check_bc_frame(bc, &mut out);
            }
        }
        for combo in &self.combinations {
            for (case, factor) in &combo.terms {
                if !self.cases.iter().any(|c| &c.name == case) {
                    out.push(Violation {
                        code: "combo-case-missing",
                        what: format!(
                            "combination {:?} references unknown case {case:?}",
                            combo.name
                        ),
                        fix: "reference only defined load cases".to_string(),
                    });
                }
                if !factor.is_finite() {
                    out.push(Violation {
                        code: "combo-factor",
                        what: format!("combination {:?} has non-finite factor", combo.name),
                        fix: "use finite combination factors".to_string(),
                    });
                }
            }
        }
        for e in &self.ensembles {
            e.check(&mut out);
        }
        for c in &self.contacts {
            if let ContactModel::Coulomb { mu_s, mu_k } = c.model
                && !(mu_k >= 0.0 && mu_s >= mu_k && mu_s.is_finite())
            {
                out.push(Violation {
                    code: "contact-coulomb-range",
                    what: format!(
                        "contact {:?}/{:?}: mu_s={mu_s}, mu_k={mu_k}",
                        c.region_a, c.region_b
                    ),
                    fix: "require 0 <= mu_k <= mu_s < inf".to_string(),
                });
            }
        }
        self.check_net_flux(&mut out);
        out
    }

    fn check_bc_frame(&self, bc: &BoundaryCondition, out: &mut Vec<Violation>) {
        if bc.frame != 0
            && !self
                .frames
                .frames
                .iter()
                .any(|f| f.id == crate::frame::FrameId(bc.frame))
        {
            out.push(Violation {
                code: "bc-frame-missing",
                what: format!(
                    "bc on {:?} references unknown frame {}",
                    bc.region, bc.frame
                ),
                fix: "reference a defined frame id (or 0 for the world)".to_string(),
            });
        }
    }

    /// Net-flux compatibility (the admission check the bead names): for
    /// every effective BC set (base alone, and base + each case), if any
    /// inlet declares `incompressible`, the declared mass flows must
    /// balance to tolerance OR a pressure outlet must exist to absorb the
    /// imbalance.
    fn check_net_flux(&self, out: &mut Vec<Violation>) {
        let base: Vec<&BoundaryCondition> = self.base_bcs.iter().collect();
        let mut sets: Vec<(String, Vec<&BoundaryCondition>)> =
            vec![("base".to_string(), base.clone())];
        for case in &self.cases {
            let mut set = base.clone();
            set.extend(case.bcs.iter());
            sets.push((format!("base+{}", case.name), set));
        }
        for (label, set) in sets {
            let declares_incompressible = set
                .iter()
                .any(|bc| bc.compatibility == Some(Compat::Incompressible));
            if !declares_incompressible {
                continue;
            }
            let has_pressure_outlet = set.iter().any(|bc| {
                bc.physics == Physics::IncompressibleFlow && bc.kind == BcKind::PressureOutlet
            });
            if has_pressure_outlet {
                continue; // the outlet absorbs whatever the inlets push in
            }
            let mut net = 0.0f64;
            let mut gross = 0.0f64;
            for bc in &set {
                if let Some(flux) = bc.mass_flow_at(0.0) {
                    net += flux;
                    gross += flux.abs();
                }
            }
            if gross > 0.0 && net.abs() > FLUX_REL_TOL * gross {
                out.push(Violation {
                    code: "flux-imbalance",
                    what: format!(
                        "set {label:?}: declared incompressible but net mass flow is \
                         {net:+.6e} kg/s over gross {gross:.6e} kg/s with no pressure outlet"
                    ),
                    fix: "balance the declared inlet/outlet mass flows or add a \
                          pressure outlet to absorb the imbalance"
                        .to_string(),
                });
            }
        }
    }
}
