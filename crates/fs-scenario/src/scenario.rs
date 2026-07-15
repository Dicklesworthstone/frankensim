//! The `Scenario` root value: frames + BCs + load cases + combinations +
//! ensembles + environment + contact laws, with whole-scenario validity
//! checking. Violations are STRUCTURED FIXES (code, what, fix) — the
//! agent-facing refusal format — never free-text panics. Environmental
//! fields are explicit constructor arguments: nothing is defaulted
//! silently.

use crate::bc::{BcKind, BcValue, BoundaryCondition, Compat, Physics};
use crate::ensemble::StochasticEnsemble;
use crate::frame::FrameTree;
use fs_qty::{Dims, QtyAny};
use std::collections::{BTreeMap, BTreeSet};

const ACCEL_DIMS: Dims = Dims([1, 0, -2, 0, 0, 0]);
const TEMP_DIMS: Dims = Dims([0, 0, 0, 1, 0, 0]);
const PRESSURE_DIMS: Dims = Dims([-1, 1, -2, 0, 0, 0]);
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
                    fix: "express gravity in m/s² (SI exponents [1,0,-2,0,0,0])".to_string(),
                });
            }
            if !g.value.is_finite() {
                out.push(Violation {
                    code: "env-gravity-nonfinite",
                    what: format!("gravity component {i} is non-finite ({})", g.value),
                    fix: "replace every gravity component with a finite acceleration".to_string(),
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
        if !(self.ambient_temperature.value.is_finite() && self.ambient_temperature.value >= 0.0) {
            out.push(Violation {
                code: "env-temperature-range",
                what: format!(
                    "ambient absolute temperature {} K is non-finite or below absolute zero",
                    self.ambient_temperature.value
                ),
                fix: "use a finite absolute temperature greater than or equal to 0 K".to_string(),
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
        if !(self.ambient_pressure.value.is_finite() && self.ambient_pressure.value >= 0.0) {
            out.push(Violation {
                code: "env-pressure-range",
                what: format!(
                    "ambient absolute pressure {} Pa is non-finite or negative",
                    self.ambient_pressure.value
                ),
                fix: "use a finite, nonnegative absolute pressure".to_string(),
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
        if self.name.is_empty() {
            out.push(Violation {
                code: "scenario-name-empty",
                what: "scenario identity is empty".to_string(),
                fix: "give the scenario a nonempty exact UTF-8 name".to_string(),
            });
        }
        self.environment.check(&mut out);
        self.frames.check(&mut out);
        let frame_ids: BTreeSet<u32> = self.frames.frames.iter().map(|frame| frame.id.0).collect();
        for bc in &self.base_bcs {
            bc.check(&mut out);
            Self::check_bc_frame(bc, &frame_ids, &mut out);
        }

        let mut first_case_by_name = BTreeMap::new();
        for (i, case) in self.cases.iter().enumerate() {
            if case.name.is_empty() {
                out.push(Violation {
                    code: "case-name-empty",
                    what: format!("load case row {i} has an empty identity"),
                    fix: "give every load case a nonempty exact UTF-8 name".to_string(),
                });
            }
            if let Some(&first) = first_case_by_name.get(case.name.as_str()) {
                out.push(Violation {
                    code: "case-name-duplicate",
                    what: format!(
                        "load case {:?} first appears at row {first} and repeats at row {i}",
                        case.name
                    ),
                    fix: "give every load case a unique name".to_string(),
                });
            } else {
                first_case_by_name.insert(case.name.as_str(), i);
            }
            for bc in &case.bcs {
                bc.check(&mut out);
                Self::check_bc_frame(bc, &frame_ids, &mut out);
            }
        }

        let mut first_combo_by_name = BTreeMap::new();
        for (combo_index, combo) in self.combinations.iter().enumerate() {
            if combo.name.is_empty() {
                out.push(Violation {
                    code: "combo-name-empty",
                    what: format!("combination row {combo_index} has an empty identity"),
                    fix: "give every combination a nonempty exact UTF-8 name".to_string(),
                });
            }
            if let Some(&first) = first_combo_by_name.get(combo.name.as_str()) {
                out.push(Violation {
                    code: "combo-name-duplicate",
                    what: format!(
                        "combination {:?} first appears at row {first} and repeats at row {combo_index}",
                        combo.name
                    ),
                    fix: "give every combination a unique name".to_string(),
                });
            } else {
                first_combo_by_name.insert(combo.name.as_str(), combo_index);
            }

            let mut first_term_by_case = BTreeMap::new();
            for (term_index, (case, factor)) in combo.terms.iter().enumerate() {
                if case.is_empty() {
                    out.push(Violation {
                        code: "combo-case-empty",
                        what: format!(
                            "combination {:?} term {term_index} has an empty case reference",
                            combo.name
                        ),
                        fix: "reference a nonempty defined load-case name".to_string(),
                    });
                }
                if let Some(&first) = first_term_by_case.get(case.as_str()) {
                    out.push(Violation {
                        code: "combo-term-duplicate",
                        what: format!(
                            "combination {:?} references case {case:?} at terms {first} and {term_index}",
                            combo.name
                        ),
                        fix: "combine repeated case factors into one term".to_string(),
                    });
                } else {
                    first_term_by_case.insert(case.as_str(), term_index);
                }
                if !first_case_by_name.contains_key(case.as_str()) {
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

        let mut first_ensemble_by_name = BTreeMap::new();
        for (ensemble_index, e) in self.ensembles.iter().enumerate() {
            e.check(&mut out);
            if let Some(&first) = first_ensemble_by_name.get(e.name.as_str()) {
                out.push(Violation {
                    code: "ensemble-name-duplicate",
                    what: format!(
                        "ensemble {:?} first appears at row {first} and repeats at row {ensemble_index}",
                        e.name
                    ),
                    fix: "give every ensemble a unique name".to_string(),
                });
            } else {
                first_ensemble_by_name.insert(e.name.as_str(), ensemble_index);
            }
        }

        let mut first_contact_by_pair = BTreeMap::new();
        for (contact_index, c) in self.contacts.iter().enumerate() {
            if c.region_a.is_empty() || c.region_b.is_empty() {
                out.push(Violation {
                    code: "contact-region-empty",
                    what: format!(
                        "contact row {contact_index} has empty region identity {:?}/{:?}",
                        c.region_a, c.region_b
                    ),
                    fix: "bind both sides to nonempty exact UTF-8 region names".to_string(),
                });
            }
            if c.region_a == c.region_b {
                out.push(Violation {
                    code: "contact-self-pair",
                    what: format!(
                        "contact row {contact_index} pairs region {:?} with itself",
                        c.region_a
                    ),
                    fix: "name two distinct contact regions".to_string(),
                });
            }
            let pair = if c.region_a <= c.region_b {
                (c.region_a.as_str(), c.region_b.as_str())
            } else {
                (c.region_b.as_str(), c.region_a.as_str())
            };
            if let Some(&(first, first_model)) = first_contact_by_pair.get(&pair) {
                let (code, fix) = if first_model == &c.model {
                    (
                        "contact-pair-duplicate",
                        "remove the repeated unordered contact pair",
                    )
                } else {
                    (
                        "contact-pair-conflict",
                        "choose one contact model for the unordered region pair",
                    )
                };
                out.push(Violation {
                    code,
                    what: format!(
                        "unordered contact pair {pair:?} first appears at row {first} and repeats at row {contact_index}"
                    ),
                    fix: fix.to_string(),
                });
            } else {
                first_contact_by_pair.insert(pair, (contact_index, &c.model));
            }
            if let ContactModel::Coulomb { mu_s, mu_k } = c.model
                && !(mu_s.is_finite() && mu_k.is_finite() && mu_k >= 0.0 && mu_s >= mu_k)
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

    fn check_bc_frame(bc: &BoundaryCondition, frame_ids: &BTreeSet<u32>, out: &mut Vec<Violation>) {
        if bc.frame != 0 && !frame_ids.contains(&bc.frame) {
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
    /// balance to tolerance at every deterministic signal checkpoint OR a
    /// pressure outlet must exist to absorb the imbalance.
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
            let mut validation_times = vec![0.0f64];
            for bc in &set {
                if bc.kind == BcKind::MassFlowInlet
                    && let Some(BcValue::Signal(signal)) = &bc.value
                {
                    signal.append_net_flux_validation_times(&mut validation_times);
                }
            }
            validation_times.sort_by(|a, b| a.total_cmp(b));
            validation_times.dedup_by(|a, b| *a == *b);

            let mut first_imbalance = None;
            let mut compatibility_failed = false;
            for time in validation_times {
                let mut net = 0.0f64;
                let mut gross = 0.0f64;
                let mut aggregation_finite = true;
                let mut evaluation_failed = false;
                for bc in &set {
                    match bc.mass_flow_at(time) {
                        Ok(Some(flux)) => {
                            net += flux;
                            gross += flux.abs();
                            if !flux.is_finite() || !net.is_finite() || !gross.is_finite() {
                                aggregation_finite = false;
                                break;
                            }
                        }
                        Ok(None) => {}
                        Err(error) => {
                            out.push(Violation {
                                code: "flux-evaluation",
                                what: format!(
                                    "set {label:?} at t={time:.6e} s: declared mass flow could not be evaluated: {error}"
                                ),
                                fix: "repair the total-flow value or use a uniform/time-signal kg/s declaration that is finite at every compatibility checkpoint"
                                    .to_string(),
                            });
                            evaluation_failed = true;
                            break;
                        }
                    }
                }
                if evaluation_failed {
                    compatibility_failed = true;
                    break;
                }
                if !aggregation_finite {
                    out.push(Violation {
                        code: "flux-aggregation-nonfinite",
                        what: format!(
                            "set {label:?} at t={time:.6e} s: declared mass-flow aggregation overflowed or contained a non-finite value"
                        ),
                        fix: "rescale or partition the declared mass flows so their finite net/gross balance can be certified"
                            .to_string(),
                    });
                    compatibility_failed = true;
                    break;
                }
                if !has_pressure_outlet
                    && first_imbalance.is_none()
                    && gross > 0.0
                    && net.abs() > FLUX_REL_TOL * gross
                {
                    first_imbalance = Some((time, net, gross));
                }
            }
            if compatibility_failed || has_pressure_outlet {
                // The outlet may absorb every finite per-instant imbalance,
                // but it cannot make malformed or unevaluable declarations
                // admissible.
                continue;
            }
            if let Some((time, net, gross)) = first_imbalance {
                out.push(Violation {
                    code: "flux-imbalance",
                    what: format!(
                        "set {label:?} at t={time:.6e} s: declared incompressible but net mass flow is \
                         {net:+.6e} kg/s over gross {gross:.6e} kg/s with no pressure outlet"
                    ),
                    fix: "balance the declared inlet/outlet mass flows at every instant or add a \
                          pressure outlet to absorb the imbalance"
                        .to_string(),
                });
            }
        }
    }
}
