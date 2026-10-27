//! The existing FSMODPK v1 parameter block is sufficient for these laws.
//! Required dimensionless selectors are explicit: `volumetric_kind` (0 =
//! quadratic, 1 = log squared), and HGO `fibre_invariants` (0 = isochoric,
//! 1 = full). No tissue properties or fallback parameters are inferred.
//! Graph inputs are H=F-I, so the graph's zero-input admission probe is the
//! stress-free reference state. Outputs are work-conjugate Piola components.

use super::{FibreInvariants, Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};
use crate::MaterialError;
use crate::graph::{
    Differentiability, EnergyBehavior, GraphError, LawNode, LawRegistry, NodeDeclaration,
    NodeOutput, NodeRole, Port, TimeParity,
};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, QueryPoint};
use fs_qty::Dims;
use std::collections::BTreeSet;

const ONE: Dims = Dims([0, 0, 0, 0, 0, 0]);
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);

/// Immutable access to the existing hyperelastic energy/stress/tangent owner.
///
/// A spatial solver may borrow a native law or a source-bound law without
/// extracting bare parameters or introducing another constitutive update.
pub trait HyperelasticMaterial {
    /// The underlying law owns constitutive evaluation and its elastic domain.
    fn hyperelastic(&self) -> &Hyperelastic;

    /// Stored energy per original reference volume.
    ///
    /// A reference-deformation adapter must transform this together with stress
    /// and tangent. Spatial consumers must retain this material seam rather
    /// than extracting a bare law and losing its reference configuration.
    fn energy_density(&self, deformation: &[f64; 9]) -> f64 {
        self.hyperelastic().energy(deformation)
    }

    /// First Piola stress conjugate to the caller's deformation gradient.
    fn first_piola(&self, deformation: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        self.hyperelastic().piola(deformation)
    }

    /// Exact derivative of first Piola stress in the same reference frame.
    fn consistent_tangent(
        &self,
        deformation: &[f64; 9],
    ) -> Result<[[f64; 9]; 9], MaterialError> {
        self.hyperelastic().tangent(deformation)
    }
}

impl HyperelasticMaterial for Hyperelastic {
    fn hyperelastic(&self) -> &Hyperelastic {
        self
    }
}

/// A native hyperelastic law retaining the exact admitted model card.
///
/// Parameters and their physical context stay fixed over the mechanical load
/// path. This does not evolve temperature or grant additional calibration,
/// uncertainty or source-domain authority beyond the supplied model card.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedHyperelastic {
    card: ConstitutiveModelCard,
    context: QueryPoint,
    law: Hyperelastic,
}

impl ResolvedHyperelastic {
    /// Admit a model with unconstrained external source conditions.
    /// Bounded cards require [`Self::try_from_model_card_at`]; an empty context
    /// cannot satisfy a temperature or other source-domain constraint.
    pub fn try_from_model_card(card: &ConstitutiveModelCard) -> Result<Self, MaterialError> {
        Self::try_from_model_card_at(card, &QueryPoint::new())
    }

    /// Admit a model at fixed external conditions in its exact typed domain.
    ///
    /// Source bounds may constrain static temperature, pressure, composition,
    /// moisture or concentration. Legacy `T` retains its kelvin convention.
    /// Mechanical strain, stress, stretch and rate axes require an evolving
    /// coordinate owner and are refused here, even when supplied in `context`.
    /// The native law's separate principal-stretch bound still gates actual F.
    /// Extra context coordinates not constrained by this card are retained but
    /// make no additional source-coverage claim.
    pub fn try_from_model_card_at(
        card: &ConstitutiveModelCard,
        context: &QueryPoint,
    ) -> Result<Self, MaterialError> {
        let law = Hyperelastic::from_model_card(card)?;
        if !card
            .validity
            .contains_typed(context.axes(), context.axis_quantities())
        {
            return Err(refusal(
                "hyperelastic source context is missing, mistyped or outside the model domain",
            ));
        }
        for axis in card.validity.bounds().keys() {
            if !crate::state_point::is_fixed_external_coordinate(
                axis,
                context.axis_quantities().get(axis).copied(),
                context.axes().get(axis).copied(),
            ) {
                return Err(refusal(format!(
                    "hyperelastic source axis {axis} is not a supported fixed external condition; evolving mechanical coordinates require an explicit owner"
                )));
            }
        }
        Ok(Self {
            card: card.clone(),
            context: context.clone(),
            law,
        })
    }

    /// Complete immutable source, including validity, provenance and identity.
    pub const fn model_card(&self) -> &ConstitutiveModelCard {
        &self.card
    }

    /// Exact external query admitted for this fixed-condition mechanical path.
    pub const fn context(&self) -> &QueryPoint {
        &self.context
    }

    /// Existing native energy, exact Piola stress and consistent tangent.
    pub const fn law(&self) -> &Hyperelastic {
        &self.law
    }
}

impl HyperelasticMaterial for ResolvedHyperelastic {
    fn hyperelastic(&self) -> &Hyperelastic {
        &self.law
    }
}

struct Parameters<'a> {
    card: &'a ConstitutiveModelCard,
    consumed: BTreeSet<String>,
}

fn refusal(what: impl Into<String>) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}

impl Parameters<'_> {
    fn get(&mut self, name: &str, dims: Dims) -> Result<f64, MaterialError> {
        let value = self
            .card
            .parameters
            .get(name)
            .ok_or_else(|| refusal(format!("{}: missing parameter {name}", self.card.law.0)))?;
        if value.dims != dims || !value.value.is_finite() {
            return Err(refusal(format!(
                "parameter {name}: expected finite {dims:?}, got {value:?}"
            )));
        }
        self.consumed.insert(name.to_string());
        Ok(value.value)
    }

    fn count(&mut self, name: &str, cap: usize) -> Result<usize, MaterialError> {
        let value = self.get(name, ONE)?;
        if value < 0.0 || value > cap as f64 || value.fract() != 0.0 {
            return Err(refusal(format!(
                "{name}: expected integer in [0,{cap}], got {value}"
            )));
        }
        Ok(value as usize)
    }
}

impl Hyperelastic {
    /// Instantiate a strict version-1 Ogden/HGO parameter block. Counts and
    /// selectors must be dimensionless, moduli Pa, directions dimensionless.
    /// Ogden uses `mu_0, alpha_0, ...`; HGO `fibre_0_x/y/z, ...`.
    /// Unknown parameters refuse rather than silently changing card meaning.
    pub fn from_model_card(card: &ConstitutiveModelCard) -> Result<Self, MaterialError> {
        card.validate().map_err(|e| refusal(e.to_string()))?;
        if card.law_version != 1
            || card.state_schema_version != 1
            || card.initial_state != InitialStatePolicy::ZeroInternalState
        {
            return Err(refusal(
                "hyperelastic cards require law/state v1 and zero internal state",
            ));
        }
        let mut p = Parameters {
            card,
            consumed: BTreeSet::new(),
        };
        let limit = p.get("stretch_limit", ONE)?;
        let kappa = p.get("kappa", PA)?;
        let volumetric = match p.count("volumetric_kind", 1)? {
            0 => Volumetric::Quadratic { kappa },
            _ => Volumetric::LogSquared { kappa },
        };
        let model = match card.law.0.as_str() {
            "ogden-decoupled" => {
                let count = p.count("term_count", super::MAX_OGDEN_TERMS)?;
                let mut terms = Vec::with_capacity(count);
                for i in 0..count {
                    terms.push(OgdenTerm {
                        mu_pa: p.get(&format!("mu_{i}"), PA)?,
                        alpha: p.get(&format!("alpha_{i}"), ONE)?,
                    });
                }
                HyperelasticModel::Ogden { terms, volumetric }
            }
            "hgo-dispersed" => {
                let c10 = p.get("c10", PA)?;
                let k1 = p.get("k1", PA)?;
                let k2 = p.get("k2", ONE)?;
                let dispersion = p.get("dispersion", ONE)?;
                let invariants = match p.count("fibre_invariants", 1)? {
                    0 => FibreInvariants::Isochoric,
                    _ => FibreInvariants::Full,
                };
                let count = p.count("fibre_count", super::MAX_FIBRE_FAMILIES)?;
                let mut fibres = Vec::with_capacity(count);
                for i in 0..count {
                    fibres.push([
                        p.get(&format!("fibre_{i}_x"), ONE)?,
                        p.get(&format!("fibre_{i}_y"), ONE)?,
                        p.get(&format!("fibre_{i}_z"), ONE)?,
                    ]);
                }
                HyperelasticModel::HolzapfelGasserOgden {
                    c10,
                    k1,
                    k2,
                    dispersion,
                    fibres,
                    invariants,
                    volumetric,
                }
            }
            _ => {
                return Err(refusal(format!(
                    "unsupported hyperelastic law {}",
                    card.law.0
                )));
            }
        };
        if p.consumed.len() != card.parameters.len() {
            return Err(refusal("unrecognized hyperelastic card parameters"));
        }
        Self::new(model, limit)
    }
}

struct HyperelasticNode {
    declaration: NodeDeclaration,
    material: Hyperelastic,
}

impl HyperelasticNode {
    fn from_card(card: &ConstitutiveModelCard) -> Result<Box<dyn LawNode>, GraphError> {
        let material = Hyperelastic::from_model_card(card).map_err(|_| GraphError::CardMismatch {
            law: card.law.0.clone(), obligation: "inadmissible hyperelastic parameters; use Hyperelastic::from_model_card for parameter diagnostics",
        })?;
        let ports = |prefix: &str, dims| {
            (0..9)
                .map(|i| Port {
                    name: format!("{prefix}_{}{}", i / 3, i % 3),
                    dims,
                    parity: TimeParity::Even,
                })
                .collect()
        };
        Ok(Box::new(Self {
            declaration: NodeDeclaration {
                law: card.law.clone(),
                law_version: 1,
                // The graph's existing constitutive storage seam (also used
                // by ElasticEnergyNode); no internal history for this law.
                role: NodeRole::InternalMemory,
                inputs: ports("H", ONE),
                outputs: ports("P", PA),
                state_slots: Vec::new(),
                state_schema_version: 1,
                calibration: card.validity.clone(),
                differentiability: if matches!(&material.model, HyperelasticModel::Ogden { .. }) {
                    Differentiability::Smooth
                } else {
                    Differentiability::PiecewiseSmooth
                },
                energy: EnergyBehavior::FreeEnergyStorage,
                tangent_claimed: true,
            },
            material,
        }))
    }

    fn deformation(&self, state: &[f64], h: &[f64]) -> Result<[f64; 9], GraphError> {
        for (which, found, expected) in [("state", state.len(), 0), ("inputs", h.len(), 9)] {
            if found != expected {
                return Err(GraphError::ArityMismatch {
                    node: self.declaration.law.0.clone(),
                    which,
                    expected,
                    found,
                });
            }
        }
        let f = core::array::from_fn(|i| h[i] + if i % 4 == 0 { 1.0 } else { 0.0 });
        self.material
            .check_state(&f)
            .map_err(|_| GraphError::CardMismatch {
                law: self.declaration.law.0.clone(),
                obligation: "F=I+H outside positive-determinant/principal-stretch domain",
            })?;
        Ok(f)
    }
}

impl LawNode for HyperelasticNode {
    fn declaration(&self) -> &NodeDeclaration {
        &self.declaration
    }

    fn evaluate(&self, state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        let f = self.deformation(state, inputs)?;
        let p = self
            .material
            .piola(&f)
            .map_err(|_| GraphError::CardMismatch {
                law: self.declaration.law.0.clone(),
                obligation: "nonfinite hyperelastic response",
            })?;
        Ok(NodeOutput {
            outputs: p.to_vec(),
            next_state: Vec::new(),
            dissipation_rate: None,
        })
    }

    fn tangent(&self, state: &[f64], inputs: &[f64]) -> Option<Vec<f64>> {
        let f = self.deformation(state, inputs).ok()?;
        Some(
            self.material
                .tangent(&f)
                .ok()?
                .into_iter()
                .flatten()
                .collect(),
        )
    }

    fn free_energy(&self, state: &[f64], inputs: &[f64]) -> Option<f64> {
        let f = self.deformation(state, inputs).ok()?;
        let w = self.material.energy(&f);
        w.is_finite().then_some(w)
    }
}

/// Add the two shipped (law,version) factories to an existing registry.
pub fn register_hyperelastic_laws(registry: &mut LawRegistry) {
    for name in ["ogden-decoupled", "hgo-dispersed"] {
        registry.register(&LawId(name.to_string()), 1, HyperelasticNode::from_card);
    }
}
