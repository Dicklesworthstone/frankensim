//! Mineral, aqueous and gas reaction thermodynamics at a shared source state.
//!
//! `Delta G° = Delta H° - T Delta S°` and `ln K = -Delta G°/(R T)` use
//! the original enthalpy/entropy claims, including their source rounding.
//! See IUPAC <https://goldbook.iupac.org/terms/view/G02629> and
//! <https://goldbook.iupac.org/terms/view/S05915>. These are dimensionless
//! standard-state activity constants, not concentration constants or a solved
//! equilibrium composition. No pH, activity coefficients, temperature
//! extrapolation, kinetics or mineral conversion is inferred.
//!
//! Charged aqueous species use conventional ionic entropy on the explicit
//! hydrogen-ion-zero reference. Neutral species use standard molar entropy.
//! Charge conservation cancels the common ionic reference. The caller still
//! associates species with source packs; pack names do not certify formulas.

use super::super::{
    ConservedReaction, FormationEnthalpySource, ReactionBalanceError, ReactionEnthalpy,
    ReactionTerm, Result, product, representable,
};
use crate::{MolarEnergyQuantityV1, MolarEntropyV1, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};
use fs_matdb::{MaterialAnswer, PropertyKey, SelectionPolicy};
use fs_qty::{Dimensionless, Dims, Pressure, QuantitySpec, Temperature};
use std::collections::BTreeMap;

/// Entropy reference selected from the admitted species charge and phase.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EntropyConvention {
    /// Standard molar entropy of a neutral or nonaqueous species.
    Absolute,
    /// Conventional aqueous-ion entropy relative to H+(aq) = zero.
    AqueousHydrogenZero,
}

/// One original entropy sample; its uncertainty and receipt are unchanged.
#[derive(Clone, Debug)]
pub struct ReactionEntropyTerm {
    /// Exactly conserved species identity and signed coefficient.
    pub term: ReactionTerm,
    /// Original normalized source identity.
    pub pack_id: String,
    /// Admitted absolute or conventional ionic entropy reference.
    pub convention: EntropyConvention,
    /// Selected entropy evidence and verified source receipt.
    pub answer: MaterialAnswer,
}

/// Standard reaction properties at the common, source-admitted temperature/pressure.
#[derive(Clone, Debug)]
pub struct ReactionEquilibrium {
    enthalpy: ReactionEnthalpy,
    entropy: MolarEntropyV1,
    gibbs: MolarEnergyQuantityV1,
    log_k: Dimensionless,
    temperature: Temperature,
    pressure: Pressure,
    entropy_terms: Vec<ReactionEntropyTerm>,
}

/// Chemical driving force for explicitly supplied dimensionless log activities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactionDrivingForce {
    /// Natural logarithm of Q = product(a_i ^ nu_i), including every participant.
    pub log_reaction_quotient: Dimensionless,
    /// Delta G° + R T ln Q [J/mol written reaction]; negative favors forward progress.
    pub molar_gibbs: MolarEnergyQuantityV1,
}

impl ConservedReaction {
    /// Resolve formation enthalpy and entropy from the same original source queries.
    ///
    /// All active species require matching typed temperature, pressure and
    /// elemental references. Both selected properties must declare the same
    /// phase, polymorph and solution convention. A charged aqueous participant
    /// additionally requires source-declared unit-activity/infinite-dilution and
    /// H+-zero qualifiers; caller flags cannot replace absent source authority.
    /// Uncertainty is retained per input, without a propagated confidence claim.
    pub fn standard_equilibrium(
        &self,
        sources: &[FormationEnthalpySource<'_>],
    ) -> Result<ReactionEquilibrium> {
        let enthalpy = self.standard_enthalpy_change(sources)?;
        let by_species: BTreeMap<_, _> = sources.iter().map(|s| (s.species, s)).collect();
        let temperature = sources[0].at.axes()["temperature"];
        let pressure = sources[0].at.axes()["pressure"];
        let mut entropy_terms = Vec::with_capacity(self.terms.len());
        let mut total = Sum::default();
        for (heat, definition) in enthalpy.terms().iter().zip(&self.definitions) {
            let source = by_species[&heat.term.species];
            let phases: Vec<_> = [
                "source-phase-solid",
                "source-phase-liquid",
                "source-phase-gas",
                "source-phase-aqueous",
            ]
            .into_iter()
            .filter(|axis| source.at.axes().get(*axis) == Some(&1.0))
            .collect();
            if phases.len() != 1 {
                return Err(ReactionBalanceError::Invalid(
                    "entropy query needs one unambiguous phase",
                ));
            }
            let ionic = definition.charge != 0 && phases[0] == "source-phase-aqueous";
            let (key, convention) = if ionic {
                (
                    "conventional-standard-molar-entropy",
                    EntropyConvention::AqueousHydrogenZero,
                )
            } else {
                ("standard-molar-entropy", EntropyConvention::Absolute)
            };
            let answer = source.pack.claims().query_typed(
                &PropertyKey::with_quantity(
                    key,
                    QuantitySpec::dimensional(Dims([2, 1, -2, -1, 0, -1])),
                ),
                source.at,
                SelectionPolicy::SingleClaimOnly,
            )?;
            source.pack.claims().verify_receipt(&answer.receipt)?;
            source.validate_selected(&answer)?;
            same_conventions(source, &heat.answer, &answer, ionic)?;
            total.add(product(
                heat.term.coefficient as f64,
                answer.evidence.value.value,
                "reaction entropy term",
            )?)?;
            entropy_terms.push(ReactionEntropyTerm {
                term: heat.term.clone(),
                pack_id: source.pack.pack_id().to_owned(),
                convention,
                answer,
            });
        }
        let entropy = total.value()?;
        product(temperature, entropy, "temperature times reaction entropy")?;
        let gibbs = (-temperature).mul_add(entropy, enthalpy.molar().value());
        representable(gibbs, false, "standard reaction Gibbs energy")?;
        let rt = product(
            UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
            temperature,
            "R times temperature",
        )?;
        let log_k = -gibbs / rt;
        representable(log_k, gibbs != 0.0, "log equilibrium constant")?;
        Ok(ReactionEquilibrium {
            enthalpy,
            entropy: MolarEntropyV1::new(entropy),
            gibbs: MolarEnergyQuantityV1::new(gibbs),
            log_k: Dimensionless::new(log_k),
            temperature: Temperature::new(temperature),
            pressure: Pressure::new(pressure),
            entropy_terms,
        })
    }
}

impl ReactionEquilibrium {
    /// Original standard formation-enthalpy result and source receipts.
    pub const fn enthalpy(&self) -> &ReactionEnthalpy {
        &self.enthalpy
    }
    /// Reaction entropy [J/(mol written reaction K)], products minus reactants.
    pub const fn molar_entropy(&self) -> MolarEntropyV1 {
        self.entropy
    }
    /// Standard reaction Gibbs energy [J/mol written reaction].
    pub const fn molar_gibbs(&self) -> MolarEnergyQuantityV1 {
        self.gibbs
    }
    /// Natural logarithm of the dimensionless standard equilibrium constant.
    /// This remains usable when exponentiating K would overflow or underflow.
    pub const fn log_equilibrium_constant(&self) -> Dimensionless {
        self.log_k
    }
    /// Common source-admitted absolute temperature; no extrapolation is performed.
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    /// Common source-admitted pressure.
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    /// Canonical original entropy samples, including conventional negative ionic values.
    pub fn entropy_terms(&self) -> &[ReactionEntropyTerm] {
        &self.entropy_terms
    }

    /// Positive dimensionless K, refusing a range-losing exponentiation.
    pub fn equilibrium_constant(&self) -> Result<Dimensionless> {
        let value = fs_math::det::exp(self.log_k.value());
        representable(value, true, "equilibrium constant exponential")?;
        Ok(Dimensionless::new(value))
    }

    /// Evaluate Delta G° + R T sum(nu_i ln(a_i)) at explicit log activities.
    ///
    /// Supply every active participant exactly once, including zero log activity
    /// for a chosen unit-activity pure phase. Activities must use the selected
    /// sources' standard states. Finite logs represent strictly positive
    /// activities; zero-activity singular limits are not silently approximated.
    /// This computes a driving force, not equilibrium speciation or a rate.
    pub fn at_log_activities(
        &self,
        activities: &[(SpeciesId, Dimensionless)],
    ) -> Result<ReactionDrivingForce> {
        if activities.len() != self.entropy_terms.len() {
            return Err(ReactionBalanceError::Invalid(
                "supply one log activity per active species",
            ));
        }
        let mut by_species = BTreeMap::new();
        for (species, value) in activities {
            if !value.value().is_finite() {
                return Err(ReactionBalanceError::Invalid("log activity must be finite"));
            }
            if by_species.insert(species, value.value()).is_some() {
                return Err(ReactionBalanceError::Species {
                    species: species.clone(),
                    reason: "duplicate log activity",
                });
            }
        }
        let mut sum = Sum::default();
        for term in &self.entropy_terms {
            let value = by_species.get(&term.term.species).ok_or_else(|| {
                ReactionBalanceError::Species {
                    species: term.term.species.clone(),
                    reason: "missing log activity",
                }
            })?;
            sum.add(product(
                term.term.coefficient as f64,
                *value,
                "log reaction quotient term",
            )?)?;
        }
        let log_q = sum.value()?;
        let rt = product(
            UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
            self.temperature.value(),
            "R times temperature",
        )?;
        product(rt, log_q, "activity chemical potential")?;
        let value = rt.mul_add(log_q, self.gibbs.value());
        representable(value, false, "reaction Gibbs energy at activities")?;
        Ok(ReactionDrivingForce {
            log_reaction_quotient: Dimensionless::new(log_q),
            molar_gibbs: MolarEnergyQuantityV1::new(value),
        })
    }
}

fn same_conventions(
    source: &FormationEnthalpySource<'_>,
    heat: &MaterialAnswer,
    entropy: &MaterialAnswer,
    ionic: bool,
) -> Result<()> {
    let claims = source.pack.claims();
    let h = claims
        .claim(heat.receipt.selected)
        .ok_or(ReactionBalanceError::Invalid("formation claim is absent"))?;
    let s = claims
        .claim(entropy.receipt.selected)
        .ok_or(ReactionBalanceError::Invalid("entropy claim is absent"))?;
    for axis in h
        .validity
        .bounds()
        .keys()
        .chain(s.validity.bounds().keys())
        .filter(|axis| {
            axis.starts_with("source-phase-")
                || axis.starts_with("source-aqueous-")
                || axis.starts_with("source-ionic-reference-")
                || axis.as_str() == "source-pure-species"
                || axis.as_str() == "source-individual-ion"
        })
    {
        if h.validity.bound(axis) != s.validity.bound(axis)
            || h.validity.axis_quantities().get(axis) != s.validity.axis_quantities().get(axis)
        {
            return Err(ReactionBalanceError::Invalid(
                "enthalpy and entropy phase or solution conventions differ",
            ));
        }
    }
    if ionic {
        for axis in [
            "source-ionic-reference-hydrogen-zero",
            "source-aqueous-unit-activity-infinite-dilution",
            "source-individual-ion",
        ] {
            if h.validity.bound(axis) != Some((1.0, 1.0))
                || s.validity.bound(axis) != Some((1.0, 1.0))
                || source.at.axes().get(axis) != Some(&1.0)
                || h.validity.axis_quantities().get(axis)
                    != Some(&QuantitySpec::dimensional(Dims::NONE))
            {
                return Err(ReactionBalanceError::Invalid(
                    "charged aqueous entropy needs a common source-declared H+-zero reference",
                ));
            }
        }
    }
    Ok(())
}

#[derive(Default)]
struct Sum {
    total: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<()> {
        let next = self.total + value;
        representable(next, false, "reaction thermodynamic sum")?;
        self.correction += if self.total.abs() >= value.abs() {
            (self.total - next) + value
        } else {
            (value - next) + self.total
        };
        self.total = next;
        self.value().map(|_| ())
    }
    fn value(&self) -> Result<f64> {
        let value = self.total + self.correction;
        representable(value, false, "reaction thermodynamic sum")?;
        Ok(value)
    }
}
