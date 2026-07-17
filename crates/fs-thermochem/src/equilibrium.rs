//! Bounded ideal-gas standard-state reaction equilibrium.
//!
//! This module combines one exactly conserved [`StoichiometricMatrix`] column
//! with the existing NASA-9 standard-state models. It evaluates standard
//! reaction Gibbs energy and the dimensionless pressure-activity equilibrium
//! constant under one exact gas/ideal-gas/reference-pressure convention. It
//! does not define activities away from that standard-state model, evaluate a
//! reacting mixture, integrate kinetics, infer a reverse rate, or prove that a
//! named reaction is physically meaningful.

use core::fmt;

use fs_qty::{Dimensionless, Temperature};

use crate::{
    ConservationCertificate, ElementalReferenceIdV1, MolarEnergyQuantityV1,
    Nasa9EvaluationReceiptV1, Nasa9StandardStateModelV1, ReactionId, ReferenceEquationOfStateV1,
    SpeciesId, StandardStateConventionV1, StandardStatePhaseV1, StoichiometricMatrix,
    ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
};

/// Version of the fixed standard-state reaction-equilibrium operation tree.
pub const REACTION_EQUILIBRIUM_EVALUATOR_VERSION_V1: u32 = 1;
/// Maximum species rows retained or scanned by one reaction model.
pub const MAX_REACTION_SPECIES_V1: usize = 128;
/// Maximum reaction columns whose matrix identity may be bound by one model.
pub const MAX_REACTION_COLUMNS_V1: usize = 128;
/// Maximum total stoichiometric cells hashed during admission.
pub const MAX_REACTION_MATRIX_CELLS_V1: usize = MAX_REACTION_SPECIES_V1 * MAX_REACTION_COLUMNS_V1;

const MAX_EXACT_STOICHIOMETRIC_COEFFICIENT_V1: u128 = 1u128 << 53;

/// Convention field named by an exact cross-species mismatch.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReactionConventionFieldV1 {
    /// Standard-state phase.
    Phase,
    /// Reference equation of state.
    EquationOfState,
    /// Exact reference-pressure bits.
    ReferencePressure,
    /// Opaque elemental-reference identity.
    ElementalReference,
}

/// Allocation stage named by a bounded fallible reservation refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReactionEquilibriumAllocationStageV1 {
    /// Active nonzero stoichiometric terms retained at admission.
    ActiveTerms,
    /// Fully bound model terms retained at admission.
    BoundTerms,
    /// Per-species evaluation receipts retained for replay.
    EvaluationReceipts,
}

/// Fixed-order arithmetic value named by a non-finite refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReactionEquilibriumArithmeticValueV1 {
    /// One integer stoichiometric coefficient multiplied by molar Gibbs energy.
    SpeciesGibbsContribution,
    /// Canonical partial sum of standard reaction Gibbs energy.
    StandardReactionGibbsEnergy,
    /// Positive molar scale `R T`.
    GasConstantTemperatureProduct,
    /// Dimensionless logarithmic equilibrium constant `-delta_g / (R T)`.
    LogEquilibriumConstant,
}

/// Construction or evaluation refusal for standard-state reaction equilibrium.
#[derive(Debug, Clone, PartialEq)]
pub enum ReactionEquilibriumErrorV1 {
    /// The stoichiometric species axis exceeded the admission bound.
    TooManySpecies {
        /// Offered species rows.
        offered: usize,
        /// Hard limit.
        limit: usize,
    },
    /// The stoichiometric reaction axis exceeded the admission bound.
    TooManyReactions {
        /// Offered reaction columns.
        offered: usize,
        /// Hard limit.
        limit: usize,
    },
    /// The bounded matrix-cell product overflowed or exceeded its cap.
    TooManyMatrixCells {
        /// Offered cells, or `usize::MAX` after multiplication overflow.
        offered: usize,
        /// Hard limit.
        limit: usize,
    },
    /// The exact conservation certificate names another stoichiometric matrix.
    ConservationCertificateMismatch {
        /// Matrix identity required by the supplied certificate.
        certified: [u8; 32],
        /// Identity of the supplied stoichiometric matrix.
        supplied: [u8; 32],
    },
    /// The requested reaction is absent from the canonical reaction axis.
    UnknownReaction {
        /// Missing reaction identity.
        reaction: ReactionId,
    },
    /// A zero reaction column has no equilibrium semantics.
    EmptyReaction {
        /// Zero-column reaction identity.
        reaction: ReactionId,
    },
    /// A nonzero column omitted either all reactants or all products.
    MissingReactionSide {
        /// Reaction identity.
        reaction: ReactionId,
        /// `reactant` or `product`.
        side: &'static str,
    },
    /// A matrix entry was unexpectedly unavailable after shape admission.
    MissingStoichiometricEntry {
        /// Canonical species row.
        row: usize,
        /// Canonical reaction column.
        column: usize,
    },
    /// An exact `i128` coefficient cannot be converted to binary64 unchanged.
    CoefficientNotExactlyRepresentable {
        /// Affected canonical species.
        species: SpeciesId,
        /// Exact refused coefficient.
        coefficient: i128,
    },
    /// More NASA-9 models were offered than the active-species bound.
    TooManyModels {
        /// Offered model count.
        offered: usize,
        /// Hard limit.
        limit: usize,
    },
    /// Two supplied standard-state models name the same species.
    DuplicateModelSpecies {
        /// Repeated species identity.
        species: SpeciesId,
    },
    /// A nonzero stoichiometric term has no supplied standard-state model.
    MissingSpeciesModel {
        /// Unbound active species.
        species: SpeciesId,
    },
    /// A supplied standard-state model is not an active term in the reaction.
    UnexpectedSpeciesModel {
        /// Extraneous model species.
        species: SpeciesId,
    },
    /// Active species disagree on one exact standard-state convention field.
    ConventionMismatch {
        /// Affected species.
        species: SpeciesId,
        /// First mismatched convention field.
        field: ReactionConventionFieldV1,
    },
    /// A bounded vector reservation failed before publication.
    AllocationRefused {
        /// Failed allocation stage.
        stage: ReactionEquilibriumAllocationStageV1,
        /// Requested element capacity.
        requested: usize,
    },
    /// One nested NASA-9 standard-state evaluation refused.
    SpeciesEvaluation {
        /// Affected species.
        species: SpeciesId,
        /// Exact nested refusal.
        source: ThermochemErrorV1,
    },
    /// Fixed-order arithmetic produced a non-finite value.
    UnrepresentableArithmeticValue {
        /// Active species for a term-local failure.
        species: Option<SpeciesId>,
        /// First failed operation-tree value.
        value: ReactionEquilibriumArithmeticValueV1,
        /// Exact rejected IEEE-754 bits.
        bits: u64,
    },
}

impl fmt::Display for ReactionEquilibriumErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TooManySpecies { offered, limit } => write!(
                formatter,
                "reaction equilibrium offered {offered} species, exceeding limit {limit}"
            ),
            Self::TooManyReactions { offered, limit } => write!(
                formatter,
                "reaction equilibrium offered {offered} reactions, exceeding limit {limit}"
            ),
            Self::TooManyMatrixCells { offered, limit } => write!(
                formatter,
                "reaction equilibrium matrix has {offered} cells, exceeding limit {limit}"
            ),
            Self::ConservationCertificateMismatch { .. } => formatter.write_str(
                "reaction equilibrium conservation certificate names another stoichiometric matrix",
            ),
            Self::UnknownReaction { reaction } => {
                write!(formatter, "reaction equilibrium has no reaction {reaction}")
            }
            Self::EmptyReaction { reaction } => {
                write!(formatter, "reaction equilibrium column {reaction} is zero")
            }
            Self::MissingReactionSide { reaction, side } => write!(
                formatter,
                "reaction equilibrium {reaction} has no {side} side"
            ),
            Self::MissingStoichiometricEntry { row, column } => write!(
                formatter,
                "reaction equilibrium matrix entry ({row}, {column}) is unavailable"
            ),
            Self::CoefficientNotExactlyRepresentable {
                species,
                coefficient,
            } => write!(
                formatter,
                "reaction coefficient {coefficient} for {species} is not exactly representable as binary64"
            ),
            Self::TooManyModels { offered, limit } => write!(
                formatter,
                "reaction equilibrium offered {offered} species models, exceeding limit {limit}"
            ),
            Self::DuplicateModelSpecies { species } => {
                write!(
                    formatter,
                    "reaction equilibrium repeats model species {species}"
                )
            }
            Self::MissingSpeciesModel { species } => {
                write!(
                    formatter,
                    "reaction equilibrium lacks a model for {species}"
                )
            }
            Self::UnexpectedSpeciesModel { species } => write!(
                formatter,
                "reaction equilibrium received inactive species model {species}"
            ),
            Self::ConventionMismatch { species, field } => write!(
                formatter,
                "reaction equilibrium species {species} mismatches {field:?}"
            ),
            Self::AllocationRefused { stage, requested } => write!(
                formatter,
                "reaction equilibrium {stage:?} allocation refused for {requested} elements"
            ),
            Self::SpeciesEvaluation { species, source } => write!(
                formatter,
                "reaction equilibrium standard-state evaluation for {species} refused: {source}"
            ),
            Self::UnrepresentableArithmeticValue {
                species,
                value,
                bits,
            } => write!(
                formatter,
                "reaction equilibrium {value:?} for {species:?} is non-finite (bits {bits:#018x})"
            ),
        }
    }
}

impl core::error::Error for ReactionEquilibriumErrorV1 {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::SpeciesEvaluation { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// One canonical nonzero stoichiometric term bound to a NASA-9 model.
#[derive(Debug, Clone, PartialEq)]
pub struct IdealGasReactionTermV1 {
    coefficient: i128,
    model: Nasa9StandardStateModelV1,
}

impl IdealGasReactionTermV1 {
    /// Exact negative-reactant/positive-product coefficient.
    #[must_use]
    pub const fn coefficient(&self) -> i128 {
        self.coefficient
    }

    /// Exact active species identity.
    #[must_use]
    pub const fn species(&self) -> &SpeciesId {
        self.model.species()
    }

    /// Bound standard-state model and its source cards.
    #[must_use]
    pub const fn model(&self) -> &Nasa9StandardStateModelV1 {
        &self.model
    }
}

/// One exactly conserved reaction column bound to standard-state models.
#[derive(Debug, Clone, PartialEq)]
pub struct IdealGasReactionEquilibriumV1 {
    stoichiometric: StoichiometricMatrix,
    conservation: ConservationCertificate,
    reaction: ReactionId,
    reaction_index: usize,
    convention: StandardStateConventionV1,
    terms: Vec<IdealGasReactionTermV1>,
    stoichiometric_identity: [u8; 32],
    conservation_identity: [u8; 32],
}

impl IdealGasReactionEquilibriumV1 {
    /// Bind one exactly conserved reaction to all and only its active species
    /// standard-state models.
    ///
    /// Model input order is irrelevant. Admission canonicalizes it to the
    /// `fs-qty` species axis, refuses any absent or inactive model, requires a
    /// reactant and product side, and requires every active model to use the
    /// exact same standard-state convention. The conservation certificate is
    /// consumed as existing exact `A N = 0` and `z^T N = 0` authority; this
    /// crate does not recompute or replace that authority.
    ///
    /// # Errors
    /// Refuses resource bounds, certificate/matrix drift, unknown or degenerate
    /// reactions, inexact coefficient conversion, model coverage drift, or
    /// standard-state convention drift.
    pub fn new(
        stoichiometric: StoichiometricMatrix,
        conservation: ConservationCertificate,
        reaction: ReactionId,
        mut models: Vec<Nasa9StandardStateModelV1>,
    ) -> Result<Self, ReactionEquilibriumErrorV1> {
        let species_count = stoichiometric.row_count();
        if species_count > MAX_REACTION_SPECIES_V1 {
            return Err(ReactionEquilibriumErrorV1::TooManySpecies {
                offered: species_count,
                limit: MAX_REACTION_SPECIES_V1,
            });
        }
        let reaction_count = stoichiometric.column_count();
        if reaction_count > MAX_REACTION_COLUMNS_V1 {
            return Err(ReactionEquilibriumErrorV1::TooManyReactions {
                offered: reaction_count,
                limit: MAX_REACTION_COLUMNS_V1,
            });
        }
        let matrix_cells = match species_count.checked_mul(reaction_count) {
            Some(cells) => cells,
            None => usize::MAX,
        };
        if matrix_cells > MAX_REACTION_MATRIX_CELLS_V1 {
            return Err(ReactionEquilibriumErrorV1::TooManyMatrixCells {
                offered: matrix_cells,
                limit: MAX_REACTION_MATRIX_CELLS_V1,
            });
        }
        if models.len() > MAX_REACTION_SPECIES_V1 {
            return Err(ReactionEquilibriumErrorV1::TooManyModels {
                offered: models.len(),
                limit: MAX_REACTION_SPECIES_V1,
            });
        }

        let stoichiometric_hash = stoichiometric.identity();
        let stoichiometric_identity = *stoichiometric_hash.as_bytes();
        if conservation.stoichiometric_matrix() != stoichiometric_hash {
            return Err(
                ReactionEquilibriumErrorV1::ConservationCertificateMismatch {
                    certified: *conservation.stoichiometric_matrix().as_bytes(),
                    supplied: stoichiometric_identity,
                },
            );
        }
        let reaction_index = stoichiometric
            .reactions()
            .binary_search(&reaction)
            .map_err(|_| ReactionEquilibriumErrorV1::UnknownReaction {
                reaction: reaction.clone(),
            })?;

        let mut active = Vec::new();
        active.try_reserve_exact(species_count).map_err(|_| {
            ReactionEquilibriumErrorV1::AllocationRefused {
                stage: ReactionEquilibriumAllocationStageV1::ActiveTerms,
                requested: species_count,
            }
        })?;
        let mut has_reactant = false;
        let mut has_product = false;
        for (row, species) in stoichiometric.species().iter().enumerate() {
            let coefficient = stoichiometric.get(row, reaction_index).ok_or(
                ReactionEquilibriumErrorV1::MissingStoichiometricEntry {
                    row,
                    column: reaction_index,
                },
            )?;
            if coefficient == 0 {
                continue;
            }
            if coefficient.unsigned_abs() > MAX_EXACT_STOICHIOMETRIC_COEFFICIENT_V1 {
                return Err(
                    ReactionEquilibriumErrorV1::CoefficientNotExactlyRepresentable {
                        species: species.clone(),
                        coefficient,
                    },
                );
            }
            has_reactant |= coefficient < 0;
            has_product |= coefficient > 0;
            active.push((species.clone(), coefficient));
        }
        if active.is_empty() {
            return Err(ReactionEquilibriumErrorV1::EmptyReaction { reaction });
        }
        if !has_reactant {
            return Err(ReactionEquilibriumErrorV1::MissingReactionSide {
                reaction,
                side: "reactant",
            });
        }
        if !has_product {
            return Err(ReactionEquilibriumErrorV1::MissingReactionSide {
                reaction,
                side: "product",
            });
        }

        models.sort_by(|left, right| left.species().cmp(right.species()));
        for pair in models.windows(2) {
            let [first, second] = pair else {
                continue;
            };
            if first.species() == second.species() {
                return Err(ReactionEquilibriumErrorV1::DuplicateModelSpecies {
                    species: first.species().clone(),
                });
            }
        }
        let mut active_iter = active.iter();
        let mut model_iter = models.iter();
        loop {
            match (active_iter.next(), model_iter.next()) {
                (Some((active_species, _)), Some(model)) => {
                    match active_species.cmp(model.species()) {
                        core::cmp::Ordering::Less => {
                            return Err(ReactionEquilibriumErrorV1::MissingSpeciesModel {
                                species: active_species.clone(),
                            });
                        }
                        core::cmp::Ordering::Greater => {
                            return Err(ReactionEquilibriumErrorV1::UnexpectedSpeciesModel {
                                species: model.species().clone(),
                            });
                        }
                        core::cmp::Ordering::Equal => {}
                    }
                }
                (Some((active_species, _)), None) => {
                    return Err(ReactionEquilibriumErrorV1::MissingSpeciesModel {
                        species: active_species.clone(),
                    });
                }
                (None, Some(model)) => {
                    return Err(ReactionEquilibriumErrorV1::UnexpectedSpeciesModel {
                        species: model.species().clone(),
                    });
                }
                (None, None) => break,
            }
        }

        let first_active_species = active
            .first()
            .map(|(species, _)| species.clone())
            .ok_or_else(|| ReactionEquilibriumErrorV1::EmptyReaction {
                reaction: reaction.clone(),
            })?;
        let convention = models
            .first()
            .map(|model| model.convention().clone())
            .ok_or_else(|| ReactionEquilibriumErrorV1::MissingSpeciesModel {
                species: first_active_species,
            })?;
        let mut terms = Vec::new();
        terms.try_reserve_exact(active.len()).map_err(|_| {
            ReactionEquilibriumErrorV1::AllocationRefused {
                stage: ReactionEquilibriumAllocationStageV1::BoundTerms,
                requested: active.len(),
            }
        })?;
        for ((species, coefficient), model) in active.into_iter().zip(models) {
            let field = first_convention_mismatch(&convention, model.convention());
            if let Some(field) = field {
                return Err(ReactionEquilibriumErrorV1::ConventionMismatch { species, field });
            }
            terms.push(IdealGasReactionTermV1 { coefficient, model });
        }

        Ok(Self {
            stoichiometric,
            conservation,
            reaction,
            reaction_index,
            convention,
            terms,
            stoichiometric_identity,
            conservation_identity: *conservation.identity().as_bytes(),
        })
    }

    /// Exact reaction identity selected from the matrix.
    #[must_use]
    pub const fn reaction(&self) -> &ReactionId {
        &self.reaction
    }

    /// Canonical reaction-column index.
    #[must_use]
    pub const fn reaction_index(&self) -> usize {
        self.reaction_index
    }

    /// Complete immutable stoichiometric matrix bound by the certificate.
    #[must_use]
    pub const fn stoichiometric_matrix(&self) -> &StoichiometricMatrix {
        &self.stoichiometric
    }

    /// Existing exact conservation certificate consumed by this model.
    #[must_use]
    pub const fn conservation_certificate(&self) -> ConservationCertificate {
        self.conservation
    }

    /// Common exact standard-state convention.
    #[must_use]
    pub const fn convention(&self) -> &StandardStateConventionV1 {
        &self.convention
    }

    /// Canonically ordered nonzero species terms.
    #[must_use]
    pub fn terms(&self) -> &[IdealGasReactionTermV1] {
        &self.terms
    }

    /// Evaluate standard reaction Gibbs energy and its dimensionless
    /// pressure-activity equilibrium constant.
    ///
    /// Version 1 evaluates every NASA-9 model once in canonical species order,
    /// accumulates `delta_g = sum_i(nu_i g_i^0)`, then computes
    /// `ln(K_p) = -delta_g / (R T)`. `K_p` uses dimensionless ideal-gas
    /// activities `p_i / p0`; the common `p0` is retained in every receipt.
    /// If deterministic exponentiation overflows or underflows, the finite log
    /// value remains successful and the direct constant is explicitly absent.
    ///
    /// # Errors
    /// Refuses the first nested species evaluation or non-finite fixed-order
    /// intermediate. No partial evaluation escapes.
    pub fn evaluate(
        &self,
        temperature: Temperature,
    ) -> Result<IdealGasReactionEquilibriumEvaluationV1, ReactionEquilibriumErrorV1> {
        let mut term_receipts = Vec::new();
        term_receipts
            .try_reserve_exact(self.terms.len())
            .map_err(|_| ReactionEquilibriumErrorV1::AllocationRefused {
                stage: ReactionEquilibriumAllocationStageV1::EvaluationReceipts,
                requested: self.terms.len(),
            })?;
        let mut delta_g = 0.0f64;
        for term in &self.terms {
            let evaluation = term.model.evaluate(temperature).map_err(|source| {
                ReactionEquilibriumErrorV1::SpeciesEvaluation {
                    species: term.species().clone(),
                    source,
                }
            })?;
            let gibbs = evaluation.properties().g().value();
            let contribution = (term.coefficient as f64) * gibbs;
            if !contribution.is_finite() {
                return Err(ReactionEquilibriumErrorV1::UnrepresentableArithmeticValue {
                    species: Some(term.species().clone()),
                    value: ReactionEquilibriumArithmeticValueV1::SpeciesGibbsContribution,
                    bits: contribution.to_bits(),
                });
            }
            delta_g += contribution;
            if !delta_g.is_finite() {
                return Err(ReactionEquilibriumErrorV1::UnrepresentableArithmeticValue {
                    species: Some(term.species().clone()),
                    value: ReactionEquilibriumArithmeticValueV1::StandardReactionGibbsEnergy,
                    bits: delta_g.to_bits(),
                });
            }
            term_receipts.push(ReactionEquilibriumTermReceiptV1 {
                species: term.species().clone(),
                coefficient: term.coefficient,
                standard_gibbs_bits: gibbs.to_bits(),
                contribution_bits: contribution.to_bits(),
                standard_state_receipt: evaluation.receipt().clone(),
            });
        }

        let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * temperature.value();
        if !rt.is_finite() || rt <= 0.0 {
            return Err(ReactionEquilibriumErrorV1::UnrepresentableArithmeticValue {
                species: None,
                value: ReactionEquilibriumArithmeticValueV1::GasConstantTemperatureProduct,
                bits: rt.to_bits(),
            });
        }
        let log_equilibrium_constant = -delta_g / rt;
        if !log_equilibrium_constant.is_finite() {
            return Err(ReactionEquilibriumErrorV1::UnrepresentableArithmeticValue {
                species: None,
                value: ReactionEquilibriumArithmeticValueV1::LogEquilibriumConstant,
                bits: log_equilibrium_constant.to_bits(),
            });
        }
        let raw_equilibrium_constant = fs_math::det::exp(log_equilibrium_constant);
        let (status, equilibrium_constant) = if raw_equilibrium_constant.is_infinite() {
            (EquilibriumConstantStatusV1::LogOnlyOverflow, None)
        } else if raw_equilibrium_constant == 0.0 {
            (EquilibriumConstantStatusV1::LogOnlyUnderflow, None)
        } else {
            (
                EquilibriumConstantStatusV1::FinitePositive,
                Some(ReactionEquilibriumConstantV1::new(raw_equilibrium_constant)),
            )
        };
        let receipt = IdealGasReactionEquilibriumReceiptV1 {
            evaluator_version: REACTION_EQUILIBRIUM_EVALUATOR_VERSION_V1,
            fs_math_version: fs_math::VERSION,
            fs_qty_version: fs_qty::VERSION,
            gas_constant_bits: UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K.to_bits(),
            reaction: self.reaction.clone(),
            reaction_index: self.reaction_index,
            stoichiometric_identity: self.stoichiometric_identity,
            conservation_identity: self.conservation_identity,
            temperature_bits: temperature.value().to_bits(),
            phase: self.convention.phase(),
            eos: self.convention.eos(),
            reference_pressure_bits: self.convention.reference_pressure().value().to_bits(),
            elemental_reference: self.convention.elemental_reference().clone(),
            term_receipts,
            standard_reaction_gibbs_bits: delta_g.to_bits(),
            log_equilibrium_constant_bits: log_equilibrium_constant.to_bits(),
            equilibrium_constant_status: status,
            equilibrium_constant_bits: equilibrium_constant.map(|value| value.value().to_bits()),
        };
        Ok(IdealGasReactionEquilibriumEvaluationV1 {
            temperature,
            standard_reaction_gibbs: StandardReactionGibbsEnergyV1::new(delta_g),
            log_equilibrium_constant: LogReactionEquilibriumConstantV1::new(
                log_equilibrium_constant,
            ),
            equilibrium_constant,
            status,
            receipt,
        })
    }
}

fn first_convention_mismatch(
    expected: &StandardStateConventionV1,
    found: &StandardStateConventionV1,
) -> Option<ReactionConventionFieldV1> {
    if expected.phase() != found.phase() {
        Some(ReactionConventionFieldV1::Phase)
    } else if expected.eos() != found.eos() {
        Some(ReactionConventionFieldV1::EquationOfState)
    } else if expected.reference_pressure().value().to_bits()
        != found.reference_pressure().value().to_bits()
    {
        Some(ReactionConventionFieldV1::ReferencePressure)
    } else if expected.elemental_reference() != found.elemental_reference() {
        Some(ReactionConventionFieldV1::ElementalReference)
    } else {
        None
    }
}

/// Standard reaction Gibbs energy `sum_i(nu_i g_i^0)`, J/mol-reaction.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct StandardReactionGibbsEnergyV1(MolarEnergyQuantityV1);

impl StandardReactionGibbsEnergyV1 {
    const fn new(value: f64) -> Self {
        Self(MolarEnergyQuantityV1::new(value))
    }

    /// Raw coherent-SI scalar.
    #[must_use]
    pub const fn value(self) -> f64 {
        self.0.value()
    }

    /// Dimensioned coherent-SI quantity.
    #[must_use]
    pub const fn quantity(self) -> MolarEnergyQuantityV1 {
        self.0
    }
}

/// Dimensionless `ln(K_p)` under the retained standard-state convention.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct LogReactionEquilibriumConstantV1(Dimensionless);

impl LogReactionEquilibriumConstantV1 {
    const fn new(value: f64) -> Self {
        Self(Dimensionless::new(value))
    }

    /// Raw dimensionless scalar.
    #[must_use]
    pub const fn value(self) -> f64 {
        self.0.value()
    }

    /// Dimensionless quantity wrapper.
    #[must_use]
    pub const fn quantity(self) -> Dimensionless {
        self.0
    }
}

/// Positive finite dimensionless `K_p` when representable directly.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct ReactionEquilibriumConstantV1(Dimensionless);

impl ReactionEquilibriumConstantV1 {
    const fn new(value: f64) -> Self {
        Self(Dimensionless::new(value))
    }

    /// Raw positive finite dimensionless scalar.
    #[must_use]
    pub const fn value(self) -> f64 {
        self.0.value()
    }

    /// Dimensionless quantity wrapper.
    #[must_use]
    pub const fn quantity(self) -> Dimensionless {
        self.0
    }
}

/// Direct equilibrium-constant representation status.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EquilibriumConstantStatusV1 {
    /// `K_p = exp(ln K_p)` is positive and finite.
    FinitePositive,
    /// `ln K_p` is finite but direct exponentiation overflowed.
    LogOnlyOverflow,
    /// `ln K_p` is finite but direct exponentiation underflowed to zero.
    LogOnlyUnderflow,
}

/// Exact per-species replay record for one reaction evaluation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReactionEquilibriumTermReceiptV1 {
    species: SpeciesId,
    coefficient: i128,
    standard_gibbs_bits: u64,
    contribution_bits: u64,
    standard_state_receipt: Nasa9EvaluationReceiptV1,
}

impl ReactionEquilibriumTermReceiptV1 {
    /// Canonical active species.
    #[must_use]
    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }

    /// Exact stoichiometric coefficient.
    #[must_use]
    pub const fn coefficient(&self) -> i128 {
        self.coefficient
    }

    /// Exact evaluated standard-state molar Gibbs-energy bits.
    #[must_use]
    pub const fn standard_gibbs_bits(&self) -> u64 {
        self.standard_gibbs_bits
    }

    /// Exact coefficient-times-Gibbs contribution bits.
    #[must_use]
    pub const fn contribution_bits(&self) -> u64 {
        self.contribution_bits
    }

    /// Complete nested NASA-9 source/convention receipt.
    #[must_use]
    pub const fn standard_state_receipt(&self) -> &Nasa9EvaluationReceiptV1 {
        &self.standard_state_receipt
    }
}

/// Immutable exact-field receipt for standard-state reaction equilibrium.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IdealGasReactionEquilibriumReceiptV1 {
    evaluator_version: u32,
    fs_math_version: &'static str,
    fs_qty_version: &'static str,
    gas_constant_bits: u64,
    reaction: ReactionId,
    reaction_index: usize,
    stoichiometric_identity: [u8; 32],
    conservation_identity: [u8; 32],
    temperature_bits: u64,
    phase: StandardStatePhaseV1,
    eos: ReferenceEquationOfStateV1,
    reference_pressure_bits: u64,
    elemental_reference: ElementalReferenceIdV1,
    term_receipts: Vec<ReactionEquilibriumTermReceiptV1>,
    standard_reaction_gibbs_bits: u64,
    log_equilibrium_constant_bits: u64,
    equilibrium_constant_status: EquilibriumConstantStatusV1,
    equilibrium_constant_bits: Option<u64>,
}

impl IdealGasReactionEquilibriumReceiptV1 {
    /// Evaluator/operation-tree version.
    #[must_use]
    pub const fn evaluator_version(&self) -> u32 {
        self.evaluator_version
    }

    /// Deterministic elementary-math crate version.
    #[must_use]
    pub const fn fs_math_version(&self) -> &'static str {
        self.fs_math_version
    }

    /// Quantity/chemistry authority crate version.
    #[must_use]
    pub const fn fs_qty_version(&self) -> &'static str {
        self.fs_qty_version
    }

    /// Exact universal-gas-constant bits.
    #[must_use]
    pub const fn gas_constant_bits(&self) -> u64 {
        self.gas_constant_bits
    }

    /// Exact reaction identity.
    #[must_use]
    pub const fn reaction(&self) -> &ReactionId {
        &self.reaction
    }

    /// Canonical reaction-column index.
    #[must_use]
    pub const fn reaction_index(&self) -> usize {
        self.reaction_index
    }

    /// Exact stoichiometric-matrix identity.
    #[must_use]
    pub const fn stoichiometric_identity(&self) -> &[u8; 32] {
        &self.stoichiometric_identity
    }

    /// Exact prior conservation-certificate identity.
    #[must_use]
    pub const fn conservation_identity(&self) -> &[u8; 32] {
        &self.conservation_identity
    }

    /// Exact evaluation-temperature bits.
    #[must_use]
    pub const fn temperature_bits(&self) -> u64 {
        self.temperature_bits
    }

    /// Common standard-state phase.
    #[must_use]
    pub const fn phase(&self) -> StandardStatePhaseV1 {
        self.phase
    }

    /// Common reference EOS.
    #[must_use]
    pub const fn eos(&self) -> ReferenceEquationOfStateV1 {
        self.eos
    }

    /// Exact common reference-pressure bits.
    #[must_use]
    pub const fn reference_pressure_bits(&self) -> u64 {
        self.reference_pressure_bits
    }

    /// Common opaque elemental-reference identity.
    #[must_use]
    pub const fn elemental_reference(&self) -> &ElementalReferenceIdV1 {
        &self.elemental_reference
    }

    /// Canonical exact-field per-species receipts.
    #[must_use]
    pub fn terms(&self) -> &[ReactionEquilibriumTermReceiptV1] {
        &self.term_receipts
    }

    /// Exact standard reaction Gibbs-energy bits.
    #[must_use]
    pub const fn standard_reaction_gibbs_bits(&self) -> u64 {
        self.standard_reaction_gibbs_bits
    }

    /// Exact finite `ln(K_p)` bits.
    #[must_use]
    pub const fn log_equilibrium_constant_bits(&self) -> u64 {
        self.log_equilibrium_constant_bits
    }

    /// Direct-constant representation status.
    #[must_use]
    pub const fn equilibrium_constant_status(&self) -> EquilibriumConstantStatusV1 {
        self.equilibrium_constant_status
    }

    /// Exact positive finite direct-constant bits, when representable.
    #[must_use]
    pub const fn equilibrium_constant_bits(&self) -> Option<u64> {
        self.equilibrium_constant_bits
    }
}

/// Successful standard-state reaction-equilibrium evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct IdealGasReactionEquilibriumEvaluationV1 {
    temperature: Temperature,
    standard_reaction_gibbs: StandardReactionGibbsEnergyV1,
    log_equilibrium_constant: LogReactionEquilibriumConstantV1,
    equilibrium_constant: Option<ReactionEquilibriumConstantV1>,
    status: EquilibriumConstantStatusV1,
    receipt: IdealGasReactionEquilibriumReceiptV1,
}

impl IdealGasReactionEquilibriumEvaluationV1 {
    /// Evaluated absolute temperature.
    #[must_use]
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }

    /// Standard reaction Gibbs energy.
    #[must_use]
    pub const fn standard_reaction_gibbs(&self) -> StandardReactionGibbsEnergyV1 {
        self.standard_reaction_gibbs
    }

    /// Always-retained finite logarithmic equilibrium constant.
    #[must_use]
    pub const fn log_equilibrium_constant(&self) -> LogReactionEquilibriumConstantV1 {
        self.log_equilibrium_constant
    }

    /// Direct positive finite equilibrium constant, when representable.
    #[must_use]
    pub const fn equilibrium_constant(&self) -> Option<ReactionEquilibriumConstantV1> {
        self.equilibrium_constant
    }

    /// Whether the direct constant is finite or the log representation is
    /// authoritative because exponentiation left the representable range.
    #[must_use]
    pub const fn status(&self) -> EquilibriumConstantStatusV1 {
        self.status
    }

    /// Exact convention, source, arithmetic, and result receipt.
    #[must_use]
    pub const fn receipt(&self) -> &IdealGasReactionEquilibriumReceiptV1 {
        &self.receipt
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
    use fs_qty::chemistry::{ChargeVector, ElementId, ElementalMatrix, verify_conservation};
    use fs_qty::{Dims, MolarMass, Pressure};
    use std::collections::BTreeMap;

    fn species(value: &str) -> SpeciesId {
        SpeciesId::new(value).expect("canonical species")
    }

    fn reaction(value: &str) -> ReactionId {
        ReactionId::new(value).expect("canonical reaction")
    }

    fn convention(reference_pressure: f64) -> StandardStateConventionV1 {
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(reference_pressure),
            ElementalReferenceIdV1::new("nasa-reference-elements").expect("canonical reference"),
        )
        .expect("valid convention")
    }

    fn card(coefficients: [f64; 9], reference_pressure: f64) -> ConstitutiveModelCard {
        let mut parameters = BTreeMap::new();
        for (index, (name, dims)) in [
            ("a0", Dims([0, 0, 0, 2, 0, 0])),
            ("a1", Dims([0, 0, 0, 1, 0, 0])),
            ("a2", Dims::NONE),
            ("a3", Dims([0, 0, 0, -1, 0, 0])),
            ("a4", Dims([0, 0, 0, -2, 0, 0])),
            ("a5", Dims([0, 0, 0, -3, 0, 0])),
            ("a6", Dims([0, 0, 0, -4, 0, 0])),
            ("a7", Dims([0, 0, 0, 1, 0, 0])),
            ("a8", Dims::NONE),
        ]
        .into_iter()
        .enumerate()
        {
            parameters.insert(
                name.to_string(),
                LawParameter {
                    value: coefficients[index],
                    dims,
                },
            );
        }
        parameters.insert(
            "reference_pressure".to_string(),
            LawParameter {
                value: reference_pressure,
                dims: Pressure::DIMS,
            },
        );
        ConstitutiveModelCard {
            law: LawId(crate::NASA9_LAW_ID_V1.to_string()),
            law_version: crate::NASA9_LAW_VERSION_V1,
            parameters,
            state_schema_version: crate::NASA9_STATE_SCHEMA_VERSION_V1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", 200.0, 6_000.0),
            sources: Vec::new(),
            provenance: Provenance {
                source: "reaction-equilibrium synthetic NASA-9 fixture".to_string(),
                license: "test fixture".to_string(),
                artifact: None,
            },
        }
    }

    fn model(
        species_id: &str,
        molar_mass: f64,
        a2: f64,
        a8: f64,
        reference_pressure: f64,
    ) -> Nasa9StandardStateModelV1 {
        let mut coefficients = [0.0; 9];
        coefficients[2] = a2;
        coefficients[8] = a8;
        Nasa9StandardStateModelV1::new(
            species(species_id),
            MolarMass::new(molar_mass),
            convention(reference_pressure),
            vec![
                crate::Nasa9RegionV1::from_card(card(coefficients, reference_pressure))
                    .expect("valid region"),
            ],
        )
        .expect("valid model")
    }

    fn water_models() -> Vec<Nasa9StandardStateModelV1> {
        vec![
            model("O2", 0.031_998_8, 3.5, 1.2, 100_000.0),
            model("H2O", 0.018_015_28, 4.0, 2.0, 100_000.0),
            model("H2", 0.002_015_88, 3.0, 1.0, 100_000.0),
        ]
    }

    fn water_stoichiometry(
        reversed: bool,
    ) -> (StoichiometricMatrix, ConservationCertificate, ReactionId) {
        let species_axis = vec![species("H2"), species("O2"), species("H2O")];
        let elemental = ElementalMatrix::new(
            vec![
                ElementId::new("H").expect("H"),
                ElementId::new("O").expect("O"),
            ],
            species_axis.clone(),
            vec![vec![2, 0, 2], vec![0, 2, 1]],
        )
        .expect("elemental matrix");
        let reaction_id = if reversed {
            reaction("water-dissociation")
        } else {
            reaction("water-formation")
        };
        let sign = if reversed { -1 } else { 1 };
        let stoichiometric = StoichiometricMatrix::new(
            species_axis.clone(),
            vec![reaction_id.clone()],
            vec![vec![-2 * sign], vec![-sign], vec![2 * sign]],
        )
        .expect("stoichiometric matrix");
        let charge = ChargeVector::new(species_axis, vec![0, 0, 0]).expect("charge vector");
        let certificate =
            verify_conservation(&elemental, &stoichiometric, &charge).expect("conserved reaction");
        (stoichiometric, certificate, reaction_id)
    }

    fn water_equilibrium(reversed: bool) -> IdealGasReactionEquilibriumV1 {
        let (stoichiometric, certificate, reaction_id) = water_stoichiometry(reversed);
        IdealGasReactionEquilibriumV1::new(stoichiometric, certificate, reaction_id, water_models())
            .expect("admitted water equilibrium")
    }

    #[test]
    fn g0_conserved_reaction_uses_canonical_fixed_order_and_exact_receipts() {
        let equilibrium = water_equilibrium(false);
        assert_eq!(
            equilibrium
                .terms()
                .iter()
                .map(|term| (term.species().as_str(), term.coefficient()))
                .collect::<Vec<_>>(),
            vec![("H2", -2), ("H2O", 2), ("O2", -1)]
        );
        assert_eq!(
            equilibrium
                .conservation_certificate()
                .stoichiometric_matrix(),
            equilibrium.stoichiometric_matrix().identity()
        );

        let temperature = Temperature::new(800.0);
        let evaluation = equilibrium
            .evaluate(temperature)
            .expect("equilibrium evaluation");
        let expected_delta_g = equilibrium.terms().iter().fold(0.0, |sum, term| {
            let gibbs = term
                .model()
                .evaluate(temperature)
                .expect("species evaluation")
                .properties()
                .g()
                .value();
            sum + (term.coefficient() as f64) * gibbs
        });
        assert_eq!(
            evaluation.standard_reaction_gibbs().value().to_bits(),
            expected_delta_g.to_bits()
        );
        let expected_log =
            -expected_delta_g / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * temperature.value());
        assert_eq!(
            evaluation.log_equilibrium_constant().value().to_bits(),
            expected_log.to_bits()
        );
        assert_eq!(
            evaluation.status(),
            EquilibriumConstantStatusV1::FinitePositive
        );
        assert!(evaluation.equilibrium_constant().is_some());
        assert_eq!(evaluation.receipt().terms().len(), 3);
        assert_eq!(
            evaluation.receipt().stoichiometric_identity(),
            equilibrium.stoichiometric_matrix().identity().as_bytes()
        );
        assert_eq!(
            evaluation.receipt().conservation_identity(),
            equilibrium.conservation_certificate().identity().as_bytes()
        );
    }

    #[test]
    fn g5_model_permutation_replays_and_reaction_reversal_inverts_log_constant() {
        let (stoichiometric, certificate, reaction_id) = water_stoichiometry(false);
        let mut reversed_models = water_models();
        reversed_models.reverse();
        let first = IdealGasReactionEquilibriumV1::new(
            stoichiometric.clone(),
            certificate,
            reaction_id.clone(),
            water_models(),
        )
        .expect("first model");
        let replay = IdealGasReactionEquilibriumV1::new(
            stoichiometric,
            certificate,
            reaction_id,
            reversed_models,
        )
        .expect("permuted model");
        assert_eq!(first, replay);

        let temperature = Temperature::new(1_200.0);
        let forward = first.evaluate(temperature).expect("forward equilibrium");
        let reverse = water_equilibrium(true)
            .evaluate(temperature)
            .expect("reverse equilibrium");
        assert_eq!(
            reverse.standard_reaction_gibbs().value().to_bits(),
            (-forward.standard_reaction_gibbs().value()).to_bits()
        );
        assert_eq!(
            reverse.log_equilibrium_constant().value().to_bits(),
            (-forward.log_equilibrium_constant().value()).to_bits()
        );
        let product = forward
            .equilibrium_constant()
            .expect("finite forward K")
            .value()
            * reverse
                .equilibrium_constant()
                .expect("finite reverse K")
                .value();
        assert!((product - 1.0).abs() <= 16.0 * f64::EPSILON);
    }

    #[test]
    fn g3_certificate_model_and_convention_drift_refuse_before_evaluation() {
        let (forward_stoichiometric, _, forward_reaction) = water_stoichiometry(false);
        let (_, reverse_certificate, _) = water_stoichiometry(true);
        assert!(matches!(
            IdealGasReactionEquilibriumV1::new(
                forward_stoichiometric.clone(),
                reverse_certificate,
                forward_reaction.clone(),
                water_models(),
            ),
            Err(ReactionEquilibriumErrorV1::ConservationCertificateMismatch { .. })
        ));

        let (_, forward_certificate, _) = water_stoichiometry(false);
        let mut missing = water_models();
        missing.retain(|model| model.species().as_str() != "O2");
        assert!(matches!(
            IdealGasReactionEquilibriumV1::new(
                forward_stoichiometric.clone(),
                forward_certificate,
                forward_reaction.clone(),
                missing,
            ),
            Err(ReactionEquilibriumErrorV1::MissingSpeciesModel { ref species })
                if species.as_str() == "O2"
        ));

        let (_, forward_certificate, _) = water_stoichiometry(false);
        let mut mismatched = water_models();
        mismatched[0] = model("O2", 0.031_998_8, 3.5, 1.2, 101_325.0);
        assert!(matches!(
            IdealGasReactionEquilibriumV1::new(
                forward_stoichiometric,
                forward_certificate,
                forward_reaction,
                mismatched,
            ),
            Err(ReactionEquilibriumErrorV1::ConventionMismatch {
                ref species,
                field: ReactionConventionFieldV1::ReferencePressure,
            }) if species.as_str() == "O2"
        ));
    }

    fn isomer_stoichiometry(
        coefficient: i128,
    ) -> (StoichiometricMatrix, ConservationCertificate, ReactionId) {
        let species_axis = vec![species("A"), species("B")];
        let elemental = ElementalMatrix::new(
            vec![ElementId::new("H").expect("H")],
            species_axis.clone(),
            vec![vec![1, 1]],
        )
        .expect("elemental matrix");
        let reaction_id = reaction("isomerization");
        let stoichiometric = StoichiometricMatrix::new(
            species_axis.clone(),
            vec![reaction_id.clone()],
            vec![vec![-coefficient], vec![coefficient]],
        )
        .expect("stoichiometric matrix");
        let charge = ChargeVector::new(species_axis, vec![0, 0]).expect("charge vector");
        let certificate = verify_conservation(&elemental, &stoichiometric, &charge)
            .expect("conserved isomerization");
        (stoichiometric, certificate, reaction_id)
    }

    fn degenerate_stoichiometry(
        coefficient: i128,
    ) -> (StoichiometricMatrix, ConservationCertificate, ReactionId) {
        let species_axis = vec![species("X")];
        let elemental = ElementalMatrix::new(
            vec![ElementId::new("H").expect("H")],
            species_axis.clone(),
            vec![vec![0]],
        )
        .expect("zero-count elemental matrix");
        let reaction_id = reaction("degenerate");
        let stoichiometric = StoichiometricMatrix::new(
            species_axis.clone(),
            vec![reaction_id.clone()],
            vec![vec![coefficient]],
        )
        .expect("degenerate stoichiometric matrix");
        let charge = ChargeVector::new(species_axis, vec![0]).expect("charge vector");
        let certificate = verify_conservation(&elemental, &stoichiometric, &charge)
            .expect("bookkeeping-conserved degenerate column");
        (stoichiometric, certificate, reaction_id)
    }

    #[test]
    fn g3_unknown_zero_and_one_sided_columns_have_no_equilibrium_semantics() {
        let (stoichiometric, certificate, _) = water_stoichiometry(false);
        assert!(matches!(
            IdealGasReactionEquilibriumV1::new(
                stoichiometric,
                certificate,
                reaction("not-in-matrix"),
                water_models(),
            ),
            Err(ReactionEquilibriumErrorV1::UnknownReaction { .. })
        ));

        let (stoichiometric, certificate, reaction_id) = degenerate_stoichiometry(0);
        assert!(matches!(
            IdealGasReactionEquilibriumV1::new(
                stoichiometric,
                certificate,
                reaction_id,
                Vec::new(),
            ),
            Err(ReactionEquilibriumErrorV1::EmptyReaction { .. })
        ));

        for (coefficient, side) in [(1, "reactant"), (-1, "product")] {
            let (stoichiometric, certificate, reaction_id) = degenerate_stoichiometry(coefficient);
            assert!(matches!(
                IdealGasReactionEquilibriumV1::new(
                    stoichiometric,
                    certificate,
                    reaction_id,
                    vec![model("X", 0.01, 1.0, 0.0, 100_000.0)],
                ),
                Err(ReactionEquilibriumErrorV1::MissingReactionSide {
                    side: found,
                    ..
                }) if found == side
            ));
        }
    }

    #[test]
    fn g3_large_exact_coefficients_keep_log_authority_and_inexact_ones_refuse() {
        let exact = 1i128 << 53;
        let (stoichiometric, certificate, reaction_id) = isomer_stoichiometry(exact);
        let equilibrium = IdealGasReactionEquilibriumV1::new(
            stoichiometric,
            certificate,
            reaction_id,
            vec![
                model("A", 0.01, 1.0, 0.0, 100_000.0),
                model("B", 0.01, 2.0, 0.0, 100_000.0),
            ],
        )
        .expect("maximum exact coefficient");
        let evaluation = equilibrium
            .evaluate(Temperature::new(500.0))
            .expect("log-only equilibrium");
        assert_eq!(
            evaluation.status(),
            EquilibriumConstantStatusV1::LogOnlyOverflow
        );
        assert!(evaluation.log_equilibrium_constant().value().is_finite());
        assert!(evaluation.equilibrium_constant().is_none());
        assert_eq!(evaluation.receipt().equilibrium_constant_bits(), None);

        let (stoichiometric, certificate, reaction_id) = isomer_stoichiometry(exact);
        let underflow = IdealGasReactionEquilibriumV1::new(
            stoichiometric,
            certificate,
            reaction_id,
            vec![
                model("A", 0.01, 2.0, 0.0, 100_000.0),
                model("B", 0.01, 1.0, 0.0, 100_000.0),
            ],
        )
        .expect("maximum exact coefficient with reverse Gibbs polarity")
        .evaluate(Temperature::new(500.0))
        .expect("underflow keeps log authority");
        assert_eq!(
            underflow.status(),
            EquilibriumConstantStatusV1::LogOnlyUnderflow
        );
        assert!(underflow.log_equilibrium_constant().value().is_finite());
        assert!(underflow.equilibrium_constant().is_none());

        let inexact = exact + 1;
        let (stoichiometric, certificate, reaction_id) = isomer_stoichiometry(inexact);
        assert!(matches!(
            IdealGasReactionEquilibriumV1::new(
                stoichiometric,
                certificate,
                reaction_id,
                vec![
                    model("A", 0.01, 1.0, 0.0, 100_000.0),
                    model("B", 0.01, 2.0, 0.0, 100_000.0),
                ],
            ),
            Err(ReactionEquilibriumErrorV1::CoefficientNotExactlyRepresentable {
                coefficient,
                ..
            }) if coefficient == -inexact
        ));
    }
}
