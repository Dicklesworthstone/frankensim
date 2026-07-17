//! Bounded ideal-gas standard-state reaction equilibrium.
//!
//! This module combines one exactly conserved [`StoichiometricMatrix`] column
//! with the existing NASA-9 standard-state models. It evaluates standard
//! reaction Gibbs energy and the dimensionless pressure-activity equilibrium
//! constant under one exact gas/ideal-gas/reference-pressure convention. It
//! also supports an explicitly caller-declared stoichiometric mass-action
//! closure over positive dimensionless ideal-gas activities. It does not
//! establish empirical kinetic orders, integrate kinetics, or prove that a
//! named reaction is physically meaningful.

use core::fmt;

use fs_qty::{Dimensionless, Qty, Temperature};

use crate::{
    ConservationCertificate, ElementalReferenceIdV1, MolarEnergyQuantityV1,
    Nasa9EvaluationReceiptV1, Nasa9StandardStateModelV1, ReactionId, ReferenceEquationOfStateV1,
    SpeciesId, StandardStateConventionV1, StandardStatePhaseV1, StoichiometricMatrix,
    ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
};

/// Version of the fixed standard-state reaction-equilibrium operation tree.
pub const REACTION_EQUILIBRIUM_EVALUATOR_VERSION_V1: u32 = 1;
/// Version of the reverse-progress-rate consistency operation tree.
pub const REVERSE_RATE_CLOSURE_EVALUATOR_VERSION_V1: u32 = 1;
/// Version of the declared stoichiometric mass-action operation tree.
pub const STOICHIOMETRIC_MASS_ACTION_EVALUATOR_VERSION_V1: u32 = 1;
/// Maximum species rows retained or scanned by one reaction model.
pub const MAX_REACTION_SPECIES_V1: usize = 128;
/// Maximum reaction columns whose matrix identity may be bound by one model.
pub const MAX_REACTION_COLUMNS_V1: usize = 128;
/// Maximum total stoichiometric cells hashed during admission.
pub const MAX_REACTION_MATRIX_CELLS_V1: usize = MAX_REACTION_SPECIES_V1 * MAX_REACTION_COLUMNS_V1;

const MAX_EXACT_STOICHIOMETRIC_COEFFICIENT_V1: u128 = 1u128 << 53;

/// Coherent-SI reaction-progress-rate scale, mol m^-3 s^-1.
///
/// Forward and reverse scales share these dimensions because this module's
/// mass-action activities are dimensionless ratios to the retained standard
/// state. This alias does not define kinetic orders or evaluate a net rate.
pub type ReactionProgressRateQuantityV1 = Qty<-3, 0, -1, 0, 0, 1>;

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

/// Construction or evaluation refusal for a reverse-progress-rate closure.
#[derive(Debug, Clone, PartialEq)]
pub enum ReverseRateClosureErrorV1 {
    /// The forward progress-rate scale was zero, negative, or non-finite.
    InvalidForwardProgressRateScale {
        /// Exact rejected IEEE-754 bits.
        bits: u64,
    },
    /// The nested standard-state equilibrium evaluation refused.
    Equilibrium {
        /// Exact nested refusal.
        source: ReactionEquilibriumErrorV1,
    },
    /// Positive finite inputs produced an unexpected NaN reverse scale.
    UnrepresentableReverseProgressRateScale {
        /// Exact rejected IEEE-754 bits.
        bits: u64,
    },
}

impl fmt::Display for ReverseRateClosureErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidForwardProgressRateScale { bits } => write!(
                formatter,
                "forward progress-rate scale must be positive and finite (bits {bits:#018x})"
            ),
            Self::Equilibrium { source } => {
                write!(
                    formatter,
                    "reverse-rate equilibrium evaluation refused: {source}"
                )
            }
            Self::UnrepresentableReverseProgressRateScale { bits } => write!(
                formatter,
                "reverse progress-rate scale is unrepresentable (bits {bits:#018x})"
            ),
        }
    }
}

impl core::error::Error for ReverseRateClosureErrorV1 {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Equilibrium { source } => Some(source),
            _ => None,
        }
    }
}

/// Positive finite progress-rate scale for dimensionless mass-action
/// activities.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct ReactionProgressRateScaleV1(ReactionProgressRateQuantityV1);

impl ReactionProgressRateScaleV1 {
    /// Admit a positive finite reaction-progress-rate scale.
    ///
    /// # Errors
    /// Refuses zero, negative, NaN, or infinite values.
    pub fn new(
        quantity: ReactionProgressRateQuantityV1,
    ) -> Result<Self, ReverseRateClosureErrorV1> {
        let value = quantity.value();
        if !value.is_finite() || value <= 0.0 {
            return Err(ReverseRateClosureErrorV1::InvalidForwardProgressRateScale {
                bits: value.to_bits(),
            });
        }
        Ok(Self(quantity))
    }

    /// Raw coherent-SI scalar in mol m^-3 s^-1.
    #[must_use]
    pub const fn value(self) -> f64 {
        self.0.value()
    }

    /// Dimensioned coherent-SI quantity.
    #[must_use]
    pub const fn quantity(self) -> ReactionProgressRateQuantityV1 {
        self.0
    }
}

/// Dimensionless logarithmic ratio `ln(k_reverse / k_forward)`.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct LogReverseToForwardProgressRateRatioV1(Dimensionless);

impl LogReverseToForwardProgressRateRatioV1 {
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

/// Positive finite direct ratio `k_reverse / k_forward`, when representable.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct ReverseToForwardProgressRateRatioV1(Dimensionless);

impl ReverseToForwardProgressRateRatioV1 {
    const fn new(value: f64) -> Self {
        Self(Dimensionless::new(value))
    }

    /// Raw positive finite ratio.
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

/// Direct reverse-progress-rate representation status.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReverseProgressRateStatusV1 {
    /// The direct ratio and reverse scale are both positive and finite.
    FinitePositive,
    /// `exp(ln(k_reverse / k_forward))` overflowed.
    LogOnlyRatioOverflow,
    /// `exp(ln(k_reverse / k_forward))` underflowed to zero.
    LogOnlyRatioUnderflow,
    /// The direct ratio is finite but multiplying the forward scale overflowed.
    LogOnlyScaleOverflow,
    /// The direct ratio is finite but multiplying the forward scale underflowed.
    LogOnlyScaleUnderflow,
}

/// One forward progress-rate scale closed against a standard-state
/// equilibrium model.
///
/// The algebra assumes a mass-action formulation whose activities are the
/// same dimensionless ideal-gas `p_i / p0` ratios named by the equilibrium
/// model. Under that convention, thermodynamic consistency requires
/// `k_reverse / k_forward = 1 / K_p`. This type binds and evaluates that ratio;
/// it does not choose kinetic orders, evaluate activities, or integrate a
/// mechanism.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermodynamicReverseRateClosureV1 {
    equilibrium: IdealGasReactionEquilibriumV1,
    forward_progress_rate_scale: ReactionProgressRateScaleV1,
}

impl ThermodynamicReverseRateClosureV1 {
    /// Bind a validated forward progress-rate scale to one reaction's exact
    /// standard-state equilibrium model.
    #[must_use]
    pub const fn new(
        equilibrium: IdealGasReactionEquilibriumV1,
        forward_progress_rate_scale: ReactionProgressRateScaleV1,
    ) -> Self {
        Self {
            equilibrium,
            forward_progress_rate_scale,
        }
    }

    /// Bound standard-state equilibrium authority.
    #[must_use]
    pub const fn equilibrium(&self) -> &IdealGasReactionEquilibriumV1 {
        &self.equilibrium
    }

    /// Positive finite forward progress-rate scale.
    #[must_use]
    pub const fn forward_progress_rate_scale(&self) -> ReactionProgressRateScaleV1 {
        self.forward_progress_rate_scale
    }

    /// Evaluate the thermodynamically consistent reverse-to-forward scale
    /// relation at one temperature.
    ///
    /// The logarithmic relation `ln(k_reverse / k_forward) = -ln(K_p)` is
    /// always retained after a successful nested equilibrium evaluation.
    /// Direct ratio or scale range loss is a successful explicit log-only
    /// status, never a zero or infinite physical rate publication.
    ///
    /// # Errors
    /// Propagates the exact nested equilibrium refusal or an impossible NaN
    /// arithmetic result. No partial receipt escapes.
    pub fn evaluate(
        &self,
        temperature: Temperature,
    ) -> Result<ThermodynamicReverseRateEvaluationV1, ReverseRateClosureErrorV1> {
        let equilibrium = self
            .equilibrium
            .evaluate(temperature)
            .map_err(|source| ReverseRateClosureErrorV1::Equilibrium { source })?;
        let log_ratio = -equilibrium.log_equilibrium_constant().value();
        let raw_ratio = fs_math::det::exp(log_ratio);
        let (status, ratio, reverse_progress_rate_scale) = if raw_ratio.is_infinite() {
            (
                ReverseProgressRateStatusV1::LogOnlyRatioOverflow,
                None,
                None,
            )
        } else if raw_ratio == 0.0 {
            (
                ReverseProgressRateStatusV1::LogOnlyRatioUnderflow,
                None,
                None,
            )
        } else {
            let ratio = ReverseToForwardProgressRateRatioV1::new(raw_ratio);
            let raw_reverse = self.forward_progress_rate_scale.value() * raw_ratio;
            if raw_reverse.is_infinite() {
                (
                    ReverseProgressRateStatusV1::LogOnlyScaleOverflow,
                    Some(ratio),
                    None,
                )
            } else if raw_reverse == 0.0 {
                (
                    ReverseProgressRateStatusV1::LogOnlyScaleUnderflow,
                    Some(ratio),
                    None,
                )
            } else if raw_reverse.is_finite() {
                (
                    ReverseProgressRateStatusV1::FinitePositive,
                    Some(ratio),
                    Some(ReactionProgressRateScaleV1(
                        ReactionProgressRateQuantityV1::new(raw_reverse),
                    )),
                )
            } else {
                return Err(
                    ReverseRateClosureErrorV1::UnrepresentableReverseProgressRateScale {
                        bits: raw_reverse.to_bits(),
                    },
                );
            }
        };
        let receipt = ThermodynamicReverseRateReceiptV1 {
            evaluator_version: REVERSE_RATE_CLOSURE_EVALUATOR_VERSION_V1,
            fs_math_version: fs_math::VERSION,
            equilibrium_receipt: equilibrium.receipt().clone(),
            forward_progress_rate_scale_bits: self.forward_progress_rate_scale.value().to_bits(),
            log_reverse_to_forward_ratio_bits: log_ratio.to_bits(),
            status,
            reverse_to_forward_ratio_bits: ratio.map(|value| value.value().to_bits()),
            reverse_progress_rate_scale_bits: reverse_progress_rate_scale
                .map(|value| value.value().to_bits()),
        };
        Ok(ThermodynamicReverseRateEvaluationV1 {
            equilibrium,
            forward_progress_rate_scale: self.forward_progress_rate_scale,
            log_reverse_to_forward_ratio: LogReverseToForwardProgressRateRatioV1::new(log_ratio),
            reverse_to_forward_ratio: ratio,
            reverse_progress_rate_scale,
            status,
            receipt,
        })
    }
}

/// Exact-field replay receipt for one reverse-progress-rate closure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ThermodynamicReverseRateReceiptV1 {
    evaluator_version: u32,
    fs_math_version: &'static str,
    equilibrium_receipt: IdealGasReactionEquilibriumReceiptV1,
    forward_progress_rate_scale_bits: u64,
    log_reverse_to_forward_ratio_bits: u64,
    status: ReverseProgressRateStatusV1,
    reverse_to_forward_ratio_bits: Option<u64>,
    reverse_progress_rate_scale_bits: Option<u64>,
}

impl ThermodynamicReverseRateReceiptV1 {
    /// Reverse-rate operation-tree version.
    #[must_use]
    pub const fn evaluator_version(&self) -> u32 {
        self.evaluator_version
    }

    /// Deterministic elementary-math crate version.
    #[must_use]
    pub const fn fs_math_version(&self) -> &'static str {
        self.fs_math_version
    }

    /// Complete nested standard-state equilibrium receipt.
    #[must_use]
    pub const fn equilibrium_receipt(&self) -> &IdealGasReactionEquilibriumReceiptV1 {
        &self.equilibrium_receipt
    }

    /// Exact forward progress-rate-scale bits.
    #[must_use]
    pub const fn forward_progress_rate_scale_bits(&self) -> u64 {
        self.forward_progress_rate_scale_bits
    }

    /// Exact finite `ln(k_reverse / k_forward)` bits.
    #[must_use]
    pub const fn log_reverse_to_forward_ratio_bits(&self) -> u64 {
        self.log_reverse_to_forward_ratio_bits
    }

    /// Direct reverse representation status.
    #[must_use]
    pub const fn status(&self) -> ReverseProgressRateStatusV1 {
        self.status
    }

    /// Exact direct ratio bits, when representable.
    #[must_use]
    pub const fn reverse_to_forward_ratio_bits(&self) -> Option<u64> {
        self.reverse_to_forward_ratio_bits
    }

    /// Exact direct reverse progress-rate-scale bits, when representable.
    #[must_use]
    pub const fn reverse_progress_rate_scale_bits(&self) -> Option<u64> {
        self.reverse_progress_rate_scale_bits
    }
}

/// Successful thermodynamic reverse-progress-rate closure evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermodynamicReverseRateEvaluationV1 {
    equilibrium: IdealGasReactionEquilibriumEvaluationV1,
    forward_progress_rate_scale: ReactionProgressRateScaleV1,
    log_reverse_to_forward_ratio: LogReverseToForwardProgressRateRatioV1,
    reverse_to_forward_ratio: Option<ReverseToForwardProgressRateRatioV1>,
    reverse_progress_rate_scale: Option<ReactionProgressRateScaleV1>,
    status: ReverseProgressRateStatusV1,
    receipt: ThermodynamicReverseRateReceiptV1,
}

impl ThermodynamicReverseRateEvaluationV1 {
    /// Nested standard-state equilibrium evaluation.
    #[must_use]
    pub const fn equilibrium(&self) -> &IdealGasReactionEquilibriumEvaluationV1 {
        &self.equilibrium
    }

    /// Positive finite forward progress-rate scale.
    #[must_use]
    pub const fn forward_progress_rate_scale(&self) -> ReactionProgressRateScaleV1 {
        self.forward_progress_rate_scale
    }

    /// Always-retained finite `ln(k_reverse / k_forward)`.
    #[must_use]
    pub const fn log_reverse_to_forward_ratio(&self) -> LogReverseToForwardProgressRateRatioV1 {
        self.log_reverse_to_forward_ratio
    }

    /// Direct positive finite reverse-to-forward ratio, when representable.
    #[must_use]
    pub const fn reverse_to_forward_ratio(&self) -> Option<ReverseToForwardProgressRateRatioV1> {
        self.reverse_to_forward_ratio
    }

    /// Direct positive finite reverse progress-rate scale, when representable.
    #[must_use]
    pub const fn reverse_progress_rate_scale(&self) -> Option<ReactionProgressRateScaleV1> {
        self.reverse_progress_rate_scale
    }

    /// Whether direct ratio and reverse-scale representations survived range.
    #[must_use]
    pub const fn status(&self) -> ReverseProgressRateStatusV1 {
        self.status
    }

    /// Exact nested-authority and arithmetic receipt.
    #[must_use]
    pub const fn receipt(&self) -> &ThermodynamicReverseRateReceiptV1 {
        &self.receipt
    }
}

/// Explicit source of the forward and reverse mass-action exponents.
///
/// Version 1 supports only a caller declaration that the exact
/// stoichiometric coefficients are the kinetic orders. The declaration is
/// retained because exact reaction bookkeeping does not establish this
/// empirical kinetic-law choice.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MassActionKineticOrderConventionV1 {
    /// Negative coefficients define forward orders and positive coefficients
    /// define reverse orders.
    CallerDeclaredStoichiometricCoefficients,
}

/// Direction named by a mass-action arithmetic refusal or status.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MassActionDirectionV1 {
    /// Reactant-to-product direction.
    Forward,
    /// Product-to-reactant direction.
    Reverse,
}

/// Allocation stage named by a bounded mass-action reservation refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MassActionAllocationStageV1 {
    /// Canonical per-species term receipts.
    TermReceipts,
}

/// Fixed-order value named by a non-finite mass-action refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MassActionArithmeticValueV1 {
    /// Natural logarithm of one positive dimensionless activity.
    LogActivity,
    /// Exact stoichiometric exponent multiplied by one log activity.
    LogActivityContribution,
    /// Canonical partial sum of one directional log activity product.
    LogActivityProduct,
    /// Deterministic exponential of one log activity product.
    ActivityProduct,
    /// Directional scale multiplied by its direct activity product.
    DirectionalProgressRate,
    /// Forward progress rate minus reverse progress rate.
    NetProgressRate,
}

/// Direct representation status for one directional mass-action rate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MassActionDirectionalRateStatusV1 {
    /// Activity product and directional progress rate are positive and finite.
    FinitePositive,
    /// The log product is finite but direct exponentiation overflowed.
    LogOnlyActivityProductOverflow,
    /// The log product is finite but direct exponentiation underflowed.
    LogOnlyActivityProductUnderflow,
    /// The direct activity product exists, but the thermodynamic reverse scale
    /// is log-only and cannot be multiplied directly.
    ProgressRateScaleUnavailable,
    /// Direct activity product exists, but multiplication by the scale
    /// overflowed.
    ProgressRateOverflow,
    /// Direct activity product exists, but multiplication by the scale
    /// underflowed.
    ProgressRateUnderflow,
}

/// Construction or evaluation refusal for declared stoichiometric
/// mass-action progress.
#[derive(Debug, Clone, PartialEq)]
pub enum StoichiometricMassActionErrorV1 {
    /// One activity was zero, negative, or non-finite.
    InvalidActivity {
        /// Affected species.
        species: SpeciesId,
        /// Exact rejected IEEE-754 bits.
        bits: u64,
    },
    /// More activity entries were offered than the hard reaction bound.
    TooManyActivities {
        /// Offered activity entries.
        offered: usize,
        /// Hard limit.
        limit: usize,
    },
    /// Activity count differs from the exact active reaction-term count.
    ActivityCountMismatch {
        /// Required active terms.
        expected: usize,
        /// Offered activity entries.
        found: usize,
    },
    /// Two activity entries name the same species.
    DuplicateActivitySpecies {
        /// Repeated species identity.
        species: SpeciesId,
    },
    /// One active reaction term has no activity.
    MissingSpeciesActivity {
        /// Missing active species.
        species: SpeciesId,
    },
    /// One supplied activity is not active in the selected reaction.
    UnexpectedSpeciesActivity {
        /// Extraneous species identity.
        species: SpeciesId,
    },
    /// A bounded vector reservation failed before publication.
    AllocationRefused {
        /// Failed allocation stage.
        stage: MassActionAllocationStageV1,
        /// Requested element capacity.
        requested: usize,
    },
    /// The nested thermodynamic reverse-rate closure refused.
    ReverseRateClosure {
        /// Exact nested refusal.
        source: ReverseRateClosureErrorV1,
    },
    /// Fixed-order arithmetic produced an unexpected non-finite value.
    UnrepresentableArithmeticValue {
        /// Species for a term-local refusal.
        species: Option<SpeciesId>,
        /// Direction for a directional refusal.
        direction: Option<MassActionDirectionV1>,
        /// First failed operation-tree value.
        value: MassActionArithmeticValueV1,
        /// Exact rejected IEEE-754 bits.
        bits: u64,
    },
}

impl fmt::Display for StoichiometricMassActionErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidActivity { species, bits } => write!(
                formatter,
                "ideal-gas activity for {species} must be positive and finite (bits {bits:#018x})"
            ),
            Self::TooManyActivities { offered, limit } => write!(
                formatter,
                "mass-action evaluation offered {offered} activities, exceeding limit {limit}"
            ),
            Self::ActivityCountMismatch { expected, found } => write!(
                formatter,
                "mass-action evaluation requires {expected} activities but received {found}"
            ),
            Self::DuplicateActivitySpecies { species } => {
                write!(formatter, "mass-action activities repeat species {species}")
            }
            Self::MissingSpeciesActivity { species } => {
                write!(
                    formatter,
                    "mass-action evaluation lacks activity for {species}"
                )
            }
            Self::UnexpectedSpeciesActivity { species } => write!(
                formatter,
                "mass-action evaluation received inactive species activity {species}"
            ),
            Self::AllocationRefused { stage, requested } => write!(
                formatter,
                "mass-action {stage:?} allocation refused for {requested} elements"
            ),
            Self::ReverseRateClosure { source } => {
                write!(
                    formatter,
                    "mass-action reverse-rate closure refused: {source}"
                )
            }
            Self::UnrepresentableArithmeticValue {
                species,
                direction,
                value,
                bits,
            } => write!(
                formatter,
                "mass-action {value:?} for {species:?} in {direction:?} is non-finite (bits {bits:#018x})"
            ),
        }
    }
}

impl core::error::Error for StoichiometricMassActionErrorV1 {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::ReverseRateClosure { source } => Some(source),
            _ => None,
        }
    }
}

/// One positive finite dimensionless ideal-gas activity `p_i / p0`.
#[derive(Debug, Clone, PartialEq)]
pub struct IdealGasSpeciesActivityV1 {
    species: SpeciesId,
    activity: Dimensionless,
}

impl IdealGasSpeciesActivityV1 {
    /// Admit one positive finite dimensionless activity.
    ///
    /// # Errors
    /// Refuses zero, negative, NaN, or infinite values.
    pub fn new(
        species: SpeciesId,
        activity: Dimensionless,
    ) -> Result<Self, StoichiometricMassActionErrorV1> {
        let value = activity.value();
        if !value.is_finite() || value <= 0.0 {
            return Err(StoichiometricMassActionErrorV1::InvalidActivity {
                species,
                bits: value.to_bits(),
            });
        }
        Ok(Self { species, activity })
    }

    /// Canonical species identity.
    #[must_use]
    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }

    /// Positive finite dimensionless activity.
    #[must_use]
    pub const fn activity(&self) -> Dimensionless {
        self.activity
    }

    /// Raw dimensionless scalar.
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.activity.value()
    }
}

/// Finite logarithm of one directional dimensionless activity product.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct LogMassActionActivityProductV1(Dimensionless);

impl LogMassActionActivityProductV1 {
    const fn new(value: f64) -> Self {
        Self(Dimensionless::new(value))
    }

    /// Raw dimensionless logarithm.
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

/// Positive finite direct directional activity product, when representable.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct MassActionActivityProductV1(Dimensionless);

impl MassActionActivityProductV1 {
    const fn new(value: f64) -> Self {
        Self(Dimensionless::new(value))
    }

    /// Raw positive finite dimensionless product.
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

/// Positive finite forward or reverse reaction-progress rate.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct DirectionalReactionProgressRateV1(ReactionProgressRateQuantityV1);

impl DirectionalReactionProgressRateV1 {
    const fn new(value: f64) -> Self {
        Self(ReactionProgressRateQuantityV1::new(value))
    }

    /// Raw coherent-SI scalar in mol m^-3 s^-1.
    #[must_use]
    pub const fn value(self) -> f64 {
        self.0.value()
    }

    /// Dimensioned coherent-SI quantity.
    #[must_use]
    pub const fn quantity(self) -> ReactionProgressRateQuantityV1 {
        self.0
    }
}

/// Finite signed net reaction-progress rate, forward minus reverse.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct NetReactionProgressRateV1(ReactionProgressRateQuantityV1);

impl NetReactionProgressRateV1 {
    const fn new(value: f64) -> Self {
        Self(ReactionProgressRateQuantityV1::new(value))
    }

    /// Raw coherent-SI scalar in mol m^-3 s^-1.
    #[must_use]
    pub const fn value(self) -> f64 {
        self.0.value()
    }

    /// Dimensioned coherent-SI quantity.
    #[must_use]
    pub const fn quantity(self) -> ReactionProgressRateQuantityV1 {
        self.0
    }
}

/// One caller-declared stoichiometric mass-action law closed against exact
/// standard-state equilibrium authority.
///
/// This type deliberately requires an explicit kinetic-order convention. It
/// evaluates only one already-admitted reaction and does not infer that the
/// reaction is elementary, validate a forward scale against source data, or
/// own mechanism evolution.
#[derive(Debug, Clone, PartialEq)]
pub struct DeclaredStoichiometricMassActionV1 {
    reverse_rate_closure: ThermodynamicReverseRateClosureV1,
    kinetic_order_convention: MassActionKineticOrderConventionV1,
}

impl DeclaredStoichiometricMassActionV1 {
    /// Bind one thermodynamically closed progress-rate scale to an explicit
    /// kinetic-order declaration.
    #[must_use]
    pub const fn new(
        reverse_rate_closure: ThermodynamicReverseRateClosureV1,
        kinetic_order_convention: MassActionKineticOrderConventionV1,
    ) -> Self {
        Self {
            reverse_rate_closure,
            kinetic_order_convention,
        }
    }

    /// Bound thermodynamic reverse-rate closure.
    #[must_use]
    pub const fn reverse_rate_closure(&self) -> &ThermodynamicReverseRateClosureV1 {
        &self.reverse_rate_closure
    }

    /// Explicit source of the kinetic exponents.
    #[must_use]
    pub const fn kinetic_order_convention(&self) -> MassActionKineticOrderConventionV1 {
        self.kinetic_order_convention
    }

    /// Evaluate forward, reverse, and net progress for one activity state.
    ///
    /// Activities are canonicalized by `SpeciesId` and must cover all and only
    /// the selected reaction's active terms. Version 1 evaluates directional
    /// products in log space using the exact coefficient magnitudes, retains
    /// finite log products across direct range loss, and publishes a signed
    /// net rate only when both directional rates are directly representable.
    ///
    /// # Errors
    /// Refuses activity count/identity drift, bounded receipt allocation,
    /// nested thermodynamic closure failure, or an impossible non-finite
    /// fixed-order intermediate. Expected exponential and multiplication
    /// range loss is a successful explicit status.
    pub fn evaluate(
        &self,
        temperature: Temperature,
        mut activities: Vec<IdealGasSpeciesActivityV1>,
    ) -> Result<DeclaredStoichiometricMassActionEvaluationV1, StoichiometricMassActionErrorV1> {
        if activities.len() > MAX_REACTION_SPECIES_V1 {
            return Err(StoichiometricMassActionErrorV1::TooManyActivities {
                offered: activities.len(),
                limit: MAX_REACTION_SPECIES_V1,
            });
        }
        let terms = self.reverse_rate_closure.equilibrium().terms();
        if activities.len() != terms.len() {
            return Err(StoichiometricMassActionErrorV1::ActivityCountMismatch {
                expected: terms.len(),
                found: activities.len(),
            });
        }

        activities.sort_unstable_by(|left, right| left.species().cmp(right.species()));
        for pair in activities.windows(2) {
            let [first, second] = pair else {
                continue;
            };
            if first.species() == second.species() {
                return Err(StoichiometricMassActionErrorV1::DuplicateActivitySpecies {
                    species: first.species().clone(),
                });
            }
        }
        for (term, activity) in terms.iter().zip(&activities) {
            match term.species().cmp(activity.species()) {
                core::cmp::Ordering::Less => {
                    return Err(StoichiometricMassActionErrorV1::MissingSpeciesActivity {
                        species: term.species().clone(),
                    });
                }
                core::cmp::Ordering::Greater => {
                    return Err(StoichiometricMassActionErrorV1::UnexpectedSpeciesActivity {
                        species: activity.species().clone(),
                    });
                }
                core::cmp::Ordering::Equal => {}
            }
        }

        let mut term_receipts = Vec::new();
        term_receipts.try_reserve_exact(terms.len()).map_err(|_| {
            StoichiometricMassActionErrorV1::AllocationRefused {
                stage: MassActionAllocationStageV1::TermReceipts,
                requested: terms.len(),
            }
        })?;
        let mut log_forward_activity_product = 0.0f64;
        let mut log_reverse_activity_product = 0.0f64;
        for (term, activity) in terms.iter().zip(&activities) {
            let log_activity = fs_math::det::ln(activity.value());
            if !log_activity.is_finite() {
                return Err(
                    StoichiometricMassActionErrorV1::UnrepresentableArithmeticValue {
                        species: Some(activity.species().clone()),
                        direction: None,
                        value: MassActionArithmeticValueV1::LogActivity,
                        bits: log_activity.to_bits(),
                    },
                );
            }
            let exponent = term.coefficient().unsigned_abs() as f64;
            let contribution = exponent * log_activity;
            let direction = if term.coefficient() < 0 {
                MassActionDirectionV1::Forward
            } else {
                MassActionDirectionV1::Reverse
            };
            if !contribution.is_finite() {
                return Err(
                    StoichiometricMassActionErrorV1::UnrepresentableArithmeticValue {
                        species: Some(activity.species().clone()),
                        direction: Some(direction),
                        value: MassActionArithmeticValueV1::LogActivityContribution,
                        bits: contribution.to_bits(),
                    },
                );
            }
            let directional_sum = if direction == MassActionDirectionV1::Forward {
                &mut log_forward_activity_product
            } else {
                &mut log_reverse_activity_product
            };
            *directional_sum += contribution;
            if !directional_sum.is_finite() {
                return Err(
                    StoichiometricMassActionErrorV1::UnrepresentableArithmeticValue {
                        species: Some(activity.species().clone()),
                        direction: Some(direction),
                        value: MassActionArithmeticValueV1::LogActivityProduct,
                        bits: directional_sum.to_bits(),
                    },
                );
            }
            term_receipts.push(StoichiometricMassActionTermReceiptV1 {
                species: activity.species().clone(),
                coefficient: term.coefficient(),
                activity_bits: activity.value().to_bits(),
                log_activity_bits: log_activity.to_bits(),
                direction,
                log_contribution_bits: contribution.to_bits(),
            });
        }

        let reverse_rate = self
            .reverse_rate_closure
            .evaluate(temperature)
            .map_err(|source| StoichiometricMassActionErrorV1::ReverseRateClosure { source })?;
        let (forward_status, forward_activity_product, forward_progress_rate) =
            evaluate_mass_action_direction(
                log_forward_activity_product,
                Some(reverse_rate.forward_progress_rate_scale()),
                MassActionDirectionV1::Forward,
            )?;
        let (reverse_status, reverse_activity_product, reverse_progress_rate) =
            evaluate_mass_action_direction(
                log_reverse_activity_product,
                reverse_rate.reverse_progress_rate_scale(),
                MassActionDirectionV1::Reverse,
            )?;
        let net_progress_rate = match (forward_progress_rate, reverse_progress_rate) {
            (Some(forward), Some(reverse)) => {
                let value = forward.value() - reverse.value();
                if !value.is_finite() {
                    return Err(
                        StoichiometricMassActionErrorV1::UnrepresentableArithmeticValue {
                            species: None,
                            direction: None,
                            value: MassActionArithmeticValueV1::NetProgressRate,
                            bits: value.to_bits(),
                        },
                    );
                }
                Some(NetReactionProgressRateV1::new(value))
            }
            _ => None,
        };
        let receipt = DeclaredStoichiometricMassActionReceiptV1 {
            evaluator_version: STOICHIOMETRIC_MASS_ACTION_EVALUATOR_VERSION_V1,
            fs_math_version: fs_math::VERSION,
            fs_qty_version: fs_qty::VERSION,
            kinetic_order_convention: self.kinetic_order_convention,
            reverse_rate_receipt: reverse_rate.receipt().clone(),
            term_receipts,
            log_forward_activity_product_bits: log_forward_activity_product.to_bits(),
            log_reverse_activity_product_bits: log_reverse_activity_product.to_bits(),
            forward_status,
            reverse_status,
            forward_activity_product_bits: forward_activity_product
                .map(|value| value.value().to_bits()),
            reverse_activity_product_bits: reverse_activity_product
                .map(|value| value.value().to_bits()),
            forward_progress_rate_bits: forward_progress_rate.map(|value| value.value().to_bits()),
            reverse_progress_rate_bits: reverse_progress_rate.map(|value| value.value().to_bits()),
            net_progress_rate_bits: net_progress_rate.map(|value| value.value().to_bits()),
        };
        Ok(DeclaredStoichiometricMassActionEvaluationV1 {
            reverse_rate,
            log_forward_activity_product: LogMassActionActivityProductV1::new(
                log_forward_activity_product,
            ),
            log_reverse_activity_product: LogMassActionActivityProductV1::new(
                log_reverse_activity_product,
            ),
            forward_activity_product,
            reverse_activity_product,
            forward_progress_rate,
            reverse_progress_rate,
            net_progress_rate,
            forward_status,
            reverse_status,
            receipt,
        })
    }
}

fn evaluate_mass_action_direction(
    log_activity_product: f64,
    progress_rate_scale: Option<ReactionProgressRateScaleV1>,
    direction: MassActionDirectionV1,
) -> Result<
    (
        MassActionDirectionalRateStatusV1,
        Option<MassActionActivityProductV1>,
        Option<DirectionalReactionProgressRateV1>,
    ),
    StoichiometricMassActionErrorV1,
> {
    let raw_product = fs_math::det::exp(log_activity_product);
    if raw_product.is_infinite() {
        return Ok((
            MassActionDirectionalRateStatusV1::LogOnlyActivityProductOverflow,
            None,
            None,
        ));
    }
    if raw_product == 0.0 {
        return Ok((
            MassActionDirectionalRateStatusV1::LogOnlyActivityProductUnderflow,
            None,
            None,
        ));
    }
    if !raw_product.is_finite() || raw_product < 0.0 {
        return Err(
            StoichiometricMassActionErrorV1::UnrepresentableArithmeticValue {
                species: None,
                direction: Some(direction),
                value: MassActionArithmeticValueV1::ActivityProduct,
                bits: raw_product.to_bits(),
            },
        );
    }
    let product = MassActionActivityProductV1::new(raw_product);
    let Some(scale) = progress_rate_scale else {
        return Ok((
            MassActionDirectionalRateStatusV1::ProgressRateScaleUnavailable,
            Some(product),
            None,
        ));
    };
    let raw_rate = scale.value() * raw_product;
    if raw_rate.is_infinite() {
        return Ok((
            MassActionDirectionalRateStatusV1::ProgressRateOverflow,
            Some(product),
            None,
        ));
    }
    if raw_rate == 0.0 {
        return Ok((
            MassActionDirectionalRateStatusV1::ProgressRateUnderflow,
            Some(product),
            None,
        ));
    }
    if !raw_rate.is_finite() || raw_rate < 0.0 {
        return Err(
            StoichiometricMassActionErrorV1::UnrepresentableArithmeticValue {
                species: None,
                direction: Some(direction),
                value: MassActionArithmeticValueV1::DirectionalProgressRate,
                bits: raw_rate.to_bits(),
            },
        );
    }
    Ok((
        MassActionDirectionalRateStatusV1::FinitePositive,
        Some(product),
        Some(DirectionalReactionProgressRateV1::new(raw_rate)),
    ))
}

/// Canonical exact-field receipt for one active mass-action term.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StoichiometricMassActionTermReceiptV1 {
    species: SpeciesId,
    coefficient: i128,
    activity_bits: u64,
    log_activity_bits: u64,
    direction: MassActionDirectionV1,
    log_contribution_bits: u64,
}

impl StoichiometricMassActionTermReceiptV1 {
    /// Canonical active species.
    #[must_use]
    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }

    /// Exact signed stoichiometric coefficient.
    #[must_use]
    pub const fn coefficient(&self) -> i128 {
        self.coefficient
    }

    /// Exact positive dimensionless activity bits.
    #[must_use]
    pub const fn activity_bits(&self) -> u64 {
        self.activity_bits
    }

    /// Exact deterministic log-activity bits.
    #[must_use]
    pub const fn log_activity_bits(&self) -> u64 {
        self.log_activity_bits
    }

    /// Direction selected by the coefficient sign.
    #[must_use]
    pub const fn direction(&self) -> MassActionDirectionV1 {
        self.direction
    }

    /// Exact exponent-times-log contribution bits.
    #[must_use]
    pub const fn log_contribution_bits(&self) -> u64 {
        self.log_contribution_bits
    }
}

/// Complete replay receipt for one declared stoichiometric mass-action
/// evaluation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DeclaredStoichiometricMassActionReceiptV1 {
    evaluator_version: u32,
    fs_math_version: &'static str,
    fs_qty_version: &'static str,
    kinetic_order_convention: MassActionKineticOrderConventionV1,
    reverse_rate_receipt: ThermodynamicReverseRateReceiptV1,
    term_receipts: Vec<StoichiometricMassActionTermReceiptV1>,
    log_forward_activity_product_bits: u64,
    log_reverse_activity_product_bits: u64,
    forward_status: MassActionDirectionalRateStatusV1,
    reverse_status: MassActionDirectionalRateStatusV1,
    forward_activity_product_bits: Option<u64>,
    reverse_activity_product_bits: Option<u64>,
    forward_progress_rate_bits: Option<u64>,
    reverse_progress_rate_bits: Option<u64>,
    net_progress_rate_bits: Option<u64>,
}

impl DeclaredStoichiometricMassActionReceiptV1 {
    /// Mass-action operation-tree version.
    #[must_use]
    pub const fn evaluator_version(&self) -> u32 {
        self.evaluator_version
    }

    /// Deterministic elementary-math crate version.
    #[must_use]
    pub const fn fs_math_version(&self) -> &'static str {
        self.fs_math_version
    }

    /// Quantity-system crate version.
    #[must_use]
    pub const fn fs_qty_version(&self) -> &'static str {
        self.fs_qty_version
    }

    /// Explicit kinetic-order declaration.
    #[must_use]
    pub const fn kinetic_order_convention(&self) -> MassActionKineticOrderConventionV1 {
        self.kinetic_order_convention
    }

    /// Complete nested thermodynamic reverse-rate receipt.
    #[must_use]
    pub const fn reverse_rate_receipt(&self) -> &ThermodynamicReverseRateReceiptV1 {
        &self.reverse_rate_receipt
    }

    /// Canonical per-species activity and contribution receipts.
    #[must_use]
    pub fn terms(&self) -> &[StoichiometricMassActionTermReceiptV1] {
        &self.term_receipts
    }

    /// Exact finite forward log-product bits.
    #[must_use]
    pub const fn log_forward_activity_product_bits(&self) -> u64 {
        self.log_forward_activity_product_bits
    }

    /// Exact finite reverse log-product bits.
    #[must_use]
    pub const fn log_reverse_activity_product_bits(&self) -> u64 {
        self.log_reverse_activity_product_bits
    }

    /// Forward direct representation status.
    #[must_use]
    pub const fn forward_status(&self) -> MassActionDirectionalRateStatusV1 {
        self.forward_status
    }

    /// Reverse direct representation status.
    #[must_use]
    pub const fn reverse_status(&self) -> MassActionDirectionalRateStatusV1 {
        self.reverse_status
    }

    /// Direct forward activity-product bits, when representable.
    #[must_use]
    pub const fn forward_activity_product_bits(&self) -> Option<u64> {
        self.forward_activity_product_bits
    }

    /// Direct reverse activity-product bits, when representable.
    #[must_use]
    pub const fn reverse_activity_product_bits(&self) -> Option<u64> {
        self.reverse_activity_product_bits
    }

    /// Direct forward progress-rate bits, when representable.
    #[must_use]
    pub const fn forward_progress_rate_bits(&self) -> Option<u64> {
        self.forward_progress_rate_bits
    }

    /// Direct reverse progress-rate bits, when representable.
    #[must_use]
    pub const fn reverse_progress_rate_bits(&self) -> Option<u64> {
        self.reverse_progress_rate_bits
    }

    /// Signed net progress-rate bits when both directions are representable.
    #[must_use]
    pub const fn net_progress_rate_bits(&self) -> Option<u64> {
        self.net_progress_rate_bits
    }
}

/// Successful evaluation of one declared stoichiometric mass-action law.
#[derive(Debug, Clone, PartialEq)]
pub struct DeclaredStoichiometricMassActionEvaluationV1 {
    reverse_rate: ThermodynamicReverseRateEvaluationV1,
    log_forward_activity_product: LogMassActionActivityProductV1,
    log_reverse_activity_product: LogMassActionActivityProductV1,
    forward_activity_product: Option<MassActionActivityProductV1>,
    reverse_activity_product: Option<MassActionActivityProductV1>,
    forward_progress_rate: Option<DirectionalReactionProgressRateV1>,
    reverse_progress_rate: Option<DirectionalReactionProgressRateV1>,
    net_progress_rate: Option<NetReactionProgressRateV1>,
    forward_status: MassActionDirectionalRateStatusV1,
    reverse_status: MassActionDirectionalRateStatusV1,
    receipt: DeclaredStoichiometricMassActionReceiptV1,
}

impl DeclaredStoichiometricMassActionEvaluationV1 {
    /// Nested equilibrium and reverse-scale evaluation.
    #[must_use]
    pub const fn reverse_rate(&self) -> &ThermodynamicReverseRateEvaluationV1 {
        &self.reverse_rate
    }

    /// Always-retained finite forward log activity product.
    #[must_use]
    pub const fn log_forward_activity_product(&self) -> LogMassActionActivityProductV1 {
        self.log_forward_activity_product
    }

    /// Always-retained finite reverse log activity product.
    #[must_use]
    pub const fn log_reverse_activity_product(&self) -> LogMassActionActivityProductV1 {
        self.log_reverse_activity_product
    }

    /// Direct positive finite forward activity product, when representable.
    #[must_use]
    pub const fn forward_activity_product(&self) -> Option<MassActionActivityProductV1> {
        self.forward_activity_product
    }

    /// Direct positive finite reverse activity product, when representable.
    #[must_use]
    pub const fn reverse_activity_product(&self) -> Option<MassActionActivityProductV1> {
        self.reverse_activity_product
    }

    /// Direct positive finite forward progress rate, when representable.
    #[must_use]
    pub const fn forward_progress_rate(&self) -> Option<DirectionalReactionProgressRateV1> {
        self.forward_progress_rate
    }

    /// Direct positive finite reverse progress rate, when representable.
    #[must_use]
    pub const fn reverse_progress_rate(&self) -> Option<DirectionalReactionProgressRateV1> {
        self.reverse_progress_rate
    }

    /// Signed net progress rate when both directions are directly
    /// representable.
    #[must_use]
    pub const fn net_progress_rate(&self) -> Option<NetReactionProgressRateV1> {
        self.net_progress_rate
    }

    /// Forward direct representation status.
    #[must_use]
    pub const fn forward_status(&self) -> MassActionDirectionalRateStatusV1 {
        self.forward_status
    }

    /// Reverse direct representation status.
    #[must_use]
    pub const fn reverse_status(&self) -> MassActionDirectionalRateStatusV1 {
        self.reverse_status
    }

    /// Complete nested-authority, declaration, activity, arithmetic, and
    /// result receipt.
    #[must_use]
    pub const fn receipt(&self) -> &DeclaredStoichiometricMassActionReceiptV1 {
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

    fn progress_rate_scale(value: f64) -> ReactionProgressRateScaleV1 {
        ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(value))
            .expect("positive finite progress-rate scale")
    }

    fn activity(species_id: &str, value: f64) -> IdealGasSpeciesActivityV1 {
        IdealGasSpeciesActivityV1::new(species(species_id), Dimensionless::new(value))
            .expect("positive finite activity")
    }

    fn water_mass_action(scale: f64) -> DeclaredStoichiometricMassActionV1 {
        DeclaredStoichiometricMassActionV1::new(
            ThermodynamicReverseRateClosureV1::new(
                water_equilibrium(false),
                progress_rate_scale(scale),
            ),
            MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients,
        )
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

    #[test]
    fn g0_reverse_progress_rate_closes_against_the_exact_equilibrium_receipt() {
        let forward_scale = progress_rate_scale(25.0);
        let closure =
            ThermodynamicReverseRateClosureV1::new(water_equilibrium(false), forward_scale);
        let evaluation = closure
            .evaluate(Temperature::new(1_200.0))
            .expect("finite reverse-rate closure");

        assert_eq!(
            evaluation.status(),
            ReverseProgressRateStatusV1::FinitePositive
        );
        assert_eq!(
            evaluation.log_reverse_to_forward_ratio().value().to_bits(),
            (-evaluation.equilibrium().log_equilibrium_constant().value()).to_bits()
        );
        let ratio = evaluation.reverse_to_forward_ratio().expect("direct ratio");
        let reverse_scale = evaluation
            .reverse_progress_rate_scale()
            .expect("direct reverse scale");
        assert_eq!(
            reverse_scale.value().to_bits(),
            (forward_scale.value() * ratio.value()).to_bits()
        );
        let equilibrium_constant = evaluation
            .equilibrium()
            .equilibrium_constant()
            .expect("direct equilibrium constant");
        assert!((ratio.value() * equilibrium_constant.value() - 1.0).abs() <= 16.0 * f64::EPSILON);
        assert_eq!(
            evaluation.receipt().equilibrium_receipt(),
            evaluation.equilibrium().receipt()
        );
        assert_eq!(
            evaluation.receipt().forward_progress_rate_scale_bits(),
            forward_scale.value().to_bits()
        );
        assert_eq!(
            evaluation.receipt().log_reverse_to_forward_ratio_bits(),
            evaluation.log_reverse_to_forward_ratio().value().to_bits()
        );
        assert_eq!(
            evaluation.receipt().reverse_to_forward_ratio_bits(),
            Some(ratio.value().to_bits())
        );
        assert_eq!(
            evaluation.receipt().reverse_progress_rate_scale_bits(),
            Some(reverse_scale.value().to_bits())
        );
    }

    #[test]
    fn g5_reverse_rate_replays_and_reaction_reversal_inverts_the_ratio() {
        let temperature = Temperature::new(1_200.0);
        let scale = progress_rate_scale(3.0);
        let closure = ThermodynamicReverseRateClosureV1::new(water_equilibrium(false), scale);
        let first = closure
            .evaluate(temperature)
            .expect("first reverse-rate evaluation");
        let replay = closure
            .evaluate(temperature)
            .expect("replayed reverse-rate evaluation");
        assert_eq!(first, replay);

        let reversed = ThermodynamicReverseRateClosureV1::new(water_equilibrium(true), scale)
            .evaluate(temperature)
            .expect("reversed reaction closure");
        assert_eq!(
            reversed.log_reverse_to_forward_ratio().value().to_bits(),
            (-first.log_reverse_to_forward_ratio().value()).to_bits()
        );
        let ratio_product = first
            .reverse_to_forward_ratio()
            .expect("forward direct ratio")
            .value()
            * reversed
                .reverse_to_forward_ratio()
                .expect("reversed direct ratio")
                .value();
        assert!((ratio_product - 1.0).abs() <= 16.0 * f64::EPSILON);
    }

    #[test]
    fn g3_reverse_rate_retains_log_authority_across_ratio_and_scale_range_loss() {
        let exact = 1i128 << 53;
        let (stoichiometric, certificate, reaction_id) = isomer_stoichiometry(exact);
        let ratio_underflow = ThermodynamicReverseRateClosureV1::new(
            IdealGasReactionEquilibriumV1::new(
                stoichiometric,
                certificate,
                reaction_id,
                vec![
                    model("A", 0.01, 1.0, 0.0, 100_000.0),
                    model("B", 0.01, 2.0, 0.0, 100_000.0),
                ],
            )
            .expect("large positive log equilibrium"),
            progress_rate_scale(1.0),
        )
        .evaluate(Temperature::new(500.0))
        .expect("log-only ratio underflow");
        assert_eq!(
            ratio_underflow.status(),
            ReverseProgressRateStatusV1::LogOnlyRatioUnderflow
        );
        assert!(ratio_underflow.reverse_to_forward_ratio().is_none());
        assert!(ratio_underflow.reverse_progress_rate_scale().is_none());

        let (stoichiometric, certificate, reaction_id) = isomer_stoichiometry(exact);
        let ratio_overflow = ThermodynamicReverseRateClosureV1::new(
            IdealGasReactionEquilibriumV1::new(
                stoichiometric,
                certificate,
                reaction_id,
                vec![
                    model("A", 0.01, 2.0, 0.0, 100_000.0),
                    model("B", 0.01, 1.0, 0.0, 100_000.0),
                ],
            )
            .expect("large negative log equilibrium"),
            progress_rate_scale(1.0),
        )
        .evaluate(Temperature::new(500.0))
        .expect("log-only ratio overflow");
        assert_eq!(
            ratio_overflow.status(),
            ReverseProgressRateStatusV1::LogOnlyRatioOverflow
        );
        assert!(ratio_overflow.reverse_to_forward_ratio().is_none());
        assert!(ratio_overflow.reverse_progress_rate_scale().is_none());

        let (stoichiometric, certificate, reaction_id) = isomer_stoichiometry(1);
        let scale_overflow = ThermodynamicReverseRateClosureV1::new(
            IdealGasReactionEquilibriumV1::new(
                stoichiometric,
                certificate,
                reaction_id,
                vec![
                    model("A", 0.01, 2.0, 0.0, 100_000.0),
                    model("B", 0.01, 1.0, 0.0, 100_000.0),
                ],
            )
            .expect("finite reverse-to-forward ratio above one"),
            progress_rate_scale(f64::MAX),
        )
        .evaluate(Temperature::new(500.0))
        .expect("log-only scale overflow");
        assert_eq!(
            scale_overflow.status(),
            ReverseProgressRateStatusV1::LogOnlyScaleOverflow
        );
        assert!(scale_overflow.reverse_to_forward_ratio().is_some());
        assert!(scale_overflow.reverse_progress_rate_scale().is_none());

        let (stoichiometric, certificate, reaction_id) = isomer_stoichiometry(1);
        let scale_underflow = ThermodynamicReverseRateClosureV1::new(
            IdealGasReactionEquilibriumV1::new(
                stoichiometric,
                certificate,
                reaction_id,
                vec![
                    model("A", 0.01, 1.0, 0.0, 100_000.0),
                    model("B", 0.01, 2.0, 0.0, 100_000.0),
                ],
            )
            .expect("finite reverse-to-forward ratio below one"),
            progress_rate_scale(f64::from_bits(1)),
        )
        .evaluate(Temperature::new(500.0))
        .expect("log-only scale underflow");
        assert_eq!(
            scale_underflow.status(),
            ReverseProgressRateStatusV1::LogOnlyScaleUnderflow
        );
        assert!(scale_underflow.reverse_to_forward_ratio().is_some());
        assert!(scale_underflow.reverse_progress_rate_scale().is_none());
    }

    #[test]
    fn g3_reverse_rate_refuses_invalid_forward_scale_and_nested_temperature() {
        for value in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(matches!(
                ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(value)),
                Err(ReverseRateClosureErrorV1::InvalidForwardProgressRateScale { bits })
                    if bits == value.to_bits()
            ));
        }

        let error = ThermodynamicReverseRateClosureV1::new(
            water_equilibrium(false),
            progress_rate_scale(1.0),
        )
        .evaluate(Temperature::new(0.0))
        .expect_err("nested zero-temperature equilibrium must refuse");
        assert!(matches!(
            error,
            ReverseRateClosureErrorV1::Equilibrium { .. }
        ));
    }

    #[test]
    fn g0_declared_mass_action_evaluates_canonical_forward_reverse_and_net_progress() {
        let model = water_mass_action(7.0);
        let evaluation = model
            .evaluate(
                Temperature::new(1_200.0),
                vec![
                    activity("O2", 3.0),
                    activity("H2", 2.0),
                    activity("H2O", 5.0),
                ],
            )
            .expect("finite declared mass-action evaluation");

        let expected_forward_log = 2.0 * fs_math::det::ln(2.0) + fs_math::det::ln(3.0);
        let expected_reverse_log = 2.0 * fs_math::det::ln(5.0);
        assert_eq!(
            evaluation.log_forward_activity_product().value().to_bits(),
            expected_forward_log.to_bits()
        );
        assert_eq!(
            evaluation.log_reverse_activity_product().value().to_bits(),
            expected_reverse_log.to_bits()
        );
        assert_eq!(
            evaluation.forward_status(),
            MassActionDirectionalRateStatusV1::FinitePositive
        );
        assert_eq!(
            evaluation.reverse_status(),
            MassActionDirectionalRateStatusV1::FinitePositive
        );

        let forward_product = evaluation
            .forward_activity_product()
            .expect("direct forward activity product");
        let reverse_product = evaluation
            .reverse_activity_product()
            .expect("direct reverse activity product");
        assert_eq!(
            forward_product.value().to_bits(),
            fs_math::det::exp(expected_forward_log).to_bits()
        );
        assert_eq!(
            reverse_product.value().to_bits(),
            fs_math::det::exp(expected_reverse_log).to_bits()
        );

        let forward_rate = evaluation
            .forward_progress_rate()
            .expect("direct forward rate");
        let reverse_rate = evaluation
            .reverse_progress_rate()
            .expect("direct reverse rate");
        let reverse_scale = evaluation
            .reverse_rate()
            .reverse_progress_rate_scale()
            .expect("direct reverse scale");
        assert_eq!(
            forward_rate.value().to_bits(),
            (7.0 * forward_product.value()).to_bits()
        );
        assert_eq!(
            reverse_rate.value().to_bits(),
            (reverse_scale.value() * reverse_product.value()).to_bits()
        );
        assert_eq!(
            evaluation
                .net_progress_rate()
                .expect("direct net rate")
                .value()
                .to_bits(),
            (forward_rate.value() - reverse_rate.value()).to_bits()
        );

        assert_eq!(
            evaluation.receipt().kinetic_order_convention(),
            MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients
        );
        assert_eq!(
            evaluation.receipt().reverse_rate_receipt(),
            evaluation.reverse_rate().receipt()
        );
        assert_eq!(
            evaluation
                .receipt()
                .terms()
                .iter()
                .map(|term| (
                    term.species().as_str(),
                    term.coefficient(),
                    term.direction()
                ))
                .collect::<Vec<_>>(),
            vec![
                ("H2", -2, MassActionDirectionV1::Forward),
                ("H2O", 2, MassActionDirectionV1::Reverse),
                ("O2", -1, MassActionDirectionV1::Forward),
            ]
        );
    }

    #[test]
    fn g5_declared_mass_action_permutation_replays_exactly() {
        let model = water_mass_action(2.5);
        let temperature = Temperature::new(900.0);
        let first = model
            .evaluate(
                temperature,
                vec![
                    activity("H2", 0.75),
                    activity("O2", 1.25),
                    activity("H2O", 0.5),
                ],
            )
            .expect("first mass-action evaluation");
        let replay = model
            .evaluate(
                temperature,
                vec![
                    activity("H2O", 0.5),
                    activity("H2", 0.75),
                    activity("O2", 1.25),
                ],
            )
            .expect("permuted mass-action evaluation");
        assert_eq!(first, replay);
        assert_eq!(first.receipt(), replay.receipt());
    }

    #[test]
    fn g3_declared_mass_action_refuses_activity_and_coverage_drift() {
        for value in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(matches!(
                IdealGasSpeciesActivityV1::new(species("H2"), Dimensionless::new(value)),
                Err(StoichiometricMassActionErrorV1::InvalidActivity { bits, .. })
                    if bits == value.to_bits()
            ));
        }

        let model = water_mass_action(1.0);
        assert!(matches!(
            model.evaluate(
                Temperature::new(1_000.0),
                vec![activity("H2", 1.0), activity("O2", 1.0)],
            ),
            Err(StoichiometricMassActionErrorV1::ActivityCountMismatch {
                expected: 3,
                found: 2,
            })
        ));
        assert!(matches!(
            model.evaluate(
                Temperature::new(1_000.0),
                vec![
                    activity("H2", 1.0),
                    activity("H2", 2.0),
                    activity("H2O", 1.0),
                ],
            ),
            Err(StoichiometricMassActionErrorV1::DuplicateActivitySpecies { ref species })
                if species.as_str() == "H2"
        ));
        assert!(matches!(
            model.evaluate(
                Temperature::new(1_000.0),
                vec![
                    activity("A", 1.0),
                    activity("H2O", 1.0),
                    activity("O2", 1.0),
                ],
            ),
            Err(StoichiometricMassActionErrorV1::UnexpectedSpeciesActivity { ref species })
                if species.as_str() == "A"
        ));
        assert!(matches!(
            model.evaluate(
                Temperature::new(1_000.0),
                vec![
                    activity("H2", 1.0),
                    activity("H2O", 1.0),
                    activity("X", 1.0),
                ],
            ),
            Err(StoichiometricMassActionErrorV1::MissingSpeciesActivity { ref species })
                if species.as_str() == "O2"
        ));
        assert!(matches!(
            model.evaluate(
                Temperature::new(0.0),
                vec![
                    activity("H2", 1.0),
                    activity("H2O", 1.0),
                    activity("O2", 1.0),
                ],
            ),
            Err(StoichiometricMassActionErrorV1::ReverseRateClosure {
                source: ReverseRateClosureErrorV1::Equilibrium { .. },
            })
        ));
    }

    fn isomer_mass_action(coefficient: i128, scale: f64) -> DeclaredStoichiometricMassActionV1 {
        let (stoichiometric, certificate, reaction_id) = isomer_stoichiometry(coefficient);
        let equilibrium = IdealGasReactionEquilibriumV1::new(
            stoichiometric,
            certificate,
            reaction_id,
            vec![
                model("A", 0.01, 1.0, 0.0, 100_000.0),
                model("B", 0.01, 1.0, 0.0, 100_000.0),
            ],
        )
        .expect("thermoneutral isomer equilibrium");
        DeclaredStoichiometricMassActionV1::new(
            ThermodynamicReverseRateClosureV1::new(equilibrium, progress_rate_scale(scale)),
            MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients,
        )
    }

    #[test]
    fn g3_declared_mass_action_retains_log_products_across_direct_range_loss() {
        let exact = 1i128 << 53;
        let overflow = isomer_mass_action(exact, 1.0)
            .evaluate(
                Temperature::new(500.0),
                vec![activity("A", 2.0), activity("B", 1.0)],
            )
            .expect("log-only forward product overflow");
        assert_eq!(
            overflow.forward_status(),
            MassActionDirectionalRateStatusV1::LogOnlyActivityProductOverflow
        );
        assert!(overflow.log_forward_activity_product().value().is_finite());
        assert!(overflow.forward_activity_product().is_none());
        assert!(overflow.forward_progress_rate().is_none());
        assert!(overflow.net_progress_rate().is_none());
        assert_eq!(
            overflow.reverse_status(),
            MassActionDirectionalRateStatusV1::FinitePositive
        );

        let underflow = isomer_mass_action(exact, 1.0)
            .evaluate(
                Temperature::new(500.0),
                vec![activity("A", f64::MIN_POSITIVE), activity("B", 1.0)],
            )
            .expect("log-only forward product underflow");
        assert_eq!(
            underflow.forward_status(),
            MassActionDirectionalRateStatusV1::LogOnlyActivityProductUnderflow
        );
        assert!(underflow.log_forward_activity_product().value().is_finite());
        assert!(underflow.forward_activity_product().is_none());

        let rate_overflow = isomer_mass_action(1, f64::MAX)
            .evaluate(
                Temperature::new(500.0),
                vec![activity("A", 2.0), activity("B", 1.0)],
            )
            .expect("direct forward-rate overflow status");
        assert_eq!(
            rate_overflow.forward_status(),
            MassActionDirectionalRateStatusV1::ProgressRateOverflow
        );
        assert!(rate_overflow.forward_activity_product().is_some());
        assert!(rate_overflow.forward_progress_rate().is_none());

        let rate_underflow = isomer_mass_action(1, f64::from_bits(1))
            .evaluate(
                Temperature::new(500.0),
                vec![activity("A", 0.5), activity("B", 1.0)],
            )
            .expect("direct forward-rate underflow status");
        assert_eq!(
            rate_underflow.forward_status(),
            MassActionDirectionalRateStatusV1::ProgressRateUnderflow
        );
        assert!(rate_underflow.forward_activity_product().is_some());
        assert!(rate_underflow.forward_progress_rate().is_none());

        let (stoichiometric, certificate, reaction_id) = isomer_stoichiometry(exact);
        let reverse_scale_unavailable = DeclaredStoichiometricMassActionV1::new(
            ThermodynamicReverseRateClosureV1::new(
                IdealGasReactionEquilibriumV1::new(
                    stoichiometric,
                    certificate,
                    reaction_id,
                    vec![
                        model("A", 0.01, 1.0, 0.0, 100_000.0),
                        model("B", 0.01, 2.0, 0.0, 100_000.0),
                    ],
                )
                .expect("large-log equilibrium"),
                progress_rate_scale(1.0),
            ),
            MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients,
        )
        .evaluate(
            Temperature::new(500.0),
            vec![activity("A", 1.0), activity("B", 1.0)],
        )
        .expect("direct activities with log-only reverse scale");
        assert_eq!(
            reverse_scale_unavailable.reverse_status(),
            MassActionDirectionalRateStatusV1::ProgressRateScaleUnavailable
        );
        assert!(
            reverse_scale_unavailable
                .reverse_activity_product()
                .is_some()
        );
        assert!(reverse_scale_unavailable.reverse_progress_rate().is_none());
        assert!(reverse_scale_unavailable.net_progress_rate().is_none());
    }
}
