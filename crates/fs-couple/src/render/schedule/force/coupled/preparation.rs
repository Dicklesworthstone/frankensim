//! Explicit construction purpose for the same physical modal operator.
//! Harmonic analysis retains sources and coefficients without a time solver.
use super::*;

/// Select construction work without changing any retained physical coordinate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ModalPreparation {
    /// Existing ordinary exact propagation, with its original island limits.
    Ordinary,
    /// Existing prepared held-force propagation and integrated loss channels.
    Prepared(PreparedModalBudget),
    /// Continuous harmonic queries only; no time exponential or contact
    /// compliance is constructed. The declared budget bounds mode count,
    /// coefficient setup and dense harmonic workspace. Time stepping refuses.
    HarmonicOnly(PreparedModalBudget),
    /// The complete bank of [`Self::HarmonicOnly`], admitted as [`Self::Prepared`] admits it. Every construction
    /// and rebuild runs what a prepared construction runs before it computes its operator (the held-load pressure
    /// preparation, its reserve, and the operator's mode, byte and work admission) and refuses as it
    /// does. No time exponential is computed: a cold system built so gains it once, from
    /// `CoupledModalSystem::with_time_propagation`.
    StructureFirst(PreparedModalBudget),
}

impl From<Option<PreparedModalBudget>> for ModalPreparation {
    fn from(value: Option<PreparedModalBudget>) -> Self {
        value.map_or(Self::Ordinary, Self::Prepared)
    }
}

impl ModalPreparation {
    #[must_use]
    pub const fn prepared_budget(self) -> Option<PreparedModalBudget> {
        match self {
            Self::Ordinary => None,
            Self::Prepared(budget) | Self::HarmonicOnly(budget) | Self::StructureFirst(budget) => {
                Some(budget)
            }
        }
    }

    #[must_use]
    pub const fn is_harmonic_only(self) -> bool {
        matches!(self, Self::HarmonicOnly(_) | Self::StructureFirst(_))
    }
}

/// One shared preparation choice used by constructors and every coefficient,
/// material-history and exterior-load rebuild. Analysis never creates a fake
/// or empty propagator that a time path could accidentally execute.
pub(super) struct PreparedOperator {
    pub(super) exact: Option<exact::ExactNetwork>,
    pub(super) harmonic_damping_columns: Vec<Option<Vec<f64>>>,
    pub(super) pressure: equilibrium::pressure::EquilibriumPressure,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare(
    models: &[ModalAcousticTimeModel],
    columns: &[Vec<f64>],
    links: &[ConnectionLaw],
    count: usize,
    dt: f64,
    choice: ModalPreparation,
    branches: &[fs_phs::RelaxationBranch],
    radiation: Option<&radiation::Binding>,
    config: ModalCouplingConfig,
    sums_by_rows: bool,
    gate: &CancelGate,
) -> Result<PreparedOperator, ModalCouplingError> {
    poll(Some(gate))?;
    if let ModalPreparation::StructureFirst(budget) = choice {
        // A prepared construction's refusals before its operator is computed, in its order.
        let held = equilibrium::pressure::EquilibriumPressure::prepare;
        let held = held(models, columns, links, config, Some(budget), gate)?;
        let operator = held.reserve(budget)?;
        operator.admit_structure(models, links, count, branches.len(), radiation)?;
    }
    let pressure = if choice.is_harmonic_only() {
        equilibrium::pressure::EquilibriumPressure::default()
    } else {
        equilibrium::pressure::EquilibriumPressure::prepare(
            models, columns, links, config, choice.prepared_budget(), gate,
        )?
    };
    let mut exact = match choice {
        ModalPreparation::HarmonicOnly(budget) | ModalPreparation::StructureFirst(budget) => {
            admit_harmonic(
                budget,
                count,
                columns.len(),
                branches.len(),
                links,
                radiation.map_or(0, |r| r.model.poles.len()),
            )?;
            let harmonic_damping_columns =
                exact::ExactNetwork::harmonic_damping_columns(models, columns, links, count, gate)?;
            poll(Some(gate))?;
            return Ok(PreparedOperator {
                exact: None,
                harmonic_damping_columns,
                pressure,
            });
        }
        ModalPreparation::Prepared(budget) => {
            Some(exact::ExactNetwork::new_prepared_with_relaxation(
                models, columns, links, count, dt, pressure.reserve(budget)?, branches, radiation, gate,
            )?)
        }
        ModalPreparation::Ordinary
            if links
                .iter()
                .any(|link| link.stiffness_n_m != 0.0 || link.damping_n_s_m != 0.0)
                || !branches.is_empty()
                || radiation.is_some() =>
        {
            Some(exact::ExactNetwork::new_with_extensions(
                models, columns, links, count, branches, radiation, gate,
            )?)
        }
        ModalPreparation::Ordinary => None,
    };
    // A system that sums by rows keeps that order in every operator prepared for it, and each one's loss Gramians
    // are held to their transposes again.
    if sums_by_rows {
        exact
            .as_mut()
            .ok_or_else(|| invalid("sums by rows need a time-domain operator"))?
            .sum_by_rows()?;
    }
    Ok(PreparedOperator {
        exact,
        harmonic_damping_columns: Vec::new(),
        pressure,
    })
}

fn admit_harmonic(
    budget: PreparedModalBudget,
    count: usize,
    columns: usize,
    memories: usize,
    links: &[ConnectionLaw],
    poles: usize,
) -> Result<(), ModalCouplingError> {
    budget.validate()?;
    let dimension = count
        .checked_add(poles)
        .ok_or_else(|| invalid("harmonic construction dimension overflow"))?;
    if count == 0 || count > budget.max_modes || dimension > harmonic::MAX_HARMONIC_UNKNOWNS {
        return Err(invalid(
            "complete harmonic construction exceeds its mode/air capacity",
        ));
    }
    let anchored = links
        .iter()
        .filter(|link| link.constant_loss_reference_omega_rad_s.is_some())
        .count();
    let square = dimension
        .checked_mul(dimension)
        .ok_or_else(|| invalid("harmonic construction matrix extent overflow"))?;
    let mechanical_square = count
        .checked_mul(count)
        .ok_or_else(|| invalid("harmonic coefficient matrix extent overflow"))?;
    // Four complex squares conservatively cover original/scaled/factored
    // operators and solve scratch. Anchored loss additionally retains the
    // original stiffness eigensolve's four real squares and damping columns.
    // Histories stay physical scalars; they are not exponentiated time states.
    let entries = square
        .checked_mul(8)
        .and_then(|v| v.checked_add(dimension.checked_mul(16)?))
        .and_then(|v| v.checked_add(memories.checked_mul(8)?))
        .and_then(|v| {
            v.checked_add(if anchored == 0 {
                0
            } else {
                mechanical_square.checked_mul(4)?
            })
        })
        .and_then(|v| v.checked_add(anchored.checked_mul(count)?))
        .ok_or_else(|| invalid("harmonic construction byte extent overflow"))?;
    let bytes = entries
        .checked_mul(core::mem::size_of::<f64>())
        .and_then(|v| v.checked_add(columns.checked_mul(core::mem::size_of::<Option<Vec<f64>>>())?))
        .ok_or_else(|| invalid("harmonic construction byte overflow"))?;
    let count64 = count as u64;
    let terms = count64
        .checked_mul(columns as u64)
        .and_then(|v| v.checked_add(count64.checked_mul(memories as u64)?))
        .and_then(|v| {
            v.checked_add(if anchored == 0 {
                0
            } else {
                count64.checked_pow(3)?.checked_mul(9)?.checked_add(
                    count64
                        .checked_pow(2)?
                        .checked_mul((columns + 2 * anchored) as u64)?,
                )?
            })
        })
        .ok_or_else(|| invalid("harmonic coefficient setup work overflow"))?;
    if bytes > budget.max_dense_bytes || terms > budget.max_multiply_terms {
        return Err(invalid(
            "harmonic coefficient/workspace preparation exceeds its explicit byte/work budget",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
