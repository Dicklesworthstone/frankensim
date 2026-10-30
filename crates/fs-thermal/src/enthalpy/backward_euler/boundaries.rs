//! Endpoint coefficients for the original prescribed thermal boundaries.
use super::*;

/// Total nonnegative conductances in W/K at one proposed physical state.
/// Omitted cells retain their authored coefficients. These maps carry no
/// temperature, surface coverage, storage or heat source; those original owners
/// remain authoritative. In particular, a finite thermal body is an internal
/// link, never converted into a prescribed infinite reservoir by this ingress.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BoundaryConductances {
    /// Total free-surface Robin G, replacing reference exposed area * authored h.
    pub ambient_w_k: BTreeMap<CellId, f64>,
    /// Total G of an already-present prescribed reservoir; its temperature and
    /// covered reference area do not change. An absent reservoir refuses.
    pub reservoir_w_k: BTreeMap<CellId, f64>,
}

/// Ordinary thermal coefficients and explicit physical-boundary replacements.
/// Kept separate so existing source/link callbacks preserve their meaning.
#[derive(Debug)]
pub struct BoundaryStepCoefficients<T> {
    pub bulk: StepCoefficients<T>,
    pub boundaries: BoundaryConductances,
}
impl<T> From<StepCoefficients<T>> for BoundaryStepCoefficients<T> {
    fn from(bulk: StepCoefficients<T>) -> Self {
        Self {
            bulk,
            boundaries: BoundaryConductances::default(),
        }
    }
}

impl BoundaryConductances {
    pub(crate) fn admit(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        if self.ambient_w_k.len() > thermal.cells.len()
            || self.reservoir_w_k.len() > thermal.cells.len()
        {
            return Err(DepositionError::InvalidInput(
                "too many boundary conductance replacements",
            ));
        }
        for (i, (&id, &g)) in self.ambient_w_k.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            if !thermal.cells.contains_key(&id) || !nonnegative(g) {
                return Err(DepositionError::InvalidInput(
                    "ambient replacement needs an active cell and finite nonnegative G",
                ));
            }
        }
        for (i, (&id, &g)) in self.reservoir_w_k.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            if thermal
                .cells
                .get(&id)
                .and_then(|cell| cell.boundary.reservoir)
                .is_none()
                || !nonnegative(g)
            {
                return Err(DepositionError::InvalidInput(
                    "reservoir replacement needs an existing reservoir and finite nonnegative G",
                ));
            }
        }
        Ok(())
    }

    pub(crate) fn ambient(
        &self,
        id: CellId,
        reference_area_m2: f64,
        authored_htc: f64,
    ) -> Result<f64, DepositionError> {
        if let Some(&g) = self.ambient_w_k.get(&id) {
            if reference_area_m2 == 0.0 && g != 0.0 {
                return Err(DepositionError::InvalidInput(
                    "nonzero Robin conductance needs an exposed reference patch",
                ));
            }
            Ok(g)
        } else {
            finite(reference_area_m2 * authored_htc)
        }
    }
}
