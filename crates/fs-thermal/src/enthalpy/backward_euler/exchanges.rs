//! Explicit bath heat on an existing physical boundary channel.
use super::*;

/// Signed heat leaving a cell toward one prescribed reservoir [W], evaluated
/// at that reservoir's actual absolute temperature. Opposing contributions
/// remain distinct even when their net heat is zero.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrescribedBathHeat {
    pub outward_power_w: f64,
    pub temperature_k: f64,
}

impl PrescribedBathHeat {
    pub(crate) fn admit(self) -> Result<(), DepositionError> {
        if !self.outward_power_w.is_finite() || !positive(self.temperature_k) {
            return Err(DepositionError::InvalidInput(
                "bath exchange needs finite signed power and positive absolute temperature",
            ));
        }
        Ok(())
    }
}

/// Two explicit reservoir contributions replacing one existing face channel.
/// This is endpoint heat, not an extra cell source or an explicit-step bound.
/// A selected channel cannot also carry a conductance replacement. Surface
/// coverage, the original caloric law, mass and accepted boundary data remain
/// with the timeline. A finite body still belongs on an internal link.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BoundaryExchanges {
    pub ambient: BTreeMap<CellId, [PrescribedBathHeat; 2]>,
    pub reservoir: BTreeMap<CellId, [PrescribedBathHeat; 2]>,
}

/// Original thermal/transport coefficients plus resolved physical bath heat.
/// Existing coefficient structs keep their original fields and API.
#[derive(Debug)]
pub struct BoundaryExchangeStepCoefficients<T> {
    pub transport: TransportStepCoefficients<T>,
    pub exchanges: BoundaryExchanges,
}

impl<T> From<TransportStepCoefficients<T>> for BoundaryExchangeStepCoefficients<T> {
    fn from(transport: TransportStepCoefficients<T>) -> Self {
        Self {
            transport,
            exchanges: BoundaryExchanges::default(),
        }
    }
}

impl BoundaryExchanges {
    pub(crate) fn is_empty(&self) -> bool {
        self.ambient.is_empty() && self.reservoir.is_empty()
    }

    pub(crate) fn outward_power(baths: &[PrescribedBathHeat; 2]) -> Result<f64, DepositionError> {
        finite(baths[0].outward_power_w + baths[1].outward_power_w)
    }

    pub(crate) fn admit(
        &self,
        thermal: &DepositionSimulation<'_>,
        conductances: &BoundaryConductances,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        for (reservoir, entries) in [(false, &self.ambient), (true, &self.reservoir)] {
            if entries.len() > thermal.cells.len() {
                return Err(DepositionError::InvalidInput("too many boundary exchanges"));
            }
            for (i, (&id, baths)) in entries.iter().enumerate() {
                if i % 256 == 0 {
                    checkpoint()?;
                }
                let cell = thermal.cells.get(&id).ok_or(DepositionError::InvalidInput(
                    "boundary exchange needs an active cell",
                ))?;
                let overlaps = if reservoir {
                    if cell.boundary.reservoir.is_none() {
                        return Err(DepositionError::InvalidInput(
                            "exchange replacement needs an existing reservoir channel",
                        ));
                    }
                    conductances.reservoir_w_k.contains_key(&id)
                } else {
                    conductances.ambient_w_k.contains_key(&id)
                };
                if overlaps {
                    return Err(DepositionError::InvalidInput(
                        "one boundary channel cannot own both resolved heat and conductance",
                    ));
                }
                for bath in baths {
                    bath.admit()?;
                }
                Self::outward_power(baths)?;
            }
        }
        Ok(())
    }
}

impl<'s, 'm> BackwardEulerStep<'s, 'm> {
    /// Evaluate the original joule residual using explicit signed heat to two
    /// actual prescribed baths on each selected channel. Resolve a common
    /// physical surface before returning its separate reservoir contributions.
    pub fn evaluate_with_boundary_exchanges<T, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
        ) -> Result<BoundaryExchangeStepCoefficients<T>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, T>, E> {
        self.evaluate_exchange_option(energy_increments_j, resolve, false, None, checkpoint)
    }

    /// The same original residual with reservoir entropy from each actual
    /// signed Q/T. Net-zero heat can still transfer entropy between baths.
    pub fn evaluate_with_boundary_exchanges_and_entropy<T, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
        ) -> Result<BoundaryExchangeStepCoefficients<T>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, T>, E> {
        self.evaluate_exchange_option(energy_increments_j, resolve, true, None, checkpoint)
    }

    pub(in crate::enthalpy) fn evaluate_exchange_option<T, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
        ) -> Result<BoundaryExchangeStepCoefficients<T>, E>,
        with_entropy: bool,
        quadrature: Option<StageQuadrature<'_>>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, T>, E> {
        self.evaluate_endpoint_option(
            energy_increments_j,
            |state| resolve(state).map(|value| (value.transport, None, value.exchanges)),
            false,
            with_entropy,
            quadrature,
            checkpoint,
        )
    }
}
