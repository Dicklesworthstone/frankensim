//! Closed-component thermal response without grounding an energy zero mode.
//!
//! The tangent, caloric slopes, complex LU and original cell-power residuals
//! are the native harmonic owner's. One redundant power row per original
//! conducting component is replaced by its extensive stored-energy equation.
//! At nonzero frequency this is exactly sum(E)=i*sum(P)/omega; at DC it is the
//! caller's inventory increment (zero by default), and net heat must vanish.
//! Thus a balanced heat transfer has a finite DC response while net heating
//! retains its physical storage pole. No reservoir or time step is added.
//!
//! Energy means the original declared storage coordinate, including enthalpy
//! where that is the material's convention. This is a local linear response:
//! callers retain the admitted differentiable branch and small amplitudes.
//! Several latent parcels leave undetermined DC phase allocations and refuse;
//! their finite-frequency stored-energy response is still supported.

use super::conservation::ConservedThermalEnergy;
pub use super::conservation::{
    IsolatedThermalImpedanceControls, ThermalComponentResponse, ThermalEnergyConstraintError,
};
use super::*;

#[derive(Clone, Debug, PartialEq)]
pub enum IsolatedThermalImpedanceError {
    Original(ThermalImpedanceError),
    Conservation(ThermalEnergyConstraintError<DepositionError>),
    CellResidual { maximum_w: f64, tolerance_w: f64 },
}
impl From<ThermalImpedanceError> for IsolatedThermalImpedanceError {
    fn from(value: ThermalImpedanceError) -> Self {
        Self::Original(value)
    }
}
impl From<ThermalEnergyConstraintError<DepositionError>> for IsolatedThermalImpedanceError {
    fn from(value: ThermalEnergyConstraintError<DepositionError>) -> Self {
        Self::Conservation(value)
    }
}
impl From<DepositionError> for IsolatedThermalImpedanceError {
    fn from(value: DepositionError) -> Self {
        Self::Original(value.into())
    }
}
impl core::fmt::Display for IsolatedThermalImpedanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "isolated thermal impedance: {self:?}")
    }
}
impl std::error::Error for IsolatedThermalImpedanceError {}
type IsolatedResult<T> = core::result::Result<T, IsolatedThermalImpedanceError>;

/// Prepared only from a source-free isolated equilibrium and its original
/// component partition; the ordinary impedance API keeps its DC refusals.
pub struct IsolatedThermalImpedanceModel {
    original: ThermalImpedanceModel,
    energy: ConservedThermalEnergy,
    components: Vec<Vec<CellId>>,
    cell_power_tolerance_w: f64,
}

pub struct IsolatedThermalImpedanceFactor<'a> {
    original: ThermalImpedanceFactor<'a>,
    energy: &'a ConservedThermalEnergy,
    cell_power_tolerance_w: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct IsolatedThermalImpedanceResponse {
    /// Same original cell order, temperatures, energy and full power residuals.
    pub field: ThermalImpedanceResponse,
    /// The model's original component order; no inter-component inventory flow.
    pub components: Vec<ThermalComponentResponse>,
}

impl IsolatedThermalImpedanceModel {
    pub(crate) fn new(
        original: ThermalImpedanceModel,
        components: Vec<Vec<CellId>>,
        controls: IsolatedThermalImpedanceControls,
        checkpoint: &mut impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> IsolatedResult<Self> {
        checkpoint()?;
        let n = original.ids.len() as u128;
        // The ordinary bound already covers the shared dense factor. Include
        // component indices, compensated sums, RHS and response metadata here.
        if 64 * n * n + 512 * n > original.budget.maximum_dense_bytes as u128 {
            return Err(invalid("closed thermal response needs 64*n*n+512*n bytes").into());
        }
        let mut indices = Vec::with_capacity(components.len());
        for cells in &components {
            checkpoint()?;
            indices.push(
                cells
                    .iter()
                    .map(|id| {
                        original.ids.binary_search(id).map_err(|_| {
                            invalid("closed thermal component contains an unknown original cell")
                        })
                    })
                    .collect::<Result<Vec<_>>>()?,
            );
        }
        let energy = ConservedThermalEnergy::new(
            indices,
            &original.temperature_per_energy_k_j,
            &original.outward_tangent_per_s,
            controls,
            checkpoint,
        )?;
        Ok(Self {
            original, energy, components,
            cell_power_tolerance_w: controls.cell_power_tolerance_w,
        })
    }

    pub fn original(&self) -> &ThermalImpedanceModel {
        &self.original
    }
    pub fn components(&self) -> &[Vec<CellId>] {
        &self.components
    }

    /// Reuses the native LU, with exact component storage rows at every
    /// frequency. Multiple latent DC freedoms refuse instead of choosing phases.
    pub fn factor(
        &self,
        omega_rad_s: f64,
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> IsolatedResult<IsolatedThermalImpedanceFactor<'_>> {
        let matrix = self.energy.matrix(
            &self.original.outward_tangent_per_s,
            omega_rad_s,
            &mut checkpoint,
        )?;
        let original = self.original.factor_scaled_matrix(
            omega_rad_s,
            matrix.row_scales,
            matrix.values,
            &mut checkpoint,
        )?;
        Ok(IsolatedThermalImpedanceFactor {
            original, energy: &self.energy,
            cell_power_tolerance_w: self.cell_power_tolerance_w,
        })
    }
}

impl IsolatedThermalImpedanceFactor<'_> {
    pub fn omega_rad_s(&self) -> f64 {
        self.original.omega_rad_s()
    }

    /// Original cell heat phasors. At DC every component must have zero net
    /// heating and retains its original inventory. At finite frequency its
    /// total energy follows the physical -i*omega storage equation.
    pub fn response(
        &self,
        heat_w: &[C64],
        checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> IsolatedResult<IsolatedThermalImpedanceResponse> {
        self.solve(heat_w, None, checkpoint)
    }

    /// A DC tangent between explicitly different component inventories. These
    /// joule increments do not advance a trajectory or credit applied heat/work.
    /// Nonzero frequency refuses an independently prescribed energy increment.
    pub fn response_with_component_energy(
        &self,
        heat_w: &[C64],
        component_energy_j: &[C64],
        checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> IsolatedResult<IsolatedThermalImpedanceResponse> {
        self.solve(heat_w, Some(component_energy_j), checkpoint)
    }

    fn solve(
        &self,
        heat_w: &[C64],
        component_energy_j: Option<&[C64]>,
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> IsolatedResult<IsolatedThermalImpedanceResponse> {
        let drive = self.energy.drive(
            self.original.omega_rad_s(),
            heat_w,
            component_energy_j,
            &mut checkpoint,
        )?;
        let field = self.original.response_with_rhs(heat_w, &drive.rhs, &mut checkpoint)?;
        // A large common-mode energy near the storage pole must not hide a
        // finite watt defect behind a small relative backward error.
        if field.maximum_residual_w > self.cell_power_tolerance_w {
            return Err(IsolatedThermalImpedanceError::CellResidual {
                maximum_w: field.maximum_residual_w,
                tolerance_w: self.cell_power_tolerance_w,
            });
        }
        let components = self.energy.admit(&drive, &field.energy_j, &mut checkpoint)?;
        checkpoint()?;
        Ok(IsolatedThermalImpedanceResponse { field, components })
    }

    /// Unit heat at one original cell, K/W. At DC a nonzero net component drive
    /// has no steady response and refuses. Finite frequencies retain storage.
    pub fn impedance_column(
        &self,
        source: CellId,
        checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> IsolatedResult<Vec<C64>> {
        let column = self.original.model.ids.binary_search(&source)
            .map_err(|_| invalid("unknown thermal impedance source cell"))?;
        let mut heat = vec![C64::ZERO; self.original.model.ids.len()];
        heat[column] = C64::ONE;
        Ok(self.response(&heat, checkpoint)?.field.temperature_k)
    }

    /// Unit heat transfer from sink to source, K/W. A DC transfer must remain
    /// inside one closed component; disconnected inventories cannot compensate.
    pub fn heat_transfer_impedance(
        &self,
        source: CellId,
        sink: CellId,
        checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> IsolatedResult<Vec<C64>> {
        if source == sink {
            return Err(invalid("heat transfer needs two distinct original cells").into());
        }
        let locate = |id| self.original.model.ids.binary_search(&id)
            .map_err(|_| invalid("unknown thermal impedance transfer cell"));
        let mut heat = vec![C64::ZERO; self.original.model.ids.len()];
        heat[locate(source)?] = C64::ONE;
        heat[locate(sink)?] = -C64::ONE;
        Ok(self.response(&heat, checkpoint)?.field.temperature_k)
    }
}
