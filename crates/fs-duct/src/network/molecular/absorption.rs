//! Positive finite-volume absorption of original molecular network cylinders.

use super::*;
use crate::network::field::NetworkCylinderAbsorptionCell;
use crate::network::field::absorption::{CylinderLossOperator, admit_cells, observe_channels};
use fs_math::c64::C64;

/// Cycle-mean physical losses integrated over one original cylindrical cell.
/// Coherent peak p/U amplitudes are combined before the shared positive forms.
/// No surface load, numerical stabilization or signed power residual is heat.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MolecularNetworkCylinderAbsorption {
    /// Original wall-shear dissipation [W]; a direct caloric heat source.
    pub viscous_power_w: f64,
    /// Fourier thermal availability loss [W], not an extra caloric source.
    pub thermal_power_w: f64,
    /// Original molecular-reservoir availability loss [W]. It does not itself
    /// assign mean caloric deposition or transport between gas and walls.
    pub molecular_power_w: f64,
    /// Sum of the three independently positive channels [W].
    pub absorbed_power_w: f64,
}

/// Bound to the same source-bearing network, frequency, geometry probes and
/// complete coherent transfer bank. Preparation uses the original cached Cp,
/// local wave and radial spectra. Warm queries perform source contractions
/// and positive quadratic forms only; no state, mesh or material is refreshed.
#[derive(Debug)]
pub struct PreparedMolecularNetworkCylinderAbsorption<'fields, 'network> {
    fields: &'fields PreparedMolecularNetworkFields<'network>,
    cells: Vec<NetworkCylinderAbsorptionCell>,
    operators: Vec<CylinderLossOperator<3>>,
}

impl<'network> PreparedMolecularNetworkFields<'network> {
    /// Prepare exact cycle-mean loss integrals between existing probes in the
    /// same original cylinder, with physical cell volumes retained. Equal
    /// endpoints are zero; overlapping cells are never automatically summed.
    ///
    /// The original circular Foster resistance supplies viscosity. Its thermal
    /// spectral measure acts on the SAME complex Cp to separate Fourier and
    /// reservoir availability losses. No equilibrium LossModel, extra bulk
    /// attenuation, power subtraction or radial/axial resampling is used.
    /// Native positive two-wave integration retains standing-wave pressure and
    /// velocity nodes and remains stable for opaque intervals.
    ///
    /// # Errors
    /// Original-cylinder identity/order, nonfinite volume or positive loss
    /// forms, work/cell bounds and cancellation. Cone, profile, wall-mechanical
    /// and exterior heat need their own volume laws and refuse here. Optional
    /// absorption refusal does not invalidate an otherwise admitted p/U field.
    pub fn prepare_cylinder_absorption(
        &self,
        probe_pairs: &[[usize; 2]],
        maximum_cells: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularNetworkCylinderAbsorption<'_, 'network>, DuctError> {
        poll(&mut cancelled, 0)?;
        admit_cells(probe_pairs.len(), maximum_cells)?;
        let mut cells = Vec::with_capacity(probe_pairs.len());
        let mut operators = Vec::with_capacity(probe_pairs.len());
        // At most one coefficient evaluation per selected original cylinder,
        // bounded by the same caller cell cap. No material query occurs here.
        let mut coefficients = BTreeMap::new();
        for (index, &pair) in probe_pairs.iter().enumerate() {
            poll(&mut cancelled, index)?;
            let cell = self.fields().cylinder_absorption_cell(pair)?;
            let length = cell.end_m - cell.start_m;
            let operator = if length == 0.0 {
                CylinderLossOperator::ZERO
            } else {
                let id = self.fields().locations()[pair[0]].probe.edge_id.as_str();
                let (wave, impedance, weights) =
                    if let Some(&value) = coefficients.get(&(id, cell.segment)) {
                        value
                    } else {
                        let original = self.cylinder(id, cell.segment)?;
                        let area = core::f64::consts::PI * original.radius_m * original.radius_m;
                        // Match the actual-area normalization in TwoPort exactly.
                        let impedance = original.wave.specific_impedance / C64::from_re(area);
                        let value = (original.wave, impedance, original.loss_weights()?);
                        coefficients.insert((id, cell.segment), value);
                        value
                    };
                CylinderLossOperator::new(wave.wavenumber, impedance, length, weights)?
            };
            cells.push(cell);
            operators.push(operator);
        }
        poll(&mut cancelled, cells.len())?;
        Ok(PreparedMolecularNetworkCylinderAbsorption {
            fields: self,
            cells,
            operators,
        })
    }
}

impl<'fields, 'network> PreparedMolecularNetworkCylinderAbsorption<'fields, 'network> {
    pub fn fields(&self) -> &'fields PreparedMolecularNetworkFields<'network> {
        self.fields
    }
    pub fn cells(&self) -> &[NetworkCylinderAbsorptionCell] {
        &self.cells
    }

    /// Combine a complete coherent source vector before taking loss squares.
    pub fn sample(
        &self,
        excitation: &[C64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<MolecularNetworkCylinderAbsorption>, DuctError> {
        self.observe_coherence(Some(excitation), &[], cancelled)
    }

    /// Cycle-mean losses for an original source mean and covariance factors.
    /// Each factor retains its complete internal phase/interference; distinct
    /// factors add positive powers. This needs neither a Gaussian nor a proper
    /// process assumption. No cross-frequency or finite-shutter heat average
    /// is inferred. Source units, order, factor caps and cancellation belong to
    /// the existing network contraction owner, shared with equilibrium cells.
    pub fn observe_coherence(
        &self,
        mean: Option<&[C64]>,
        factors: &[&[C64]],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<MolecularNetworkCylinderAbsorption>, DuctError> {
        observe_channels(
            self.fields.fields(),
            &self.cells,
            &self.operators,
            mean,
            factors,
            cancelled,
            |watts, total| MolecularNetworkCylinderAbsorption {
                viscous_power_w: watts[0],
                thermal_power_w: watts[1],
                molecular_power_w: watts[2],
                absorbed_power_w: total,
            },
        )
    }
}
