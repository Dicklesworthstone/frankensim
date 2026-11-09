//! Common original-layer evolution with Maxwell and Kelvin–Maxwell constraints.

use super::{
    C64, Cx, Dd, MaterialPlane, TensorPronyHistory, TensorPronyShell,
    TensorShellHarmonicError, TensorShellMaterialCausalBudget,
    TensorShellMaterialCausalField, TensorShellMaterialCausalHistory,
    TensorShellMaterialCausalObserver, TensorShellMaterialCausalPoint,
    TensorShellMaterialCausalStep, checkpoint, invalid,
};
use fs_material::visco::TensorPronyHarmonicWindowPoint;
use fs_material::visco_tensor::{
    TensorPronyError,
    kelvin_control::{
        PreparedTensorKelvinPlaneStress, TensorKelvinPlaneStressHarmonicWindow,
        TensorKelvinStressControlBudget,
    },
    stress_control::{
        PreparedTensorMixedControl, TensorPlaneStressHarmonicBudget,
        TensorPlaneStressHarmonicWindow,
    },
};

/// Original physical material energy densities [J/m³]. The combined dissipation
/// is independently integrated even when the Kelvin/Maxwell partition is not
/// available from the coupled realization. An absent channel is never zero heat.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorShellViscoelasticEnergyDensity {
    pub start_stored_energy_j_m3: f64,
    pub end_stored_energy_j_m3: f64,
    pub stored_energy_change_j_m3: f64,
    pub material_heat_j_m3: f64,
    pub stress_work_j_m3: f64,
    pub energy_balance_residual_j_m3: f64,
    pub kelvin_heat_j_m3: Option<f64>,
    pub maxwell_heat_j_m3: Option<f64>,
}

impl From<TensorPronyHarmonicWindowPoint> for TensorShellViscoelasticEnergyDensity {
    fn from(value: TensorPronyHarmonicWindowPoint) -> Self {
        Self {
            start_stored_energy_j_m3: value.start_stored_energy_j_m3,
            end_stored_energy_j_m3: value.end_stored_energy_j_m3,
            stored_energy_change_j_m3: value.stored_energy_change_j_m3,
            material_heat_j_m3: value.material_heat_j_m3,
            stress_work_j_m3: value.stress_work_j_m3,
            energy_balance_residual_j_m3: value.energy_balance_residual_j_m3,
            kelvin_heat_j_m3: Some(value.kelvin_heat_j_m3),
            maxwell_heat_j_m3: Some(value.maxwell_heat_j_m3),
        }
    }
}

/// One original positive shell integration volume and its complete material state.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorShellViscoelasticCausalPoint {
    pub element: usize,
    pub layer: usize,
    pub triangle_quadrature_point: usize,
    pub thickness_quadrature_point: usize,
    pub triangle_barycentric: [f64; 3],
    pub thickness_position_m: f64,
    pub volume_m3: f64,
    /// Ordered engineering [epsilon_aa, epsilon_bb, gamma_ab] in layer axes.
    pub in_plane_strain: Vec<[C64; 3]>,
    pub in_plane_bias: [f64; 3],
    /// Original tensor-shear Voigt strain and stress, including transverse response.
    pub end_strain: [f64; 6],
    pub end_stress: [f64; 6],
    /// The Kelvin constraint supplies its actual full endpoint strain rate [1/s].
    /// The unchanged Maxwell-only controller does not publish this quantity.
    pub end_strain_rate: Option<[f64; 6]>,
    pub strain_bounds: [[f64; 6]; 2],
    pub material: TensorShellViscoelasticEnergyDensity,
}

/// Independent interval accounts on original layer volumes. Channel partitions
/// are present only when every physical point supplies that complete partition.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorShellViscoelasticCausalField {
    pub start_time_s_dd: Dd,
    pub end_time_s_dd: Dd,
    pub duration_s: f64,
    pub angular_frequencies_rad_s: Vec<f64>,
    pub points: Vec<TensorShellViscoelasticCausalPoint>,
    pub start_stored_energy_j: f64,
    pub end_stored_energy_j: f64,
    pub stored_energy_change_j: f64,
    pub kelvin_heat_j: Option<f64>,
    pub maxwell_heat_j: Option<f64>,
    pub material_heat_j: f64,
    pub stress_work_j: f64,
    pub energy_balance_residual_j: f64,
}

/// Both the original retained memory and its physical receipt are unpublished
/// candidates until the caller accepts the complete result.
pub struct TensorShellViscoelasticCausalStep<'a> {
    pub history: TensorShellMaterialCausalHistory<'a>,
    pub field: TensorShellViscoelasticCausalField,
}

/// One shared CST/DKT causal engine for layers with and without Kelvin viscosity.
/// It retains the same private twelve-history spatial bank and complete nodal
/// continuity as the Maxwell-only observer, with no alternative shell theory.
pub struct TensorShellViscoelasticCausalObserver<'a> {
    inner: TensorShellMaterialCausalObserver<'a>,
}

impl TensorPronyShell {
    /// Prepare causal harmonic material evolution including original Kelvin loss.
    /// Each nonzero Kelvin layer must have a numerically resolved positive-definite
    /// transverse viscosity block. The full original viscosity may be semidefinite;
    /// neither null Maxwell directions nor anisotropic cross terms are discarded.
    ///
    /// The two existing controller caps apply to their respective layer types.
    /// The field cap covers the common engine and both accepted/candidate banks.
    /// Kelvin–Maxwell heat is reported as its actual combined physical quantity;
    /// no unsupported channel partition or temperature feedback is inferred.
    pub fn prepare_viscoelastic_material_causal_observer(
        &self,
        budget: TensorShellMaterialCausalBudget,
        kelvin_budget: TensorKelvinStressControlBudget,
        cx: &Cx<'_>,
    ) -> Result<TensorShellViscoelasticCausalObserver<'_>, TensorShellHarmonicError> {
        Ok(TensorShellViscoelasticCausalObserver {
            inner: self.prepare_causal_material(budget, Some(kelvin_budget), cx)?,
        })
    }
}

impl<'a> TensorShellViscoelasticCausalObserver<'a> {
    pub const fn source(&self) -> &'a TensorPronyShell {
        self.inner.source()
    }

    pub fn start_relaxed(
        &self,
        start: Dd,
        cx: &Cx<'_>,
    ) -> Result<TensorShellMaterialCausalHistory<'a>, TensorShellHarmonicError> {
        self.inner.start_relaxed(start, cx)
    }

    /// Evolve the original layer memories under a continuous prescribed field.
    /// Carriers may change between intervals. An empty bank holds displacement
    /// while Kelvin/Maxwell relaxation continues under zero transverse stress.
    /// Supports remain stationary unless their coordinates were retained in the
    /// original source. No force-equilibrium or thermal-recoupling claim is made.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_window(
        &self,
        accepted: &TensorShellMaterialCausalHistory<'a>,
        angular_frequencies_rad_s: &[f64],
        displacements: &[&[C64]],
        bias: &[f64],
        duration_s: f64,
        cx: &Cx<'_>,
    ) -> Result<TensorShellViscoelasticCausalStep<'a>, TensorShellHarmonicError> {
        self.inner.advance_viscoelastic_window(
            accepted, angular_frequencies_rad_s, displacements, bias, duration_s, cx,
        )
    }
}

pub(super) enum LayerControl {
    Maxwell(PreparedTensorMixedControl),
    Kelvin(PreparedTensorKelvinPlaneStress),
}

pub(super) enum LayerWindow<'a> {
    Maxwell(TensorPlaneStressHarmonicWindow<'a>),
    Kelvin(TensorKelvinPlaneStressHarmonicWindow<'a>),
}

pub(super) struct PointCandidate {
    pub history: TensorPronyHistory,
    pub material: TensorShellViscoelasticEnergyDensity,
    pub strain_bounds: [[f64; 6]; 2],
    pub strain_rate: Option<[f64; 6]>,
}

impl LayerControl {
    pub(super) fn prepare_window(
        &self,
        frequencies: &[f64],
        start: Dd,
        duration: f64,
        budget: TensorPlaneStressHarmonicBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<LayerWindow<'_>, TensorPronyError> {
        match self {
            Self::Maxwell(control) => control.prepare_plane_stress_harmonic_window(
                MaterialPlane::Axes12, frequencies, start, duration, budget, cancelled,
            ).map(LayerWindow::Maxwell),
            Self::Kelvin(control) => control.prepare_harmonic_window(
                frequencies, start, duration, budget, cancelled,
            ).map(LayerWindow::Kelvin),
        }
    }
}

impl LayerWindow<'_> {
    pub(super) fn evaluate(
        &self,
        strains: &[[C64; 3]],
        bias: &[f64; 3],
        history: &TensorPronyHistory,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PointCandidate, TensorPronyError> {
        match self {
            Self::Maxwell(window) => {
                let point = window.evaluate(strains, bias, history, cancelled)?;
                Ok(PointCandidate {
                    history: point.history,
                    material: point.material.into(),
                    strain_bounds: point.strain_bounds,
                    strain_rate: None,
                })
            }
            Self::Kelvin(window) => {
                let point = window.evaluate(strains, bias, history, cancelled)?;
                let account = point.account;
                Ok(PointCandidate {
                    history: point.history,
                    material: TensorShellViscoelasticEnergyDensity {
                        start_stored_energy_j_m3: account.start_stored_energy_j_m3,
                        end_stored_energy_j_m3: account.end_stored_energy_j_m3,
                        stored_energy_change_j_m3: account.stored_energy_change_j_m3,
                        material_heat_j_m3: account.material_heat_j_m3,
                        stress_work_j_m3: account.stress_work_j_m3,
                        energy_balance_residual_j_m3: account.energy_balance_residual_j_m3,
                        kelvin_heat_j_m3: None,
                        maxwell_heat_j_m3: None,
                    },
                    strain_bounds: point.strain_bounds,
                    strain_rate: Some(point.strain_rate),
                })
            }
        }
    }
}

impl<'a> TensorShellViscoelasticCausalStep<'a> {
    // Keep the existing Maxwell-only public result and its complete channel
    // fields unchanged. The common engine performs no second material evaluation.
    pub(super) fn into_maxwell(self, cx: &Cx<'_>) -> Result<TensorShellMaterialCausalStep<'a>, TensorShellHarmonicError> {
        let channel = |value: Option<f64>| value.ok_or_else(|| invalid(
            "Maxwell-only shell receipt requires a complete heat partition",
        ));
        let field = self.field;
        let points = field.points.into_iter().map(|point| {
            checkpoint(cx)?;
            let value = point.material;
            Ok(TensorShellMaterialCausalPoint {
                element: point.element,
                layer: point.layer,
                triangle_quadrature_point: point.triangle_quadrature_point,
                thickness_quadrature_point: point.thickness_quadrature_point,
                triangle_barycentric: point.triangle_barycentric,
                thickness_position_m: point.thickness_position_m,
                volume_m3: point.volume_m3,
                in_plane_strain: point.in_plane_strain,
                in_plane_bias: point.in_plane_bias,
                end_strain: point.end_strain,
                end_stress: point.end_stress,
                strain_bounds: point.strain_bounds,
                material: TensorPronyHarmonicWindowPoint {
                    start_stored_energy_j_m3: value.start_stored_energy_j_m3,
                    end_stored_energy_j_m3: value.end_stored_energy_j_m3,
                    stored_energy_change_j_m3: value.stored_energy_change_j_m3,
                    kelvin_heat_j_m3: channel(value.kelvin_heat_j_m3)?,
                    maxwell_heat_j_m3: channel(value.maxwell_heat_j_m3)?,
                    material_heat_j_m3: value.material_heat_j_m3,
                    stress_work_j_m3: value.stress_work_j_m3,
                    energy_balance_residual_j_m3: value.energy_balance_residual_j_m3,
                },
            })
        }).collect::<Result<Vec<_>, TensorShellHarmonicError>>()?;
        checkpoint(cx)?;
        Ok(TensorShellMaterialCausalStep {
            history: self.history,
            field: TensorShellMaterialCausalField {
                start_time_s_dd: field.start_time_s_dd,
                end_time_s_dd: field.end_time_s_dd,
                duration_s: field.duration_s,
                angular_frequencies_rad_s: field.angular_frequencies_rad_s,
                points,
                start_stored_energy_j: field.start_stored_energy_j,
                end_stored_energy_j: field.end_stored_energy_j,
                stored_energy_change_j: field.stored_energy_change_j,
                kelvin_heat_j: channel(field.kelvin_heat_j)?,
                maxwell_heat_j: channel(field.maxwell_heat_j)?,
                material_heat_j: field.material_heat_j,
                stress_work_j: field.stress_work_j,
                energy_balance_residual_j: field.energy_balance_residual_j,
            },
        })
    }
}
