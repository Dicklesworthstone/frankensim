//! Prescribed thermal-strain history on the original bonded laminate.
//!
//! Each Prony layer must share one plane-stress strain extension across its
//! equilibrium, every Maxwell tensor and optional Kelvin viscosity. This
//! includes fixed-Poisson isotropic spectra. Different bulk/shear relaxation
//! spectra generally refuse: condensing their branches independently is wrong.
//! Original full-3D histories integrate affine strain ramps exactly, including
//! work and dissipation. Two points integrate affine thickness fields exactly.
//!
//! Material tensors and relaxation times remain fixed. Profiles are already
//! integrated free strains, affine in time between accepted endpoints. This
//! is not a varying-temperature Prony thermodynamic potential or a heat solver.

use super::{LaminaThermalStrain, LayerThermalStress, dot, matvec, rotate_strain};
use crate::Tangent6;
use crate::laminate::{LaminateSection, LayerLaw};
use crate::visco_tensor::{TensorProny, TensorPronyError, TensorPronyHistory};
use fs_la::factor::cholesky;
use std::sync::Arc;

const GAUSS: [f64; 2] = [0.211_324_865_405_187_13, 0.788_675_134_594_812_9];
const IN: [usize; 3] = [0, 1, 3];
const OUT: [usize; 3] = [2, 4, 5];
type Extension = [[f64; 3]; 6];

/// Failure leaves the entire accepted laminate history unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum LaminateHistoryError {
    Invalid(&'static str),
    /// The original tensors do not share a time-independent plane-stress map.
    IncompatiblePlaneStress {
        layer: usize,
    },
    /// Original native material refusal, retaining its layer index.
    Material {
        layer: usize,
        source: TensorPronyError,
    },
    Cancelled,
}
impl core::fmt::Display for LaminateHistoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MATERIAL-LAMINATE-HISTORY: {self:?}")
    }
}
impl std::error::Error for LaminateHistoryError {}
type Result<T> = core::result::Result<T, LaminateHistoryError>;

/// Resultants and independently integrated work per original midsurface area.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LaminateThermalStep {
    pub membrane_force_n_m: [f64; 3],
    pub bending_moment_n: [f64; 3],
    /// Equilibrium plus all retained Maxwell spring storage [J/m²].
    pub stored_energy_j_m2: f64,
    /// Recoverable Maxwell spring storage alone [J/m²].
    pub branch_stored_energy_j_m2: f64,
    /// Integral of stress against the actual total-strain increment [J/m²].
    pub mechanical_work_j_m2: f64,
    /// Minus the integral of stress against the free-strain increment [J/m²].
    /// This is an eigenstrain work port, not additional dissipative heat.
    pub thermal_strain_work_j_m2: f64,
    /// Independently integrated nonnegative Maxwell/Kelvin loss [J/m²].
    pub dissipated_j_m2: f64,
    /// Both work ports minus storage change and dissipation; no balancing term.
    pub balance_residual_j_m2: f64,
}

#[derive(Clone, Debug, PartialEq)]
enum Memory {
    Elastic([f64; 9]),
    Prony {
        extension: Extension,
        points: [TensorPronyHistory; 2],
    },
}
#[derive(Clone, Debug, PartialEq)]
struct LayerHistory {
    memory: Memory,
    rotation: [[f64; 3]; 3],
    z: [f64; 2],
    weight: f64,
    total: [[f64; 3]; 2],
    free: [[f64; 3]; 2],
    stress: [[f64; 3]; 2],
}

/// Virgin-origin, source-preserving section history. Every layer retains its
/// original spatial position, density and law; no reference-frequency stiffness
/// is used for a Prony step. Preparation and stepping are linear in branch count.
#[derive(Clone, Debug, PartialEq)]
pub struct LaminateThermalHistory {
    section: Arc<LaminateSection>,
    layers: Vec<LayerHistory>,
    time_s: f64,
    stored: f64,
    dissipated: f64,
}

impl LaminateThermalHistory {
    /// Start relaxed with zero total and free strain, using the original laws.
    /// A later zero-duration Maxwell jump is admitted by the native owner;
    /// nonzero Kelvin viscosity requires positive duration for strain changes.
    pub fn try_new(
        section: Arc<LaminateSection>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let mut layers = Vec::with_capacity(section.layers().len());
        for (i, (layer, geometry)) in section
            .layers()
            .iter()
            .zip(section.layer_integrals())
            .enumerate()
        {
            poll(&mut cancelled)?;
            let memory = match &layer.law {
                LayerLaw::Prony(law) => Memory::Prony {
                    extension: common_extension(law, i, &mut cancelled)?,
                    points: [law.history(), law.history()],
                },
                _ => Memory::Elastic(
                    crate::laminate::principal_stiffness(&layer.law, i)
                        .map_err(|_| {
                            LaminateHistoryError::Invalid("original elastic layer refused")
                        })?
                        .0,
                ),
            };
            // Local engineering strain = T * section engineering strain;
            // section physical stress = T^T * local physical stress.
            let columns = core::array::from_fn::<_, 3, _>(|j| {
                let mut unit = [0.0; 3];
                unit[j] = 1.0;
                rotate_strain(unit, -layer.orientation_rad)
            });
            let weight = 0.5 * (geometry.z_top_m - geometry.z_bottom_m);
            if !weight.is_finite() || weight <= 0.0 {
                return Err(LaminateHistoryError::Invalid(
                    "unrepresentable layer quadrature weight",
                ));
            }
            layers.push(LayerHistory {
                memory,
                rotation: core::array::from_fn(|i| core::array::from_fn(|j| columns[j][i])),
                z: GAUSS.map(|u| (1.0 - u) * geometry.z_bottom_m + u * geometry.z_top_m),
                weight,
                total: [[0.0; 3]; 2],
                free: [[0.0; 3]; 2],
                stress: [[0.0; 3]; 2],
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            section,
            layers,
            time_s: 0.0,
            stored: 0.0,
            dissipated: 0.0,
        })
    }

    pub fn section(&self) -> &LaminateSection {
        &self.section
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn stored_energy_j_m2(&self) -> f64 {
        self.stored
    }
    pub const fn dissipated_j_m2(&self) -> f64 {
        self.dissipated
    }

    /// Advance prescribed `epsilon(z)=membrane+z*gradient_per_m` and one
    /// principal-frame free-strain profile per original layer. All endpoint
    /// fields interpolate linearly in time; branch memory is never reset.
    /// Complete proposal cloning makes failures/cancellation transactional.
    pub fn advance(
        &mut self,
        membrane: [f64; 3],
        gradient_per_m: [f64; 3],
        profiles: &[LaminaThermalStrain],
        dt_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateThermalStep> {
        poll(&mut cancelled)?;
        let time = self.time_s + dt_s;
        if profiles.len() != self.layers.len()
            || !dt_s.is_finite()
            || dt_s < 0.0
            || !time.is_finite()
            || (dt_s > 0.0 && time <= self.time_s)
            || membrane
                .iter()
                .chain(&gradient_per_m)
                .any(|v| !v.is_finite())
        {
            return Err(LaminateHistoryError::Invalid(
                "complete finite laminate strain profiles and representable time required",
            ));
        }
        let mut next = self.clone();
        let mut report = LaminateThermalStep::default();
        for (i, (layer, profile)) in next.layers.iter_mut().zip(profiles).enumerate() {
            poll(&mut cancelled)?;
            let faces = profile.faces();
            let local_total = |z| {
                apply3(
                    &layer.rotation,
                    core::array::from_fn(|j| membrane[j] + z * gradient_per_m[j]),
                )
            };
            if let (Memory::Prony { extension, .. }, LayerLaw::Prony(law)) =
                (&layer.memory, &self.section.layers()[i].law)
            {
                let geometry = &self.section.layer_integrals()[i];
                // Norm convexity bounds every spatial/time point by its face
                // endpoints. Gauss points alone do not admit the layer faces.
                for (z, free) in [geometry.z_bottom_m, geometry.z_top_m]
                    .into_iter()
                    .zip(faces)
                {
                    let total = local_total(z);
                    let in_plane = [
                        total[0],
                        total[1],
                        0.0,
                        total[2] * core::f64::consts::FRAC_1_SQRT_2,
                        0.0,
                        0.0,
                    ];
                    law.validate_mandel_strain(&in_plane)
                        .and_then(|_| {
                            law.validate_mandel_strain(&extend(
                                extension,
                                core::array::from_fn(|j| total[j] - free[j]),
                            ))
                        })
                        .map_err(|_| {
                            LaminateHistoryError::Invalid(
                                "layer face exceeds original strain domain",
                            )
                        })?;
                }
            }
            for p in 0..2 {
                let total = local_total(layer.z[p]);
                let free = core::array::from_fn(|j| {
                    (1.0 - GAUSS[p]) * faces[0][j] + GAUSS[p] * faces[1][j]
                });
                let elastic = core::array::from_fn(|j| total[j] - free[j]);
                let old: [f64; 3] = core::array::from_fn(|j| layer.total[p][j] - layer.free[p][j]);
                let (stress, mean, stored, branch, heat) = match &mut layer.memory {
                    Memory::Elastic(q) => {
                        let stress = matvec(q, elastic);
                        let mean = matvec(q, core::array::from_fn(|j| 0.5 * (old[j] + elastic[j])));
                        (stress, mean, 0.5 * dot(elastic, stress), 0.0, 0.0)
                    }
                    Memory::Prony { extension, points } => {
                        let mandel = extend(extension, elastic);
                        let tensor = core::array::from_fn(|j| {
                            if j < 3 {
                                mandel[j]
                            } else {
                                mandel[j] * core::f64::consts::FRAC_1_SQRT_2
                            }
                        });
                        let step = points[p]
                            .advance(tensor, dt_s, &mut cancelled)
                            .map_err(|source| material_error(i, source))?;
                        (
                            IN.map(|j| step.stress[j]),
                            IN.map(|j| step.mean_stress[j]),
                            step.stored_energy,
                            points[p].branch_stored_energy(),
                            step.dissipated,
                        )
                    }
                };
                let global = core::array::from_fn(|j| {
                    (0..3).map(|k| layer.rotation[k][j] * stress[k]).sum()
                });
                for j in 0..3 {
                    report.membrane_force_n_m[j] += layer.weight * global[j];
                    report.bending_moment_n[j] += layer.weight * layer.z[p] * global[j];
                }
                report.stored_energy_j_m2 += layer.weight * stored;
                report.branch_stored_energy_j_m2 += layer.weight * branch;
                report.mechanical_work_j_m2 += layer.weight
                    * dot(mean, core::array::from_fn(|j| total[j] - layer.total[p][j]));
                report.thermal_strain_work_j_m2 -=
                    layer.weight * dot(mean, core::array::from_fn(|j| free[j] - layer.free[p][j]));
                report.dissipated_j_m2 += layer.weight * heat;
                layer.total[p] = total;
                layer.free[p] = free;
                layer.stress[p] = global;
            }
        }
        report.balance_residual_j_m2 = report.mechanical_work_j_m2
            + report.thermal_strain_work_j_m2
            - (report.stored_energy_j_m2 - self.stored)
            - report.dissipated_j_m2;
        next.stored = report.stored_energy_j_m2;
        next.dissipated += report.dissipated_j_m2;
        next.time_s = time;
        if report
            .membrane_force_n_m
            .iter()
            .chain(&report.bending_moment_n)
            .chain([
                &report.stored_energy_j_m2,
                &report.branch_stored_energy_j_m2,
                &report.mechanical_work_j_m2,
                &report.thermal_strain_work_j_m2,
                &report.dissipated_j_m2,
                &report.balance_residual_j_m2,
                &next.dissipated,
            ])
            .any(|v| !v.is_finite())
            || report.stored_energy_j_m2 < 0.0
        {
            return Err(LaminateHistoryError::Invalid(
                "unrepresentable laminate history response",
            ));
        }
        poll(&mut cancelled)?;
        *self = next;
        Ok(report)
    }

    /// Recover each original layer's physical face stress from its exact
    /// affine stress field. Interfaces may retain incompatible residual stress.
    pub fn layer_stress(&self, layer: usize) -> Result<LayerThermalStress> {
        let layer = self
            .layers
            .get(layer)
            .ok_or(LaminateHistoryError::Invalid("unknown recovery layer"))?;
        let slope: [f64; 3] = core::array::from_fn(|j| {
            (layer.stress[1][j] - layer.stress[0][j]) / (GAUSS[1] - GAUSS[0])
        });
        let stress = LayerThermalStress {
            bottom_pa: core::array::from_fn(|j| layer.stress[0][j] - GAUSS[0] * slope[j]),
            top_pa: core::array::from_fn(|j| layer.stress[1][j] + (1.0 - GAUSS[1]) * slope[j]),
        };
        if stress
            .bottom_pa
            .iter()
            .chain(&stress.top_pa)
            .any(|v| !v.is_finite())
        {
            return Err(LaminateHistoryError::Invalid(
                "unrepresentable recovered layer stress",
            ));
        }
        Ok(stress)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(LaminateHistoryError::Cancelled)
    } else {
        Ok(())
    }
}
fn material_error(layer: usize, source: TensorPronyError) -> LaminateHistoryError {
    if source == TensorPronyError::Cancelled {
        LaminateHistoryError::Cancelled
    } else {
        LaminateHistoryError::Material { layer, source }
    }
}
fn apply3(matrix: &[[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| dot(matrix[i], vector))
}
fn extend(extension: &Extension, strain: [f64; 3]) -> [f64; 6] {
    core::array::from_fn(|i| dot(extension[i], strain))
}

fn common_extension(
    law: &TensorProny,
    layer: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Extension> {
    let refused = || LaminateHistoryError::IncompatiblePlaneStress { layer };
    let c = law.equilibrium_mandel();
    let scale = c.iter().flatten().fold(0.0_f64, |m, v| m.max(v.abs()));
    if !scale.is_finite() || scale <= 0.0 {
        return Err(refused());
    }
    let normal: Vec<_> = OUT
        .iter()
        .flat_map(|&i| OUT.map(|j| c[i][j] / scale))
        .collect();
    let factor = cholesky(&normal, 3).map_err(|_| refused())?;
    let mut extension = [[0.0; 3]; 6];
    for (j, &axis) in IN.iter().enumerate() {
        let weight = if j == 2 {
            core::f64::consts::FRAC_1_SQRT_2
        } else {
            1.0
        };
        extension[axis][j] = weight;
        let mut rhs = OUT.map(|i| -(c[i][axis] / scale) * weight);
        factor.solve(&mut rhs);
        for (k, &axis) in OUT.iter().enumerate() {
            extension[axis][j] = rhs[k];
        }
    }
    if extension.iter().flatten().any(|v| !v.is_finite()) {
        return Err(refused());
    }
    // Each original tensor must annihilate the SAME out-of-plane stress map.
    // This is an admission test, never a branchwise Schur complement or refit.
    for tensor in core::iter::once(c)
        .chain(law.branches().iter().map(|b| &b.c_pa))
        .chain(law.kv_pa_s())
    {
        poll(cancelled)?;
        let scale = tensor.iter().flatten().fold(0.0_f64, |m, v| m.max(v.abs()));
        if scale == 0.0 {
            continue;
        }
        if !common_map(tensor, &extension, scale) {
            return Err(refused());
        }
    }
    Ok(extension)
}
fn common_map(c: &Tangent6, extension: &Extension, scale: f64) -> bool {
    for i in OUT {
        for j in 0..3 {
            let mut residual = 0.0;
            let mut magnitude = 0.0;
            for k in 0..6 {
                let term = (c[i][k] / scale) * extension[k][j];
                residual += term;
                magnitude += term.abs();
            }
            if !residual.is_finite()
                || !magnitude.is_finite()
                || residual.abs() > 512.0 * f64::EPSILON * magnitude
            {
                return false;
            }
        }
    }
    true
}
