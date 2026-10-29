//! Material-derived loss in the existing three-dimensional valve reduction.
//!
//! A single storage pencil K'(reference_omega) supplies ALL modes and pressure
//! ports. Its symmetric viscous image is C = K''(reference_omega)/reference_omega.
//! Both matrices come from the ordinary tensor-Prony element owner, including
//! anisotropic, oriented, regional and Kelvin-Voigt material data.
//!
//! This is an explicitly admitted frozen-frequency image, not the hereditary
//! solid. The existing fs-modal damping owner checks whether the retained
//! matrix can use independent decay rates. At every retained modal frequency,
//! full projected storage and loss matrices are compared with the frozen image.
//! The reported row-norm bounds apply to the retained subspace at those samples;
//! they do not bound between-sample dispersion, omitted modes or nonlinear motion.
//! Static compliance is compared against the original EQUILIBRIUM stiffness,
//! never relabeled as the compliance of the frozen storage pencil.

use super::{
    PreparedValve, PressureLoadedValve, ReduceError, SingleModeValveAdmission, ValveCardRequest,
    ValveElementOrder,
};
use crate::linear3::{
    TetAssemblyBudget,
    visco::{
        ComplexStiffness, SolidMeshRef, TetViscoField, TetViscoMaterial, ViscoElements,
        ViscoStateBudget,
    },
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_modal::{
    CertifiedModalDamping, ModalDampingBlock, ModalDampingBudget, ModePair,
    certify_modal_damping,
};
use fs_sparse::Csr;

/// Complete source-bound modal images, retaining nonproportional loss.
pub mod modal;
/// Equilibrium-mode images with the original tensor-Prony internal histories.
pub mod hereditary;
/// Original exterior geometry for complete source-owned surface interactions.
pub mod surface;
pub use surface::MaterialValveSurface;
pub use hereditary::{HereditaryModalValve, HereditaryValveControls, reduce_hereditary_pressure_loaded_valve};
pub use modal::{MaterialModalValve, ModalStaticCompliance, reduce_modal_pressure_loaded_valve};

/// Explicit finite-dimensional accuracy and allocation policies. No default
/// chooses a material reference frequency or promotes a decay model to truth.
#[derive(Clone, Copy, Debug)]
pub struct Controls {
    /// Common positive frequency at which BOTH storage and viscosity are frozen.
    pub reference_omega_rad_s: f64,
    /// Dense projection cap, in 1..=64. The existing sparse modal solve is reused.
    pub maximum_modes: usize,
    pub assembly: TetAssemblyBudget,
    pub internal_variables: ViscoStateBudget,
    /// Maximum relative storage quadratic-form error at each modal sample.
    /// Finite in [0, 1); the reference metric is the frozen storage energy.
    pub maximum_storage_sample_error: f64,
    /// Maximum absolute loss-factor operator error at each modal sample.
    /// This bounds omega*C_frozen-K''(omega), normalized by frozen storage.
    /// Finite in [0, 1); it is NOT a relative error divided by a vanishing loss.
    pub maximum_loss_sample_error: f64,
    /// Largest admitted frozen-image modal loss factor 2*zeta, in (0, 1).
    pub maximum_modal_loss_factor: f64,
    /// Existing retained-matrix pole accuracy policy. The legacy MaterialValve
    /// refuses a required full block; the complete material image retains it.
    pub damping: ModalDampingBudget,
}

/// Report for ALL modes selected by the requested window, including a mode
/// whose mean opening vanishes and is therefore absent from the valve card.
#[derive(Clone, Debug, PartialEq)]
pub struct Approximation {
    pub reference_omega_rad_s: f64,
    pub modal_omega_rad_s: Vec<f64>,
    /// Diagonal entries of the projected mass-normalized C [1/s].
    pub modal_damping_per_s: Vec<f64>,
    /// Complete symmetric mass-normalized viscosity C [1/s], row-major in
    /// modal_omega_rad_s order. Off-diagonal entries are never replaced by fits.
    pub modal_viscosity_per_s: Vec<f64>,
    /// Largest symmetric row-norm bound at the retained modal frequencies.
    pub maximum_storage_sample_error: f64,
    pub maximum_loss_sample_error: f64,
    /// Certificate of the RETAINED constant-coefficient damping matrix only.
    pub damping: CertifiedModalDamping,
}

/// The original material field remains attached to its admitted playing image.
/// No history is reset or transplanted: this object is an offline reduction.
pub struct MaterialValve<'a> {
    valve: PressureLoadedValve,
    materials: TetViscoField<'a>,
    approximation: Approximation,
}
impl<'a> MaterialValve<'a> {
    #[must_use]
    pub const fn valve(&self) -> &PressureLoadedValve {
        &self.valve
    }
    #[must_use]
    pub const fn materials(&self) -> TetViscoField<'a> {
        self.materials
    }
    #[must_use]
    pub const fn approximation(&self) -> &Approximation {
        &self.approximation
    }

    /// Lower the actual material-derived first mode through the existing
    /// pressure/compliance admission and BeatingReed runtime image.
    pub fn beating_reed_checked(
        &self,
        pressure_pa: f64,
        attack_s: f64,
        admission: SingleModeValveAdmission,
    ) -> Result<fs_scenario::BeatingReed, ReduceError> {
        self.valve
            .beating_reed_checked(pressure_pa, attack_s, admission)
    }
}

impl fs_scenario::reed_source::BeatingReedSource for MaterialValve<'_> {
    type Admission = SingleModeValveAdmission;
    type Error = ReduceError;
    fn compile_beating_reed(
        &self,
        blowing_pressure_pa: f64,
        attack_s: f64,
        admission: &Self::Admission,
    ) -> Result<fs_scenario::BeatingReed, Self::Error> {
        self.beating_reed_checked(blowing_pressure_pa, attack_s, *admission)
    }
}

/// Reduce a source tensor-Prony/Kelvin-Voigt field with a work-conjugate
/// distributed pressure port. Authored loss must be zero to avoid double count.
/// Every source band must admit equilibrium (omega=0), the reference frequency,
/// and every retained modal sample. Geometry/order/supports/projections and the
/// exact static-pressure comparisons are the SAME owners as elastic reduction.
pub fn reduce_pressure_loaded_valve<'a>(
    request: &ValveCardRequest<'_, TetViscoField<'a>>,
    nodal_load_per_pa_m2: &[[f64; 3]],
    controls: &Controls,
    cx: &Cx<'_>,
) -> Result<MaterialValve<'a>, ReduceError> {
    let (_, mut reduced) = reduce_material(
        request, nodal_load_per_pa_m2, controls, false, cx,
    )?;
    let approximation = reduced
        .material_loss
        .take()
        .ok_or_else(|| invalid("material valve reduction did not retain its damping admission"))?;
    let valve = super::pressure_loaded_from_reduction(reduced, nodal_load_per_pa_m2)?;
    Ok(MaterialValve {
        valve,
        materials: request.material,
        approximation,
    })
}

fn reduce_material(
    request: &ValveCardRequest<'_, TetViscoField<'_>>,
    nodal_load_per_pa_m2: &[[f64; 3]],
    controls: &Controls,
    complete_basis: bool,
    cx: &Cx<'_>,
) -> Result<(ViscoElements, super::ReducedValve), ReduceError> {
    controls.validate()?;
    if request.loss_factor != 0.0 {
        return Err(invalid(
            "source viscoelastic valve requires zero authored loss",
        ));
    }
    if nodal_load_per_pa_m2.len() != request.nodes_m.len()
        || nodal_load_per_pa_m2
            .iter()
            .flatten()
            .any(|v| !v.is_finite())
    {
        return Err(invalid(
            "pressure load must contain one finite vector per mesh node",
        ));
    }
    super::validate_element_topology(request)?;
    let mesh = match request.element_order {
        ValveElementOrder::Tet4 => SolidMeshRef::Tet4 {
            nodes_m: request.nodes_m,
            elements: request.tetrahedra,
        },
        ValveElementOrder::Tet10 {
            elements,
            quadrature,
        } => SolidMeshRef::Tet10 {
            nodes_m: request.nodes_m,
            elements,
            quadrature,
        },
    };
    let system = ViscoElements::new(
        mesh,
        request.material,
        request.fixed_dofs,
        controls.assembly,
        controls.internal_variables,
        cx,
    )?;
    // This query admits the DC comparator against each original source band.
    // Its real matrix is the existing equilibrium assembly, not a new model.
    let equilibrium = system.complex_stiffness(0.0, cx)?;
    let reference = system.complex_stiffness(controls.reference_omega_rad_s, cx)?;
    let mut assembly = system.equilibrium().clone();
    assembly.stiffness = reference.real.clone();
    let projection = Projection {
        system: &system,
        reference: &reference,
        controls,
        complete_basis,
    };
    let reduced = super::reduce_valve_impl(
        request,
        None,
        field_identity(
            request.material,
            request.tetrahedra.len(),
            controls.reference_omega_rad_s,
            cx,
        )?,
        cx,
        Some(nodal_load_per_pa_m2),
        Some(PreparedValve {
            assembly: &assembly,
            observation_nodes_m: request.nodes_m,
            static_stiffness: Some(&equilibrium.real),
            material_loss: Some(&projection),
            complete_basis,
        }),
    )?;
    Ok((system, reduced))
}

fn invalid(what: &'static str) -> ReduceError {
    ReduceError::Invalid { what }
}
fn poll(cx: &Cx<'_>) -> Result<(), ReduceError> {
    cx.checkpoint().map_err(|_| {
        crate::linear3::visco::ViscoElementError::Elastic(
            crate::linear3::TetElasticError::Cancelled,
        )
        .into()
    })
}
impl Controls {
    fn validate(&self) -> Result<(), ReduceError> {
        if !self.reference_omega_rad_s.is_finite()
            || self.reference_omega_rad_s <= 0.0
            || self.maximum_modes == 0
            || self.maximum_modes > 64
            || [
                self.maximum_storage_sample_error,
                self.maximum_loss_sample_error,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0 || *x >= 1.0)
            || !self.maximum_modal_loss_factor.is_finite()
            || self.maximum_modal_loss_factor <= 0.0
            || self.maximum_modal_loss_factor >= 1.0
        {
            return Err(invalid(
                "invalid frozen material-valve frequency or approximation policy",
            ));
        }
        Ok(())
    }
}

pub(super) struct Projection<'a> {
    system: &'a ViscoElements,
    reference: &'a ComplexStiffness,
    controls: &'a Controls,
    pub(super) complete_basis: bool,
}
impl Projection<'_> {
    pub(super) fn project(
        &self,
        modes: &[ModePair],
        cx: &Cx<'_>,
    ) -> Result<Approximation, ReduceError> {
        let n = modes.len();
        if n == 0 || n > self.controls.maximum_modes {
            return Err(invalid(
                "material valve exceeds its admitted modal projection size",
            ));
        }
        let omega: Vec<f64> = modes.iter().map(|m| m.lambda.sqrt()).collect();
        let mut damping = project_matrix(&self.reference.imag, modes, cx)?;
        for value in &mut damping {
            *value /= self.reference.omega_rad_s;
        }
        let diag: Vec<f64> = (0..n).map(|i| damping[i * n + i]).collect();
        if diag.iter().zip(&omega).any(|(&c, &w)| {
            !c.is_finite()
                || c < 0.0
                || !w.is_finite()
                || w <= 0.0
                || c / w > self.controls.maximum_modal_loss_factor
        }) {
            return Err(invalid(
                "material valve exceeds its admitted modal loss factor",
            ));
        }
        // This is the existing coupled-pole certificate, not an independently
        // invented small-zeta or eigenvalue-gap heuristic.
        let certificate = certify_modal_damping(&omega, &damping, &self.controls.damping)?;
        if !self.complete_basis && certificate
            .blocks
            .iter()
            .any(|b| matches!(b, ModalDampingBlock::Full { .. }))
        {
            return Err(invalid(
                "material valve needs coupled damping blocks absent from its runtime",
            ));
        }
        let mut storage_error: f64 = 0.0;
        let mut loss_error: f64 = 0.0;
        for &sample in &omega {
            poll(cx)?;
            let actual = self.system.complex_stiffness(sample, cx)?;
            let real = project_matrix(&actual.real, modes, cx)?;
            let imag = project_matrix(&actual.imag, modes, cx)?;
            for i in 0..n {
                let mut storage_row = 0.0;
                let mut loss_row = 0.0;
                for j in 0..n {
                    let metric = omega[i] * omega[j];
                    let frozen = if i == j { modes[i].lambda } else { 0.0 };
                    storage_row += ((real[i * n + j] - frozen) / metric).abs();
                    loss_row += ((imag[i * n + j] - sample * damping[i * n + j]) / metric).abs();
                }
                if !storage_row.is_finite() || !loss_row.is_finite() {
                    return Err(invalid(
                        "unrepresentable material-valve approximation comparison",
                    ));
                }
                storage_error = storage_error.max(storage_row);
                loss_error = loss_error.max(loss_row);
            }
        }
        if storage_error > self.controls.maximum_storage_sample_error {
            return Err(invalid(
                "material storage dispersion exceeds the admitted frozen valve image",
            ));
        }
        if loss_error > self.controls.maximum_loss_sample_error {
            return Err(invalid(
                "material loss dispersion exceeds the admitted frozen valve image",
            ));
        }
        Ok(Approximation {
            reference_omega_rad_s: self.reference.omega_rad_s,
            modal_omega_rad_s: omega,
            modal_damping_per_s: diag,
            modal_viscosity_per_s: damping,
            maximum_storage_sample_error: storage_error,
            maximum_loss_sample_error: loss_error,
            damping: certificate,
        })
    }
}

// Form the symmetric Galerkin projection with one work-conjugate application
// per column. Never average material tensors or mix shapes from different K.
fn project_matrix(matrix: &Csr, modes: &[ModePair], cx: &Cx<'_>) -> Result<Vec<f64>, ReduceError> {
    let n = modes.len();
    let mut projected = vec![0.0; n * n];
    let mut product = vec![0.0; matrix.nrows()];
    for j in 0..n {
        poll(cx)?;
        matrix.spmv(&modes[j].phi, &mut product);
        for i in 0..=j {
            let value = modes[i]
                .phi
                .iter()
                .zip(&product)
                .map(|(a, b)| a * b)
                .sum::<f64>();
            if !value.is_finite() {
                return Err(invalid("nonfinite material modal projection"));
            }
            projected[i * n + j] = value;
            projected[j * n + i] = value;
        }
    }
    Ok(projected)
}

fn material_identity(material: &TetViscoMaterial) -> ContentHash {
    let mut h = DomainHasher::new("org.frankensim.fs-solid.valve-visco-material.v1");
    h.update(&material.material_state_identity);
    h.update(&material.density_kg_m3.to_bits().to_le_bytes());
    for x in material.law.c_eq().iter().flatten() {
        h.update(&x.to_bits().to_le_bytes());
    }
    match material.law.kv_pa_s() {
        None => h.update(&[0]),
        Some(kv) => {
            h.update(&[1]);
            for x in kv.iter().flatten() {
                h.update(&x.to_bits().to_le_bytes());
            }
        }
    }
    h.update(&(material.law.branches().len() as u64).to_le_bytes());
    for branch in material.law.branches() {
        h.update(&branch.tau_s.to_bits().to_le_bytes());
        for x in branch.c_pa.iter().flatten() {
            h.update(&x.to_bits().to_le_bytes());
        }
    }
    let (low, high) = material.law.band_rad_s();
    h.update(&low.to_bits().to_le_bytes());
    h.update(&high.to_bits().to_le_bytes());
    h.finalize()
}
fn field_identity(
    field: TetViscoField<'_>,
    count: usize,
    reference_omega_rad_s: f64,
    cx: &Cx<'_>,
) -> Result<ContentHash, ReduceError> {
    let mut h = DomainHasher::new("org.frankensim.fs-solid.valve-visco-field.v1");
    h.update(&(count as u64).to_le_bytes());
    h.update(&reference_omega_rad_s.to_bits().to_le_bytes());
    match field {
        TetViscoField::Uniform(material) => {
            let identity = material_identity(material);
            for _ in 0..count {
                poll(cx)?;
                h.update(identity.as_bytes());
            }
        }
        TetViscoField::PerElement(materials) => {
            if materials.len() != count {
                return Err(invalid("visco material count mismatch"));
            }
            for material in materials {
                poll(cx)?;
                h.update(material_identity(material).as_bytes());
            }
        }
    }
    Ok(h.finalize())
}
