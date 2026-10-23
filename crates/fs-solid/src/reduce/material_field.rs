//! Valve reduction over the same oriented/heterogeneous material field as the
//! three-dimensional solid operator. This is an ingress adapter, not a second
//! elastic model, modal solver, pressure projection or acoustic integrator.
//!
//! Resolve each material with `fs-material`, lower its tensor and orientation
//! through `TetElasticMaterial`, then assign it uniformly or per tetrahedron.
//! Cane axes, regional moisture states and synthetic laminates are therefore
//! data rather than an isotropic effective-E override. The caller must provide
//! actual admissible constitutive data: no transverse modulus, orientation or
//! dissipation is inferred from a material name. The scalar loss factor remains
//! an explicitly authored reduction input, not a sourced tensor loss law.
//!
//! The compatibility functions in the parent module retain their original
//! isotropic material identity. Here a field identity binds the constitutive
//! bits AND upstream identities in element order. A uniform field and identical
//! per-element expansion have the same identity. Reordering different materials
//! cannot reuse a valve card for a different assignment.

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;

use crate::linear3::{TetElasticError, TetElasticMaterial, TetMaterialField};
use super::{PressureLoadedValve, ReduceError, ValveCard, ValveCardRequest};

/// Run the ordinary certified valve reduction with a full solid material field.
///
/// The same geometry, support, mode-window, compliance-floor and orifice
/// admissions apply as on the isotropic path. No physics is selected by an
/// instrument or material name.
///
/// # Errors
/// Propagates the common assembler/reducer's material, geometry and compliance
/// refusals; in particular every tetrahedron must have exactly one material.
pub fn reduce_valve(
    request: &ValveCardRequest<'_, TetMaterialField<'_>>,
    cx: &Cx<'_>,
) -> Result<ValveCard, ReduceError> {
    let identity = field_identity(request.material, request.tetrahedra.len())?;
    super::reduce_valve_impl(request, request.material, identity, cx, None, None)
        .map(|reduced| reduced.card)
}

/// Reduce a material field and its work-conjugate distributed pressure port.
///
/// Pressure-force and swept-volume projections use the SAME anisotropic
/// modes, normalization, mesh and material assignment. A material swap cannot
/// leave a pressure area from the previous specimen attached to the new mass.
/// Signed higher-mode pressure areas remain signed.
///
/// # Errors
/// The common pressure-load and structural reduction refusals propagate.
pub fn reduce_pressure_loaded_valve(
    request: &ValveCardRequest<'_, TetMaterialField<'_>>,
    nodal_load_per_pa_m2: &[[f64; 3]],
    cx: &Cx<'_>,
) -> Result<PressureLoadedValve, ReduceError> {
    let identity = field_identity(request.material, request.tetrahedra.len())?;
    super::reduce_pressure_loaded_valve_impl(
        request, request.material, identity, nodal_load_per_pa_m2, cx, None,
    )
}

fn material_identity(material: &TetElasticMaterial) -> ContentHash {
    let mut hash = DomainHasher::new("org.frankensim.fs-solid.valve-constitutive-state.v1");
    hash.update(&material.material_state_identity);
    hash.update(&material.density_kg_m3.to_bits().to_le_bytes());
    for value in material.stiffness_mandel_pa().iter().flatten() {
        hash.update(&value.to_bits().to_le_bytes());
    }
    hash.finalize()
}

pub(super) fn field_identity(
    field: TetMaterialField<'_>,
    element_count: usize,
) -> Result<ContentHash, ReduceError> {
    let mut hash = DomainHasher::new("org.frankensim.fs-solid.valve-material-field.v1");
    hash.update(&(element_count as u64).to_le_bytes());
    match field {
        TetMaterialField::Uniform(material) => {
            // Hash the tensor once, not once per tetrahedron. The element
            // sequence below is identical to the explicit per-element form.
            let identity = material_identity(material);
            for _ in 0..element_count {
                hash.update(identity.as_bytes());
            }
        }
        TetMaterialField::PerElement(materials) => {
            if materials.len() != element_count {
                return Err(TetElasticError::MaterialCountMismatch {
                    expected: element_count,
                    got: materials.len(),
                }.into());
            }
            for material in materials {
                hash.update(material_identity(material).as_bytes());
            }
        }
    }
    // Do not duplicate or weaken SPD/density/identity admission. The common
    // TetLinearElasticProblem validates these values before minting a card.
    Ok(hash.finalize())
}

// Dependency inversion: the acoustic renderer consumes the scenario image,
// not fs-solid. In particular this does not add an fs-couple -> fs-solid cycle.
impl fs_scenario::reed_source::BeatingReedSource for super::PressureLoadedValve {
    type Admission = super::SingleModeValveAdmission;
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

