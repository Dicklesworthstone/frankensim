//! Native electrostatic mechanical loads and resolved guided motion.
//!
//! The electric Lorentz term is `rho_free E` [N/m3]. Charges must be bound to
//! the selected material, or transfer their momentum to it on the timescale
//! being resolved. See Haus and Melcher, *Electromagnetic Fields and Energy*,
//! [11.9, Macroscopic Force Densities](https://web.mit.edu/6.013_book/www/chapter11/11.9.html).
//! This is only the force on free VOLUME charge: electrode and internal-sheet
//! charge, polarization traction, electrostriction and dielectric-interface
//! forces are separate physical contributions, never inferred from this term.
//!
//! The free-volume-charge adapter transfers one accepted field into the native
//! solid load integral. It is not a charge-transport or moving-mesh integrator;
//! it supplies neither electromagnetic feedback from the resulting deformation
//! nor a closed electric/mechanical energy trajectory. Do not combine it with
//! an enclosing Maxwell-stress load that already includes these same charges.

/// Resolved moving capacitor geometry and reciprocal guided-body work.
pub mod motion;

use fs_em::{
    electrostatic::nonlinear::spatial::SpatialNonlinearElectricField,
    spatial::electrostatic::SpatialElectricField,
};
use fs_exec::Cx;
use fs_solid::loads3::{LoadError3, TetLoadCase3, TetLoadGeometry3, VolumeLoad3};

fn poll(cx: &Cx<'_>) -> Result<(), LoadError3> {
    cx.checkpoint().map_err(|_| LoadError3::Cancelled)
}

/// Integrate the solved `rho_free E` body force over explicitly selected cells.
///
/// `charged_cells` is a strictly increasing list of original tetrahedron IDs;
/// unselected cells receive zero force, even if they carry electric charge in
/// the field problem. An empty selection gives a zero load. Each selected
/// charge belongs mechanically to this solid; choosing a subset does not
/// change or re-solve the original electric field.
///
/// Positions, connectivity, vertex order and element order must exactly match
/// the mechanical load geometry. No nearest-cell guess, coordinate scaling,
/// material-interface averaging or remeshing occurs. The source field retains
/// its complete material laws and all contributions to E, while ONLY its
/// original free volume charge enters the force density.
///
/// The existing P1 load owner integrates `V rho E/4` per tetrahedron node,
/// preserving resultant, first moment and virtual work for the represented
/// cell-constant force and affine displacement. Forces at fixed DOFs remain in
/// the returned native load, so support reactions include them. Its `solve`
/// method uses the existing solid displacement/stress solver on the same mesh.
/// Bounds follow the already-admitted source and load meshes; work and storage
/// are linear in their size. Cancellation or failure publishes no partial load.
pub fn spatial_coulomb_volume_load<'g, 'm>(
    field: &SpatialElectricField<'_>,
    geometry: &'g TetLoadGeometry3<'m>,
    charged_cells: &[usize],
    cx: &Cx<'_>,
) -> Result<TetLoadCase3<'g, 'm>, LoadError3> {
    let source = field.problem();
    assemble_coulomb_volume_load(
        SpatialCoulombData {
            positions_m: source.positions_m(),
            tetrahedra: &source.mesh().tets,
            free_charge_density_c_m3: field.free_charge_density_c_m3(),
            electric_field_v_m: field.electric_field_v_m(),
        },
        geometry,
        charged_cells,
        cx,
    )
}

/// Integrate `rho_free E` from an accepted native nonlinear tetrahedral field.
///
/// This uses the actual field at its accepted material operating point and the
/// original free volume charge. Constitutive polarization, electrode charge and
/// internal-interface forces are not reclassified as volume charge. The exact
/// mesh/frame checks, cell selection and native solid P1 integration are the
/// same as [`spatial_coulomb_volume_load`], including forces on fixed nodes.
///
/// The result is an instantaneous load on mechanically bound free charge. It
/// does not advance the electric state or its geometry, infer charge transport,
/// or close a deformation-dependent energy trajectory. An enclosing Maxwell
/// load that already includes these charges must not be added a second time.
pub fn spatial_nonlinear_coulomb_volume_load<'g, 'm>(
    field: &SpatialNonlinearElectricField<'_>,
    geometry: &'g TetLoadGeometry3<'m>,
    charged_cells: &[usize],
    cx: &Cx<'_>,
) -> Result<TetLoadCase3<'g, 'm>, LoadError3> {
    let source = field.problem();
    assemble_coulomb_volume_load(
        SpatialCoulombData {
            positions_m: source.positions_m(),
            tetrahedra: &source.mesh().tets,
            free_charge_density_c_m3: field.free_charge_density_c_m3(),
            electric_field_v_m: field.electric_field_v_m(),
        },
        geometry,
        charged_cells,
        cx,
    )
}

// Only accepted native field owners construct this view. The common adapter
// cannot accept independent caller-authored E/rho arrays or a substitute mesh.
struct SpatialCoulombData<'a> {
    positions_m: &'a [[f64; 3]],
    tetrahedra: &'a [[u32; 4]],
    free_charge_density_c_m3: &'a [f64],
    electric_field_v_m: &'a [[f64; 3]],
}

fn assemble_coulomb_volume_load<'g, 'm>(
    source: SpatialCoulombData<'_>,
    geometry: &'g TetLoadGeometry3<'m>,
    charged_cells: &[usize],
    cx: &Cx<'_>,
) -> Result<TetLoadCase3<'g, 'm>, LoadError3> {
    poll(cx)?;
    let source_tets = source.tetrahedra;
    let nodes = geometry.nodes_m();
    let tets = geometry.tetrahedra();
    if source.positions_m.len() != nodes.len() || source_tets.len() != tets.len() {
        return Err(LoadError3::Invalid(
            "electrostatic body loads require the original reference volume mesh",
        ));
    }
    for (electric, mechanical) in source.positions_m.iter().zip(nodes) {
        poll(cx)?;
        if electric != mechanical {
            return Err(LoadError3::Invalid(
                "electrostatic body-load positions or frame differ from the solid",
            ));
        }
    }
    for (electric, mechanical) in source_tets.iter().zip(tets) {
        poll(cx)?;
        if electric.map(|v| v as usize) != *mechanical {
            return Err(LoadError3::Invalid(
                "electrostatic body-load connectivity or cell order differs from the solid",
            ));
        }
    }
    if charged_cells.len() > tets.len() {
        return Err(LoadError3::Budget);
    }
    let mut previous = None;
    let mut force_density_n_m3 = vec![[0.0; 3]; tets.len()];
    for &cell in charged_cells {
        poll(cx)?;
        if cell >= tets.len() || previous.is_some_and(|p| cell <= p) {
            return Err(LoadError3::Invalid(
                "charged body cells must be sorted, unique original tetrahedron IDs",
            ));
        }
        previous = Some(cell);
        let rho = source.free_charge_density_c_m3[cell];
        let electric = source.electric_field_v_m[cell];
        let density = electric.map(|component| rho * component);
        let nodal_volume = geometry.element_volumes_m3()[cell] / 4.0;
        for axis in 0..3 {
            let nodal_force = nodal_volume * density[axis];
            if !density[axis].is_finite()
                || !nodal_force.is_finite()
                || (rho != 0.0 && electric[axis] != 0.0 && density[axis] == 0.0)
                || (density[axis] != 0.0 && nodal_force == 0.0)
            {
                return Err(LoadError3::Invalid(
                    "electrostatic body force is not representable in SI load units",
                ));
            }
        }
        force_density_n_m3[cell] = density;
    }
    poll(cx)?;
    geometry.assemble(&[], VolumeLoad3::PerElement(&force_density_n_m3), &[], cx)
}
