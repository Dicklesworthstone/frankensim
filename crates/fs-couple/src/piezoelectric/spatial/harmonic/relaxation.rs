//! Constitutive increments on the original tetrahedral work-conjugate maps.
//!
//! No element geometry, reference inertia or electrical coupling is rebuilt.
//! Frequency tensors are evaluated once per input material. Native gradients
//! and the original Mandel strain map then assemble symmetric complex blocks.
//! Recovery independently integrates physical stress/polarization increments;
//! it does not use the assembled matrix as its force/charge oracle.

use super::super::{Result, TetPiezoelectricModel, invalid, poll, strain_matrix};
use super::finite;
use fs_exec::Cx;
use fs_material::piezoelectric::tensor::relaxation::{
    PiezoelectricFrequencyCoefficients, RelaxingPiezoelectricHarmonicResponse,
    RelaxingTensorPiezoelectric,
};
use fs_math::c64::C64;

pub(super) fn prepare(
    model: &TetPiezoelectricModel<'_>,
    materials: &[RelaxingTensorPiezoelectric],
    omega: f64,
    cx: &Cx<'_>,
) -> Result<Vec<PiezoelectricFrequencyCoefficients>> {
    if materials.len() != 1 && materials.len() != model.tetrahedra.len() {
        return Err(invalid(
            "intrinsic relaxation requires one uniform law or one per original cell",
        ));
    }
    for cell in 0..model.tetrahedra.len() {
        poll(cx)?;
        let original = &model.materials[if model.materials.len() == 1 { 0 } else { cell }];
        let material = &materials[if materials.len() == 1 { 0 } else { cell }];
        if material.base().identity() != original.identity()
            || material.base().frame_identity() != model.mesh_frame
        {
            return Err(invalid(
                "relaxation must retain the original cell tensor source and mesh frame",
            ));
        }
    }
    materials
        .iter()
        .map(|material| {
            poll(cx)?;
            Ok(material.frequency_coefficients(omega)?)
        })
        .collect()
}

pub(super) fn stamp(
    model: &TetPiezoelectricModel<'_>,
    coefficients: &[PiezoelectricFrequencyCoefficients],
    add: &mut impl FnMut(usize, usize, C64, f64) -> Result<()>,
    cx: &Cx<'_>,
) -> Result<()> {
    let nu = 3 * model.nodes_m.len();
    for (cell, tet) in model.tetrahedra.iter().enumerate() {
        poll(cx)?;
        let coefficient = &coefficients[if coefficients.len() == 1 { 0 } else { cell }];
        let gradients = &model.dielectric.gradients()[cell];
        let b = strain_matrix(gradients);
        let volume = model.dielectric.volumes_m3()[cell];
        // Form C_delta B once; then only the upper triangle of B^T C_delta B.
        // Mirror the exact same rounded value so native MINRES sees symmetry.
        let cb: [[C64; 12]; 6] = core::array::from_fn(|row| {
            core::array::from_fn(|column| {
                (0..6).fold(C64::ZERO, |sum, k| {
                    sum + coefficient.stiffness_increment_mandel_pa[row][k].scale(b[k][column])
                })
            })
        });
        for row in 0..12 {
            let global_row = 3 * tet[row / 3] + row % 3;
            for column in row..12 {
                let global_column = 3 * tet[column / 3] + column % 3;
                let value = (0..6)
                    .fold(C64::ZERO, |sum, k| sum + cb[k][column].scale(b[k][row]))
                    .scale(volume);
                add(global_row, global_column, value, value.abs())?;
                if row != column {
                    add(global_column, global_row, value, value.abs())?;
                }
            }
        }
        let eg: [[C64; 3]; 4] = core::array::from_fn(|node| {
            core::array::from_fn(|axis| {
                (0..3).fold(C64::ZERO, |sum, k| {
                    sum + coefficient.permittivity_increment_f_m[axis][k].scale(gradients[node][k])
                })
            })
        });
        for row in 0..4 {
            for column in row..4 {
                let value = (0..3)
                    .fold(C64::ZERO, |sum, k| {
                        sum + eg[column][k].scale(gradients[row][k])
                    })
                    .scale(-volume);
                add(nu + tet[row], nu + tet[column], value, value.abs())?;
                if row != column {
                    add(nu + tet[column], nu + tet[row], value, value.abs())?;
                }
            }
        }
    }
    Ok(())
}

pub(super) struct Recovery {
    pub(super) force_n: Vec<[C64; 3]>,
    pub(super) charge_c: Vec<C64>,
    pub(super) cells: Vec<RelaxingPiezoelectricHarmonicResponse>,
}

pub(super) fn recover(
    model: &TetPiezoelectricModel<'_>,
    materials: &[RelaxingTensorPiezoelectric],
    omega: f64,
    full: &[C64],
    cx: &Cx<'_>,
) -> Result<Recovery> {
    let n = model.nodes_m.len();
    let nu = 3 * n;
    let mut force_n = vec![[C64::ZERO; 3]; n];
    let mut charge_c = vec![C64::ZERO; n];
    let mut cells = Vec::with_capacity(model.tetrahedra.len());
    for (cell, tet) in model.tetrahedra.iter().enumerate() {
        poll(cx)?;
        let material = &materials[if materials.len() == 1 { 0 } else { cell }];
        let gradients = &model.dielectric.gradients()[cell];
        let b = strain_matrix(gradients);
        let volume = model.dielectric.volumes_m3()[cell];
        let mut strain = [C64::ZERO; 6];
        let mut field = [C64::ZERO; 3];
        for axis in 0..6 {
            for local in 0..12 {
                strain[axis] = strain[axis]
                    + (full[3 * tet[local / 3] + local % 3] - full[3 * tet[0] + local % 3])
                        .scale(b[axis][local]);
            }
        }
        for local in 1..4 {
            let potential = full[nu + tet[local]] - full[nu + tet[0]];
            for axis in 0..3 {
                field[axis] = field[axis] - potential.scale(gradients[local][axis]);
            }
        }
        let response = material.harmonic_response(omega, strain, field)?;
        for local in 0..12 {
            let value = (0..6)
                .fold(C64::ZERO, |sum, axis| {
                    sum + response.additional_stress_mandel_pa[axis].scale(b[axis][local])
                })
                .scale(volume);
            let node = tet[local / 3];
            force_n[node][local % 3] = force_n[node][local % 3] + value;
        }
        for local in 0..4 {
            let value = (0..3)
                .fold(C64::ZERO, |sum, axis| {
                    sum + response.polarization_c_m2[axis].scale(gradients[local][axis])
                })
                .scale(-volume);
            charge_c[tet[local]] = charge_c[tet[local]] + value;
        }
        cells.push(response);
    }
    if force_n
        .iter()
        .flatten()
        .chain(&charge_c)
        .any(|&z| !finite(z))
    {
        return Err(invalid(
            "intrinsic material force or polarization charge is not representable",
        ));
    }
    poll(cx)?;
    Ok(Recovery {
        force_n,
        charge_c,
        cells,
    })
}
