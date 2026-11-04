//! Exact P1 charge moments on the original extruded planar complex.
//!
//! D=epsilon E+P, div(D)=rho and E=-grad(V) give
//! K V = integral(N rho) + integral_interface(N sigma)
//!       - integral_boundary(N D.n) + integral(grad(N).P).
//! Impressed P is prescribed in the mesh frame and independent of E; induced
//! polarization is already part of the sourced permittivity and must not be
//! supplied again. Reference: <https://mfem.org/electromagnetics/#electrostatics>.

use std::collections::BTreeMap;

use fs_rep_mesh::TriComplex2;

use super::{ElectricField, ElectrostaticError, accumulate, finite, geometry, poll};
use crate::planar::gradients;

/// Constant free surface charge on one original INTERNAL edge, extruded through
/// the mesh thickness. Vertex orientation does not change its physical sign.
/// With n from the first dielectric to the second, (D_second-D_first).n=sigma.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectrostaticChargeSheet {
    pub vertices: [usize; 2],
    pub free_charge_density_c_m2: f64,
}

/// Constant OUTWARD D.n [C/m2] on one original exterior edge. This specifies
/// displacement flux, not an equipotential floating conductor or its charge.
/// The corresponding weak nodal load has the opposite sign.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectrostaticDisplacementTrace {
    pub vertices: [usize; 2],
    pub outward_displacement_c_m2: f64,
}

/// Prescribed sources for one frozen material/geometry solve. All volume
/// quantities use original cell order; surface quantities name original edges.
/// A missing polarization slice means zero; an explicit slice must be complete.
#[derive(Clone, Copy, Debug)]
pub struct PlanarElectrostaticSources<'a> {
    pub free_charge_density_c_m3: &'a [f64],
    pub impressed_polarization_c_m2: Option<&'a [[f64; 2]]>,
    pub charge_sheets: &'a [ElectrostaticChargeSheet],
    pub displacement_traces: &'a [ElectrostaticDisplacementTrace],
}
impl<'a> PlanarElectrostaticSources<'a> {
    /// The original volume-charge-only problem, with zero impressed P and
    /// homogeneous natural displacement on unprescribed exterior edges.
    pub fn volume_charge(free_charge_density_c_m3: &'a [f64]) -> Self {
        Self {
            free_charge_density_c_m3,
            impressed_polarization_c_m2: None,
            charge_sheets: &[],
            displacement_traces: &[],
        }
    }
}

pub(super) struct SourceLoads {
    pub rhs: Vec<f64>,
    pub volume_charge_c: f64,
    pub sheet_charge_c: f64,
    pub outward_displacement_c: f64,
    sheets: Vec<ElectrostaticChargeSheet>,
    traces: Vec<ElectrostaticDisplacementTrace>,
    loaded_cells: Vec<bool>,
}

pub(super) fn assemble(
    mesh: &TriComplex2,
    sources: PlanarElectrostaticSources<'_>,
    prescribed: &[Option<f64>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SourceLoads, ElectrostaticError> {
    poll(cancelled)?;
    let cells = mesh.faces().len();
    if sources.free_charge_density_c_m3.len() != cells
        || sources
            .free_charge_density_c_m3
            .iter()
            .any(|v| !v.is_finite())
        || sources
            .impressed_polarization_c_m2
            .is_some_and(|p| p.len() != cells || p.iter().flatten().any(|v| !v.is_finite()))
        || sources.charge_sheets.len() > mesh.edges().len()
        || sources.displacement_traces.len() > mesh.edges().len()
    {
        return Err(ElectrostaticError::Invalid(
            "complete finite cell sources and bounded edge source counts required",
        ));
    }
    let mut rhs = vec![0.0; mesh.vertices().len()];
    let mut corrections = vec![0.0; rhs.len()];
    let mut volume_charge_c = 0.0;
    let mut volume_correction = 0.0;
    let has_sources = sources.impressed_polarization_c_m2.is_some()
        || !sources.charge_sheets.is_empty()
        || !sources.displacement_traces.is_empty();
    let mut loaded_cells = if has_sources {
        vec![false; cells]
    } else {
        Vec::new()
    };
    for (cell, face) in mesh.faces().iter().enumerate() {
        poll(cancelled)?;
        let volume = mesh.face_measure(cell).ok_or(ElectrostaticError::Invalid(
            "missing dielectric cell measure",
        ))?;
        let charge = finite(sources.free_charge_density_c_m3[cell] * volume)?;
        accumulate(&mut volume_charge_c, &mut volume_correction, charge)?;
        for &vertex in face {
            accumulate(
                &mut rhs[vertex as usize],
                &mut corrections[vertex as usize],
                charge / 3.0,
            )?;
        }
        if let Some(polarization) = sources.impressed_polarization_c_m2 {
            let p = polarization[cell];
            if p != [0.0; 2] {
                loaded_cells[cell] = true;
                let g = gradients(face.map(|v| mesh.vertices()[v as usize])).map_err(geometry)?;
                let b1 = finite(volume * g[1][0].mul_add(p[0], g[1][1] * p[1]))?;
                let b2 = finite(volume * g[2][0].mul_add(p[0], g[2][1] * p[1]))?;
                // Enforce the constant-potential nullspace locally: impressed
                // polarization creates bound charge, never net free charge.
                for (vertex, load) in face.iter().zip([-finite(b1 + b2)?, b1, b2]) {
                    accumulate(
                        &mut rhs[*vertex as usize],
                        &mut corrections[*vertex as usize],
                        load,
                    )?;
                }
            }
        }
    }
    let mut sheets = BTreeMap::new();
    let mut traces = BTreeMap::new();
    let mut edge_lookup = BTreeMap::new();
    if !sources.charge_sheets.is_empty() || !sources.displacement_traces.is_empty() {
        // Reuse the mesh's existing edge index and retain incidence only for
        // authored sources; do not rebuild a second full-mesh edge map.
        for mut vertices in sources
            .charge_sheets
            .iter()
            .map(|s| s.vertices)
            .chain(sources.displacement_traces.iter().map(|s| s.vertices))
        {
            poll(cancelled)?;
            vertices.sort_unstable();
            let index = u32::try_from(vertices[0])
                .ok()
                .zip(u32::try_from(vertices[1]).ok())
                .and_then(|(a, b)| mesh.edge_index(a, b))
                .ok_or(ElectrostaticError::Invalid(
                    "charge source must name an original mesh edge",
                ))?;
            edge_lookup
                .entry(vertices)
                .or_insert_with(|| (index, Vec::with_capacity(2)));
        }
        for (cell, face) in mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            for [a, b] in [[0, 1], [1, 2], [2, 0]] {
                let mut vertices = [face[a] as usize, face[b] as usize];
                vertices.sort_unstable();
                if let Some((_, owners)) = edge_lookup.get_mut(&vertices) {
                    owners.push(cell);
                }
            }
        }
    }
    for sheet in sources.charge_sheets {
        poll(cancelled)?;
        let mut vertices = sheet.vertices;
        vertices.sort_unstable();
        if !sheet.free_charge_density_c_m2.is_finite()
            || !edge_lookup
                .get(&vertices)
                .is_some_and(|(_, owners)| owners.len() == 2)
            || sheets
                .insert(vertices, sheet.free_charge_density_c_m2)
                .is_some()
        {
            return Err(ElectrostaticError::Invalid(
                "charge sheets require unique internal edges and finite surface densities",
            ));
        }
    }
    for trace in sources.displacement_traces {
        poll(cancelled)?;
        let mut vertices = trace.vertices;
        vertices.sort_unstable();
        if !trace.outward_displacement_c_m2.is_finite()
            || !edge_lookup
                .get(&vertices)
                .is_some_and(|(_, owners)| owners.len() == 1)
            || (prescribed[vertices[0]].is_some() && prescribed[vertices[1]].is_some())
            || traces
                .insert(vertices, trace.outward_displacement_c_m2)
                .is_some()
        {
            return Err(ElectrostaticError::Invalid(
                "displacement traces require unique exterior edges without a prescribed voltage trace",
            ));
        }
    }
    let mut sheet_charge_c = 0.0;
    let mut sheet_correction = 0.0;
    let mut outward_displacement_c = 0.0;
    let mut outward_correction = 0.0;
    // Sorted edges make permutation and endpoint reversal physically identical.
    for (is_sheet, entries) in [(true, &sheets), (false, &traces)] {
        for (&vertices, &density) in entries {
            poll(cancelled)?;
            let (edge, owners) = &edge_lookup[&vertices];
            let area = mesh.edge_measure(*edge).ok_or(ElectrostaticError::Invalid(
                "charge source lacks physical edge area",
            ))?;
            let charge = finite(area * density)?;
            let half = 0.5 * charge;
            if density != 0.0 && half == 0.0 {
                return Err(ElectrostaticError::Invalid(
                    "edge source charge is not representable",
                ));
            }
            if is_sheet {
                accumulate(&mut sheet_charge_c, &mut sheet_correction, charge)?;
            } else {
                accumulate(&mut outward_displacement_c, &mut outward_correction, charge)?;
            }
            for vertex in vertices {
                accumulate(
                    &mut rhs[vertex],
                    &mut corrections[vertex],
                    if is_sheet { half } else { -half },
                )?;
            }
            if density != 0.0 {
                for &cell in owners {
                    loaded_cells[cell] = true;
                }
            }
        }
    }
    for (value, correction) in rhs.iter_mut().zip(corrections) {
        poll(cancelled)?;
        *value = finite(*value + correction)?;
    }
    Ok(SourceLoads {
        rhs,
        volume_charge_c: finite(volume_charge_c + volume_correction)?,
        sheet_charge_c: finite(sheet_charge_c + sheet_correction)?,
        outward_displacement_c: finite(outward_displacement_c + outward_correction)?,
        sheets: sheets
            .into_iter()
            .map(
                |(vertices, free_charge_density_c_m2)| ElectrostaticChargeSheet {
                    vertices,
                    free_charge_density_c_m2,
                },
            )
            .collect(),
        traces: traces
            .into_iter()
            .map(
                |(vertices, outward_displacement_c_m2)| ElectrostaticDisplacementTrace {
                    vertices,
                    outward_displacement_c_m2,
                },
            )
            .collect(),
        loaded_cells,
    })
}

pub(super) fn retain(
    field: &mut ElectricField<'_>,
    sources: PlanarElectrostaticSources<'_>,
    loads: SourceLoads,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ElectrostaticError> {
    let mut pairing = 0.0;
    let mut correction = 0.0;
    if let Some(polarization) = sources.impressed_polarization_c_m2 {
        for (cell, &p) in polarization.iter().enumerate() {
            poll(cancelled)?;
            for (d, p) in field.displacement_c_m2[cell].iter_mut().zip(p) {
                *d = finite(*d + p)?;
            }
            let e = field.electric_field_v_m[cell];
            let volume =
                field
                    .problem
                    .mesh
                    .face_measure(cell)
                    .ok_or(ElectrostaticError::Invalid(
                        "missing polarization cell measure",
                    ))?;
            accumulate(
                &mut pairing,
                &mut correction,
                finite(volume * e[0].mul_add(p[0], e[1] * p[1]))?,
            )?;
        }
        field.impressed_polarization_c_m2 = Some(polarization.to_vec());
    }
    field.polarization_pairing_j = finite(pairing + correction)?;
    field.volume_charge_c = loads.volume_charge_c;
    field.sheet_charge_c = loads.sheet_charge_c;
    field.outward_displacement_c = loads.outward_displacement_c;
    field.free_charge_sheets = loads.sheets;
    field.displacement_traces = loads.traces;
    field.source_loaded_cells = loads.loaded_cells;
    Ok(())
}

impl ElectricField<'_> {
    /// Fixed impressed P in the original mesh frame [C/m2], or zero when absent.
    /// This excludes the induced polarization already contained in epsilon.
    pub fn impressed_polarization_c_m2(&self) -> Option<&[[f64; 2]]> {
        self.impressed_polarization_c_m2.as_deref()
    }
    pub fn free_charge_sheets(&self) -> &[ElectrostaticChargeSheet] {
        &self.free_charge_sheets
    }
    pub fn displacement_traces(&self) -> &[ElectrostaticDisplacementTrace] {
        &self.displacement_traces
    }
    /// Total authored INTERNAL free surface charge [C].
    pub fn sheet_charge_c(&self) -> f64 {
        self.sheet_charge_c
    }
    /// Total prescribed OUTWARD displacement flux [C]. Electrode reactions
    /// therefore obey Q_electrode+Q_volume+Q_sheet-Q_outward=0.
    pub fn outward_displacement_c(&self) -> f64 {
        self.outward_displacement_c
    }
    /// Integral E.P_impressed [J]. Its negative is the voltage pairing of the
    /// polarization load. It is not heat or an inferred polarization energy.
    pub fn polarization_pairing_j(&self) -> f64 {
        self.polarization_pairing_j
    }
}
