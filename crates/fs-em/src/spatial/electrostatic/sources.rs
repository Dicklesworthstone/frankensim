//! Original tetrahedral free-charge, polarization and boundary source moments.
//! Preparation depends on geometry and physical traces, never a linear Hodge.
use super::*;

pub(crate) struct PreparedSources {
    pub(super) rho: Vec<f64>,
    pub(super) polarization: Option<Vec<[f64; 3]>>,
    pub(super) sheets: Vec<ElectrostaticChargeSheet>,
    pub(super) volume_load: Vec<f64>,
    pub(super) sheet_load: Vec<f64>,
    pub(super) polarization_load: Vec<f64>,
    pub(super) total_load: Vec<f64>,
    pub(super) sheet_charge_by_face: Vec<f64>,
    pub(super) charged_sheet_vertices: Vec<bool>,
}

impl PreparedSources {
    pub(crate) fn weak_load_c(&self) -> &[f64] {
        if self.total_load.is_empty() {
            &self.volume_load
        } else {
            &self.total_load
        }
    }
}

impl PreparedSources {
    pub(crate) fn volume_load_c(&self) -> &[f64] {
        &self.volume_load
    }
    pub(crate) fn sheet_load_c(&self) -> &[f64] {
        &self.sheet_load
    }
    pub(crate) fn polarization_load_c(&self) -> &[f64] {
        &self.polarization_load
    }
    pub(crate) fn polarization_c_m2(&self) -> Option<&[[f64; 3]]> {
        self.polarization.as_deref()
    }
    pub(crate) fn sheet_charge_by_face_c(&self) -> &[f64] {
        &self.sheet_charge_by_face
    }
    pub(crate) fn into_declarations(
        self,
    ) -> (
        Vec<f64>,
        Option<Vec<[f64; 3]>>,
        Vec<ElectrostaticChargeSheet>,
    ) {
        (self.rho, self.polarization, self.sheets)
    }
}

pub(crate) fn prepare(
    geometry: &SpatialGeometry,
    sources: SpatialElectrostaticSources<'_>,
    boundary: &SpatialElectrostaticBoundary,
    controls: ElectrostaticControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(PreparedSources, ScalarBoundary)> {
    check(cancelled)?;
    let rho = sources.free_charge_density_c_m3;
    if rho.len() != geometry.mesh.tets.len()
        || rho.iter().any(|v| !v.is_finite())
        || sources.impressed_polarization_c_m2.is_some_and(|p| {
            p.len() != geometry.mesh.tets.len() || p.iter().flatten().any(|v| !v.is_finite())
        })
        || sources.sheets.len() > geometry.faces.len()
        || boundary.voltage.len() > geometry.faces.len()
        || boundary.displacement.len() > geometry.faces.len()
        || !positive(controls.linear.tolerance)
        || controls.linear.tolerance >= 1.0
        || !positive(controls.charge_tolerance_c)
        || !positive(controls.energy_tolerance_j)
    {
        return Err(invalid(
            "complete finite cell sources, bounded sheets/traces and positive tolerances required",
        ));
    }
    let g = geometry;
    let mut volume_load = vec![0.0; g.mesh.vertex_count];
    let mut corrections = vec![0.0; g.mesh.vertex_count];
    for (cell, tet) in g.mesh.tets.iter().enumerate() {
        check(cancelled)?;
        let load = finite(rho[cell] * g.geometry.vol_signed[cell].abs() / 4.0)?;
        for &v in tet {
            add(
                &mut volume_load[v as usize],
                &mut corrections[v as usize],
                load,
            )?;
        }
    }
    let mut sheets = sources.sheets.to_vec();
    sheets.sort_unstable_by_key(|s| s.face);
    let mut sheet_load = Vec::new();
    let mut sheet_charge_by_face = Vec::new();
    let mut charged_sheet_vertices = Vec::new();
    if !sheets.is_empty() {
        sheet_load.resize(g.mesh.vertex_count, 0.0);
        sheet_charge_by_face.resize(g.faces.len(), 0.0);
        charged_sheet_vertices.resize(g.mesh.vertex_count, false);
        corrections.fill(0.0);
    }
    for (i, sheet) in sheets.iter().enumerate() {
        check(cancelled)?;
        if g.faces.get(sheet.face).is_none_or(|f| f.cells.len() != 2)
            || i > 0 && sheets[i - 1].face == sheet.face
            || !sheet.free_charge_density_c_m2.is_finite()
        {
            return Err(invalid(
                "finite free sheets require distinct original internal faces",
            ));
        }
        let charge = finite(sheet.free_charge_density_c_m2 * g.faces[sheet.face].area_m2)?;
        sheet_charge_by_face[sheet.face] = charge;
        for v in g.mesh.faces[sheet.face] {
            let v = v as usize;
            add(&mut sheet_load[v], &mut corrections[v], charge / 3.0)?;
            charged_sheet_vertices[v] |= sheet.free_charge_density_c_m2 != 0.0;
        }
    }
    // No extra nodal/face arrays are needed for the original volume-only
    // path. When sheets are present, reuse their complete free-charge
    // buffer for the eventual polarized weak load after gauge admission.
    let mut total_load = if sheets.is_empty() {
        Vec::new()
    } else {
        volume_load
            .iter()
            .zip(&sheet_load)
            .map(|(&volume, &sheet)| finite(volume + sheet))
            .collect::<Result<Vec<_>>>()?
    };
    let traces = prepare_scalar_boundary_traces_with_sources(
        &g.mesh,
        &g.faces,
        &g.components,
        boundary.voltage.iter().map(|b| (b.face, b.voltage_v)),
        boundary
            .displacement
            .iter()
            .map(|b| (b.face, -b.outward_displacement_c_m2)),
        // Polarization creates no net free charge. Check the physical
        // compatibility of volume + sheet charge against outward D.n;
        // never let a numerical gauge hide uncompensated free charge.
        Some(if total_load.is_empty() {
            &volume_load
        } else {
            &total_load
        }),
        controls.charge_tolerance_c,
        cancelled,
    )?;
    let mut polarization_load = Vec::new();
    if let Some(polarization) = sources.impressed_polarization_c_m2 {
        polarization_load.resize(g.mesh.vertex_count, 0.0);
        corrections.fill(0.0);
        for (cell, tet) in g.mesh.tets.iter().enumerate() {
            check(cancelled)?;
            let (mut first, mut first_correction) = (0.0, 0.0);
            // The constant test function has zero gradient. Form its
            // zeroth contribution from the other three, matching the
            // voltage-difference gradient used by field recovery.
            for local in 1..4 {
                let term = finite(
                    g.geometry.vol_signed[cell].abs()
                        * dot(g.geometry.grads[cell][local], polarization[cell]),
                )?;
                let v = tet[local] as usize;
                add(&mut polarization_load[v], &mut corrections[v], term)?;
                add(&mut first, &mut first_correction, -term)?;
            }
            let v = tet[0] as usize;
            add(&mut polarization_load[v], &mut corrections[v], first)?;
        }
        if total_load.is_empty() {
            total_load.clone_from(&volume_load);
        }
        for (load, &polarized) in total_load.iter_mut().zip(&polarization_load) {
            check(cancelled)?;
            *load = finite(*load + polarized)?;
        }
    }
    check(cancelled)?;
    Ok((
        PreparedSources {
            rho: rho.to_vec(),
            polarization: sources.impressed_polarization_c_m2.map(<[_]>::to_vec),
            sheets,
            volume_load,
            sheet_load,
            polarization_load,
            total_load,
            sheet_charge_by_face,
            charged_sheet_vertices,
        },
        traces,
    ))
}
