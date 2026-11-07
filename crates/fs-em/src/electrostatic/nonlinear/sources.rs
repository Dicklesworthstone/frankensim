//! Prescribed nonlinear sources use the original exact P1 charge integrator.
//! Source P is affine in the voltage objective and changes no material tangent.
use super::*;
use crate::electrostatic::{
    ElectrostaticChargeSheet, ElectrostaticDisplacementTrace, PlanarElectrostaticSources,
    sources::{SourceLoads, assemble, volume_nodal_loads},
};

pub(super) struct PreparedSources {
    pub loads: SourceLoads,
    pub component_charge_c: Vec<f64>,
    pub component_correction_c: Vec<f64>,
}

pub(super) struct RetainedSources {
    impressed_polarization_c_m2: Option<Vec<[f64; 2]>>,
    pub total_displacement_c_m2: Option<Vec<[f64; 2]>>,
    sheets: Vec<ElectrostaticChargeSheet>,
    traces: Vec<ElectrostaticDisplacementTrace>,
    loaded_cells: Vec<bool>,
    sheet_charge_c: f64,
    outward_displacement_c: f64,
    polarization_pairing_j: f64,
}

impl RetainedSources {
    pub fn cell_has_source(&self, cell: usize) -> bool {
        self.loaded_cells.get(cell).copied().unwrap_or(false)
    }
    pub fn has_geometry_sources(&self) -> bool {
        self.loaded_cells.iter().any(|&loaded| loaded)
    }
}

pub(super) fn prepare(
    electrodes: &NonlinearElectrodes<'_>,
    references: &electrodes::VoltageReferences,
    sources: PlanarElectrostaticSources<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedSources, Error> {
    let problem = electrodes.problem;
    // Both voltage and TOTAL-charge electrodes already determine their weak
    // surface reaction. A separate normal-displacement trace there would
    // prescribe a second, potentially incompatible local boundary condition.
    let prescribed: Vec<_> = electrodes
        .layout
        .node_terminal
        .iter()
        .map(|terminal| terminal.map(|_| 0.0))
        .collect();
    let loads = assemble(&problem.mesh, sources, &prescribed, cancelled)?;
    let mut component_charge_c = vec![0.0; electrodes.layout.terminals.len()];
    let mut component_correction_c = vec![0.0; component_charge_c.len()];
    // Integrate physical totals independently of the assembled nodal load.
    // P contributes no net free charge, even across material/source jumps.
    for (cell, (element, &density)) in problem
        .elements
        .iter()
        .zip(sources.free_charge_density_c_m3)
        .enumerate()
    {
        poll(cancelled)?;
        let charge = finite(element.volume_m3 * density)?;
        if density != 0.0
            && volume_nodal_loads(&problem.mesh, problem.mesh.faces()[cell], charge).contains(&0.0)
        {
            return Err(invalid("free charge load underflowed"));
        }
        let component = references.terminal_by_vertex[element.vertices[0]];
        add(
            &mut component_charge_c[component],
            &mut component_correction_c[component],
            charge,
        )?;
    }
    for (vertices, density) in loads
        .sheets
        .iter()
        .map(|sheet| (sheet.vertices, sheet.free_charge_density_c_m2))
        .chain(
            loads
                .traces
                .iter()
                .map(|trace| (trace.vertices, -trace.outward_displacement_c_m2)),
        )
    {
        poll(cancelled)?;
        // The shared source assembler has already admitted and canonicalized
        // these original edges, including their incidence and physical area.
        let edge = problem
            .mesh
            .edge_index(vertices[0] as u32, vertices[1] as u32)
            .ok_or_else(|| invalid("original source edge correspondence changed"))?;
        let area = problem
            .mesh
            .edge_measure(edge)
            .ok_or_else(|| invalid("original source surface measure is missing"))?;
        let component = references.terminal_by_vertex[vertices[0]];
        add(
            &mut component_charge_c[component],
            &mut component_correction_c[component],
            finite(area * density)?,
        )?;
    }
    Ok(PreparedSources {
        loads,
        component_charge_c,
        component_correction_c,
    })
}

pub(super) fn retain(
    problem: &PlanarNonlinearElectrostatics,
    evaluation: &Evaluation,
    sources: PlanarElectrostaticSources<'_>,
    loads: SourceLoads,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RetainedSources, Error> {
    let mut total_displacement_c_m2 = None;
    let (mut pairing, mut correction) = (0.0, 0.0);
    if let Some(polarization) = sources.impressed_polarization_c_m2 {
        let mut total = Vec::with_capacity(polarization.len());
        for (cell, &p) in polarization.iter().enumerate() {
            poll(cancelled)?;
            let d = evaluation.displacement[cell];
            total.push([finite(d[0] + p[0])?, finite(d[1] + p[1])?]);
            let e = evaluation.electric[cell];
            add(
                &mut pairing,
                &mut correction,
                problem.elements[cell].volume_m3 * e[0].mul_add(p[0], e[1] * p[1]),
            )?;
        }
        total_displacement_c_m2 = Some(total);
    }
    Ok(RetainedSources {
        impressed_polarization_c_m2: sources.impressed_polarization_c_m2.map(<[_]>::to_vec),
        total_displacement_c_m2,
        sheets: loads.sheets,
        traces: loads.traces,
        loaded_cells: loads.loaded_cells,
        sheet_charge_c: loads.sheet_charge_c,
        outward_displacement_c: loads.outward_displacement_c,
        polarization_pairing_j: finite(pairing + correction)?,
    })
}

impl NonlinearElectricField<'_> {
    /// Independent fixed P in the mesh frame, excluding constitutive response.
    pub fn impressed_polarization_c_m2(&self) -> Option<&[[f64; 2]]> {
        self.additional_sources
            .impressed_polarization_c_m2
            .as_deref()
    }
    pub fn free_charge_sheets(&self) -> &[ElectrostaticChargeSheet] {
        &self.additional_sources.sheets
    }
    pub fn displacement_traces(&self) -> &[ElectrostaticDisplacementTrace] {
        &self.additional_sources.traces
    }
    pub fn sheet_charge_c(&self) -> f64 {
        self.additional_sources.sheet_charge_c
    }
    /// Prescribed outward D.n integrated on the original physical boundary.
    /// Gauss balance is Q_electrode+Q_volume+Q_sheet-Q_outward=0.
    pub fn outward_displacement_c(&self) -> f64 {
        self.additional_sources.outward_displacement_c
    }
    /// Integral E.P_impressed [J], not inferred material storage or heat.
    pub fn polarization_pairing_j(&self) -> f64 {
        self.additional_sources.polarization_pairing_j
    }
}
