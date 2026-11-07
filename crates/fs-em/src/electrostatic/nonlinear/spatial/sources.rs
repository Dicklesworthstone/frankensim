//! Native tetrahedral source preparation and independent impressed-source work.
use super::*;
use crate::electrostatic::ElectrostaticControls;
use crate::spatial::electrostatic::{
    ElectrostaticVoltageFace, PreparedSources as NativeSources, SpatialElectrostaticBoundary,
    prepare_sources_on_geometry,
};

pub(super) struct PreparedSources {
    pub rhs: Vec<f64>,
    pub component_charge_c: Vec<f64>,
    pub component_correction_c: Vec<f64>,
    pub volume_charge_c: f64,
    pub sheet_charge_c: f64,
    pub outward_displacement_c: f64,
    native: NativeSources,
    surface_load_c: Vec<f64>,
    traces: Vec<ElectrostaticDisplacementFace>,
}

pub(super) struct Response {
    pub component_polarization_pairing_j: Vec<f64>,
    pub polarization_pairing_j: f64,
    pub polarization_work_residual_j: f64,
    pub worst_component_polarization_work_residual_j: f64,
    total_displacement_c_m2: Option<Vec<[f64; 3]>>,
}

pub(super) struct RetainedSources {
    polarization: Option<Vec<[f64; 3]>>,
    pub total_displacement_c_m2: Option<Vec<[f64; 3]>>,
    sheets: Vec<ElectrostaticChargeSheet>,
    traces: Vec<ElectrostaticDisplacementFace>,
    sheet_charge_c: f64,
    outward_displacement_c: f64,
    polarization_pairing_j: f64,
    polarization_work_residual_j: f64,
    maximum_component_polarization_work_residual_j: f64,
    has_nonzero: bool,
}

pub(super) fn prepare(
    electrodes: &SpatialNonlinearElectrodes<'_>,
    references: &VoltageReferences,
    sources: SpatialElectrostaticSources<'_>,
    displacement_faces: &[ElectrostaticDisplacementFace],
    settings: NonlinearElectrostaticControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedSources, Error> {
    poll(cancelled)?;
    let problem = electrodes.problem;
    if displacement_faces.len() > problem.mesh().faces.len() {
        return Err(invalid("bounded original displacement faces required"));
    }
    for trace in displacement_faces {
        poll(cancelled)?;
        let vertices = problem
            .mesh()
            .faces
            .get(trace.face)
            .ok_or_else(|| invalid("displacement face lies outside the original mesh"))?;
        if vertices
            .iter()
            .all(|&v| electrodes.layout.node_terminal[v as usize].is_some())
        {
            return Err(invalid(
                "a displacement face cannot also carry a complete voltage or floating electrode trace",
            ));
        }
    }
    let mut traces = displacement_faces.to_vec();
    traces.sort_unstable_by_key(|trace| trace.face);
    // These zero traces admit original geometry only. Actual voltages, total
    // floating charges and electrical gauges belong to Layout and Gauss. A
    // floating metal is also a physical trace, never an independent D.n source.
    let boundary = SpatialElectrostaticBoundary {
        voltage: electrodes
            .faces
            .iter()
            .flatten()
            .map(|&face| ElectrostaticVoltageFace {
                face,
                voltage_v: [0.0; 3],
            })
            .collect(),
        displacement: traces.clone(),
    };
    let (native, boundary) = prepare_sources_on_geometry(
        &problem.geometry,
        sources,
        &boundary,
        ElectrostaticControls {
            linear: settings.linear,
            charge_tolerance_c: settings.charge_tolerance_c,
            energy_tolerance_j: settings.energy_tolerance_j,
        },
        cancelled,
    )
    .map_err(trace_error)?;
    let surface_load_c = boundary.into_flux_load();
    let rhs = native
        .weak_load_c()
        .iter()
        .zip(&surface_load_c)
        .map(|(&volume, &surface)| finite(volume + surface))
        .collect::<Result<Vec<_>, _>>()?;
    let mut component_charge_c = vec![0.0; electrodes.layout.terminals.len()];
    let mut component_correction_c = vec![0.0; component_charge_c.len()];
    let (mut volume_charge_c, mut volume_correction) = (0.0, 0.0);
    let (mut sheet_charge_c, mut sheet_correction) = (0.0, 0.0);
    let (mut outward_displacement_c, mut outward_correction) = (0.0, 0.0);
    // Physical totals are independently integrated, not obtained by summing
    // the assembled RHS (which also contains polarization and outward signs).
    for (cell, &density) in sources.free_charge_density_c_m3.iter().enumerate() {
        poll(cancelled)?;
        let element = problem.element(cell);
        let charge = finite(element.volume_m3 * density)?;
        if density != 0.0 && 0.25 * charge == 0.0 {
            return Err(invalid("tetrahedral free-charge load underflowed"));
        }
        let reference = references.terminal_by_vertex[element.vertices[0] as usize];
        add(
            &mut component_charge_c[reference],
            &mut component_correction_c[reference],
            charge,
        )?;
        add(&mut volume_charge_c, &mut volume_correction, charge)?;
    }
    for sheet in sources.sheets {
        poll(cancelled)?;
        let charge = *native
            .sheet_charge_by_face_c()
            .get(sheet.face)
            .ok_or_else(|| invalid("original free-sheet correspondence changed"))?;
        if sheet.free_charge_density_c_m2 != 0.0 && charge / 3.0 == 0.0 {
            return Err(invalid("free-sheet face moment underflowed"));
        }
        let reference = references.terminal_by_vertex[problem.mesh().faces[sheet.face][0] as usize];
        add(
            &mut component_charge_c[reference],
            &mut component_correction_c[reference],
            charge,
        )?;
        add(&mut sheet_charge_c, &mut sheet_correction, charge)?;
    }
    for trace in &traces {
        poll(cancelled)?;
        let area = problem
            .geometry
            .face_area_m2(trace.face)
            .ok_or_else(|| invalid("original displacement surface measure is missing"))?;
        let charge = finite(area * trace.outward_displacement_c_m2)?;
        if trace.outward_displacement_c_m2 != 0.0 && charge / 3.0 == 0.0 {
            return Err(invalid("outward displacement face moment underflowed"));
        }
        let reference = references.terminal_by_vertex[problem.mesh().faces[trace.face][0] as usize];
        add(
            &mut component_charge_c[reference],
            &mut component_correction_c[reference],
            -charge,
        )?;
        add(&mut outward_displacement_c, &mut outward_correction, charge)?;
    }
    poll(cancelled)?;
    Ok(PreparedSources {
        rhs,
        component_charge_c,
        component_correction_c,
        volume_charge_c: finite(volume_charge_c + volume_correction)?,
        sheet_charge_c: finite(sheet_charge_c + sheet_correction)?,
        outward_displacement_c: finite(outward_displacement_c + outward_correction)?,
        native,
        surface_load_c,
        traces,
    })
}

impl PreparedSources {
    /// Free-charge and outward-displacement work, excluding impressed P.
    /// Retaining the original separate arrays avoids subtracting a large P
    /// load from a much smaller physical free-charge source after rounding.
    pub fn physical_load_c(&self, vertex: usize) -> Result<f64, Error> {
        let (mut sum, mut correction) = (self.native.volume_load_c()[vertex], 0.0);
        add(
            &mut sum,
            &mut correction,
            self.native
                .sheet_load_c()
                .get(vertex)
                .copied()
                .unwrap_or(0.0),
        )?;
        add(&mut sum, &mut correction, self.surface_load_c[vertex])?;
        Ok(finite(sum + correction)?)
    }

    pub fn response(
        &self,
        problem: &SpatialNonlinearElectrostatics,
        evaluation: &Evaluation,
        potential: &[f64],
        references: &VoltageReferences,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Response, Error> {
        let count = self.component_charge_c.len();
        let mut component_pairing = vec![0.0; count];
        let mut component_correction = vec![0.0; count];
        let (mut pairing, mut correction) = (0.0, 0.0);
        let mut total_displacement_c_m2 = None;
        if let Some(polarization) = self.native.polarization_c_m2() {
            let mut total = Vec::with_capacity(polarization.len());
            for (cell, &p) in polarization.iter().enumerate() {
                poll(cancelled)?;
                let element = problem.element(cell);
                let d = evaluation.displacement[cell];
                total.push([
                    finite(d[0] + p[0])?,
                    finite(d[1] + p[1])?,
                    finite(d[2] + p[2])?,
                ]);
                let work = finite(element.volume_m3 * dot(evaluation.electric[cell], p))?;
                let reference = references.terminal_by_vertex[element.vertices[0] as usize];
                add(
                    &mut component_pairing[reference],
                    &mut component_correction[reference],
                    work,
                )?;
                add(&mut pairing, &mut correction, work)?;
            }
            total_displacement_c_m2 = Some(total);
        }
        for (value, correction) in component_pairing.iter_mut().zip(&component_correction) {
            *value = finite(*value + *correction)?;
        }
        let polarization_pairing_j = finite(pairing + correction)?;
        let mut work = component_pairing.clone();
        component_correction.fill(0.0);
        for (vertex, &load) in self.native.polarization_load_c().iter().enumerate() {
            poll(cancelled)?;
            let reference = references.terminal_by_vertex[vertex];
            add(
                &mut work[reference],
                &mut component_correction[reference],
                potential[vertex] * load,
            )?;
        }
        let (mut total_work, mut work_correction) = (0.0, 0.0);
        let mut worst = 0.0_f64;
        for (&work, &correction) in work.iter().zip(&component_correction) {
            poll(cancelled)?;
            let residual = finite(work + correction)?;
            if residual.abs() > worst.abs() {
                worst = residual;
            }
            add(&mut total_work, &mut work_correction, residual)?;
        }
        Ok(Response {
            component_polarization_pairing_j: component_pairing,
            polarization_pairing_j,
            polarization_work_residual_j: finite(total_work + work_correction)?,
            worst_component_polarization_work_residual_j: worst,
            total_displacement_c_m2,
        })
    }

    pub fn retain(self, response: Response) -> (Vec<f64>, RetainedSources) {
        let (rho, polarization, sheets) = self.native.into_declarations();
        let has_nonzero = polarization
            .as_ref()
            .is_some_and(|p| p.iter().flatten().any(|&v| v != 0.0))
            || sheets
                .iter()
                .any(|sheet| sheet.free_charge_density_c_m2 != 0.0)
            || self
                .traces
                .iter()
                .any(|trace| trace.outward_displacement_c_m2 != 0.0);
        (
            rho,
            RetainedSources {
                polarization,
                total_displacement_c_m2: response.total_displacement_c_m2,
                sheets,
                traces: self.traces,
                sheet_charge_c: self.sheet_charge_c,
                outward_displacement_c: self.outward_displacement_c,
                polarization_pairing_j: response.polarization_pairing_j,
                polarization_work_residual_j: response.polarization_work_residual_j,
                maximum_component_polarization_work_residual_j: response
                    .worst_component_polarization_work_residual_j
                    .abs(),
                has_nonzero,
            },
        )
    }
}

impl SpatialNonlinearElectricField<'_> {
    /// Fixed mesh-frame P, separate from all polarization in the material law.
    pub fn impressed_polarization_c_m2(&self) -> Option<&[[f64; 3]]> {
        self.additional_sources.polarization.as_deref()
    }
    pub fn charge_sheets(&self) -> &[ElectrostaticChargeSheet] {
        &self.additional_sources.sheets
    }
    pub fn displacement_faces(&self) -> &[ElectrostaticDisplacementFace] {
        &self.additional_sources.traces
    }
    pub fn sheet_charge_c(&self) -> f64 {
        self.additional_sources.sheet_charge_c
    }
    /// Integral of prescribed outward D.n. Physical Gauss balance is
    /// Q_electrode + Q_volume + Q_sheet - Q_outward = 0.
    pub fn outward_displacement_c(&self) -> f64 {
        self.additional_sources.outward_displacement_c
    }
    /// Integral E.P_impressed [J], distinct from material storage or heat.
    pub fn polarization_pairing_j(&self) -> f64 {
        self.additional_sources.polarization_pairing_j
    }
    /// Direct E.P plus voltage paired with the original weak polarization load.
    pub fn polarization_work_residual_j(&self) -> f64 {
        self.additional_sources.polarization_work_residual_j
    }
    pub fn maximum_component_polarization_work_residual_j(&self) -> f64 {
        self.additional_sources
            .maximum_component_polarization_work_residual_j
    }
    /// Nonzero fixed P, sheets or outward displacement require a separate
    /// transport law for shape derivatives. Conserved volume charge is excluded.
    pub fn has_additional_sources(&self) -> bool {
        self.additional_sources.has_nonzero
    }
}
