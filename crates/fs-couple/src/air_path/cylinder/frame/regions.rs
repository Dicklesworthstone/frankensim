//! Piecewise exterior source coverage on the original structural mesh.
use super::*;

/// One constant-fluid span in the oriented reference coordinate of an original
/// member. Require 0 <= start < end <= 1; overlapping spans of the same member
/// refuse, while shared endpoints are allowed. Gaps omit exterior loading:
/// they do not select an inferred vacuum or air material. Declaration order
/// controls source, memory and per-region loss ordering.
pub struct FrameFluidRegion<F> {
    pub element: usize,
    pub span: [f64; 2],
    pub fluid: F,
}
impl<F> From<FrameFluidBinding<F>> for FrameFluidRegion<F> {
    fn from(binding: FrameFluidBinding<F>) -> Self {
        Self { element: binding.element, span: [0.0, 1.0], fluid: binding.fluid }
    }
}

impl<'a, F: FluidTransport> MaterialFluidFrame<'a, F> {
    /// Wholly immersed source-defined members: the [0,1] case of the SAME
    /// regional assembly. Duplicate whole-member bindings refuse as overlaps.
    pub fn new(frame: &'a Frame2d, retained_modes: usize, bindings: Vec<FrameFluidBinding<F>>,
        resolution: CylinderFlowResolution, frame_budget: FrameVibrationBudget,
        budget: FrameFluidBudget, gate: &CancelGate) -> Result<Self> {
        Self::from_regions(frame, retained_modes, bindings.into_iter().map(FrameFluidRegion::from),
            resolution, frame_budget, budget, gate)
    }

    /// Partial immersion or heterogeneous exterior data without remeshing the
    /// solid, averaging transport coefficients, or multiplying whole-member
    /// damping by a coverage fraction. Four positive Gauss points integrate the
    /// original cubic-trace products on EACH span. Every spectral momentum and
    /// full added-mass cross term is retained by the original assembly/stepper.
    ///
    /// Counts are admitted before regional allocation or source evaluation.
    /// Overlap checking is O(regions log regions), sorting indices only. Inputs
    /// are frozen reference conditions. This local-cylinder approximation does
    /// not resolve free surfaces, menisci, interface waves, axial coupling,
    /// buoyancy or interactions between adjacent fluid regions. No fluid phase
    /// is inferred. Thermal and electrical consumers use this SAME prepared model.
    pub fn new_with_regions(frame: &'a Frame2d, retained_modes: usize, bindings: Vec<FrameFluidRegion<F>>,
        resolution: CylinderFlowResolution, frame_budget: FrameVibrationBudget,
        budget: FrameFluidBudget, gate: &CancelGate) -> Result<Self> {
        Self::from_regions(frame, retained_modes, bindings.into_iter(),
            resolution, frame_budget, budget, gate)
    }

    fn from_regions(frame: &'a Frame2d, retained_modes: usize,
        bindings: impl ExactSizeIterator<Item = FrameFluidRegion<F>>,
        resolution: CylinderFlowResolution, frame_budget: FrameVibrationBudget,
        budget: FrameFluidBudget, gate: &CancelGate) -> Result<Self> {
        poll(gate)?;
        resolution.validate().map_err(FrameFluidError::Fluid)?;
        if retained_modes == 0 || !budget.mass_inverse_tolerance.is_finite()
            || budget.mass_inverse_tolerance <= 0.0 || budget.mass_inverse_tolerance >= 1.0 {
            return Err(invalid("retained modes and a finite inverse residual allowance in (0,1) required"));
        }
        let memory_count = bindings.len().checked_mul(4)
            .and_then(|n| n.checked_mul(resolution.intervals))
            .ok_or_else(|| invalid("complete frame fluid-memory count overflow"))?;
        let n = retained_modes.checked_mul(2).and_then(|n| n.checked_add(memory_count))
            .ok_or_else(|| invalid("complete frame/fluid state count overflow"))?;
        let square = n.checked_mul(n).ok_or_else(|| invalid("frame/fluid matrix size overflow"))?;
        if n > budget.max_states || square > budget.max_dense_entries {
            return Err(invalid("complete frame/fluid state or dense matrix budget exceeded"));
        }
        let bindings: Vec<_> = bindings.collect();
        let lengths = admit_regions(frame, &bindings, gate)?;
        let basis = FrameVibrationBasis::from_frame(frame, retained_modes, frame_budget,
            || gate.is_requested())?;
        let m = retained_modes;
        let port_count = basis.material_damping().map_or(0, |d| d.ports().len());
        if bindings.len().checked_mul(4).and_then(|n| n.checked_add(port_count))
            .and_then(|n| n.checked_mul(m)).and_then(|n| n.checked_mul(6))
            .is_none_or(|n| n > budget.max_projection_terms) {
            return Err(invalid("frame/fluid physical projection budget exceeded"));
        }
        let mut members = Vec::with_capacity(bindings.len());
        let mut stations = Vec::with_capacity(bindings.len()*4);
        for (binding, length_m) in bindings.into_iter().zip(lengths) {
            poll(gate)?;
            let section = frame.material_section(binding.element)
                .ok_or_else(|| invalid("original circular material section is absent"))?;
            let flow = CylinderFlow::from_fluid(section.outer_radius_m(), binding.fluid, resolution, gate)
                .map_err(|error| match error {
                    AcousticRealizeError::Cancelled => FrameFluidError::Cancelled,
                    other => FrameFluidError::Fluid(other),
                })?;
            let member = members.len();
            for (xi, weight) in quadrature(binding.span)? {
                poll(gate)?;
                let (dofs, weights, _) = frame.reference_transverse_trace(binding.element, xi)?;
                let projection = project(&basis, &dofs, &weights)?;
                let weight_m = positive(weight*length_m)?;
                let first_memory = stations.len()*resolution.intervals;
                stations.push(Station { member, weight_m, projection, first_memory });
            }
            members.push(LoadedFrameMember { element: binding.element, span: binding.span, length_m, flow });
        }
        let mut material_ports = Vec::new();
        if let Some(damping) = basis.material_damping() {
            for port in damping.ports() {
                poll(gate)?;
                let resistance = finite(port.damping_n_s_m)?;
                if resistance < 0.0 { return Err(invalid("negative original material resistance")); }
                material_ports.push(MaterialPort {
                    projection: project(&basis, &port.dofs, &port.displacement_weights)?, resistance,
                });
            }
        }
        let FrameFluidOperators { modal_mass, memory_masses, inverse_residual, system } =
            assemble(&basis, &members, &stations, &material_ports, budget.mass_inverse_tolerance, gate)?;
        poll(gate)?;
        Ok(Self { frame, basis, members, stations, material_ports, modal_mass,
            memory_masses, inverse_residual, mass_inverse_tolerance: budget.mass_inverse_tolerance, system })
    }
}

/// The SAME original physical assembly serves preparation and a prescribed
/// material intervention. The retained basis, station traces and coordinate
/// topology are inputs; this does not compute replacement modes or poles.
pub(super) struct FrameFluidOperators {
    pub modal_mass: Vec<f64>,
    pub memory_masses: Vec<f64>,
    pub inverse_residual: f64,
    pub system: PortHamiltonian,
}

pub(super) fn assemble<F>(basis: &FrameVibrationBasis, members: &[LoadedFrameMember<F>],
    stations: &[Station], material_ports: &[MaterialPort], mass_inverse_tolerance: f64,
    gate: &CancelGate) -> Result<FrameFluidOperators> {
    poll(gate)?;
    let m = basis.modes().len();
    let mut modal_mass = vec![0.0; m*m];
    for i in 0..m { modal_mass[i*m+i] = 1.0; }
    let mut memory_masses = Vec::new();
    for station in stations {
        poll(gate)?;
        let flow = &members[station.member].flow;
        if station.first_memory != memory_masses.len() {
            return Err(invalid("original fluid spectral coordinate ordering changed"));
        }
        add_outer(&mut modal_mass, m, 0, &station.projection,
            positive(station.weight_m*flow.added_mass_per_length_kg_m())?)?;
        for branch in flow.branches() {
            let mass = positive(station.weight_m*branch.mass_per_length_kg_m)?;
            positive(1.0/mass)?;
            memory_masses.push(mass);
        }
    }
    let n = 2*m+memory_masses.len();
    let square = n.checked_mul(n).ok_or_else(|| invalid("frame/fluid matrix size overflow"))?;
    let (inverse, inverse_residual) = inverse_mass(&modal_mass, m, mass_inverse_tolerance, gate)?;
    let mut q = vec![0.0; square];
    let mut j = vec![0.0; square];
    let mut r = vec![0.0; square];
    let mut g = vec![0.0; n*m];
    for i in 0..m {
        q[i*n+i] = basis.modes()[i].eigenvalue_per_s2;
        j[i*n+m+i] = 1.0;
        j[(m+i)*n+i] = -1.0;
        g[(m+i)*m+i] = 1.0;
        for k in 0..m { q[(m+i)*n+m+k] = inverse[i*m+k]; }
    }
    for port in material_ports {
        poll(gate)?;
        add_outer(&mut r, n, m, &port.projection, port.resistance)?;
    }
    for station in stations {
        poll(gate)?;
        let flow = &members[station.member].flow;
        add_outer(&mut r, n, m, &station.projection,
            positive(station.weight_m*flow.low_tail_resistance_per_length_n_s_m2())?)?;
        for (k, branch) in flow.branches().iter().enumerate() {
            poll(gate)?;
            let memory = station.first_memory+k;
            let index = 2*m+memory;
            let c = positive(station.weight_m*branch.resistance_per_length_n_s_m2)?;
            q[index*n+index] = 1.0/memory_masses[memory];
            r[index*n+index] = c;
            add_outer(&mut r, n, m, &station.projection, c)?;
            for i in 0..m {
                let cross = finite(-c*station.projection[i])?;
                r[(m+i)*n+index] = cross;
                r[index*n+m+i] = cross;
            }
        }
    }
    let storage = QuadraticStorage::new(q, n)?;
    poll(gate)?;
    let system = PortHamiltonian::new(n, m, j, r, g, Box::new(storage))?;
    poll(gate)?;
    Ok(FrameFluidOperators { modal_mass, memory_masses, inverse_residual, system })
}

/// Map the quadrature inside the original coordinate, NOT a shorter new beam.
fn quadrature(span: [f64; 2]) -> Result<[(f64, f64); 4]> {
    let [start, end] = span;
    if !start.is_finite() || !end.is_finite() || start < 0.0 || end > 1.0 || start >= end {
        return Err(invalid("fluid region requires 0 <= start < end <= 1"));
    }
    let mapped = GAUSS.map(|(xi, weight)| ((end-start).mul_add(xi, start), weight));
    let mut previous = start;
    for (xi, _) in mapped {
        if xi <= previous || xi >= end {
            return Err(invalid("fluid region is too narrow for distinct interior quadrature points"));
        }
        previous = xi;
    }
    Ok(mapped)
}

fn admit_regions<F>(frame: &Frame2d, regions: &[FrameFluidRegion<F>], gate: &CancelGate) -> Result<Vec<f64>> {
    let mut lengths = Vec::with_capacity(regions.len());
    for region in regions {
        poll(gate)?;
        let points = quadrature(region.span)?;
        let (_, _, length) = frame.reference_transverse_trace(region.element, points[0].0)?;
        if frame.material_section(region.element).is_none() {
            return Err(invalid("fluid loading needs the original member material section and circular geometry"));
        }
        lengths.push(positive(length*(region.span[1]-region.span[0]))?);
    }
    // Only indices are sorted: retain original source objects and state order.
    let mut order: Vec<_> = (0..regions.len()).collect();
    order.sort_unstable_by(|&a, &b| regions[a].element.cmp(&regions[b].element)
        .then_with(|| regions[a].span[0].total_cmp(&regions[b].span[0]))
        .then_with(|| a.cmp(&b)));
    for pair in order.windows(2) {
        poll(gate)?;
        let (a, b) = (&regions[pair[0]], &regions[pair[1]]);
        if a.element == b.element && a.span[1] > b.span[0] {
            return Err(invalid("exterior fluid regions overlap on the same original member"));
        }
    }
    poll(gate)?;
    Ok(lengths)
}
