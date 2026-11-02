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
        let FrameFluidOperators { modal_mass, memory_masses, inverse_residual, system, witness_refusal } =
            assemble(&basis, &members, &stations, &material_ports, budget.mass_inverse_tolerance, gate)?;
        poll(gate)?;
        Ok(Self { frame, basis, members, stations, material_ports, modal_mass,
            memory_masses, inverse_residual, mass_inverse_tolerance: budget.mass_inverse_tolerance, system,
            witness_refusal })
    }
}

/// The fs-phs reconstruction gate's refusal text for a source factor that cannot rebuild `R` entry by entry.
const GATE_REFUSAL: &str = "physical resistance factor exceeds componentwise roundoff";

/// Why an assembled system carries no `ResistanceWitness`: the refusal fs-phs's componentwise reconstruction gate
/// returned for this assembly's own source squares, kept as the value it returned, with the uncompressed source rows it
/// was given (row-major, `n` entries per row). The system is then exactly the one built without a witness: no substitute
/// factor, no loosened gate. An exact flow over it factors `R` itself; fs-phs decides whether that builds.
#[derive(Clone, Debug)]
pub(super) struct WitnessRefusal {
    pub error: PhsError,
    #[cfg_attr(not(test), allow(dead_code))]
    pub rows: Vec<f64>,
}

/// The SAME original physical assembly serves preparation and a prescribed
/// material intervention. The retained basis, station traces and coordinate
/// topology are inputs; this does not compute replacement modes or poles.
pub(super) struct FrameFluidOperators {
    pub modal_mass: Vec<f64>,
    pub memory_masses: Vec<f64>,
    pub inverse_residual: f64,
    pub system: PortHamiltonian,
    pub witness_refusal: Option<WitnessRefusal>,
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
    // Every loss below is a physical square: the solid port and the low tail dissipate r (projection . qdot)^2, a memory
    // branch c (projection . qdot - v_branch)^2. R is their sum, so each square is also one row of a source factor
    // (sqrt of its coefficient times the same projection / slip vector). Retaining those rows in a checked
    // ResistanceWitness lets exact flows use the factor the assembly already holds instead of re-deriving one from the
    // summed R, which is not diagonally dominant here (a station projection above 1 makes c(1-projection) negative) and
    // so falls to an eigen factor that cannot meet fs-phs's componentwise reconstruction gate. The gate is unchanged.
    let mut source_rows: Vec<f64> = Vec::new();
    for port in material_ports {
        poll(gate)?;
        add_outer(&mut r, n, m, &port.projection, port.resistance)?;
        source_rows.extend(source_row(n, m, port.resistance, &port.projection, None)?);
    }
    for station in stations {
        poll(gate)?;
        let flow = &members[station.member].flow;
        let tail = positive(station.weight_m*flow.low_tail_resistance_per_length_n_s_m2())?;
        add_outer(&mut r, n, m, &station.projection, tail)?;
        source_rows.extend(source_row(n, m, tail, &station.projection, None)?);
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
            source_rows.extend(source_row(n, m, c, &station.projection, Some(index))?);
        }
    }
    let storage = QuadraticStorage::new(q, n)?;
    poll(gate)?;
    // The witness is evidence the assembly offers, not a precondition of the system: attach it when the fs-phs gate
    // admits these rows. When the gate's componentwise reconstruction check refuses them, build the system as before,
    // without a witness, and keep the refusal (see `gate_report` for what the gate objects to in the pinned fixtures).
    // Any other refusal (malformed rows, asymmetry) is a fault in this assembly and propagates.
    let attempt = fs_phs::ResistanceWitness::try_from_source_rows(n, r.clone(), source_rows.clone());
    let system = PortHamiltonian::new(n, m, j, r, g, Box::new(storage))?;
    let (system, witness_refusal) = match attempt {
        Ok(witness) => (system.with_resistance_witness(witness)?, None),
        Err(error @ PhsError::LinearSplitBudget { what }) if what == GATE_REFUSAL =>
            (system, Some(WitnessRefusal { error, rows: source_rows })),
        Err(other) => return Err(other.into()),
    };
    poll(gate)?;
    Ok(FrameFluidOperators { modal_mass, memory_masses, inverse_residual, system, witness_refusal })
}

/// One row of the source factor of a loss `coefficient * (projection . qdot - v_branch)^2`: sqrt(coefficient) times
/// the modal-momentum projection (state columns m..2m) and, for a memory branch, -1 on its own momentum column.
fn source_row(n: usize, m: usize, coefficient: f64, projection: &[f64], branch: Option<usize>)
    -> Result<Vec<f64>> {
    let weight = finite(fs_math::det::sqrt(coefficient))?;
    let mut row = vec![0.0; n];
    for (i, &p) in projection.iter().enumerate() { row[m+i] = finite(weight*p)?; }
    if let Some(index) = branch { row[index] = -weight; }
    Ok(row)
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

/// Test-only measurement of why fs-phs's reconstruction gate refuses an assembly's source rows. It repeats, step for
/// step, what `ResistanceWitness::try_from_source_rows` does (QR-compress when there are more rows than states, else
/// pad the rows to `n x n`; then the componentwise gate `|gram - R_ij| <= 128 n eps magnitude`, with the local product
/// sum as magnitude at an exactly zero `R_ij`), and reports the entries instead of refusing. It changes no gate.
#[cfg(test)]
pub(super) mod gate_report {
    use super::GATE_REFUSAL;
    use crate::air_path::cylinder::frame::MaterialFluidFrame;
    use fs_la::factor::qr;
    use fs_phs::PhsError;

    /// One entry of the reconstruction check, with the source terms that make it up.
    #[derive(Clone, Debug)]
    pub struct Entry {
        pub i: usize,
        pub j: usize,
        pub gram: f64,
        pub resistance: f64,
        pub magnitude: f64,
        pub error: f64,
        pub limit: f64,
        /// Largest and total absolute single source product `rows[k][i]*rows[k][j]` (the `c*p_i*p_j` terms).
        pub largest_term: f64,
        pub sum_abs_terms: f64,
    }

    #[derive(Debug)]
    pub struct Report {
        pub n: usize,
        pub height: usize,
        pub qr_path: bool,
        pub failing: usize,
        pub first_failure: Option<Entry>,
        pub worst: Option<Entry>,
    }

    pub fn report(n: usize, resistance: &[f64], rows: &[f64]) -> Report {
        assert!(n > 0 && rows.len() % n == 0 && resistance.len() == n*n);
        let height = rows.len()/n;
        let qr_path = height > n;
        let scale = rows.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
        let root: Vec<f64> = if qr_path {
            if scale == 0.0 { vec![0.0; n*n] } else {
                let normalized: Vec<f64> = rows.iter().map(|x| x/scale).collect();
                let decomposition = qr(&normalized, height, n);
                let mut root = vec![0.0; n*n];
                for i in 0..n { for j in i..n { root[i*n+j] = scale*decomposition.r(i, j); } }
                root
            }
        } else {
            let mut padded = rows.to_vec();
            padded.resize(n*n, 0.0);
            padded
        };
        let mut failing = 0;
        let (mut first_failure, mut worst): (Option<Entry>, Option<Entry>) = (None, None);
        let mut worst_violation = -1.0_f64;
        for i in 0..n {
            for j in i..n {
                let (mut gram, mut products) = (0.0, 0.0);
                for k in 0..n {
                    let product = root[k*n+i]*root[k*n+j];
                    gram += product;
                    products += product.abs();
                }
                let r = resistance[i*n+j];
                let magnitude = if r != 0.0 { r.abs() } else { products };
                let error = (gram - r).abs();
                let limit = 128.0*n as f64*f64::EPSILON*magnitude;
                let fails = !gram.is_finite() || !magnitude.is_finite() || error > limit;
                let (mut largest_term, mut sum_abs_terms) = (0.0_f64, 0.0_f64);
                for k in 0..height {
                    let term = (rows[k*n+i]*rows[k*n+j]).abs();
                    largest_term = largest_term.max(term);
                    sum_abs_terms += term;
                }
                let entry = Entry { i, j, gram, resistance: r, magnitude, error, limit, largest_term, sum_abs_terms };
                let violation = if limit > 0.0 { error/limit } else if error > 0.0 { f64::INFINITY } else { 0.0 };
                if fails {
                    failing += 1;
                    if first_failure.is_none() { first_failure = Some(entry.clone()); }
                }
                if violation > worst_violation { worst_violation = violation; worst = Some(entry); }
            }
        }
        Report { n, height, qr_path, failing, first_failure, worst }
    }

    fn describe(label: &str, which: &str, m: usize, report: &Report, e: &Entry) -> String {
        let block = |x: usize| if x < m { "position" } else if x < 2*m { "momentum" } else { "memory" };
        format!("GATE-REPORT {label} | {which}: m={m} n={} rows={} qr_path={} failing_entries={} \
            entry=({},{}) [{}x{}] gram={:e} R_ij={:e} |gram-R_ij|={:e} limit(128 n eps magnitude)={:e} \
            largest_term={:e} sum_abs_terms={:e} |R_ij|/largest_term={:e} |R_ij|/sum_abs_terms={:e}",
            report.n, report.height, report.qr_path, report.failing, e.i, e.j, block(e.i), block(e.j), e.gram,
            e.resistance, e.error, e.limit, e.largest_term, e.sum_abs_terms, e.resistance.abs()/e.largest_term,
            e.resistance.abs()/e.sum_abs_terms)
    }

    /// The pinned behaviour of an assembly whose gate refused its source rows: no witness is attached, the retained
    /// refusal is fs-phs's own value, and the replicated gate reproduces the refusal. Prints the measured entries.
    pub fn pin_refusal<F>(label: &str, model: &MaterialFluidFrame<'_, F>) -> Report {
        assert!(model.system().resistance_witness().is_none(), "{label}: a refused gate must leave no witness");
        let refusal = model.witness_refusal.as_ref().unwrap_or_else(|| panic!("{label}: the refusal must be retained"));
        assert_eq!(refusal.error, PhsError::LinearSplitBudget { what: GATE_REFUSAL });
        assert_eq!(model.resistance_witness_refusal(), Some(&refusal.error));
        let n = model.system().state_dim();
        let m = model.basis().modes().len();
        let (_, r, _) = model.system().structure();
        let report = report(n, r, &refusal.rows);
        eprintln!("GATE-REPORT {label} | refusal kept: {:?}", refusal.error);
        let first = report.first_failure.as_ref()
            .unwrap_or_else(|| panic!("{label}: the replicated gate must fail where fs-phs's did"));
        eprintln!("{}", describe(label, "first failing entry in gate order", m, &report, first));
        let worst = report.worst.as_ref().expect("worst entry");
        eprintln!("{}", describe(label, "worst entry (largest error/limit)", m, &report, worst));
        report
    }
}
