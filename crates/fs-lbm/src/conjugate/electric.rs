//! DC electric conduction and Joule heating on the voxel domain.
//!
//! Solid materials may conduct electricity ([`Conductor`]: a resistivity
//! with a linear temperature coefficient). [`ElectricSetup::solve_with_resistivity`]
//! admits a caller's state-resolved constitutive law instead, without replacing
//! the current-balance or Joule-heating discretization. Electrodes are boxes of
//! cells treated as perfect conductors, one potential each, held at a voltage or
//! carrying a declared total current. The potential solves the discrete
//! current balance `sum_f G_f (phi_P - phi_N) = 0` on every conductor cell,
//! with `G_f = A_f / (rho_P w_P / 2 + rho_N w_N / 2 + R_contact)` between
//! conductor cells; an electrode's bulk half resistance is zero. Contacts are
//! optional explicit [`ElectricalContact`] data. Other faces retain the original
//! perfect-contact rule; unassigned cells and the domain boundary carry no current.
//! The legacy paths insulate all fluids; [`ElectricSetup::solve_with_cell_resistivity`]
//! instead lets explicit material data admit either fluid or solid cells without
//! changing their mechanical occupancy or introducing a second electrical solver.
//! Each face dissipates `G_f (phi_P - phi_N)^2`. Bulk heat stays in its own
//! half cells; interface heat follows the declared surface partition. Thus the
//! returned Joule heat equals electrode power `sum_e V_e I_e` to solver precision.
//!
//! [`solve_electrothermal`] couples the two fields: resistivity at each
//! cell's temperature, the potential, the Joule heat as a source, the
//! temperature, repeated to a fixed point (one pass when no conductor
//! depends on temperature).
//!
//! No-claims: DC only (no skin effect, inductance or capacitance); ohmic,
//! isotropic bulk conductors and explicit zero-thickness resistive contacts;
//! electrodes are ideal equipotentials; no contact-area or heat-partition
//! model is inferred. The Picard iteration does not converge past a thermal
//! runaway (current-driven conductors whose resistivity rises with temperature
//! can have no steady state, and that refusal is the answer).

use fs_exec::CancelGate;
use fs_sparse::Coo;
use fs_sparse::precond::{SaAmg, pcg};

use super::domain::{Voxel, VoxelDomain};
use super::energy::{EnergySolution, ThermalSetup};
use super::flow::scale_rows;
use super::krylov::bicgstab_ilu0;
use super::{ChtError, finite, finite_positive, poll};

mod contact;
pub mod motion;
pub use contact::ElectricalContact;

/// An ohmic conductor: `rho(T) = rho_0 (1 + alpha (T - T_ref))`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Conductor {
    /// Resistivity at the reference temperature, ohm m.
    pub resistivity_ohm_m: f64,
    /// Linear temperature coefficient `alpha`, 1/K (copper: 3.9e-3).
    pub temperature_coefficient_per_k: f64,
    /// Reference temperature of `resistivity_ohm_m`, K.
    pub reference_temperature_k: f64,
}

impl Conductor {
    /// A temperature-independent conductor.
    #[must_use]
    pub const fn new(resistivity_ohm_m: f64) -> Self {
        Self {
            resistivity_ohm_m,
            temperature_coefficient_per_k: 0.0,
            reference_temperature_k: 293.15,
        }
    }

    /// The same conductor with a linear temperature coefficient about
    /// `reference_temperature_k`.
    #[must_use]
    pub const fn with_temperature_coefficient(
        mut self,
        per_k: f64,
        reference_temperature_k: f64,
    ) -> Self {
        self.temperature_coefficient_per_k = per_k;
        self.reference_temperature_k = reference_temperature_k;
        self
    }

    /// Resistivity at `t`, ohm m.
    #[must_use]
    pub fn resistivity(&self, t: f64) -> f64 {
        self.resistivity_ohm_m
            * self
                .temperature_coefficient_per_k
                .mul_add(t - self.reference_temperature_k, 1.0)
    }
}

/// What an electrode imposes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ElectrodeDrive {
    /// Held at this potential, V.
    Voltage(f64),
    /// Delivers this total current into the conductors, A
    /// (negative: draws it); its potential is solved.
    Current(f64),
}

/// A perfectly conducting box of cells `lo..hi`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Electrode {
    /// First cells of the box (inclusive).
    pub lo: [usize; 3],
    /// Last cells of the box (exclusive).
    pub hi: [usize; 3],
    /// Its voltage or current.
    pub drive: ElectrodeDrive,
}

/// Conductors and electrodes of a domain.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ElectricSetup {
    /// Per solid material index: its conductor model, or `None` (an
    /// insulator). Missing trailing entries are insulators. Ignored by the
    /// explicit per-cell constitutive path, including its topology decisions.
    pub conductors: Vec<Option<Conductor>>,
    /// Electrodes; grounded paths require a voltage reference.
    /// The opt-in isolated-motion path requires this list to be empty.
    pub electrodes: Vec<Electrode>,
}

impl ElectricSetup {
    /// Whether any conductor's resistivity depends on temperature.
    #[must_use]
    pub fn temperature_dependent(&self) -> bool {
        self.conductors
            .iter()
            .flatten()
            .any(|c| c.temperature_coefficient_per_k != 0.0)
    }

    /// Solve with a state-resolved resistivity law in the existing spatial kernel.
    ///
    /// `resolve(cell, material_index, temperature_k)` is called once per
    /// non-electrode conductor cell, in increasing cell order. The caller can
    /// resolve its material card at that cell's full state (including pressure,
    /// composition and processing coordinates captured by the closure), retain
    /// the resulting property-use receipts, and return resistivity in ohm metres.
    /// No material name, temperature coefficient, fallback value or extrapolation
    /// is inferred here. An error from the constitutive owner aborts the solve.
    ///
    /// `conductors` still declares electrical topology and valid reference data;
    /// its `Some`/`None` entries determine which solid materials conduct. On this
    /// path its reference/linear coefficients do not replace the callback result.
    /// Temperature needs one entry per domain cell and must be finite and positive
    /// at every cell where the material law is evaluated. Electrode and insulating
    /// cells do not demand an electrical material property.
    ///
    /// This is a frozen-state electrical solve, not a thermal evolution or an
    /// experimental qualification. The callback owns its evidence; numerical
    /// coefficients alone do not acquire material-card authority.
    ///
    /// # Errors
    /// Shape/temperature/resistivity errors, constitutive refusals, the ordinary
    /// electrical solver's refusals, or cancellation. No field is published on
    /// failure; externally captured callback state remains the caller's concern.
    pub fn solve_with_resistivity(
        &self,
        domain: &VoxelDomain,
        temperature: &[f64],
        mut resolve: impl FnMut(usize, u16, f64) -> Result<f64, ChtError>,
        gate: &CancelGate,
    ) -> Result<ElectricSolution, ChtError> {
        solve_electric_impl(
            domain,
            self,
            Some(temperature),
            |cell, material, point, _| {
                // The shared kernel checks the complete shape before evaluation.
                resolve(cell, material, point.expect("explicit temperature field"))
            },
            &[],
            gate,
        )
    }

    /// Solve DC conduction in material-assigned fluid AND solid cells.
    ///
    /// `resolve(cell, occupancy, temperature_k)` is called once per non-electrode
    /// cell, in increasing index order. `Some(rho)` declares a finite positive
    /// isotropic resistivity [ohm m]; `None` explicitly declares insulation.
    /// Missing/unsupported material data must return an error, not `None`.
    /// The caller owns cell-to-card assignments, full query coordinates and
    /// retained receipts. Neither phase nor a material name implies conductivity.
    /// `conductors` is ignored on this path: no dummy solid preset is required.
    ///
    /// Electrode geometry/drive validation, grounding, matrix assembly, current
    /// solve and half-resistance Joule deposition are shared with the legacy
    /// path. Ideal electrodes need no resistivity/temperature evaluation.
    /// Mechanical occupancy is unchanged: a conducting fluid still flows and
    /// advects heat in the existing thermal owner. An electrical electrode
    /// designation alone does not create a mechanical wall.
    ///
    /// This is bulk ohmic conduction at a frozen state, not electrochemistry,
    /// ion transport, plasma, phase change, motional EMF or magnetic forces.
    ///
    /// # Errors
    /// Invalid electrode data, shape or non-electrode temperatures refuse before
    /// invoking the law. Invalid/missing constitutive data, grounding/solver
    /// failures, unrepresentable derived fields and cancellation abort the solve.
    /// No partial field is returned; captured callback state remains caller-owned.
    pub fn solve_with_cell_resistivity(
        &self,
        domain: &VoxelDomain,
        temperature: &[f64],
        resolve: impl FnMut(usize, Voxel, f64) -> Result<Option<f64>, ChtError>,
        gate: &CancelGate,
    ) -> Result<ElectricSolution, ChtError> {
        self.solve_with_cell_contacts(domain, temperature, &[], resolve, gate)
    }
}

/// One electrode's solved state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElectrodeSolution {
    /// Potential, V.
    pub voltage_v: f64,
    /// Current delivered into the conductors, A.
    pub current_a: f64,
}

/// The potential and its Joule heating.
#[derive(Debug, Clone, PartialEq)]
pub struct ElectricSolution {
    /// Potential per cell, V (NaN for electrical insulators and conductors
    /// with no path to a voltage electrode on legacy paths). Isolated-motion
    /// potentials use a numerical gauge, not an absolute voltage measurement.
    pub potential_v: Vec<f64>,
    /// Joule heat per cell, W, including explicitly partitioned contact heat.
    pub joule_w: Vec<f64>,
    /// Sum of `joule_w`, W.
    pub total_joule_w: f64,
    /// Power the electrodes deliver, `sum_e V_e I_e`, W (equal to
    /// `total_joule_w` up to the solver residual in stationary solves; the
    /// moving-conductor path additionally includes Lorentz mechanical power).
    pub delivered_w: f64,
    /// One entry per electrode, in order.
    pub electrodes: Vec<ElectrodeSolution>,
    /// Krylov iterations of the potential solve.
    pub iterations: usize,
}

/// Evidence of the electro-thermal fixed point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElectrothermalReport {
    /// Potential/temperature passes.
    pub passes: usize,
    /// Largest temperature change of the last pass, K.
    pub temperature_change_k: f64,
}

/// Where a cell sits in the potential system.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Node {
    /// Carries no current.
    Off,
    /// Inside electrode `e`.
    Electrode(usize),
    /// A free conductor cell (row assigned later).
    Conductor,
}

fn admit(
    domain: &VoxelDomain,
    setup: &ElectricSetup,
    use_solid_presets: bool,
) -> Result<Vec<Node>, ChtError> {
    for conductor in setup.conductors.iter().flatten().filter(|_| use_solid_presets) {
        finite_positive("electric.resistivity_ohm_m", conductor.resistivity_ohm_m)?;
        finite(
            "electric.temperature_coefficient_per_k",
            conductor.temperature_coefficient_per_k,
        )?;
        finite(
            "electric.reference_temperature_k",
            conductor.reference_temperature_k,
        )?;
    }
    if !setup
        .electrodes
        .iter()
        .any(|e| matches!(e.drive, ElectrodeDrive::Voltage(_)))
    {
        return Err(ChtError::InvalidInput {
            field: "electric.electrodes",
            reason: "at least one electrode must hold a voltage (the potential reference)".into(),
        });
    }
    let dims = domain.dims();
    let mut nodes: Vec<Node> = (0..domain.cell_count())
        .map(|c| match domain.voxel_at(c) {
            Voxel::Solid(m)
                if use_solid_presets && setup
                    .conductors
                    .get(m as usize)
                    .copied()
                    .flatten()
                    .is_some() =>
            {
                Node::Conductor
            }
            _ => Node::Off,
        })
        .collect();
    for (e, electrode) in setup.electrodes.iter().enumerate() {
        match electrode.drive {
            ElectrodeDrive::Voltage(v) | ElectrodeDrive::Current(v) => {
                finite("electric.electrodes.drive", v)?;
            }
        }
        if (0..3).any(|a| electrode.lo[a] >= electrode.hi[a] || electrode.hi[a] > dims[a]) {
            return Err(ChtError::InvalidInput {
                field: "electric.electrodes",
                reason: format!("electrode {e} is empty or leaves the domain"),
            });
        }
        for z in electrode.lo[2]..electrode.hi[2] {
            for y in electrode.lo[1]..electrode.hi[1] {
                for x in electrode.lo[0]..electrode.hi[0] {
                    let c = domain.index(x, y, z);
                    if let Node::Electrode(other) = nodes[c] {
                        return Err(ChtError::InvalidInput {
                            field: "electric.electrodes",
                            reason: format!("electrodes {other} and {e} overlap"),
                        });
                    }
                    nodes[c] = Node::Electrode(e);
                }
            }
        }
    }
    // Two electrodes in contact would short through an infinite conductance.
    for c in 0..domain.cell_count() {
        if let Node::Electrode(e) = nodes[c] {
            for face in 0..6 {
                if let Some(n) = domain.neighbor(c, face)
                    && let Node::Electrode(other) = nodes[n]
                    && other != e
                {
                    return Err(ChtError::InvalidInput {
                        field: "electric.electrodes",
                        reason: format!("electrodes {e} and {other} touch (a short circuit)"),
                    });
                }
            }
        }
    }
    Ok(nodes)
}

/// Solve the potential with each conductor's resistivity at `temperature`
/// (its reference resistivity without one).
///
/// # Errors
/// [`ChtError::InvalidInput`] for invalid conductors or electrodes, a
/// current electrode with no path to a voltage electrode, or a resistivity
/// that is not positive at the given temperature; solver refusals;
/// [`ChtError::Cancelled`].
pub fn solve_electric(
    domain: &VoxelDomain,
    setup: &ElectricSetup,
    temperature: Option<&[f64]>,
    gate: &CancelGate,
) -> Result<ElectricSolution, ChtError> {
    solve_electric_impl(
        domain,
        setup,
        temperature,
        |_, _, t, conductor| {
            Ok(t.map_or(conductor.resistivity_ohm_m, |t| conductor.resistivity(t)))
        },
        &[],
        gate,
    )
}

#[allow(clippy::too_many_lines)] // admission, assembly, solve, post-processing
fn solve_electric_impl(
    domain: &VoxelDomain,
    setup: &ElectricSetup,
    temperature: Option<&[f64]>,
    mut resolve: impl FnMut(usize, u16, Option<f64>, Conductor) -> Result<f64, ChtError>,
    contacts: &[ElectricalContact],
    gate: &CancelGate,
) -> Result<ElectricSolution, ChtError> {
    poll(gate)?;
    let nodes = admit(domain, setup, true)?;
    if temperature.is_some_and(|t| t.len() != domain.cell_count()) {
        return Err(ChtError::InvalidInput {
            field: "electric.temperature",
            reason: "needs one temperature per cell".into(),
        });
    }
    let cells = domain.cell_count();
    // Resolve once per cell, before matrix assembly. Both entry points use
    // this same field for current balance AND resistive heat deposition.
    let mut rho = vec![f64::NAN; cells];
    for c in 0..cells {
        poll_every(gate, c)?;
        if nodes[c] == Node::Conductor
            && let Voxel::Solid(m) = domain.voxel_at(c)
        {
            let conductor = setup.conductors[m as usize].expect("conductor node");
            let point = temperature.map(|t| t[c]);
            if let Some(t) = point {
                finite_positive("electric.temperature", t)?;
            }
            let value = resolve(c, m, point, conductor)?;
            if !(value.is_finite() && value > 0.0) {
                return Err(ChtError::InvalidInput {
                    field: "electric.resistivity_ohm_m",
                    reason: format!(
                        "resistivity {value} ohm m at cell {c} is not positive (temperature {} K)",
                        point.unwrap_or(f64::NAN)
                    ),
                });
            }
            rho[c] = value;
        }
    }
    solve_electric_field(domain, setup, nodes, &rho, contacts, gate)
}

// One numerical owner for preset solids, source-resolved solids and fluids.
// All topology and coefficients have been admitted before this boundary.
#[allow(clippy::too_many_lines)] // reachability, assembly, solve and conservative heat
fn solve_electric_field(
    domain: &VoxelDomain,
    setup: &ElectricSetup,
    nodes: Vec<Node>,
    rho: &[f64],
    contacts: &[ElectricalContact],
    gate: &CancelGate,
) -> Result<ElectricSolution, ChtError> {
    solve_electric_field_with_motion(domain, setup, nodes, rho, contacts, None, gate)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn solve_electric_field_with_motion(
    domain: &VoxelDomain,
    setup: &ElectricSetup,
    nodes: Vec<Node>,
    rho: &[f64],
    contacts: &[ElectricalContact],
    motion: Option<&mut motion::MotionCoupling<'_>>,
    gate: &CancelGate,
) -> Result<ElectricSolution, ChtError> {
    poll(gate)?;
    let contact_field = contact::ContactField::new(domain, &nodes, contacts, gate)?;
    solve_electric_field_with_contact_field(domain, setup, nodes, rho, contact_field, motion, gate)
}

// A single numerical owner consumes both resistive and explicitly open faces.
// Moving-conductor gauges use this same face topology when finding components.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn solve_electric_field_with_contact_field(
    domain: &VoxelDomain,
    setup: &ElectricSetup,
    mut nodes: Vec<Node>,
    rho: &[f64],
    contact_field: contact::ContactField,
    mut motion: Option<&mut motion::MotionCoupling<'_>>,
    gate: &CancelGate,
) -> Result<ElectricSolution, ChtError> {
    poll(gate)?;
    let cells = domain.cell_count();
    let half = |c: usize, axis: usize| 0.5 * rho[c] * domain.widths(c)[axis];
    // Only isolated-motion calls supply gauges. Gauge cells remain conductors
    // in all constitutive, contact, current, heating and force operations.
    let potential_gauges = motion.as_deref().and_then(|m| m.potential_gauges());
    // Reachability from the voltage electrodes over conducting links.
    let electrodes = setup.electrodes.len();
    let mut grounded = vec![false; cells];
    let mut electrode_grounded = vec![false; electrodes];
    let mut stack: Vec<usize> = Vec::new();
    for (c, node) in nodes.iter().enumerate() {
        if potential_gauges.is_some() && *node == Node::Conductor {
            grounded[c] = true; // solved with a component gauge, not discarded
        }
        if let Node::Electrode(e) = *node
            && matches!(setup.electrodes[e].drive, ElectrodeDrive::Voltage(_))
        {
            electrode_grounded[e] = true;
            grounded[c] = true;
            stack.push(c);
        }
    }
    let mut visit_electrode = |e: usize, grounded: &mut [bool], stack: &mut Vec<usize>| {
        if !electrode_grounded[e] {
            electrode_grounded[e] = true;
            let electrode = setup.electrodes[e];
            for z in electrode.lo[2]..electrode.hi[2] {
                for y in electrode.lo[1]..electrode.hi[1] {
                    for x in electrode.lo[0]..electrode.hi[0] {
                        let c = domain.index(x, y, z);
                        grounded[c] = true;
                        stack.push(c);
                    }
                }
            }
        }
    };
    while let Some(c) = stack.pop() {
        for face in 0..6 {
            let Some(n) = domain.neighbor(c, face) else {
                continue;
            };
            if grounded[n] || contact_field.is_open(c, n) {
                continue;
            }
            match (nodes[c], nodes[n]) {
                (Node::Off, _) | (_, Node::Off) => {}
                (_, Node::Conductor) => {
                    grounded[n] = true;
                    stack.push(n);
                }
                (Node::Conductor, Node::Electrode(e)) => {
                    visit_electrode(e, &mut grounded, &mut stack);
                }
                (Node::Electrode(_), Node::Electrode(_)) => {}
            }
        }
    }
    for (e, electrode) in setup.electrodes.iter().enumerate() {
        if let ElectrodeDrive::Current(i) = electrode.drive
            && !electrode_grounded[e]
            && i != 0.0
        {
            return Err(ChtError::InvalidInput {
                field: "electric.electrodes",
                reason: format!(
                    "current electrode {e} has no conducting path to a voltage electrode"
                ),
            });
        }
    }
    for (c, node) in nodes.iter_mut().enumerate() {
        if *node == Node::Conductor && !grounded[c] {
            if motion.is_some() {
                return Err(ChtError::InvalidInput {
                    field: "electric.motion",
                    reason: format!("conducting component at cell {c} needs a voltage reference; motional currents cannot be discarded"),
                });
            }
            *node = Node::Off;
        }
    }
    // Unknowns: free conductor cells, then the grounded current electrodes.
    let mut row_of = vec![usize::MAX; cells];
    let mut rows = 0;
    for c in 0..cells {
        if nodes[c] == Node::Conductor {
            row_of[c] = rows;
            rows += 1;
        }
    }
    let mut electrode_row = vec![usize::MAX; electrodes];
    for (e, electrode) in setup.electrodes.iter().enumerate() {
        if matches!(electrode.drive, ElectrodeDrive::Current(_)) && electrode_grounded[e] {
            electrode_row[e] = rows;
            rows += 1;
        }
    }
    let mut coo = Coo::new(rows, rows);
    let mut b = vec![0.0f64; rows];
    for (e, electrode) in setup.electrodes.iter().enumerate() {
        if let ElectrodeDrive::Current(i) = electrode.drive
            && electrode_row[e] != usize::MAX
        {
            b[electrode_row[e]] = i;
        }
    }
    let mut linked = vec![false; electrodes];
    for c in 0..cells {
        if nodes[c] != Node::Conductor {
            continue;
        }
        poll_every(gate, c)?;
        let r = row_of[c];
        if potential_gauges.is_some_and(|pins| pins[c]) {
            // Symmetric row/column elimination of phi_c = 0. This replaces a
            // redundant balance equation, NOT the cell's physical resistance.
            coo.push(r, r, 1.0);
            continue;
        }
        for face in 0..6 {
            let Some(n) = domain.neighbor(c, face) else {
                continue;
            };
            if contact_field.is_open(c, n) {
                continue;
            }
            let axis = face / 2;
            let area = domain.face_area(c, axis);
            let (contact_r, _) = contact_field.at(c, n);
            match nodes[n] {
                Node::Off => {}
                Node::Conductor => {
                    let g = area / (half(c, axis) + half(n, axis) + contact_r);
                    finite_positive("electric.face_conductance_s", g)?;
                    coo.push(r, r, g);
                    if !potential_gauges.is_some_and(|pins| pins[n]) {
                        coo.push(r, row_of[n], -g);
                    }
                    if let Some(m) = motion.as_deref() {
                        b[r] -= g * m.emf(domain, c, n, face, nodes[n])?;
                    }
                }
                Node::Electrode(e) => {
                    linked[e] = true;
                    let g = area / (half(c, axis) + contact_r);
                    finite_positive("electric.face_conductance_s", g)?;
                    coo.push(r, r, g);
                    if let Some(m) = motion.as_deref() {
                        let induced = g * m.emf(domain, c, n, face, nodes[n])?;
                        b[r] -= induced;
                        if matches!(setup.electrodes[e].drive, ElectrodeDrive::Current(_)) {
                            b[electrode_row[e]] += induced;
                        }
                    }
                    match setup.electrodes[e].drive {
                        ElectrodeDrive::Voltage(v) => b[r] += g * v,
                        ElectrodeDrive::Current(_) => {
                            let re = electrode_row[e];
                            coo.push(r, re, -g);
                            coo.push(re, re, g);
                            coo.push(re, r, -g);
                        }
                    }
                }
            }
        }
    }
    for (e, &row) in electrode_row.iter().enumerate() {
        if row != usize::MAX && !linked[e] {
            return Err(ChtError::InvalidInput {
                field: "electric.electrodes",
                reason: format!("current electrode {e} touches no conductor cell"),
            });
        }
    }
    if motion.is_some() && b.iter().any(|value| !value.is_finite()) {
        return Err(ChtError::InvalidInput {
            field: "electric.motion", reason: "motional current right-hand side overflows".into(),
        });
    }
    let mut x = vec![0.0f64; rows];
    let mut iterations = 0;
    if rows > 0 && b.iter().any(|v| *v != 0.0) {
        poll(gate)?;
        let matrix = coo.assemble();
        let amg = SaAmg::new(&matrix, 0.08, 3);
        let report = pcg(&matrix, &b, &mut x, &amg, 1e-12, 2000);
        iterations = report.iters;
        if !(report.converged && x.iter().all(|v| v.is_finite())) {
            x.fill(0.0);
            let scaled = scale_rows(&coo, &mut b);
            let outcome = bicgstab_ilu0(
                "electric potential",
                &scaled,
                &b,
                &mut x,
                1e-12,
                20_000,
                gate,
            )?;
            iterations += outcome.iterations;
        }
    }
    let electrode_potential: Vec<f64> = setup
        .electrodes
        .iter()
        .enumerate()
        .map(|(e, electrode)| match electrode.drive {
            ElectrodeDrive::Voltage(v) => v,
            ElectrodeDrive::Current(_) if electrode_row[e] != usize::MAX => x[electrode_row[e]],
            ElectrodeDrive::Current(_) => f64::NAN,
        })
        .collect();
    let potential_v: Vec<f64> = (0..cells)
        .map(|c| match nodes[c] {
            Node::Conductor => x[row_of[c]],
            Node::Electrode(e) => electrode_potential[e],
            Node::Off => f64::NAN,
        })
        .collect();
    // Dissipation per face: bulk heat stays local; contact heat uses its
    // explicit surface partition. Both use the exact assembly conductance.
    let mut joule_w = vec![0.0f64; cells];
    let mut current = vec![0.0f64; electrodes];
    for c in 0..cells {
        poll_every(gate, c)?;
        if nodes[c] != Node::Conductor {
            continue;
        }
        for face in 0..6 {
            let Some(n) = domain.neighbor(c, face) else {
                continue;
            };
            if contact_field.is_open(c, n) {
                continue;
            }
            let axis = face / 2;
            let area = domain.face_area(c, axis);
            let (contact_r, fraction) = contact_field.at(c, n);
            match nodes[n] {
                Node::Off => {}
                // Each conductor pair once, from its lower cell.
                Node::Conductor if face % 2 == 1 => {
                    let (hc, hn) = (half(c, axis), half(n, axis));
                    let resistance = hc + hn + contact_r;
                    let mut drop = potential_v[c] - potential_v[n];
                    if let Some(m) = motion.as_deref_mut() {
                        drop += m.emf(domain, c, n, face, nodes[n])?;
                        m.deposit(domain, c, n, face, nodes[n], area / resistance * drop)?;
                    }
                    let heat = area / resistance * drop * drop;
                    if contact_r == 0.0 {
                        // Preserve the established perfect-contact arithmetic.
                        joule_w[c] += heat * hc / (hc + hn);
                        joule_w[n] += heat * hn / (hc + hn);
                    } else {
                        let contact_heat = heat * (contact_r / resistance);
                        joule_w[c] += heat * (hc / resistance) + fraction * contact_heat;
                        joule_w[n] += heat * (hn / resistance) + (1.0 - fraction) * contact_heat;
                    }
                }
                Node::Conductor => {}
                Node::Electrode(e) => {
                    let hc = half(c, axis);
                    let resistance = hc + contact_r;
                    let g = area / resistance;
                    let mut drop = electrode_potential[e] - potential_v[c];
                    if let Some(m) = motion.as_deref_mut() {
                        // drop is electrode -> c; the face helper is c -> electrode.
                        drop -= m.emf(domain, c, n, face, nodes[n])?;
                        m.deposit(domain, c, n, face, nodes[n], -g * drop)?;
                    }
                    let heat = g * drop * drop;
                    if contact_r == 0.0 {
                        joule_w[c] += heat;
                    } else {
                        let contact_heat = heat * (contact_r / resistance);
                        joule_w[c] += heat * (hc / resistance) + fraction * contact_heat;
                        joule_w[n] += (1.0 - fraction) * contact_heat;
                    }
                    current[e] += g * drop;
                }
            }
        }
    }
    let total_joule_w: f64 = joule_w.iter().sum();
    let delivered_w: f64 = electrode_potential
        .iter()
        .zip(&current)
        .filter(|(v, _)| v.is_finite())
        .map(|(v, i)| v * i)
        .sum();
    if joule_w.iter().any(|q| !q.is_finite() || *q < 0.0)
        || current.iter().any(|i| !i.is_finite())
        || !total_joule_w.is_finite() || !delivered_w.is_finite()
    {
        return Err(ChtError::InvalidInput {
            field: "electric.result",
            reason: "current or Joule heating is unrepresentable".into(),
        });
    }
    poll(gate)?;
    Ok(ElectricSolution {
        potential_v,
        joule_w,
        total_joule_w,
        delivered_w,
        electrodes: electrode_potential
            .iter()
            .zip(&current)
            .map(|(&voltage_v, &current_a)| ElectrodeSolution {
                voltage_v,
                current_a,
            })
            .collect(),
        iterations,
    })
}

fn poll_every(gate: &CancelGate, i: usize) -> Result<(), ChtError> {
    if i % 4096 == 0 { poll(gate) } else { Ok(()) }
}

/// Couple the potential and the temperature: Joule heat from the potential
/// at the current temperature joins `thermal`'s sources, `energy` solves
/// the temperature, and the two repeat until the largest temperature
/// change falls below `1e-9` of the temperature span (one pass when no
/// conductor depends on temperature).
///
/// # Errors
/// The refusals of [`solve_electric`] and of `energy`;
/// [`ChtError::SolverNotConverged`] (`"electrothermal"`) after 200 passes,
/// which a thermal runaway produces.
pub fn solve_electrothermal(
    domain: &VoxelDomain,
    electric: &ElectricSetup,
    thermal: &ThermalSetup,
    mut energy: impl FnMut(&ThermalSetup) -> Result<EnergySolution, ChtError>,
    gate: &CancelGate,
) -> Result<(EnergySolution, ElectricSolution, ElectrothermalReport), ChtError> {
    const PASSES: usize = 200;
    let cells = domain.cell_count();
    let coupled = electric.temperature_dependent();
    let mut heated = thermal.clone();
    let mut temperature: Option<Vec<f64>> = None;
    let mut change = f64::INFINITY;
    for pass in 1..=PASSES {
        poll(gate)?;
        let field = solve_electric(domain, electric, temperature.as_deref(), gate)?;
        heated.power_w = if thermal.power_w.is_empty() {
            field.joule_w.clone()
        } else {
            thermal
                .power_w
                .iter()
                .zip(&field.joule_w)
                .map(|(p, q)| p + q)
                .collect()
        };
        debug_assert_eq!(heated.power_w.len(), cells);
        let solution = energy(&heated)?;
        if let Some(previous) = &temperature {
            change = previous
                .iter()
                .zip(&solution.temperature)
                .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
        }
        let (lo, hi) = solution
            .temperature
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), t| {
                (lo.min(*t), hi.max(*t))
            });
        if !coupled || change <= 1e-9 * (hi - lo).max(1.0) {
            return Ok((
                solution,
                field,
                ElectrothermalReport {
                    passes: pass,
                    temperature_change_k: if coupled { change } else { 0.0 },
                },
            ));
        }
        if !(hi.is_finite() && lo.is_finite()) {
            break;
        }
        temperature = Some(solution.temperature);
    }
    Err(ChtError::SolverNotConverged {
        system: "electrothermal (thermal runaway?)",
        iterations: PASSES,
        relative_residual: change,
        tolerance: 1e-9,
    })
}

#[cfg(test)]
mod material_tests {
    use super::*;

    fn bar() -> (VoxelDomain, ElectricSetup) {
        let domain = VoxelDomain::from_fn(5, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
        let setup = ElectricSetup {
            conductors: vec![Some(Conductor::new(2.0))],
            electrodes: vec![
                Electrode {
                    lo: [0, 0, 0],
                    hi: [1, 1, 1],
                    drive: ElectrodeDrive::Voltage(0.0),
                },
                Electrode {
                    lo: [4, 0, 0],
                    hi: [5, 1, 1],
                    drive: ElectrodeDrive::Voltage(9.0),
                },
            ],
        };
        (domain, setup)
    }

    #[test]
    fn g1_g3_state_resistivity_drives_series_current_and_local_joule_heat() {
        let (domain, setup) = bar();
        let temperature = [300.0, 300.0, 310.0, 320.0, 300.0];
        let gate = CancelGate::new();
        for scale in [1.0, 2.0] {
            let mut visited = Vec::new();
            let field = setup
                .solve_with_resistivity(
                    &domain,
                    &temperature,
                    |cell, material, t| {
                        visited.push((cell, material, t));
                        Ok(scale * [2.0, 4.0, 8.0][cell - 1])
                    },
                    &gate,
                )
                .unwrap();
            assert_eq!(visited, vec![(1, 0, 300.0), (2, 0, 310.0), (3, 0, 320.0)]);
            // Three unit-length, unit-area resistors between electrode surfaces.
            let current = 9.0 / (scale * 14.0);
            assert!((field.electrodes[1].current_a - current).abs() < 1e-9);
            assert!((field.electrodes[0].current_a + current).abs() < 1e-9);
            for (cell, rho) in [(1, 2.0), (2, 4.0), (3, 8.0)] {
                assert!((field.joule_w[cell] - current * current * scale * rho).abs() < 1e-9);
            }
            assert!((field.total_joule_w - 9.0 * current).abs() < 1e-9);
            assert!((field.total_joule_w - field.delivered_w).abs() < 1e-9);
        }
    }

    #[test]
    fn g5_material_callback_and_reference_law_share_identical_kernel() {
        let (domain, setup) = bar();
        let gate = CancelGate::new();
        let temperature = [300.0; 5];
        let reference = solve_electric(&domain, &setup, Some(&temperature), &gate).unwrap();
        let resolved = setup
            .solve_with_resistivity(&domain, &temperature, |_, _, _| Ok(2.0), &gate)
            .unwrap();
        assert_eq!(reference, resolved);
    }

    #[test]
    fn g0_material_refusal_and_invalid_values_never_publish_an_electrical_field() {
        let (domain, setup) = bar();
        let gate = CancelGate::new();
        let refusal = ChtError::InvalidInput {
            field: "material.state",
            reason: "outside the selected resistivity claim's temperature domain".into(),
        };
        assert_eq!(
            setup.solve_with_resistivity(&domain, &[300.0; 5], |_, _, _| Err(refusal.clone()), &gate),
            Err(refusal.clone())
        );
        for rho in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                setup.solve_with_resistivity(&domain, &[300.0; 5], |_, _, _| Ok(rho), &gate),
                Err(ChtError::InvalidInput { field: "electric.resistivity_ohm_m", .. })
            ));
        }
        for temperature in [vec![300.0; 4], vec![300.0; 6], vec![f64::NAN; 5], vec![0.0; 5]] {
            let mut calls = 0;
            assert!(matches!(
                setup.solve_with_resistivity(&domain, &temperature, |_, _, _| {
                    calls += 1;
                    Ok(2.0)
                }, &gate),
                Err(ChtError::InvalidInput { field: "electric.temperature", .. })
            ));
            assert_eq!(calls, 0);
        }
    }

    #[test]
    fn fluid_material_law_conducts_without_relabeling_mechanical_cells() {
        let domain = VoxelDomain::new(5, 1, 1, 1.0).unwrap();
        let (_, mut setup) = bar();
        setup.conductors.clear();
        let gate = CancelGate::new();
        let old = solve_electric(&domain, &setup, None, &gate).unwrap();
        assert_eq!(old.total_joule_w, 0.0); // existing behavior remains opt-in
        let mut visited = Vec::new();
        let result = setup.solve_with_cell_resistivity(&domain,
            &[f64::NAN, 300.0, 310.0, 320.0, f64::NAN],
            |cell, occupancy, t| {
                visited.push((cell, occupancy, t));
                Ok(Some([2.0, 4.0, 8.0][cell - 1]))
            }, &gate).unwrap();
        assert_eq!(visited, vec![(1, Voxel::Fluid, 300.0), (2, Voxel::Fluid, 310.0), (3, Voxel::Fluid, 320.0)]);
        let current = 9.0 / 14.0;
        assert!((result.electrodes[1].current_a - current).abs() < 1e-9);
        assert!((result.electrodes[0].current_a + current).abs() < 1e-9);
        for (cell, rho) in [(1, 2.0), (2, 4.0), (3, 8.0)] {
            assert!((result.joule_w[cell] - current * current * rho).abs() < 1e-9);
        }
        assert!((result.delivered_w - result.total_joule_w).abs() < 1e-9);
        assert!((0..5).all(|cell| domain.voxel_at(cell) == Voxel::Fluid));
    }

    #[test]
    fn mixed_solid_liquid_graded_path_conserves_current_and_power_for_both_drives() {
        let domain = VoxelDomain::graded_from_fn(
            [vec![1.0, 0.5, 2.0, 1.5, 1.0], vec![2.0], vec![3.0]],
            |point| if point[0] < 1.5 { Voxel::Solid(7) } else { Voxel::Fluid },
        ).unwrap();
        let (_, mut setup) = bar();
        setup.conductors.clear(); // electrical topology is NOT a solid table
        let gate = CancelGate::new();
        // Interior R = sum(rho*dx)/A = (2*.5 + 4*2 + 8*1.5)/6 = 3.5 ohm.
        let resistance = 3.5;
        for drive in [ElectrodeDrive::Voltage(9.0), ElectrodeDrive::Current(-2.0)] {
            setup.electrodes[1].drive = drive;
            let result = setup.solve_with_cell_resistivity(&domain, &[300.0; 5],
                |cell, _, _| Ok(Some([2.0, 4.0, 8.0][cell - 1])), &gate).unwrap();
            let current = match drive {
                ElectrodeDrive::Voltage(v) => v / resistance,
                ElectrodeDrive::Current(i) => i,
            };
            assert!((result.electrodes[1].voltage_v - current * resistance).abs() < 1e-9);
            assert!((result.electrodes[1].current_a - current).abs() < 1e-9);
            assert!((result.electrodes[0].current_a + current).abs() < 1e-9);
            for (cell, r) in [(1, 1.0/6.0), (2, 8.0/6.0), (3, 2.0)] {
                assert!((result.joule_w[cell] - current*current*r).abs() < 1e-9);
            }
            assert!((result.total_joule_w - current*current*resistance).abs() < 1e-9);
            assert!((result.delivered_w - result.total_joule_w).abs() < 1e-9);
        }
    }

    #[test]
    fn per_cell_and_legacy_solid_paths_share_bit_identical_numerics() {
        let (domain, setup) = bar();
        let gate = CancelGate::new();
        let legacy = solve_electric(&domain, &setup, Some(&[300.0; 5]), &gate).unwrap();
        let cells = setup.solve_with_cell_resistivity(&domain, &[300.0; 5],
            |_, _, _| Ok(Some(2.0)), &gate).unwrap();
        assert_eq!(legacy, cells);
    }

    #[test]
    fn explicit_insulation_is_not_a_fallback_for_a_missing_material_claim() {
        let domain = VoxelDomain::new(5, 1, 1, 1.0).unwrap();
        let (_, mut setup) = bar();
        let gate = CancelGate::new();
        let result = setup.solve_with_cell_resistivity(&domain, &[300.0; 5],
            |cell, _, _| Ok(if cell == 2 { None } else { Some(2.0) }), &gate).unwrap();
        assert!(result.potential_v[2].is_nan());
        assert!(result.total_joule_w.abs() < 1e-20);
        let refusal = ChtError::InvalidInput { field: "material.resistivity", reason: "missing source claim".into() };
        assert_eq!(setup.solve_with_cell_resistivity(&domain, &[300.0; 5],
            |_, _, _| Err(refusal.clone()), &gate), Err(refusal));
        setup.electrodes[1].drive = ElectrodeDrive::Current(1.0);
        assert!(setup.solve_with_cell_resistivity(&domain, &[300.0; 5],
            |cell, _, _| Ok(if cell == 2 { None } else { Some(2.0) }), &gate).is_err());
    }

    #[test]
    fn per_cell_queries_preflight_temperatures_and_reject_unrepresentable_laws() {
        let domain = VoxelDomain::new(5, 1, 1, 1.0).unwrap();
        let (_, setup) = bar();
        let gate = CancelGate::new();
        for temperature in [vec![300.0; 4], vec![300.0, 300.0, 300.0, f64::NAN, 300.0]] {
            assert!(setup.solve_with_cell_resistivity(&domain, &temperature,
                |_, _, _| panic!("whole field must be admitted before evaluating material data"), &gate).is_err());
        }
        for rho in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
            assert!(setup.solve_with_cell_resistivity(&domain, &[300.0; 5],
                |_, _, _| Ok(Some(rho)), &gate).is_err());
        }
    }

}
