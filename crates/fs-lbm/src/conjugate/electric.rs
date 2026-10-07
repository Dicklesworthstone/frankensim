//! DC electric conduction and Joule heating on the voxel domain.
//!
//! Solid materials may conduct electricity ([`Conductor`]: a resistivity
//! with a linear temperature coefficient). Electrodes are boxes of cells
//! treated as perfect conductors, one potential each, held at a voltage or
//! carrying a declared total current. The potential solves the discrete
//! current balance `sum_f G_f (phi_P - phi_N) = 0` on every conductor cell,
//! with `G_f = A_f / (rho_P w_P / 2 + rho_N w_N / 2)` between conductor
//! cells and `A_f / (rho_P w_P / 2)` to an electrode's surface; every other
//! face (fluid, insulating solids, the domain boundary) carries no current.
//! Each face dissipates `G_f (phi_P - phi_N)^2`, split between its two
//! cells in proportion to their half resistances (an electrode's half is
//! zero), so the Joule heat equals the power the electrodes deliver,
//! `sum_e V_e I_e` (the discrete Tellegen theorem), to solver precision.
//!
//! [`solve_electrothermal`] couples the two fields: resistivity at each
//! cell's temperature, the potential, the Joule heat as a source, the
//! temperature, repeated to a fixed point (one pass when no conductor
//! depends on temperature).
//!
//! No-claims: DC only (no skin effect, inductance or capacitance); ohmic,
//! isotropic conductors; no electrical contact resistance between
//! conductors; electrodes are ideal equipotentials; the Picard iteration
//! does not converge past a thermal runaway (current-driven conductors
//! whose resistivity rises with temperature can have no steady state, and
//! that refusal is the answer).

use fs_exec::CancelGate;
use fs_sparse::Coo;
use fs_sparse::precond::{SaAmg, pcg};

use super::domain::{Voxel, VoxelDomain};
use super::energy::{EnergySolution, ThermalSetup};
use super::flow::scale_rows;
use super::krylov::bicgstab_ilu0;
use super::{ChtError, finite, finite_positive, poll};

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
    /// Delivers this total current into the conductors it touches, A
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
    /// insulator). Missing trailing entries are insulators.
    pub conductors: Vec<Option<Conductor>>,
    /// Electrodes; at least one must hold a voltage.
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
    /// Potential per cell, V (NaN where no current can flow: fluid,
    /// insulators, conductors with no path to a voltage electrode).
    pub potential_v: Vec<f64>,
    /// Joule heat per cell, W.
    pub joule_w: Vec<f64>,
    /// Sum of `joule_w`, W.
    pub total_joule_w: f64,
    /// Power the electrodes deliver, `sum_e V_e I_e`, W (equal to
    /// `total_joule_w` up to the solver residual).
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

fn admit(domain: &VoxelDomain, setup: &ElectricSetup) -> Result<Vec<Node>, ChtError> {
    for conductor in setup.conductors.iter().flatten() {
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
                if setup
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
    // Two electrodes in contact would short through an infinite
    // conductance.
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
#[allow(clippy::too_many_lines)] // admission, assembly, solve, post-processing
pub fn solve_electric(
    domain: &VoxelDomain,
    setup: &ElectricSetup,
    temperature: Option<&[f64]>,
    gate: &CancelGate,
) -> Result<ElectricSolution, ChtError> {
    let mut nodes = admit(domain, setup)?;
    if temperature.is_some_and(|t| t.len() != domain.cell_count()) {
        return Err(ChtError::InvalidInput {
            field: "electric.temperature",
            reason: "needs one temperature per cell".into(),
        });
    }
    let cells = domain.cell_count();
    // Resistivity per conductor cell.
    let mut rho = vec![f64::NAN; cells];
    for c in 0..cells {
        if nodes[c] == Node::Conductor
            && let Voxel::Solid(m) = domain.voxel_at(c)
        {
            let conductor = setup.conductors[m as usize].expect("conductor node");
            let value =
                temperature.map_or(conductor.resistivity_ohm_m, |t| conductor.resistivity(t[c]));
            if !(value.is_finite() && value > 0.0) {
                return Err(ChtError::InvalidInput {
                    field: "electric.resistivity_ohm_m",
                    reason: format!(
                        "resistivity {value} ohm m at cell {c} is not positive (temperature {} K)",
                        temperature.map_or(f64::NAN, |t| t[c])
                    ),
                });
            }
            rho[c] = value;
        }
    }
    let half = |c: usize, axis: usize| 0.5 * rho[c] * domain.widths(c)[axis];
    // Reachability from the voltage electrodes over conducting links.
    let electrodes = setup.electrodes.len();
    let mut grounded = vec![false; cells];
    let mut electrode_grounded = vec![false; electrodes];
    let mut stack: Vec<usize> = Vec::new();
    for (c, node) in nodes.iter().enumerate() {
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
            if grounded[n] {
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
        for face in 0..6 {
            let Some(n) = domain.neighbor(c, face) else {
                continue;
            };
            let axis = face / 2;
            let area = domain.face_area(c, axis);
            match nodes[n] {
                Node::Off => {}
                Node::Conductor => {
                    let g = area / (half(c, axis) + half(n, axis));
                    coo.push(r, r, g);
                    coo.push(r, row_of[n], -g);
                }
                Node::Electrode(e) => {
                    linked[e] = true;
                    let g = area / half(c, axis);
                    coo.push(r, r, g);
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
    // Dissipation per face, split by half resistance; electrode currents.
    let mut joule_w = vec![0.0f64; cells];
    let mut current = vec![0.0f64; electrodes];
    for c in 0..cells {
        if nodes[c] != Node::Conductor {
            continue;
        }
        for face in 0..6 {
            let Some(n) = domain.neighbor(c, face) else {
                continue;
            };
            let axis = face / 2;
            let area = domain.face_area(c, axis);
            match nodes[n] {
                Node::Off => {}
                // Each conductor pair once, from its lower cell.
                Node::Conductor if face % 2 == 1 => {
                    let (hc, hn) = (half(c, axis), half(n, axis));
                    let drop = potential_v[c] - potential_v[n];
                    let heat = area / (hc + hn) * drop * drop;
                    joule_w[c] += heat * hc / (hc + hn);
                    joule_w[n] += heat * hn / (hc + hn);
                }
                Node::Conductor => {}
                Node::Electrode(e) => {
                    let g = area / half(c, axis);
                    let drop = electrode_potential[e] - potential_v[c];
                    joule_w[c] += g * drop * drop;
                    current[e] += g * drop;
                }
            }
        }
    }
    let total_joule_w = joule_w.iter().sum();
    let delivered_w = electrode_potential
        .iter()
        .zip(&current)
        .filter(|(v, _)| v.is_finite())
        .map(|(v, i)| v * i)
        .sum();
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
