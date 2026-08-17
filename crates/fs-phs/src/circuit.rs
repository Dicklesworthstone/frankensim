//! Kirchhoff DAE interconnection (music bead
//! `frankensim-music-v8-root-3ez8g.9.1`): wire circuit elements into
//! graphs where KCL/KVL are ALGEBRAIC constraints of a
//! [`DescriptorPortHamiltonian`], not soft couplings. THIS MODULE
//! CONSUMES THE RECORDED DEFERRAL TRIGGER — the crate CONTRACT
//! deferred constrained Dirac structures "until the first consumer
//! needing constraints"; the electric-guitar track is that consumer,
//! by name.
//!
//! Formulation (the displacement argument for descriptor-pHS over an
//! MNA library is passivity BY CONSTRUCTION): extended state
//! `[φ_L…, q_C… | v_node…, i_C…, i_V…, i_T…]` — fluxes and charges
//! are differential and carry ALL the energy
//! (`H = Σ φ²/2L + Σ q²/2C`); node potentials and constraint
//! currents are algebraic multipliers with zero storage. Every
//! constraint enters the composite `J` ONCE, and Dirac ANTISYMMETRY
//! manufactures its dual for free: the inductor's KVL row
//! (`φ̇ = v_a − v_b`) transposes into the node rows as exactly the
//! inductor's KCL current; the capacitor's branch row
//! (`v_a − v_b − q/C = 0`) transposes into both the charge dynamics
//! (`q̇ = i_C`) and the node currents; the voltage-source row pins
//! the branch voltage to the port while its transpose draws the
//! source current; the ideal transformer's single row
//! (`v₁ − N v₂ = 0`) transposes into `±i` and `∓N i` at the two
//! branch's nodes — power-conserving by antisymmetry, never by
//! bookkeeping. Resistors are pure `R`-matrix dissipation
//! (`(v_a − v_b)²/R`, PSD stamps).
//!
//! Admission is LOUD: floating nodes, shorted or parallel voltage
//! sources, and non-physical element values refuse BY NAME before any
//! stepping; consistent initial conditions are solved by a `dt = 0`
//! descriptor step (the algebraic subsystem alone), and a stall there
//! refuses as inconsistent ICs rather than NaN-ing later.
//!
//! NO-CLAIMS (v1): no event switching — diodes/triodes arrive as
//! SMOOTH device laws in the device bead, and switching topologies
//! are a later, NAMED extension; index analysis is structural
//! (the admission catches the classic pathologies; an exotic graph
//! that stalls Newton refuses loudly rather than silently drifting).

use crate::{
    DescriptorPortHamiltonian, PhsError, QuadraticStorage, StepRecord, step_descriptor,
};

/// One circuit branch between two node indices (node 0 is GROUND).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Branch {
    /// Resistor [ohm].
    Resistor {
        /// Resistance.
        ohms: f64,
    },
    /// Inductor [henry].
    Inductor {
        /// Inductance.
        henries: f64,
    },
    /// Capacitor [farad].
    Capacitor {
        /// Capacitance.
        farads: f64,
    },
    /// Ideal voltage source driven by external port `port`.
    VoltageSource {
        /// Port index.
        port: usize,
    },
    /// Ideal current source driven by external port `port`.
    CurrentSource {
        /// Port index.
        port: usize,
    },
}

/// An ideal transformer coupling two node pairs (`v1 = ratio * v2`,
/// `i2 = -ratio * i1` — power-conserving by construction).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransformerLink {
    /// Primary terminals (a, b).
    pub primary: (usize, usize),
    /// Secondary terminals (a, b).
    pub secondary: (usize, usize),
    /// Turns ratio `v1 / v2`.
    pub ratio: f64,
}

/// A circuit graph: `node_count` nodes (0 = ground), branches, and
/// ideal transformer links.
#[derive(Debug, Clone, PartialEq)]
pub struct CircuitGraph {
    /// Total node count including ground.
    pub node_count: usize,
    /// Branches as (node_a, node_b, element).
    pub branches: Vec<(usize, usize, Branch)>,
    /// Ideal transformers.
    pub transformers: Vec<TransformerLink>,
}

/// Typed circuit refusals.
#[derive(Debug)]
pub enum CircuitError {
    /// Structural problem, by name.
    Invalid {
        /// What.
        what: &'static str,
    },
    /// A node with no incident branch (its potential would be free).
    FloatingNode {
        /// Which node.
        node: usize,
    },
    /// A voltage source with both terminals on the same node.
    ShortedSource {
        /// Branch index.
        branch: usize,
    },
    /// Two ideal voltage sources across the same node pair.
    ParallelSources {
        /// The two branch indices.
        branches: (usize, usize),
    },
    /// The `dt = 0` consistency solve stalled: the supplied initial
    /// differential state contradicts the algebraic constraints.
    InconsistentInitialConditions {
        /// Newton residual at the stall.
        residual: f64,
    },
    /// Underlying pHS admission/step refusal.
    Phs(PhsError),
}

impl From<PhsError> for CircuitError {
    fn from(e: PhsError) -> Self {
        CircuitError::Phs(e)
    }
}

/// The assembled circuit: the descriptor system plus the index maps a
/// consumer needs to read states and audit the ledger.
pub struct CircuitDae {
    /// The composed descriptor-pHS.
    pub system: DescriptorPortHamiltonian,
    /// Extended-state offsets of each inductor flux, in branch order.
    pub flux_index: Vec<usize>,
    /// Extended-state offsets of each capacitor charge, in branch order.
    pub charge_index: Vec<usize>,
    /// Extended-state offset of each non-ground node's potential.
    pub node_potential_index: Vec<usize>,
    /// Port count.
    pub ports: usize,
}

/// Assemble a circuit graph into a [`CircuitDae`].
///
/// # Errors
/// [`CircuitError`] refusals as documented on the module.
#[allow(clippy::too_many_lines)] // one assembly pass, kept whole
pub fn assemble_circuit(graph: &CircuitGraph) -> Result<CircuitDae, CircuitError> {
    if graph.node_count < 2 {
        return Err(CircuitError::Invalid {
            what: "a circuit needs ground plus at least one node",
        });
    }
    if graph.branches.is_empty() {
        return Err(CircuitError::Invalid {
            what: "a circuit needs at least one branch",
        });
    }
    let mut incident = vec![0usize; graph.node_count];
    let mut ports = 0usize;
    for (idx, &(a, b, el)) in graph.branches.iter().enumerate() {
        if a >= graph.node_count || b >= graph.node_count {
            return Err(CircuitError::Invalid {
                what: "branch terminal out of range",
            });
        }
        incident[a] += 1;
        incident[b] += 1;
        match el {
            Branch::Resistor { ohms } => {
                if !(ohms.is_finite() && ohms > 0.0) {
                    return Err(CircuitError::Invalid {
                        what: "resistance must be positive and finite",
                    });
                }
            }
            Branch::Inductor { henries } => {
                if !(henries.is_finite() && henries > 0.0) {
                    return Err(CircuitError::Invalid {
                        what: "inductance must be positive and finite",
                    });
                }
            }
            Branch::Capacitor { farads } => {
                if !(farads.is_finite() && farads > 0.0) {
                    return Err(CircuitError::Invalid {
                        what: "capacitance must be positive and finite",
                    });
                }
            }
            Branch::VoltageSource { port } | Branch::CurrentSource { port } => {
                if a == b {
                    return Err(CircuitError::ShortedSource { branch: idx });
                }
                ports = ports.max(port + 1);
            }
        }
        if a == b && !matches!(el, Branch::Resistor { .. }) {
            return Err(CircuitError::Invalid {
                what: "a storage branch cannot loop one node",
            });
        }
    }
    for t in &graph.transformers {
        for &n in &[t.primary.0, t.primary.1, t.secondary.0, t.secondary.1] {
            if n >= graph.node_count {
                return Err(CircuitError::Invalid {
                    what: "transformer terminal out of range",
                });
            }
            incident[n] += 1;
        }
        if !(t.ratio.is_finite() && t.ratio != 0.0) {
            return Err(CircuitError::Invalid {
                what: "transformer ratio must be finite and nonzero",
            });
        }
    }
    for (node, &count) in incident.iter().enumerate().skip(1) {
        if count == 0 {
            return Err(CircuitError::FloatingNode { node });
        }
    }
    // Parallel ideal voltage sources over-constrain the same pair.
    for i in 0..graph.branches.len() {
        for jdx in i + 1..graph.branches.len() {
            let (a1, b1, e1) = graph.branches[i];
            let (a2, b2, e2) = graph.branches[jdx];
            if matches!(e1, Branch::VoltageSource { .. })
                && matches!(e2, Branch::VoltageSource { .. })
                && ((a1, b1) == (a2, b2) || (a1, b1) == (b2, a2))
            {
                return Err(CircuitError::ParallelSources { branches: (i, jdx) });
            }
        }
    }

    // ---- extended-state layout ----------------------------------
    let n_l = graph
        .branches
        .iter()
        .filter(|(_, _, e)| matches!(e, Branch::Inductor { .. }))
        .count();
    let n_c = graph
        .branches
        .iter()
        .filter(|(_, _, e)| matches!(e, Branch::Capacitor { .. }))
        .count();
    let n_v = graph
        .branches
        .iter()
        .filter(|(_, _, e)| matches!(e, Branch::VoltageSource { .. }))
        .count();
    let n_nodes = graph.node_count - 1; // non-ground potentials
    let n_t = graph.transformers.len();
    let n_diff = n_l + n_c;
    let n = n_diff + n_nodes + n_c + n_v + n_t;
    // Offsets.
    let off_flux = 0;
    let off_charge = n_l;
    let off_v = n_diff;
    let off_ic = off_v + n_nodes;
    let off_iv = off_ic + n_c;
    let off_it = off_iv + n_v;
    // Node potential coordinate for node k (k >= 1).
    let vcoord = |node: usize| -> Option<usize> {
        if node == 0 { None } else { Some(off_v + node - 1) }
    };

    let mut j = vec![0.0; n * n];
    let mut r = vec![0.0; n * n];
    let mut g = vec![0.0; n * ports.max(1)];
    let mut q_storage = vec![0.0; n * n];
    let mut flux_index = Vec::new();
    let mut charge_index = Vec::new();
    let (mut l_seen, mut c_seen, mut v_seen) = (0usize, 0usize, 0usize);
    let mut set_j = |j: &mut Vec<f64>, row: usize, col: usize, val: f64| {
        j[row * n + col] += val;
        j[col * n + row] -= val;
    };
    for &(a, b, el) in &graph.branches {
        match el {
            Branch::Inductor { henries } => {
                let x = off_flux + l_seen;
                l_seen += 1;
                flux_index.push(x);
                q_storage[x * n + x] = 1.0 / henries;
                // KVL row: dφ/dt = v_a − v_b; the transpose is the
                // inductor current (∂H/∂φ) leaving/entering the nodes.
                if let Some(va) = vcoord(a) {
                    set_j(&mut j, x, va, 1.0);
                }
                if let Some(vb) = vcoord(b) {
                    set_j(&mut j, x, vb, -1.0);
                }
            }
            Branch::Capacitor { farads } => {
                let xq = off_charge + c_seen;
                let xi = off_ic + c_seen;
                c_seen += 1;
                charge_index.push(xq);
                q_storage[xq * n + xq] = 1.0 / farads;
                // dq/dt = i_C (the transpose puts −∂H/∂q into the i_C
                // row) and the branch row v_a − v_b − q/C = 0 with the
                // node transposes carrying i_C into KCL.
                set_j(&mut j, xq, xi, 1.0);
                if let Some(va) = vcoord(a) {
                    set_j(&mut j, xi, va, 1.0);
                }
                if let Some(vb) = vcoord(b) {
                    set_j(&mut j, xi, vb, -1.0);
                }
            }
            Branch::Resistor { ohms } => {
                let cond = 1.0 / ohms;
                let stamp = |r: &mut Vec<f64>, i: Option<usize>, k: Option<usize>| {
                    if let (Some(i), Some(k)) = (i, k) {
                        r[i * n + k] += if i == k { cond } else { -cond };
                    }
                };
                // PSD conductance stamp on the node-potential block.
                stamp(&mut r, vcoord(a), vcoord(a));
                stamp(&mut r, vcoord(b), vcoord(b));
                if a != b {
                    stamp(&mut r, vcoord(a), vcoord(b));
                    stamp(&mut r, vcoord(b), vcoord(a));
                }
            }
            Branch::VoltageSource { port } => {
                let xi = off_iv + v_seen;
                v_seen += 1;
                // Branch row: v_a − v_b − u = 0; the transpose draws
                // the source current i_V from the nodes; the port
                // output is −i_V (current delivered INTO the source's
                // + terminal, power-conjugate to u).
                if let Some(va) = vcoord(a) {
                    set_j(&mut j, xi, va, 1.0);
                }
                if let Some(vb) = vcoord(b) {
                    set_j(&mut j, xi, vb, -1.0);
                }
                g[xi * ports + port] -= 1.0;
            }
            Branch::CurrentSource { port } => {
                // External current injected at the nodes; the output
                // is the branch voltage (power-conjugate).
                if let Some(va) = vcoord(a) {
                    g[va * ports + port] += 1.0;
                }
                if let Some(vb) = vcoord(b) {
                    g[vb * ports + port] -= 1.0;
                }
            }
        }
    }
    for (t_idx, t) in graph.transformers.iter().enumerate() {
        let xi = off_it + t_idx;
        // Single algebraic row: (v1a − v1b) − ratio (v2a − v2b) = 0;
        // antisymmetry delivers ±i at the primary and ∓ratio·i at the
        // secondary — the ideal transformer, power-exact.
        if let Some(v) = vcoord(t.primary.0) {
            set_j(&mut j, xi, v, 1.0);
        }
        if let Some(v) = vcoord(t.primary.1) {
            set_j(&mut j, xi, v, -1.0);
        }
        if let Some(v) = vcoord(t.secondary.0) {
            set_j(&mut j, xi, v, -t.ratio);
        }
        if let Some(v) = vcoord(t.secondary.1) {
            set_j(&mut j, xi, v, t.ratio);
        }
    }

    let storage = Box::new(QuadraticStorage::new(q_storage, n)?);
    let system = DescriptorPortHamiltonian {
        n,
        n_diff,
        m: ports,
        j,
        r,
        g,
        storage,
    };
    let node_potential_index = (1..graph.node_count)
        .map(|k| off_v + k - 1)
        .collect();
    Ok(CircuitDae {
        system,
        flux_index,
        charge_index,
        node_potential_index,
        ports,
    })
}

impl CircuitDae {
    /// Solve consistent initial conditions from the given differential
    /// state (fluxes/charges kept; algebraic coordinates solved by a
    /// `dt = 0` descriptor step).
    ///
    /// # Errors
    /// [`CircuitError::InconsistentInitialConditions`] when the
    /// algebraic solve stalls; dimension refusals.
    pub fn consistent_initial_state(
        &self,
        x0: &[f64],
        u0: &[f64],
    ) -> Result<Vec<f64>, CircuitError> {
        match step_descriptor(&self.system, x0, u0, 0.0) {
            Ok(rec) => Ok(rec.x),
            Err(PhsError::NewtonStalled { residual }) => {
                Err(CircuitError::InconsistentInitialConditions { residual })
            }
            Err(e) => Err(CircuitError::Phs(e)),
        }
    }

    /// One audited step: advances the DAE and returns the record plus
    /// the SUPPLY-DEFECT residual `|ΔH + dissipated − supplied|` — the
    /// independent ledger check (balance alone is tautological under
    /// sign mutations; the supply rate is computed from the port
    /// variables and catches a smuggled sign).
    ///
    /// # Errors
    /// Step refusals.
    pub fn step_audited(
        &self,
        x0: &[f64],
        u: &[f64],
        dt: f64,
    ) -> Result<(StepRecord, f64), CircuitError> {
        let rec = step_descriptor(&self.system, x0, u, dt)?;
        let defect = (rec.delta_h + rec.dissipated - rec.supplied).abs();
        Ok((rec, defect))
    }
}

#[cfg(test)]
mod circuit_tests {
    use super::*;
    use crate::det;

    fn series_rlc(r_ohm: f64, l_h: f64, c_f: f64) -> CircuitGraph {
        // node 0 = ground; V source 0->1; L 1->2; C 2->3; R 3->0.
        CircuitGraph {
            node_count: 4,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: l_h }),
                (2, 3, Branch::Capacitor { farads: c_f }),
                (3, 0, Branch::Resistor { ohms: r_ohm }),
            ],
            transformers: vec![],
        }
    }

    #[test]
    fn kd_001_series_rlc_matches_the_analytic_impedance() {
        // Drive the series RLC with a sinusoid and compare the
        // steady-state current amplitude against |Z| = |R + j(ωL −
        // 1/ωC)| at three frequencies straddling resonance — the
        // analytic oracle, computed in-test from the same three
        // numbers, never from the stepped system.
        let (r_ohm, l_h, c_f) = (220.0f64, 0.1f64, 2.2e-7f64);
        let f_res = 1.0 / (core::f64::consts::TAU * (l_h * c_f).sqrt());
        let dae = assemble_circuit(&series_rlc(r_ohm, l_h, c_f)).expect("assemble");
        let dt = 1.0 / 480_000.0;
        for f_ratio in [0.5f64, 1.0, 1.8] {
            let f = f_ratio * f_res;
            let omega = core::f64::consts::TAU * f;
            let z = (r_ohm * r_ohm + (omega * l_h - 1.0 / (omega * c_f)).powi(2)).sqrt();
            let mut x = dae
                .consistent_initial_state(&vec![0.0; dae.system.state_dim()], &[0.0])
                .expect("ics");
            // 14 cycles: the transient (tau = 2L/R ~ 0.9 ms) must be
            // gone before measuring — the 6-cycle window read a 13%
            // transient beat at 1.8x resonance on the first run.
            let n = (14.0 / f / dt) as usize;
            let mut peak = 0.0f64;
            let mut worst_defect = 0.0f64;
            for k in 0..n {
                let u = [10.0 * det::sin(omega * k as f64 * dt)];
                let (rec, defect) = dae.step_audited(&x, &u, dt).expect("step");
                x = rec.x;
                worst_defect = worst_defect.max(defect);
                if k > 2 * n / 3 {
                    // Inductor current IS the series current.
                    let i_l = x[dae.flux_index[0]] / l_h;
                    peak = peak.max(i_l.abs());
                }
            }
            let expected = 10.0 / z;
            let rel = (peak - expected).abs() / expected;
            assert!(
                rel < 0.03,
                "series RLC at {f_ratio}x resonance: |I| {peak:.5e} vs analytic {expected:.5e} \
                 (rel {rel:.4})"
            );
            assert!(
                worst_defect < 1.0e-9,
                "supply-defect audit at {f_ratio}x: {worst_defect:.3e}"
            );
            println!(
                "{{\"suite\":\"fs-phs\",\"case\":\"kd-001-rlc\",\"f_ratio\":{f_ratio},\
                 \"i_peak\":{peak:.5e},\"i_analytic\":{expected:.5e},\"defect\":{worst_defect:.2e}}}"
            );
        }
    }

    #[test]
    fn kd_002_pure_lc_loop_steps_energy_exactly() {
        // The index-relevant case: L parallel C, no damping, no
        // source. The discrete-gradient descriptor step must hold the
        // energy to machine class over a LONG horizon (1e5 steps).
        let graph = CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Inductor { henries: 0.05 }),
                (1, 0, Branch::Capacitor { farads: 1.0e-6 }),
            ],
            transformers: vec![],
        };
        let dae = assemble_circuit(&graph).expect("assemble");
        let mut x0 = vec![0.0; dae.system.state_dim()];
        x0[dae.charge_index[0]] = 1.0e-5; // 10 uC on the cap
        let mut x = dae.consistent_initial_state(&x0, &[]).expect("ics");
        let h0 = dae.system.hamiltonian(&x);
        assert!(h0 > 0.0);
        let dt = 1.0 / 480_000.0;
        let mut worst_drift = 0.0f64;
        let mut worst_defect = 0.0f64;
        for _ in 0..100_000 {
            let (rec, defect) = dae.step_audited(&x, &[], dt).expect("step");
            x = rec.x;
            worst_defect = worst_defect.max(defect);
            worst_drift = worst_drift.max((dae.system.hamiltonian(&x) - h0).abs());
        }
        assert!(
            worst_drift < 1.0e-10 * h0.max(1e-30) * 1.0e4 || worst_drift < 1.0e-12,
            "LC energy drift {worst_drift:.3e} over 1e5 steps (H0 {h0:.3e})"
        );
        assert!(worst_defect < 1.0e-12, "defect {worst_defect:.3e}");
        println!(
            "{{\"suite\":\"fs-phs\",\"case\":\"kd-002-lc\",\"h0\":{h0:.6e},\
             \"worst_drift\":{worst_drift:.3e},\"worst_defect\":{worst_defect:.3e}}}"
        );
    }

    #[test]
    fn kd_003_transformer_reflects_the_load() {
        // V source + R1 on the primary; R2 on the secondary through an
        // ideal N:1 transformer. DC steady state: the reflected load is
        // N^2 R2, so I1 = U / (R1 + N^2 R2) — analytic; and the ideal
        // transformer must pass power EXACTLY (the supply audit is the
        // proof: the only dissipation is the two resistors).
        let n_ratio = 3.0;
        let (r1, r2) = (100.0, 50.0);
        let graph = CircuitGraph {
            node_count: 4,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Resistor { ohms: r1 }),
                (3, 0, Branch::Resistor { ohms: r2 }),
            ],
            transformers: vec![TransformerLink {
                primary: (2, 0),
                secondary: (3, 0),
                ratio: n_ratio,
            }],
        };
        let dae = assemble_circuit(&graph).expect("assemble");
        let mut x = dae
            .consistent_initial_state(&vec![0.0; dae.system.state_dim()], &[0.0])
            .expect("ics");
        let dt = 1.0e-5;
        let mut worst_defect = 0.0f64;
        let mut last_y = vec![0.0];
        for _ in 0..200 {
            let (rec, defect) = dae.step_audited(&x, &[12.0], dt).expect("step");
            x = rec.x;
            last_y = rec.y.clone();
            worst_defect = worst_defect.max(defect);
        }
        // Primary current from the STEP record's port output (the
        // step's effort carries the multipliers; the static output()
        // uses grad H alone, which is ZERO on multipliers — reading it
        // for a source current returns exactly 0, measured).
        // Sign convention, measured: y[0] = -i_V is already the current
        // the source DELIVERS into its + terminal (the first read
        // negated it again and matched the analytic value to the digit
        // with the wrong sign).
        let i1 = last_y[0];
        let expected = 12.0 / (r1 + n_ratio * n_ratio * r2);
        let rel = (i1 - expected).abs() / expected.abs();
        assert!(
            rel < 1.0e-6,
            "reflected-load current {i1:.6e} vs analytic {expected:.6e} (rel {rel:.2e})"
        );
        // Secondary node voltage = N-th fraction: v3 = -? magnitude check
        let v2 = x[dae.node_potential_index[1]];
        let v3 = x[dae.node_potential_index[2]];
        assert!(
            (v2 - n_ratio * v3).abs() < 1.0e-9 * v2.abs().max(1.0),
            "transformer voltage law v2 {v2:.4e} vs N*v3 {:.4e}",
            n_ratio * v3
        );
        assert!(worst_defect < 1.0e-9, "defect {worst_defect:.3e}");
        println!(
            "{{\"suite\":\"fs-phs\",\"case\":\"kd-003-transformer\",\"i1\":{i1:.6e},\
             \"expected\":{expected:.6e},\"defect\":{worst_defect:.2e}}}"
        );
    }

    #[test]
    fn kd_004_refusals_fire_by_name() {
        // Floating node.
        let floating = CircuitGraph {
            node_count: 3,
            branches: vec![(1, 0, Branch::Resistor { ohms: 10.0 })],
            transformers: vec![],
        };
        assert!(matches!(
            assemble_circuit(&floating),
            Err(CircuitError::FloatingNode { node: 2 })
        ));
        // Shorted source.
        let shorted = CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 1, Branch::VoltageSource { port: 0 }),
                (1, 0, Branch::Resistor { ohms: 10.0 }),
            ],
            transformers: vec![],
        };
        assert!(matches!(
            assemble_circuit(&shorted),
            Err(CircuitError::ShortedSource { branch: 0 })
        ));
        // Parallel ideal sources.
        let parallel = CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (0, 1, Branch::VoltageSource { port: 1 }),
                (1, 0, Branch::Resistor { ohms: 10.0 }),
            ],
            transformers: vec![],
        };
        assert!(matches!(
            assemble_circuit(&parallel),
            Err(CircuitError::ParallelSources { .. })
        ));
        // Non-physical element.
        let bad = CircuitGraph {
            node_count: 2,
            branches: vec![(1, 0, Branch::Inductor { henries: -1.0 })],
            transformers: vec![],
        };
        assert!(matches!(
            assemble_circuit(&bad),
            Err(CircuitError::Invalid { .. })
        ));
        // Inconsistent ICs: a V-source directly across a capacitor
        // whose charge disagrees — the dt=0 solve cannot repair a
        // differential coordinate and must refuse.
        let vc = CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 0, Branch::Capacitor { farads: 1.0e-6 }),
            ],
            transformers: vec![],
        };
        let dae = assemble_circuit(&vc).expect("assemble");
        let mut x0 = vec![0.0; dae.system.state_dim()];
        x0[dae.charge_index[0]] = 5.0e-6; // 5 V on the cap
        let clash = dae.consistent_initial_state(&x0, &[12.0]); // source says 12 V
        assert!(
            matches!(
                clash,
                Err(CircuitError::InconsistentInitialConditions { .. })
            ),
            "V-across-C with a contradicting charge must refuse ({clash:?})"
        );
        println!("{{\"suite\":\"fs-phs\",\"case\":\"kd-004-refusals\",\"verdict\":\"pass\"}}");
    }

    #[test]
    fn kd_005_supply_audit_is_not_vacuous() {
        // The independent-check falsifier: with the SUPPLY measurement
        // sign-flipped, the audit residual on a driven damped circuit
        // must be LARGE — proving the audit actually constrains the
        // ledger instead of comparing a quantity with itself.
        let dae = assemble_circuit(&series_rlc(220.0, 0.1, 2.2e-7)).expect("assemble");
        let dt = 1.0 / 480_000.0;
        let mut x = dae
            .consistent_initial_state(&vec![0.0; dae.system.state_dim()], &[0.0])
            .expect("ics");
        let omega = core::f64::consts::TAU * 1000.0;
        let mut honest = 0.0f64;
        let mut flipped = 0.0f64;
        for k in 0..2000 {
            let u = [10.0 * det::sin(omega * k as f64 * dt)];
            let (rec, defect) = dae.step_audited(&x, &u, dt).expect("step");
            x = rec.x;
            honest = honest.max(defect);
            flipped = flipped.max((rec.delta_h + rec.dissipated + rec.supplied).abs());
        }
        assert!(honest < 1.0e-9, "honest audit {honest:.3e}");
        assert!(
            flipped > 1.0e3 * honest.max(1.0e-15),
            "a sign-flipped supply must fail the audit loudly \
             ({flipped:.3e} vs honest {honest:.3e})"
        );
        println!(
            "{{\"suite\":\"fs-phs\",\"case\":\"kd-005-audit\",\"honest\":{honest:.3e},\
             \"flipped\":{flipped:.3e}}}"
        );
    }
}
