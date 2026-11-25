//! Integration tests for the static structural planner (bead `.5.7`):
//! - Plan determinism (`plan_hash` bitwise-identical across independent builds
//!   and invariant under stepping; changes on topology, class, or solver policy).
//! - Block classification (`ExactLinear`, `Smooth`, `Monotone1d`, `NonSmooth`, `Mixed`).
//! - Algebraic loops, tearing variables, follower substitution order, and
//!   load-line Woodbury effort tearing (`LoadLineWoodbury`) on a drum-like
//!   modal/mass contact graph.
//! - Parallel execution waves over independent subgraphs.
//! - Time-scale grouping from component `limit_dt_s` probes.
//! - `explain()`, canonical JSON export, and `CausalizationAudit`.

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::{
    Class, Component, EndOfStep, LoadLine, PortRole, PortSpec, PortTangent, Proposal, Receipt,
};
use fs_ports::{PortKind, StableId};
use fs_system::elements::{
    Damper, Force, HeatCapacity, Mass, Motion, Slider, Slip, SourceAccount, Spring, Waveform,
};
use fs_system::{
    AitkenConfig, BlockClass, BlockSolveMode, BlockSolverPolicy, CausalDetermination,
    EXECUTION_PLAN_SCHEMA, System, SystemBuilder, TearStrategy, TimeScaleTier,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 57,
            kernel_id: 7,
            tile: 0,
            iteration: 0,
        };
        f(&Cx::new(
            &gate,
            arena,
            key,
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

/// Minimal monotone scalar component (`Class::Monotone1d`) with a configurable
/// `limit_dt_s` bound to exercise block classification and time-scale grouping.
struct MonotoneResistive {
    ports: [PortSpec; 1],
    conductance: f64,
    dt_cap_s: f64,
}

impl MonotoneResistive {
    fn new(conductance: f64, dt_cap_s: f64) -> Self {
        Self {
            ports: [PortSpec::power(
                StableId::new("p").unwrap(),
                PortKind::MechanicalForceVelocity,
                PortRole::EffortIn,
                NonZeroUsize::MIN,
            )
            .unwrap()],
            conductance,
            dt_cap_s,
        }
    }
}

#[derive(Debug)]
struct MonotoneErr;

impl Component for MonotoneResistive {
    type State = ();
    type Error = MonotoneErr;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Monotone1d
    }

    fn stored_energy_j(&self, (): &()) -> Result<f64, MonotoneErr> {
        Ok(0.0)
    }

    fn limit_dt_s(&self, (): &(), proposed_s: f64) -> Result<f64, MonotoneErr> {
        Ok(proposed_s.min(self.dt_cap_s))
    }

    fn propose(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        (): &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, MonotoneErr> {
        cx.checkpoint().map_err(|_| MonotoneErr)?;
        let e = inputs[0];
        let f = self.conductance * e;
        outputs[0] = f;
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j: dt_s * e * f,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: 1.0e-12,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

fn build_damped_chain(n: usize, damped: bool) -> System {
    let mut b = SystemBuilder::new();
    let masses: Vec<_> = (0..n)
        .map(|i| {
            b.add_with_tangent(
                &format!("mass{i}"),
                Mass::new(0.4).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s: 0.0,
                },
            )
            .unwrap()
        })
        .collect();
    let push = b
        .add_with_tangent(
            "push",
            Force::new(Waveform::constant(2.0).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();
    let mut ends: Vec<Vec<_>> = vec![Vec::new(); n];
    ends[0].push(b.port(push, "p").unwrap());
    let mut heats = Vec::new();
    for i in 0..n {
        let spring = b
            .add_with_tangent(&format!("spring{i}"), Spring::new(900.0, 0.0).unwrap(), 0.0)
            .unwrap();
        let mut between = vec![spring];
        if damped {
            let damper = b
                .add_with_tangent(&format!("damper{i}"), Damper::new(1.5).unwrap(), ())
                .unwrap();
            between.push(damper);
            heats.push(b.port(damper, "heat").unwrap());
        }
        for elem in between {
            ends[i].push(b.port(elem, "a").unwrap());
            let far = b.port(elem, "b").unwrap();
            if i + 1 < n {
                ends[i + 1].push(far);
            } else {
                b.boundary(far, Box::new(|_, _, v| v[0] = 0.0)).unwrap();
            }
        }
    }
    for (mass, members) in masses.iter().zip(&ends) {
        b.common_flow(b.port(*mass, "p").unwrap(), members).unwrap();
    }
    if damped {
        let body = b
            .add_with_tangent("body", HeatCapacity::new(50.0).unwrap(), 300.0)
            .unwrap();
        b.common_effort(b.port(body, "heat").unwrap(), &heats)
            .unwrap();
    }
    let mut sys = b.seal().unwrap();
    for h in heats {
        assert!(sys.guess(h, &[300.0]));
    }
    sys
}

#[test]
fn plan_hash_is_deterministic_invariant_across_steps_and_sensitive_to_graph_and_policy() {
    let mut sys_a = build_damped_chain(3, false);
    let sys_b = build_damped_chain(3, false);

    let plan_a0 = sys_a.execution_plan();
    let plan_b = sys_b.execution_plan();
    assert_eq!(plan_a0.plan_hash, plan_b.plan_hash);
    assert_eq!(plan_a0.plan_hash_hex(), plan_b.plan_hash_hex());
    assert_eq!(plan_a0.to_canonical_json(), plan_b.to_canonical_json());
    assert_eq!(plan_a0.explain(), plan_b.explain());
    assert!(plan_a0.static_zero_alloc_workspaces);

    // Stepping must NOT perturb the static structural plan hash.
    with_cx(|cx| {
        for _ in 0..5 {
            sys_a.step(1.0e-3, cx).unwrap();
        }
    });
    let plan_a1 = sys_a.execution_plan();
    assert_eq!(plan_a0.plan_hash, plan_a1.plan_hash);

    // Changing solver policy must produce a distinct plan hash.
    sys_a.set_solver_policy(BlockSolverPolicy::FixedPointAitken(AitkenConfig::default()));
    let plan_fp = sys_a.execution_plan();
    assert_ne!(plan_a0.plan_hash, plan_fp.plan_hash);
    assert_eq!(
        plan_fp.blocks[0].solve_mode,
        BlockSolveMode::FixedPointAitken
    );
    assert_eq!(
        plan_fp.blocks[0].tear_strategy,
        TearStrategy::FixedPointAccelerated
    );

    // Changing the graph (adding dampers -> Smooth instead of ExactLinear) changes plan_hash.
    let damped = build_damped_chain(3, true);
    let plan_damped = damped.execution_plan();
    assert_ne!(plan_a0.plan_hash, plan_damped.plan_hash);

    // ExactLinear affine chain is cached; damped chain is Smooth and uncached.
    assert_eq!(plan_a0.block_count, 1);
    assert_eq!(plan_a0.loop_count, 1);
    assert_eq!(plan_a0.blocks[0].class, BlockClass::ExactLinear);
    assert!(plan_a0.blocks[0].affine_cached);
    assert_eq!(
        plan_a0.blocks[0].tear_strategy,
        TearStrategy::EliminationCondensed
    );
    assert_eq!(plan_a0.blocks[0].total_unknowns, 9);
    assert_eq!(plan_a0.blocks[0].kept_unknowns, 3);
    assert_eq!(plan_a0.blocks[0].substituted_count, 6);
    assert_eq!(
        plan_a0.blocks[0].tear_variables,
        vec!["mass0.p", "mass1.p", "mass2.p"]
    );
    assert_eq!(plan_a0.blocks[0].cost.estimated_dense_flops_per_factor, 0);

    assert_eq!(plan_damped.blocks[0].class, BlockClass::Smooth);
    assert!(!plan_damped.blocks[0].affine_cached);
    assert_eq!(plan_damped.blocks[0].total_unknowns, 18);
    assert_eq!(plan_damped.blocks[0].kept_unknowns, 4);
    assert_eq!(plan_damped.blocks[0].tear_variables.len(), 4);
    assert!(
        plan_damped.blocks[0]
            .tear_variables
            .contains(&"body.heat".to_owned())
    );
}

/// Exact-linear mass that publishes its end-of-step velocity on a companion
/// signal output (the drum-strike `ModalStructure` / `ExactPortHamiltonian`
/// contract), used to verify Woodbury load-line tearing on 6-unknown contacts.
struct ToldMass {
    mass_kg: f64,
    ports: [PortSpec; 2],
}

impl ToldMass {
    fn new(mass_kg: f64) -> Self {
        Self {
            mass_kg,
            ports: [
                PortSpec::power(
                    StableId::new("p").unwrap(),
                    PortKind::MechanicalForceVelocity,
                    PortRole::EffortIn,
                    NonZeroUsize::MIN,
                )
                .unwrap(),
                PortSpec::signal_out(StableId::new("end").unwrap(), NonZeroUsize::MIN),
            ],
        }
    }
}

impl Component for ToldMass {
    type State = f64;
    type Error = MonotoneErr;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn end_of_step(&self, port: usize) -> EndOfStep {
        if port == 0 {
            EndOfStep::Told(1)
        } else {
            EndOfStep::Midpoint
        }
    }

    fn stored_energy_j(&self, v: &f64) -> Result<f64, MonotoneErr> {
        Ok(0.5 * self.mass_kg * v * v)
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, MonotoneErr> {
        cx.checkpoint().map_err(|_| MonotoneErr)?;
        let f = inputs[0];
        let v_new = old + dt_s * f / self.mass_kg;
        let v_mid = 0.5 * (old + v_new);
        *new = v_new;
        outputs[0] = v_mid;
        outputs[1] = v_new;
        let d_e = 0.5 * self.mass_kg * (v_new * v_new - old * old);
        let scale = 0.5 * self.mass_kg * (old * old + v_new * v_new) + (dt_s * f * v_mid).abs();
        Ok(Proposal::Step(Receipt::lossless(
            d_e,
            64.0 * f64::EPSILON * scale.max(1.0),
        )))
    }
}

impl PortTangent for ToldMass {
    fn jvp(
        &self,
        _: &f64,
        dt_s: f64,
        _: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), MonotoneErr> {
        cx.checkpoint().map_err(|_| MonotoneErr)?;
        let dv = dt_s * d_inputs[0] / self.mass_kg;
        d_outputs[0] = 0.5 * dv;
        d_outputs[1] = dv;
        Ok(())
    }

    fn vjp(
        &self,
        _: &f64,
        dt_s: f64,
        _: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), MonotoneErr> {
        cx.checkpoint().map_err(|_| MonotoneErr)?;
        bar_inputs[0] = (0.5 * bar_outputs[0] + bar_outputs[1]) * dt_s / self.mass_kg;
        Ok(())
    }
}

/// Two-sided non-smooth contact on a load line with `line_ends` (exact shape of
/// `ContactInterface` in the drum-strike graph: 2 flow ports + 2 end-of-step
/// signal inputs -> 6 loop unknowns torn to 2 active effort outputs).
struct DrumContact {
    ports: [PortSpec; 4],
    stiffness: f64,
}

impl DrumContact {
    fn new(stiffness: f64) -> Self {
        Self {
            ports: [
                PortSpec::power(
                    StableId::new("a").unwrap(),
                    PortKind::MechanicalForceVelocity,
                    PortRole::FlowIn,
                    NonZeroUsize::MIN,
                )
                .unwrap(),
                PortSpec::power(
                    StableId::new("b").unwrap(),
                    PortKind::MechanicalForceVelocity,
                    PortRole::FlowIn,
                    NonZeroUsize::MIN,
                )
                .unwrap(),
                PortSpec::signal_in(StableId::new("end_a").unwrap(), NonZeroUsize::MIN),
                PortSpec::signal_in(StableId::new("end_b").unwrap(), NonZeroUsize::MIN),
            ],
            stiffness,
        }
    }
}

impl Component for DrumContact {
    type State = f64;
    type Error = MonotoneErr;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, gap: &f64) -> Result<f64, MonotoneErr> {
        let pen = (-*gap).max(0.0);
        Ok(0.5 * self.stiffness * pen * pen)
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, MonotoneErr> {
        self.propose_on_line(old, dt_s, inputs, &[0.0; 8], new, outputs, cx)
    }
}

impl LoadLine for DrumContact {
    fn line_ports(&self) -> &[usize] {
        &[0, 1]
    }

    fn line_ends(&self) -> &[Option<usize>] {
        &[Some(2), Some(3)]
    }

    fn propose_on_line(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, MonotoneErr> {
        cx.checkpoint().map_err(|_| MonotoneErr)?;
        // Released relative closing rate: v_rel_free = inputs[0] - inputs[1].
        // Effective scalar compliance for equal-and-opposite reaction R:
        // c_eff = C00 - C01 - C10 + C11.
        let c_eff = (compliance[0] - compliance[1] - compliance[2] + compliance[3]).max(0.0);
        let v_rel_free = inputs[0] - inputs[1];
        let gap_free = *old + dt_s * v_rel_free;
        // Linear unilateral spring reaction solved directly on the compliance line.
        let reaction = if gap_free < 0.0 {
            (-self.stiffness * (*old + 0.5 * dt_s * v_rel_free))
                / (1.0 + 0.5 * self.stiffness * dt_s * c_eff)
        } else {
            0.0
        }
        .max(0.0);
        outputs[0] = reaction;
        outputs[1] = -reaction;
        let flow_a = inputs[0] - (compliance[0] * reaction - compliance[1] * reaction);
        let flow_b = inputs[1] - (compliance[2] * reaction - compliance[3] * reaction);
        let new_gap = *old + dt_s * (flow_a - flow_b);
        *new = new_gap;
        let e0 = self.stored_energy_j(old)?;
        let e1 = self.stored_energy_j(&new_gap)?;
        let work = dt_s * (reaction * flow_a - reaction * flow_b);
        Ok(Proposal::Step(Receipt {
            stored_change_j: e1 - e0,
            dissipated_j: 0.0,
            numerical_loss_j: work - (e1 - e0),
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: 1.0e-9 * (e0.abs() + e1.abs() + work.abs() + 1.0),
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

#[test]
fn load_line_contact_compiles_to_woodbury_effort_tear_and_reports_exact_tear_set() {
    // 1. Plain 4-unknown Slider: Condensed already achieves 2 kept unknowns (optimal),
    //    with affine_line_cached = true and LoadLineNewton solve mode.
    let mut b_slider = SystemBuilder::new();
    let ma = b_slider
        .add_with_tangent(
            "a",
            Mass::new(0.04).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.8,
            },
        )
        .unwrap();
    let mb = b_slider
        .add_with_tangent(
            "b",
            Mass::new(0.25).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let slider = b_slider
        .add_on_line(
            "slider",
            Slider::new(20.0, 0.4, 0.25).unwrap(),
            Slip::sliding(-0.8),
        )
        .unwrap();
    b_slider
        .bond(
            b_slider.port(ma, "p").unwrap(),
            b_slider.port(slider, "a").unwrap(),
        )
        .unwrap();
    b_slider
        .bond(
            b_slider.port(mb, "p").unwrap(),
            b_slider.port(slider, "b").unwrap(),
        )
        .unwrap();
    b_slider
        .boundary(
            b_slider.port(slider, "heat").unwrap(),
            Box::new(|_, _, v| v[0] = 293.15),
        )
        .unwrap();
    let slider_plan = b_slider.seal().unwrap().execution_plan();
    assert_eq!(slider_plan.blocks[0].class, BlockClass::NonSmooth);
    assert_eq!(
        slider_plan.blocks[0].solve_mode,
        BlockSolveMode::LoadLineNewton
    );
    assert!(slider_plan.blocks[0].affine_line_cached);
    assert_eq!(slider_plan.blocks[0].total_unknowns, 4);
    assert_eq!(slider_plan.blocks[0].kept_unknowns, 2);

    // 2. Drum-strike shape: 2 ToldMass neighbours + 2-port Contact with line_ends
    //    => 6 loop unknowns (2 efforts + 2 flows + 2 end-of-step flows), where
    //    plain Condensed would keep 4 and Woodbury LoadLineWoodbury tears down to 2!
    let mut b = SystemBuilder::new();
    let stick = b
        .add_with_tangent("stick", ToldMass::new(0.04), 0.8)
        .unwrap();
    let head = b
        .add_with_tangent("head", ToldMass::new(0.25), 0.0)
        .unwrap();
    let contact = b
        .add_on_line("contact", DrumContact::new(1.0e5), 1.0e-4)
        .unwrap();
    b.bond(b.port(stick, "p").unwrap(), b.port(contact, "a").unwrap())
        .unwrap();
    b.bond(b.port(head, "p").unwrap(), b.port(contact, "b").unwrap())
        .unwrap();

    let sys = b.seal().unwrap();
    let plan = sys.execution_plan();
    assert_eq!(plan.block_count, 1);
    assert_eq!(plan.loop_count, 1);
    let blk = &plan.blocks[0];
    assert_eq!(blk.class, BlockClass::NonSmooth);
    assert_eq!(blk.solve_mode, BlockSolveMode::LoadLineNewton);
    assert_eq!(blk.tear_strategy, TearStrategy::LoadLineWoodbury);
    assert!(blk.affine_line_cached);
    assert_eq!(blk.load_line_node.as_deref(), Some("contact"));
    assert_eq!(blk.total_unknowns, 6);
    assert_eq!(blk.kept_unknowns, 2);
    assert_eq!(blk.tear_variables, vec!["contact.a", "contact.b"]);

    let explanation = plan.explain();
    assert!(explanation.contains(EXECUTION_PLAN_SCHEMA));
    assert!(explanation.contains("load_line_woodbury"));
    assert!(explanation.contains("load_line_node=contact"));

    let audit = sys.causalization_audit();
    assert_eq!(audit.determination, CausalDetermination::WellDetermined);
    assert_eq!(audit.unmatched_equations, 0);
    assert_eq!(audit.unmatched_variables, 0);
    assert_eq!(audit.algebraic_loops.len(), 1);
    assert_eq!(audit.algebraic_loops[0].unknown_count, 6);
    assert_eq!(audit.algebraic_loops[0].kept_count, 2);
    assert_eq!(
        audit.algebraic_loops[0].tear_variables,
        vec!["contact.a", "contact.b"]
    );

    let audit_json = audit.to_canonical_json();
    assert!(audit_json.contains("\"determination\":\"well_determined\""));
    assert!(audit_json.contains("\"tear_strategy\":\"load_line_woodbury\""));
}

#[test]
fn independent_subgraphs_group_into_parallel_waves_and_timescale_tiers() {
    let mut b = SystemBuilder::new();
    // Independent fast monotone block (dt cap = 1e-4 s).
    let fast_r = b
        .add("fast_r", MonotoneResistive::new(2.0, 1.0e-4), ())
        .unwrap();
    b.boundary(b.port(fast_r, "p").unwrap(), Box::new(|_, _, v| v[0] = 1.0))
        .unwrap();

    // Independent medium monotone block (dt cap = 2e-2 s).
    let med_r = b
        .add("med_r", MonotoneResistive::new(1.0, 2.0e-2), ())
        .unwrap();
    b.boundary(b.port(med_r, "p").unwrap(), Box::new(|_, _, v| v[0] = 1.0))
        .unwrap();

    // Independent slow/unconstrained ExactLinear mass-spring loop feeding a downstream heat sink.
    let mass = b
        .add_with_tangent(
            "mass",
            Mass::new(1.0).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.5,
            },
        )
        .unwrap();
    let spring = b
        .add_with_tangent("spring", Spring::new(10.0, 0.0).unwrap(), 0.0)
        .unwrap();
    b.bond(b.port(mass, "p").unwrap(), b.port(spring, "a").unwrap())
        .unwrap();
    b.boundary(b.port(spring, "b").unwrap(), Box::new(|_, _, v| v[0] = 0.0))
        .unwrap();

    let sys = b.seal().unwrap();
    let plan = sys.execution_plan_for_dt(1.0);
    assert_eq!(plan.block_count, 3);
    assert_eq!(plan.loop_count, 1);
    // All 3 subgraphs are independent -> single parallel wave 0!
    assert_eq!(plan.waves.len(), 1);
    assert_eq!(plan.waves[0].len(), 3);

    let fast_blk = plan
        .blocks
        .iter()
        .find(|blk| blk.nodes == ["fast_r"])
        .unwrap();
    let med_blk = plan
        .blocks
        .iter()
        .find(|blk| blk.nodes == ["med_r"])
        .unwrap();
    let mech_blk = plan
        .blocks
        .iter()
        .find(|blk| blk.nodes.contains(&"mass".to_owned()))
        .unwrap();
    assert_eq!(fast_blk.class, BlockClass::Monotone1d);
    assert_eq!(fast_blk.tear_strategy, TearStrategy::Acyclic);
    assert_eq!(fast_blk.time_scale_tier, TimeScaleTier::Fast);
    assert_eq!(med_blk.class, BlockClass::Monotone1d);
    assert_eq!(med_blk.time_scale_tier, TimeScaleTier::Medium);
    assert_eq!(mech_blk.class, BlockClass::ExactLinear);
    assert_eq!(mech_blk.time_scale_tier, TimeScaleTier::Unconstrained);

    assert_eq!(plan.time_scale_groups.len(), 3);
    assert_eq!(plan.time_scale_groups[0].tier, TimeScaleTier::Fast);
    assert!(plan.time_scale_groups[0].stiffness_ratio >= 9_999.0);
    assert_eq!(plan.time_scale_groups[1].tier, TimeScaleTier::Medium);
    assert_eq!(plan.time_scale_groups[2].tier, TimeScaleTier::Unconstrained);

    let json = plan.to_canonical_json();
    assert!(json.starts_with('{') && json.ends_with('}'));
    assert!(json.contains("\"waves\":[[0,1,2]]"));
    assert!(json.contains("\"tier\":\"fast\""));

    // Regression: when every block is constrained below probe_dt (1e-2 s and
    // 1e-4 s under probe_dt = 1.0 s), relative stiffness is measured against
    // the slowest block (1e-2 s), not probe_dt (1.0 s), so the slower block is
    // Slow (1x) and the faster block is Fast (100x).
    let mut b2 = SystemBuilder::new();
    let slow_r = b2
        .add("slow_r", MonotoneResistive::new(1.0, 1.0e-2), ())
        .unwrap();
    b2.boundary(
        b2.port(slow_r, "p").unwrap(),
        Box::new(|_, _, v| v[0] = 1.0),
    )
    .unwrap();
    let fast_r2 = b2
        .add("fast_r2", MonotoneResistive::new(1.0, 1.0e-4), ())
        .unwrap();
    b2.boundary(
        b2.port(fast_r2, "p").unwrap(),
        Box::new(|_, _, v| v[0] = 1.0),
    )
    .unwrap();
    let plan2 = b2.seal().unwrap().execution_plan_for_dt(1.0);
    let b_slow = plan2
        .blocks
        .iter()
        .find(|blk| blk.nodes == ["slow_r"])
        .unwrap();
    let b_fast = plan2
        .blocks
        .iter()
        .find(|blk| blk.nodes == ["fast_r2"])
        .unwrap();
    assert_eq!(b_slow.time_scale_tier, TimeScaleTier::Slow);
    assert_eq!(b_fast.time_scale_tier, TimeScaleTier::Fast);
}
