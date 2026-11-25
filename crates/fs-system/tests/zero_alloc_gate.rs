//! Counting-allocator zero-allocation gate across gallery assemblies, solve
//! modes, and the conformance-kit `check_zero_alloc_propose` check (bead `.15.1`).
//!
//! Runs in its own test binary so the thread-local counting `GlobalAlloc`
//! measures only this thread's accepted-step hot path.

use core::num::NonZeroUsize;
use std::alloc::{GlobalAlloc, Layout, System as SystemAllocator};
use std::cell::Cell;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Failure, check_zero_alloc_propose};
use fs_ports::component::{
    Class, Component, EndOfStep, LoadLine, PortRole, PortSpec, PortTangent, Proposal, Receipt,
};
use fs_ports::{PortKind, StableId};
use fs_system::elements::{
    Damper, Force, HeatCapacity, Mass, Motion, Slider, Slip, SourceAccount, Spring, Waveform,
};
use fs_system::{
    AitkenConfig, BlockSolverPolicy, IqnIlsConfig, RouteTemperature, System, SystemBuilder,
};

struct Counting;

thread_local! {
    static TRACKING: Cell<bool> = const { Cell::new(false) };
    static ALLOCATIONS: Cell<usize> = const { Cell::new(0) };
}

fn record_alloc() {
    if TRACKING.try_with(Cell::get).unwrap_or(false) {
        let _ = ALLOCATIONS.try_with(|c| c.set(c.get() + 1));
    }
}

fn allocations_of(f: impl FnOnce()) -> usize {
    ALLOCATIONS.with(|c| c.set(0));
    TRACKING.with(|t| t.set(true));
    f();
    TRACKING.with(|t| t.set(false));
    ALLOCATIONS.with(Cell::get)
}

// SAFETY: pure forwarding wrapper over `std::alloc::System` with lock-free
// thread-local counting that never allocates or inspects delegated memory.
#[allow(unsafe_code)]
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        record_alloc();
        // SAFETY: forwarded unchanged to the system allocator.
        unsafe { SystemAllocator.alloc(layout) }
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        record_alloc();
        // SAFETY: forwarded unchanged to the system allocator.
        unsafe { SystemAllocator.alloc_zeroed(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: pointer and layout came from the system allocator.
        unsafe { SystemAllocator.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        record_alloc();
        // SAFETY: pointer, layout, and new_size forwarded unchanged.
        unsafe { SystemAllocator.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOCATOR: Counting = Counting;

const DT_S: f64 = 1.0e-3;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 151,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// Deliberately allocating fixture component used to prove the counting
/// allocator and conformance-kit gate catch per-step heap allocations.
struct AllocatingDamper {
    ports: [PortSpec; 2],
    damping: f64,
}

impl AllocatingDamper {
    fn new(damping: f64) -> Self {
        Self {
            ports: [
                PortSpec::power(
                    StableId::new("a").unwrap(),
                    PortKind::MechanicalForceVelocity,
                    PortRole::FlowIn,
                    NonZeroUsize::MIN,
                )
                .unwrap(),
                PortSpec::heat(StableId::new("heat").unwrap(), PortRole::EffortIn).unwrap(),
            ],
            damping,
        }
    }
}

#[derive(Debug)]
struct FixtureErr;

impl Component for AllocatingDamper {
    type State = ();
    type Error = FixtureErr;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, (): &()) -> Result<f64, FixtureErr> {
        Ok(0.0)
    }

    fn propose(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        (): &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, FixtureErr> {
        cx.checkpoint().map_err(|_| FixtureErr)?;
        // Deliberate heap allocation inside `propose` (black-boxed via Vec).
        let scratch = Vec::from([inputs[0], inputs[1]]);
        let v = scratch[0];
        let t_k = scratch[1];
        let force = self.damping * v;
        let heat_w = force * v;
        outputs[0] = force;
        outputs[1] = -heat_w;
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j: dt_s * heat_w,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dt_s * heat_w / t_k,
            energy_tolerance_j: 64.0 * f64::EPSILON * (dt_s * heat_w.abs() + 1.0),
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

/// Exact-linear mass with companion end-of-step velocity port for Woodbury tear.
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
    type Error = FixtureErr;

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

    fn stored_energy_j(&self, v: &f64) -> Result<f64, FixtureErr> {
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
    ) -> Result<Proposal, FixtureErr> {
        cx.checkpoint().map_err(|_| FixtureErr)?;
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
    ) -> Result<(), FixtureErr> {
        cx.checkpoint().map_err(|_| FixtureErr)?;
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
    ) -> Result<(), FixtureErr> {
        cx.checkpoint().map_err(|_| FixtureErr)?;
        bar_inputs[0] = (0.5 * bar_outputs[0] + bar_outputs[1]) * dt_s / self.mass_kg;
        Ok(())
    }
}

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
    type Error = FixtureErr;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, gap: &f64) -> Result<f64, FixtureErr> {
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
    ) -> Result<Proposal, FixtureErr> {
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
    ) -> Result<Proposal, FixtureErr> {
        cx.checkpoint().map_err(|_| FixtureErr)?;
        let c_eff = (compliance[0] - compliance[1] - compliance[2] + compliance[3]).max(0.0);
        let v_rel_free = inputs[0] - inputs[1];
        let gap_free = *old + dt_s * v_rel_free;
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

fn assert_zero_alloc_steps(label: &str, mut sys: System, dt_s: f64) {
    with_cx(|cx| {
        for _ in 0..5 {
            sys.step(dt_s, cx)
                .unwrap_or_else(|e| panic!("{label} warmup: {e}"));
        }
        let allocs = allocations_of(|| {
            for step in 0..80 {
                sys.step(dt_s, cx)
                    .unwrap_or_else(|e| panic!("{label} step {step}: {e}"));
            }
        });
        assert_eq!(
            allocs, 0,
            "{label}: allocated {allocs} times over 80 warm steps"
        );
    });
}

fn build_routed_chain(n: usize) -> System {
    let mut b = SystemBuilder::new();
    let masses: Vec<_> = (0..n)
        .map(|i| {
            b.add_with_tangent(
                &format!("m{i}"),
                Mass::new(0.4).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s: 0.25,
                },
            )
            .unwrap()
        })
        .collect();
    let body_a = b
        .add_with_tangent("body_a", HeatCapacity::new(50.0).unwrap(), 300.0)
        .unwrap();
    let body_b = b
        .add_with_tangent("body_b", HeatCapacity::new(25.0).unwrap(), 310.0)
        .unwrap();
    let push = b
        .add_with_tangent(
            "push",
            Force::new(Waveform::constant(1.5).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();
    let mut ends = vec![Vec::new(); n];
    ends[0].push(b.port(push, "p").unwrap());
    let mut heats = Vec::new();
    for i in 0..n {
        let sp = b
            .add_with_tangent(&format!("sp{i}"), Spring::new(400.0, 0.0).unwrap(), 0.0)
            .unwrap();
        let dm = b
            .add_with_tangent(&format!("dm{i}"), Damper::new(1.5).unwrap(), ())
            .unwrap();
        heats.push(b.port(dm, "heat").unwrap());
        for elem in [sp, dm] {
            ends[i].push(b.port(elem, "a").unwrap());
            let far = b.port(elem, "b").unwrap();
            if i + 1 < n {
                ends[i + 1].push(far);
            } else {
                b.boundary(far, Box::new(|_, _, v| v[0] = 0.0)).unwrap();
            }
        }
    }
    for (m, members) in masses.iter().zip(&ends) {
        b.common_flow(b.port(*m, "p").unwrap(), members).unwrap();
    }
    let shares = [
        (b.port(body_a, "heat").unwrap(), 0.6),
        (b.port(body_b, "heat").unwrap(), 0.4),
    ];
    for &h in &heats {
        b.heat_route(h, &shares, None, RouteTemperature::SinkEquivalent)
            .unwrap();
    }
    let mut sys = b.seal().unwrap();
    for h in heats {
        assert!(sys.guess(h, &[303.0]));
    }
    sys
}

#[test]
fn gallery_assemblies_and_solver_modes_allocate_zero_bytes_per_warm_step() {
    // 1. Routed thermo-mechanical chain (nonlinear EquivalentTemperature driver).
    assert_zero_alloc_steps("routed_chain_interface_newton", build_routed_chain(4), DT_S);

    // 2. Same thermo-mechanical chain under FixedPointAitken and FixedPointIqnIls.
    let mut aitken_sys = build_routed_chain(2);
    aitken_sys.set_solver_policy(BlockSolverPolicy::FixedPointAitken(AitkenConfig::default()));
    assert_zero_alloc_steps("routed_chain_fixed_point_aitken", aitken_sys, 2.0e-4);

    let mut iqn_sys = build_routed_chain(2);
    iqn_sys.set_solver_policy(BlockSolverPolicy::FixedPointIqnIls(IqnIlsConfig::default()));
    assert_zero_alloc_steps("routed_chain_fixed_point_iqn_ils", iqn_sys, 2.0e-4);

    // 3. Load-line Slider (NonSmooth + affine_line).
    let mut b_slider = SystemBuilder::new();
    let a = b_slider
        .add_with_tangent(
            "a",
            Mass::new(2.0).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 1.5,
            },
        )
        .unwrap();
    let b_mass = b_slider
        .add_with_tangent(
            "b",
            Mass::new(3.0).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: -0.5,
            },
        )
        .unwrap();
    let slider = b_slider
        .add_on_line(
            "slider",
            Slider::new(20.0, 0.4, 0.25).unwrap(),
            Slip::sliding(-2.0),
        )
        .unwrap();
    b_slider
        .bond(
            b_slider.port(a, "p").unwrap(),
            b_slider.port(slider, "a").unwrap(),
        )
        .unwrap();
    b_slider
        .bond(
            b_slider.port(b_mass, "p").unwrap(),
            b_slider.port(slider, "b").unwrap(),
        )
        .unwrap();
    b_slider
        .boundary(
            b_slider.port(slider, "heat").unwrap(),
            Box::new(|_, _, v| v[0] = 300.0),
        )
        .unwrap();
    assert_zero_alloc_steps("slider_load_line", b_slider.seal().unwrap(), DT_S);

    // 4. Drum-strike shape (ToldMass + Contact with line_ends -> LoadLineWoodbury).
    let mut b_drum = SystemBuilder::new();
    let stick = b_drum
        .add_with_tangent("stick", ToldMass::new(0.04), -0.8)
        .unwrap();
    let head = b_drum
        .add_with_tangent("head", ToldMass::new(0.25), 0.0)
        .unwrap();
    let contact = b_drum
        .add_on_line("contact", DrumContact::new(1.0e5), 2.0e-4)
        .unwrap();
    b_drum
        .bond(
            b_drum.port(stick, "p").unwrap(),
            b_drum.port(contact, "a").unwrap(),
        )
        .unwrap();
    b_drum
        .bond(
            b_drum.port(head, "p").unwrap(),
            b_drum.port(contact, "b").unwrap(),
        )
        .unwrap();
    assert_zero_alloc_steps("drum_woodbury_contact", b_drum.seal().unwrap(), 2.0e-5);
}

#[test]
fn conformance_zero_alloc_gate_passes_clean_elements_and_catches_allocating_fixture() {
    with_cx(|cx| {
        // Core library elements pass the conformance-kit zero-allocation check.
        let mass = Mass::new(1.0).unwrap();
        let m0 = Motion {
            position_m: 0.1,
            velocity_m_s: 0.5,
        };
        let m_stale = Motion {
            position_m: -9.0,
            velocity_m_s: 7.0,
        };
        check_zero_alloc_propose(
            &mass,
            &Case {
                state: &m0,
                stale: &m_stale,
                dt_s: DT_S,
                inputs: &[2.0],
            },
            2,
            50,
            |f| allocations_of(f),
            cx,
        )
        .expect("Mass::propose allocates zero bytes");

        let damper = Damper::new(1.5).unwrap();
        check_zero_alloc_propose(
            &damper,
            &Case {
                state: &(),
                stale: &(),
                dt_s: DT_S,
                inputs: &[1.0, -0.5, 300.0],
            },
            2,
            50,
            |f| allocations_of(f),
            cx,
        )
        .expect("Damper::propose allocates zero bytes");

        // Deliberately allocating fixture MUST fail the conformance-kit gate.
        let leaky = AllocatingDamper::new(1.5);
        let err = check_zero_alloc_propose(
            &leaky,
            &Case {
                state: &(),
                stale: &(),
                dt_s: DT_S,
                inputs: &[1.0, 300.0],
            },
            2,
            10,
            |f| allocations_of(f),
            cx,
        )
        .expect_err("AllocatingDamper must fail zero-alloc conformance gate");
        assert!(matches!(
            err,
            Failure::Law {
                law: "propose allocates no heap memory after warm-up",
                ..
            }
        ));

        // And inside a sealed System, the counting allocator catches the same fixture.
        let mut b = SystemBuilder::new();
        let m = b
            .add_with_tangent(
                "m",
                Mass::new(1.0).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s: 1.0,
                },
            )
            .unwrap();
        let d = b.add("leaky_d", AllocatingDamper::new(1.5), ()).unwrap();
        b.bond(b.port(m, "p").unwrap(), b.port(d, "a").unwrap())
            .unwrap();
        b.boundary(b.port(d, "heat").unwrap(), Box::new(|_, _, v| v[0] = 300.0))
            .unwrap();
        let mut sys = b.seal().unwrap();
        for _ in 0..3 {
            sys.step(DT_S, cx).unwrap();
        }
        let allocs = allocations_of(|| {
            for _ in 0..10 {
                sys.step(DT_S, cx).unwrap();
            }
        });
        assert!(
            allocs >= 10,
            "counting allocator must detect per-step allocations from fixture (got {allocs})"
        );
    });
}
