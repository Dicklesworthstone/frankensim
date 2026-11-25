//! Bead `frankensim-rc-1109-root-3gc12.5.6`: Window ledger from the driver —
//! energy, entropy, mass, momentum, charge, evidence coloring, leak attribution,
//! 1e6-step lossless closure, and export structure.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::{Color, validate_color_payload};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::StableId;
use fs_ports::component::{
    Class, Component, PortRole, PortSpec, Proposal, Receipt, StreamBasis, StreamQuantity,
};
use fs_qty::{Current, Length, Mass as MassQty, Time};
use fs_system::elements::{Damper, GasVolume, HeatCapacity, Mass, Motion, Restriction, Spring};
use fs_system::{System, SystemBuilder};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let gate = CancelGate::new_clock_free();
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 56,
                kernel_id: 56,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// A portless component that declares an intentional `leak_j` energy mismatch
/// inside a wide claimed tolerance so the per-step audit passes it into the
/// window ledger's defect accounts, letting the ledger attribute the leak to
/// the exact owning component.
struct LeakyElement {
    leak_j: f64,
    leak_step: u64,
}

impl Component for LeakyElement {
    type State = (f64, u64);
    type Error = &'static str;

    fn ports(&self) -> &[PortSpec] {
        &[]
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, Self::Error> {
        Ok(state.0)
    }

    fn propose(
        &self,
        old: &Self::State,
        _dt_s: f64,
        _inputs: &[f64],
        new: &mut Self::State,
        _outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        cx.checkpoint().map_err(|_| "cancelled")?;
        let change = if old.1 == self.leak_step {
            self.leak_j
        } else {
            0.0
        };
        *new = (old.0 + change, old.1 + 1);
        Ok(Proposal::Step(Receipt::lossless(
            change,
            self.leak_j.abs() * 10.0 + 1.0e-12,
        )))
    }
}

/// A component whose integrator loses a small known amount of numerical energy
/// per step without physically dissipating anything.
struct NumericalLossElement {
    loss_per_step_j: f64,
}

impl Component for NumericalLossElement {
    type State = f64;
    type Error = &'static str;

    fn ports(&self) -> &[PortSpec] {
        &[]
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, Self::Error> {
        Ok(*state)
    }

    fn propose(
        &self,
        old: &Self::State,
        _dt_s: f64,
        _inputs: &[f64],
        new: &mut Self::State,
        _outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        cx.checkpoint().map_err(|_| "cancelled")?;
        let change = -self.loss_per_step_j;
        *new = old + change;
        Ok(Proposal::Step(Receipt {
            stored_change_j: change,
            dissipated_j: 0.0,
            numerical_loss_j: self.loss_per_step_j,
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: 1.0e-14,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

fn build_lossless_oscillator() -> System {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(1.0).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 1.0,
            },
        )
        .unwrap();
    let spring = builder
        .add_with_tangent("spring", Spring::new(100.0, 0.0).unwrap(), 0.0)
        .unwrap();
    builder
        .bond(
            builder.port(mass, "p").unwrap(),
            builder.port(spring, "a").unwrap(),
        )
        .unwrap();
    builder
        .boundary(
            builder.port(spring, "b").unwrap(),
            Box::new(|_, _, out| out[0] = 0.0),
        )
        .unwrap();
    builder.seal().unwrap()
}

#[test]
fn lossless_graph_closes_to_roundoff_over_one_million_steps_and_verifies_color() {
    let mut sys = build_lossless_oscillator();
    let dt_s = 1.0e-3;
    with_cx(|cx| {
        for _ in 0..1_000_000 {
            sys.step(dt_s, cx).unwrap();
        }
    });
    let ledger = sys.ledger();
    assert_eq!(ledger.steps, 1_000_000);
    assert!(
        (ledger.duration_s.value() - 1_000.0).abs() <= 1.0e-9,
        "duration {}",
        ledger.duration_s.value()
    );
    assert!(
        ledger.energy_residual_j().abs() <= 1.0e-11,
        "1e6-step lossless energy residual too large: {:e} J",
        ledger.energy_residual_j()
    );
    assert_eq!(ledger.total_dissipated_j(), 0.0);
    assert_eq!(ledger.numerical_loss_j(), 0.0);

    let outcome = ledger.outcome();
    assert!(validate_color_payload(&outcome.overall).is_ok());
    assert!(validate_color_payload(&outcome.energy).is_ok());
    assert!(validate_color_payload(&outcome.entropy).is_ok());
    assert!(
        matches!(outcome.overall, Color::Verified { lo, hi } if lo <= 0.0 && hi >= 0.0),
        "expected Verified window outcome for 1e6-step lossless closure, got {:?}",
        outcome.overall
    );
    assert!(matches!(outcome.energy, Color::Verified { .. }));
    assert!(matches!(outcome.entropy, Color::Verified { .. }));
    assert!(outcome.mass.is_none());
    assert_eq!(outcome.carried, Vec::new());

    let export = ledger.export();
    assert_eq!(export.outcome, outcome);
    assert!(
        export
            .canonical_json
            .contains("\"schema\":\"frankensim.window_ledger.v1\"")
    );
    assert!(export.canonical_json.contains("\"color\":\"verified\""));
    assert!(export.report_markdown.contains("overall: `verified`"));
    assert!(
        export
            .metrics
            .iter()
            .any(|(k, v)| k == "window.steps" && (*v - 1_000_000.0).abs() == 0.0)
    );
}

#[test]
fn injected_1e_9_j_leak_is_caught_attributed_to_owning_component_and_demotes_color() {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "clean_mass",
            Mass::new(1.0).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 1.0,
            },
        )
        .unwrap();
    let spring = builder
        .add_with_tangent("clean_spring", Spring::new(25.0, 0.0).unwrap(), 0.0)
        .unwrap();
    builder
        .bond(
            builder.port(mass, "p").unwrap(),
            builder.port(spring, "a").unwrap(),
        )
        .unwrap();
    builder
        .boundary(
            builder.port(spring, "b").unwrap(),
            Box::new(|_, _, out| out[0] = 0.0),
        )
        .unwrap();
    builder
        .add(
            "suspect_valve",
            LeakyElement {
                leak_j: 1.0e-9,
                leak_step: 4,
            },
            (10.0, 0),
        )
        .unwrap();
    let mut sys = builder.seal().unwrap();

    with_cx(|cx| {
        for _ in 0..10 {
            sys.step(1.0e-3, cx).unwrap();
        }
    });

    let ledger = sys.ledger();
    assert!(
        (ledger.energy_residual_j() - 1.0e-9).abs() <= 1.0e-14,
        "residual = {:e}",
        ledger.energy_residual_j()
    );
    let leak = ledger
        .worst_energy_leak_node()
        .expect("injected 1e-9 J leak must be attributed");
    assert_eq!(leak.node, "suspect_valve");
    assert!((leak.balance_defect_j + 1.0e-9).abs() <= 1.0e-15);
    assert!((leak.absolute_defect_j - 1.0e-9).abs() <= 1.0e-15);

    let outcome = ledger.outcome();
    assert!(validate_color_payload(&outcome.overall).is_ok());
    assert!(validate_color_payload(&outcome.energy).is_ok());
    assert!(
        matches!(&outcome.energy, Color::Estimated { estimator, dispersion } if estimator == "window-ledger:energy-defect" && *dispersion > 0.0),
        "1e-9 J leak must demote energy color to Estimated, got {:?}",
        outcome.energy
    );
    assert!(matches!(outcome.overall, Color::Estimated { .. }));
    assert_eq!(
        outcome.worst_energy_leak_node.as_deref(),
        Some("suspect_valve")
    );
    assert!((outcome.worst_energy_leak_j - 1.0e-9).abs() <= 1.0e-15);

    let export = ledger.export();
    assert!(export.canonical_json.contains("\"node\":\"suspect_valve\""));
    assert!(export.report_markdown.contains("`suspect_valve`"));
    assert!(export.report_markdown.contains("overall: `estimated`"));
}

#[test]
fn numerical_loss_is_never_booked_as_physical_dissipation_and_stays_estimated() {
    let mut builder = SystemBuilder::new();
    builder
        .add(
            "dissipative_integrator",
            NumericalLossElement {
                loss_per_step_j: 2.5e-6,
            },
            1.0,
        )
        .unwrap();
    let mut sys = builder.seal().unwrap();

    with_cx(|cx| {
        for _ in 0..4 {
            sys.step(1.0e-2, cx).unwrap();
        }
    });

    let ledger = sys.ledger();
    assert_eq!(
        ledger.total_dissipated_j(),
        0.0,
        "numerical loss must NEVER be booked as physical dissipation"
    );
    assert_eq!(ledger.nodes[0].dissipated_j.value(), 0.0);
    assert!((ledger.numerical_loss_j() - 1.0e-5).abs() <= 1.0e-18);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-18);

    let outcome = ledger.outcome();
    assert!(validate_color_payload(&outcome.energy).is_ok());
    assert!(validate_color_payload(&outcome.overall).is_ok());
    assert!(
        matches!(&outcome.energy, Color::Estimated { estimator, .. } if estimator == "window-ledger:numerical-loss"),
        "nonzero numerical loss must remain Estimated, got {:?}",
        outcome.energy
    );
    assert!(matches!(outcome.overall, Color::Estimated { .. }));

    // Contrast: physical dissipation routed into a heat capacity closes First Law
    // with zero numerical loss and positive entropy generation -> Verified.
    let mut b2 = SystemBuilder::new();
    let m = b2
        .add_with_tangent(
            "m",
            Mass::new(1.0).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 2.0,
            },
        )
        .unwrap();
    let d = b2
        .add_with_tangent("d", Damper::new(1.5).unwrap(), ())
        .unwrap();
    let c = b2
        .add_with_tangent("c", HeatCapacity::new(50.0).unwrap(), 300.0)
        .unwrap();
    b2.bond(b2.port(m, "p").unwrap(), b2.port(d, "a").unwrap())
        .unwrap();
    b2.boundary(b2.port(d, "b").unwrap(), Box::new(|_, _, o| o[0] = 0.0))
        .unwrap();
    b2.bond(b2.port(c, "heat").unwrap(), b2.port(d, "heat").unwrap())
        .unwrap();
    let mut damped = b2.seal().unwrap();
    assert!(damped.guess(damped.port(d, "heat").unwrap(), &[300.0]));
    with_cx(|cx| {
        for _ in 0..20 {
            damped.step(1.0e-3, cx).unwrap();
        }
    });
    assert!(damped.ledger().total_dissipated_j() > 0.0);
    assert_eq!(damped.ledger().numerical_loss_j(), 0.0);
    assert!(damped.ledger().entropy_generated_j_per_k() > 0.0);
    assert!(matches!(
        damped.ledger().outcome().overall,
        Color::Verified { .. }
    ));
}

/// Minimal storage carrying identified stream coordinates for momentum and charge
/// alongside mass, allowing verification of conserved inventory closure and
/// incomplete-measurement color demotion.
struct ConservedCarriedTank {
    ports: [PortSpec; 1],
    measure_inventories: bool,
}

#[derive(Clone, Debug)]
struct TankState {
    mass_kg: f64,
    energy_j: f64,
    momentum_kg_m_s: f64,
    charge_c: f64,
}

impl ConservedCarriedTank {
    fn basis() -> StreamBasis {
        let momentum_dims = MassQty::DIMS
            .checked_plus(Length::DIMS)
            .unwrap()
            .checked_minus(Time::DIMS)
            .unwrap();
        let charge_dims = Current::DIMS.checked_plus(Time::DIMS).unwrap();
        StreamBasis::new(
            StableId::new("charged-fluid-v1").unwrap(),
            StableId::new("zero-datum").unwrap(),
            vec![
                StreamQuantity::Named {
                    id: StableId::new("linear-momentum-x").unwrap(),
                    dimensions: momentum_dims,
                },
                StreamQuantity::Named {
                    id: StableId::new("electric-charge").unwrap(),
                    dimensions: charge_dims,
                },
            ],
        )
        .unwrap()
    }

    fn new(measure_inventories: bool) -> Self {
        let port = PortSpec::identified_stream(
            StableId::new("stream").unwrap(),
            PortRole::FlowIn,
            Self::basis(),
        )
        .unwrap();
        Self {
            ports: [port],
            measure_inventories,
        }
    }
}

impl Component for ConservedCarriedTank {
    type State = TankState;
    type Error = &'static str;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, Self::Error> {
        Ok(state.energy_j)
    }

    fn stored_mass_kg(&self, state: &Self::State) -> Result<Option<f64>, Self::Error> {
        Ok(self.measure_inventories.then_some(state.mass_kg))
    }

    fn stored_carried_quantity(
        &self,
        state: &Self::State,
        carried: usize,
    ) -> Result<Option<f64>, Self::Error> {
        if !self.measure_inventories {
            return Ok(None);
        }
        match carried {
            0 => Ok(Some(state.momentum_kg_m_s)),
            1 => Ok(Some(state.charge_c)),
            _ => Ok(None),
        }
    }

    fn propose(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        cx.checkpoint().map_err(|_| "cancelled")?;
        let dm = dt_s * inputs[0];
        let de = dt_s * inputs[1];
        let dp = dt_s * inputs[2];
        let dq = dt_s * inputs[3];
        *new = TankState {
            mass_kg: old.mass_kg + dm,
            energy_j: old.energy_j + de,
            momentum_kg_m_s: old.momentum_kg_m_s + dp,
            charge_c: old.charge_c + dq,
        };
        outputs.copy_from_slice(&[1.0e5, 1.0e3, 2.0, -0.5]);
        Ok(Proposal::Step(Receipt::lossless(de, 1.0e-12)))
    }
}

#[test]
fn mass_momentum_and_charge_inventories_verify_when_measured_and_demote_when_unmeasured() {
    // 1. Measured tank: mass, momentum, and charge close to roundoff -> Verified.
    let mut b = SystemBuilder::new();
    let tank = b
        .add(
            "measured_tank",
            ConservedCarriedTank::new(true),
            TankState {
                mass_kg: 2.0,
                energy_j: 2000.0,
                momentum_kg_m_s: 4.0,
                charge_c: -1.0,
            },
        )
        .unwrap();
    b.boundary(
        b.port(tank, "stream").unwrap(),
        Box::new(|_, _, vals| vals.copy_from_slice(&[0.5, 500.0, 1.25, -0.25])),
    )
    .unwrap();
    let mut sys = b.seal().unwrap();
    with_cx(|cx| {
        for _ in 0..25 {
            sys.step(0.01, cx).unwrap();
        }
    });
    let outcome = sys.ledger().outcome();
    assert!(validate_color_payload(&outcome.overall).is_ok());
    assert!(
        matches!(outcome.overall, Color::Verified { .. }),
        "{outcome:?}"
    );
    assert!(matches!(outcome.mass, Some(Color::Verified { .. })));
    assert_eq!(outcome.carried.len(), 2);
    assert!(matches!(outcome.carried[0].color, Color::Verified { .. }));
    assert!(matches!(outcome.carried[1].color, Color::Verified { .. }));
    assert!(sys.ledger().mass_residual_kg().unwrap().abs() <= 1.0e-14);
    assert!(sys.ledger().carried_residual(0).unwrap().abs() <= 1.0e-14);
    assert!(sys.ledger().carried_residual(1).unwrap().abs() <= 1.0e-14);

    // 2. Unmeasured tank: missing inventory measurement NEVER launders to Verified.
    let mut b2 = SystemBuilder::new();
    let tank2 = b2
        .add(
            "unmeasured_tank",
            ConservedCarriedTank::new(false),
            TankState {
                mass_kg: 2.0,
                energy_j: 2000.0,
                momentum_kg_m_s: 4.0,
                charge_c: -1.0,
            },
        )
        .unwrap();
    b2.boundary(
        b2.port(tank2, "stream").unwrap(),
        Box::new(|_, _, vals| vals.copy_from_slice(&[0.5, 500.0, 1.25, -0.25])),
    )
    .unwrap();
    let mut sys2 = b2.seal().unwrap();
    with_cx(|cx| {
        for _ in 0..5 {
            sys2.step(0.01, cx).unwrap();
        }
    });
    let outcome2 = sys2.ledger().outcome();
    assert!(validate_color_payload(&outcome2.overall).is_ok());
    assert!(matches!(outcome2.energy, Color::Verified { .. }));
    assert!(matches!(
        outcome2.mass,
        Some(Color::Estimated { ref estimator, .. }) if estimator == "window-ledger:mass-unmeasured"
    ));
    assert!(matches!(outcome2.overall, Color::Estimated { .. }));
}

#[test]
fn coupled_gas_inventories_verify_and_entropy_or_dissipation_corruption_demotes() {
    // 3. Coupled two-vessel GasVolume + Restriction stream exchange closes mass and species.
    let gas = GasVolume::new(287.0, 1.4, 1).unwrap();
    let mut b3 = SystemBuilder::new();
    let pipe = b3
        .add_with_tangent("pipe", Restriction::for_gas(1.0e-8, &gas).unwrap(), ())
        .unwrap();
    let a = b3
        .add_with_tangent(
            "a",
            gas.clone(),
            gas.contents(2.0e5, 320.0, 0.02, &[0.7]).unwrap(),
        )
        .unwrap();
    let c = b3
        .add_with_tangent(
            "b",
            gas.clone(),
            gas.contents(1.5e5, 300.0, 0.02, &[0.2]).unwrap(),
        )
        .unwrap();
    b3.bond(b3.port(a, "stream").unwrap(), b3.port(pipe, "a").unwrap())
        .unwrap();
    b3.bond(b3.port(c, "stream").unwrap(), b3.port(pipe, "b").unwrap())
        .unwrap();
    b3.boundary(b3.port(a, "heat").unwrap(), Box::new(|_, _, v| v[0] = 0.0))
        .unwrap();
    b3.boundary(b3.port(c, "heat").unwrap(), Box::new(|_, _, v| v[0] = 0.0))
        .unwrap();
    let mut coupled_gas = b3.seal().unwrap();
    with_cx(|cx| {
        for _ in 0..40 {
            coupled_gas.step(0.005, cx).unwrap();
        }
    });
    let gas_outcome = coupled_gas.ledger().outcome();
    assert!(
        matches!(gas_outcome.overall, Color::Verified { .. }),
        "{gas_outcome:?}"
    );
    let gas_export = coupled_gas.ledger().export();
    assert!(
        gas_export
            .metrics
            .iter()
            .any(|(k, _)| k == "window.mass.residual_kg")
    );
    assert!(
        gas_export
            .metrics
            .iter()
            .any(|(k, _)| k == "window.carried.0.residual")
    );

    // 4. Negative connection entropy or negative node dissipation must never verify.
    let mut corrupted = coupled_gas.ledger().clone();
    corrupted.connections[0]
        .entropy_generated_j_per_k
        .add(-1.0e-2);
    let bad_ent = corrupted.outcome();
    assert!(validate_color_payload(&bad_ent.entropy).is_ok());
    assert!(
        matches!(&bad_ent.entropy, Color::Estimated { estimator, .. } if estimator == "window-ledger:entropy-defect"),
        "negative connection entropy must demote entropy color, got {:?}",
        bad_ent.entropy
    );
    assert!(matches!(bad_ent.overall, Color::Estimated { .. }));

    let mut neg_diss = coupled_gas.ledger().clone();
    neg_diss.nodes[0].dissipated_j.add(-1.0e-3);
    let bad_diss = neg_diss.outcome();
    assert!(validate_color_payload(&bad_diss.energy).is_ok());
    assert!(
        matches!(&bad_diss.energy, Color::Estimated { estimator, .. } if estimator == "window-ledger:negative-dissipation"),
        "negative node dissipation must demote energy color, got {:?}",
        bad_diss.energy
    );
    assert!(matches!(bad_diss.overall, Color::Estimated { .. }));
}
