//! G0/G2: ideal transformers and rectangular power maps, with trajectories
//! checked against independent monolithic solutions rather than their ledger.

#![allow(clippy::float_cmp)]

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_ports::component::{
    Carrier, Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt,
};
use fs_ports::{PortKind, StableId};
use fs_system::elements::{Force, Mass, Motion, SourceAccount, Spring, Waveform};
use fs_system::{AdmissionError, PortRef, SolveConfig, SystemBuilder};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 11,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx);
    });
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} versus {expected}"
    );
}

#[test]
fn a_geared_oscillator_has_reflected_stiffness_and_the_midpoint_phase() {
    for ratio in [2.75_f64, -1.8] {
        let (mass_kg, stiffness, x0, v0) = (2.0, 70.0, 0.13, -0.4);
        let mut builder = SystemBuilder::new();
        let mass = builder
            .add_with_tangent(
                "mass",
                Mass::new(mass_kg).unwrap(),
                Motion {
                    position_m: x0,
                    velocity_m_s: v0,
                },
            )
            .unwrap();
        let spring = builder
            .add_with_tangent("spring", Spring::new(stiffness, 0.0).unwrap(), -ratio * x0)
            .unwrap();
        let a = builder.port(mass, "p").unwrap();
        let b = builder.port(spring, "a").unwrap();
        let anchor = builder.port(spring, "b").unwrap();
        builder.transformer(a, b, ratio).unwrap();
        builder
            .boundary(anchor, Box::new(|_, _, values| values.fill(0.0)))
            .unwrap();
        let mut system = builder.seal().unwrap();
        // The affine loop has to converge with one exact Newton update.
        system.set_config(SolveConfig {
            max_iterations: 1,
            ..SolveConfig::default()
        });

        let (dt, steps) = (0.015, 2_500);
        let omega = (stiffness * ratio * ratio / mass_kg).sqrt();
        let energy0 = 0.5 * mass_kg * v0 * v0 + 0.5 * stiffness * (ratio * x0).powi(2);
        with_cx(|cx| {
            for i in 1..=steps {
                system.step(dt, cx).unwrap();
                let angle = f64::from(i) * 2.0 * (0.5 * omega * dt).atan();
                let (sin, cos) = angle.sin_cos();
                let x = x0 * cos + v0 / omega * sin;
                let v = v0 * cos - omega * x0 * sin;
                let motion = system.state::<Motion>(mass).unwrap();
                let q = *system.state::<f64>(spring).unwrap();
                close(motion.position_m, x, 2.0e-11);
                close(motion.velocity_m_s, v, 3.0e-10);
                close(q, -ratio * x, 6.0e-11);
                let energy = 0.5 * mass_kg * motion.velocity_m_s.powi(2) + 0.5 * stiffness * q * q;
                close(energy, energy0, 2.0e-11 * energy0);
            }
        });
        close(
            system.ledger().connections[0].energy_defect_j.value(),
            0.0,
            1.0e-11 * energy0,
        );
        close(system.ledger().energy_residual_j(), 0.0, 1.0e-11 * energy0);
    }
}

#[test]
fn a_pressure_area_transformer_drives_a_piston_and_accounts_for_source_work() {
    let (area, pressure, mass_kg, dt, steps) = (0.0025, 32_000.0, 4.0, 0.001, 1_000);
    let mut builder = SystemBuilder::new();
    let piston = builder
        .add_with_tangent(
            "piston",
            Mass::new(mass_kg).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let reservoir = builder
        .add_with_tangent(
            "reservoir",
            Force::in_domain(
                PortKind::FluidPressureFlux,
                Waveform::constant(pressure).unwrap(),
            ),
            SourceAccount::default(),
        )
        .unwrap();
    let a = builder.port(piston, "p").unwrap();
    let b = builder.port(reservoir, "p").unwrap();
    assert!(matches!(
        builder.bond(a, b),
        Err(AdmissionError::Mismatch { .. })
    ));
    // An explicit area, in m², bridges pressure/volume flux and force/velocity.
    builder.transformer(a, b, area).unwrap();
    let mut system = builder.seal().unwrap();
    with_cx(|cx| {
        for _ in 0..steps {
            system.step(dt, cx).unwrap();
        }
    });
    let time = f64::from(steps) * dt;
    let acceleration = area * pressure / mass_kg;
    let motion = system.state::<Motion>(piston).unwrap();
    close(motion.position_m, 0.5 * acceleration * time * time, 1.0e-10);
    close(motion.velocity_m_s, acceleration * time, 1.0e-10);
    let kinetic = 0.5 * mass_kg * motion.velocity_m_s.powi(2);
    let delivered = -system
        .state::<SourceAccount>(reservoir)
        .unwrap()
        .energy_in_j;
    close(delivered, area * pressure * motion.position_m, 1.0e-10);
    close(delivered, kinetic, 1.0e-10);
    close(system.ledger().energy_residual_j(), 0.0, 1.0e-10);
}

fn power(role: PortRole, width: usize) -> PortSpec {
    PortSpec::power(
        StableId::new("p").unwrap(),
        PortKind::MechanicalForceVelocity,
        role,
        NonZeroUsize::new(width).unwrap(),
    )
    .unwrap()
}

/// Independent diagonal storage: an inertia bank when effort is the input,
/// a spring bank when flow is the input. State is velocity or extension.
struct Storage {
    ports: [PortSpec; 1],
    coefficients: Vec<f64>,
}

impl Storage {
    fn new(role: PortRole, coefficients: &[f64]) -> Self {
        Self {
            ports: [power(role, coefficients.len())],
            coefficients: coefficients.to_vec(),
        }
    }

    fn inertial(&self) -> bool {
        self.ports[0].role() == PortRole::EffortIn
    }
}

impl Component for Storage {
    type State = Vec<f64>;
    type Error = Cancelled;
    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }
    fn class(&self) -> Class {
        Class::ExactLinear
    }
    fn stored_energy_j(&self, state: &Vec<f64>) -> Result<f64, Cancelled> {
        Ok(self
            .coefficients
            .iter()
            .zip(state)
            .map(|(c, x)| 0.5 * c * x * x)
            .sum())
    }
    fn propose(
        &self,
        old: &Vec<f64>,
        dt: f64,
        input: &[f64],
        new: &mut Vec<f64>,
        output: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        let mut change = 0.0;
        for i in 0..old.len() {
            let c = self.coefficients[i];
            new[i] = old[i] + dt * input[i] / if self.inertial() { c } else { 1.0 };
            let mean = 0.5 * (new[i] + old[i]);
            output[i] = mean * if self.inertial() { 1.0 } else { c };
            change += c * mean * (new[i] - old[i]);
        }
        let scale = self.stored_energy_j(old)? + self.stored_energy_j(new)? + change.abs();
        Ok(Proposal::Step(Receipt::lossless(
            change,
            128.0 * f64::EPSILON * scale,
        )))
    }
}

impl PortTangent for Storage {
    fn jvp(
        &self,
        _old: &Vec<f64>,
        dt: f64,
        _input: &[f64],
        direction: &[f64],
        output: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        for i in 0..direction.len() {
            let c = self.coefficients[i];
            output[i] = 0.5 * dt * direction[i] * if self.inertial() { 1.0 / c } else { c };
        }
        Ok(())
    }
    fn vjp(
        &self,
        old: &Vec<f64>,
        dt: f64,
        input: &[f64],
        direction: &[f64],
        output: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        self.jvp(old, dt, input, direction, output, cx)
    }
}

#[test]
fn rectangular_maps_match_the_monolithic_solve_including_unconstrained_modes() {
    for map in [
        [[1.0, 2.0], [0.5, -1.0], [-2.0, 0.25]],
        [[1.0, 2.0], [2.0, 4.0], [-0.5, -1.0]], // rank one
        [[0.0; 2]; 3],
    ] {
        let (m, k) = ([1.5, 2.5], [8.0, 12.0, 5.0]);
        let (mut velocity, mut extension) = ([0.35, -0.7], [0.1, -0.2, 0.4]);
        let energy = |v: &[f64], q: &[f64]| -> f64 {
            m.iter()
                .zip(v)
                .chain(k.iter().zip(q))
                .map(|(c, x)| 0.5 * c * x * x)
                .sum()
        };
        let energy0 = energy(&velocity, &extension);
        let mut builder = SystemBuilder::new();
        let inertias = builder
            .add_with_tangent(
                "inertias",
                Storage::new(PortRole::EffortIn, &m),
                velocity.to_vec(),
            )
            .unwrap();
        let springs = builder
            .add_with_tangent(
                "springs",
                Storage::new(PortRole::FlowIn, &k),
                extension.to_vec(),
            )
            .unwrap();
        let a = builder.port(inertias, "p").unwrap();
        let b = builder.port(springs, "p").unwrap();
        builder
            .power_map(a, b, map.into_iter().flatten().collect())
            .unwrap();
        let mut system = builder.seal().unwrap();
        system.set_config(SolveConfig {
            max_iterations: 1,
            ..SolveConfig::default()
        });
        let dt = 0.03;
        // Condense the spring variables analytically. The two-dimensional
        // Schur matrix is M + dt²/4 Pᵀ K P; invert it by its determinant.
        let mut schur = [[m[0], 0.0], [0.0, m[1]]];
        for i in 0..2 {
            for j in 0..2 {
                schur[i][j] +=
                    0.25 * dt * dt * (0..3).map(|r| map[r][i] * k[r] * map[r][j]).sum::<f64>();
            }
        }
        let det = schur[0][0] * schur[1][1] - schur[0][1] * schur[1][0];
        with_cx(|cx| {
            for _ in 0..1_500 {
                let rhs = std::array::from_fn::<_, 2, _>(|i| {
                    m[i] * velocity[i]
                        + 0.5 * dt * (0..3).map(|r| map[r][i] * k[r] * extension[r]).sum::<f64>()
                });
                let mean = [
                    (schur[1][1] * rhs[0] - schur[0][1] * rhs[1]) / det,
                    (schur[0][0] * rhs[1] - schur[1][0] * rhs[0]) / det,
                ];
                for i in 0..2 {
                    velocity[i] = 2.0 * mean[i] - velocity[i];
                }
                for r in 0..3 {
                    extension[r] -= dt * (map[r][0] * mean[0] + map[r][1] * mean[1]);
                }
                system.step(dt, cx).unwrap();
                let v = system.state::<Vec<f64>>(inertias).unwrap();
                let q = system.state::<Vec<f64>>(springs).unwrap();
                for i in 0..2 {
                    close(v[i], velocity[i], 2.0e-10);
                }
                for r in 0..3 {
                    close(q[r], extension[r], 2.0e-10);
                }
                close(energy(v, q), energy0, 2.0e-11 * energy0);
            }
        });
        close(
            system.ledger().connections[0].energy_defect_j.value(),
            0.0,
            1.0e-11 * energy0,
        );
        close(system.ledger().energy_residual_j(), 0.0, 1.0e-11 * energy0);
    }
}

/// Zero exchange, except that a heat receiver reports 300 K. Used only to
/// exercise description admission and prove a repaired builder can step.
struct TerminatedPort([PortSpec; 1]);
impl Component for TerminatedPort {
    type State = ();
    type Error = Cancelled;
    fn ports(&self) -> &[PortSpec] {
        &self.0
    }
    fn class(&self) -> Class {
        Class::ExactLinear
    }
    fn stored_energy_j(&self, _: &()) -> Result<f64, Cancelled> {
        Ok(0.0)
    }
    fn propose(
        &self,
        _: &(),
        _: f64,
        _: &[f64],
        _: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        outputs.fill(
            if self.0[0].carrier() == Carrier::Heat && self.0[0].role() == PortRole::FlowIn {
                300.0
            } else {
                0.0
            },
        );
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

fn terminated(builder: &mut SystemBuilder, name: &str, spec: PortSpec) -> PortRef {
    let node = builder.add(name, TerminatedPort([spec]), ()).unwrap();
    builder.port(node, "p").unwrap()
}

#[test]
fn invalid_transformers_and_maps_leave_every_unclaimed_port_usable() {
    let mut builder = SystemBuilder::new();
    let a = terminated(&mut builder, "a", power(PortRole::EffortIn, 1));
    let b = terminated(&mut builder, "b", power(PortRole::FlowIn, 1));
    let wide_a = terminated(&mut builder, "wide_a", power(PortRole::EffortIn, 2));
    let wide_b = terminated(&mut builder, "wide_b", power(PortRole::FlowIn, 3));
    let hot = terminated(
        &mut builder,
        "hot",
        PortSpec::heat(StableId::new("p").unwrap(), PortRole::EffortIn).unwrap(),
    );
    let cold = terminated(
        &mut builder,
        "cold",
        PortSpec::heat(StableId::new("p").unwrap(), PortRole::FlowIn).unwrap(),
    );
    let signal = terminated(
        &mut builder,
        "signal",
        PortSpec::signal_in(StableId::new("p").unwrap(), NonZeroUsize::MIN),
    );

    for ratio in [0.0, -0.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            builder.transformer(a, b, ratio),
            Err(AdmissionError::TransformerRatio { ratio: rejected }) if rejected.to_bits() == ratio.to_bits()
        ));
    }
    for (left, right) in [
        (b, a),
        (a, a),
        (a, hot),
        (cold, b),
        (hot, cold),
        (a, signal),
        (signal, b),
    ] {
        assert!(matches!(
            builder.transformer(left, right, 2.0),
            Err(AdmissionError::Role { .. })
        ));
        assert!(matches!(
            builder.power_map(left, right, vec![2.0]),
            Err(AdmissionError::Role { .. })
        ));
    }
    assert!(matches!(
        builder.transformer(wide_a, wide_b, 2.0),
        Err(AdmissionError::Mismatch { .. })
    ));
    for entries in [0, 5, 7] {
        assert!(matches!(
            builder.power_map(wide_a, wide_b, vec![1.0; entries]),
            Err(AdmissionError::MapShape { rows: 3, columns: 2, entries: rejected }) if rejected == entries
        ));
    }
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            builder.power_map(wide_a, wide_b, vec![value; 6]),
            Err(AdmissionError::MapCoefficient { index: 0, value: rejected }) if rejected.to_bits() == value.to_bits()
        ));
    }
    builder.transformer(a, b, -2.0).unwrap();
    // A refusal on the second endpoint must not claim the first endpoint.
    assert!(matches!(
        builder.power_map(wide_a, b, vec![1.0, 2.0]),
        Err(AdmissionError::AlreadyConnected { .. })
    ));
    assert!(matches!(
        builder.transformer(a, b, 1.0),
        Err(AdmissionError::AlreadyConnected { .. })
    ));
    builder
        .power_map(wide_a, wide_b, vec![1.0, 2.0, 0.0, 1.0, -1.0, 0.0])
        .unwrap();
    builder
        .boundary(hot, Box::new(|_, _, values| values.fill(300.0)))
        .unwrap();
    builder
        .boundary(cold, Box::new(|_, _, values| values.fill(0.0)))
        .unwrap();
    builder
        .boundary(signal, Box::new(|_, _, values| values.fill(0.0)))
        .unwrap();
    let mut system = builder.seal().unwrap();
    with_cx(|cx| {
        system.step(0.01, cx).unwrap();
    });
    assert_eq!(system.steps(), 1);
    assert_eq!(system.ledger().connections.len(), 2);
    assert_eq!(system.ledger().energy_residual_j(), 0.0);
}

/// The geared oscillator of the first test after 400 steps, wired by
/// `transformer` or as a junction whose one member has that ratio.
fn geared_oscillator(as_junction: bool, ratio: f64) -> (Motion, f64) {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(2.0).unwrap(),
            Motion {
                position_m: 0.13,
                velocity_m_s: -0.4,
            },
        )
        .unwrap();
    let spring = builder
        .add_with_tangent("spring", Spring::new(70.0, 0.0).unwrap(), -ratio * 0.13)
        .unwrap();
    let a = builder.port(mass, "p").unwrap();
    let b = builder.port(spring, "a").unwrap();
    let anchor = builder.port(spring, "b").unwrap();
    if as_junction {
        builder.common_flow_ratios(a, &[(b, ratio)]).unwrap();
    } else {
        builder.transformer(a, b, ratio).unwrap();
    }
    builder
        .boundary(anchor, Box::new(|_, _, values| values.fill(0.0)))
        .unwrap();
    let mut system = builder.seal().unwrap();
    with_cx(|cx| {
        for _ in 0..400 {
            system.step(0.015, cx).unwrap();
        }
    });
    (
        *system.state::<Motion>(mass).unwrap(),
        *system.state::<f64>(spring).unwrap(),
    )
}

#[test]
fn a_junction_with_one_ratio_is_the_transformer_to_the_bit() {
    for ratio in [2.75_f64, -1.8] {
        let (by_transformer, spring_by_transformer) = geared_oscillator(false, ratio);
        let (by_junction, spring_by_junction) = geared_oscillator(true, ratio);
        assert_ne!(by_transformer.position_m, 0.13);
        assert_eq!(
            by_junction.position_m.to_bits(),
            by_transformer.position_m.to_bits()
        );
        assert_eq!(
            by_junction.velocity_m_s.to_bits(),
            by_transformer.velocity_m_s.to_bits()
        );
        assert_eq!(
            spring_by_junction.to_bits(),
            spring_by_transformer.to_bits()
        );
    }
}

#[test]
fn a_mass_between_two_levered_springs_and_a_force_is_the_reflected_oscillator() {
    // m x'' = -(r1² k1 + r2² k2) x + s F: the mass sees each spring through
    // its own lever and the force through another. The midpoint rule is exact
    // for the constant offset and turns the phase by 2 atan(w dt / 2) a step.
    let (mass_kg, x0, v0) = (1.5, 0.08, 0.3);
    let springs = [(40.0, 1.7), (25.0, -0.6)];
    let (force_n, lever) = (3.0, -2.5);
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(mass_kg).unwrap(),
            Motion {
                position_m: x0,
                velocity_m_s: v0,
            },
        )
        .unwrap();
    let mut members = Vec::new();
    let mut spring_nodes = Vec::new();
    for (name, (stiffness, ratio)) in ["near", "far"].into_iter().zip(springs) {
        let spring = builder
            .add_with_tangent(name, Spring::new(stiffness, 0.0).unwrap(), -ratio * x0)
            .unwrap();
        members.push((builder.port(spring, "a").unwrap(), ratio));
        let anchor = builder.port(spring, "b").unwrap();
        builder
            .boundary(anchor, Box::new(|_, _, values| values.fill(0.0)))
            .unwrap();
        spring_nodes.push(spring);
    }
    let push = builder
        .add_with_tangent(
            "push",
            Force::new(Waveform::constant(force_n).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();
    members.push((builder.port(push, "p").unwrap(), lever));
    let wall = builder.port(mass, "p").unwrap();
    builder.common_flow_ratios(wall, &members).unwrap();
    let mut system = builder.seal().unwrap();
    // The loop is affine, so one exact Newton update has to close it: the
    // ratios are in the interface Jacobian as well as in the residual.
    system.set_config(SolveConfig {
        max_iterations: 1,
        ..SolveConfig::default()
    });

    let reflected: f64 = springs.iter().map(|(k, r)| r * r * k).sum();
    let omega = (reflected / mass_kg).sqrt();
    let rest = lever * force_n / reflected;
    let energy0 = 0.5 * mass_kg * v0 * v0
        + springs
            .iter()
            .map(|(k, r)| 0.5 * k * (r * x0).powi(2))
            .sum::<f64>();
    let (dt, steps) = (0.01, 2_000);
    with_cx(|cx| {
        for i in 1..=steps {
            system.step(dt, cx).unwrap();
            let (sin, cos) = (f64::from(i) * 2.0 * (0.5 * omega * dt).atan()).sin_cos();
            let x = rest + (x0 - rest) * cos + v0 / omega * sin;
            let v = v0 * cos - omega * (x0 - rest) * sin;
            let motion = system.state::<Motion>(mass).unwrap();
            close(motion.position_m, x, 2.0e-12);
            close(motion.velocity_m_s, v, 2.0e-11);
            let mut energy = 0.5 * mass_kg * motion.velocity_m_s.powi(2);
            for (node, (stiffness, ratio)) in spring_nodes.iter().zip(springs) {
                let q = *system.state::<f64>(*node).unwrap();
                close(q, -ratio * x, 4.0e-12);
                energy += 0.5 * stiffness * q * q;
            }
            // The source's account holds the work it received: minus the
            // levered force times the mass's travel.
            let account = system.state::<SourceAccount>(push).unwrap().energy_in_j;
            close(
                account,
                -lever * force_n * (motion.position_m - x0),
                1.0e-12,
            );
            close(energy + account, energy0, 4.0e-12 * energy0);
        }
    });
    assert_eq!(system.ledger().connections.len(), 1);
    close(
        system.ledger().connections[0].energy_defect_j.value(),
        0.0,
        1.0e-11 * energy0,
    );
    close(system.ledger().energy_residual_j(), 0.0, 1.0e-11 * energy0);
}

#[test]
fn a_capacitor_feeding_two_coils_through_turns_ratios_rings_at_the_reflected_frequency() {
    // C e' = -(n1 f1 + n2 f2) and L_i f_i' = n_i e. The junction current
    // g = n1 f1 + n2 f2 and the voltage e are one oscillator with
    // w² = (n1²/L1 + n2²/L2) / C, and each coil follows g by its own share.
    let kind = PortKind::ElectricalVoltageCurrent;
    let (capacitance, charge0) = (2.0e-3, 0.004);
    let coils = [(0.5, 3.0, 0.2), (0.125, -0.5, -0.1)];
    let mut builder = SystemBuilder::new();
    let capacitor = builder
        .add_with_tangent(
            "capacitor",
            fs_system::elements::Storage::in_domain(kind, capacitance).unwrap(),
            charge0,
        )
        .unwrap();
    let mut members = Vec::new();
    let mut coil_nodes = Vec::new();
    for (name, (inductance, turns, current0)) in ["primary", "secondary"].into_iter().zip(coils) {
        let coil = builder
            .add_with_tangent(
                name,
                Mass::in_domain(kind, inductance).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s: current0,
                },
            )
            .unwrap();
        members.push((builder.port(coil, "p").unwrap(), turns));
        coil_nodes.push(coil);
    }
    let node = builder.port(capacitor, "p").unwrap();
    builder.common_effort_ratios(node, &members).unwrap();
    let mut system = builder.seal().unwrap();
    system.set_config(SolveConfig {
        max_iterations: 1,
        ..SolveConfig::default()
    });

    let share: f64 = coils.iter().map(|(l, n, _)| n * n / l).sum();
    let omega = (share / capacitance).sqrt();
    let e0 = charge0 / capacitance;
    let g0: f64 = coils.iter().map(|(_, n, f)| n * f).sum();
    let energy0 = 0.5 * charge0 * e0
        + coils
            .iter()
            .map(|(l, _, f)| 0.5 * l * f * f)
            .sum::<f64>();
    let (dt, steps) = (2.0e-4, 3_000);
    with_cx(|cx| {
        for i in 1..=steps {
            system.step(dt, cx).unwrap();
            let (sin, cos) = (f64::from(i) * 2.0 * (0.5 * omega * dt).atan()).sin_cos();
            let e = e0 * cos - g0 / (capacitance * omega) * sin;
            let g = g0 * cos + capacitance * omega * e0 * sin;
            let charge = *system.state::<f64>(capacitor).unwrap();
            close(charge / capacitance, e, 2.0e-11 * e0);
            let mut energy = 0.5 * charge * charge / capacitance;
            for (coil, (inductance, turns, current0)) in coil_nodes.iter().zip(coils) {
                let current = system.state::<Motion>(*coil).unwrap().velocity_m_s;
                close(
                    current,
                    current0 + turns / inductance / share * (g - g0),
                    2.0e-11,
                );
                energy += 0.5 * inductance * current * current;
            }
            close(energy, energy0, 2.0e-11 * energy0);
        }
    });
    close(
        system.ledger().connections[0].energy_defect_j.value(),
        0.0,
        1.0e-11 * energy0,
    );
    close(system.ledger().energy_residual_j(), 0.0, 1.0e-11 * energy0);
}

#[test]
fn refused_ratio_junctions_leave_every_port_usable() {
    let id = || StableId::new("p").unwrap();
    let mut builder = SystemBuilder::new();
    let wall = terminated(&mut builder, "wall", power(PortRole::EffortIn, 1));
    let near = terminated(&mut builder, "near", power(PortRole::FlowIn, 1));
    let far = terminated(&mut builder, "far", power(PortRole::FlowIn, 1));
    let tank = terminated(&mut builder, "tank", power(PortRole::FlowIn, 1));
    let pipe = terminated(&mut builder, "pipe", power(PortRole::EffortIn, 1));
    let drain = terminated(&mut builder, "drain", power(PortRole::EffortIn, 1));
    let wide = terminated(&mut builder, "wide", power(PortRole::FlowIn, 2));
    let hot = terminated(
        &mut builder,
        "hot",
        PortSpec::heat(id(), PortRole::EffortIn).unwrap(),
    );
    let cold = terminated(
        &mut builder,
        "cold",
        PortSpec::heat(id(), PortRole::FlowIn).unwrap(),
    );
    let signal = terminated(
        &mut builder,
        "signal",
        PortSpec::signal_in(id(), NonZeroUsize::MIN),
    );

    assert!(matches!(
        builder.common_flow_ratios(wall, &[]),
        Err(AdmissionError::Empty)
    ));
    assert!(matches!(
        builder.common_effort_ratios(tank, &[]),
        Err(AdmissionError::Empty)
    ));
    for ratio in [0.0, -0.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            builder.common_flow_ratios(wall, &[(near, 2.0), (far, ratio)]),
            Err(AdmissionError::TransformerRatio { ratio: rejected }) if rejected.to_bits() == ratio.to_bits()
        ));
        assert!(matches!(
            builder.common_effort_ratios(tank, &[(pipe, ratio), (drain, 2.0)]),
            Err(AdmissionError::TransformerRatio { ratio: rejected }) if rejected.to_bits() == ratio.to_bits()
        ));
    }
    // A common flow is set by a port that receives effort and shared by
    // ports that receive flow; heat and signals are not power pairs.
    for (setter, member) in [
        (near, far),
        (wall, pipe),
        (hot, near),
        (wall, cold),
        (wall, signal),
    ] {
        assert!(matches!(
            builder.common_flow_ratios(setter, &[(member, 2.0)]),
            Err(AdmissionError::Role { .. })
        ));
    }
    // The dual: set by a port that receives flow, shared by ports that
    // receive effort.
    for (setter, member) in [
        (pipe, drain),
        (tank, near),
        (cold, pipe),
        (tank, hot),
        (tank, signal),
    ] {
        assert!(matches!(
            builder.common_effort_ratios(setter, &[(member, 2.0)]),
            Err(AdmissionError::Role { .. })
        ));
    }
    assert!(matches!(
        builder.common_flow_ratios(wall, &[(near, 1.0), (wide, 1.0)]),
        Err(AdmissionError::Mismatch { .. })
    ));
    assert!(matches!(
        builder.common_flow_ratios(wall, &[(near, 1.0), (near, 2.0)]),
        Err(AdmissionError::AlreadyConnected { .. })
    ));

    // Every refusal above named these ports, and none was claimed.
    builder
        .common_flow_ratios(wall, &[(near, 2.0), (far, -3.0)])
        .unwrap();
    builder
        .common_effort_ratios(tank, &[(pipe, 0.5), (drain, -4.0)])
        .unwrap();
    assert!(matches!(
        builder.common_flow_ratios(wall, &[(near, 2.0)]),
        Err(AdmissionError::AlreadyConnected { .. })
    ));
    for (port, value) in [(wide, 0.0), (hot, 300.0), (cold, 0.0), (signal, 0.0)] {
        builder
            .boundary(port, Box::new(move |_, _, values| values.fill(value)))
            .unwrap();
    }
    let mut system = builder.seal().unwrap();
    with_cx(|cx| {
        system.step(0.01, cx).unwrap();
    });
    assert_eq!(system.ledger().connections.len(), 2);
    assert_eq!(system.ledger().energy_residual_j(), 0.0);
}
