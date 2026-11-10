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
