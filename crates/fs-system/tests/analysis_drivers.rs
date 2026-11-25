//! Direct verification of bead `frankensim-rc-1109-root-3gc12.5.10`:
//! analysis drivers over a sealed `System` — stationary (`StageResidual`),
//! linearised harmonic response & modal analysis, periodic steady state
//! (shooting + `fs-orbit` harmonic balance), forward-tangent & matrix-free IFT
//! adjoint sensitivities (`fs-adjoint::ift_gradient_matfree` + `verify_gradient`),
//! and strict capability gating naming the component and missing tier.

#![allow(clippy::float_cmp)]

use core::f64::consts::TAU;
use core::num::NonZeroUsize;

use fs_adjoint::verify_gradient;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_orbit::{HbAnchor, HbBudget, HbScales};
use fs_ports::component::{
    Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt, StageResidual,
};
use fs_ports::{PortKind, StableId};
use fs_system::elements::{
    Conductance, ControlledSource, Damper, HeatCapacity, Mass, Motion, SourceAccount, Spring,
};
use fs_system::{
    AnalysisError, AnalysisKind, AnalysisOutcome, AnalysisRequest, BoundaryChannel,
    LinearizeConfig, ModalConfig, NodeId, OutputChannel, PeriodicHbConfig, PeriodicShootingConfig,
    RequiredTier, RouteTemperature, SensitivityConfig, StateMap, StationaryConfig, StepObjective,
    System, SystemBuilder, TrajectoryStepPartials,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 2026,
                kernel_id: 510,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn id(s: &str) -> StableId {
    StableId::new(s).expect("valid id")
}

fn mech_port(name: &str, role: PortRole) -> PortSpec {
    PortSpec::power(
        id(name),
        PortKind::MechanicalForceVelocity,
        role,
        NonZeroUsize::MIN,
    )
    .expect("mech port")
}

fn roundoff(scale: f64) -> f64 {
    64.0 * f64::EPSILON * scale.max(1.0)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct StageErr;

/// Point mass implementing Component + PortTangent + StageResidual (exact midpoint).
struct StageMass {
    mass_kg: f64,
    ports: [PortSpec; 1],
}

impl StageMass {
    fn new(mass_kg: f64) -> Self {
        Self {
            mass_kg,
            ports: [mech_port("p", PortRole::EffortIn)],
        }
    }
}

impl Component for StageMass {
    type State = f64;
    type Error = StageErr;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, v: &f64) -> Result<f64, StageErr> {
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
    ) -> Result<Proposal, StageErr> {
        cx.checkpoint().map_err(|_| StageErr)?;
        let f = inputs[0];
        *new = old + dt_s * f / self.mass_kg;
        let v_mid = 0.5 * (old + *new);
        outputs[0] = v_mid;
        let de = 0.5 * self.mass_kg * (*new * *new - old * old);
        let scale = 0.5 * self.mass_kg * (old * old + *new * *new) + (dt_s * f * v_mid).abs();
        Ok(Proposal::Step(Receipt::lossless(de, roundoff(scale))))
    }
}

impl PortTangent for StageMass {
    fn jvp(
        &self,
        _old: &f64,
        dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), StageErr> {
        d_outputs[0] = 0.5 * dt_s / self.mass_kg * d_inputs[0];
        Ok(())
    }

    fn vjp(
        &self,
        _old: &f64,
        dt_s: f64,
        _inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), StageErr> {
        bar_inputs[0] = 0.5 * dt_s / self.mass_kg * bar_outputs[0];
        Ok(())
    }
}

impl StageResidual for StageMass {
    fn unknowns(&self) -> usize {
        1
    }

    fn initial_guess(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        z: &mut [f64],
    ) -> Result<(), StageErr> {
        z[0] = old + dt_s * inputs[0] / self.mass_kg;
        Ok(())
    }

    fn residual(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        r: &mut [f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), StageErr> {
        cx.checkpoint().map_err(|_| StageErr)?;
        r[0] = z[0] - old - dt_s * inputs[0] / self.mass_kg;
        outputs[0] = 0.5 * (old + z[0]);
        Ok(())
    }

    fn jvp_z(
        &self,
        _old: &f64,
        _dt_s: f64,
        _inputs: &[f64],
        _z: &[f64],
        dz: &[f64],
        dr: &mut [f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), StageErr> {
        dr[0] = dz[0];
        d_outputs[0] = 0.5 * dz[0];
        Ok(())
    }

    fn vjp_z(
        &self,
        _old: &f64,
        _dt_s: f64,
        _inputs: &[f64],
        _z: &[f64],
        bar_r: &[f64],
        bar_outputs: &[f64],
        bar_z: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), StageErr> {
        bar_z[0] = bar_r[0] + 0.5 * bar_outputs[0];
        Ok(())
    }

    fn accept(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        new: &mut f64,
    ) -> Result<Receipt, StageErr> {
        *new = z[0];
        let v_mid = 0.5 * (old + *new);
        let de = 0.5 * self.mass_kg * (*new * *new - old * old);
        let scale =
            0.5 * self.mass_kg * (old * old + *new * *new) + (dt_s * inputs[0] * v_mid).abs();
        Ok(Receipt::lossless(de, roundoff(scale)))
    }
}

/// Nonlinear cubic spring (`k x + c x^3`) with exact discrete-gradient
/// Component + PortTangent + StageResidual.
struct StageSpring {
    k: f64,
    cubic: f64,
    ports: [PortSpec; 2],
}

impl StageSpring {
    fn new(k: f64, cubic: f64) -> Self {
        Self {
            k,
            cubic,
            ports: [
                mech_port("a", PortRole::FlowIn),
                mech_port("b", PortRole::FlowIn),
            ],
        }
    }

    fn discrete_force(&self, a: f64, b: f64) -> f64 {
        0.5 * self.k * (a + b) + 0.25 * self.cubic * (a + b) * (a * a + b * b)
    }

    fn discrete_slope(&self, a: f64, b: f64) -> f64 {
        0.5 * self.k + 0.25 * self.cubic * (a * a + 2.0 * a * b + 3.0 * b * b)
    }
}

impl Component for StageSpring {
    type State = f64;
    type Error = StageErr;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        if self.cubic == 0.0 {
            Class::ExactLinear
        } else {
            Class::Smooth
        }
    }

    fn stored_energy_j(&self, x: &f64) -> Result<f64, StageErr> {
        Ok(0.5 * self.k * x * x + 0.25 * self.cubic * x * x * x * x)
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, StageErr> {
        cx.checkpoint().map_err(|_| StageErr)?;
        *new = old + dt_s * (inputs[0] - inputs[1]);
        let f = self.discrete_force(*old, *new);
        outputs[0] = f;
        outputs[1] = -f;
        let de = f * (*new - old);
        let scale = self.stored_energy_j(old)? + self.stored_energy_j(new)? + de.abs();
        Ok(Proposal::Step(Receipt::lossless(de, roundoff(scale))))
    }
}

impl PortTangent for StageSpring {
    fn jvp(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), StageErr> {
        let new = old + dt_s * (inputs[0] - inputs[1]);
        let df = self.discrete_slope(*old, new) * dt_s * (d_inputs[0] - d_inputs[1]);
        d_outputs[0] = df;
        d_outputs[1] = -df;
        Ok(())
    }

    fn vjp(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), StageErr> {
        self.jvp(old, dt_s, inputs, bar_outputs, bar_inputs, cx)
    }
}

impl StageResidual for StageSpring {
    fn unknowns(&self) -> usize {
        1
    }

    fn initial_guess(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        z: &mut [f64],
    ) -> Result<(), StageErr> {
        z[0] = old + dt_s * (inputs[0] - inputs[1]);
        Ok(())
    }

    fn residual(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        r: &mut [f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), StageErr> {
        cx.checkpoint().map_err(|_| StageErr)?;
        r[0] = z[0] - old - dt_s * (inputs[0] - inputs[1]);
        let f = self.discrete_force(*old, z[0]);
        outputs[0] = f;
        outputs[1] = -f;
        Ok(())
    }

    fn jvp_z(
        &self,
        old: &f64,
        _dt_s: f64,
        _inputs: &[f64],
        z: &[f64],
        dz: &[f64],
        dr: &mut [f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), StageErr> {
        dr[0] = dz[0];
        let df = self.discrete_slope(*old, z[0]) * dz[0];
        d_outputs[0] = df;
        d_outputs[1] = -df;
        Ok(())
    }

    fn vjp_z(
        &self,
        old: &f64,
        _dt_s: f64,
        _inputs: &[f64],
        z: &[f64],
        bar_r: &[f64],
        bar_outputs: &[f64],
        bar_z: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), StageErr> {
        let s = self.discrete_slope(*old, z[0]);
        bar_z[0] = bar_r[0] + s * (bar_outputs[0] - bar_outputs[1]);
        Ok(())
    }

    fn accept(
        &self,
        old: &f64,
        _dt_s: f64,
        _inputs: &[f64],
        z: &[f64],
        new: &mut f64,
    ) -> Result<Receipt, StageErr> {
        *new = z[0];
        let f = self.discrete_force(*old, *new);
        let de = f * (*new - old);
        let scale = self.stored_energy_j(old)? + self.stored_energy_j(new)? + de.abs();
        Ok(Receipt::lossless(de, roundoff(scale)))
    }
}

/// State map for `[extension_m, velocity_m_s]` on `(StageSpring, StageMass)`.
struct StageOscMap {
    spring: NodeId,
    mass: NodeId,
}

impl StateMap for StageOscMap {
    fn dim(&self) -> usize {
        2
    }

    fn read(&self, system: &System, out: &mut [f64]) -> Result<(), AnalysisError> {
        out[0] = *system
            .state::<f64>(self.spring)
            .ok_or(AnalysisError::Invalid {
                what: "missing StageSpring state",
            })?;
        out[1] = *system
            .state::<f64>(self.mass)
            .ok_or(AnalysisError::Invalid {
                what: "missing StageMass state",
            })?;
        Ok(())
    }

    fn write(&self, system: &mut System, state: &[f64]) -> Result<(), AnalysisError> {
        system.set_accepted_state(self.spring, state[0])?;
        system.set_accepted_state(self.mass, state[1])?;
        Ok(())
    }
}

/// State map for `[spring_x, mass_v, body_T]` on `(Spring, Mass, HeatCapacity)`.
struct ThermoMechMap {
    spring: NodeId,
    mass: NodeId,
    body: NodeId,
    damper: NodeId,
}

impl StateMap for ThermoMechMap {
    fn dim(&self) -> usize {
        3
    }

    fn read(&self, system: &System, out: &mut [f64]) -> Result<(), AnalysisError> {
        out[0] = *system
            .state::<f64>(self.spring)
            .ok_or(AnalysisError::Invalid {
                what: "missing Spring state",
            })?;
        out[1] = system
            .state::<Motion>(self.mass)
            .ok_or(AnalysisError::Invalid {
                what: "missing Mass state",
            })?
            .velocity_m_s;
        out[2] = *system
            .state::<f64>(self.body)
            .ok_or(AnalysisError::Invalid {
                what: "missing HeatCapacity state",
            })?;
        Ok(())
    }

    fn write(&self, system: &mut System, state: &[f64]) -> Result<(), AnalysisError> {
        if !(state[2].is_finite() && state[2] > 0.0) {
            return Err(AnalysisError::Invalid {
                what: "temperature state must stay positive",
            });
        }
        system.set_accepted_state(self.spring, state[0])?;
        let mut motion = *system
            .state::<Motion>(self.mass)
            .ok_or(AnalysisError::Invalid {
                what: "missing Mass state",
            })?;
        motion.velocity_m_s = state[1];
        system.set_accepted_state(self.mass, motion)?;
        system.set_accepted_state(self.body, state[2])?;
        Ok(())
    }
}

/// State map for mechanical-only `[spring_x, mass_v]` (2D).
struct Mech2Map {
    spring: NodeId,
    mass: NodeId,
}

impl StateMap for Mech2Map {
    fn dim(&self) -> usize {
        2
    }

    fn read(&self, system: &System, out: &mut [f64]) -> Result<(), AnalysisError> {
        out[0] = *system
            .state::<f64>(self.spring)
            .ok_or(AnalysisError::Invalid {
                what: "missing Spring state",
            })?;
        out[1] = system
            .state::<Motion>(self.mass)
            .ok_or(AnalysisError::Invalid {
                what: "missing Mass state",
            })?
            .velocity_m_s;
        Ok(())
    }

    fn write(&self, system: &mut System, state: &[f64]) -> Result<(), AnalysisError> {
        system.set_accepted_state(self.spring, state[0])?;
        let mut motion = *system
            .state::<Motion>(self.mass)
            .ok_or(AnalysisError::Invalid {
                what: "missing Mass state",
            })?;
        motion.velocity_m_s = state[1];
        system.set_accepted_state(self.mass, motion)?;
        Ok(())
    }
}

/// Build a composed mass-spring-damper-thermal system with external force
/// boundary and ambient temperature boundary.
#[allow(clippy::too_many_arguments)]
fn build_thermo_mech_oscillator(
    mass_kg: f64,
    stiffness_n_m: f64,
    cubic_n_m3: f64,
    damping_n_s_m: f64,
    heat_cap_j_k: f64,
    conductance_w_k: f64,
    force_fn: impl Fn(f64, f64) -> f64 + 'static,
    t_amb_k: f64,
    initial_x: f64,
    initial_v: f64,
    initial_t: f64,
) -> (
    System,
    ThermoMechMap,
    fs_system::PortRef,
    fs_system::PortRef,
) {
    let mut b = SystemBuilder::new();
    let mass = b
        .add_with_tangent(
            "mass",
            Mass::new(mass_kg).expect("mass"),
            Motion {
                position_m: 0.0,
                velocity_m_s: initial_v,
            },
        )
        .expect("add mass");
    let spring = b
        .add_with_tangent(
            "spring",
            Spring::new(stiffness_n_m, cubic_n_m3).expect("spring"),
            initial_x,
        )
        .expect("add spring");
    let damper = b
        .add_with_tangent("damper", Damper::new(damping_n_s_m).expect("damper"), ())
        .expect("add damper");
    let force_src = b
        .add_with_tangent(
            "ext_force",
            ControlledSource::effort(PortKind::MechanicalForceVelocity),
            SourceAccount::default(),
        )
        .expect("add ext_force");
    let body = b
        .add_with_tangent(
            "body",
            HeatCapacity::new(heat_cap_j_k).expect("heat cap"),
            initial_t,
        )
        .expect("add body");
    let cond = b
        .add_with_tangent("cond", Conductance::new(conductance_w_k).expect("cond"), ())
        .expect("add cond");

    let m_p = b.port(mass, "p").unwrap();
    let s_a = b.port(spring, "a").unwrap();
    let s_b = b.port(spring, "b").unwrap();
    let d_a = b.port(damper, "a").unwrap();
    let d_b = b.port(damper, "b").unwrap();
    let d_h = b.port(damper, "heat").unwrap();
    let f_p = b.port(force_src, "p").unwrap();
    let f_cmd = b.port(force_src, "command").unwrap();
    let b_h = b.port(body, "heat").unwrap();
    let c_a = b.port(cond, "a").unwrap();
    let c_b = b.port(cond, "b").unwrap();

    // Common flow on mass: mass returns v, receives F_spring + F_damper + F_ext.
    // Spring & damper receive -v at port `a` and 0 at anchored port `b`, so
    // extension grows as (-v - 0) = -v? To make x_dot = +v (positive extension
    // with positive velocity), put port `b` on the mass and anchor port `a` at 0!
    // Check Spring signs: new = old + dt * (inputs[0] - inputs[1]).
    // With `b` on mass (`inputs[1] = -v`) and `a` anchored (`inputs[0] = 0`):
    // new = old + dt * (0 - (-v)) = old + dt * v!
    // And output[1] = -force_n = -k*x, which goes into mass: m dv/dt = -k*x - c*v + F_ext!
    b.common_flow(m_p, &[s_b, d_b, f_p]).unwrap();
    b.boundary(s_a, Box::new(|_, _, out| out[0] = 0.0)).unwrap();
    b.boundary(d_a, Box::new(|_, _, out| out[0] = 0.0)).unwrap();
    b.boundary(
        f_cmd,
        Box::new(move |t0, dt, out| out[0] = force_fn(t0, dt)),
    )
    .unwrap();

    // Thermal loop: damper heat -> body heat via common effort (temperature) + conductance to ambient.
    b.common_effort(b_h, &[d_h, c_a]).unwrap();
    b.boundary(c_b, Box::new(move |_, _, out| out[0] = t_amb_k))
        .unwrap();

    let mut sys = b.seal().expect("seal thermo-mech");
    sys.guess(d_h, &[initial_t]);
    sys.guess(c_a, &[initial_t]);
    (
        sys,
        ThermoMechMap {
            spring,
            mass,
            body,
            damper,
        },
        f_cmd,
        c_b,
    )
}

#[test]
#[allow(clippy::too_many_lines)]
fn capability_gating_names_component_and_missing_tier() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        // 1. Plain Mass (mandatory Component only: no PortTangent, no StageResidual).
        let mut b = SystemBuilder::new();
        let m = b
            .add(
                "plain_rotor",
                Mass::new(2.0).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s: 1.0,
                },
            )
            .unwrap();
        let s = b
            .add_with_tangent("tangent_spring", Spring::new(10.0, 0.0).unwrap(), 0.1)
            .unwrap();
        let m_p = b.port(m, "p").unwrap();
        let s_a = b.port(s, "a").unwrap();
        let s_b = b.port(s, "b").unwrap();
        b.common_flow(m_p, &[s_b]).unwrap();
        b.boundary(s_a, Box::new(|_, _, out| out[0] = 0.0)).unwrap();
        let mut sys = b.seal().unwrap();
        let map = Mech2Map { spring: s, mass: m };

        // Every PortTangent-dependent driver must refuse naming `plain_rotor` and `PortTangent`.
        for kind in [
            AnalysisKind::Harmonic,
            AnalysisKind::Modal,
            AnalysisKind::PeriodicShooting,
            AnalysisKind::PeriodicHarmonicBalance,
            AnalysisKind::TangentSensitivity,
            AnalysisKind::AdjointSensitivity,
        ] {
            let err = sys
                .check_analysis_capabilities(kind)
                .expect_err("must refuse missing PortTangent");
            match &err {
                AnalysisError::MissingTier {
                    analysis,
                    node,
                    required,
                } => {
                    assert_eq!(*analysis, kind);
                    assert_eq!(node, "plain_rotor");
                    assert_eq!(*required, RequiredTier::PortTangent);
                }
                other => panic!("expected MissingTier, got {other:?}"),
            }
            let msg = err.to_string();
            assert!(
                msg.contains("plain_rotor"),
                "message missing node name: {msg}"
            );
            assert!(
                msg.contains("PortTangent"),
                "message missing required tier: {msg}"
            );
        }

        // Calling the actual drivers returns the exact same typed refusal.
        let modal_err = sys
            .modal_analysis(&map, ModalConfig::default(), cx)
            .expect_err("modal must refuse plain_rotor");
        assert!(modal_err.to_string().contains("plain_rotor"));
        assert!(modal_err.to_string().contains("PortTangent"));

        let adj_err = sys
            .step_sensitivity_adjoint(
                1.0e-3,
                &StepObjective {
                    dj_du: Vec::new(),
                    dj_dy: Vec::new(),
                    dj_db: Vec::new(),
                },
                SensitivityConfig::default(),
                cx,
            )
            .expect_err("adjoint must refuse plain_rotor");
        assert!(adj_err.to_string().contains("plain_rotor"));
        assert!(adj_err.to_string().contains("PortTangent"));

        // 2. Tangent-only system (has PortTangent, lacks StageResidual) passes Modal
        //    but refuses Stationary naming the first non-StageResidual component & `StageResidual`.
        let mut b2 = SystemBuilder::new();
        let m2 = b2
            .add_with_tangent(
                "tangent_only_mass",
                Mass::new(2.0).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s: 0.0,
                },
            )
            .unwrap();
        let s2 = b2
            .add_with_tangent_and_stage_residual("stage_spring", StageSpring::new(10.0, 0.0), 0.0)
            .unwrap();
        let m2_p = b2.port(m2, "p").unwrap();
        let s2_a = b2.port(s2, "a").unwrap();
        let s2_b = b2.port(s2, "b").unwrap();
        b2.common_flow(m2_p, &[s2_b]).unwrap();
        b2.boundary(s2_a, Box::new(|_, _, out| out[0] = 0.0))
            .unwrap();
        let mut sys2 = b2.seal().unwrap();
        let map2 = Mech2Map {
            spring: s2,
            mass: m2,
        };

        let stat_err = sys2
            .solve_stationary(&map2, StationaryConfig::default(), cx)
            .expect_err("stationary must refuse tangent_only_mass");
        match &stat_err {
            AnalysisError::MissingTier {
                analysis,
                node,
                required,
            } => {
                assert_eq!(*analysis, AnalysisKind::Stationary);
                assert_eq!(node, "tangent_only_mass");
                assert_eq!(*required, RequiredTier::StageResidual);
            }
            other => panic!("expected MissingTier(StageResidual), got {other:?}"),
        }
        let msg = stat_err.to_string();
        assert!(msg.contains("tangent_only_mass"), "{msg}");
        assert!(msg.contains("StageResidual"), "{msg}");
    });
}

#[test]
#[allow(clippy::too_many_lines, clippy::items_after_statements)]
fn stationary_equilibrium_stage_residual_and_unified_analysis_switch() {
    // Undamped nonlinear cubic oscillator (k x + c x^3 = F_ext, v = 0) built
    // with StageMass + StageSpring (both PortTangent + StageResidual).
    // Relaxation-to-rest would oscillate forever; stationary Newton-Krylov
    // converges directly to the exact root, and then `run_analysis` switches
    // to Modal, Harmonic, Transient, and Adjoint on the UNCHANGED assembly.
    let mass_kg = 2.0;
    let k = 25.0;
    let cubic = 40.0;
    let f_ext = 15.0;

    // Conservative energy-reservoir actuator (`State = reservoir_energy_j`,
    // `stored_energy_j = E`, exact lossless work balance) offering
    // `Component + PortTangent + StageResidual` (`unknowns() == 1`).
    struct StageActuator {
        ports: [PortSpec; 2],
    }
    impl StageActuator {
        fn new() -> Self {
            Self {
                ports: [
                    mech_port("p", PortRole::FlowIn),
                    PortSpec::signal_in(id("cmd"), NonZeroUsize::MIN),
                ],
            }
        }
    }
    impl Component for StageActuator {
        type State = f64; // reservoir energy [J]
        type Error = StageErr;
        fn ports(&self) -> &[PortSpec] {
            &self.ports
        }
        fn class(&self) -> Class {
            Class::ExactLinear
        }
        fn stored_energy_j(&self, e: &f64) -> Result<f64, StageErr> {
            Ok(*e)
        }
        fn propose(
            &self,
            old: &f64,
            dt_s: f64,
            inputs: &[f64],
            new: &mut f64,
            outputs: &mut [f64],
            cx: &Cx<'_>,
        ) -> Result<Proposal, StageErr> {
            cx.checkpoint().map_err(|_| StageErr)?;
            let flow = inputs[0];
            let force = inputs[1];
            outputs[0] = force;
            let de = dt_s * force * flow;
            *new = old + de;
            let scale = old.abs() + new.abs() + de.abs();
            Ok(Proposal::Step(Receipt::lossless(de, roundoff(scale))))
        }
    }
    impl PortTangent for StageActuator {
        fn jvp(
            &self,
            _old: &f64,
            _dt_s: f64,
            _inputs: &[f64],
            d_inputs: &[f64],
            d_outputs: &mut [f64],
            _cx: &Cx<'_>,
        ) -> Result<(), StageErr> {
            d_outputs[0] = d_inputs[1];
            Ok(())
        }
        fn vjp(
            &self,
            _old: &f64,
            _dt_s: f64,
            _inputs: &[f64],
            bar_outputs: &[f64],
            bar_inputs: &mut [f64],
            _cx: &Cx<'_>,
        ) -> Result<(), StageErr> {
            bar_inputs[0] = 0.0;
            bar_inputs[1] = bar_outputs[0];
            Ok(())
        }
    }
    impl StageResidual for StageActuator {
        fn unknowns(&self) -> usize {
            1
        }
        fn initial_guess(
            &self,
            _old: &f64,
            _dt_s: f64,
            inputs: &[f64],
            z: &mut [f64],
        ) -> Result<(), StageErr> {
            z[0] = inputs[1];
            Ok(())
        }
        fn residual(
            &self,
            _old: &f64,
            _dt_s: f64,
            inputs: &[f64],
            z: &[f64],
            r: &mut [f64],
            outputs: &mut [f64],
            cx: &Cx<'_>,
        ) -> Result<(), StageErr> {
            cx.checkpoint().map_err(|_| StageErr)?;
            r[0] = z[0] - inputs[1];
            outputs[0] = z[0];
            Ok(())
        }
        fn jvp_z(
            &self,
            _old: &f64,
            _dt_s: f64,
            _inputs: &[f64],
            _z: &[f64],
            dz: &[f64],
            dr: &mut [f64],
            d_outputs: &mut [f64],
            _cx: &Cx<'_>,
        ) -> Result<(), StageErr> {
            dr[0] = dz[0];
            d_outputs[0] = dz[0];
            Ok(())
        }
        fn vjp_z(
            &self,
            _old: &f64,
            _dt_s: f64,
            _inputs: &[f64],
            _z: &[f64],
            bar_r: &[f64],
            bar_outputs: &[f64],
            bar_z: &mut [f64],
            _cx: &Cx<'_>,
        ) -> Result<(), StageErr> {
            bar_z[0] = bar_r[0] + bar_outputs[0];
            Ok(())
        }
        fn accept(
            &self,
            old: &f64,
            dt_s: f64,
            inputs: &[f64],
            z: &[f64],
            new: &mut f64,
        ) -> Result<Receipt, StageErr> {
            let de = dt_s * z[0] * inputs[0];
            *new = old + de;
            let scale = old.abs() + new.abs() + de.abs();
            Ok(Receipt::lossless(de, roundoff(scale)))
        }
    }

    let mut b = SystemBuilder::new();
    let mass = b
        .add_with_tangent_and_stage_residual("mass", StageMass::new(mass_kg), 0.3)
        .unwrap();
    let spring = b
        .add_with_tangent_and_stage_residual("spring", StageSpring::new(k, cubic), 0.05)
        .unwrap();
    let act = b
        .add_with_tangent_and_stage_residual("actuator", StageActuator::new(), 1_000.0)
        .unwrap();

    let m_p = b.port(mass, "p").unwrap();
    let s_a = b.port(spring, "a").unwrap();
    let s_b = b.port(spring, "b").unwrap();
    let a_p = b.port(act, "p").unwrap();
    let a_cmd = b.port(act, "cmd").unwrap();

    b.common_flow(m_p, &[s_b, a_p]).unwrap();
    b.boundary(s_a, Box::new(|_, _, out| out[0] = 0.0)).unwrap();
    b.boundary(a_cmd, Box::new(move |_, _, out| out[0] = f_ext))
        .unwrap();

    let mut sys = b.seal().unwrap();
    let map = StageOscMap { spring, mass };
    let gate = CancelGate::new();

    with_cx(&gate, |cx| {
        // 1. Stationary solve with commit = false leaves accepted state untouched.
        let outcome = sys
            .run_analysis(
                AnalysisRequest::Stationary {
                    state_map: &map,
                    config: StationaryConfig {
                        probe_dt_s: 1.0e-2,
                        absolute_tolerance: 1.0e-11,
                        commit: false,
                        ..StationaryConfig::default()
                    },
                },
                cx,
            )
            .expect("stationary solve");
        let stat = match outcome {
            AnalysisOutcome::Stationary(s) => s,
            other => panic!("expected Stationary outcome, got {other:?}"),
        };
        let x_eq = stat.state[0];
        let v_eq = stat.state[1];
        assert!(v_eq.abs() < 1.0e-11, "v_eq = {v_eq}");
        let balance = k * x_eq + cubic * x_eq.powi(3) - f_ext;
        assert!(balance.abs() < 1.0e-10, "nonlinear balance = {balance}");
        // Accepted state remained at initial (0.05, 0.3) because commit = false.
        assert!((*sys.state::<f64>(spring).unwrap() - 0.05).abs() < 1.0e-15);
        assert!((*sys.state::<f64>(mass).unwrap() - 0.3).abs() < 1.0e-15);

        // 2. Stationary solve with commit = true installs x_eq into the SAME System.
        let _ = sys
            .run_analysis(
                AnalysisRequest::Stationary {
                    state_map: &map,
                    config: StationaryConfig {
                        probe_dt_s: 1.0e-2,
                        absolute_tolerance: 1.0e-11,
                        commit: true,
                        ..StationaryConfig::default()
                    },
                },
                cx,
            )
            .expect("stationary commit");
        assert!((*sys.state::<f64>(spring).unwrap() - x_eq).abs() < 1.0e-12);

        // 3. Switch directly to Modal on the UNCHANGED System about the nonlinear
        //    equilibrium: tangent stiffness is k_eff = k + 3 * cubic * x_eq^2!
        let k_eff = k + 3.0 * cubic * x_eq * x_eq;
        let omega_expected = (k_eff / mass_kg).sqrt();
        let modal_out = sys
            .run_analysis(
                AnalysisRequest::Modal {
                    state_map: &map,
                    config: ModalConfig {
                        linearize: LinearizeConfig {
                            dt_s: 1.0e-3,
                            fd_eps: 1.0e-6,
                        },
                        max_qr_sweeps: 160,
                    },
                },
                cx,
            )
            .expect("modal at nonlinear equilibrium");
        let modal = match modal_out {
            AnalysisOutcome::Modal(m) => m,
            other => panic!("expected Modal, got {other:?}"),
        };
        assert_eq!(modal.modes.len(), 2);
        for mode in &modal.modes {
            assert!(
                (mode.natural_frequency_rad_s - omega_expected).abs() < 1.0e-6 * omega_expected,
                "got {}, expected {}",
                mode.natural_frequency_rad_s,
                omega_expected
            );
            assert!(
                mode.decay_rate_per_s.abs() < 1.0e-6,
                "lossless decay {}",
                mode.decay_rate_per_s
            );
        }

        // 4. Switch to Harmonic, AdjointSensitivity, and Transient via the SAME switch!
        let in_ch = [BoundaryChannel {
            port: a_cmd,
            component: 0,
        }];
        let out_ch = [OutputChannel::State(0)];
        let harm_out = sys
            .run_analysis(
                AnalysisRequest::Harmonic {
                    state_map: &map,
                    inputs: &in_ch,
                    outputs: &out_ch,
                    frequencies_rad_s: &[0.0],
                    config: LinearizeConfig {
                        dt_s: 1.0e-3,
                        fd_eps: 1.0e-6,
                    },
                },
                cx,
            )
            .expect("DC compliance at equilibrium");
        let harm = match harm_out {
            AnalysisOutcome::Harmonic(h) => h,
            other => panic!("expected Harmonic, got {other:?}"),
        };
        // DC compliance dx/dF = 1 / k_eff!
        let dc_gain = harm.continuous(0, 0, 0);
        assert!(
            (dc_gain.re - 1.0 / k_eff).abs() < 1.0e-7,
            "DC compliance {} vs {}",
            dc_gain.re,
            1.0 / k_eff
        );
        assert!(dc_gain.im.abs() < 1.0e-9);

        let trans_out = sys
            .run_analysis(
                AnalysisRequest::Transient {
                    dt_s: 1.0e-3,
                    steps: 5,
                },
                cx,
            )
            .expect("transient at equilibrium stays stationary");
        match trans_out {
            AnalysisOutcome::Transient(reps) => assert_eq!(reps.len(), 5),
            other => panic!("expected Transient, got {other:?}"),
        }
        assert!((*sys.state::<f64>(spring).unwrap() - x_eq).abs() < 1.0e-11);
        assert!(sys.state::<f64>(mass).unwrap().abs() < 1.0e-11);
    });

    // 5. Cancellation check: pre-cancelled gate refuses cleanly without mutating state.
    let cancelled = CancelGate::new();
    cancelled.request();
    with_cx(&cancelled, |cx| {
        let before_x = *sys.state::<f64>(spring).unwrap();
        let err = sys
            .solve_stationary(&map, StationaryConfig::default(), cx)
            .expect_err("cancelled stationary must refuse");
        assert!(matches!(err, AnalysisError::Cancelled));
        assert_eq!(
            sys.state::<f64>(spring).unwrap().to_bits(),
            before_x.to_bits()
        );
    });
}

#[test]
#[allow(clippy::too_many_lines)]
fn harmonic_response_and_modal_match_analytic_and_transient_steady_state() {
    let m = 1.5;
    let k = 60.0;
    let c = 1.8;
    let c_th = 200.0;
    let g_th = 10.0;
    let t_amb = 300.0;
    let f_amp = 4.0;
    let dt_step = 2.0e-3;

    // 1. Modal & Harmonic about zero-force rest state.
    let (mut sys_lin, map, f_cmd, _c_b) =
        build_thermo_mech_oscillator(m, k, 0.0, c, c_th, g_th, |_, _| 0.0, t_amb, 0.0, 0.0, t_amb);
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let modal = sys_lin
            .modal_analysis(
                &map,
                ModalConfig {
                    linearize: LinearizeConfig {
                        dt_s: dt_step,
                        fd_eps: 1.0e-5,
                    },
                    max_qr_sweeps: 160,
                },
                cx,
            )
            .expect("modal analysis");
        assert_eq!(modal.modes.len(), 3);
        // Mode 0: real thermal pole at lambda = -G / C_th = -0.05 rad/s.
        let thermal_pole = -g_th / c_th;
        assert!(
            (modal.modes[0].eigenvalue.re - thermal_pole).abs() < 1.0e-8,
            "thermal pole {:?} vs {thermal_pole}",
            modal.modes[0].eigenvalue
        );
        assert!(modal.modes[0].eigenvalue.im.abs() < 1.0e-10);
        // Modes 1 & 2: conjugate mechanical pair -c/(2m) +/- i sqrt(k/m - (c/(2m))^2).
        let sigma_ex = -c / (2.0 * m);
        let omega_d_ex = (k / m - sigma_ex * sigma_ex).sqrt();
        for mode in &modal.modes[1..3] {
            assert!(
                (mode.eigenvalue.re - sigma_ex).abs() < 1.0e-8,
                "sigma {} vs {sigma_ex}",
                mode.eigenvalue.re
            );
            assert!(
                (mode.angular_frequency_rad_s - omega_d_ex).abs() < 1.0e-8,
                "omega_d {} vs {omega_d_ex}",
                mode.angular_frequency_rad_s
            );
        }

        // Harmonic response at 3 frequencies: sub-resonant, near-resonant, super-resonant.
        let omegas = [3.5_f64, 6.2_f64, 11.0_f64];
        let in_ch = [BoundaryChannel {
            port: f_cmd,
            component: 0,
        }];
        let mass_node = map.mass;
        let mass_port = sys_lin.port(mass_node, "p").unwrap();
        let out_ch = [
            OutputChannel::State(0), // displacement x
            OutputChannel::PortOutput {
                port: mass_port,
                component: 0,
            }, // step-mean velocity v_mid
        ];
        let harm = sys_lin
            .harmonic_response(
                &map,
                &in_ch,
                &out_ch,
                &omegas,
                LinearizeConfig {
                    dt_s: dt_step,
                    fd_eps: 1.0e-6,
                },
                cx,
            )
            .expect("harmonic response");

        // Compare continuous H(i omega) against exact analytic formula:
        // H_x(i w) = 1 / (k - m w^2 + i c w),  H_v(i w) = i w * H_x(i w).
        for (idx, &w) in omegas.iter().enumerate() {
            let denom = C64::new(k - m * w * w, c * w);
            let hx_analytic = C64::ONE / denom;
            let hv_analytic = C64::new(0.0, w) * hx_analytic;
            let hx_got = harm.continuous(idx, 0, 0);
            let hv_got = harm.continuous(idx, 1, 0);
            assert!(
                (hx_got - hx_analytic).abs() < 1.0e-9 * hx_analytic.abs(),
                "w={w}: H_x got {hx_got:?}, want {hx_analytic:?}"
            );
            assert!(
                (hv_got - hv_analytic).abs() < 1.0e-9 * hv_analytic.abs(),
                "w={w}: H_v got {hv_got:?}, want {hv_analytic:?}"
            );
        }

        // Now verify discrete H_d(e^{i omega dt}) directly against transient
        // sinusoidal steady-state marching on the EXACT SAME linear graph!
        for (idx, &w) in omegas.iter().enumerate() {
            let period = TAU / w;
            let steps_per_period = (period / dt_step).round() as usize;
            let dt_exact = period / steps_per_period as f64;

            let (mut sys_t, map_t, f_cmd_t, _) = build_thermo_mech_oscillator(
                m,
                k,
                0.0,
                c,
                c_th,
                g_th,
                move |t0, dt| f_amp * (w * (t0 + 0.5 * dt)).cos(),
                t_amb,
                0.0,
                0.0,
                t_amb,
            );
            // Re-evaluate harmonic response at the exact integer-period step `dt_exact`.
            let m_port_t = sys_t.port(map_t.mass, "p").unwrap();
            let harm_exact = sys_t
                .harmonic_response(
                    &map_t,
                    &[BoundaryChannel {
                        port: f_cmd_t,
                        component: 0,
                    }],
                    &[OutputChannel::PortOutput {
                        port: m_port_t,
                        component: 0,
                    }],
                    &[w],
                    LinearizeConfig {
                        dt_s: dt_exact,
                        fd_eps: 1.0e-6,
                    },
                    cx,
                )
                .expect("harmonic at dt_exact");

            // March 25 periods (> 12 damping time constants: e^{-0.6 * 25 * T} < 1e-7)
            // then demodulate the fundamental phasor of step-mean velocity over 1 period.
            let warmup_periods = ((20.0 / (-sigma_ex)) / period).ceil() as usize;
            for _ in 0..(warmup_periods * steps_per_period) {
                sys_t.step(dt_exact, cx).expect("transient warmup");
            }
            let mut re_sum = 0.0_f64;
            let mut im_sum = 0.0_f64;
            for step_i in 0..steps_per_period {
                let t0 = sys_t.time_s();
                sys_t.step(dt_exact, cx).expect("transient sample");
                let v_mid = sys_t.last_output(m_port_t)[0];
                let phase = w * (t0 + 0.5 * dt_exact);
                // For v(t) = Re(V e^{i w t}) = |V| cos(w t + phi):
                // (2/N) sum v_k e^{-i w t_k} = V.
                re_sum += v_mid * phase.cos();
                im_sum -= v_mid * phase.sin();
                let _ = step_i;
            }
            let scale = 2.0 / (steps_per_period as f64 * f_amp);
            let h_transient = C64::new(re_sum * scale, im_sum * scale);
            let h_discrete = harm_exact.discrete(0, 0, 0);
            let err = (h_transient - h_discrete).abs();
            assert!(
                err < 1.0e-6 * h_discrete.abs(),
                "freq index {idx} (w={w}): transient {h_transient:?} vs H_d {h_discrete:?} (err {err:.3e})"
            );
        }
    });
}

#[test]
fn periodic_steady_state_shooting_and_fs_orbit_hb_match_on_nonlinear_duffing_system() {
    // Nonlinear Duffing oscillator (cubic stiffness) driven by sinusoidal force:
    // Both `solve_periodic_shooting` and `solve_periodic_hb` (through `fs-orbit`)
    // run directly on the EXACT SAME `System` and agree on fundamental + 3rd harmonic!
    let m = 1.0;
    let k = 16.0;
    let cubic = 30.0;
    let c = 1.2;
    let w = 4.5;
    let period = TAU / w;
    let f_amp = 3.0;

    let (mut sys, map_tm, _f_cmd, _c_b) = build_thermo_mech_oscillator(
        m,
        k,
        cubic,
        c,
        1_000.0,
        10.0,
        move |t0, dt| f_amp * (w * (t0 + 0.5 * dt)).cos(),
        300.0,
        0.0,
        0.0,
        300.0,
    );
    // Mechanical state map `[x, v]` (2D). Note thermal body warms minimally; fix
    // ambient conductance or solve mechanical 2D periodic orbit? To keep the full
    // 3D state strictly periodic (including thermal equilibrium T* > T_amb!),
    // let's solve on the FULL 3D `[x, v, T]` state map!
    // Check: in steady state, damper dissipates average power P_avg > 0, so body
    // reaches periodic temperature `T(t) = T_amb + P_avg / G + 2*omega ripple`!
    // Both shooting AND fs-orbit HB capture BOTH mechanical 3rd harmonic AND
    // thermal 2nd harmonic + DC temperature rise!
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let shoot = sys
            .solve_periodic_shooting(
                &map_tm,
                None,
                PeriodicShootingConfig {
                    period_s: period,
                    steps_per_period: 120,
                    tolerance: 1.0e-10,
                    max_iterations: 20,
                    fd_eps: 1.0e-6,
                    max_qr_sweeps: 160,
                    commit: false,
                },
                cx,
            )
            .expect("periodic shooting on nonlinear thermo-mech Duffing system");
        assert!(shoot.closure_residual < 1.0e-10);
        assert!(
            shoot.spectral_radius < 1.0,
            "damped orbit must be Floquet-stable, rho = {}",
            shoot.spectral_radius
        );
        // Thermal body must have heated above 300 K from damper dissipation!
        assert!(
            shoot.initial_state[2] > 300.01,
            "expected steady thermal rise above 300 K, got {}",
            shoot.initial_state[2]
        );

        // Now solve the EXACT SAME unchanged System via `fs-orbit` harmonic balance!
        let hb = sys
            .solve_periodic_hb(
                &map_tm,
                PeriodicHbConfig {
                    anchor: HbAnchor::Forced { omega: w },
                    budget: HbBudget {
                        harmonics: 5,
                        max_newton: 30,
                        tolerance: 1.0e-7,
                        deflate_equilibrium: false,
                    },
                    scales: HbScales {
                        states: vec![0.2, 1.0, 300.0],
                        equations: vec![1.0, 5.0, 1.0],
                        angular_frequency: w,
                    },
                    seed: None,
                    probe_dt_s: 1.0e-5,
                },
                cx,
            )
            .expect("fs-orbit harmonic balance on nonlinear thermo-mech Duffing system");

        // Compare HB waveform at t = 0 (`theta = 0`) against shooting initial_state `x0*`!
        let x_hb_0 = hb.orbit.orbit.sample(0, 0.0);
        let v_hb_0 = hb.orbit.orbit.sample(1, 0.0);
        let t_hb_0 = hb.orbit.orbit.sample(2, 0.0);
        assert!(
            (x_hb_0 - shoot.initial_state[0]).abs() < 2.0e-3,
            "x(0): HB {x_hb_0} vs shooting {}",
            shoot.initial_state[0]
        );
        assert!(
            (v_hb_0 - shoot.initial_state[1]).abs() < 5.0e-3,
            "v(0): HB {v_hb_0} vs shooting {}",
            shoot.initial_state[1]
        );
        assert!(
            (t_hb_0 - shoot.initial_state[2]).abs() < 5.0e-3,
            "T(0): HB {t_hb_0} vs shooting {}",
            shoot.initial_state[2]
        );
        // Confirm cubic Duffing nonlinearity actually excited a non-trivial 3rd harmonic!
        let amp3 = 2.0 * hb.orbit.orbit.coeffs[0][3].abs();
        assert!(
            amp3 > 1.0e-4,
            "expected nonzero 3rd harmonic from cubic spring, got {amp3:.3e}"
        );
    });
}

#[test]
#[allow(clippy::too_many_lines)]
fn adjoint_verified_via_verify_gradient_case1_thermo_mechanical_and_trajectory() {
    // Composed Case 1: Nonlinear cubic spring + mass + damper (quadratic heat)
    // + heat capacity + conductance with 2 active boundaries (force + T_amb).
    // Verifies:
    //  (a) single-step IFT adjoint via `fs_adjoint::verify_gradient`
    //  (b) tangent-adjoint dot-product identity `<J^{-1} db, dj> == <db, g>`
    //  (c) card-derived parameter gradient (stiffness k, damping c) via pullback + `verify_gradient`
    //  (d) multi-step checkpointed `trajectory_adjoint` via `verify_gradient`!
    let m = 1.4;
    let k = 32.0;
    let cubic = 45.0;
    let c = 2.2;
    let c_th = 80.0;
    let g_th = 6.0;
    let f0 = 7.5;
    let t_amb = 305.0;
    let dt = 0.02;

    let (mut sys, map, f_cmd, _c_b) = build_thermo_mech_oscillator(
        m,
        k,
        cubic,
        c,
        c_th,
        g_th,
        move |_, _| f0,
        t_amb,
        0.18,
        0.75,
        320.0,
    );

    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let base_cp = sys.checkpoint().unwrap();

        // Find indices for a nonlinear composite step objective J(u, y):
        // J = 1.3 * v_mass + 0.7 * F_spring_b - 0.4 * q_damper + 0.25 * T_body + 0.1 * u_f_cmd.
        let mass_p = sys.port(map.mass, "p").unwrap();
        let spring_b = sys.port(map.spring, "b").unwrap();
        let body_h = sys.port(map.body, "heat").unwrap();
        let idx_v = sys.output_index(mass_p, 0).unwrap();
        let idx_fs = sys.output_index(spring_b, 0).unwrap();
        let idx_tb = sys.output_index(body_h, 0).unwrap();
        let idx_ucmd = sys.input_index(f_cmd, 0).unwrap();

        let mut dj_dy = vec![0.0_f64; sys.output_len()];
        dj_dy[idx_v] = 1.3;
        dj_dy[idx_fs] = 0.7;
        dj_dy[idx_tb] = 0.25;
        let mut dj_du = vec![0.0_f64; sys.input_len()];
        dj_du[idx_ucmd] = 0.1;
        let obj = StepObjective {
            dj_du: dj_du.clone(),
            dj_dy: dj_dy.clone(),
            dj_db: Vec::new(),
        };

        // (a) Single-step IFT adjoint via `fs_adjoint::ift_gradient_matfree`.
        let adj = sys
            .step_sensitivity_adjoint(dt, &obj, SensitivityConfig::default(), cx)
            .expect("single-step IFT adjoint");
        assert!(adj.report.converged);
        // Checkpoint must be restored (commit_step = false).
        assert_eq!(sys.steps(), 0);

        // (b) Tangent-adjoint dot-product identity across all boundary scalars!
        let db_vec = vec![0.3, -0.7, 1.2, -0.5]; // boundary order: s_a, d_a, f_cmd, c_b
        let tan = sys
            .step_sensitivity_tangent(dt, &db_vec, SensitivityConfig::default(), cx)
            .expect("single-step forward tangent");
        let lhs: f64 = dj_du
            .iter()
            .zip(&tan.d_inputs)
            .map(|(a, b)| a * b)
            .sum::<f64>()
            + dj_dy
                .iter()
                .zip(&tan.d_outputs)
                .map(|(a, b)| a * b)
                .sum::<f64>();
        let rhs: f64 = adj
            .boundary_gradient
            .iter()
            .zip(&db_vec)
            .map(|(a, b)| a * b)
            .sum();
        assert!(
            (lhs - rhs).abs() < 1.0e-11 * lhs.abs().max(1.0),
            "tangent-adjoint mismatch: lhs={lhs:.12e} vs rhs={rhs:.12e}"
        );

        // Verify dJ / db against central FD using `fs_adjoint::verify_gradient`!
        let channels = sys.boundary_channels();
        assert_eq!(channels.len(), 4);
        let eval_j = |p: &[f64]| -> f64 {
            let (mut s_eval, map_e, f_cmd_e, _) = build_thermo_mech_oscillator(
                m,
                k,
                cubic,
                c,
                c_th,
                g_th,
                move |_, _| f0,
                t_amb,
                0.18,
                0.75,
                320.0,
            );
            let chs = s_eval.boundary_channels();
            for (ch, &bias) in chs.iter().zip(p) {
                s_eval
                    .set_boundary_bias(ch.port, ch.component, bias)
                    .unwrap();
            }
            with_cx(&CancelGate::new(), |cx2| {
                s_eval.step(dt, cx2).unwrap();
            });
            let mp = s_eval.port(map_e.mass, "p").unwrap();
            let sb = s_eval.port(map_e.spring, "b").unwrap();
            let bh = s_eval.port(map_e.body, "heat").unwrap();
            1.3 * s_eval.last_output(mp)[0]
                + 0.7 * s_eval.last_output(sb)[0]
                + 0.25 * s_eval.last_output(bh)[0]
                + 0.1 * s_eval.last_input(f_cmd_e)[0]
        };
        let p0 = vec![0.0_f64; 4];
        let dirs = vec![
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],
            vec![0.4, -0.6, 0.8, -0.3],
        ];
        let verdict = verify_gradient(&eval_j, &p0, &adj.boundary_gradient, &dirs, 1.0e-4, 1.0e-6);
        assert!(
            verdict.pass,
            "case1 single-step IFT gradient failed verify_gradient: {verdict:?}"
        );
        assert_eq!(verdict.informative_directions, dirs.len());

        // (c) Card-derived parameter gradient for [stiffness k, damping c] using
        //     `adjoint_outputs` returned by `step_sensitivity_adjoint_with_pullback`!
        let mut sys_p = build_thermo_mech_oscillator(
            m,
            k,
            cubic,
            c,
            c_th,
            g_th,
            move |_, _| f0,
            t_amb,
            0.18,
            0.75,
            320.0,
        )
        .0;
        let step_ref = {
            let cp = sys_p.checkpoint().unwrap();
            sys_p.step(dt, cx).unwrap();
            let s_a_in = sys_p.last_input(sys_p.port(map.spring, "a").unwrap())[0];
            let s_b_in = sys_p.last_input(sys_p.port(map.spring, "b").unwrap())[0];
            let d_a_in = sys_p.last_input(sys_p.port(map.damper, "a").unwrap())[0];
            let d_b_in = sys_p.last_input(sys_p.port(map.damper, "b").unwrap())[0];
            sys_p.restore(&cp).unwrap();
            (s_a_in, s_b_in, d_a_in, d_b_in)
        };
        let x_old = 0.18;
        let x_new = x_old + dt * (step_ref.0 - step_ref.1);
        let rate = step_ref.2 - step_ref.3;
        let spring_out_a = sys_p
            .output_index(sys_p.port(map.spring, "a").unwrap(), 0)
            .unwrap();
        let spring_out_b = sys_p
            .output_index(sys_p.port(map.spring, "b").unwrap(), 0)
            .unwrap();
        let damper_out_a = sys_p
            .output_index(sys_p.port(map.damper, "a").unwrap(), 0)
            .unwrap();
        let damper_out_b = sys_p
            .output_index(sys_p.port(map.damper, "b").unwrap(), 0)
            .unwrap();
        let damper_out_h = sys_p
            .output_index(sys_p.port(map.damper, "heat").unwrap(), 0)
            .unwrap();

        let param_pullback = |_lambda: &[f64], bar_y: &[f64]| -> Vec<f64> {
            // Sign convention: `dJ/dp = dj_dp - pullback`, so `pullback = -(dy/dp)^T bar_y`.
            // Spring discrete force: F = 0.5*k*(x_old + x_new) + ..., so dF/dk = 0.5*(x_old + x_new).
            let df_dk = 0.5 * (x_old + x_new);
            let dj_dk = bar_y[spring_out_a] * df_dk + bar_y[spring_out_b] * (-df_dk);
            // Damper: out_a = c*rate, out_b = -c*rate, out_h = -c*rate^2.
            let dj_dc = bar_y[damper_out_a] * rate
                + bar_y[damper_out_b] * (-rate)
                + bar_y[damper_out_h] * (-rate * rate);
            vec![-dj_dk, -dj_dc]
        };
        let adj_p = sys_p
            .step_sensitivity_adjoint_with_pullback(
                dt,
                &obj,
                &[],
                Some(&param_pullback),
                SensitivityConfig::default(),
                cx,
            )
            .expect("card-parameter pullback adjoint");
        let eval_j_params = |params: &[f64]| -> f64 {
            let (mut se, me, fce, _) = build_thermo_mech_oscillator(
                m,
                params[0],
                cubic,
                params[1],
                c_th,
                g_th,
                move |_, _| f0,
                t_amb,
                0.18,
                0.75,
                320.0,
            );
            with_cx(&CancelGate::new(), |cx2| {
                se.step(dt, cx2).unwrap();
            });
            let mp = se.port(me.mass, "p").unwrap();
            let sb = se.port(me.spring, "b").unwrap();
            let bh = se.port(me.body, "heat").unwrap();
            1.3 * se.last_output(mp)[0]
                + 0.7 * se.last_output(sb)[0]
                + 0.25 * se.last_output(bh)[0]
                + 0.1 * se.last_input(fce)[0]
        };
        let p_verdict = verify_gradient(
            &eval_j_params,
            &[k, c],
            &adj_p.parameter_gradient,
            &[vec![1.0, 0.0], vec![0.0, 1.0], vec![0.7, -0.5]],
            1.0e-5,
            1.0e-6,
        );
        assert!(
            p_verdict.pass,
            "card-parameter IFT gradient failed verify_gradient: {p_verdict:?}"
        );

        // (d) Multi-step checkpointed `trajectory_adjoint` over 6 steps!
        let n_steps = 6;
        let term_weights = vec![2.0, -1.5, 0.4];
        let mut step_parts = vec![TrajectoryStepPartials::default(); n_steps];
        step_parts[2].dj_dy = dj_dy.clone();
        step_parts[4].dj_ds = vec![0.5, 0.3, -0.1];

        let traj_adj = sys
            .trajectory_adjoint(
                &map,
                dt,
                n_steps,
                &term_weights,
                &step_parts,
                SensitivityConfig::default(),
                cx,
            )
            .expect("multi-step checkpointed trajectory_adjoint");
        assert_eq!(traj_adj.steps, n_steps);
        assert_eq!(sys.steps(), 0);
        // Verify initial checkpoint restored bit-for-bit.
        assert_eq!(
            sys.state::<f64>(map.spring).unwrap().to_bits(),
            base_cp.state::<f64>(map.spring).unwrap().to_bits()
        );

        // Joint vector `[x0, v0, T0, bias_f_cmd, bias_t_amb]` verified via `verify_gradient`!
        let combined_grad = vec![
            traj_adj.grad_initial_state[0],
            traj_adj.grad_initial_state[1],
            traj_adj.grad_initial_state[2],
            traj_adj.grad_boundary_total[2], // f_cmd
            traj_adj.grad_boundary_total[3], // c_b (T_amb)
        ];
        let eval_traj = |q: &[f64]| -> f64 {
            let (mut se, me, fce, cbe) = build_thermo_mech_oscillator(
                m,
                k,
                cubic,
                c,
                c_th,
                g_th,
                move |_, _| f0,
                t_amb,
                q[0],
                q[1],
                q[2],
            );
            se.set_boundary_bias(fce, 0, q[3]).unwrap();
            se.set_boundary_bias(cbe, 0, q[4]).unwrap();
            let mut total = 0.0_f64;
            let mut st = [0.0_f64; 3];
            with_cx(&CancelGate::new(), |cx2| {
                for step_k in 0..n_steps {
                    if step_k == 4 {
                        me.read(&se, &mut st).unwrap();
                        total += 0.5 * st[0] + 0.3 * st[1] - 0.1 * st[2];
                    }
                    se.step(dt, cx2).unwrap();
                    if step_k == 2 {
                        let mp = se.port(me.mass, "p").unwrap();
                        let sb = se.port(me.spring, "b").unwrap();
                        let bh = se.port(me.body, "heat").unwrap();
                        total += 1.3 * se.last_output(mp)[0]
                            + 0.7 * se.last_output(sb)[0]
                            + 0.25 * se.last_output(bh)[0];
                    }
                }
                me.read(&se, &mut st).unwrap();
                total +=
                    term_weights[0] * st[0] + term_weights[1] * st[1] + term_weights[2] * st[2];
            });
            total
        };
        let q0 = vec![0.18, 0.75, 320.0, 0.0, 0.0];
        let q_dirs = vec![
            vec![1.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0, 1.0],
            vec![0.3, -0.5, 0.8, 0.4, -0.6],
        ];
        let t_verdict = verify_gradient(&eval_traj, &q0, &combined_grad, &q_dirs, 1.0e-5, 1.0e-5);
        assert!(
            t_verdict.pass,
            "case1 trajectory_adjoint failed verify_gradient: {t_verdict:?}"
        );
        assert_eq!(t_verdict.informative_directions, q_dirs.len());
    });
}

#[test]
#[allow(clippy::too_many_lines)]
fn adjoint_verified_via_verify_gradient_case2_equivalent_temperature_heat_route() {
    // Composed Case 2: Two-body thermal-mechanical coupled system with nonlinear
    // `RouteTemperature::Equivalent` heat routing (55% -> body_a at 310 K,
    // 35% -> body_b at 365 K, 10% -> ambient sink at 295 K), verifying exact
    // `Driver::EquivalentTemperature` transpose scatter + `Damper::vjp` +
    // `fs_adjoint::ift_gradient_matfree` against `verify_gradient`.
    let build_clean = || {
        let mut b = SystemBuilder::new();
        let mass = b
            .add_with_tangent(
                "mass",
                Mass::new(1.8).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s: 1.4,
                },
            )
            .unwrap();
        let damper = b
            .add_with_tangent("damper", Damper::new(3.5).unwrap(), ())
            .unwrap();
        let f_src = b
            .add_with_tangent(
                "drive",
                ControlledSource::effort(PortKind::MechanicalForceVelocity),
                SourceAccount::default(),
            )
            .unwrap();
        let body_a = b
            .add_with_tangent("body_a", HeatCapacity::new(45.0).unwrap(), 310.0)
            .unwrap();
        let body_b = b
            .add_with_tangent("body_b", HeatCapacity::new(65.0).unwrap(), 365.0)
            .unwrap();

        let m_p = b.port(mass, "p").unwrap();
        let d_a = b.port(damper, "a").unwrap();
        let d_b = b.port(damper, "b").unwrap();
        let d_h = b.port(damper, "heat").unwrap();
        let f_p = b.port(f_src, "p").unwrap();
        let f_cmd = b.port(f_src, "command").unwrap();
        let ba_h = b.port(body_a, "heat").unwrap();
        let bb_h = b.port(body_b, "heat").unwrap();

        b.common_flow(m_p, &[d_b, f_p]).unwrap();
        b.boundary(d_a, Box::new(|_, _, out| out[0] = 0.25))
            .unwrap();
        b.boundary(f_cmd, Box::new(|_, _, out| out[0] = 9.0))
            .unwrap();
        b.heat_route(
            d_h,
            &[(ba_h, 0.55), (bb_h, 0.35)],
            Some(fs_system::heat_sink("ambient", 295.0).unwrap()),
            RouteTemperature::SinkEquivalent,
        )
        .unwrap();

        let mut sys = b.seal().unwrap();
        sys.guess(d_h, &[325.0]);
        (sys, d_h, ba_h, bb_h, m_p)
    };

    let dt = 0.05;
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let (mut sys, d_h, ba_h, bb_h, m_p) = build_clean();
        // Objective explicitly includes `d_h` INPUT (the nonlinear
        // `EquivalentTemperature` driver value!) + `body_a`/`body_b` output
        // temperatures + `mass` velocity!
        let idx_t_eq = sys.input_index(d_h, 0).unwrap();
        let idx_ta = sys.output_index(ba_h, 0).unwrap();
        let idx_tb = sys.output_index(bb_h, 0).unwrap();
        let idx_v = sys.output_index(m_p, 0).unwrap();

        let mut dj_du = vec![0.0_f64; sys.input_len()];
        dj_du[idx_t_eq] = 2.5;
        let mut dj_dy = vec![0.0_f64; sys.output_len()];
        dj_dy[idx_ta] = 0.8;
        dj_dy[idx_tb] = -0.6;
        dj_dy[idx_v] = 1.1;

        let adj = sys
            .step_sensitivity_adjoint(
                dt,
                &StepObjective {
                    dj_du,
                    dj_dy,
                    dj_db: Vec::new(),
                },
                SensitivityConfig::default(),
                cx,
            )
            .expect("case2 EquivalentTemperature IFT adjoint");
        assert!(adj.report.converged);

        let eval_j = |p: &[f64]| -> f64 {
            let (mut se, dh, bah, bbh, mp) = build_clean();
            let chs = se.boundary_channels();
            for (ch, &bias) in chs.iter().zip(p) {
                se.set_boundary_bias(ch.port, ch.component, bias).unwrap();
            }
            with_cx(&CancelGate::new(), |cx2| {
                se.step(dt, cx2).unwrap();
            });
            2.5 * se.last_input(dh)[0] + 0.8 * se.last_output(bah)[0] - 0.6 * se.last_output(bbh)[0]
                + 1.1 * se.last_output(mp)[0]
        };

        let dirs = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![0.6, -0.8],
            vec![-0.4, 0.9],
        ];
        let verdict = verify_gradient(
            &eval_j,
            &[0.0, 0.0],
            &adj.boundary_gradient,
            &dirs,
            1.0e-4,
            1.0e-6,
        );
        assert!(
            verdict.pass,
            "case2 EquivalentTemperature IFT gradient failed verify_gradient: {verdict:?}"
        );
        assert_eq!(verdict.informative_directions, dirs.len());

        // Regression check: non-finite parameter pullback must fail closed as Singular.
        let nan_pullback = |_lambda: &[f64], _bar_y: &[f64]| vec![f64::NAN];
        let nan_err = sys
            .step_sensitivity_adjoint_with_pullback(
                dt,
                &StepObjective {
                    dj_du: Vec::new(),
                    dj_dy: Vec::new(),
                    dj_db: Vec::new(),
                },
                &[],
                Some(&nan_pullback),
                SensitivityConfig::default(),
                cx,
            )
            .expect_err("NaN parameter_pullback must refuse");
        assert!(matches!(
            nan_err,
            AnalysisError::Singular {
                analysis: AnalysisKind::AdjointSensitivity,
                ..
            }
        ));

        // Regression check: TangentSensitivity works through unified `run_analysis` switch.
        let tan_out = sys
            .run_analysis(
                AnalysisRequest::TangentSensitivity {
                    dt_s: dt,
                    d_boundaries: &[0.5, -0.25],
                    config: SensitivityConfig::default(),
                },
                cx,
            )
            .expect("TangentSensitivity via run_analysis");
        match tan_out {
            AnalysisOutcome::TangentSensitivity(t) => {
                assert!(t.d_outputs.iter().all(|v| v.is_finite()));
            }
            other => panic!("expected TangentSensitivity outcome, got {other:?}"),
        }
    });
}
