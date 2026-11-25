//! Direct verification of `StageResidual` on default `fs-system` lumped
//! elements (`Mass`, `Spring`, `Damper`, `Force`, `HeatCapacity`, `Conductance`,
//! `Heater`), derivative identities (`jvp_z` FD + `vjp_z` adjoint), hand-derived
//! steady-state oracles (`x* = F/k`, `T* = T_amb + Q/G`), monolithic parity,
//! and explicit `MissingTier(StageResidual)` refusal on elements without
//! `StageResidual` (`Storage`).

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::{Component, Proposal, StageResidual};
use fs_system::elements::{
    Conductance, Damper, DissipationLaw, ElementError, Force, HeatCapacity, Heater, Mass, Motion,
    SourceAccount, Spring, Storage, Waveform,
};
use fs_system::{
    AnalysisError, AnalysisKind, BlockSolveMode, BlockSolverPolicy, NodeId, RequiredTier,
    RouteTemperature, StateMap, StationaryConfig, System, SystemBuilder,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 2026,
                kernel_id: 2117,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn assert_close(actual: f64, expected: f64, tol: f64, label: &str) {
    assert!(
        (actual - expected).abs() <= tol,
        "{label}: actual={actual:.16e}, expected={expected:.16e}, diff={:.3e} > tol={tol:.3e}",
        (actual - expected).abs()
    );
}

#[derive(Debug)]
struct ThermalFlowLaw;
impl DissipationLaw for ThermalFlowLaw {
    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(0x5448_4552_4d31);
    }
    fn factor(&self, temperature_k: f64, flow: f64) -> Result<(f64, [f64; 2]), String> {
        let f_t = 1.0 + 1.0e-3 * (temperature_k - 300.0);
        let f_u = 1.0 + 0.05 * flow * flow;
        Ok((f_t * f_u, [1.0e-3 * f_u, f_t * 0.1 * flow]))
    }
}

#[allow(clippy::too_many_arguments)]
fn check_jvp_vjp_identities<C: StageResidual>(
    label: &str,
    comp: &C,
    old: &C::State,
    dt_s: f64,
    inputs: &[f64],
    z: &[f64],
    dz: &[f64],
    bar_r: &[f64],
    bar_outputs: &[f64],
    cx: &Cx<'_>,
) {
    let nz = comp.unknowns();
    let n_out = comp
        .ports()
        .iter()
        .map(fs_ports::component::PortSpec::width)
        .sum::<usize>();
    let mut dr = vec![f64::NAN; nz];
    let mut d_out = vec![f64::NAN; n_out];
    comp.jvp_z(old, dt_s, inputs, z, dz, &mut dr, &mut d_out, cx)
        .unwrap_or_else(|e| panic!("{label} jvp_z failed: {e:?}"));

    if nz == 0 {
        assert!(
            d_out.iter().all(|&v| v.abs() <= 0.0),
            "{label}: 0-unknown jvp_z must zero d_outputs, got {d_out:?}"
        );
        let mut bar_z = Vec::new();
        comp.vjp_z(old, dt_s, inputs, z, bar_r, bar_outputs, &mut bar_z, cx)
            .unwrap_or_else(|e| panic!("{label} vjp_z failed: {e:?}"));
        return;
    }

    // Central finite differences on (r(z), outputs(z)).
    let h = 1.0e-6;
    let zp: Vec<f64> = z.iter().zip(dz).map(|(a, b)| a + h * b).collect();
    let zm: Vec<f64> = z.iter().zip(dz).map(|(a, b)| a - h * b).collect();
    let (mut rp, mut rm) = (vec![0.0; nz], vec![0.0; nz]);
    let (mut op, mut om) = (vec![0.0; n_out], vec![0.0; n_out]);
    comp.residual(old, dt_s, inputs, &zp, &mut rp, &mut op, cx)
        .unwrap();
    comp.residual(old, dt_s, inputs, &zm, &mut rm, &mut om, cx)
        .unwrap();
    for i in 0..nz {
        let fd = (rp[i] - rm[i]) / (2.0 * h);
        let tol = 1.0e-7 * fd.abs().max(1.0);
        assert_close(dr[i], fd, tol, &format!("{label} dr[{i}] vs FD"));
    }
    for i in 0..n_out {
        let fd = (op[i] - om[i]) / (2.0 * h);
        let tol = 1.0e-7 * fd.abs().max(1.0);
        assert_close(d_out[i], fd, tol, &format!("{label} d_out[{i}] vs FD"));
    }

    // Adjoint dot-product identity: <dr, bar_r> + <d_out, bar_out> == <dz, bar_z>.
    let mut bar_z = vec![f64::NAN; nz];
    comp.vjp_z(old, dt_s, inputs, z, bar_r, bar_outputs, &mut bar_z, cx)
        .unwrap_or_else(|e| panic!("{label} vjp_z failed: {e:?}"));
    let lhs: f64 = dr
        .iter()
        .zip(bar_r)
        .map(|(a, b)| a * b)
        .chain(d_out.iter().zip(bar_outputs).map(|(a, b)| a * b))
        .sum();
    let rhs: f64 = dz.iter().zip(&bar_z).map(|(a, b)| a * b).sum();
    let tol = 1.0e-12 * lhs.abs().max(rhs.abs()).max(1.0);
    assert_close(lhs, rhs, tol, &format!("{label} adjoint dot identity"));
}

#[test]
#[allow(clippy::too_many_lines)]
fn residual_and_accept_match_propose_to_roundoff_for_all_seven_elements() {
    with_cx(|cx| {
        let dt = 2.5e-3;

        // 1. Mass (plain + publishing)
        for (tag, mass) in [
            ("Mass", Mass::new(2.5).unwrap()),
            ("Mass::publishing", Mass::new(2.5).unwrap().publishing()),
        ] {
            let old = Motion {
                position_m: 0.4,
                velocity_m_s: -1.3,
            };
            let inputs = [14.7];
            let n_out = mass.ports().len();
            let mut new_p = old;
            let mut out_p = vec![0.0; n_out];
            let Proposal::Step(rec_p) = mass
                .propose(&old, dt, &inputs, &mut new_p, &mut out_p, cx)
                .unwrap()
            else {
                panic!("{tag} must step");
            };
            let mut z = [0.0];
            mass.initial_guess(&old, dt, &inputs, &mut z).unwrap();
            let mut r = [f64::NAN];
            let mut out_r = vec![0.0; n_out];
            mass.residual(&old, dt, &inputs, &z, &mut r, &mut out_r, cx)
                .unwrap();
            assert_close(r[0], 0.0, 1.0e-13, &format!("{tag} residual"));
            for i in 0..n_out {
                assert_close(out_r[i], out_p[i], 1.0e-14, &format!("{tag} out[{i}]"));
            }
            let mut new_a = old;
            let rec_a = mass.accept(&old, dt, &inputs, &z, &mut new_a).unwrap();
            assert_close(new_a.position_m, new_p.position_m, 1.0e-14, "Mass pos");
            assert_close(new_a.velocity_m_s, new_p.velocity_m_s, 1.0e-14, "Mass vel");
            assert_close(
                rec_a.stored_change_j,
                rec_p.stored_change_j,
                1.0e-13,
                "Mass dE",
            );
        }

        // 2. Spring (linear + cubic)
        for (tag, spring) in [
            ("Spring(linear)", Spring::new(150.0, 0.0).unwrap()),
            ("Spring(cubic)", Spring::new(150.0, 42.0).unwrap()),
        ] {
            let old = 0.18;
            let inputs = [1.2, -0.7];
            let (mut new_p, mut out_p) = (0.0, [0.0; 2]);
            let Proposal::Step(rec_p) = spring
                .propose(&old, dt, &inputs, &mut new_p, &mut out_p, cx)
                .unwrap()
            else {
                panic!("{tag} must step");
            };
            let mut z = [0.0];
            spring.initial_guess(&old, dt, &inputs, &mut z).unwrap();
            let (mut r, mut out_r) = ([f64::NAN], [0.0; 2]);
            spring
                .residual(&old, dt, &inputs, &z, &mut r, &mut out_r, cx)
                .unwrap();
            assert_close(r[0], 0.0, 1.0e-14, &format!("{tag} residual"));
            for i in 0..2 {
                assert_close(out_r[i], out_p[i], 1.0e-13, &format!("{tag} out[{i}]"));
            }
            let mut new_a = 0.0;
            let rec_a = spring.accept(&old, dt, &inputs, &z, &mut new_a).unwrap();
            assert_close(new_a, new_p, 1.0e-14, "Spring ext");
            assert_close(
                rec_a.stored_change_j,
                rec_p.stored_change_j,
                1.0e-13,
                "Spring dE",
            );
        }

        // 3. HeatCapacity (plain + publishing + nonpositive refusal)
        for (tag, hc) in [
            ("HeatCapacity", HeatCapacity::new(420.0).unwrap()),
            (
                "HeatCapacity::publishing",
                HeatCapacity::new(420.0).unwrap().publishing(),
            ),
        ] {
            let old = 310.0;
            let inputs = [85.0];
            let n_out = hc.ports().len();
            let mut new_p = 0.0;
            let mut out_p = vec![0.0; n_out];
            let Proposal::Step(rec_p) = hc
                .propose(&old, dt, &inputs, &mut new_p, &mut out_p, cx)
                .unwrap()
            else {
                panic!("{tag} must step");
            };
            let mut z = [0.0];
            hc.initial_guess(&old, dt, &inputs, &mut z).unwrap();
            let mut r = [f64::NAN];
            let mut out_r = vec![0.0; n_out];
            hc.residual(&old, dt, &inputs, &z, &mut r, &mut out_r, cx)
                .unwrap();
            assert_close(r[0], 0.0, 1.0e-11, &format!("{tag} residual"));
            for i in 0..n_out {
                assert_close(out_r[i], out_p[i], 1.0e-13, &format!("{tag} out[{i}]"));
            }
            let mut new_a = 0.0;
            let rec_a = hc.accept(&old, dt, &inputs, &z, &mut new_a).unwrap();
            assert_close(new_a, new_p, 1.0e-13, "HeatCapacity T");
            assert_close(
                rec_a.stored_change_j,
                rec_p.stored_change_j,
                1.0e-11,
                "HeatCapacity dE",
            );
            let bad = hc.residual(&old, dt, &inputs, &[-1.0], &mut r, &mut out_r, cx);
            assert!(matches!(bad, Err(ElementError::Temperature { .. })));
        }

        // 4. Damper (plain + DissipationLaw)
        for (tag, damper) in [
            ("Damper(plain)", Damper::new(8.5).unwrap()),
            (
                "Damper(law)",
                Damper::new(8.5)
                    .unwrap()
                    .following(std::sync::Arc::new(ThermalFlowLaw)),
            ),
        ] {
            assert_eq!(damper.unknowns(), 0);
            let inputs = [1.4, -0.3, 325.0];
            let (mut new_p, mut out_p) = ((), [0.0; 3]);
            let Proposal::Step(rec_p) = damper
                .propose(&(), dt, &inputs, &mut new_p, &mut out_p, cx)
                .unwrap()
            else {
                panic!("{tag} must step");
            };
            let mut out_r = [0.0; 3];
            damper
                .residual(&(), dt, &inputs, &[], &mut [], &mut out_r, cx)
                .unwrap();
            for i in 0..3 {
                assert_close(out_r[i], out_p[i], 1.0e-14, &format!("{tag} out[{i}]"));
            }
            let rec_a = damper.accept(&(), dt, &inputs, &[], &mut ()).unwrap();
            assert_close(rec_a.dissipated_j, rec_p.dissipated_j, 1.0e-14, tag);
            assert_close(
                rec_a.entropy_generated_j_per_k,
                rec_p.entropy_generated_j_per_k,
                1.0e-14,
                tag,
            );
        }

        // 5. Force
        let force = Force::new(Waveform::sine(12.0, 3.0, 0.2, 5.0).unwrap());
        assert_eq!(force.unknowns(), 0);
        let old_f = SourceAccount::default();
        let (mut new_fp, mut out_fp) = (old_f, [0.0]);
        let Proposal::Step(rec_fp) = force
            .propose(&old_f, dt, &[0.9], &mut new_fp, &mut out_fp, cx)
            .unwrap()
        else {
            panic!("Force must step");
        };
        let mut out_fr = [0.0];
        force
            .residual(&old_f, dt, &[0.9], &[], &mut [], &mut out_fr, cx)
            .unwrap();
        assert_close(out_fr[0], out_fp[0], 1.0e-14, "Force out");
        let mut new_fa = old_f;
        let rec_fa = force.accept(&old_f, dt, &[0.9], &[], &mut new_fa).unwrap();
        assert_close(
            rec_fa.stored_change_j,
            rec_fp.stored_change_j,
            1.0e-14,
            "Force dE",
        );
        assert_close(
            new_fa.energy_in_j,
            new_fp.energy_in_j,
            1.0e-14,
            "Force account",
        );

        // 6. Conductance
        let cond = Conductance::new(6.4).unwrap();
        assert_eq!(cond.unknowns(), 0);
        let (mut new_cp, mut out_cp) = ((), [0.0; 2]);
        let Proposal::Step(rec_cp) = cond
            .propose(&(), dt, &[350.0, 295.0], &mut new_cp, &mut out_cp, cx)
            .unwrap()
        else {
            panic!("Conductance must step");
        };
        let mut out_cr = [0.0; 2];
        cond.residual(&(), dt, &[350.0, 295.0], &[], &mut [], &mut out_cr, cx)
            .unwrap();
        for i in 0..2 {
            assert_close(out_cr[i], out_cp[i], 1.0e-14, "Conductance out");
        }
        let rec_ca = cond.accept(&(), dt, &[350.0, 295.0], &[], &mut ()).unwrap();
        assert_close(
            rec_ca.entropy_generated_j_per_k,
            rec_cp.entropy_generated_j_per_k,
            1.0e-14,
            "Conductance dS",
        );

        // 7. Heater
        let heater = Heater::new(Waveform::constant(45.0).unwrap()).unwrap();
        assert_eq!(heater.unknowns(), 0);
        let old_h = SourceAccount::default();
        let (mut new_hp, mut out_hp) = (old_h, [0.0]);
        let Proposal::Step(rec_hp) = heater
            .propose(&old_h, dt, &[312.0], &mut new_hp, &mut out_hp, cx)
            .unwrap()
        else {
            panic!("Heater must step");
        };
        let mut out_hr = [0.0];
        heater
            .residual(&old_h, dt, &[312.0], &[], &mut [], &mut out_hr, cx)
            .unwrap();
        assert_close(out_hr[0], out_hp[0], 1.0e-14, "Heater out");
        let mut new_ha = old_h;
        let rec_ha = heater
            .accept(&old_h, dt, &[312.0], &[], &mut new_ha)
            .unwrap();
        assert_close(
            rec_ha.dissipated_j,
            rec_hp.dissipated_j,
            1.0e-14,
            "Heater diss",
        );
        assert_close(
            rec_ha.entropy_generated_j_per_k,
            rec_hp.entropy_generated_j_per_k,
            1.0e-14,
            "Heater dS",
        );
    });
}

#[test]
#[allow(clippy::too_many_lines)]
fn jvp_z_matches_finite_differences_and_vjp_z_satisfies_adjoint_identity() {
    with_cx(|cx| {
        let dt = 3.0e-3;
        check_jvp_vjp_identities(
            "Mass",
            &Mass::new(3.2).unwrap(),
            &Motion {
                position_m: 0.25,
                velocity_m_s: -0.9,
            },
            dt,
            &[11.0],
            &[-0.82],
            &[0.73],
            &[1.4],
            &[-2.1],
            cx,
        );
        check_jvp_vjp_identities(
            "Mass::publishing",
            &Mass::new(3.2).unwrap().publishing(),
            &Motion {
                position_m: 0.25,
                velocity_m_s: -0.9,
            },
            dt,
            &[11.0],
            &[-0.82],
            &[0.73],
            &[1.4],
            &[-2.1, 3.3, -4.4],
            cx,
        );
        check_jvp_vjp_identities(
            "Spring(cubic)",
            &Spring::new(95.0, 38.0).unwrap(),
            &0.21,
            dt,
            &[0.8, -0.5],
            &[0.34],
            &[-0.65],
            &[1.7],
            &[2.3, -1.9],
            cx,
        );
        check_jvp_vjp_identities(
            "HeatCapacity::publishing",
            &HeatCapacity::new(620.0).unwrap().publishing(),
            &305.0,
            dt,
            &[60.0],
            &[309.5],
            &[1.25],
            &[-0.8],
            &[2.7, 5.1],
            cx,
        );
        check_jvp_vjp_identities(
            "Damper",
            &Damper::new(5.0).unwrap(),
            &(),
            dt,
            &[1.0, -0.5, 300.0],
            &[],
            &[],
            &[],
            &[1.0, -1.0, 0.5],
            cx,
        );
        check_jvp_vjp_identities(
            "Force",
            &Force::new(Waveform::constant(7.0).unwrap()),
            &SourceAccount::default(),
            dt,
            &[0.5],
            &[],
            &[],
            &[],
            &[1.0],
            cx,
        );
        check_jvp_vjp_identities(
            "Conductance",
            &Conductance::new(4.0).unwrap(),
            &(),
            dt,
            &[330.0, 295.0],
            &[],
            &[],
            &[],
            &[1.0, -1.0],
            cx,
        );
        check_jvp_vjp_identities(
            "Heater",
            &Heater::new(Waveform::constant(20.0).unwrap()).unwrap(),
            &SourceAccount::default(),
            dt,
            &[300.0],
            &[],
            &[],
            &[],
            &[1.0],
            cx,
        );
    });
}

struct SpringMassMap {
    spring: NodeId,
    mass: NodeId,
}
impl StateMap for SpringMassMap {
    fn dim(&self) -> usize {
        2
    }
    fn read(&self, system: &System, out: &mut [f64]) -> Result<(), AnalysisError> {
        out[0] = *system
            .state::<f64>(self.spring)
            .ok_or(AnalysisError::Invalid { what: "spring" })?;
        out[1] = system
            .state::<Motion>(self.mass)
            .ok_or(AnalysisError::Invalid { what: "mass" })?
            .velocity_m_s;
        Ok(())
    }
    fn write(&self, system: &mut System, state: &[f64]) -> Result<(), AnalysisError> {
        system.set_accepted_state(self.spring, state[0])?;
        let mut m = *system
            .state::<Motion>(self.mass)
            .ok_or(AnalysisError::Invalid { what: "mass" })?;
        m.velocity_m_s = state[1];
        system.set_accepted_state(self.mass, m)?;
        Ok(())
    }
}

struct ThermalScalarMap {
    body: NodeId,
}
impl StateMap for ThermalScalarMap {
    fn dim(&self) -> usize {
        1
    }
    fn read(&self, system: &System, out: &mut [f64]) -> Result<(), AnalysisError> {
        out[0] = *system
            .state::<f64>(self.body)
            .ok_or(AnalysisError::Invalid { what: "body" })?;
        Ok(())
    }
    fn write(&self, system: &mut System, state: &[f64]) -> Result<(), AnalysisError> {
        if !(state[0].is_finite() && state[0] > 0.0) {
            return Err(AnalysisError::Invalid {
                what: "non-positive temperature",
            });
        }
        system.set_accepted_state(self.body, state[0])?;
        Ok(())
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn steady_state_oracles_forced_spring_mass_and_heater_conductance_sink() {
    with_cx(|cx| {
        // Oracle 1: forced spring-mass-damper (x* = F/k = 10 / 200 = 0.05 m, v* = 0).
        let (k, force_n) = (200.0, 10.0);
        let mut bm = SystemBuilder::new();
        let mass = bm
            .add_with_tangent_and_stage_residual(
                "mass",
                Mass::new(2.0).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s: -0.7,
                },
            )
            .unwrap();
        let spring = bm
            .add_with_tangent_and_stage_residual("spring", Spring::new(k, 0.0).unwrap(), 0.2)
            .unwrap();
        let damper = bm
            .add_with_tangent_and_stage_residual("damper", Damper::new(4.0).unwrap(), ())
            .unwrap();
        let drive = bm
            .add_with_tangent_and_stage_residual(
                "drive",
                Force::new(Waveform::constant(force_n).unwrap()),
                SourceAccount::default(),
            )
            .unwrap();
        let (mp, sa, sb, da, db, dh, fp) = (
            bm.port(mass, "p").unwrap(),
            bm.port(spring, "a").unwrap(),
            bm.port(spring, "b").unwrap(),
            bm.port(damper, "a").unwrap(),
            bm.port(damper, "b").unwrap(),
            bm.port(damper, "heat").unwrap(),
            bm.port(drive, "p").unwrap(),
        );
        // Port `b` on mass, `a` anchored at 0 => x_dot = +v, F_on_mass = -k*x - c*v + F.
        bm.common_flow(mp, &[sb, db, fp]).unwrap();
        bm.boundary(sa, Box::new(|_, _, out| out[0] = 0.0)).unwrap();
        bm.boundary(da, Box::new(|_, _, out| out[0] = 0.0)).unwrap();
        bm.heat_route(
            dh,
            &[],
            Some(fs_system::heat_sink("amb", 300.0).unwrap()),
            RouteTemperature::Fixed(300.0),
        )
        .unwrap();
        let mut sys_m = bm.seal().unwrap();
        let sol_m = sys_m
            .solve_stationary(
                &SpringMassMap { spring, mass },
                StationaryConfig::default(),
                cx,
            )
            .expect("forced spring-mass stationary converges");
        assert!(sol_m.residual_norm < 1.0e-10, "res={}", sol_m.residual_norm);
        assert_close(sol_m.state[0], force_n / k, 1.0e-10, "x* = F/k");
        assert_close(sol_m.state[1], 0.0, 1.0e-10, "v* = 0");

        // Oracle 1b: purely stateful undamped cubic Mass+Spring under MonolithicNewtonKrylov
        // agrees with InterfaceNewton to roundoff.
        let build_undamped = || {
            let mut b = SystemBuilder::new();
            let m = b
                .add_with_tangent_and_stage_residual(
                    "m",
                    Mass::new(1.5).unwrap(),
                    Motion {
                        position_m: 0.0,
                        velocity_m_s: 0.8,
                    },
                )
                .unwrap();
            let s = b
                .add_with_tangent_and_stage_residual("s", Spring::new(80.0, 25.0).unwrap(), 0.25)
                .unwrap();
            let (mp, sa, sb) = (
                b.port(m, "p").unwrap(),
                b.port(s, "a").unwrap(),
                b.port(s, "b").unwrap(),
            );
            b.common_flow(mp, &[sb]).unwrap();
            b.boundary(sa, Box::new(|_, _, out| out[0] = 0.0)).unwrap();
            (b.seal().unwrap(), m, s)
        };
        let (mut sys_mono, m_mono, s_mono) = build_undamped();
        let (mut sys_iface, m_iface, s_iface) = build_undamped();
        assert_eq!(
            sys_mono.execution_plan_for_dt(1.0e-3).blocks[0].solve_mode,
            BlockSolveMode::MonolithicNewtonKrylov
        );
        sys_iface.set_solver_policy(BlockSolverPolicy::InterfaceNewton);
        for _ in 0..20 {
            sys_mono.step(1.0e-3, cx).unwrap();
            sys_iface.step(1.0e-3, cx).unwrap();
        }
        assert_close(
            sys_mono.state::<Motion>(m_mono).unwrap().velocity_m_s,
            sys_iface.state::<Motion>(m_iface).unwrap().velocity_m_s,
            1.0e-12,
            "monolithic vs interface-newton velocity",
        );
        assert_close(
            *sys_mono.state::<f64>(s_mono).unwrap(),
            *sys_iface.state::<f64>(s_iface).unwrap(),
            1.0e-12,
            "monolithic vs interface-newton extension",
        );

        // Oracle 2: heater-conductance-sink thermal chain (T* = T_amb + Q / G = 315 K).
        let (t_amb, q_w, g_w_k) = (300.0, 75.0, 5.0);
        let mut bt = SystemBuilder::new();
        let body = bt
            .add_with_tangent_and_stage_residual("body", HeatCapacity::new(500.0).unwrap(), 290.0)
            .unwrap();
        let cond = bt
            .add_with_tangent_and_stage_residual("cond", Conductance::new(g_w_k).unwrap(), ())
            .unwrap();
        let heater = bt
            .add_with_tangent_and_stage_residual(
                "heater",
                Heater::new(Waveform::constant(q_w).unwrap()).unwrap(),
                SourceAccount::default(),
            )
            .unwrap();
        let (bh, ca, cb, hh) = (
            bt.port(body, "heat").unwrap(),
            bt.port(cond, "a").unwrap(),
            bt.port(cond, "b").unwrap(),
            bt.port(heater, "heat").unwrap(),
        );
        bt.common_effort(bh, &[ca, hh]).unwrap();
        bt.boundary(cb, Box::new(move |_, _, out| out[0] = t_amb))
            .unwrap();
        let mut sys_t = bt.seal().unwrap();
        assert!(sys_t.guess(ca, &[290.0]));
        assert!(sys_t.guess(hh, &[290.0]));
        let sol_t = sys_t
            .solve_stationary(&ThermalScalarMap { body }, StationaryConfig::default(), cx)
            .expect("heater-conductance-sink stationary converges");
        assert!(sol_t.residual_norm < 1.0e-10, "res={}", sol_t.residual_norm);
        // |dT/dt| <= 1e-10 K/s with tau = C/G = 100 s implies |T* - T| <= 1e-8 K (< 1e-10 rel).
        assert_close(
            sol_t.state[0],
            t_amb + q_w / g_w_k,
            1.0e-8,
            "T* = T_amb + Q/G",
        );
    });
}

#[test]
fn missing_tier_stage_residual_refuses_storage_element() {
    with_cx(|cx| {
        // `Storage` intentionally does not implement `StageResidual` after S2-B.
        let mut b = SystemBuilder::new();
        let m = b
            .add_with_tangent_and_stage_residual(
                "mass",
                Mass::new(1.0).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s: 0.0,
                },
            )
            .unwrap();
        let s = b
            .add_with_tangent("cap_storage", Storage::new(0.05).unwrap(), 0.0)
            .unwrap();
        let (mp, sp) = (b.port(m, "p").unwrap(), b.port(s, "p").unwrap());
        b.common_flow(mp, &[sp]).unwrap();
        let mut sys = b.seal().unwrap();
        let err = sys
            .solve_stationary(
                &SpringMassMap { spring: s, mass: m },
                StationaryConfig::default(),
                cx,
            )
            .expect_err("Storage must trigger MissingTier(StageResidual)");
        assert!(matches!(
            err,
            AnalysisError::MissingTier {
                analysis: AnalysisKind::Stationary,
                ref node,
                required: RequiredTier::StageResidual,
            } if node == "cap_storage"
        ));
    });
}
