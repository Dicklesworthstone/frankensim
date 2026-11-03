//! G1/G3/G4 actual original MIMO sources on full hereditary solid equations.
use fs_couple::solid_rational::{
    PreparedSolidMimoResonance, SolidMimoResonanceBudget, SolidMimoResonanceError,
};
use fs_exec::Cx;
use fs_la::nonlinear_eig::PoleOptions;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::visco::*;
use fs_vfit::{
    mimo::{MimoModel, load::{MimoLoadLimits, PreparedMimoLoad}},
    model::PoleTerm,
};

fn with_cx(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 419,
                    kernel_id: 96,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}
fn diagonal(value: f64) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { value } else { 0.0 }))
}
// The one free displacement has original M=1, K=2.5. With memory it also
// has Kelvin C=0.2 and six original Maxwell strains with k=6.5, tau=0.25.
fn source(memory: bool, cx: &Cx<'_>) -> ViscoElements {
    let branches = if memory {
        vec![PronyBranch::new(diagonal(39.0), 0.25).unwrap()]
    } else {
        vec![]
    };
    let mut law = TensorProny::new_mandel(diagonal(15.0), branches, 0.01).unwrap();
    if memory {
        law = law.with_kelvin_viscosity_mandel(diagonal(1.2)).unwrap();
    }
    let material = TetViscoMaterial {
        law,
        density_kg_m3: 60.0,
        material_state_identity: [96; 32],
    };
    let fixed: Vec<_> = (0..12).filter(|&i| i != 3).collect();
    ViscoElements::new(
        SolidMeshRef::Tet4 {
            nodes_m: &[
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            elements: &[[0, 1, 2, 3]],
        },
        TetViscoField::Uniform(&material),
        &fixed,
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 1000,
        },
        cx,
    )
    .unwrap()
}
fn preparation() -> ViscoPoleBudget {
    ViscoPoleBudget {
        maximum_coordinates: 100,
        maximum_payload_bytes: 8 << 20,
        maximum_assembly_work: 1_000_000,
        maximum_evaluations: 1024,
        reference_rate_per_s: 2.0,
        physical_residual_tolerance: 2e-11,
    }
}
fn dynamic_budget() -> ViscoPoleDynamicBudget {
    ViscoPoleDynamicBudget {
        maximum_states: 20,
        maximum_ports: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_preparation_work: 1_000_000,
    }
}
fn response_budget() -> ViscoJointHarmonicBudget {
    ViscoJointHarmonicBudget {
        maximum_coordinates: 100,
        maximum_payload_bytes: 8 << 20,
        maximum_work: 100_000_000,
    }
}
fn pole_options() -> PoleOptions {
    PoleOptions {
        spectral_scale: 2.0,
        step_tolerance: 2e-12,
        residual_tolerance: 2e-12,
        max_iterations: 40,
        max_dimension: 100,
        max_dense_bytes: 8 << 20,
        max_work: 100_000_000,
    }
}


fn coupled_budget() -> SolidMimoResonanceBudget {
    SolidMimoResonanceBudget {
        maximum_realization_bytes: 1 << 20,
        poles: preparation(),
        load: dynamic_budget(),
    }
}
fn original_model() -> MimoModel {
    MimoModel {
        ports: 2,
        poles: vec![
            PoleTerm::Real { pole: -2.0, residue: 0.0 },
            PoleTerm::Pair { pole: C64::new(-0.4, 1.8), residue: C64::ZERO },
        ],
        residues: vec![
            [0.8, -0.25, -0.25, 0.6].map(C64::from_re).to_vec(),
            [0.4, 0.12, 0.12, 0.3].map(C64::from_re).to_vec(),
        ],
        d: vec![0.9, 0.15, 0.15, 0.7],
        e: vec![0.0; 4],
        reciprocal: true,
    }
}
fn original_load(model: &MimoModel, dt: f64) -> PreparedMimoLoad {
    PreparedMimoLoad::new(
        model,
        dt,
        MimoLoadLimits {
            band_rad_s: [0.1, 100.0],
            sampled_passivity_points: 33,
            max_states: 32,
            max_memory_bytes: 1 << 20,
        },
    ).unwrap()
}
fn near(a: C64, b: C64, scale: f64) {
    assert!((a-b).abs() < 5e-10*scale.max(a.abs()).max(b.abs()), "{a:?} != {b:?}");
}

#[test]
fn g1_original_real_and_pair_states_mutual_impedance_and_history_survive_attachment() {
    with_cx(|cx, _| {
        let body = source(true, cx);
        let model = original_model();
        let mut load = original_load(&model, 0.05);
        let control = original_load(&model, 0.1);
        for input in [[1.0, -0.4], [0.2, 0.3], [-0.5, 0.1]] {
            load.step_mean_into(&input, &mut [0.0; 2]).unwrap();
        }
        let accepted = load.pole_states().to_vec();
        assert!(accepted.iter().any(|z| z.abs() > 0.0));
        assert_ne!(load.mean_direct_term(), model.d.as_slice());
        let map = [1.0, -0.4];
        let prepared = PreparedSolidMimoResonance::new(
            &body, &load, &map, coupled_budget(), cx,
        ).unwrap();
        let cold = PreparedSolidMimoResonance::new(
            &body, &control, &map, coupled_budget(), cx,
        ).unwrap();
        assert!(std::ptr::eq(prepared.source(), &load));
        assert!(std::ptr::eq(prepared.reports().0, load.reports().0));
        assert!(std::ptr::eq(prepared.reports().1, load.reports().1));
        assert_eq!(prepared.source().source_model(), &model);
        assert_eq!(prepared.source().time_step_s(), 0.05);
        assert_eq!(prepared.model().dimension(), 13);
        assert_eq!(prepared.model().state_space_load_range(), Some(7..13));
        let matrices = load.continuous_operator()
            .state_space(1 << 20, &mut || false).unwrap();
        let force = C64::new(1e-4, -2e-4);
        for omega in [0.2, 0.9, 1.8, 2.5, 12.0] {
            let s = C64::new(0.0, -omega);
            // Independent original pole-residue transfer; the production
            // joint solve uses all original real state-space coordinates.
            let impedance = model.eval(s);
            let mut projected = C64::ZERO;
            for i in 0..2 {
                for j in 0..2 {
                    projected = projected + impedance[2*i+j].scale(map[i]*map[j]);
                }
            }
            let st = s.scale(0.25);
            let host = C64::from_re(2.5)+s*s+s.scale(0.2)
                +st.scale(6.5)/(C64::ONE+st);
            let expected = force/(host+s*projected);
            let factor = prepared.harmonic_factor(omega, response_budget(), cx).unwrap();
            let r = factor.solve(&[force], cx).unwrap();
            let q = st*expected/(C64::ONE+st);
            let boundary = r.state_space_load.as_ref().unwrap();
            near(r.displacement_m[0], expected, force.abs());
            near(r.physical_vector[1], q, force.abs());
            assert_eq!(boundary.state.len(), 6);
            for i in 0..matrices.n {
                let mut residual = s*boundary.state[i];
                for j in 0..matrices.n {
                    residual = residual - boundary.state[j].scale(matrices.a[i*matrices.n+j]);
                }
                for p in 0..2 {
                    residual = residual - boundary.flow[p].scale(matrices.b[i*2+p]);
                }
                near(residual, C64::ZERO, force.abs());
            }
            for i in 0..2 {
                near(boundary.flow[i], (s*expected).scale(map[i]), force.abs());
                let effort = (0..2).fold(C64::ZERO, |sum, j| {
                    sum + impedance[2*i+j]*boundary.flow[j]
                });
                near(boundary.effort[i], effort, force.abs());
            }
            let velocity = s*expected;
            let heat = 0.5*(0.2*velocity.norm_sq()+(6.5/0.25)*q.norm_sq());
            let input = 0.5*(force.conj()*velocity).re;
            let exchange = 0.5*boundary.flow.iter().zip(&boundary.effort)
                .map(|(v,e)| (v.conj()* *e).re).sum::<f64>();
            near(C64::from_re(r.average_material_heat_w), C64::from_re(heat), 1e-25);
            near(C64::from_re(input-heat), C64::from_re(exchange), 1e-25);
            assert_eq!(r.average_total_storage_j, None);
            assert_eq!(r.total_power_balance_residual_w, None);
            let other = cold.harmonic_factor(omega, response_budget(), cx).unwrap()
                .solve(&[force], cx).unwrap();
            // Source clocks and accepted transient history do not warp the
            // continuous response or become a second periodic drive.
            assert_eq!(r.physical_vector, other.physical_vector);
            assert_eq!(load.pole_states(), accepted);
        }
        assert_eq!(load.source_model(), &model);
        assert_eq!(load.time_step_s(), 0.05);
        assert_eq!(load.pole_states(), accepted);
    });
}

#[test]
fn g1_original_rational_source_has_the_analytic_loaded_finite_poles() {
    with_cx(|cx, _| {
        let body = source(false, cx);
        let real_root = 1.8_f64;
        let decay = 0.2_f64;
        let frequency_squared = 5.0/real_root;
        let direct = real_root+2.0*decay-2.0;
        let residue = frequency_squared+2.0*real_root*decay-(2.5+2.0*direct);
        // Original joint cubic factors as
        // (s+1.8)*(s^2+0.4s+5/1.8).
        let source_model = MimoModel {
            ports: 1,
            poles: vec![PoleTerm::Real { pole: -2.0, residue: 0.0 }],
            residues: vec![vec![C64::from_re(residue)]],
            d: vec![direct], e: vec![0.0], reciprocal: true,
        };
        let load = original_load(&source_model, 0.05);
        let prepared = PreparedSolidMimoResonance::new(
            &body, &load, &[1.0], coupled_budget(), cx,
        ).unwrap();
        for expected in [
            C64::from_re(-real_root),
            C64::new(-decay, (frequency_squared-decay*decay).sqrt()),
        ] {
            let seed = expected+C64::new(-0.01, 0.02);
            let pole = prepared.model().resonance_pole_from_state(
                seed, &[C64::ONE, seed/(seed+C64::from_re(2.0))],
                &pole_options(), cx,
            ).unwrap();
            near(pole.rate_per_s, expected, 1.0);
            let response = pole.state_space_load.as_ref().unwrap();
            near(response.state[0],
                pole.rate_per_s*pole.physical_vector[0]/(pole.rate_per_s+C64::from_re(2.0)),
                1.0);
            let scale = 2.0*pole.rate_per_s.abs()*response.body_storage_j
                +response.absorbed_power_w.abs()+1e-30;
            assert!(response.body_decay_power_residual_w.abs() < 5e-10*scale);
        }
        assert!(load.pole_states().iter().all(|z| *z == C64::ZERO));
    });
}

#[test]
fn g4_original_admission_budgets_band_and_cancellation_preserve_source() {
    with_cx(|cx, gate| {
        let body = source(true, cx);
        let model = original_model();
        let mut load = original_load(&model, 0.05);
        load.step_mean_into(&[1.0, -0.5], &mut [0.0; 2]).unwrap();
        let before = load.pole_states().to_vec();
        for map in [&[][..], &[1.0][..], &[1.0,f64::NAN][..]] {
            assert!(PreparedSolidMimoResonance::new(
                &body, &load, map, coupled_budget(), cx,
            ).is_err());
        }
        for budget in [
            SolidMimoResonanceBudget { maximum_realization_bytes: 1, ..coupled_budget() },
            SolidMimoResonanceBudget {
                poles: ViscoPoleBudget { maximum_coordinates: 12, ..preparation() },
                ..coupled_budget()
            },
            SolidMimoResonanceBudget {
                load: ViscoPoleDynamicBudget { maximum_states: 5, ..dynamic_budget() },
                ..coupled_budget()
            },
            SolidMimoResonanceBudget {
                load: ViscoPoleDynamicBudget { maximum_ports: 1, ..dynamic_budget() },
                ..coupled_budget()
            },
        ] {
            assert!(PreparedSolidMimoResonance::new(&body,&load,&[1.0,-0.4],budget,cx).is_err());
        }
        let prepared = PreparedSolidMimoResonance::new(
            &body, &load, &[1.0,-0.4], coupled_budget(), cx,
        ).unwrap();
        for omega in [0.05, 100.1, f64::NAN] {
            assert!(prepared.harmonic_factor(omega,response_budget(),cx).is_err());
        }
        gate.request();
        assert!(matches!(
            PreparedSolidMimoResonance::new(&body,&load,&[],coupled_budget(),cx),
            Err(SolidMimoResonanceError::Cancelled),
        ));
        assert!(prepared.harmonic_factor(1.0,response_budget(),cx).is_err());
        assert_eq!(load.pole_states(), before);
    });
}
