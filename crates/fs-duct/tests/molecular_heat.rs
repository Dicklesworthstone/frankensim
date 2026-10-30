//! Positive molecular heat on the original held-port DG trajectory.
use fs_duct::{
    Duct, DuctError, Segment,
    dg::{DgBoundary, DgDuct, DgMesh, DgState, DgWallResolution, PortId, PortInput},
};
use fs_material::{
    gas::{GasSpec, GasState},
    relaxation::{
        InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
        RelaxingMedium, SpeciesRelaxation, VibrationalMode,
    },
};
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_thermochem::SpeciesId;

fn budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 128 * 1024 * 1024,
        max_multiply_terms: 100_000_000_000,
    }
}
fn duct(diffusion: usize, rate_scale: f64, active: bool) -> DgDuct {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let id = SpeciesId::new("synthetic").unwrap();
    let law = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: id.clone(),
        processes: if active {
            vec![RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1300.0,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 2.0 * rate_scale,
                    t_ref_k: gas.temperature,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 500.0,
                },
            }]
        } else {
            Vec::new()
        },
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(id, 1.0)]).unwrap())
    .unwrap();
    let mut dg = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.001,
                length: 0.05,
            }],
        },
        &gas,
        DgMesh {
            elements: 1,
            order: 1,
        },
        DgBoundary::Port(PortId(7)),
        DgBoundary::Closed,
    )
    .unwrap();
    dg = match diffusion {
        0 => dg,
        1 => dg.with_viscothermal(8, [50.0, 2000.0], 0.02).unwrap(),
        _ => dg
            .with_spatial_viscothermal(DgWallResolution {
                poles_per_field: 8,
                quadrature_points: 2,
                band_rad_s: [50.0, 2000.0],
                maximum_relative_error: 0.02,
            })
            .unwrap(),
    };
    dg.with_molecular_relaxation(law).unwrap()
}
fn quadrature(dg: &DgDuct, initial: &DgState, ports: &[PortInput], dt: f64, count: usize) -> f64 {
    let step = dg.propagator(dt / count as f64, 1e-12).unwrap();
    let mut state = initial.clone();
    let mut sum = dg.molecular_dissipated_power_w(&state).unwrap();
    for k in 1..=count {
        step.step(&mut state, ports).unwrap();
        let weight = if k == count {
            1.0
        } else if k % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * dg.molecular_dissipated_power_w(&state).unwrap();
    }
    sum * dt / (3.0 * count as f64)
}

#[test]
fn molecular_heat_integrates_native_forcing_and_local_diffusion_without_altering_motion() {
    for diffusion in 0..=2 {
        let dg = duct(diffusion, 1.0, true);
        let step = dg.propagator(2e-5, 1e-12).unwrap();
        let mut initial = dg.project_initial_state(|x| 200.0 + 100.0 * x, |_| 1e-8);
        for _ in 0..3 {
            step.step(&mut initial, &[]).unwrap();
        }
        let frozen = initial.clone();
        let generator = dg.continuous_operator().generator_per_s().to_vec();
        let ports = [PortInput {
            id: PortId(7),
            effort: 3e-8,
        }];
        let heat = dg
            .prepare_molecular_heat(&step, budget(), &mut || false)
            .unwrap();
        let actual = heat
            .integrated_heat_j(&initial, &ports, &mut || false)
            .unwrap();
        let coarse = quadrature(&dg, &initial, &ports, step.dt(), 128);
        let reference = quadrature(&dg, &initial, &ports, step.dt(), 256);
        assert!(actual > 0.0);
        assert!(
            (actual - reference).abs() < 2e-9 * actual,
            "diffusion={diffusion}: actual={actual:e}, independent={reference:e}"
        );
        assert!((actual - reference).abs() <= (actual - coarse).abs() + 3e-13 * actual);
        assert_eq!(initial, frozen);
        assert_eq!(generator, dg.continuous_operator().generator_per_s());
        assert_eq!(heat.operator_identity(), dg.operator_identity());
        assert_eq!(heat.duration(), step.dt());
        let mut a = initial.clone();
        let mut b = initial.clone();
        assert_eq!(
            step.step(&mut a, &ports).unwrap(),
            step.step(&mut b, &ports).unwrap()
        );
        assert_eq!(a, b);
        let total = a.cumulative_dissipated_j() - initial.cumulative_dissipated_j();
        if diffusion == 0 {
            assert!((total - actual).abs() < 2e-10 * initial.stored_energy_j());
        } else {
            assert!(total > actual); // Wall transport keeps its independent channel.
        }
    }
}

#[test]
fn source_cache_rejects_foreign_history_budgets_and_cancelled_observations() {
    let dg = duct(0, 1.0, true);
    let step = dg.propagator(1e-6, 1e-12).unwrap();
    let initial = dg.project_initial_state(|_| 100.0, |_| 0.0);
    let heat = dg
        .prepare_molecular_heat(&step, budget(), &mut || false)
        .unwrap();
    let foreign = duct(0, 2.0, true);
    assert!(
        dg.prepare_molecular_heat(
            &foreign.propagator(1e-6, 1e-12).unwrap(),
            budget(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        heat.integrated_heat_j(
            &foreign.project_initial_state(|_| 100.0, |_| 0.0),
            &[],
            &mut || false
        )
        .is_err()
    );
    assert!(
        heat.integrated_heat_j(
            &initial,
            &[PortInput {
                id: PortId(99),
                effort: 0.0
            }],
            &mut || false
        )
        .is_err()
    );
    assert!(
        dg.prepare_molecular_heat(
            &step,
            DenseExponentialBudget {
                max_multiply_terms: 0,
                ..budget()
            },
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        dg.prepare_molecular_heat(&step, budget(), &mut || true),
        Err(DuctError::Cancelled { .. })
    ));
    assert!(matches!(
        heat.integrated_heat_j(&initial, &[], &mut || true),
        Err(DuctError::Cancelled { .. })
    ));
    assert!(
        heat.integrated_heat_j(&initial, &[], &mut || false)
            .unwrap()
            > 0.0
    );
    let empty = duct(0, 1.0, false);
    let empty_step = empty.propagator(1e-6, 1e-12).unwrap();
    let empty_heat = empty
        .prepare_molecular_heat(&empty_step, budget(), &mut || false)
        .unwrap();
    assert_eq!(
        empty_heat
            .integrated_heat_j(
                &empty.project_initial_state(|_| 100.0, |_| 0.0),
                &[],
                &mut || false
            )
            .unwrap(),
        0.0
    );
}
