use super::*;
use fs_couple::cohesive_contact::patch::CohesivePatchConfig;
use fs_couple::render::schedule::force::coupled::contact::cohesive::patch::{
    CohesiveModalSite, CohesivePatchModalConfig, CohesivePatchModalSystem,
};

fn patch_config() -> CohesivePatchModalConfig {
    let mut contact = config().contact;
    contact.max_material_evaluations = 3000;
    CohesivePatchModalConfig {
        patch: CohesivePatchConfig {
            max_points: 8,
            contact,
        },
        max_setup_terms: 100_000,
        force_absolute_tolerance_n: config().force_absolute_tolerance_n,
        force_relative_tolerance: config().force_relative_tolerance,
    }
}

#[test]
fn coincident_quadrature_preserves_original_bond_response_and_reference_area() {
    let law = law().with_viscosity(CohesiveViscosity::new(2e6, 7e5, 100.0).unwrap());
    let relative = [0.05, 0.04, -0.025];
    let mut one = system(&law, relative);
    let sites = [0.2, 0.3, 0.5]
        .into_iter()
        .map(|fraction| CohesiveModalSite::new(&law, AREA * fraction, columns()).unwrap())
        .collect();
    let mut patch = CohesivePatchModalSystem::new(
        network(relative, 1.0, false),
        sites,
        patch_config(),
        &CancelGate::new(),
    )
    .unwrap();
    let initial = patch.total_energy_j().unwrap();
    let mut dissipated = 0.0;
    for _ in 0..150 {
        let a = *one.step(&[0.0; 6], [0.0; 3]).unwrap();
        let b = *patch.step(&[0.0; 6], &[[0.0; 3]; 3]).unwrap();
        close(a.cohesive_energy_j, b.cohesive_energy_j, 2e-12);
        close(a.fracture_dissipation_j, b.fracture_dissipation_j, 2e-12);
        close(a.viscous_dissipation_j, b.viscous_dissipation_j, 2e-12);
        for (left, right) in one.components().iter().zip(patch.components()) {
            for (x, y) in left.states().iter().zip(right.states()) {
                close(x.displacement_m_sqrt_kg, y.displacement_m_sqrt_kg, 2e-12);
                close(
                    x.velocity_m_sqrt_kg_per_s,
                    y.velocity_m_sqrt_kg_per_s,
                    2e-11,
                );
            }
        }
        for site in patch.last_point_frames().unwrap() {
            close(site.damage, a.damage, 2e-10);
            for i in 0..3 {
                close(site.jump_m[i], a.jump_m[i], 2e-12);
                close(
                    site.mean_force_n[i],
                    a.mean_force_n[i] * site.area_m2 / AREA,
                    3e-9,
                );
            }
        }
        dissipated += b.fracture_dissipation_j + b.viscous_dissipation_j + b.network_dissipation_j;
        close(patch.total_energy_j().unwrap() + dissipated, initial, 3e-11);
        assert!(b.material_evaluations <= patch_config().patch.contact.max_material_evaluations);
    }
}

// Two rigid bars retain translation in three axes and pitch about y. In this
// mass-normalized basis q_pitch=sqrt(I)*theta. A point at x has z motion -x*theta;
// the SAME lever arm maps its normal force into physical torque on both bodies.
const IA: f64 = 0.003;
const IB: f64 = 0.004;
const LOCATIONS: [f64; 3] = [-0.02, 0.0, 0.02];

fn bar_columns(x: f64) -> [Vec<f64>; 3] {
    let mut columns: [Vec<f64>; 3] = core::array::from_fn(|_| vec![0.0; 8]);
    for i in 0..3 {
        columns[i][i] = -1.0 / det::sqrt(MA);
        columns[i][4 + i] = 1.0 / det::sqrt(MB);
    }
    columns[0][3] = x / det::sqrt(IA);
    columns[0][7] = -x / det::sqrt(IB);
    columns
}

fn bars(relative_velocity: [f64; 4]) -> CoupledModalSystem {
    let models = [(MA, IA), (MB, IB)]
        .into_iter()
        .enumerate()
        .map(|(side, (mass, inertia))| {
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                RATE,
                vec![
                    ModalAcousticMode {
                        angular_frequency_rad_s: 0.0,
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: C64::ZERO
                    };
                    4
                ],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
            let states: Vec<_> = (0..4)
                .map(|i| {
                    let (own, other) = if i == 3 {
                        (inertia, if side == 0 { IB } else { IA })
                    } else {
                        (mass, if side == 0 { MB } else { MA })
                    };
                    ModalAcousticState {
                        displacement_m_sqrt_kg: 0.0,
                        velocity_m_sqrt_kg_per_s: det::sqrt(own)
                            * (if side == 0 {
                                -other / (own + other)
                            } else {
                                other / (own + other)
                            })
                            * relative_velocity[i],
                    }
                })
                .collect();
            model.restore_states(&states).unwrap();
            model
        })
        .collect();
    CoupledModalSystem::new(
        models,
        vec![],
        ModalCouplingConfig {
            max_modes: 8,
            max_connections: 0,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1.0,
            maximum_abs_pressure_pa: 1000.0,
            maximum_abs_connection_force_n: 1e5,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-10,
        },
        &CancelGate::new(),
    )
    .unwrap()
}

fn bar_patch<'a>(law: &'a BilinearCohesive, velocity: [f64; 4]) -> CohesivePatchModalSystem<'a> {
    let sites = LOCATIONS
        .into_iter()
        .map(|x| CohesiveModalSite::new(law, AREA / 3.0, bar_columns(x)).unwrap())
        .collect();
    CohesivePatchModalSystem::new(bars(velocity), sites, patch_config(), &CancelGate::new())
        .unwrap()
}

#[test]
fn finite_bond_peels_at_the_opening_edge_and_transmits_equal_opposite_force_and_torque() {
    let law = law();
    let mut patch = bar_patch(&law, [0.008, 0.01, -0.005, -2.0]);
    let initial = patch.total_energy_j().unwrap();
    let mut loss = 0.0;
    let mut edge_ahead = false;
    let mut old_velocity = [0.0; 8];
    for _ in 0..250 {
        for (i, state) in patch
            .components()
            .iter()
            .flat_map(|x| x.states())
            .enumerate()
        {
            old_velocity[i] = state.velocity_m_sqrt_kg_per_s;
        }
        let frame = *patch.step(&[0.0; 8], &[[0.0; 3]; 3]).unwrap();
        let responses = patch.last_point_frames().unwrap();
        let total_force: [f64; 3] =
            core::array::from_fn(|i| responses.iter().map(|r| r.mean_force_n[i]).sum());
        let moment: f64 = responses
            .iter()
            .zip(LOCATIONS)
            .map(|(r, x)| x * r.mean_force_n[0])
            .sum();
        let dt = patch.sample_period_s();
        for i in 0..3 {
            let va = patch.components()[0].states()[i].velocity_m_sqrt_kg_per_s;
            let vb = patch.components()[1].states()[i].velocity_m_sqrt_kg_per_s;
            close(det::sqrt(MA) * va + det::sqrt(MB) * vb, 0.0, 2e-13);
            close(
                det::sqrt(MA) * (va - old_velocity[i]),
                dt * total_force[i],
                2e-13,
            );
            close(
                det::sqrt(MB) * (vb - old_velocity[4 + i]),
                -dt * total_force[i],
                2e-13,
            );
        }
        let wa = patch.components()[0].states()[3].velocity_m_sqrt_kg_per_s;
        let wb = patch.components()[1].states()[3].velocity_m_sqrt_kg_per_s;
        close(det::sqrt(IA) * wa + det::sqrt(IB) * wb, 0.0, 2e-14);
        close(det::sqrt(IA) * (wa - old_velocity[3]), -dt * moment, 2e-14);
        close(det::sqrt(IB) * (wb - old_velocity[7]), dt * moment, 2e-14);
        edge_ahead |= responses[2].damage > responses[0].damage + 0.05
            && responses[2].damage > responses[1].damage + 0.02;
        loss += frame.network_dissipation_j
            + frame.fracture_dissipation_j
            + frame.viscous_dissipation_j;
        close(patch.total_energy_j().unwrap() + loss, initial, 3e-11);
    }
    assert!(
        edge_ahead,
        "the opening edge must peel before the compressed/central sites"
    );
    let s = patch.compliance_m_per_n();
    // Distinct contact locations remain mechanically coupled through translation
    // and pitch; replacing this by three independent solves is a different model.
    assert!(s[6].abs() > 0.0 && s[3].abs() > 0.0);
    assert_ne!(s[0], s[3 * 9 + 3]);
}

#[test]
fn late_local_refusal_and_cancellation_preserve_the_entire_patch_and_replay() {
    let ordinary = law();
    let viscous = law().with_viscosity(CohesiveViscosity::new(2e6, 7e5, 0.02).unwrap());
    let build = || {
        let sites = vec![
            CohesiveModalSite::new(&ordinary, AREA, columns()).unwrap(),
            CohesiveModalSite::new(&viscous, AREA, columns()).unwrap(),
        ];
        CohesivePatchModalSystem::new(
            network([0.0; 3], 1.0, false),
            sites,
            patch_config(),
            &CancelGate::new(),
        )
        .unwrap()
    };
    let mut owner = build();
    let mut reference = build();
    let zero = [[0.0; 3]; 2];
    owner.step(&[0.0; 6], &zero).unwrap();
    reference.step(&[0.0; 6], &zero).unwrap();
    let before_points = owner.points().to_vec();
    let before_observations = owner.last_point_frames().unwrap().to_vec();
    let before_frame = *owner.last_frame().unwrap();
    let before_states: Vec<_> = owner
        .components()
        .iter()
        .flat_map(|x| x.states().iter().copied())
        .collect();
    assert!(owner.step(&[0.0; 6], &[[0.0; 3], [1.0, 0.0, 0.0]]).is_err());
    assert_eq!(owner.points(), before_points);
    assert_eq!(owner.last_point_frames().unwrap(), before_observations);
    assert_eq!(*owner.last_frame().unwrap(), before_frame);
    assert_eq!(
        owner
            .components()
            .iter()
            .flat_map(|x| x.states().iter().copied())
            .collect::<Vec<_>>(),
        before_states
    );
    let gate = CancelGate::new();
    gate.request();
    assert!(owner.step_under_gate(&[0.0; 6], &zero, &gate).is_err());
    assert_eq!(owner.points(), before_points);
    let drive = [[0.001, 0.0, 0.0], [0.001, 0.0, 0.0]];
    assert_eq!(
        owner.step(&[0.0; 6], &drive).unwrap(),
        reference.step(&[0.0; 6], &drive).unwrap()
    );
    assert_eq!(owner.last_point_frames(), reference.last_point_frames());
}
