//! Exact P1 recoil relaxation and a stationary magnet driving a moving conductor.
use super::*;
use fs_em::planar::eddy::{EddyControls, EddyError, PlanarEddy, moving::MovingEddyControls};
#[path = "../winding_ohmic/common.rs"]
mod ohmic;

fn controls() -> EddyControls {
    EddyControls {
        linear: common::config(),
        energy_tolerance_j: 1e-9,
        current_balance_tolerance_a: 1e-9,
    }
}
fn moving_controls() -> MovingEddyControls {
    MovingEddyControls {
        relative_tolerance: 1e-10,
        absolute_residual_a_m: 1e-12,
        restart: 4,
        max_cycles: 20,
        max_krylov_scalars: 1000,
        max_cell_magnetic_peclet: 20.0,
        energy_tolerance_j: 1e-9,
        current_balance_tolerance_a: 1e-9,
    }
}
fn conductors(source: &ResolvedConductorState, count: usize) -> Vec<EddyConductor> {
    (0..count)
        .map(|triangle| EddyConductor {
            triangle,
            region: 7,
            conductor: source.clone(),
        })
        .collect()
}

#[test]
fn midpoint_recoil_retains_interface_and_exterior_sources_with_second_order_decay() {
    let (mu, br, sigma) = (1.0, 0.5, 18.0);
    for natural_boundary in [false, true] {
        let (law, conductor) = recoil(mu, br, [-5.0, 5.0], 5.0, sigma);
        let directions = if natural_boundary {
            vec![[1.0, 0.0]]
        } else {
            vec![[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]]
        };
        let models = directions
            .iter()
            .map(|&direction| {
                PlanarMagneticMaterial::Recoil(
                    PlanarRecoilMaterial::new(law.clone(), direction).unwrap(),
                )
            })
            .collect();
        let problem = PlanarMagnetostatics::new_with_material_models(
            if natural_boundary {
                triangle()
            } else {
                induction::mesh(0.5)
            },
            models,
            common::limits(),
            &mut || false,
        )
        .unwrap();
        let boundary = if natural_boundary {
            vec![(0, 0.0), (1, 0.0)]
        } else {
            induction::boundary(0.0)
        };
        let count = directions.len();
        let initial = problem
            .solve(
                &vec![if natural_boundary { 3.0 } else { 6.0 }; count],
                &boundary,
                common::config(),
                &mut || false,
            )
            .unwrap();
        // Elementary triangle/fan integrals, independent of FEEC assembly:
        // C_sigma is the consistent mass with its floating mean removed.
        let (free, equilibrium, k, c) = if natural_boundary {
            (2, br, 0.1 / mu, sigma * 0.1 / 18.0)
        } else {
            (4, br / 2.0, 4.0 * 0.5 / mu, sigma * 0.5 / 18.0)
        };
        let a0 = initial.potential_wb_m()[free];
        let original = PlanarEddy::new(
            &initial,
            conductors(&conductor, count),
            &[],
            &boundary,
            0.0,
            controls(),
            &mut || false,
        )
        .unwrap();
        let equilibrium_field = problem
            .solve(&vec![0.0; count], &boundary, common::config(), &mut || {
                false
            })
            .unwrap();
        common::near(equilibrium_field.potential_wb_m()[free], equilibrium);
        let mut errors = Vec::new();
        for steps in [4, 8] {
            let mut state = original.clone();
            let dt = 0.5 / f64::from(steps);
            let factor = (c - 0.5 * dt * k) / (c + 0.5 * dt * k);
            let mut expected = a0;
            for index in 1..=steps {
                let before = expected;
                expected = equilibrium + factor * (before - equilibrium);
                let candidate = state
                    .step_to_candidate(
                        f64::from(index) * dt,
                        state.conductors().to_vec(),
                        &vec![0.0; count],
                        controls(),
                        &mut || false,
                    )
                    .unwrap();
                let step = candidate.step();
                common::near(candidate.state().potential_wb_m()[free], expected);
                common::near(
                    step.stored_before_j,
                    0.5 * k * (before - equilibrium).powi(2),
                );
                common::near(
                    step.stored_after_j,
                    0.5 * k * (expected - equilibrium).powi(2),
                );
                common::near(step.heat_j, c * (expected - before).powi(2) / dt);
                common::near(step.parameter_work_j, 0.0);
                common::near(step.supplied_work_j, 0.0);
                common::near(
                    step.stored_after_j - step.stored_before_j + step.heat_j,
                    0.0,
                );
                for cell in &step.magnetic_cells {
                    let remanence = directions[cell.triangle].map(|x| br * x);
                    for endpoint in 0..2 {
                        let b = cell.flux_density_t[endpoint];
                        common::near(
                            cell.free_energy_j[endpoint],
                            cell.volume_m3 / (2.0 * mu)
                                * ((b[0] - remanence[0]).powi(2) + (b[1] - remanence[1]).powi(2)),
                        );
                    }
                }
                common::near(step.region_currents_a[0].1, 0.0);
                state = candidate.into_parts().0;
            }
            let exact = equilibrium + (a0 - equilibrium) * (-k * 0.5 / c).exp();
            errors.push((state.potential_wb_m()[free] - exact).abs());
        }
        assert!(errors[0] / errors[1] > 3.8 && errors[0] / errors[1] < 4.4);

        // Replacing the conductive source or erasing its magnetic role cannot
        // be passed off as a fixed-law continuation or as physical Joule heat.
        let (_, foreign) = recoil(mu, br, [-5.0, 5.0], 5.0, 2.0 * sigma);
        assert!(
            original
                .step_to_candidate(
                    0.1,
                    conductors(&foreign, count),
                    &vec![0.0; count],
                    controls(),
                    &mut || false,
                )
                .is_err()
        );
        assert!(matches!(
            original.step_with_materials_candidate(
                0.1,
                vec![common::material(mu); count],
                original.conductors().to_vec(),
                &vec![0.0; count],
                common::limits(),
                controls(),
                &mut || false,
            ),
            Err(EddyError::Invalid(
                "fixed-remanence induction cannot replace magnetic materials"
            ))
        ));
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.potential_wb_m(), initial.potential_wb_m());
    }
}

#[test]
fn fixed_magnet_drives_moving_conductor_drag_with_original_recoil_energy() {
    let (mu, br, sigma, depth) = (1.0, 0.5, 18.0, 0.5);
    let (law, _) = recoil(mu, br, [-5.0, 5.0], 5.0, sigma);
    let basis = ElectricalPropertyBasis::Conductivity;
    let source = ohmic::card(1.0 / sigma, 0.0, mu, basis);
    let material = ohmic::bundle(&source, 300.0, basis);
    let conductor =
        ResolvedConductorState::try_from_material_with_basis(material.clone(), basis).unwrap();
    // Bottom magnet and top conductor only meet at the fixed-zero-velocity
    // centre vertex. The side triangles are insulating induced material.
    let problem = PlanarMagnetostatics::new_with_material_models(
        induction::mesh(depth),
        vec![
            PlanarMagneticMaterial::Recoil(PlanarRecoilMaterial::new(law, [1.0, 0.0]).unwrap()),
            PlanarMagneticMaterial::Linear(common::material(mu)),
            PlanarMagneticMaterial::Linear(material),
            PlanarMagneticMaterial::Linear(common::material(mu)),
        ],
        common::limits(),
        &mut || false,
    )
    .unwrap();
    let boundary = induction::boundary(0.0);
    let initial = problem
        .solve(&[0.0; 4], &boundary, common::config(), &mut || false)
        .unwrap();
    let original = PlanarEddy::new(
        &initial,
        vec![EddyConductor {
            triangle: 2,
            region: 7,
            conductor,
        }],
        &[0, 1, 3],
        &boundary,
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let (a0, k, c) = (br / 8.0, 4.0 * depth / mu, sigma * depth / 72.0);
    let bias_energy = 3.0 * depth * br * br / (32.0 * mu);
    common::near(initial.potential_wb_m()[4], a0);
    common::near(initial.energy_j(), bias_energy);
    let dt = 0.125;
    for speed in [0.0, 0.2, -0.2] {
        let mut velocity = vec![[0.0; 2]; 5];
        velocity[2][1] = speed;
        velocity[3][1] = speed;
        // On the conductor, v_y = speed*(1-phi_c), B_x = -2*A_c:
        // its projected advection is exactly 2*speed*C_sigma*A_c.
        let advection = 2.0 * speed * c;
        let delta = -dt * advection * a0 / (c + 0.5 * dt * (k + advection));
        let stage = a0 + 0.5 * delta;
        let effective_rate = delta / dt + 2.0 * speed * stage;
        let candidate = original
            .step_with_eulerian_velocity_candidate(
                dt,
                &velocity,
                &[0.0; 4],
                moving_controls(),
                &mut || false,
            )
            .unwrap();
        let step = candidate.step();
        common::near(candidate.state().potential_wb_m()[4], a0 + delta);
        common::near(step.stored_before_j, bias_energy);
        common::near(step.stored_after_j, bias_energy + 0.5 * k * delta * delta);
        common::near(step.heat_j, dt * c * effective_rate * effective_rate);
        common::near(
            step.mechanical_work_j,
            -dt * advection * stage * effective_rate,
        );
        common::near(step.supplied_work_j, 0.0);
        common::near(step.energy_residual_j, 0.0);
        common::near(step.region_currents_a[0].1, 0.0);
        for (actual, multiplier) in step.cells[0]
            .joule
            .electric_field_v_m
            .iter()
            .zip([1.0, 1.0, -2.0])
        {
            common::near(*actual, multiplier * effective_rate / 3.0);
        }
        if speed == 0.0 {
            let stationary = original
                .step_to_candidate(
                    dt,
                    original.conductors().to_vec(),
                    &[0.0; 4],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            assert_eq!(
                candidate.state().potential_wb_m(),
                stationary.state().potential_wb_m()
            );
            common::near(step.heat_j, 0.0);
            common::near(step.mechanical_work_j, 0.0);
        } else {
            assert!(step.heat_j > 0.0 && step.mechanical_work_j < 0.0);
        }
    }
}

#[test]
fn midpoint_cannot_overshoot_the_original_recoil_envelope() {
    let (law, conductor) = recoil(1.0, 0.5, [0.0, 5.0], 5.0, 18.0);
    let problem = PlanarMagnetostatics::new_with_material_models(
        triangle(),
        vec![PlanarMagneticMaterial::Recoil(
            PlanarRecoilMaterial::new(law, [1.0, 0.0]).unwrap(),
        )],
        common::limits(),
        &mut || false,
    )
    .unwrap();
    let boundary = [(0, 0.0), (1, 0.0)];
    let initial = problem
        .solve(&[3.0], &boundary, common::config(), &mut || false)
        .unwrap();
    let original = PlanarEddy::new(
        &initial,
        conductors(&conductor, 1),
        &[],
        &boundary,
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    // tau = mu*sigma/18 = 1 s. Midpoint at dt=3 s changes H=1 to H=-0.2,
    // outside this authored [0,5] A/m recoil branch. No clamping is allowed.
    assert!(
        original
            .step_to_candidate(
                3.0,
                original.conductors().to_vec(),
                &[0.0],
                controls(),
                &mut || false,
            )
            .is_err()
    );
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(original.potential_wb_m(), initial.potential_wb_m());
}
