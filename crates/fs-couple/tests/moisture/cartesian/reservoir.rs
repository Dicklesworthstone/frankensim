//! Finite reservoir depletion/concentration reaches the actual 3D solid.
//! Authored source laws supply independent conservative reference solutions.
use super::super::{atmosphere as air, solution as bath, surface as emission};
use super::*;
use fs_couple::moisture::atmosphere::{CartesianMoistureParticipant, DryBasisSorption};

fn participant<'a>(
    id: u64,
    moisture: &'a MoistureCartesian,
    faces: &'a [Option<DryBasisSorption>; 6],
) -> CartesianMoistureParticipant<'a> {
    CartesianMoistureParticipant {
        id,
        moisture,
        faces,
        transport: config(),
        nonlinear: nonlinear(),
    }
}

fn solid(field: &MoistureCartesian, cx: &Cx<'_>) -> MoistureSolid {
    let geometry = field.geometry();
    let counts = geometry.cells.map(|n| n + 1);
    let nodes: usize = counts.iter().product();
    let fixed = (0..nodes)
        .filter(|i| i % counts[0] == 0)
        .flat_map(|i| [3 * i, 3 * i + 1, 3 * i + 2])
        .collect();
    MoistureSolid::from_cartesian(
        field,
        SolidStiffnessModel::Isotropic,
        MaterialPropertySelection::SingleClaimOnly,
        MaterialPropertySelection::SingleClaimOnly,
        HygroscopicAxes {
            principal_to_mesh: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            principal_frame: ContentHash([3; 32]),
            mesh_frame: geometry.frame,
        },
        fixed,
        MoistureSolidControls {
            assembly: TetAssemblyBudget {
                maximum_nodes: 64,
                maximum_tetrahedra: 96,
                maximum_free_dofs: 192,
                minimum_scaled_jacobian: 1e-8,
            },
            equilibrium: TetStaticSolveConfig {
                relative_residual_tolerance: 1e-10,
                maximum_iterations: 4096,
            },
            maximum_absolute_mandel_strain: 0.02,
            maximum_nodal_displacement_m: 0.01,
            mass_tolerance_kg: 1e-12,
        },
        cx,
    )
    .unwrap()
}

#[test]
fn finite_vapor_drives_three_dimensional_gradient_swelling_and_atomic_late_refusal() {
    with_cx(|cx| {
        let atmosphere = air::atmosphere(0.8, 0.01);
        let geometry = grid_geometry([0.02; 3], [2; 3]);
        let initial = grid(card([0.01; 2], 0.3), geometry, vec![0.1; 8], cx);
        let mechanical = solid(&initial, cx);
        let mut faces: [Option<DryBasisSorption>; 6] = std::array::from_fn(|_| None);
        faces[5] = Some(air::sorption(0.25));
        let (candidate, response) = mechanical
            .advance_cartesian_finite_atmosphere(
                &atmosphere,
                &initial,
                0.005,
                &faces,
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap();
        let report = candidate.step();
        let r0 = atmosphere.relative_humidity().unwrap();
        let capacity = report.vapor_mass_kg[0] / r0;
        // Solve both original z layers and finite vapor independently:
        // M_low=.1+(s*r-.1)/7; M_high=.1+3*(s*r-.1)/7.
        // Thus C*(r-r0)+m_dry*(2/7)*(s*r-.1)=0.
        let weight = initial.dry_mass_kg() * 2.0 / 7.0;
        let expected_rh = (capacity * r0 + weight * 0.1) / (capacity + weight * 0.25);
        near(
            candidate.atmosphere().relative_humidity().unwrap(),
            expected_rh,
            4e-10,
        );
        for i in 0..8 {
            let expected = 0.1 + if i < 4 { 1.0 } else { 3.0 } * (0.25 * expected_rh - 0.1) / 7.0;
            near(candidate.moisture().moisture()[i], expected, 5e-11);
            near(
                response.expansion_paths()[i].free_principal_strains()[2],
                0.02 * (expected - 0.1),
                2e-12,
            );
            near(
                report.transport.transport_start[i].query().axes()["moisture_content"],
                0.1,
                0.0,
            );
            near(
                report.transport.transport_end[i].query().axes()["moisture_content"],
                candidate.moisture().moisture()[i],
                0.0,
            );
        }
        assert_eq!(report.transport.surface_moisture[5].len(), 4);
        assert!(
            report.transport.boundary_water_kg[..5]
                .iter()
                .all(|v| *v == 0.0)
        );
        near(
            report.transport.boundary_water_kg[5],
            report.water_to_solid_kg,
            1e-13,
        );
        near(
            report.vapor_mass_kg[0] - report.vapor_mass_kg[1],
            report.water_to_solid_kg,
            1e-13,
        );
        near(
            air::carrier_mass(candidate.atmosphere()),
            air::carrier_mass(&atmosphere),
            1e-13,
        );
        near(
            response.wet_mass_kg(),
            candidate.moisture().wet_mass_kg(),
            1e-13,
        );
        assert!(
            response
                .solution()
                .displacement_m()
                .iter()
                .flatten()
                .any(|u| u.abs() > 1e-7)
        );
        assert!(candidate.moisture().moisture()[4] > candidate.moisture().moisture()[0]);

        let limited = grid(card([0.01; 2], 0.105), geometry, vec![0.1; 8], cx);
        let limited_solid = solid(&limited, cx);
        let water_only = atmosphere
            .advance_cartesian(
                &limited,
                0.005,
                &faces,
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap();
        assert!(water_only.moisture().moisture()[4] > 0.105);
        let before_air = atmosphere.clone();
        assert!(
            limited_solid
                .advance_cartesian_finite_atmosphere(
                    &atmosphere,
                    &limited,
                    0.005,
                    &faces,
                    config(),
                    nonlinear(),
                    air::controls(),
                    cx
                )
                .is_err()
        );
        assert_eq!(atmosphere, before_air);
        assert_eq!(limited.time_s(), 0.0);
        assert_eq!(limited.moisture(), &[0.1; 8]);
    });
}

#[test]
fn six_face_nonlinear_uptake_uses_one_vapor_debit_and_actual_endpoint_diffusivity() {
    with_cx(|cx| {
        let atmosphere = air::atmosphere(0.8, 1.0);
        let initial = grid(
            card([0.01, 0.04], 0.3),
            grid_geometry([0.1; 3], [1; 3]),
            vec![0.1],
            cx,
        );
        let faces = std::array::from_fn(|_| Some(air::sorption(0.25)));
        let candidate = atmosphere
            .advance_cartesian(
                &initial,
                1.0 / 30.0,
                &faces,
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap();
        let report = candidate.step();
        let r0 = atmosphere.relative_humidity().unwrap();
        let capacity = report.vapor_mass_kg[0] / r0;
        let dry = initial.dry_mass_kg();
        // M-.1 = 40*(.01+.1*M)*(A-B*M), with vapor eliminated exactly.
        let a = 0.25 * r0 + 0.25 * dry * 0.1 / capacity;
        let b = 1.0 + 0.25 * dry / capacity;
        let qa = 40.0 * 0.1 * b;
        let qb = 1.0 - 40.0 * (0.1 * a - 0.01 * b);
        let qc = -0.1 - 40.0 * 0.01 * a;
        let expected = -2.0 * qc / (qb + (qb * qb - 4.0 * qa * qc).sqrt());
        near(candidate.moisture().moisture()[0], expected, 2e-11);
        near(
            candidate.atmosphere().relative_humidity().unwrap(),
            r0 - dry * (expected - 0.1) / capacity,
            3e-10,
        );
        assert!(report.transport.nonlinear_iterations > 1);
        near(
            report.transport.transport_end[0]
                .property("moisture_diffusivity")
                .unwrap()
                .value_si(),
            0.01 + 0.1 * expected,
            3e-12,
        );
        for water in report.transport.boundary_water_kg {
            near(water, dry * (expected - 0.1) / 6.0, 2e-12);
        }
        near(report.water_balance_residual_kg, 0.0, 1e-13);
        assert_eq!(initial.time_s(), 0.0);
    });
}

#[test]
fn shared_cartesian_fields_adsorb_and_desorb_with_source_surface_resistance_without_order_bias() {
    with_cx(|cx| {
        let atmosphere = air::atmosphere(0.5, 1.0);
        let ga = grid_geometry([0.1, 0.05, 0.02], [1; 3]);
        let gb = grid_geometry([0.05, 0.02, 0.02], [1; 3]);
        let a = grid(card([0.01; 2], 0.3), ga, vec![0.02], cx);
        let b = grid(card([0.003; 2], 0.3), gb, vec![0.2], cx);
        let mut fa = std::array::from_fn(|_| Some(air::sorption(0.25)));
        fa[0] = Some(emission::law([0.0; 2]));
        let h = 0.05;
        let fb = std::array::from_fn(|_| Some(emission::law([h; 2])));
        let dt = 0.01;
        let run = |reverse| {
            let entries = if reverse {
                [participant(3, &b, &fb), participant(7, &a, &fa)]
            } else {
                [participant(7, &a, &fa), participant(3, &b, &fb)]
            };
            atmosphere
                .advance_shared_cartesian(&entries, dt, 2, air::controls(), cx)
                .unwrap()
        };
        let accepted = run(false);
        let replay = run(true);
        let capacity = accepted.step().vapor_mass_kg[0] / atmosphere.relative_humidity().unwrap();
        let alpha_a = (0..6)
            .filter(|side| *side != 0)
            .map(|side| 2.0 * dt * 0.01 / ga.lengths_m[side / 2].powi(2))
            .sum::<f64>();
        let alpha_b = (0..6)
            .map(|side| {
                let dx = gb.lengths_m[side / 2];
                dt / dx / (1.0 / h + dx / (2.0 * 0.003))
            })
            .sum::<f64>();
        let weights = [
            a.dry_mass_kg() * alpha_a / (1.0 + alpha_a),
            b.dry_mass_kg() * alpha_b / (1.0 + alpha_b),
        ];
        let expected_rh = (capacity * atmosphere.relative_humidity().unwrap()
            + weights[0] * 0.02
            + weights[1] * 0.2)
            / (capacity + 0.25 * (weights[0] + weights[1]));
        near(
            accepted.atmosphere().relative_humidity().unwrap(),
            expected_rh,
            2e-10,
        );
        assert_eq!(
            accepted
                .solids()
                .iter()
                .map(|body| body.id)
                .collect::<Vec<_>>(),
            vec![3, 7]
        );
        for (body, again) in accepted.solids().iter().zip(replay.solids()) {
            let (old, alpha) = if body.id == 7 {
                (0.02, alpha_a)
            } else {
                (0.2, alpha_b)
            };
            near(
                body.moisture.moisture()[0],
                (old + alpha * 0.25 * expected_rh) / (1.0 + alpha),
                2e-11,
            );
            assert_eq!(body.moisture.moisture(), again.moisture.moisture());
            assert_eq!(body.moisture.time_s(), dt);
        }
        assert!(accepted.solids()[0].moisture.moisture()[0] < 0.2);
        assert!(accepted.solids()[1].moisture.moisture()[0] > 0.02);
        assert_eq!(accepted.solids()[1].transport.boundary_water_kg[0], 0.0);
        assert!(
            accepted.solids()[0]
                .transport
                .surface_transfer
                .iter()
                .all(Option::is_some)
        );
        near(
            air::carrier_mass(accepted.atmosphere()),
            air::carrier_mass(&atmosphere),
            1e-13,
        );
        near(accepted.step().water_balance_residual_kg, 0.0, 1e-13);
        near(accepted.step().total_mass_residual_kg, 0.0, 1e-13);
        assert!(
            atmosphere
                .advance_shared_cartesian(
                    &[participant(7, &a, &fa), participant(3, &b, &fb)],
                    dt,
                    1,
                    air::controls(),
                    cx
                )
                .is_err()
        );
        let later = b
            .advance(dt, &sealed(), config(), nonlinear(), cx)
            .unwrap()
            .0;
        assert!(
            atmosphere
                .advance_shared_cartesian(
                    &[participant(7, &a, &fa), participant(3, &later, &fb)],
                    dt,
                    2,
                    air::controls(),
                    cx
                )
                .is_err()
        );
    });
}

#[test]
fn finite_solution_concentration_closes_all_cartesian_faces_and_native_swelling() {
    with_cx(|cx| {
        let final_volume = 5e-5;
        let geometry = grid_geometry([0.1; 3], [2; 3]);
        let dt = 0.01;
        let alpha = 24.0 * dt * 0.01 / 0.1_f64.powi(2);
        for (activity, m0) in [(0.998_f64, 0.198), (0.995_f64, 0.201)] {
            let initial = grid(card([0.01; 2], 0.3), geometry, vec![m0; 8], cx);
            let mechanical = solid(&initial, cx);
            // Every cell has three exposed half-cell faces; by symmetry the
            // internal flux is zero. Manufacture the endpoint solvent volume.
            let expected = (m0 + alpha * 0.2 * activity) / (1.0 + alpha);
            let water = initial.dry_mass_kg() * (expected - m0);
            let initial_volume = final_volume + water / 1000.0;
            let solute = -activity.ln() * final_volume / bath::WATER_VOLUME;
            let reservoir = bath::reservoir(initial_volume, vec![0.6 * solute, 0.4 * solute]);
            let faces = std::array::from_fn(|_| Some(air::sorption(0.2)));
            let (candidate, response) = mechanical
                .advance_cartesian_finite_solution(
                    &reservoir,
                    &initial,
                    dt,
                    &faces,
                    config(),
                    nonlinear(),
                    bath::controls(initial_volume),
                    cx,
                )
                .unwrap();
            near(
                candidate.reservoir().solution().solvent_volume_m3(),
                final_volume,
                2e-16,
            );
            near(
                candidate.reservoir().solution().water_activity().unwrap(),
                activity,
                1e-12,
            );
            assert_eq!(
                candidate.reservoir().solution().solute_moles(),
                reservoir.solution().solute_moles()
            );
            for &m in candidate.moisture().moisture() {
                near(m, expected, 2e-12);
            }
            near(candidate.step().water_to_solid_kg, water, 1e-13);
            near(candidate.step().water_balance_residual_kg, 0.0, 1e-13);
            near(
                response.wet_mass_kg(),
                candidate.moisture().wet_mass_kg(),
                1e-13,
            );
            for path in response.expansion_paths() {
                near(
                    path.free_principal_strains()[2],
                    0.02 * (expected - 0.1),
                    1e-13,
                );
            }
            assert_eq!(initial.moisture(), &[m0; 8]);
            assert_eq!(reservoir.solution().solvent_volume_m3(), initial_volume);
        }
    });
}
