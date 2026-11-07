//! Actual wall reflection, nonlinear source viscosity and rotational SGS agree.

use super::*;
use fs_lbm::d3q19::freesurface3::{FreeSurface3, PlanarWall3, PlanarWallFieldControls3};
use fs_lbm::material::free_surface::MaterialPlanarWallVelocityField3;

fn wall_controls() -> PlanarWallFieldControls3 {
    PlanarWallFieldControls3 {
        maximum_domain_cells: 512,
        maximum_wall_faces: 64,
        maximum_relative_mass_defect: 1e-11,
        maximum_relative_face_mass_defect: 1e-12,
    }
}

fn wall(velocities: &[[f64; 2]]) -> MaterialPlanarWallVelocityField3<'_> {
    MaterialPlanarWallVelocityField3 {
        wall: PlanarWall3::ZMin,
        moment_origin_m: [-0.003, 0.009, 0.002],
        tangential_velocity_m_s: velocities,
    }
}

fn velocities() -> Vec<[f64; 2]> {
    (0..64)
        .map(|i| {
            [
                0.02 + 0.002 * ((i % 8) as f64 - 3.5),
                -0.001 * ((i / 8) as f64 - 3.5),
            ]
        })
        .collect()
}

fn liquid<'s, 'm>(
    law: &'s IsothermalRheologyLaw<'m>,
    cells: &[MaterialFreeSurfaceInitialCell3],
    gravity: [f64; 3],
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> MaterialFreeSurfaceRheology3<'s, 'm> {
    MaterialFreeSurfaceRheology3::from_cells_with_pressure_reference(
        law,
        DX,
        DT,
        0.1,
        domain(gravity),
        cells,
        solve(),
        100000.0,
        cx,
        gate,
    )
    .unwrap()
}

// Independent population reduction, including the actual accepted Guo frame.
// These fixtures populate every tile below their surface in Morton order.
fn raw_stress(flow: &FreeSurface3, c: [usize; 3]) -> f64 {
    let [x, y, z] = c;
    let bits = flow.grid().state_bits();
    let slot = (x / 4) | ((y / 4) << 1) | ((z / 4) << 2);
    let lane = x % 4 + 4 * (y % 4 + 4 * (z % 4));
    let sample = flow.sample_at(x, y, z).unwrap();
    let force = flow.force_density_at(x, y, z).unwrap();
    let eq = equilibrium3(sample.density, sample.velocity);
    let mut p = [[0.0; 3]; 3];
    for q in 0..Q3 {
        let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
        let neq = f64::from_bits(bits[(slot * Q3 + q) * 64 + lane]) - eq[q];
        for a in 0..3 {
            for b in 0..3 {
                p[a][b] += e[a] * e[b] * neq;
            }
        }
    }
    for a in 0..3 {
        for b in 0..3 {
            p[a][b] += 0.5 * (sample.velocity[a] * force[b] + sample.velocity[b] * force[a]);
        }
    }
    let mean = (p[0][0] + p[1][1] + p[2][2]) / 3.0;
    for a in 0..3 {
        p[a][a] -= mean;
    }
    3.0 * (0.5 * p.iter().flatten().map(|v| v * v).sum::<f64>()).sqrt() / sample.density
}

// The first liquid center is half a cell from the actual wall, not one cell.
// Tangential coordinates are interior here; dx/dt=1 in this SI fixture.
fn wall_gradient(flow: &FreeSurface3, velocities: &[[f64; 2]]) -> [[f64; 3]; 3] {
    let u = |x, y, z| flow.sample_at(x, y, z).unwrap().velocity;
    let center = u(3, 3, 0);
    let next = u(3, 3, 1);
    let wall = [velocities[3 * 8 + 3][0], velocities[3 * 8 + 3][1], 0.0];
    core::array::from_fn(|a| {
        [
            (u(4, 3, 0)[a] - u(2, 3, 0)[a]) / 2.0,
            (u(3, 4, 0)[a] - u(3, 2, 0)[a]) / 2.0,
            -4.0 / 3.0 * wall[a] + center[a] + next[a] / 3.0,
        ]
    })
}

#[test]
fn quadratic_rotating_two_layer_film_binds_the_actual_wall_without_advancing() {
    let gate = CancelGate::new_clock_free();
    let source = card(0.003, 1000.0, 0.001);
    let law = law(&source);
    let walls: Vec<_> = (0..64)
        .map(|i| {
            [
                0.001 * ((i / 8) as f64 - 3.5),
                -0.002 * ((i % 8) as f64 - 3.5),
            ]
        })
        .collect();
    let cells: Vec<_> = initial(true)
        .into_iter()
        .filter(|cell| cell.coordinates[2] < 2)
        .map(|mut cell| {
            let [x, y, z] = cell.coordinates.map(|v| v as f64);
            let s = z + 0.5;
            cell.velocity_m_s = [
                0.001 * (y - 3.5) + 0.0006 * s * s,
                -0.002 * (x - 3.5) + 0.0003 * s * s,
                0.0,
            ];
            if z == 1.0 {
                cell.fill_fraction = 0.6;
            }
            cell
        })
        .collect();
    with_cx(&gate, |cx| {
        let original = liquid(&law, &cells, [0.0; 3], cx, &gate);
        for model in [SgsModel::Wale { cw: 0.6 }, SgsModel::Vreman { cv: 0.2 }] {
            let mut kinematics = cache();
            let mut workspace = FreeSurfaceSgsWorkspace3::default();
            assert!(
                original
                    .clone()
                    .with_sgs_si(
                        controls(model),
                        &gate,
                        Some(&mut kinematics),
                        &mut workspace,
                    )
                    .is_err(),
                "the liquid-only two-layer stencil must still refuse"
            );
            let bound = original
                .clone()
                .with_sgs_and_planar_wall_velocity_field_si(
                    controls(model),
                    wall(&walls),
                    wall_controls(),
                    &gate,
                    Some(&mut kinematics),
                    &mut workspace,
                )
                .unwrap();
            assert_eq!(
                bound.lattice().grid().state_bits(),
                original.lattice().grid().state_bits()
            );
            assert_eq!(
                bound.ledger_mass_kg().unwrap(),
                original.ledger_mass_kg().unwrap()
            );
            for z in 0..2 {
                let gradient = [
                    [0.0, 0.001, 0.0012 * (z as f64 + 0.5)],
                    [-0.002, 0.0, 0.0006 * (z as f64 + 0.5)],
                    [0.0; 3],
                ];
                let eddy = eddy_viscosity(model, 1.7, gradient).unwrap();
                assert!(eddy > 0.0);
                let c = [3, 3, z];
                near(
                    workspace.relaxation_at(c).unwrap(),
                    bound.rheology_at(c).unwrap().collision_tau + 3.0 * eddy,
                    3e-13,
                );
                assert_eq!(
                    bound.material_at(c).unwrap().card_identity(),
                    source.content_hash()
                );
            }
        }
    });
}

#[test]
fn changing_wall_recloses_the_source_and_bulk_before_collision_and_at_the_endpoint() {
    let gate = CancelGate::new_clock_free();
    let source = bulk_card(0.003);
    let law = bulk_law(&source);
    let walls = velocities();
    let resting = [[0.0; 2]; 64];
    let c = [3, 3, 0];
    with_cx(&gate, |cx| {
        let input = primed(&law, cx, &gate);
        for model in [
            SgsModel::Smagorinsky { cs: 0.6 },
            SgsModel::Wale { cw: 0.6 },
            SgsModel::Vreman { cv: 0.2 },
        ] {
            let mut ka = cache();
            let mut kb = cache();
            let mut wa = FreeSurfaceSgsWorkspace3::default();
            let mut wb = FreeSurfaceSgsWorkspace3::default();
            let mut a = input
                .clone()
                .with_sgs_and_planar_wall_velocity_field_si(
                    controls(model),
                    wall(&resting),
                    wall_controls(),
                    &gate,
                    Some(&mut ka),
                    &mut wa,
                )
                .unwrap();
            let mut b = input
                .clone()
                .with_sgs_and_planar_wall_velocity_field_si(
                    controls(model),
                    wall(&walls),
                    wall_controls(),
                    &gate,
                    Some(&mut kb),
                    &mut wb,
                )
                .unwrap();
            if !matches!(model, SgsModel::Smagorinsky { .. }) {
                assert!(
                    (a.rheology_at(c).unwrap().shear_rate_per_s
                        - b.rheology_at(c).unwrap().shear_rate_per_s)
                        .abs()
                        > 1e-7
                );
                assert_ne!(
                    a.material_at(c)
                        .unwrap()
                        .property("bulk_eta")
                        .unwrap()
                        .value_si(),
                    b.material_at(c)
                        .unwrap()
                        .property("bulk_eta")
                        .unwrap()
                        .value_si()
                );
            }
            // Both have identical populations but different previous source/bulk
            // epochs. Applying the SAME current wall must remove that history.
            let mass = a.ledger_mass_kg().unwrap();
            let (exchange, report) = a
                .step_with_planar_wall_velocity_field_si(
                    wall(&walls),
                    wall_controls(),
                    &gate,
                    Some(&mut ka),
                    &mut wa,
                )
                .unwrap();
            b.step_with_planar_wall_velocity_field_si(
                wall(&walls),
                wall_controls(),
                &gate,
                Some(&mut kb),
                &mut wb,
            )
            .unwrap();
            for (left, right) in a
                .lattice()
                .grid()
                .state_bits()
                .iter()
                .zip(b.lattice().grid().state_bits())
            {
                near(f64::from_bits(*left), f64::from_bits(right), 2e-12);
            }
            assert!(report.sgs_dissipation_w > 0.0);
            assert!(report.resolved_dissipation_w > 0.0);
            near(a.ledger_mass_kg().unwrap(), mass, 1e-10 * mass);
            let work: f64 = exchange
                .faces
                .iter()
                .zip(&walls)
                .map(|(face, velocity)| {
                    let impulse = face.exchange.fluid_impulse_kg_m_s;
                    impulse[0] * velocity[0] + impulse[1] * velocity[1]
                })
                .sum();
            near(exchange.exchange.fluid_work_j, work, 2e-18);
            let raw = raw_stress(a.lattice(), c);
            let (alpha, beta) = if let SgsModel::Smagorinsky { cs } = model {
                (0.8, 0.009 + 3.0 * (cs * 1.7).powi(2) * DT)
            } else {
                (
                    0.8 + 3.0
                        * eddy_viscosity(model, 1.7, wall_gradient(a.lattice(), &walls)).unwrap(),
                    0.009,
                )
            };
            let gamma = 2.0 * raw / DT / (alpha + (alpha * alpha + 4.0 * beta * raw / DT).sqrt());
            let sample = a.rheology_at(c).unwrap();
            near(sample.shear_rate_per_s, gamma, 3e-9);
            near(sample.collision_tau, 0.8 + 0.009 * gamma, 3e-11);
            near(wa.relaxation_at(c).unwrap(), alpha + beta * gamma, 3e-11);
            near(
                a.material_at(c)
                    .unwrap()
                    .property("bulk_eta")
                    .unwrap()
                    .value_si(),
                0.2 + 0.004 * gamma,
                3e-11,
            );
            assert_eq!(
                a.material_at(c).unwrap().card_identity(),
                source.content_hash()
            );
        }
    });
}

#[test]
fn constant_source_wall_sgs_matches_the_original_native_gradient_collision_and_work() {
    let gate = CancelGate::new_clock_free();
    let source = card(0.0, 1000.0, 0.001);
    let law = law(&source);
    let frozen = ResolvedLiquidSurfaceState::try_from_material(
        bundle(&source),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
    )
    .unwrap();
    let walls = velocities();
    let cells = initial(true);
    let gravity = [0.03, -0.02, 0.01];
    with_cx(&gate, |cx| {
        for model in [SgsModel::Wale { cw: 0.6 }, SgsModel::Vreman { cv: 0.2 }] {
            let mut ka = cache();
            let mut kb = cache();
            let mut wa = FreeSurfaceSgsWorkspace3::default();
            let mut wb = FreeSurfaceSgsWorkspace3::default();
            let mut actual = liquid(&law, &cells, gravity, cx, &gate)
                .with_sgs_and_planar_wall_velocity_field_si(
                    controls(model),
                    wall(&walls),
                    wall_controls(),
                    &gate,
                    Some(&mut ka),
                    &mut wa,
                )
                .unwrap();
            let mut expected = plan_free_surface_from_material(&frozen, DX, DT, 0.1)
                .unwrap()
                .with_pressure_reference(100000.0)
                .unwrap()
                .free_surface3_from_cells(domain(gravity), &cells, cx)
                .unwrap();
            for _ in 0..2 {
                let (wall_a, _) = actual
                    .step_with_planar_wall_velocity_field_si(
                        wall(&walls),
                        wall_controls(),
                        &gate,
                        Some(&mut ka),
                        &mut wa,
                    )
                    .unwrap();
                let (wall_b, _) = expected
                    .step_with_planar_wall_velocity_field_and_gradient_sgs_si(
                        wall(&walls),
                        wall_controls(),
                        &gate,
                        controls(model),
                        &mut kb,
                        &mut wb,
                    )
                    .unwrap();
                assert_eq!(wall_a, wall_b);
                assert_eq!(
                    actual.lattice().grid().state_bits(),
                    expected.lattice().grid().state_bits()
                );
                assert_eq!(
                    actual.ledger_mass_kg().unwrap(),
                    expected.ledger_mass_kg().unwrap()
                );
            }
        }
    });
}

#[test]
fn zero_sgs_preserves_molecular_wall_receipts_and_late_source_refusal_is_atomic() {
    let gate = CancelGate::new_clock_free();
    let source = bulk_card(0.003);
    let law = bulk_law(&source);
    let walls = velocities();
    with_cx(&gate, |cx| {
        let input = primed(&law, cx, &gate);
        for model in [
            SgsModel::Smagorinsky { cs: 0.0 },
            SgsModel::Wale { cw: 0.0 },
            SgsModel::Vreman { cv: 0.0 },
        ] {
            let mut ordinary = input.clone();
            let mut scratch = FreeSurfaceSgsWorkspace3::default();
            let mut zero = input
                .clone()
                .with_sgs_si(controls(model), &gate, None, &mut scratch)
                .unwrap();
            for _ in 0..2 {
                let expected = ordinary
                    .step_with_planar_wall_velocity_field_si(wall(&walls), wall_controls(), &gate)
                    .unwrap();
                let (actual, report) = zero
                    .step_with_planar_wall_velocity_field_si(
                        wall(&walls),
                        wall_controls(),
                        &gate,
                        None,
                        &mut scratch,
                    )
                    .unwrap();
                assert_eq!(actual, expected);
                assert_eq!(report.sgs_dissipation_w, 0.0);
                assert_eq!(
                    zero.lattice().grid().state_bits(),
                    ordinary.lattice().grid().state_bits()
                );
                assert_eq!(
                    zero.ledger_mass_kg().unwrap(),
                    ordinary.ledger_mass_kg().unwrap()
                );
                for c in [[3, 3, 0], CELL] {
                    assert_eq!(zero.rheology_at(c), ordinary.rheology_at(c));
                    assert_eq!(zero.material_at(c), ordinary.material_at(c));
                    assert_eq!(
                        scratch.relaxation_at(c).unwrap(),
                        ordinary.rheology_at(c).unwrap().collision_tau
                    );
                }
            }
        }
        let narrow = card(0.003, 1e-6, 0.001);
        let binding = super::super::law(&narrow);
        let cells: Vec<_> = initial(true)
            .into_iter()
            .map(|mut c| {
                c.velocity_m_s = [0.0; 3];
                c
            })
            .collect();
        let mut scratch = FreeSurfaceSgsWorkspace3::default();
        let mut closed = liquid(&binding, &cells, [0.0; 3], cx, &gate)
            .with_sgs_si(
                controls(SgsModel::Smagorinsky { cs: 0.2 }),
                &gate,
                None,
                &mut scratch,
            )
            .unwrap();
        let accepted = closed.clone();
        assert!(matches!(
            closed.step_with_planar_wall_velocity_field_si(
                wall(&walls),
                wall_controls(),
                &gate,
                None,
                &mut scratch,
            ),
            Err(FreeSurfaceRheologyError::Material(_))
        ));
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            closed
                .step_with_planar_wall_velocity_field_si(
                    wall(&walls),
                    wall_controls(),
                    &cancelled,
                    None,
                    &mut scratch,
                )
                .is_err()
        );
        assert_eq!(
            closed.lattice().grid().state_bits(),
            accepted.lattice().grid().state_bits()
        );
        assert_eq!(
            closed.ledger_mass_kg().unwrap(),
            accepted.ledger_mass_kg().unwrap()
        );
        assert_eq!(
            closed.material_at([3, 3, 0]),
            accepted.material_at([3, 3, 0])
        );
        assert_eq!(
            closed.rheology_at([3, 3, 0]),
            accepted.rheology_at([3, 3, 0])
        );
        let resting = [[0.0; 2]; 64];
        let mut direct = accepted;
        let retry = closed
            .step_with_planar_wall_velocity_field_si(
                wall(&resting),
                wall_controls(),
                &gate,
                None,
                &mut scratch,
            )
            .unwrap();
        let straight = direct
            .step_with_planar_wall_velocity_field_si(
                wall(&resting),
                wall_controls(),
                &gate,
                None,
                &mut FreeSurfaceSgsWorkspace3::default(),
            )
            .unwrap();
        assert_eq!(retry, straight);
        assert_eq!(
            closed.lattice().grid().state_bits(),
            direct.lattice().grid().state_bits()
        );
    });
}
