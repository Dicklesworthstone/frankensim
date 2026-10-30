use super::*;
use crate::compressible::mixture::molecular::{MolecularTransport2d, PlanarTransport};
use crate::compressible::spatial::mixture::molecular::{
    DissipationControl, DissipationError, MolecularReport, MolecularTransport3d, MolecularWall3d,
    SpatialTransport, WallThermalCondition,
};
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, TransportCoefficients, TransportDomain,
};
use fs_qty::{DynViscosity, Temperature};
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

fn law(
    gas: &MixtureGas,
    diffusion: f64,
    shear: f64,
    bulk: f64,
    conductivity: f64,
) -> ConstantTransport {
    let n = gas.species_count();
    let mut binary = vec![diffusion; n * n];
    for i in 0..n {
        binary[i * n + i] = 0.0;
    }
    ConstantTransport::new(
        gas.source().clone(),
        "manufactured-volumetric-molecular",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(shear),
            bulk_viscosity: DynViscosity::new(bulk),
            thermal_conductivity: ThermalConductivityQuantityV1::new(conductivity),
        },
        binary,
    )
    .unwrap()
}
fn prepare(
    flow: &MixtureEuler3d,
    law: ConstantTransport,
) -> MolecularTransport3d<ConstantTransport> {
    let source = SpatialTransport::new(flow.gas().clone(), law, || false).unwrap();
    let bytes = MolecularTransport3d::required_workspace_bytes(flow, &source).unwrap();
    let result = MolecularTransport3d::new(flow, source, bytes, || false).unwrap();
    assert_eq!(result.workspace_bytes() + flow.workspace_bytes(), bytes);
    result
}
fn check(report: MolecularReport) {
    assert!(report.maximum_parabolic_number <= 0.40000000000001);
    assert!(report.maximum_mass_frame_residual_kg_m2_s < 1e-11);
    for value in report.conservation_residual.species_masses {
        close(value, 0.0, 3e-11);
    }
    for value in report.conservation_residual.momentum {
        close(value, 0.0, 3e-11);
    }
    close(report.conservation_residual.energy, 0.0, 3e-11);
    close(report.absolute_energy_residual_j, 0.0, 3e-9);
}
fn equal_mass_gas() -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.0, 40.0),
        source("B", R, 4.5, 0.0, -50.0),
        source("absent-C", R, 3.5, 0.0, 2.0),
    ])
}

fn diffusion_mode(n: usize, diffusion: f64) -> (f64, f64) {
    let h = 1.0 / n as f64;
    let z = std::f64::consts::PI * h;
    let average = (z.sin() / z).powi(3);
    let wave = |x: [f64; 3]| (2.0 * std::f64::consts::PI * (x[0] + x[1] + x[2])).sin();
    let grid = CartesianGrid::new(n, n, n, h, h, h).unwrap();
    let mut flow = field(equal_mass_gas(), grid, PERIODIC, |x| {
        let a = 0.5 + 0.2 * average * wave(x);
        ([a, 1.0 - a, 0.0], [0.0; 3], 1.0)
    });
    let before = flow.states().to_vec();
    let mut operator = prepare(&flow, law(flow.gas(), diffusion, 0.0, 0.0, 0.0));
    let duration = 0.05;
    check(
        operator
            .advance(
                &mut flow,
                Time::new(duration),
                DissipationControl::default(),
                || false,
            )
            .unwrap(),
    );
    assert_eq!(flow.time().value(), 0.0);
    assert!(flow.last_step().is_none());
    assert_ne!(flow.states(), before);
    let decay = (-diffusion * 12.0 * std::f64::consts::PI.powi(2) * duration).exp();
    let mut error = 0.0;
    let mut projection = 0.0;
    let mut norm = 0.0;
    for index in 0..grid.len() {
        let q = flow.cell(index).unwrap();
        let p = q.primitive(flow.gas(), || false).unwrap();
        let mode = wave(grid.cell_center(index).unwrap());
        error += (q.partial_densities[0] - (0.5 + 0.2 * average * decay * mode)).abs();
        projection += (q.partial_densities[0] - 0.5) * mode;
        norm += mode * mode;
        assert_eq!(q.partial_densities[2], 0.0);
        close(p.density, 1.0, 2e-12);
        // Distinct Cv and formation offsets: enthalpy transport changes E_ref
        // while the physical diffusion solution remains isothermal/isobaric.
        close(p.temperature, 1.0, 2e-10);
        close(p.pressure, 1.0, 2e-10);
    }
    (error / grid.len() as f64, projection / norm)
}

#[test]
fn g1_three_axis_species_diffusion_refines_with_original_partial_enthalpies() {
    let coarse = diffusion_mode(8, 0.02);
    let fine = diffusion_mode(16, 0.02);
    let order = (coarse.0 / fine.0).ln() / 2.0_f64.ln();
    eprintln!(
        "3D Maxwell-Stefan oblique L1 8^3/16^3: {} / {}; order={order}",
        coarse.0, fine.0
    );
    assert!(order > 1.8);
    assert!(fine.0 < 2e-4);
    let changed_data = diffusion_mode(16, 0.04);
    assert!(changed_data.1 < fine.1, "only original binary data changed");
}

#[test]
fn g1_complete_spatial_stress_matches_independent_vector_symbols_and_heating() {
    let n = 8;
    let h = 1.0 / n as f64;
    let mu = 0.015;
    let bulk = 0.005;
    let amplitude = 0.02;
    let duration: f64 = 0.02;
    let dt: f64 = 0.001;
    let wave = |x: [f64; 3]| (2.0 * std::f64::consts::PI * (x[0] + x[1] + x[2])).sin();
    let theta = 2.0 * std::f64::consts::PI * h;
    let beta = bulk + mu / 3.0;
    let diagonal = (3.0 * mu + beta) * 4.0 * (0.5 * theta).sin().powi(2) / (h * h);
    let cross = beta * (theta.sin() / h).powi(2);
    // Independent eigenvectors of the complete three-axis stencil: the
    // longitudinal mode sees two off-diagonal terms; the transverse mode
    // sees their difference. All three physical velocity components evolve.
    for (direction, eigenvalue) in [
        ([1.0, 1.0, 1.0], diagonal + 2.0 * cross),
        ([1.0, 2.0, -3.0], diagonal - cross),
    ] {
        let grid = CartesianGrid::new(n, n, n, h, h, h).unwrap();
        let mut flow = field(equal_mass_gas(), grid, PERIODIC, |x| {
            (
                [1.0, 0.0, 0.0],
                direction.map(|d| d * amplitude * wave(x)),
                1.0,
            )
        });
        let kinetic = |flow: &MixtureEuler3d| {
            (0..flow.cell_count())
                .map(|i| {
                    let q = flow.cell(i).unwrap();
                    0.5 * q.momentum.iter().map(|p| p * p).sum::<f64>()
                })
                .sum::<f64>()
                * flow.grid.cell_volume()
        };
        let before_kinetic = kinetic(&flow);
        let before_energy = flow.totals().energy;
        let mut operator = prepare(&flow, law(flow.gas(), 0.001, mu, bulk, 0.0));
        let report = operator
            .advance(
                &mut flow,
                Time::new(duration),
                DissipationControl {
                    maximum_step_s: dt,
                    ..DissipationControl::default()
                },
                || false,
            )
            .unwrap();
        check(report);
        let decay = (1.0 - eigenvalue * dt + 0.5 * (eigenvalue * dt).powi(2))
            .powi((duration / dt).round() as i32);
        for index in 0..grid.len() {
            let p = flow
                .cell(index)
                .unwrap()
                .primitive(flow.gas(), || false)
                .unwrap();
            for axis in 0..3 {
                let expected =
                    direction[axis] * amplitude * wave(grid.cell_center(index).unwrap()) * decay;
                close(p.velocity[axis], expected, 3e-12);
            }
        }
        close(flow.totals().energy, before_energy, 2e-12);
        let after_kinetic = kinetic(&flow);
        assert!(after_kinetic < before_kinetic);
        let thermal_gain = (0..grid.len())
            .map(|i| {
                let p = flow
                    .cell(i)
                    .unwrap()
                    .primitive(flow.gas(), || false)
                    .unwrap();
                2.5 * (p.temperature - 1.0) * grid.cell_volume()
            })
            .sum::<f64>();
        close(thermal_gain, before_kinetic - after_kinetic, 3e-12);
    }
}

#[test]
fn g1_spatial_fourier_heat_matches_independent_conservative_diffusion() {
    let n = 8;
    let h = 1.0 / n as f64;
    let conductivity = 0.04;
    let duration: f64 = 0.02;
    let dt: f64 = 0.001;
    let z = std::f64::consts::PI * h;
    let average = (z.sin() / z).powi(3);
    let mode = |x: [f64; 3]| (2.0 * std::f64::consts::PI * (x[0] + x[1] + x[2])).cos();
    let grid = CartesianGrid::new(n, n, n, h, h, h).unwrap();
    let mut flow = field(equal_mass_gas(), grid, PERIODIC, |x| {
        ([1.0, 0.0, 0.0], [0.0; 3], 1.0 + 0.1 * average * mode(x))
    });
    let before = flow.totals();
    let mut operator = prepare(&flow, law(flow.gas(), 0.001, 0.0, 0.0, conductivity));
    check(
        operator
            .advance(
                &mut flow,
                Time::new(duration),
                DissipationControl {
                    maximum_step_s: dt,
                    ..DissipationControl::default()
                },
                || false,
            )
            .unwrap(),
    );
    let eigenvalue = conductivity / 2.5 * 12.0 * z.sin().powi(2) / (h * h);
    let decay = (1.0 - eigenvalue * dt + 0.5 * (eigenvalue * dt).powi(2))
        .powi((duration / dt).round() as i32);
    for i in 0..grid.len() {
        let q = flow.cell(i).unwrap();
        assert_eq!(q.partial_densities, [1.0, 0.0, 0.0]);
        assert_eq!(q.momentum, [0.0; 3]);
        let p = q.primitive(flow.gas(), || false).unwrap();
        close(
            p.temperature,
            1.0 + 0.1 * average * decay * mode(grid.cell_center(i).unwrap()),
            3e-12,
        );
    }
    close(flow.totals().energy, before.energy, 2e-12);
}

#[test]
fn g3_spatial_molecular_matches_planar_extrusion_and_three_axis_rotation() {
    let g = standard(0.0);
    let initial = |x: [f64; 3]| {
        let angle = 2.0 * std::f64::consts::PI * (x[0] + x[1]);
        let a = 0.4 + 0.1 * angle.sin();
        let b = 1.0 - a;
        let temperature = 1.0 + 0.1 * angle.cos();
        (
            [a, b, 0.0],
            [0.07 * angle.sin(), 0.05 * angle.cos(), 0.23],
            (a + 0.5 * b) * temperature,
        )
    };
    let make_spatial = |axes: [usize; 3]| {
        let mut shape = [0; 3];
        let mut spacing = [0.0; 3];
        for d in 0..3 {
            shape[axes[d]] = [6, 4, 3][d];
            spacing[axes[d]] = [1.0 / 6.0, 0.25, 0.2][d];
        }
        let grid = CartesianGrid::new(
            shape[0], shape[1], shape[2], spacing[0], spacing[1], spacing[2],
        )
        .unwrap();
        field(g.clone(), grid, PERIODIC, |x| {
            let (rho, u, p) = initial(axes.map(|axis| x[axis]));
            let mut velocity = [0.0; 3];
            for d in 0..3 {
                velocity[axes[d]] = u[d];
            }
            (rho, velocity, p)
        })
    };
    let mut horizontal = make_spatial([0, 1, 2]);
    let grid = Grid2d::new(6, 4, 1.0 / 6.0, 0.25, 0.6).unwrap();
    let mut states = vec![0.0; 6 * 4 * 6];
    for (index, row) in states.chunks_exact_mut(6).enumerate() {
        let x = grid.cell_center(index).unwrap();
        let (rho, u, p) = initial([x[0], x[1], 0.0]);
        planar::state_into(&g, &rho, [u[0], u[1]], p, row).unwrap();
    }
    let bytes = planar::MixtureEuler2d::required_workspace_bytes(&g, &grid).unwrap();
    let mut planar = planar::MixtureEuler2d::new(
        g.clone(),
        grid,
        states,
        Boundaries2d {
            x: AxisBoundaries::Periodic,
            y: AxisBoundaries::Periodic,
        },
        bytes,
        || false,
    )
    .unwrap();
    let source =
        PlanarTransport::new(g.clone(), law(&g, 0.03, 0.02, 0.01, 0.03), || false).unwrap();
    let bytes = MolecularTransport2d::required_workspace_bytes(&planar, &source).unwrap();
    let mut planar_op = MolecularTransport2d::new(&planar, source, bytes, || false).unwrap();
    let controls = DissipationControl {
        maximum_step_s: 0.00005,
        ..DissipationControl::default()
    };
    let mut horizontal_op = prepare(&horizontal, law(&g, 0.03, 0.02, 0.01, 0.03));
    check(
        horizontal_op
            .advance(&mut horizontal, Time::new(0.0002), controls, || false)
            .unwrap(),
    );
    planar_op
        .advance(&mut planar, Time::new(0.0002), controls, || false)
        .unwrap();
    for index in 0..horizontal.cell_count() {
        let spatial = horizontal.cell(index).unwrap();
        let planar = planar.cell(index % 24).unwrap();
        for d in 0..3 {
            close(
                spatial.partial_densities[d],
                planar.partial_densities[d],
                3e-11,
            );
        }
        for d in 0..2 {
            close(spatial.momentum[d], planar.momentum[d], 3e-11);
        }
        close(spatial.momentum[2], 0.23, 3e-11);
        close(
            spatial.energy - 0.5 * 0.23_f64.powi(2),
            planar.energy,
            3e-11,
        );
    }
    for axes in [[2, 0, 1], [1, 2, 0]] {
        let mut rotated = make_spatial(axes);
        let mut operator = prepare(&rotated, law(&g, 0.03, 0.02, 0.01, 0.03));
        check(
            operator
                .advance(&mut rotated, Time::new(0.0002), controls, || false)
                .unwrap(),
        );
        for index in 0..horizontal.cell_count() {
            let coordinates = horizontal.grid.coordinates(index);
            let mut other = [0; 3];
            for d in 0..3 {
                other[axes[d]] = coordinates[d];
            }
            let expected = horizontal.cell(index).unwrap();
            let actual = rotated.cell(rotated.grid.index(other)).unwrap();
            for d in 0..3 {
                close(
                    actual.partial_densities[d],
                    expected.partial_densities[d],
                    3e-11,
                );
                close(actual.momentum[axes[d]], expected.momentum[d], 3e-11);
            }
            close(actual.energy, expected.energy, 3e-11);
        }
    }
}

fn wall_field() -> MixtureEuler3d {
    field(
        standard(0.0),
        CartesianGrid::new(4, 3, 2, 0.25, 1.0 / 3.0, 0.3).unwrap(),
        Boundaries3d {
            x: WALLS,
            y: WALLS,
            z: WALLS,
        },
        |x| {
            (
                [0.6, 0.4, 0.0],
                [0.15 + 0.05 * x[0], -0.2 + 0.1 * x[1], 0.1 + 0.03 * x[2]],
                0.8 + 0.1 * x[0] + 0.05 * x[2],
            )
        },
    )
}
fn impulses(operator: &MolecularTransport3d<ConstantTransport>) -> Vec<([f64; 3], [f64; 3], f64)> {
    BoundarySide::ALL
        .into_iter()
        .flat_map(|side| operator.boundary_impulses(side))
        .map(|q| (q.species_masses.try_into().unwrap(), q.momentum, q.energy))
        .collect()
}

#[test]
fn g4_spatial_molecular_walls_close_and_requested_duration_is_atomic() {
    let create = || {
        let flow = wall_field();
        let op = prepare(&flow, law(flow.gas(), 0.02, 0.03, 0.01, 0.02));
        (flow, op)
    };
    let (mut reference, mut expected_op) = create();
    expected_op
        .advance(
            &mut reference,
            Time::new(0.001),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let mut polls = 0;
    let expected = expected_op
        .advance(
            &mut reference,
            Time::new(0.002),
            DissipationControl::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    check(expected);
    let mut outward = [0.0; 3];
    for side in BoundarySide::ALL {
        for q in expected_op.boundary_impulses(side) {
            assert_eq!(q.species_masses, [0.0; 3]);
            assert_eq!(q.energy, 0.0);
            for axis in 0..3 {
                outward[axis] += q.momentum[axis];
                if axis != side.axis() {
                    assert_eq!(q.momentum[axis], 0.0);
                }
            }
        }
    }
    for axis in 0..3 {
        close(
            expected.after.momentum[axis] - expected.before.momentum[axis] + outward[axis],
            0.0,
            2e-12,
        );
    }
    assert!(outward.iter().all(|value| value.abs() > 1e-6));
    let (mut flow, mut operator) = create();
    operator
        .advance(
            &mut flow,
            Time::new(0.001),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let before = flow.states().to_vec();
    let before_impulses = impulses(&operator);
    let before_report = *operator.last_step().unwrap();
    let mut attempted = 0;
    assert_eq!(
        operator
            .advance(
                &mut flow,
                Time::new(0.002),
                DissipationControl::default(),
                || {
                    attempted += 1;
                    attempted >= polls
                },
            )
            .unwrap_err(),
        DissipationError::Cancelled,
    );
    assert_eq!(flow.states(), before);
    assert_eq!(flow.time().value(), 0.0);
    assert!(flow.last_step().is_none());
    assert_eq!(operator.last_step(), Some(&before_report));
    assert_eq!(impulses(&operator), before_impulses);
    assert!(matches!(
        operator.advance(
            &mut flow,
            Time::new(0.002),
            DissipationControl {
                maximum_step_s: 0.0001,
                maximum_steps: 1,
                ..DissipationControl::default()
            },
            || false,
        ),
        Err(DissipationError::WorkLimit { .. })
    ));
    assert_eq!(flow.states(), before);
    assert_eq!(impulses(&operator), before_impulses);
    assert_eq!(
        operator
            .advance(
                &mut flow,
                Time::new(0.002),
                DissipationControl::default(),
                || false
            )
            .unwrap(),
        expected
    );
    assert_eq!(flow.states(), reference.states());
    assert_eq!(impulses(&operator), impulses(&expected_op));
}

#[test]
fn g0_spatial_molecular_budget_capacity_and_boundary_ownership_are_explicit() {
    let mut flow = wall_field();
    let mut states = Vec::with_capacity(flow.states.len() + 37);
    states.extend_from_slice(&flow.states);
    let extra = (states.capacity() - flow.states.capacity()) * size_of::<f64>();
    flow.states = states;
    flow.workspace_bytes += extra;
    let source = SpatialTransport::new(
        flow.gas().clone(),
        law(flow.gas(), 0.01, 0.02, 0.0, 0.03),
        || false,
    )
    .unwrap();
    let required = MolecularTransport3d::required_workspace_bytes(&flow, &source).unwrap();
    assert!(matches!(
        MolecularTransport3d::new(&flow, source, required - 1, || false),
        Err(DissipationError::Invalid(_))
    ));
    let mut op = prepare(&flow, law(flow.gas(), 0.01, 0.02, 0.0, 0.03));
    let combined = flow.workspace_bytes() + op.workspace_bytes();
    for _ in 0..3 {
        op.advance(
            &mut flow,
            Time::new(0.001),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
        assert_eq!(flow.workspace_bytes() + op.workspace_bytes(), combined);
    }
    let before = flow.states().to_vec();
    let zero = op
        .advance(
            &mut flow,
            Time::new(0.0),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(zero.elapsed.value(), 0.0);
    assert_eq!(zero.accepted_steps, 0);
    assert_eq!(flow.states(), before);
    assert!(
        impulses(&op)
            .iter()
            .all(|(s, p, e)| *s == [0.0; 3] && *p == [0.0; 3] && *e == 0.0)
    );
    flow.boundaries.z = OPEN;
    assert!(matches!(
        op.advance(
            &mut flow,
            Time::new(0.001),
            DissipationControl::default(),
            || false
        ),
        Err(DissipationError::Invalid(_))
    ));
    assert_eq!(flow.states(), before);
    let source = SpatialTransport::new(
        flow.gas().clone(),
        law(flow.gas(), 0.01, 0.02, 0.0, 0.03),
        || false,
    )
    .unwrap();
    assert!(matches!(
        MolecularTransport3d::new(&flow, source, MAX_MIXTURE_WORKSPACE_BYTES, || false),
        Err(DissipationError::Invalid(_))
    ));
}


#[test]
fn g1_no_slip_isothermal_wall_transfers_shear_and_heat_conservatively() {
    let grid = CartesianGrid::new(3, 3, 3, 0.2, 0.2, 0.2).unwrap();
    let mut flow = field(equal_mass_gas(), grid, Boundaries3d {
        x: WALLS,
        y: WALLS,
        z: WALLS,
    }, |_| ([1.0, 0.0, 0.0], [0.0, 0.2, -0.1], 1.0));
    let source = SpatialTransport::new(
        flow.gas().clone(),
        law(flow.gas(), 0.0, 0.03, 0.01, 0.04),
        || false,
    )
    .unwrap();
    let bytes = MolecularTransport3d::required_workspace_bytes(&flow, &source).unwrap();
    let mut operator = MolecularTransport3d::new(&flow, source, bytes, || false)
        .unwrap()
        .with_wall(
            BoundarySide::XMin,
            MolecularWall3d::NoSlip {
                velocity: [0.0; 3],
                thermal: WallThermalCondition::Isothermal(Temperature::new(0.8)),
            },
        )
        .unwrap();
    let report = operator
        .advance(
            &mut flow,
            Time::new(0.0002),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    check(report);
    let xmin: Vec<_> = operator.boundary_impulses(BoundarySide::XMin).collect();
    assert!(xmin.iter().all(|q| q.species_masses == [0.0; 3]));
    assert!(xmin.iter().any(|q| q.momentum[1].abs() > 0.0));
    assert!(xmin.iter().any(|q| q.momentum[2].abs() > 0.0));
    assert!(xmin.iter().all(|q| q.energy > 0.0));
    let outward_energy: f64 = BoundarySide::ALL
        .into_iter()
        .flat_map(|side| operator.boundary_impulses(side))
        .map(|q| q.energy)
        .sum();
    close(
        report.after.energy - report.before.energy + outward_energy,
        0.0,
        3e-11,
    );
    assert!(MolecularTransport3d::new(
        &flow,
        SpatialTransport::new(
            flow.gas().clone(),
            law(flow.gas(), 0.0, 0.03, 0.01, 0.04),
            || false,
        )
        .unwrap(),
        bytes,
        || false,
    )
    .unwrap()
    .with_wall(
        BoundarySide::XMin,
        MolecularWall3d::NoSlip {
            velocity: [0.1, 0.0, 0.0],
            thermal: WallThermalCondition::Adiabatic,
        },
    )
    .is_err());
}
