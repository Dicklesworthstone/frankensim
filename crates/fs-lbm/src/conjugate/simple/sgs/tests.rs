//! Manufactured gradients and the live momentum/heat coupling, without a
//! claim to a developed turbulent spectrum or a calibrated wall model.

use super::*;
use crate::conjugate::{
    EnergyConfig, FluidProperties, SimpleConfig, ThermalFace, ThermalSetup, TimeScheme,
    TransientConfig, UnsteadyConfig, Voxel, VoxelDomain, march_conjugate_with_sgs, march_energy,
    simple_unsteady, simple_unsteady_with_sgs,
};

fn fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 2.0,
        conductivity_w_m_k: 0.002,
        kinematic_viscosity_m2_s: 0.01,
    }
}

fn config(model: SgsModel) -> FvSgsConfig {
    FvSgsConfig {
        model,
        filter_width_factor: 1.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

fn fill_faces(solver: &mut Solver<'_>, field: impl Fn([f64; 3]) -> [f64; 3]) {
    for a in 0..3 {
        for f in 0..solver.vel[a].len() {
            let at = solver.comps[a].coords(f);
            let point = [0, 1, 2].map(|b| {
                let low = solver.domain.face_coord(b, at[b]);
                if a == b {
                    low
                } else {
                    low + 0.5 * solver.domain.width(b, at[b])
                }
            });
            solver.vel[a][f] = field(point)[a];
        }
    }
}

#[test]
fn nonuniform_derivatives_honor_value_and_neumann_at_the_actual_face() {
    // q(x)=3+2x+4x^2 at x=0, and unequal neighbour distances.
    let q = |x: f64| 3.0 + 2.0 * x + 4.0 * x * x;
    let (hm, hp) = (0.17, 0.41);
    let got = derivative(
        q(0.0),
        (hm, Sample::Value(q(-hm))),
        (hp, Sample::Value(q(hp))),
    );
    assert!((got - 2.0).abs() < 2e-14, "{got}");
    // q=(x+hm)^2 has zero derivative at the low boundary, not the centre.
    let got = derivative(
        hm * hm,
        (hm, Sample::ZeroGradient),
        (hp, Sample::Value((hm + hp).powi(2))),
    );
    assert!((got - 2.0 * hm).abs() < 2e-14, "{got}");
    let got = derivative(
        hp * hp,
        (hm, Sample::Value((hm + hp).powi(2))),
        (hp, Sample::ZeroGradient),
    );
    assert!((got + 2.0 * hp).abs() < 2e-14, "{got}");
}

#[test]
fn graded_couette_and_voxel_walls_feed_the_same_sgs_and_heat_viscosity() {
    let gate = CancelGate::new();
    for internal_wall in [false, true] {
        let wall_y = if internal_wall { 0.1 } else { 0.0 };
        let domain = VoxelDomain::graded_from_fn(
            [
                vec![0.2, 0.3, 0.5],
                vec![0.1, 0.12, 0.18, 0.22, 0.38],
                vec![0.4, 0.6],
            ],
            |p| {
                if p[1] < wall_y {
                    Voxel::Solid(0)
                } else {
                    Voxel::Fluid
                }
            },
        )
        .unwrap();
        let faces = [
            FvBoundary::Outlet,
            FvBoundary::Outlet,
            FvBoundary::wall(),
            FvBoundary::Wall {
                velocity: [1.0 - wall_y, 0.0, 0.0],
            },
            FvBoundary::Symmetry,
            FvBoundary::Symmetry,
        ];
        let controls = SimpleConfig::new(faces);
        let fluid = fluid();
        let mut solver = Solver::new(&domain, &fluid, &controls);
        fill_faces(&mut solver, |p| [(p[1] - wall_y).max(0.0), 0.0, 0.0]);
        solver.start_time(0.01, TimeScheme::Bdf2);
        solver
            .enable_sgs(config(SgsModel::Wale { cw: 0.5 }), &gate)
            .unwrap();
        solver.update_sgs(&gate).unwrap();
        assert!(
            solver.eddy_viscosity.iter().all(|&nu| nu == 0.0),
            "WALE pure shear"
        );
        solver
            .enable_sgs(config(SgsModel::Smagorinsky { cs: 0.2 }), &gate)
            .unwrap();
        solver.update_sgs(&gate).unwrap();
        let velocity: Vec<_> = (0..domain.cell_count())
            .map(|c| solver.cell_velocity(c))
            .collect();
        for c in 0..domain.cell_count() {
            if !domain.is_fluid(c) {
                assert_eq!(solver.eddy_viscosity[c], 0.0);
                continue;
            }
            let gradient = solver.sgs_gradient(c, &velocity);
            for (i, row) in gradient.iter().enumerate() {
                for (j, &value) in row.iter().enumerate() {
                    let exact = if i == 0 && j == 1 { 1.0 } else { 0.0 };
                    assert!(
                        (value - exact).abs() < 2e-14,
                        "cell={c}, G[{i}][{j}]={value}"
                    );
                }
            }
            let expected = 0.04 * fs_math::det::pow(domain.volume(c), 2.0 / 3.0);
            assert!((solver.eddy_viscosity[c] / expected - 1.0).abs() < 2e-13);
        }
        assert_eq!(solver.heat_eddy_viscosity(), solver.eddy_viscosity);
        let flow = solver.finish_projected(&fluid, 0, (0.0, 0.0));
        for (k_t, nu_t) in flow
            .eddy_conductivity(&fluid)
            .iter()
            .zip(&flow.eddy_viscosity_m2_s)
        {
            assert!((k_t - 2.0 * nu_t / 0.9).abs() < 1e-14);
        }
    }
}

#[test]
fn symmetry_and_mac_normal_gradients_are_not_clamped_at_adjacent_centres() {
    let domain = VoxelDomain::graded([
        vec![0.2, 0.3, 0.5],
        vec![0.1, 0.15, 0.3, 0.45],
        vec![0.4, 0.6],
    ])
    .unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    faces[3] = FvBoundary::Wall {
        velocity: [1.0, 0.0, 0.0],
    };
    let controls = SimpleConfig::new(faces);
    let mut solver = Solver::new(&domain, &fluid(), &controls);
    fill_faces(&mut solver, |p| [p[1] * p[1], 0.0, 0.0]);
    let velocities: Vec<_> = (0..domain.cell_count())
        .map(|c| solver.cell_velocity(c))
        .collect();
    for c in 0..domain.cell_count() {
        let [x, y, z] = domain.coords(c);
        let exact = 2.0 * domain.center(x, y, z)[1];
        assert!((solver.sgs_gradient(c, &velocities)[0][1] - exact).abs() < 2e-14);
    }
    // The native normal derivative uses the two owning MAC faces, exactly
    // reproducing affine extension even on differently stretched axes.
    fill_faces(&mut solver, |p| [2.0 * p[0], -3.0 * p[1], p[2]]);
    let velocities: Vec<_> = (0..domain.cell_count())
        .map(|c| solver.cell_velocity(c))
        .collect();
    for c in 0..domain.cell_count() {
        let gradient = solver.sgs_gradient(c, &velocities);
        for (a, exact) in [2.0, -3.0, 1.0].into_iter().enumerate() {
            assert!((gradient[a][a] - exact).abs() < 2e-14);
        }
    }
}

#[test]
fn sgs_refuses_competing_closures_caps_and_acceptance_use() {
    let domain = VoxelDomain::new(3, 3, 1, 0.1).unwrap();
    let mut controls = SimpleConfig::new([FvBoundary::wall(); 6]);
    let sgs = config(SgsModel::Smagorinsky { cs: 0.2 });
    let gate = CancelGate::new();
    let mut solver = Solver::new(&domain, &fluid(), &controls);
    assert!(
        solver.enable_sgs(sgs, &gate).is_err(),
        "steady SGS must refuse"
    );
    solver.start_time(0.01, TimeScheme::Bdf2);
    assert!(
        solver
            .enable_sgs(
                FvSgsConfig {
                    usage: SgsUse::AcceptanceRender,
                    ..sgs
                },
                &gate
            )
            .is_err()
    );
    solver
        .enable_sgs(
            FvSgsConfig {
                maximum_eddy_to_molecular: 0.0,
                ..sgs
            },
            &gate,
        )
        .unwrap();
    fill_faces(&mut solver, |p| [p[1], 0.0, 0.0]);
    assert!(
        solver.update_sgs(&gate).is_err(),
        "cap must refuse rather than clip"
    );
    controls.turbulence = Turbulence::Lvel;
    let mut solver = Solver::new(&domain, &fluid(), &controls);
    solver.start_time(0.01, TimeScheme::Bdf2);
    assert!(solver.enable_sgs(sgs, &gate).is_err());
    gate.request();
    assert!(matches!(
        solver.enable_sgs(sgs, &gate),
        Err(ChtError::Cancelled)
    ));
}

#[test]
fn transient_sgs_drives_heat_and_zero_coefficient_recovers_the_existing_march() {
    let gate = CancelGate::new();
    let domain = VoxelDomain::new(5, 5, 1, 0.1).unwrap();
    let fluid = fluid();
    let mut faces = [FvBoundary::wall(); 6];
    faces[3] = FvBoundary::Wall {
        velocity: [0.1, 0.0, 0.0],
    };
    faces[4] = FvBoundary::Symmetry;
    faces[5] = FvBoundary::Symmetry;
    let mut controls = SimpleConfig::new(faces);
    controls.pressure_tolerance = 1e-11;
    let mut unsteady = UnsteadyConfig::new(0.01, 1);
    unsteady.inner_iterations = 400;
    unsteady.inner_tolerance = 1e-9;
    let laminar = simple_unsteady(&domain, &fluid, &controls, &unsteady, |_| 1.0, &gate).unwrap();
    let disabled = simple_unsteady_with_sgs(
        &domain,
        &fluid,
        &controls,
        &unsteady,
        config(SgsModel::Smagorinsky { cs: 0.0 }),
        |_| 1.0,
        &gate,
    )
    .unwrap();
    assert_eq!(
        laminar, disabled,
        "disabled must recover all ordinary march arithmetic"
    );
    let sgs = config(SgsModel::Smagorinsky { cs: 0.2 });
    let initial: Vec<_> = (0..domain.cell_count())
        .map(|c| 300.0 + 10.0 * domain.coords(c)[1] as f64)
        .collect();
    let setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let coupled = march_conjugate_with_sgs(
        &domain,
        &fluid,
        &[],
        &setup,
        &controls,
        &unsteady,
        sgs,
        None,
        &initial,
        |_| 1.0,
        |_| 1.0,
        &EnergyConfig::default(),
        None,
        &gate,
    )
    .unwrap();
    assert!(
        coupled
            .flow
            .flow
            .eddy_viscosity_m2_s
            .iter()
            .any(|&nu| nu > 1e-8)
    );
    assert_ne!(
        coupled.flow.flow.field, laminar.flow.field,
        "SGS must enter momentum"
    );
    let mut explicit = setup.clone();
    explicit.eddy_conductivity_w_m_k = coupled.flow.flow.eddy_conductivity(&fluid);
    let thermal = TransientConfig {
        time_step_s: unsteady.time_step_s,
        steps: 1,
        energy: EnergyConfig::default(),
    };
    let frozen = march_energy(
        &domain,
        &fluid,
        &[],
        &coupled.flow.flow.field,
        &explicit,
        &initial,
        |_| 1.0,
        &thermal,
        None,
        &gate,
    )
    .unwrap();
    let molecular = march_energy(
        &domain,
        &fluid,
        &[],
        &coupled.flow.flow.field,
        &setup,
        &initial,
        |_| 1.0,
        &thermal,
        None,
        &gate,
    )
    .unwrap();
    let error = coupled
        .temperature
        .iter()
        .zip(&frozen.temperature)
        .fold(0.0_f64, |max, (a, b)| max.max((a - b).abs()));
    assert!(
        error < 1e-8,
        "coupled SGS heat must match explicit frozen conductivity: {error}"
    );
    let effect = coupled
        .temperature
        .iter()
        .zip(&molecular.temperature)
        .fold(0.0_f64, |max, (a, b)| max.max((a - b).abs()));
    assert!(
        effect > 1e-5,
        "SGS conductivity must affect heat transport: {effect}"
    );
    assert!(coupled.energy_records[0].closure_j.abs() < 1e-8);
    assert!(
        march_conjugate_with_sgs(
            &domain,
            &fluid,
            &[],
            &explicit,
            &controls,
            &unsteady,
            sgs,
            None,
            &initial,
            |_| 1.0,
            |_| 1.0,
            &EnergyConfig::default(),
            None,
            &gate,
        )
        .is_err(),
        "a second turbulent thermal closure must refuse"
    );
}
