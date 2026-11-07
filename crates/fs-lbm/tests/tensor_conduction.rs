//! Native full-tensor face-temperature conduction: physical patch and owner tests.

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    CellSink, ChtError, ConductivityTensor, ContactResistance, EnergyConfig,
    FaceContactResistance, FlowField, FluidProperties, SolidMaterial, ThermalFace,
    ThermalSetup, Voxel, VoxelDomain, solve_energy, solve_tensor_conduction,
};

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(x,y)| x*y).sum()
}

fn center(domain: &VoxelDomain, cell: usize) -> [f64; 3] {
    let [x,y,z] = domain.coords(cell);
    domain.center(x,y,z)
}

#[test]
fn rotated_affine_patch_on_graded_mesh_retains_heat_with_equal_neighbor_temperatures() {
    let k = [[5.0,2.0,-1.0],[2.0,4.0,0.5],[-1.0,0.5,3.0]];
    let gradient = [0.0,7.0,-2.0];
    let flux = k.map(|row| -dot(row, gradient));
    let domain = VoxelDomain::graded_from_fn(
        [vec![0.2,0.35,0.6],vec![0.4,0.7],vec![0.3,0.8]],
        |p| if p[0] < 0.2 { Voxel::Solid(0) } else { Voxel::Fluid },
    ).unwrap();
    // Uniform Neumann data on each plane admit the full affine field, with
    // one cell thermostat setting only its otherwise-free constant offset.
    let faces = core::array::from_fn(|f| ThermalFace::HeatFlux(
        if f & 1 == 0 { flux[f/2] } else { -flux[f/2] },
    ));
    let mut setup = ThermalSetup::new(faces);
    setup.fixed_temperature.push((0, 350.0 + dot(center(&domain, 0), gradient)));
    let config = EnergyConfig { tolerance: 1e-11, ..EnergyConfig::default() };
    let field = vec![ConductivityTensor::new(k).unwrap(); domain.cell_count()];
    let solved = solve_tensor_conduction(
        &domain, &[SolidMaterial::new("solid", 1.0)], &setup, &field,
        &config, &CancelGate::new(),
    ).unwrap();
    for c in 0..domain.cell_count() {
        let expected = 350.0 + dot(center(&domain,c), gradient);
        assert!((solved.temperature[c]-expected).abs() < 2e-8, "cell {c}: {solved:?}");
        for f in 0..6 {
            let sign = if f & 1 == 0 { -1.0 } else { 1.0 };
            let tf = expected + sign*0.5*domain.widths(c)[f/2]*gradient[f/2];
            let q = sign*domain.face_area(c,f/2)*flux[f/2];
            assert!((solved.face_temperature_k[c][f]-tf).abs() < 2e-8);
            assert!((solved.outward_heat_w[c][f]-q).abs() < 2e-8);
            if let Some(n) = domain.neighbor(c,f) {
                assert_eq!(solved.face_temperature_k[c][f],
                    solved.face_temperature_k[n][f ^ 1]);
                assert!((solved.outward_heat_w[c][f]+solved.outward_heat_w[n][f ^ 1])
                    .abs() < 2e-8);
            }
        }
    }
    assert!((solved.temperature[0]-solved.temperature[1]).abs() < 2e-8);
    assert!(solved.outward_heat_w[0][1].abs() > 1.0);
    assert!((solved.solid_to_fluid_heat_w()-flux[0]*1.1*1.1).abs() < 2e-8);
    assert!(solved.report.max_heat_residual <= config.tolerance);
    assert!(solved.report.relative_residual <= config.tolerance);
    assert!(solved.report.balance.relative_residual < 2e-10);
}

#[test]
fn heterogeneous_rotated_tensors_resolve_both_contact_surfaces_and_all_heat_components() {
    let base = [[4.0,1.0,0.5],[1.0,3.0,0.2],[0.5,0.2,2.0]];
    let material = [SolidMaterial::new("A",1.0), SolidMaterial::new("B",1.0)];
    let domain = VoxelDomain::graded_from_fn(
        [vec![0.3,0.8],vec![0.2,0.5],vec![0.4,0.6]],
        |p| Voxel::Solid(if p[0] < 0.3 { 0 } else { 1 }),
    ).unwrap();
    let field: Vec<_> = (0..domain.cell_count()).map(|c| {
        let factor = if domain.coords(c)[0] == 0 { 1.0 } else { 3.0 };
        ConductivityTensor::new(base.map(|row| row.map(|v| factor*v))).unwrap()
    }).collect();
    let flux = [8.0,2.0,1.0]; // -K*grad(T), gradients -2 and -2/3 along x.
    for resistance in [0.0,0.2] {
        for explicit in [false,true] {
            let surface_a = 400.0-2.0*0.3;
            let surface_b = surface_a-resistance*flux[0];
            let cold = surface_b-(2.0/3.0)*0.8;
            let faces = [
                ThermalFace::Temperature(400.0),ThermalFace::Temperature(cold),
                ThermalFace::HeatFlux(flux[1]),ThermalFace::HeatFlux(-flux[1]),
                ThermalFace::HeatFlux(flux[2]),ThermalFace::HeatFlux(-flux[2]),
            ];
            let mut setup = ThermalSetup::new(faces);
            if explicit {
                for c in 0..domain.cell_count() {
                    if domain.coords(c)[0] == 0 {
                        setup.face_contacts.push(FaceContactResistance {
                            cells:[c,c+1], resistance_m2_k_w:resistance,
                        });
                    }
                }
            } else {
                setup.contacts.push(ContactResistance {
                    materials:(1,0), resistance_m2_k_w:resistance,
                });
            }
            let solved = solve_tensor_conduction(
                &domain, &material, &setup, &field, &EnergyConfig::default(), &CancelGate::new(),
            ).unwrap();
            for c in 0..domain.cell_count() {
                let left = domain.coords(c)[0] == 0;
                let expected = if left { 400.0-2.0*0.15 }
                    else { surface_b-(2.0/3.0)*0.4 };
                assert!((solved.temperature[c]-expected).abs() < 2e-8);
                for f in 0..6 {
                    let sign = if f & 1 == 0 { -1.0 } else { 1.0 };
                    assert!((solved.outward_heat_w[c][f]
                        -sign*domain.face_area(c,f/2)*flux[f/2]).abs() < 2e-8);
                }
                if left {
                    assert!((solved.face_temperature_k[c][1]-surface_a).abs() < 2e-8);
                    assert!((solved.face_temperature_k[c+1][0]-surface_b).abs() < 2e-8);
                    assert!((solved.outward_heat_w[c][1]+solved.outward_heat_w[c+1][0])
                        .abs() < 2e-8);
                }
            }
            assert!(solved.report.balance.relative_residual < 2e-10);
        }
    }
}

#[test]
fn diagonal_tensor_recovers_native_sources_film_sinks_and_face_contact() {
    let domain = VoxelDomain::graded_from_fn(
        [vec![0.1,0.3,0.6],vec![0.2],vec![0.3]], |_| Voxel::Solid(0),
    ).unwrap();
    let material = [SolidMaterial::new("isotropic",2.0)];
    let mut faces = [ThermalFace::Adiabatic;6];
    faces[0] = ThermalFace::HeatFlux(25.0);
    faces[1] = ThermalFace::Convective { h:3.0, ambient:300.0 };
    let mut setup = ThermalSetup::new(faces);
    setup.power_w = vec![1.0,2.0,3.0];
    setup.cell_sinks.push(CellSink { cell:1, conductance_w_k:0.5, temperature_k:290.0 });
    setup.face_contacts.push(FaceContactResistance {
        cells:[1,0],resistance_m2_k_w:0.2,
    });
    let config = EnergyConfig::default();
    let tensor = ConductivityTensor::new([[2.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,2.0]]).unwrap();
    let solved = solve_tensor_conduction(
        &domain, &material, &setup, &[tensor;3], &config, &CancelGate::new(),
    ).unwrap();
    let native = solve_energy(
        &domain, &FluidProperties::dry_air_300k(), &material, &FlowField::quiescent(&domain),
        &setup, &config, &CancelGate::new(),
    ).unwrap();
    for c in 0..3 {
        assert!((solved.temperature[c]-native.temperature[c]).abs() < 2e-8);
    }
    let actual_film = 3.0*domain.face_area(2,0)*(solved.face_temperature_k[2][1]-300.0);
    assert!((solved.outward_heat_w[2][1]-actual_film).abs() < 2e-9);
    assert!((solved.report.balance.sink_outflow_w-native.report.balance.sink_outflow_w)
        .abs() < 2e-9);
    assert!(solved.report.balance.residual_w.abs() < 2e-9);
}

#[test]
fn tensor_solver_refuses_ambiguous_owners_unanchored_states_and_unfunded_work() {
    let domain = VoxelDomain::from_fn(2,1,1,1.0, |_| Voxel::Solid(0)).unwrap();
    let tensor = ConductivityTensor::new([[2.0,0.5,0.0],[0.5,2.0,0.0],[0.0,0.0,1.0]])
        .unwrap();
    let material = [SolidMaterial::new("solid",1.0)];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic;6]);
    let config = EnergyConfig::default();
    let solve = |setup:&ThermalSetup, materials:&[SolidMaterial], config:&EnergyConfig,
        gate:&CancelGate| solve_tensor_conduction(
            &domain,materials,setup,&[tensor;2],config,gate,
        );
    assert!(matches!(solve(&setup,&material,&config,&CancelGate::new()),
        Err(ChtError::InvalidInput { .. })));
    setup.faces[0] = ThermalFace::Temperature(300.0);
    setup.faces[1] = ThermalFace::Temperature(400.0);
    let unfunded = EnergyConfig { max_iterations:0, ..config };
    assert!(matches!(solve(&setup,&material,&unfunded,&CancelGate::new()),
        Err(ChtError::SolverNotConverged { .. })));
    let nonlinear = [SolidMaterial::new("nonlinear",1.0)
        .with_conductivity_table(&[(300.0,1.0),(400.0,2.0)])];
    assert!(matches!(solve(&setup,&nonlinear,&config,&CancelGate::new()),
        Err(ChtError::InvalidInput { .. })));
    setup.eddy_conductivity_w_m_k = vec![0.0;2];
    assert!(matches!(solve(&setup,&material,&config,&CancelGate::new()),
        Err(ChtError::InvalidInput { .. })));
    setup.eddy_conductivity_w_m_k.clear();
    setup.faces[0] = ThermalFace::Inflow { temperature:300.0 };
    assert!(matches!(solve(&setup,&material,&config,&CancelGate::new()),
        Err(ChtError::InvalidInput { .. })));
    setup.faces[0] = ThermalFace::Temperature(300.0);
    setup.face_contacts.push(FaceContactResistance { cells:[0,1], resistance_m2_k_w:0.1 });
    setup.face_contacts.push(FaceContactResistance { cells:[1,0], resistance_m2_k_w:0.1 });
    assert!(matches!(solve(&setup,&material,&config,&CancelGate::new()),
        Err(ChtError::InvalidInput { .. })));
    setup.face_contacts.clear();
    let gate = CancelGate::new(); gate.request();
    assert_eq!(solve(&setup,&material,&config,&gate).unwrap_err(),ChtError::Cancelled);
}
