use super::*;
use crate::conjugate::{FlowFace, FluidProperties, Voxel};

fn setup(reverse: bool) -> (VoxelDomain, FlowField, ThermalSetup) {
    let domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |_| Voxel::Fluid).unwrap();
    let mut faces = [FlowFace::Wall; 6];
    faces[0] = FlowFace::Fixed;
    faces[1] = FlowFace::Fixed;
    let flow = FlowField::from_face_velocity(&domain, faces, |axis, _| {
        if axis == 0 {
            if reverse { -1.0 } else { 1.0 }
        } else {
            0.0
        }
    })
    .unwrap();
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[usize::from(reverse)] = ThermalFace::Inflow { temperature: 300.0 };
    faces[usize::from(!reverse)] = ThermalFace::Outflow {
        backflow_temperature: 300.0,
    };
    let mut thermal = ThermalSetup::new(faces);
    thermal.power_w = vec![0.0, 4.5, 0.0];
    (domain, flow, thermal)
}

fn properties(reference_cp: f64) -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: reference_cp,
        conductivity_w_m_k: 1.0,
        kinematic_viscosity_m2_s: 1.0,
    }
}

#[test]
fn g1_g3_each_face_uses_one_donor_enthalpy_with_datum_and_direction_invariance() {
    // Independent three-cell equations in theta = T - 300:
    // 14 a - b = 0; -2 a + 24 b - c = 4.5; -3 b + 45 c = 0.
    // Different donor Cp values are deliberate: using receiver Cp in one row
    // no longer telescopes and fails both this solution and physical closure.
    let middle = 4.5 / (24.0 - 2.0 / 14.0 - 3.0 / 45.0);
    let expected = [middle / 14.0, middle, middle * 3.0 / 45.0];
    for reverse in [false, true] {
        let (domain, flow, setup) = setup(reverse);
        for datum in [0.0, 1024.0] {
            let mut cp = [1.0, 2.0, 4.0];
            if reverse {
                cp.reverse();
            }
            let cells = cp.map(|c| {
                Some(AffineSpecificEnthalpy {
                    temperature_k: 300.0,
                    enthalpy_j_kg: datum,
                    heat_capacity_j_kg_k: c,
                })
            });
            let field = AdvectiveEnthalpyField {
                cells: &cells,
                boundary_j_kg: [Some(datum), Some(datum), None, None, None, None],
            };
            // This reference Cp must not leak into the admitted caloric field.
            let config = EnergyConfig {
                scheme: ConvectionScheme::Upwind,
                ..EnergyConfig::default()
            };
            let (solution, record) = config
                .step_with_enthalpy(
                    &domain,
                    &properties(999.0),
                    &[],
                    &flow,
                    &setup,
                    &[[1.0; 3]; 3],
                    field,
                    &cp,
                    &[300.0; 3],
                    0.1,
                    &CancelGate::new(),
                )
                .unwrap();
            for cell in 0..3 {
                let index = if reverse { 2 - cell } else { cell };
                assert!((solution.temperature[cell] - 300.0 - expected[index]).abs() < 1e-9);
            }
            assert!(record.closure_j.abs() < 1e-9);
            assert!((record.source_j - 0.45).abs() < 1e-12);
            assert!((solution.report.balance.advective_outflow_w - 4.0 * expected[2]).abs() < 1e-9);
        }
    }
}

#[test]
fn g0_invalid_caloric_fields_and_nonconservative_flow_refuse_before_publication() {
    let (domain, flow, setup) = setup(false);
    let h = AffineSpecificEnthalpy {
        temperature_k: 300.0,
        enthalpy_j_kg: 0.0,
        heat_capacity_j_kg_k: 1.0,
    };
    let mut cells = [Some(h); 3];
    let boundary = [Some(0.0), Some(0.0), None, None, None, None];
    let config = EnergyConfig {
        scheme: ConvectionScheme::Upwind,
        ..EnergyConfig::default()
    };
    let solve = |config: &EnergyConfig,
                 flow: &FlowField,
                 field: AdvectiveEnthalpyField<'_>,
                 gate: &CancelGate| {
        config.step_with_enthalpy(
            &domain,
            &properties(1.0),
            &[],
            flow,
            &setup,
            &[[1.0; 3]; 3],
            field,
            &[1.0; 3],
            &[300.0; 3],
            0.1,
            gate,
        )
    };
    let field = AdvectiveEnthalpyField {
        cells: &cells,
        boundary_j_kg: boundary,
    };
    assert!(solve(&EnergyConfig::default(), &flow, field, &CancelGate::new()).is_err());
    let divergent = FlowField::from_face_velocity(
        &domain,
        [
            FlowFace::Fixed,
            FlowFace::Fixed,
            FlowFace::Wall,
            FlowFace::Wall,
            FlowFace::Wall,
            FlowFace::Wall,
        ],
        |axis, point| if axis == 0 { 1.0 + point[0] } else { 0.0 },
    )
    .unwrap();
    assert!(solve(&config, &divergent, field, &CancelGate::new()).is_err());
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        solve(&config, &flow, field, &cancelled),
        Err(ChtError::Cancelled)
    ));
    cells[1].as_mut().unwrap().heat_capacity_j_kg_k = -1.0;
    assert!(
        solve(
            &config,
            &flow,
            AdvectiveEnthalpyField {
                cells: &cells,
                boundary_j_kg: boundary
            },
            &CancelGate::new()
        )
        .is_err()
    );
}

#[test]
fn g1_backflow_uses_reservoir_enthalpy_without_inventing_diffusion() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Fluid).unwrap();
    let flow = FlowField::from_face_velocity(
        &domain,
        [
            FlowFace::Fixed,
            FlowFace::Fixed,
            FlowFace::Wall,
            FlowFace::Wall,
            FlowFace::Wall,
            FlowFace::Wall,
        ],
        |axis, _| if axis == 0 { 1.0 } else { 0.0 },
    )
    .unwrap();
    let mut faces = [ThermalFace::Adiabatic; 6];
    // The left OUTflow face now receives fluid. Its supplied h(310) is 20,
    // not the cell tangent Cp=1 times a boundary temperature.
    faces[0] = ThermalFace::Outflow {
        backflow_temperature: 310.0,
    };
    faces[1] = ThermalFace::Outflow {
        backflow_temperature: 300.0,
    };
    let cells = [Some(AffineSpecificEnthalpy {
        temperature_k: 300.0,
        enthalpy_j_kg: 0.0,
        heat_capacity_j_kg_k: 1.0,
    })];
    let field = AdvectiveEnthalpyField {
        cells: &cells,
        boundary_j_kg: [Some(20.0), Some(0.0), None, None, None, None],
    };
    let config = EnergyConfig {
        scheme: ConvectionScheme::Upwind,
        ..EnergyConfig::default()
    };
    let (solution, record) = config
        .step_with_enthalpy(
            &domain,
            &properties(999.0),
            &[],
            &flow,
            &ThermalSetup::new(faces),
            &[[100.0; 3]],
            field,
            &[1.0],
            &[300.0],
            0.1,
            &CancelGate::new(),
        )
        .unwrap();
    // Ten units of storage rate, one of outgoing mass enthalpy: 11 theta=20.
    assert!((solution.temperature[0] - 300.0 - 20.0 / 11.0).abs() < 1e-9);
    assert!(record.closure_j.abs() < 1e-9);
    assert!((record.stored_energy_change_j + record.boundary_outflow_j).abs() < 1e-9);
}
