//! G1/G3: package/TIM series laws, actual source surfaces and one storage history.
//! These manufactured material laws are numerical references, not measured data.

use fs_exec::CancelGate;
use fs_lbm::Face3;
use fs_lbm::conjugate::{
    ChtError, CompactComponent, EnergyConfig, FaceContactResistance, FlowField, FluidProperties,
    SolidMaterial, SourceFaceContact, ThermalFace, ThermalSetup, TransientConfig, Voxel,
    VoxelDomain, march_energy, solve_energy,
};

fn fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0,
        kinematic_viscosity_m2_s: 1.0,
    }
}

fn component(power_w: f64) -> CompactComponent {
    CompactComponent {
        lo: [0, 0, 1],
        hi: [1, 1, 2],
        board_face: Face3::ZMin,
        power_w,
        junction_to_case_k_w: 1.0,
        junction_to_board_k_w: 1.0,
        heat_capacity_j_k: Some(6.0),
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-7,
        "{actual:.15} != {expected:.15}"
    );
}

fn source_contact(
    cells: [usize; 2],
    case_span: [f64; 2],
    outside_span: [f64; 2],
    r_at_320: f64,
    slope: f64,
    reverse: bool,
) -> SourceFaceContact {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
        PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
        UncertaintyModel,
    };
    use fs_material::state_point::MaterialPropertySelection;
    use fs_material::thermal::interface::{
        AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as DIMS,
        AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as PROPERTY,
    };
    use fs_qty::{
        QuantitySpec, Temperature,
        semantic::{QuantityKind, SemanticType, ValueForm},
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(PROPERTY, DIMS),
            value: PropertyValue::Curve {
                abscissa: "TA".into(),
                abscissa_dims: Temperature::DIMS,
                knots: case_span
                    .map(|t| (t, r_at_320 + slope * (t - 320.0)))
                    .to_vec(),
                dims: DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with("TA", case_span[0], case_span[1])
                .with("TB", outside_span[0], outside_span[1]),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured compact/TIM law".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        texture_frame: "declared Cartesian normal".into(),
    };
    let names = if reverse {
        ["outside", "package"]
    } else {
        ["package", "outside"]
    };
    let card = InterfaceSystemCard::assemble(
        surface(names[0]),
        surface(names[1]),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "fixture".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let kelvin = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let point = QueryPoint::new()
        .with_quantity(
            "TA",
            kelvin,
            case_span[0] + 0.25 * (case_span[1] - case_span[0]),
        )
        .unwrap()
        .with_quantity("TB", kelvin, 0.5 * (outside_span[0] + outside_span[1]))
        .unwrap();
    SourceFaceContact::new(
        cells,
        std::sync::Arc::new(card),
        &point,
        if reverse { ["TB", "TA"] } else { ["TA", "TB"] },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn graded_package_contacts_match_independent_parallel_series_network() {
    let domain =
        VoxelDomain::graded_from_fn([vec![1.0, 2.0], vec![1.0], vec![0.5, 0.3, 1.5]], |_| {
            Voxel::Solid(0)
        })
        .unwrap();
    let solids = [SolidMaterial::new("grid-aligned solid", 2.0).with_orthotropic([5.0, 7.0, 2.0])];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.fixed_temperature = vec![(0, 300.0), (1, 300.0), (4, 310.0), (5, 310.0)];
    setup.compact_components.push(CompactComponent {
        hi: [2, 1, 2],
        junction_to_board_k_w: 2.0,
        junction_to_case_k_w: 5.0,
        ..component(7.0)
    });
    setup.face_contacts = vec![
        FaceContactResistance {
            cells: [0, 2],
            resistance_m2_k_w: 1.0,
        },
        FaceContactResistance {
            cells: [1, 3],
            resistance_m2_k_w: 4.0,
        },
        FaceContactResistance {
            cells: [2, 4],
            resistance_m2_k_w: 2.0,
        },
        FaceContactResistance {
            cells: [3, 5],
            resistance_m2_k_w: 6.0,
        },
    ];
    let solve = |setup: &ThermalSetup| {
        solve_energy(
            &domain,
            &fluid(),
            &solids,
            &FlowField::quiescent(&domain),
            setup,
            &EnergyConfig::default(),
            &CancelGate::new(),
        )
        .unwrap()
    };
    let result = solve(&setup);
    // Total package area is 3 m²; each R_j branch splits by its patch area.
    // TIM R'' remains area-specific, and only the OUTSIDE half cell conducts.
    let branches = [
        (0, 1.0, 2.0, 1.0, 0.25, 300.0),
        (1, 2.0, 2.0, 4.0, 0.25, 300.0),
        (4, 1.0, 5.0, 2.0, 0.75, 310.0),
        (5, 2.0, 5.0, 6.0, 0.75, 310.0),
    ];
    let conductance: Vec<_> = branches
        .iter()
        .map(|&(_, a, r, contact, half, _)| a / (3.0 * r + contact + half / 2.0))
        .collect();
    let expected = (7.0
        + conductance
            .iter()
            .zip(branches)
            .map(|(g, b)| g * b.5)
            .sum::<f64>())
        / conductance.iter().sum::<f64>();
    close(result.junctions[0].temperature_k, expected);
    let mut branch_heat = [0.0; 2];
    for ((cell, _, _, _, _, reservoir), g) in branches.into_iter().zip(conductance) {
        let actual = result.junctions[0]
            .links
            .iter()
            .find(|(c, _)| *c == cell)
            .unwrap()
            .1;
        close(actual, g);
        branch_heat[usize::from(cell >= 4)] += g * (expected - reservoir);
    }
    close(result.junctions[0].board_w, branch_heat[0]);
    close(result.junctions[0].case_w, branch_heat[1]);
    close(branch_heat.iter().sum(), 7.0);
    close(result.report.balance.fixed_cell_injection_w, -7.0);
    close(result.report.balance.residual_w, 0.0);
    assert_eq!(result.temperature[2], result.temperature[3]);
    for contact in &mut setup.face_contacts {
        contact.cells.reverse();
    }
    setup.face_contacts.reverse();
    assert_eq!(result, solve(&setup));
}

#[test]
fn nonlinear_contact_queries_package_surface_after_its_junction_resistor() {
    let domain = VoxelDomain::from_fn(1, 1, 3, 1.0, |p| {
        Voxel::Solid(if p[2] < 1.0 {
            1
        } else if p[2] < 2.0 {
            2
        } else {
            0
        })
    })
    .unwrap();
    let solids = [
        SolidMaterial::new("outside k(T)", 1.0)
            .with_conductivity_table(&[(250.0, 0.5), (450.0, 2.5)]),
        SolidMaterial::new("board", 1.0),
        SolidMaterial::new("collapsed alias", 1e-6),
    ];
    // Manufacture a 20 W case branch: T_case=320, T_outside_surface=310,
    // R_TIM=0.5, R_jc=2 -> T_junction=360. Integral(k dT) across the
    // outside half cell is q*0.5=10, giving Phi(T_cell)=Phi(310)-10=0.5.
    let outside = 300.0 + 1.0 / (1.0 + 1.01_f64.sqrt());
    let board_w = 60.0 / 3.5;
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.fixed_temperature = vec![(0, 300.0), (2, outside)];
    setup.compact_components.push(CompactComponent {
        junction_to_board_k_w: 3.0,
        junction_to_case_k_w: 2.0,
        ..component(20.0 + board_w)
    });
    let mut previous = None;
    for reverse in [false, true] {
        setup.source_face_contacts = vec![source_contact(
            if reverse { [2, 1] } else { [1, 2] },
            [318.0, 322.0],
            [308.0, 312.0],
            0.5,
            0.01,
            reverse,
        )];
        let before = setup.clone();
        let solution = solve_energy(
            &domain,
            &fluid(),
            &solids,
            &FlowField::quiescent(&domain),
            &setup,
            &EnergyConfig::default(),
            &CancelGate::new(),
        )
        .unwrap();
        close(solution.junctions[0].temperature_k, 360.0);
        close(solution.junctions[0].case_w, 20.0);
        close(solution.junctions[0].board_w, board_w);
        let state = &solution.source_contact_states[0];
        let expected = if reverse {
            [310.0, 320.0]
        } else {
            [320.0, 310.0]
        };
        for (&actual, target) in state.surface_temperature_k.iter().zip(expected) {
            close(actual, target);
        }
        close(state.source.state_coordinate("TA").unwrap(), 320.0);
        close(state.source.state_coordinate("TB").unwrap(), 310.0);
        close(
            state.heat_rate_first_to_second_w,
            if reverse { -20.0 } else { 20.0 },
        );
        close(state.resistance_m2_k_w, 0.5);
        assert!(solution.source_contact_iterations > 1);
        assert!(
            solution.temperature[1] > 322.0 && solution.temperature[2] < 308.0,
            "neither bulk/junction value is inside the interface source box"
        );
        close(solution.report.balance.residual_w, 0.0);
        if let Some(old) = &previous {
            let old: &Vec<f64> = old;
            for (&a, &b) in old.iter().zip(&solution.temperature) {
                close(a, b);
            }
        }
        previous = Some(solution.temperature.clone());
        assert_eq!(setup, before);
    }
}

#[test]
fn fixed_and_source_contacts_keep_one_closed_junction_storage_step() {
    let domain = VoxelDomain::from_fn(1, 1, 3, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("ordinary solid", 2.0).with_heat_capacity(2.0)];
    let previous = [300.0, 400.0, 300.0];
    for source in [false, true] {
        let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
        setup.compact_components.push(component(0.0));
        for cells in [[1, 0], [1, 2]] {
            if source {
                setup.source_face_contacts.push(source_contact(
                    cells,
                    [250.0, 450.0],
                    [250.0, 450.0],
                    3.0,
                    0.01,
                    false,
                ));
            } else {
                setup.face_contacts.push(FaceContactResistance {
                    cells,
                    resistance_m2_k_w: 3.0,
                });
            }
        }
        let before = setup.clone();
        // Eliminate the two outside storage rows and the one Cj=6 row at
        // dt=3: Tj=400-q, Tout=300+1.5q, Tcase=400-2q. The nonlinear
        // contact R=3+0.01(Tcase-320) gives 0.02q²-7.55q+100=0.
        let q = if source {
            200.0 / (7.55 + (7.55_f64 * 7.55 - 8.0).sqrt())
        } else {
            100.0 / 6.75
        };
        let expected = [300.0 + 1.5 * q, 400.0 - q, 300.0 + 1.5 * q];
        let (solution, record) = EnergyConfig::default()
            .step_with_conductivity(
                &domain,
                &fluid(),
                &solids,
                &FlowField::quiescent(&domain),
                &setup,
                &[[2.0; 3]; 3],
                &[2.0, 0.0, 2.0],
                &previous,
                3.0,
                &CancelGate::new(),
            )
            .unwrap();
        for (&actual, target) in solution.temperature.iter().zip(expected) {
            close(actual, target);
        }
        close(solution.junctions[0].case_w, q);
        close(solution.junctions[0].board_w, q);
        close(record.stored_energy_change_j, 0.0);
        close(record.source_j, 0.0);
        close(record.boundary_outflow_j, 0.0);
        close(record.closure_j, 0.0);
        close(record.time_s, 3.0);
        if source {
            assert!(solution.source_contact_iterations > 1);
        }
        let native = march_energy(
            &domain,
            &fluid(),
            &solids,
            &FlowField::quiescent(&domain),
            &setup,
            &previous,
            |_| 1.0,
            &TransientConfig {
                time_step_s: 3.0,
                steps: 1,
                energy: EnergyConfig::default(),
            },
            None,
            &CancelGate::new(),
        )
        .unwrap();
        for (&actual, target) in native.temperature.iter().zip(expected) {
            close(actual, target);
        }
        close(native.records[0].closure_j, 0.0);
        assert_eq!(setup, before);
        assert_eq!(previous, [300.0, 400.0, 300.0]);
    }
}

#[test]
fn contacts_do_not_open_compact_sides_internal_aliases_or_duplicate_owners() {
    let domain = VoxelDomain::from_fn(3, 1, 4, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("solid", 1.0)];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.faces[4] = ThermalFace::Temperature(300.0);
    setup.faces[5] = ThermalFace::Temperature(300.0);
    setup.compact_components.push(CompactComponent {
        lo: [1, 0, 1],
        hi: [2, 1, 3],
        ..component(4.0)
    });
    let solve = |setup: &ThermalSetup| {
        solve_energy(
            &domain,
            &fluid(),
            &solids,
            &FlowField::quiescent(&domain),
            setup,
            &EnergyConfig::default(),
            &CancelGate::new(),
        )
    };
    assert!(solve(&setup).is_ok());
    let alias = domain.index(1, 0, 1);
    for other in [domain.index(0, 0, 1), domain.index(1, 0, 2)] {
        for source in [false, true] {
            let mut bad = setup.clone();
            if source {
                bad.source_face_contacts.push(source_contact(
                    [alias, other],
                    [250.0, 450.0],
                    [250.0, 450.0],
                    0.5,
                    0.0,
                    false,
                ));
            } else {
                bad.face_contacts.push(FaceContactResistance {
                    cells: [alias, other],
                    resistance_m2_k_w: 0.5,
                });
            }
            assert!(matches!(
                solve(&bad),
                Err(ChtError::InvalidInput {
                    field: "thermal.face_contacts",
                    ..
                })
            ));
        }
    }
    let outside = domain.index(1, 0, 0);
    setup.face_contacts.push(FaceContactResistance {
        cells: [alias, outside],
        resistance_m2_k_w: 0.5,
    });
    setup.source_face_contacts.push(source_contact(
        [outside, alias],
        [250.0, 450.0],
        [250.0, 450.0],
        0.5,
        0.0,
        true,
    ));
    assert!(matches!(
        solve(&setup),
        Err(ChtError::InvalidInput {
            field: "thermal.face_contacts",
            ..
        })
    ));
}
