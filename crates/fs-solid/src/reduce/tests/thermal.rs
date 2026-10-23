//! G1/G3/G4: a real pressure-loaded valve from thermal/dead-load equilibrium.
//! The one-coordinate cube has closed-form P1 AND TET10 integrals, independent
//! of the production assembler or eigensolver. Synthetic values are deliberate.
use super::*;
use crate::linear3::{TetThermalStrainField, TetThermalStrainState, TetStaticSolveConfig};
use crate::reduce::thermal::{
    ThermalValveControls, reduce_pressure_loaded_valve as reduce_thermal,
};
use fs_scenario::reed_source::BeatingReedSource;

fn controls() -> ThermalValveControls {
    ThermalValveControls {
        solve: TetStaticSolveConfig {
            relative_residual_tolerance: 1e-10,
            maximum_iterations: 100,
        },
        mesh_frame: ContentHash([9; 32]),
        maximum_absolute_mandel_strain: 0.1,
        maximum_nodal_displacement_m: 0.1,
    }
}

fn strain(material: &TetElasticMaterial, epsilon: f64, bound: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        300.0 + 100.0 * epsilon,
        300.0,
        [epsilon, epsilon, epsilon, 0.0, 0.0, 0.0],
        material.material_state_identity,
        ContentHash([11; 32]),
        bound,
    ).unwrap()
}

fn request<'a>(
    nodes: &'a [[f64; 3]],
    tetrahedra: &'a [[usize; 4]],
    quadratic: Option<&'a [[usize; 10]]>,
    fixed: &'a [usize],
    face: &'a [usize],
    material: &'a TetElasticMaterial,
) -> ValveCardRequest<'a, TetMaterialField<'a>> {
    ValveCardRequest {
        nodes_m: nodes,
        tetrahedra,
        element_order: quadratic.map_or(ValveElementOrder::Tet4, |elements| {
            ValveElementOrder::Tet10 {
                elements,
                quadrature: TetQuadrature {
                    stiffness_points_per_dir: 3,
                    mass_points_per_dir: 4,
                },
            }
        }),
        fixed_dofs: fixed,
        material: TetMaterialField::Uniform(material),
        orifice: OrificeSpec {
            face_nodes: face,
            opening_axis: [1.0, 0.0, 0.0],
            width_axis: [0.0, 1.0, 0.0],
            opposing_plane_offset_m: 0.5,
        },
        window_hz: (0.01, 10.0),
        retained_compliance_floor: 0.999,
        loss_factor: 0.02,
        source_id: "synthetic thermal valve cube",
    }
}

fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 2e-9 * expected.abs().max(1e-5),
        "actual={actual:e}, expected={expected:e}");
}

#[test]
fn g1_thermal_valve_matches_closed_form_p1_and_tet10_pressure_response() {
    with_cx(|cx| {
        let material = TetElasticMaterial::try_new(2.0, 12.0, 0.0, [7; 32]).unwrap();
        let (corners, tets) = slab_mesh(1.0, 1.0, 1.0, 1, 1, 1);
        // Node 7 is the (1,1,1) corner. Its x component is the ONLY free DOF.
        // The observation uses four explicit face nodes, including fixed nodes.
        let face = [4, 5, 6, 7];
        for quadratic in [false, true] {
            let (nodes, elements) = if quadratic {
                let (nodes, elements) = crate::linear3::tet10::promote_p1_to_p2(
                    &corners, &tets, None, cx,
                ).unwrap();
                (nodes, Some(elements))
            } else { (corners.clone(), None) };
            let fixed: Vec<_> = (0..3 * nodes.len()).filter(|d| *d != 21).collect();
            let request = request(&nodes, &tets, elements.as_deref(), &fixed, &face, &material);
            let mut pressure = vec![[0.0; 3]; nodes.len()];
            pressure[7][0] = 1.0;
            let baseline = crate::reduce::material_field::reduce_pressure_loaded_valve(
                &request, &pressure, cx,
            ).unwrap();
            let zero = strain(&material, 0.0, 0.1);
            let unchanged = reduce_thermal(
                &request, TetThermalStrainField::Uniform(&zero), &[],
                &pressure, controls(), cx,
            ).unwrap();
            assert_eq!(unchanged.valve().card(), baseline.card());

            let epsilon = -0.02;
            let thermal = strain(&material, epsilon, 0.1);
            for force in [0.0, -0.048] {
                let reduced = reduce_thermal(
                    &request, TetThermalStrainField::Uniform(&thermal), &[(21, force)],
                    &pressure, controls(), cx,
                ).unwrap();

                // For the six Freudenthal tets: grad(lambda_7) is x, y or z,
                // twice each, and V=1/6. P1: K_e=2E/3, M=rho/10, u=eps/2+f/K_e,
                // K_g=E(u-3eps)/3. TET10 corner N=lambda(2lambda-1):
                // int N²=V/70, int(4lambda-1)²=3V/5, int(4lambda-1)³=2V/5,
                // int grad N=0. Thus K_e=2E/5 and K_g=2Eu/15-3Eeps/5.
                let (elastic, mass, displacement, geometric) = if quadratic {
                    let elastic = 2.0 * 12.0 / 5.0;
                    let displacement = force / elastic;
                    (elastic, 2.0 / 70.0, displacement,
                        2.0 * 12.0 * displacement / 15.0 - 3.0 * 12.0 * epsilon / 5.0)
                } else {
                    let elastic = 2.0 * 12.0 / 3.0;
                    let displacement = epsilon / 2.0 + force / elastic;
                    (elastic, 2.0 / 10.0, displacement,
                        12.0 * (displacement - 3.0 * epsilon) / 3.0)
                };
                let tangent = elastic + geometric;
                let card = reduced.valve().card();
                let mode = &card.modes[0];
                assert_eq!(card.modes.len(), 1);
                near(mode.frequency_hz, (tangent / mass).sqrt() / core::f64::consts::TAU);
                near(mode.mass_kg, 16.0 * mass);
                near(mode.stiffness_n_m, 16.0 * tangent);
                near(card.total_mass_kg, 2.0);
                near(card.rest_gap_m, 0.5 + displacement.min(0.0));
                near(reduced.equilibrium().displacement_m()[7][0], displacement);
                near(reduced.observation_nodes_m()[7][0], 1.0 + displacement);
                near(reduced.valve().effective_areas_m2()[0], 4.0);
                near(reduced.valve().pressure_compliance().exact_volume_m3_per_pa(), 1.0 / tangent);
                near(reduced.valve().pressure_compliance().exact_opening_m_per_pa(), 1.0 / (4.0 * tangent));
                assert_ne!(reduced.identity(), unchanged.identity());
                reduced.compile_beating_reed(0.01, 0.001, &SingleModeValveAdmission {
                    minimum_opening_compliance_fraction: 0.99,
                    minimum_pressure_compliance_fraction: 0.99,
                    maximum_static_opening_relative_error: 1e-8,
                }).unwrap();
                let repeated = reduce_thermal(
                    &request, TetThermalStrainField::Uniform(&thermal), &[(21, force)],
                    &pressure, controls(), cx,
                ).unwrap();
                assert_eq!(repeated.identity(), reduced.identity());
            }
        }
    });
}

#[test]
fn g3_thermal_valve_rejects_instability_outside_the_selected_positive_window() {
    with_cx(|cx| {
        // Two disconnected, separately restrained cubes. Only the second is
        // observed/pressure-loaded. Its positive compliance cannot conceal the
        // unstable first cube. This deliberately extreme synthetic strain
        // tests algebraic admission, not a physical small-strain validity claim.
        let material = TetElasticMaterial::try_new(2.0, 12.0, 0.0, [7; 32]).unwrap();
        let (mut nodes, mut tets) = slab_mesh(1.0, 1.0, 1.0, 1, 1, 1);
        let more_nodes: Vec<_> = nodes.iter().map(|p| [p[0] + 3.0, p[1], p[2]]).collect();
        let more_tets: Vec<_> = tets.iter().map(|t| t.map(|n| n + 8)).collect();
        nodes.extend(more_nodes);
        tets.extend(more_tets);
        let fixed: Vec<_> = (0..3 * nodes.len()).filter(|d| *d != 21 && *d != 45).collect();
        let face = [12, 13, 14, 15];
        let mut request = request(&nodes, &tets, None, &fixed, &face, &material);
        request.orifice.opposing_plane_offset_m = 3.5;
        let hot = strain(&material, 1.0, 2.0);
        let cold = strain(&material, 0.0, 2.0);
        let thermal: Vec<_> = (0..12).map(|i| if i < 6 { hot.clone() } else { cold.clone() }).collect();
        let mut pressure = vec![[0.0; 3]; nodes.len()];
        pressure[15][0] = 1.0;
        let mut limits = controls();
        limits.maximum_absolute_mandel_strain = 2.0;
        // f=-E*eps/3 leaves u=0 on the first cube. K_t=2E/3-E<0,
        // while the second cube retains K_t=2E/3>0 in the requested window.
        let result = reduce_thermal(
            &request, TetThermalStrainField::PerElement(&thermal), &[(21, -4.0)],
            &pressure, limits, cx,
        );
        assert!(matches!(result, Err(ReduceError::Invalid {
            what: "thermal prestress tangent must be positive definite before valve reduction"
        })));
    });
}

#[test]
fn g4_thermal_valve_refuses_displacement_closed_aperture_and_source_mismatch() {
    with_cx(|cx| {
        let material = TetElasticMaterial::try_new(2.0, 12.0, 0.0, [7; 32]).unwrap();
        let (nodes, tets) = slab_mesh(1.0, 1.0, 1.0, 1, 1, 1);
        let original_nodes = nodes.clone();
        let fixed: Vec<_> = (0..3 * nodes.len()).filter(|d| *d != 21).collect();
        let face = [4, 5, 6, 7];
        let mut request = request(&nodes, &tets, None, &fixed, &face, &material);
        let thermal = strain(&material, -0.02, 0.1);
        let mut pressure = vec![[0.0; 3]; nodes.len()];
        pressure[7][0] = 1.0;
        let mut limits = controls();
        limits.maximum_nodal_displacement_m = 0.005;
        assert!(matches!(reduce_thermal(
            &request, TetThermalStrainField::Uniform(&thermal), &[],
            &pressure, limits, cx,
        ), Err(ReduceError::Assembly(crate::linear3::TetElasticError::InvalidGeometryUpdate { .. }))));
        request.orifice.opposing_plane_offset_m = 0.999;
        assert!(matches!(reduce_thermal(
            &request, TetThermalStrainField::Uniform(&thermal), &[],
            &pressure, controls(), cx,
        ), Err(ReduceError::Invalid {
            what: "measured rest gap is not positive (face behind the opposing plane)"
        })));
        let other = TetElasticMaterial::try_new(2.0, 12.0, 0.0, [8; 32]).unwrap();
        let foreign = strain(&other, -0.02, 0.1);
        assert!(reduce_thermal(
            &request, TetThermalStrainField::Uniform(&foreign), &[],
            &pressure, controls(), cx,
        ).is_err());
        assert_eq!(nodes, original_nodes);
    });
}

#[test]
fn g4_thermal_valve_cancellation_publishes_no_reduction() {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate, arena,
            fs_exec::StreamKey { seed: 0, kernel_id: 78, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic,
        );
        let material = TetElasticMaterial::try_new(2.0, 12.0, 0.0, [7; 32]).unwrap();
        let (nodes, tets) = slab_mesh(1.0, 1.0, 1.0, 1, 1, 1);
        let fixed: Vec<_> = (0..3 * nodes.len()).filter(|d| *d != 21).collect();
        let face = [4, 5, 6, 7];
        let request = request(&nodes, &tets, None, &fixed, &face, &material);
        let thermal = strain(&material, -0.02, 0.1);
        gate.request();
        assert!(matches!(reduce_thermal(
            &request, TetThermalStrainField::Uniform(&thermal), &[],
            &[[0.0; 3]; 8], controls(), &cx,
        ), Err(ReduceError::Assembly(crate::linear3::TetElasticError::Cancelled))));
    });
}
