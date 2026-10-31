use super::*;
use crate::acoustics::helmholtz::boundary_layer::BoundaryLayerGeometry;
use crate::acoustics::helmholtz::boundary_layer::thermal_wall::{
    ConductiveThermalWall, WallThermalTransport,
};
use crate::acoustics::thermoviscous::wall::ConductiveWallRegion;
use crate::acoustics::thermoviscous::{ThermoviscousBoundary, ThermoviscousLimits};

const SIDE: f64 = 0.002;
const DEPTH: f64 = 0.00002;

fn specimen() -> (TetComplex, Vec<[f64; 3]>) {
    (
        TetComplex::from_tets(5, vec![[0, 1, 2, 4], [1, 3, 2, 4]]),
        vec![
            [0.0, 0.0, 0.0],
            [SIDE, 0.0, 0.0],
            [0.0, SIDE, 0.0],
            [SIDE, SIDE, 0.0],
            [0.0, 0.0, SIDE],
        ],
    )
}

fn model<'c>(
    mesh: &'c TetComplex,
    positions: &[[f64; 3]],
    faces: Vec<[u32; 3]>,
) -> ThermoviscousFe<'c> {
    let medium = ThermoviscousMedium {
        density_kg_m3: 1000.0,
        pressure_pa: 101325.0,
        temperature_k: 300.0,
        cp_j_kg_k: 4000.0,
        dynamic_viscosity_pa_s: 0.001,
        thermal_conductivity_w_m_k: 0.6,
        bulk_viscosity_pa_s: 0.002,
    };
    let geometry = HelmholtzFe::new(
        mesh,
        positions,
        2,
        Medium {
            density_kg_m3: medium.density_kg_m3,
            sound_speed_m_s: 1000.0,
        },
    )
    .unwrap();
    ThermoviscousFe::new_with_thermodynamics(
        mesh,
        positions,
        medium,
        ThermoviscousThermodynamics {
            isothermal_compressibility_pa_inv: 1e-9,
            isobaric_expansion_per_k: 0.0,
        },
        &ThermoviscousBoundary {
            no_slip_faces: geometry.boundary_faces(),
            isothermal_faces: vec![],
            ports: vec![],
        },
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap()
    .with_conductive_walls(
        &[ConductiveWallRegion {
            faces,
            wall: ConductiveThermalWall {
                transport: WallThermalTransport {
                    conductivity_w_m_k: 0.003,
                    volumetric_heat_capacity_j_m3_k: 2000.0,
                },
                depth: ThermalWallDepth::Slab {
                    thickness_m: DEPTH,
                    backing: ThermalBacking::Adiabatic,
                },
                minimum_tangential_scale_m: SIDE,
            },
            geometry: BoundaryLayerGeometry {
                minimum_curvature_radius_m: f64::INFINITY,
                minimum_half_gap_m: SIDE,
                maximum_layer_ratio: 0.05,
            },
            temperature_k: 300.0,
        }],
        &mut || false,
    )
    .unwrap()
}

fn resolution(elements: &[usize]) -> ThermoviscousWallDiscretization<'_> {
    ThermoviscousWallDiscretization {
        normal_elements: elements,
        maximum_frequency_hz: 2000.0,
    }
}

fn quadratic(a: &[f64], x: &[f64]) -> f64 {
    a.chunks_exact(x.len())
        .zip(x)
        .map(|(row, xi)| xi * row.iter().zip(x).map(|(v, xj)| v * xj).sum::<f64>())
        .sum()
}

#[test]
fn original_surface_and_depth_forms_match_the_planar_volume_heat_equation() {
    let (mesh, positions) = specimen();
    let model = model(&mesh, &positions, vec![[0, 1, 2], [1, 3, 2]]);
    let transport = [TangentialWallConductivity {
        isotropic_w_m_k: 0.05,
    }];
    for elements in [1, 2, 4] {
        let (normal, n) = WallFields::new(
            &model,
            1e-5,
            Some(resolution(&[elements])),
            None,
            model.n,
            &mut || false,
        )
        .unwrap();
        let (spatial, n2) = WallFields::new(
            &model,
            1e-5,
            Some(resolution(&[elements])),
            Some(&transport),
            model.n,
            &mut || false,
        )
        .unwrap();
        assert_eq!(n, n2);
        assert_eq!(
            spatial.layouts[0].trace_dofs().len(),
            9,
            "both triangles share their actual vertex/edge coordinates"
        );
        let mut mass = vec![0.0; n * n];
        let normal_k = normal
            .assemble(&model, &mut mass, n, &mut || false)
            .unwrap();
        let mut spatial_mass = vec![0.0; n * n];
        let spatial_k = spatial
            .assemble(&model, &mut spatial_mass, n, &mut || false)
            .unwrap();
        assert_eq!(mass, spatial_mass, "spatial conduction adds no heat store");
        let layout = &spatial.layouts[0];
        let mut field = vec![0.0; n];
        // Exact theta = .001 (x/L + .3 y/L) (1 + z/d) in the solid.
        // This polynomial has no hierarchical edge contribution.
        for depth in 0..=elements {
            for (surface, &dof) in layout.trace_dofs().iter().enumerate() {
                if dof < positions.len() {
                    let row = spatial.coordinate(&model, 0, depth, surface).unwrap();
                    field[row] = 0.001 * (positions[dof][0] + 0.3 * positions[dof][1]) / SIDE
                        * (1.0 + depth as f64 / elements as f64);
                }
            }
        }
        let lateral: Vec<_> = spatial_k
            .iter()
            .zip(&normal_k)
            .map(|(a, b)| a - b)
            .collect();
        let expected_tangent = 1e-6 * 300.0 * 0.05 * DEPTH * (7.0 / 3.0) * (1.0 + 0.3 * 0.3);
        let expected_normal =
            1e-6 * 300.0 * 0.003 * SIDE * SIDE / DEPTH * (1.0 / 3.0 + 0.3 / 2.0 + 0.09 / 3.0);
        assert!((quadratic(&lateral, &field) / expected_tangent - 1.0).abs() < 2e-10);
        assert!((quadratic(&normal_k, &field) / expected_normal - 1.0).abs() < 2e-12);
        let expected_storage = 1e-6
            * 2000.0
            * 300.0
            * SIDE
            * SIDE
            * DEPTH
            * (7.0 / 3.0)
            * (1.0 / 3.0 + 0.3 / 2.0 + 0.09 / 3.0);
        assert!((quadratic(&mass, &field) / expected_storage - 1.0).abs() < 2e-12);
    }
}

#[test]
fn shared_time_owner_spreads_heat_without_changing_storage_or_dropping_p2_fields() {
    let (mesh, positions) = specimen();
    let model = model(&mesh, &positions, vec![[0, 1, 2], [1, 3, 2]]);
    let transport = [TangentialWallConductivity {
        isotropic_w_m_k: 0.05,
    }];
    let prepare = |dt| {
        model
            .prepare_time_with_spatial_wall_fields(
                dt,
                model.medium.cp_j_kg_k,
                &[],
                resolution(&[2]),
                &transport,
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .unwrap()
    };
    let prepared = prepare(1e-5);
    let half = prepare(0.5e-5);
    let mut state = prepared.zero_state();
    let layout = &prepared.conductive_wall_layouts()[0];
    let edge = layout
        .trace_dofs()
        .iter()
        .position(|d| *d >= positions.len())
        .unwrap();
    state.conductive_wall_temperature_k[0][edge] = 0.01;
    let original = state.clone();
    let whole = prepared.step(&state, &[], &[], &mut || false).unwrap();
    let first = half.step(&state, &[], &[], &mut || false).unwrap();
    let second = half.step(&first.state, &[], &[], &mut || false).unwrap();
    assert_eq!(state, original);
    assert!(whole.work.solid_conductive_heat_j > 0.0);
    assert!(whole.work.final_energy_j < whole.work.initial_energy_j);
    assert!(whole.state.temperature_k.iter().any(|v| *v != 0.0));
    let scale = whole.work.initial_energy_j;
    assert!(whole.work.balance_residual_j().abs() < 1e-10 * scale);
    assert!(
        (whole.work.solid_conductive_heat_j
            - first.work.solid_conductive_heat_j
            - second.work.solid_conductive_heat_j)
            .abs()
            < 1e-9 * scale
    );
    for (a, b) in whole.state.conductive_wall_temperature_k[0]
        .iter()
        .zip(&second.state.conductive_wall_temperature_k[0])
    {
        assert!((a - b).abs() < 1e-10);
    }
    assert!(prepared.step(&state, &[], &[], &mut || true).is_err());
    assert_eq!(
        prepared
            .step(&state, &[], &[], &mut || false)
            .unwrap()
            .state,
        whole.state
    );
}

#[test]
fn spatial_wall_admission_checks_real_geometry_and_explicit_transport() {
    let (mesh, positions) = specimen();
    let planar = model(&mesh, &positions, vec![[0, 1, 2], [1, 3, 2]]);
    let bent = model(&mesh, &positions, vec![[0, 1, 2], [0, 1, 4]]);
    let good = [TangentialWallConductivity {
        isotropic_w_m_k: 0.05,
    }];
    assert!(admit(&planar, &good, &mut || false).is_ok());
    assert!(admit(&bent, &good, &mut || false).is_err());
    assert!(admit(&planar, &[], &mut || false).is_err());
    for k in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            admit(
                &planar,
                &[TangentialWallConductivity { isotropic_w_m_k: k }],
                &mut || false
            )
            .is_err()
        );
    }
    assert!(admit(&planar, &good, &mut || true).is_err());
}
