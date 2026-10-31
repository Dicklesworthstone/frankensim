use super::super::wall::{MaterialTangentialWallRegion, MaterialWallRegion};
use super::*;
use fs_feec::acoustics::helmholtz::boundary_layer::BoundaryLayerGeometry;
use fs_feec::acoustics::helmholtz::boundary_layer::thermal_wall::{
    ThermalBacking, ThermalWallDepth,
};

#[test]
fn declared_source_isotropy_reaches_the_original_spatial_wall_and_rejects_foreign_material() {
    let (mesh, points0) = kuhn_cube(1);
    let xyz = points(&points0);
    let gas = gas();
    let geometry = HelmholtzFe::new(
        &mesh,
        &xyz,
        2,
        Medium {
            density_kg_m3: gas.density,
            sound_speed_m_s: gas.sound_speed,
        },
    )
    .unwrap();
    let face = geometry
        .boundary_faces()
        .into_iter()
        .find(|f| f.iter().all(|&v| xyz[v as usize][1] == 0.0))
        .unwrap();
    let wall = time::wall_transport(1000.0, 2.0, 0.003);
    let foreign = time::wall_transport(1000.0, 2.0, 0.012);
    let sample = wall.at_temperature(gas.temperature).unwrap();
    let make = || {
        ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, false), gas, BAND)
            .unwrap()
            .with_conductive_walls(
                &[MaterialWallRegion {
                    faces: &[face],
                    sample: &sample,
                    depth: ThermalWallDepth::Slab {
                        thickness_m: 0.00002,
                        backing: ThermalBacking::Adiabatic,
                    },
                    minimum_tangential_scale_m: 0.004,
                    geometry: BoundaryLayerGeometry {
                        minimum_curvature_radius_m: f64::INFINITY,
                        minimum_half_gap_m: 0.0001,
                        maximum_layer_ratio: 0.05,
                    },
                }],
                &mut || false,
            )
            .unwrap()
    };
    assert!(
        make()
            .with_tangentially_isotropic_wall_transport(
                &[MaterialTangentialWallRegion {
                    normal: &wall,
                    tangential: &foreign,
                }],
                &mut || false
            )
            .is_err()
    );
    assert!(
        make()
            .with_tangentially_isotropic_wall_transport(
                &[MaterialTangentialWallRegion {
                    normal: &foreign,
                    tangential: &wall,
                }],
                &mut || false
            )
            .is_err()
    );
    let normal_source = make();
    // Isotropy is declared explicitly even though this fixture selects the same
    // scalar source in all directions; no scalar normal sample implies it.
    let spatial_source = make()
        .with_tangentially_isotropic_wall_transport(
            &[MaterialTangentialWallRegion {
                normal: &wall,
                tangential: &wall,
            }],
            &mut || false,
        )
        .unwrap();
    assert_eq!(spatial_source.wall_tangential_sources(), &[sample]);
    assert!(
        spatial_source
            .response(BAND[0], &[C64::ZERO, C64::ZERO], &mut || false)
            .is_err()
    );
    let normal = normal_source
        .prepare_time_with_wall_diffusion(
            1e-7,
            &[2],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let spatial = spatial_source
        .prepare_time_with_wall_diffusion(
            1e-7,
            &[2],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let mut initial = spatial.zero_state();
    for (value, point) in initial.temperature_k.iter_mut().zip(&xyz) {
        *value = 0.01 * point[0] / LENGTH;
    }
    let layout = &spatial.conductive_wall_layouts()[0];
    for depth in initial.conductive_wall_temperature_k[0].chunks_mut(layout.trace_dofs().len()) {
        for (value, &dof) in depth.iter_mut().zip(layout.trace_dofs()) {
            *value = initial.temperature_k[dof];
        }
    }
    assert_eq!(
        spatial.stored_energy_j(&initial).unwrap(),
        normal.stored_energy_j(&initial).unwrap()
    );
    let a = normal
        .step(&initial, &[0.0, 0.0], &[], &mut || false)
        .unwrap();
    let b = spatial
        .step(&initial, &[0.0, 0.0], &[], &mut || false)
        .unwrap();
    assert!(b.work.solid_conductive_heat_j > a.work.solid_conductive_heat_j);
    assert_ne!(
        a.state.conductive_wall_temperature_k,
        b.state.conductive_wall_temperature_k
    );
    assert!(b.work.balance_residual_j().abs() < 1e-10 * b.work.initial_energy_j);
    // The same FE pointer cannot authorize dropping the source's in-plane law
    // when the resolved field is connected to the original duct owner.
    let bare = spatial_source
        .inner()
        .prepare_time_with_thermal_fields(
            1e-7,
            spatial_source.inner().medium().cp_j_kg_k,
            &[],
            ThermoviscousWallDiscretization {
                normal_elements: &[2],
                maximum_frequency_hz: BAND[1],
            },
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let gate = fs_exec::CancelGate::new_clock_free();
    assert!(
        super::super::dg::ThermoviscousDg::new(
            &spatial_source,
            bare,
            time::time_dg_load(&gas, 1e-7, fs_duct::dg::DgBoundary::Closed),
            1,
            super::super::dg::ThermoviscousDgLimits::default(),
            &gate
        )
        .is_err()
    );
    assert!(
        super::super::dg::ThermoviscousDg::new(
            &spatial_source,
            spatial,
            time::time_dg_load(&gas, 1e-7, fs_duct::dg::DgBoundary::Closed),
            1,
            super::super::dg::ThermoviscousDgLimits::default(),
            &gate
        )
        .is_ok()
    );
}
