use super::*;
use crate::acoustics::helmholtz::boundary_layer::BoundaryLayerGeometry;
use crate::acoustics::helmholtz::boundary_layer::thermal_wall::{
    ConductiveThermalWall, ThermalBacking, ThermalWallDepth, WallThermalTransport,
};
use crate::acoustics::thermoviscous::wall::ConductiveWallRegion;

fn solid() -> ConductiveThermalWall {
    ConductiveThermalWall {
        transport: WallThermalTransport {
            conductivity_w_m_k: 0.003,
            volumetric_heat_capacity_j_m3_k: 2000.0,
        },
        depth: ThermalWallDepth::Slab {
            thickness_m: 0.00002,
            backing: ThermalBacking::Adiabatic,
        },
        minimum_tangential_scale_m: 0.004,
    }
}

fn region(faces: Vec<[u32; 3]>, wall: ConductiveThermalWall) -> ConductiveWallRegion {
    ConductiveWallRegion {
        faces,
        wall,
        geometry: BoundaryLayerGeometry {
            minimum_curvature_radius_m: f64::INFINITY,
            minimum_half_gap_m: 0.0001,
            maximum_layer_ratio: 0.05,
        },
        temperature_k: air().temperature_k,
    }
}

fn relative(actual: C64, expected: C64) -> f64 {
    (actual - expected).abs() / expected.abs().max(1e-30)
}

#[test]
fn finite_solid_regions_preserve_full_lns_work_reciprocity_and_loaded_scattering() {
    let size = [0.004, 0.0008, 0.001];
    let (mesh, pos) = box_mesh([2, 2, 2], size);
    let boundary = boundaries(&mesh, &pos, size, false);
    let mut gas = air();
    gas.bulk_viscosity_pa_s = 1.2e-5;
    let (lower, upper): (Vec<_>, Vec<_>) = boundary
        .no_slip_faces
        .iter()
        .copied()
        .partition(|f| f.iter().all(|&v| pos[v as usize][1] < 0.5 * size[1]));
    let mut upper_solid = solid();
    upper_solid.depth = ThermalWallDepth::HalfSpace;
    upper_solid.transport.conductivity_w_m_k *= 2.5;
    let regions = [region(lower, solid()), region(upper, upper_solid)];
    let model = ThermoviscousFe::new(&mesh, &pos, gas, &boundary, limits(), &mut || false)
        .unwrap()
        .with_conductive_walls(&regions, &mut || false)
        .unwrap();
    assert_eq!(model.conductive_walls(), regions);
    let cp = C64::new(970.0, -45.0);
    let frequency = 1700.0;
    let prepared = model
        .prepare_with_heat_capacity(frequency, cp, &mut || false)
        .unwrap();
    let drive = [C64::new(1.0, 0.2), C64::new(0.7, -0.1)];
    let response = prepared.response(&drive, &mut || false).unwrap();
    assert_work(&response);
    let p = response.power;
    assert!(p.solid_conductive_w > 0.0 && p.conductive_w > 0.0 && p.caloric_w > 0.0);
    assert!(p.viscous_w > 0.0 && p.bulk_w > 0.0);
    let repeat = prepared.response(&drive, &mut || false).unwrap();
    assert_eq!(response.temperature_k, repeat.temperature_k);
    assert_eq!(response.velocity_m_s, repeat.velocity_m_s);
    let y = model
        .admittance_with_heat_capacity(frequency, cp, &mut || false)
        .unwrap();
    assert!(relative(y[1], y[2]) < 1e-9);
    let refs = [5e8, 3e8];
    let s = model
        .scattering_with_heat_capacity(frequency, cp, &refs, &mut || false)
        .unwrap();
    let a = [C64::new(0.4, 0.1), C64::new(-0.3, 0.7)];
    let b = [s[0] * a[0] + s[1] * a[1], s[2] * a[0] + s[3] * a[1]];
    let pressures = std::array::from_fn::<_, 2, _>(|i| (a[i] + b[i]).scale(det::sqrt(refs[i])));
    let loaded = prepared.response(&pressures, &mut || false).unwrap();
    assert_work(&loaded);
    let wave_power = 0.5
        * (a.iter().map(|z| z.norm_sq()).sum::<f64>() - b.iter().map(|z| z.norm_sq()).sum::<f64>());
    assert!((wave_power - loaded.power.dissipated_w()).abs() < 1e-9 * wave_power);
    for i in 0..2 {
        assert!(
            relative(
                loaded.inward_flow_m3_s[i],
                (a[i] - b[i]).scale(1.0 / det::sqrt(refs[i]))
            ) < 1e-9
        );
    }
    // Restoring insulation reuses exactly the original gas operator and fields.
    let restored = model.with_conductive_walls(&[], &mut || false).unwrap();
    let original =
        ThermoviscousFe::new(&mesh, &pos, gas, &boundary, limits(), &mut || false).unwrap();
    let r = restored
        .prepare_with_heat_capacity(frequency, cp, &mut || false)
        .unwrap()
        .response(&drive, &mut || false)
        .unwrap();
    let o = original
        .prepare_with_heat_capacity(frequency, cp, &mut || false)
        .unwrap()
        .response(&drive, &mut || false)
        .unwrap();
    assert_eq!(r.temperature_k, o.temperature_k);
    assert_eq!(r.inward_flow_m3_s, o.inward_flow_m3_s);
    assert_eq!(r.power.solid_conductive_w, 0.0);
    assert!(relative(response.inward_flow_m3_s[0], o.inward_flow_m3_s[0]) > 1e-4);
}

fn tanh(z: C64) -> C64 {
    let denominator = (2.0 * z.re).cosh() + (2.0 * z.im).cos();
    C64::new(
        (2.0 * z.re).sinh() / denominator,
        (2.0 * z.im).sin() / denominator,
    )
}

#[test]
fn finite_solid_full_lns_heat_mesh_converges_to_independent_slab_compliance() {
    // Short acoustic length makes pressure uniform to the separately retained
    // finite-length/stress defect. The gas heat field varies across y; z walls
    // and both traction mouths are insulated. This is a low-kL thermal slab
    // benchmark, not infinite-slit propagation or an instrument acceptance.
    let size = [0.0001, 0.0003, 0.0004];
    let gas = air();
    let frequency = 800.0;
    let omega = std::f64::consts::TAU * frequency;
    let cp = C64::new(970.0, -45.0);
    let q = (C64::new(
        0.0,
        omega * gas.density_kg_m3 / gas.thermal_conductivity_w_m_k,
    ) * cp)
        .sqrt();
    let qa = q.scale(size[1] / 2.0);
    for backing in [ThermalBacking::Isothermal, ThermalBacking::Adiabatic] {
        let mut wall = solid();
        wall.depth = ThermalWallDepth::Slab {
            thickness_m: 0.00002,
            backing,
        };
        let qw = C64::new(
            0.0,
            omega * wall.transport.volumetric_heat_capacity_j_m3_k
                / wall.transport.conductivity_w_m_k,
        )
        .sqrt();
        // Direct complex analytic solution, independent of the production
        // solid-admittance branch/series and of the FE trace/volume matrices.
        let y = qw.scale(wall.transport.conductivity_w_m_k)
            * match backing {
                ThermalBacking::Adiabatic => tanh(qw.scale(0.00002)),
                ThermalBacking::Isothermal => tanh(qw.scale(0.00002)).recip(),
            };
        let deficit =
            (tanh(qa) / qa) * y / (y + q.scale(gas.thermal_conductivity_w_m_k) * tanh(qa));
        let beta = C64::new(1.0 / gas.pressure_pa, 0.0)
            - (C64::ONE - deficit) / cp.scale(gas.density_kg_m3 * gas.temperature_k);
        let exact = C64::new(0.0, omega * size.iter().product::<f64>()) * beta;
        let mut errors = Vec::new();
        for across in [2, 4, 8] {
            let (mesh, pos) = box_mesh([1, across, 2], size);
            let boundary = boundaries(&mesh, &pos, size, false);
            let faces = boundary
                .no_slip_faces
                .iter()
                .copied()
                .filter(|f| {
                    f.iter().all(|&v| pos[v as usize][1] == 0.0)
                        || f.iter().all(|&v| pos[v as usize][1] == size[1])
                })
                .collect();
            let model = ThermoviscousFe::new(&mesh, &pos, gas, &boundary, limits(), &mut || false)
                .unwrap()
                .with_conductive_walls(&[region(faces, wall)], &mut || false)
                .unwrap();
            let field = model
                .prepare_with_heat_capacity(frequency, cp, &mut || false)
                .unwrap()
                .response(&[C64::ONE, C64::ONE], &mut || false)
                .unwrap();
            assert_work(&field);
            assert!(field.power.solid_conductive_w > 0.0 && field.power.conductive_w > 0.0);
            let error = relative(field.inward_flow_m3_s[0] + field.inward_flow_m3_s[1], exact);
            eprintln!(
                "finite-solid full LNS backing={backing:?} across={across} dofs={} error={error:e}",
                model.ndof()
            );
            errors.push(error);
        }
        assert!(errors[1] < errors[0] && errors[2] < errors[1], "{errors:?}");
        assert!(errors[2] < 2e-4, "{errors:?}");
    }
}

#[test]
fn finite_solid_regions_refuse_conflicting_faces_temperature_profile_and_cancellation() {
    let size = [0.002, 0.001, 0.001];
    let (mesh, pos) = box_mesh([1, 2, 2], size);
    let boundary = boundaries(&mesh, &pos, size, false);
    let model =
        || ThermoviscousFe::new(&mesh, &pos, air(), &boundary, limits(), &mut || false).unwrap();
    let good = region(boundary.no_slip_faces.clone(), solid());
    assert!(matches!(
        model().with_conductive_walls(&[good.clone()], &mut || true),
        Err(ThermoviscousError::Cancelled)
    ));
    let mut bad = good.clone();
    bad.temperature_k *= 1.0 + 1e-12;
    assert!(
        model()
            .with_conductive_walls(&[bad], &mut || false)
            .is_err()
    );
    let mut bad = good.clone();
    bad.faces.push(bad.faces[0]);
    assert!(
        model()
            .with_conductive_walls(&[bad], &mut || false)
            .is_err()
    );
    let mut bad = good.clone();
    bad.faces = vec![boundary.ports[0].traces[0].0];
    assert!(
        model()
            .with_conductive_walls(&[bad], &mut || false)
            .is_err()
    );
    assert!(
        model()
            .with_conductive_walls(&[good.clone(), good.clone()], &mut || false)
            .is_err()
    );
    let mut isothermal = boundary.clone();
    isothermal.isothermal_faces = boundary.no_slip_faces.clone();
    assert!(
        ThermoviscousFe::new(&mesh, &pos, air(), &isothermal, limits(), &mut || false)
            .unwrap()
            .with_conductive_walls(&[good.clone()], &mut || false)
            .is_err()
    );
    let mut bad = good.clone();
    bad.wall.transport.conductivity_w_m_k = -1.0;
    assert!(
        model()
            .with_conductive_walls(&[bad], &mut || false)
            .unwrap()
            .prepare(800.0, &mut || false)
            .is_err()
    );
    let mut bad = good.clone();
    bad.geometry.minimum_curvature_radius_m = 1e-9;
    assert!(
        model()
            .with_conductive_walls(&[bad], &mut || false)
            .unwrap()
            .prepare(800.0, &mut || false)
            .is_err()
    );
    let model = model()
        .with_conductive_walls(&[good], &mut || false)
        .unwrap();
    assert!(model.prepare(800.0, &mut || true).is_err());
    let prepared = model.prepare(800.0, &mut || false).unwrap();
    assert!(
        prepared
            .response(&[C64::ONE, C64::ZERO], &mut || true)
            .is_err()
    );
}
