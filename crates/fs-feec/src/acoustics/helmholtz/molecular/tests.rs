use super::super::Medium;
use super::super::boundary_layer::thermal_wall::{
    ConductiveThermalWall, ThermalBacking, ThermalWallDepth, WallThermalTransport,
};
use super::*;
use crate::fixtures::{kuhn_cube, single_tet};

const CP: f64 = 1100.0;
const R: f64 = 287.0;
const TEMP: f64 = 300.0;
const RHO: f64 = 1.2;
const GEOMETRY: BoundaryLayerGeometry = BoundaryLayerGeometry {
    minimum_curvature_radius_m: f64::INFINITY,
    minimum_half_gap_m: 0.0005,
    maximum_layer_ratio: 0.05,
};

fn medium() -> Medium {
    Medium {
        density_kg_m3: RHO,
        sound_speed_m_s: (CP / (CP - R) * R * TEMP).sqrt(),
    }
}
fn gas(kappa: f64, mu: f64) -> CaloricGas {
    CaloricGas {
        pressure_pa: RHO * R * TEMP,
        temperature_k: TEMP,
        transport: BoundaryLayerTransport {
            dynamic_viscosity_pa_s: mu,
            thermal_conductivity_w_m_k: kappa,
            cp_j_kg_k: CP,
            heat_capacity_ratio: CP / (CP - R),
        },
    }
}
fn basis(fe: &HelmholtzFe<'_>, constant: bool) -> HelmholtzModes {
    let shapes = if constant {
        let mut p = vec![0.0; fe.ndof()];
        p[..fe.space.complex.vertex_count].fill(1.0);
        vec![p]
    } else {
        (0..fe.ndof())
            .map(|i| {
                let mut p = vec![0.0; fe.ndof()];
                p[i] = 1.0;
                p
            })
            .collect()
    };
    let n = shapes.len();
    HelmholtzModes {
        lambda: vec![0.0; n],
        frequency_hz: vec![0.0; n],
        residual: vec![0.0; n],
        shapes,
        below_low: 0,
        below_high: n,
    }
}
fn plane(fe: &HelmholtzFe<'_>, axis: usize, value: f64) -> Vec<[u32; 3]> {
    fe.boundary_faces()
        .into_iter()
        .filter(|face| {
            face.iter()
                .all(|&v| (fe.positions[v as usize][axis] - value).abs() < 1e-13)
        })
        .collect()
}
fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs()
}
fn exp(z: C64) -> C64 {
    C64::new(z.im.cos(), z.im.sin()).scale(z.re.exp())
}
fn tanh(z: C64) -> C64 {
    let e = exp(z.scale(2.0));
    (e - C64::ONE) / (e + C64::ONE)
}
fn assert_work(response: &CaloricResponse) {
    let input = response.acoustic.input_power_w();
    assert!(input > 0.0);
    assert!(
        (input - response.dissipation.total_w()).abs() < 2e-9 * input,
        "input={input:e}, loss={:?}",
        response.dissipation
    );
}

#[test]
fn original_volume_recovers_exact_molecular_compression_and_work() {
    let (mesh, points) = single_tet();
    let points: Vec<_> = points.iter().map(|p| p.map(|v| 0.1 * v)).collect();
    let fe = HelmholtzFe::new(&mesh, &points, 1, medium()).unwrap();
    let faces = fe.boundary_faces();
    let model = CaloricInterior::new(
        &fe,
        &basis(&fe, true),
        &[BoundaryMode::uniform(&faces[..1])],
        &faces[1..],
        gas(0.0, 0.0),
        ThermalWall::Adiabatic,
        GEOMETRY,
        &|| false,
    )
    .unwrap();
    let volume = 0.1_f64.powi(3) / 6.0;
    for f in [70.0, 700.0] {
        let omega = core::f64::consts::TAU * f;
        let cp = C64::from_re(900.0) + C64::new(1.0, omega * 0.001).recip().scale(200.0);
        let prepared = model.prepare(f, cp, &|| false).unwrap();
        let q = C64::new(1e-7, -0.3e-7);
        let response = prepared.response(&[q], &|| false).unwrap();
        let bulk = (cp / (cp - C64::from_re(R))).scale(gas(0.0, 0.0).pressure_pa);
        let exact = bulk * C64::new(0.0, -1.0 / (omega * volume));
        assert!(relative(response.acoustic.pressure_pa[0] / q, exact) < 2e-13);
        let temperature = response.acoustic.pressure_pa[0] / cp.scale(RHO);
        assert!(
            response
                .temperature_k
                .iter()
                .all(|&t| relative(t, temperature) < 2e-13)
        );
        assert_work(&response);
        assert_eq!(response.dissipation.gas_conduction_w, 0.0);
        assert_eq!(response.dissipation.solid_conduction_w, 0.0);
        assert_eq!(response.dissipation.viscous_wall_w, 0.0);
    }
}

#[test]
fn resolved_joint_heat_field_converges_to_independent_slab_with_molecular_cp() {
    let length = 0.001;
    let volume = length * 0.002 * 0.002;
    let frequency = 80.0;
    let omega = core::f64::consts::TAU * frequency;
    let cp = C64::new(1000.0, -100.0);
    let conductivity = 0.03;
    let q = (C64::new(0.0, omega * RHO / conductivity) * cp).sqrt();
    let qa = q.scale(length / 2.0);
    let wall = ConductiveThermalWall {
        transport: WallThermalTransport {
            conductivity_w_m_k: 0.01,
            volumetric_heat_capacity_j_m3_k: 1000.0,
        },
        depth: ThermalWallDepth::Slab {
            thickness_m: 0.00008,
            backing: ThermalBacking::Adiabatic,
        },
        minimum_tangential_scale_m: f64::INFINITY,
    };
    for boundary in [ThermalWall::Isothermal, ThermalWall::Conductive(wall)] {
        // Independent complex tanh slab solution, not the wall owner's series.
        let deficit = match boundary {
            ThermalWall::Isothermal => tanh(qa) / qa,
            _ => {
                let qw = C64::new(0.0, omega * 1000.0 / 0.01).sqrt();
                let y = qw.scale(0.01) * tanh(qw.scale(0.00008));
                (tanh(qa) / qa) * y / (y + q.scale(conductivity) * tanh(qa))
            }
        };
        let beta = C64::from_re(1.0 / gas(conductivity, 0.0).pressure_pa)
            - (C64::ONE - deficit) / (cp.scale(RHO * TEMP));
        let exact = (C64::new(0.0, omega * volume) * beta).recip();
        let mut errors = Vec::new();
        for n in [1, 2, 3] {
            let (mesh, points) = kuhn_cube(n);
            let points: Vec<_> = points
                .iter()
                .map(|p| [p[0] * length, p[1] * 0.002, p[2] * 0.002])
                .collect();
            let fe = HelmholtzFe::new(&mesh, &points, 2, medium()).unwrap();
            let mut walls = plane(&fe, 0, 0.0);
            walls.extend(plane(&fe, 0, length));
            let ports = [BoundaryMode::uniform(&plane(&fe, 1, 0.0))];
            let model = CaloricInterior::new(
                &fe,
                &basis(&fe, true),
                &ports,
                &walls,
                gas(conductivity, 0.0),
                boundary,
                GEOMETRY,
                &|| false,
            )
            .unwrap();
            let response = model
                .prepare(frequency, cp, &|| false)
                .unwrap()
                .response(&[C64::from_re(1e-10)], &|| false)
                .unwrap();
            assert_work(&response);
            assert!(
                response.dissipation.molecular_w > 0.0
                    && response.dissipation.gas_conduction_w > 0.0
            );
            if matches!(boundary, ThermalWall::Conductive(_)) {
                assert!(response.dissipation.solid_conduction_w > 0.0);
            }
            errors.push(relative(
                response.acoustic.pressure_pa[0].scale(1e10),
                exact,
            ));
        }
        println!("joint molecular thermal slab {boundary:?}: {errors:?}");
        assert!(errors.windows(2).all(|p| p[1] < p[0]));
        assert!(
            errors[2] < 0.002,
            "thermal field requires mesh refinement: {errors:?}"
        );
    }
}

#[test]
fn spatial_heat_and_viscous_work_close_and_power_scattering_is_reciprocal() {
    let (mesh, points) = kuhn_cube(1);
    let points: Vec<_> = points.iter().map(|p| p.map(|v| 0.05 * v)).collect();
    let fe = HelmholtzFe::new(&mesh, &points, 2, medium()).unwrap();
    let left = plane(&fe, 0, 0.0);
    let right = plane(&fe, 0, 0.05);
    let walls: Vec<_> = fe
        .boundary_faces()
        .into_iter()
        .filter(|f| !left.contains(f) && !right.contains(f))
        .collect();
    let ports = [
        BoundaryMode::uniform(&left),
        BoundaryMode {
            traces: right
                .iter()
                .map(|&f| (f, f.map(|v| 0.4 + fe.positions[v as usize][1] / 0.05)))
                .collect(),
        },
    ];
    let geometry = BoundaryLayerGeometry {
        minimum_half_gap_m: 0.025,
        ..GEOMETRY
    };
    let model = CaloricInterior::new(
        &fe,
        &basis(&fe, false),
        &ports,
        &walls,
        gas(0.03, 1.8e-5),
        ThermalWall::Isothermal,
        geometry,
        &|| false,
    )
    .unwrap();
    let cp = C64::new(1000.0, -80.0);
    let prepared = model.prepare(700.0, cp, &|| false).unwrap();
    let flows = [C64::new(1e-6, -0.1e-6), C64::new(-0.3e-6, 0.7e-6)];
    let response = prepared.response(&flows, &|| false).unwrap();
    assert_work(&response);
    assert!(response.dissipation.viscous_wall_w > 0.0);
    let replay = prepared.response(&flows, &|| false).unwrap();
    assert_eq!(response.acoustic.pressure_pa, replay.acoustic.pressure_pa);
    assert_eq!(response.temperature_k, replay.temperature_k);
    let s = model.scattering(700.0, cp, &[1e6, 1e6], &|| false).unwrap();
    assert!((s[1] - s[2]).abs() < 2e-12);
    for incident in [
        [C64::ONE, C64::ZERO],
        [C64::new(0.2, 0.8), C64::new(-0.1, 0.4)],
    ] {
        let reflected = [
            s[0] * incident[0] + s[1] * incident[1],
            s[2] * incident[0] + s[3] * incident[1],
        ];
        assert!(
            reflected.iter().map(|v| v.norm_sq()).sum::<f64>()
                < incident.iter().map(|v| v.norm_sq()).sum::<f64>()
        );
    }
}

#[test]
fn active_calorics_wrong_equilibrium_and_cancellation_refuse() {
    let (mesh, points) = single_tet();
    let fe = HelmholtzFe::new(&mesh, &points, 1, medium()).unwrap();
    let faces = fe.boundary_faces();
    let ports = [BoundaryMode::uniform(&faces[..1])];
    let model = CaloricInterior::new(
        &fe,
        &basis(&fe, true),
        &ports,
        &faces[1..],
        gas(0.03, 0.0),
        ThermalWall::Adiabatic,
        GEOMETRY,
        &|| false,
    )
    .unwrap();
    for cp in [
        C64::new(1000.0, 1.0),
        C64::from_re(R),
        C64::new(f64::NAN, 0.0),
    ] {
        assert!(model.prepare(200.0, cp, &|| false).is_err());
    }
    assert!(
        model
            .prepare(200.0, C64::new(1000.0, -10.0), &|| true)
            .is_err()
    );
    let mut bad = gas(0.03, 0.0);
    bad.pressure_pa *= 1.001;
    assert!(
        CaloricInterior::new(
            &fe,
            &basis(&fe, true),
            &ports,
            &faces[1..],
            bad,
            ThermalWall::Adiabatic,
            GEOMETRY,
            &|| false
        )
        .is_err()
    );
    let prepared = model
        .prepare(200.0, C64::new(1000.0, -10.0), &|| false)
        .unwrap();
    assert!(prepared.response(&[C64::ONE], &|| true).is_err());
    assert!(prepared.response(&[], &|| false).is_err());
}
