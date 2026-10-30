//! Independent finite-line/mutual-load closure and physical-port work checks.
use super::*;
use crate::HoleState;
use fs_material::gas::{ConductivityModel, GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation,
};
use fs_thermochem::SpeciesId;

fn state(relaxing: bool) -> RelaxationState {
    let gas = GasState::try_new(
        &GasSpec {
            molar_mass: 0.0289644,
            gamma: 1.25,
            sutherland_beta: 1.458e-6,
            sutherland_s: 110.4,
            conductivity: ConductivityModel::Eucken,
        },
        293.15,
        101325.0,
    )
    .unwrap();
    let species = SpeciesId::new("A").unwrap();
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes: if relaxing {
            vec![RelaxationProcess {
                reservoir: InternalReservoir::Rotational {
                    degrees_of_freedom: 2,
                },
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 101325.0 * 0.002,
                    t_ref_k: 293.15,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 400.0,
                },
            }]
        } else {
            vec![]
        },
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap())
    .unwrap()
}
fn model(segments: Vec<Segment>, relaxing: bool, wall: bool) -> MolecularDuct {
    MolecularDuct::new(
        &Duct { segments },
        state(relaxing),
        if wall {
            Some(MolecularWallResolution {
                poles_per_field: 32,
                band_rad_s: [10.0, 8000.0],
                maximum_relative_error: 0.2,
            })
        } else {
            None
        },
        ConeSlices::Auto,
        256,
    )
    .unwrap()
}
fn hole(opening: HoleState, height: f64) -> Segment {
    Segment::ToneHole {
        hole_radius: 0.002,
        chimney_height: height,
        bore_radius: 0.008,
        state: opening,
    }
}
fn cylinder(length: f64) -> Segment {
    Segment::Cylinder {
        radius: 0.008,
        length,
    }
}
fn labels(n: usize) -> Vec<String> {
    (0..n)
        .map(|i| format!("mouth-{i}"))
        .chain(["terminal".to_string()])
        .collect()
}
fn load(omega: f64, labels: &[String], z: Vec<C64>) -> CoupledTermination {
    CoupledTermination {
        omega,
        labels: labels.to_vec(),
        z,
        source_id: "explicit-analytic-test-matrix".into(),
    }
}
#[track_caller]
fn near(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-16),
        "{a:?} != {b:?}"
    );
}
fn line(gas: &GasState, omega: f64, radius: f64, length: f64, relaxing: bool) -> [C64; 4] {
    // Independent one-reservoir rational bulk law: cv_eq=4R, cv_frozen=3R.
    let t = omega * 0.002;
    let bulk = if relaxing {
        C64::new(
            gas.pressure * (20.0 + 12.0 * t * t) / (16.0 + 9.0 * t * t),
            -gas.pressure * t / (16.0 + 9.0 * t * t),
        )
    } else {
        C64::from_re(1.25 * gas.pressure)
    };
    let k = (C64::from_re(gas.density) / bulk).sqrt().scale(omega);
    let z = C64::from_re(omega * gas.density / (core::f64::consts::PI * radius * radius)) / k;
    let phase = k.scale(length);
    let (s, c) = det::sin_cos(phase.re);
    let e = det::exp(phase.im);
    let ch = 0.5 * (e + 1.0 / e);
    let sh = 0.5 * (e - 1.0 / e);
    let sin = C64::new(s * ch, c * sh);
    let cos = C64::new(c * ch, -s * sh);
    [
        cos,
        C64::new(0.0, -1.0) * z * sin,
        C64::new(0.0, -1.0) * sin / z,
        cos,
    ]
}

#[test]
fn diagonal_matrix_recovers_scalar_chimneys_and_physical_vent_work() {
    let geometry = vec![
        cylinder(0.1),
        hole(HoleState::Open, 0.03),
        cylinder(0.08),
        hole(HoleState::Closed, 0.02),
        cylinder(0.05),
        hole(HoleState::Vent(0.27), 0.04),
        cylinder(0.13),
    ];
    for wall in [false, true] {
        let model = model(geometry.clone(), true, wall);
        let labels = labels(2);
        for omega in [20.0, 400.0, 2500.0, 6500.0] {
            let r = model.reference_impedance;
            let mouths = [C64::new(2.0 * r, -0.5 * r), C64::new(1.2 * r, -0.3 * r)];
            let terminal = C64::new(0.8 * r, -0.15 * r);
            let scalar = model
                .response_with_holes(omega, Some(terminal), &mouths, ComposeOptions::default())
                .unwrap();
            let mut z = vec![C64::ZERO; 9];
            z[0] = mouths[0];
            z[4] = mouths[1].scale(1.0 / 0.27);
            z[8] = terminal;
            let prepared = model
                .prepare_exterior_ports(omega, &labels, ComposeOptions::default())
                .unwrap();
            let response = prepared.response(&load(omega, &labels, z)).unwrap();
            near(response.z_in, scalar.duct.z_in, 4e-10);
            near(
                response.exterior_flows_over_u_in[0],
                scalar.holes[0].u_mouth_over_u_in,
                4e-10,
            );
            near(
                response.exterior_flows_over_u_in[1],
                scalar.holes[2].u_mouth_over_u_in,
                4e-10,
            );
            near(
                response.exterior_flows_over_u_in[2],
                scalar.duct.u_out_over_u_in,
                4e-10,
            );
            near(
                response.exterior_pressures_over_u_in[1],
                scalar.holes[2].p_mouth_over_u_in,
                4e-10,
            );
            let work = scalar
                .holes
                .iter()
                .map(|h| h.mouth_power_per_squared_inlet_flow)
                .sum::<f64>()
                + 0.5 * terminal.re * scalar.duct.u_out_over_u_in.norm_sq();
            assert!((response.exterior_power_per_squared_inlet_flow - work).abs() < 4e-10 * work);
            assert!(response.dissipated_power_per_squared_inlet_flow > 0.0);
        }
    }
}

#[test]
fn mutual_radiation_matches_independent_pressure_flow_closure() {
    for sigma in [1.0, 0.31] {
        let model = model(vec![hole(HoleState::Vent(sigma), 0.035)], true, false);
        let labels = labels(1);
        let r = model.reference_impedance;
        let z = [
            C64::new(2.0 * r, -0.7 * r),
            C64::new(0.35 * r, -0.2 * r),
            C64::new(0.35 * r, -0.2 * r),
            C64::new(0.8 * r, -0.3 * r),
        ];
        for omega in [2.0, 200.0, 1200.0, 7000.0] {
            let prepared = model
                .prepare_exterior_ports(omega, &labels, ComposeOptions::default())
                .unwrap();
            let response = prepared
                .response(&load(omega, &labels, z.to_vec()))
                .unwrap();
            let gas = model.state().gas();
            let o = line(
                gas,
                omega,
                0.002,
                0.035 + fs_phs::side_hole_inner_length(0.002, 0.008),
                true,
            );
            let c = line(gas, omega, 0.002, 0.035, true);
            let yc = c[2] / c[0];
            let half = C64::new(
                0.0,
                -0.5 * omega * gas.density * fs_phs::side_hole_series_length(0.002, 0.008) * sigma
                    / (core::f64::consts::PI * 0.008 * 0.008),
            );
            // Unknowns: physical outward mouth flow q, terminal flow u.
            // Junction pressure from chimney equals terminal pressure + series.
            let p_q = o[0] * z[0] + o[1].scale(1.0 / sigma);
            let p_u = o[0] * z[1];
            let a = p_q - z[2];
            let b = p_u - z[3] - half;
            let c = (o[2] * z[0]).scale(sigma) + o[3] + (yc * p_q).scale(1.0 - sigma);
            let d = C64::ONE + (o[2] * z[1]).scale(sigma) + (yc * p_u).scale(1.0 - sigma);
            let determinant = a * d - b * c;
            let q = -b / determinant;
            let u = a / determinant;
            let junction = p_q * q + p_u * u;
            near(response.z_in, junction + half, 3e-10);
            near(response.exterior_flows_over_u_in[0], q, 3e-10);
            near(response.exterior_flows_over_u_in[1], u, 3e-10);
            near(
                response.exterior_pressures_over_u_in[0],
                z[0] * q + z[1] * u,
                3e-10,
            );
            near(
                response.exterior_pressures_over_u_in[1],
                z[2] * q + z[3] * u,
                3e-10,
            );
            assert!(response.dissipated_power_per_squared_inlet_flow > 0.0);
            let diagonal = prepared
                .response(&load(
                    omega,
                    &labels,
                    vec![z[0], C64::ZERO, C64::ZERO, z[3]],
                ))
                .unwrap();
            assert!((response.z_in - diagonal.z_in).abs() > 1e-3 * response.z_in.abs());
            // A permutation of matrix identities, rows and columns is exact physics.
            let reverse = load(
                omega,
                &[labels[1].clone(), labels[0].clone()],
                vec![z[3], z[2], z[1], z[0]],
            );
            near(
                response.z_in,
                prepared.response(&reverse).unwrap().z_in,
                1e-13,
            );
        }
    }
}

#[test]
fn exposed_mouth_operator_is_reciprocal_and_lossless_work_closes() {
    let model = model(
        vec![
            cylinder(0.1),
            hole(HoleState::Vent(0.4), 0.03),
            cylinder(0.07),
            hole(HoleState::Closed, 0.025),
            cylinder(0.1),
            hole(HoleState::Open, 0.02),
            cylinder(0.12),
        ],
        false,
        false,
    );
    let labels = labels(2);
    for omega in [20.0, 1000.0, 5000.0] {
        let prepared = model
            .prepare_exterior_ports(omega, &labels, ComposeOptions::default())
            .unwrap();
        let s = prepared.scattering();
        for i in 0..4 {
            for j in 0..4 {
                assert!((s[i * 4 + j] - s[j * 4 + i]).abs() < 2e-12);
                let mut gram = C64::ZERO;
                for k in 0..4 {
                    gram = gram + s[k * 4 + i].conj() * s[k * 4 + j];
                }
                assert!((gram - if i == j { C64::ONE } else { C64::ZERO }).abs() < 3e-12);
            }
        }
        let r = model.reference_impedance;
        let mut z = vec![C64::new(0.2 * r, -0.1 * r); 9];
        for i in 0..3 {
            z[i * 3 + i] = z[i * 3 + i] + C64::new(r, -0.2 * r);
        }
        let response = prepared.response(&load(omega, &labels, z)).unwrap();
        assert!(response.exterior_power_per_squared_inlet_flow > 0.0);
        assert!(
            response.dissipated_power_per_squared_inlet_flow.abs() < 5e-11 * response.z_in.abs()
        );
    }
}

#[test]
fn rank_deficient_load_sealed_pole_and_opaque_chimney_are_admissible() {
    let gas = state(false);
    let omega = gas.gas().sound_speed * core::f64::consts::FRAC_PI_2 / 0.03;
    let closed = model(vec![hole(HoleState::Closed, 0.03)], false, false);
    let labels = labels(0);
    let prepared = closed
        .prepare_exterior_ports(omega, &labels, ComposeOptions::default())
        .unwrap();
    let response = prepared
        .response(&load(
            omega,
            &labels,
            vec![C64::from_re(closed.reference_impedance)],
        ))
        .unwrap();
    assert!(response.z_in.abs() < 1e-10 * closed.reference_impedance);
    let opaque = model(
        vec![cylinder(0.1), hole(HoleState::Open, 20000.0), cylinder(0.1)],
        true,
        false,
    );
    let labels = super::tests::labels(1);
    let prepared = opaque
        .prepare_exterior_ports(1000.0, &labels, ComposeOptions::default())
        .unwrap();
    let r = opaque.reference_impedance;
    // Rank-one positive semidefinite load; no inverse exists or is required.
    let response = prepared
        .response(&load(1000.0, &labels, vec![C64::from_re(r); 4]))
        .unwrap();
    assert!(response.z_in.re > 0.0);
    assert!(response.dissipated_power_per_squared_inlet_flow > 0.0);
    assert!(response.exterior_flows_over_u_in.iter().all(|&z| finite(z)));
    let pressure_release = prepared
        .response(&load(1000.0, &labels, vec![C64::ZERO; 4]))
        .unwrap();
    assert!(pressure_release.z_in.re > 0.0);
    assert!(
        pressure_release
            .exterior_pressures_over_u_in
            .iter()
            .all(|p| p.abs() < 1e-10 * r)
    );
}

#[test]
fn exterior_matrix_admission_caps_and_cancellation_are_explicit() {
    let model = model(
        vec![cylinder(0.1), hole(HoleState::Open, 0.03)],
        true,
        false,
    );
    let labels = labels(1);
    let options = ComposeOptions::default();
    let prepared = model
        .prepare_exterior_ports(1000.0, &labels, options)
        .unwrap();
    let r = model.reference_impedance;
    let good = load(
        1000.0,
        &labels,
        vec![C64::from_re(r), C64::ZERO, C64::ZERO, C64::from_re(r)],
    );
    assert!(prepared.response(&good).is_ok());
    let mut stale = good.clone();
    stale.omega += 1.0;
    assert!(prepared.response(&stale).is_err());
    let mut wrong = good.clone();
    wrong.labels[1] = wrong.labels[0].clone();
    assert!(prepared.response(&wrong).is_err());
    let active = load(
        1000.0,
        &labels,
        vec![
            C64::from_re(r),
            C64::from_re(2.0 * r),
            C64::from_re(2.0 * r),
            C64::from_re(r),
        ],
    );
    assert!(prepared.response(&active).is_err());
    assert!(
        model
            .prepare_exterior_ports_controlled(1000.0, &labels, options, &mut || true)
            .is_err()
    );
    assert!(prepared.response_controlled(&good, &mut || true).is_err());
    assert!(
        model
            .prepare_exterior_ports(1000.0, &["same".into(), "same".into()], options)
            .is_err()
    );
    let mut calls = 0;
    assert!(
        prepared
            .response_controlled(&good, &mut || {
                calls += 1;
                calls > 4
            })
            .is_err()
    );
    for count in [MAX_NETWORK_PORTS - 1, MAX_NETWORK_PORTS] {
        let geometry = (0..count)
            .flat_map(|_| [cylinder(0.01), hole(HoleState::Open, 0.01)])
            .collect();
        let many = super::tests::model(geometry, false, false);
        let labels = super::tests::labels(count);
        assert_eq!(
            many.prepare_exterior_ports(100.0, &labels, options).is_ok(),
            count < MAX_NETWORK_PORTS
        );
    }
}
