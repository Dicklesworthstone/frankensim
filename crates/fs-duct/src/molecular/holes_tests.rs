//! Finite-line and volume-work oracles. Synthetic rates isolate constitutive
//! coupling; these tests do not claim experimental validation of a gas or pad.
use super::*;
use crate::HoleState;
use fs_material::gas::{ConductivityModel, GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation,
};
use fs_thermochem::SpeciesId;

fn state(relaxing: bool, pressure: f64) -> RelaxationState {
    let gas = GasState::try_new(
        &GasSpec {
            molar_mass: 0.0289644,
            gamma: 1.25,
            sutherland_beta: 1.458e-6,
            sutherland_s: 110.4,
            conductivity: ConductivityModel::Eucken,
        },
        293.15,
        pressure,
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

fn geometry(hole_state: HoleState, height: f64) -> Duct {
    Duct {
        segments: vec![Segment::ToneHole {
            hole_radius: 0.002,
            chimney_height: height,
            bore_radius: 0.008,
            state: hole_state,
        }],
    }
}

fn model(
    geometry: &Duct,
    state: RelaxationState,
    wall: Option<MolecularWallResolution>,
) -> MolecularDuct {
    MolecularDuct::new(geometry, state, wall, ConeSlices::Auto, 128).unwrap()
}

#[track_caller]
fn near(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-20),
        "{a:?} != {b:?}"
    );
}

// Independent one-reservoir rational form: cv_eq=4R, cv_frozen=3R.
fn bulk(state: &RelaxationState, omega: f64) -> C64 {
    let p = state.gas().pressure;
    if state.processes().is_empty() {
        return C64::from_re(1.25 * p);
    }
    let y = omega * 0.002 * 101325.0 / p;
    C64::new(
        p * (20.0 + 12.0 * y * y) / (16.0 + 9.0 * y * y),
        -p * y / (16.0 + 9.0 * y * y),
    )
}

fn sin_cos(z: C64) -> (C64, C64) {
    let (s, c) = det::sin_cos(z.re);
    let e = det::exp(z.im);
    let ch = (e + 1.0 / e) * 0.5;
    let sh = (e - 1.0 / e) * 0.5;
    (C64::new(s * ch, c * sh), C64::new(c * ch, -s * sh))
}

fn line(state: &RelaxationState, omega: f64, radius: f64, length: f64) -> [C64; 4] {
    let area = core::f64::consts::PI * radius * radius;
    let k = (C64::from_re(state.gas().density) / bulk(state, omega))
        .sqrt()
        .scale(omega);
    let z = C64::from_re(omega * state.gas().density / area) / k;
    let (s, c) = sin_cos(k.scale(length));
    [
        c,
        C64::new(0.0, -1.0) * z * s,
        C64::new(0.0, -1.0) * s / z,
        c,
    ]
}

fn multiply(a: [C64; 4], b: [C64; 4]) -> [C64; 4] {
    [
        a[0] * b[0] + a[1] * b[2],
        a[0] * b[1] + a[1] * b[3],
        a[2] * b[0] + a[3] * b[2],
        a[2] * b[1] + a[3] * b[3],
    ]
}

fn integrated_heat(
    state: &RelaxationState,
    omega: f64,
    radius: f64,
    length: f64,
    p_end: C64,
    u_end: C64,
) -> f64 {
    let count = 256;
    let mut integral = 0.0;
    for i in 0..=count {
        let x = length * i as f64 / count as f64;
        let m = line(state, omega, radius, length - x);
        let p = m[0] * p_end + m[1] * u_end;
        let weight = if i == 0 || i == count {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight * p.norm_sq();
    }
    0.5 * omega
        * bulk(state, omega).recip().im
        * core::f64::consts::PI
        * radius
        * radius
        * integral
        * length
        / (3.0 * count as f64)
}

#[test]
fn molecular_open_closed_and_vent_chimneys_match_finite_lines_and_volume_work() {
    let state = state(true, 101325.0);
    let radius = 0.002;
    let height = 0.03;
    let length = height + fs_phs::side_hole_inner_length(radius, 0.008);
    for (hole_state, sigma) in [
        (HoleState::Open, 1.0),
        (HoleState::Closed, 0.0),
        (HoleState::Vent(0.31), 0.31),
    ] {
        let mut geometry = geometry(hole_state, height);
        geometry.segments.insert(
            0,
            Segment::Cylinder {
                radius: 0.008,
                length: 0.08,
            },
        );
        geometry.segments.push(Segment::Cylinder {
            radius: 0.008,
            length: 0.12,
        });
        let model = model(&geometry, state.clone(), None);
        assert_eq!(model.blocks(), 5);
        for omega in [2.0, 200.0, 800.0, 5000.0] {
            let reference = model.reference_impedance;
            let end = C64::new(reference * 0.9, -reference * 0.2);
            let mouth = C64::new(reference * 3.0, -reference * 0.7);
            let load = if sigma > 0.0 { vec![mouth] } else { vec![] };
            let result = model
                .response_with_holes(omega, Some(end), &load, ComposeOptions::default())
                .unwrap();
            let hole = result.holes[0];
            assert_eq!(hole.site.axial_position_m, 0.08);
            let o = line(&state, omega, radius, length);
            let c = line(&state, omega, radius, height);
            let zo = (o[0] * mouth + o[1]) / (o[2] * mouth + o[3]);
            let zc = c[0] / c[2];
            let y = zo.recip().scale(sigma) + zc.recip().scale(1.0 - sigma);
            let half = C64::new(
                0.0,
                -0.5 * omega
                    * state.gas().density
                    * fs_phs::side_hole_series_length(radius, 0.008)
                    * sigma
                    / (core::f64::consts::PI * 0.008 * 0.008),
            );
            let series = [C64::ONE, half, C64::ZERO, C64::ONE];
            let junction = multiply(multiply(series, [C64::ONE, C64::ZERO, y, C64::ONE]), series);
            let tail = line(&state, omega, 0.008, 0.12);
            let total = multiply(multiply(line(&state, omega, 0.008, 0.08), junction), tail);
            let u_end = (total[2] * end + total[3]).recip();
            let p_end = end * u_end;
            let p_j =
                tail[0] * p_end + tail[1] * u_end + half * (tail[2] * p_end + tail[3] * u_end);
            near(result.duct.z_in, (total[0] * end + total[1]) * u_end, 5e-12);
            near(result.duct.u_out_over_u_in, u_end, 5e-12);
            near(hole.p_junction_over_u_in, p_j, 5e-12);
            let open = p_j / zo;
            let closed = p_j / zc;
            if sigma > 0.0 {
                near(hole.u_open_in_over_u_in, open.scale(sigma), 5e-12);
                let u_mouth = open / (o[2] * mouth + o[3]);
                near(hole.u_mouth_over_u_in, u_mouth.scale(sigma), 5e-12);
                near(hole.p_mouth_over_u_in, mouth * u_mouth, 5e-12);
            } else {
                assert_eq!(hole.u_mouth_over_u_in, C64::ZERO);
                assert_eq!(hole.p_mouth_over_u_in, C64::ZERO);
            }
            if sigma < 1.0 {
                near(hole.u_closed_in_over_u_in, closed.scale(1.0 - sigma), 5e-12);
            }
            let u_mouth = open / (o[2] * mouth + o[3]);
            let expected_heat = sigma
                * integrated_heat(&state, omega, radius, length, mouth * u_mouth, u_mouth)
                + (1.0 - sigma)
                    * integrated_heat(&state, omega, radius, height, p_j / c[0], C64::ZERO);
            let scale = hole
                .input_power_per_squared_inlet_flow
                .abs()
                .max(expected_heat);
            assert!(
                (hole.dissipated_power_per_squared_inlet_flow - expected_heat).abs()
                    < 2e-10 * scale
            );
            assert!(hole.dissipated_power_per_squared_inlet_flow > 0.0);
            assert!(hole.mouth_power_per_squared_inlet_flow >= 0.0);
            near(
                result.duct.scattering.scattering[1],
                result.duct.scattering.scattering[2],
                5e-12,
            );
        }
    }
}

#[test]
fn sealed_chimney_poles_remain_finite_and_recover_the_compact_gas_spring() {
    let height = 0.04;
    let state = state(false, 101325.0);
    let model = model(&geometry(HoleState::Closed, height), state.clone(), None);
    for phase in [core::f64::consts::FRAC_PI_2, core::f64::consts::PI] {
        let omega = phase * state.gas().sound_speed / height;
        let result = model
            .response_with_holes(
                omega,
                Some(C64::from_re(model.reference_impedance)),
                &[],
                ComposeOptions::default(),
            )
            .unwrap();
        let flow = result.holes[0];
        assert!(finite(result.duct.z_in));
        assert!(result.duct.condition_estimate < 100.0);
        assert_eq!(flow.u_mouth_over_u_in, C64::ZERO);
        if phase == core::f64::consts::FRAC_PI_2 {
            assert!(result.duct.z_in.abs() < 1e-12 * model.reference_impedance);
            near(flow.u_closed_in_over_u_in, C64::ONE, 1e-12);
        } else {
            near(
                result.duct.z_in,
                C64::from_re(model.reference_impedance),
                1e-12,
            );
            assert!(flow.u_closed_in_over_u_in.abs() < 1e-12);
        }
    }
    let state = self::state(true, 101325.0);
    let model = self::model(&geometry(HoleState::Closed, height), state.clone(), None);
    let omega = 0.1;
    let result = model
        .response_with_holes(omega, None, &[], ComposeOptions::default())
        .unwrap();
    let volume = core::f64::consts::PI * 0.002 * 0.002 * height;
    let expected = C64::new(0.0, 1.0) * bulk(&state, omega).scale(1.0 / (omega * volume));
    near(result.duct.z_in, expected, 2e-10);
    // The stiff gas spring normalizes flow by subtracting nearly equal
    // incident/reflected waves. Its measured condition, not an arbitrary
    // frequency-independent tolerance, bounds this cancellation's roundoff.
    let roundoff = 16.0 * f64::EPSILON * result.duct.condition_estimate;
    assert!(result.duct.condition_estimate < 1e7);
    assert!((result.holes[0].u_closed_in_over_u_in - C64::ONE).abs() < roundoff);
}

#[test]
fn attenuating_molecular_chimney_stays_stable_and_absorbs_branch_work() {
    let model = model(
        &geometry(HoleState::Open, 100_000.0),
        state(true, 101325.0),
        None,
    );
    let reference = model.reference_impedance;
    let result = model
        .response_with_holes(
            10000.0,
            Some(C64::from_re(reference)),
            &[C64::from_re(reference)],
            ComposeOptions::default(),
        )
        .unwrap();
    let hole = result.holes[0];
    assert_eq!(hole.u_mouth_over_u_in, C64::ZERO);
    assert!(hole.dissipated_power_per_squared_inlet_flow > 0.0);
    let output = 0.5 * reference * result.duct.u_out_over_u_in.norm_sq();
    let input = 0.5 * result.duct.z_in.re;
    assert!((input - output - hole.input_power_per_squared_inlet_flow).abs() < 1e-12 * input);
}

#[test]
fn molecular_chimney_wall_resolution_refines_and_both_losses_are_live() {
    let shape = geometry(HoleState::Vent(0.45), 0.025);
    let law = state(true, 101325.0);
    let solve = |poles, law: RelaxationState| {
        let wall = Some(MolecularWallResolution {
            poles_per_field: poles,
            band_rad_s: [100.0, 6000.0],
            maximum_relative_error: 0.5,
        });
        let model = model(&shape, law, wall);
        model
            .response_with_holes(
                5000.0,
                Some(C64::from_re(model.reference_impedance)),
                &[C64::from_re(2e6)],
                ComposeOptions::default(),
            )
            .unwrap()
    };
    let reference = solve(32, law.clone());
    let mut errors = Vec::new();
    for poles in [4, 8, 16] {
        let response = solve(poles, law.clone());
        errors.push((response.duct.z_in - reference.duct.z_in).abs() / reference.duct.z_in.abs());
        assert!(response.holes[0].dissipated_power_per_squared_inlet_flow > 0.0);
    }
    assert!(
        errors[1] < errors[0] && errors[2] < errors[1] && errors[2] < 1e-6,
        "wall refinement: {errors:?}"
    );
    let frozen = solve(32, state(false, 101325.0));
    assert!((reference.duct.z_in - frozen.duct.z_in).abs() > 1e-3 * reference.duct.z_in.abs());
    let adiabatic = model(&shape, law, None);
    let adiabatic = adiabatic
        .response_with_holes(
            5000.0,
            Some(C64::from_re(adiabatic.reference_impedance)),
            &[C64::from_re(2e6)],
            ComposeOptions::default(),
        )
        .unwrap();
    assert!((reference.duct.z_in - adiabatic.duct.z_in).abs() > 1e-3 * reference.duct.z_in.abs());
}

#[test]
fn molecular_holes_require_explicit_loads_and_preserve_work_and_cancellation_caps() {
    let law = state(true, 101325.0);
    let shape = geometry(HoleState::Open, 0.01);
    assert!(MolecularDuct::new(&shape, law.clone(), None, ConeSlices::Auto, 2).is_err());
    let model = model(&shape, law.clone(), None);
    assert_eq!(model.blocks(), 3);
    assert_eq!(model.material_points(), 1);
    assert_eq!(model.open_hole_sites().len(), 1);
    for maximum_condition in [f64::NAN, 0.5] {
        let options = ComposeOptions {
            maximum_condition,
            ..ComposeOptions::default()
        };
        assert!(
            model
                .prepared_chain_with_hole_loads_controlled(
                    1000.0,
                    &[C64::ZERO],
                    options,
                    &mut || false,
                )
                .is_err()
        );
    }
    assert!(
        model
            .response(1000.0, None, ComposeOptions::default())
            .is_err()
    );
    for loads in [
        vec![],
        vec![C64::new(-1.0, 0.0)],
        vec![C64::new(f64::NAN, 0.0)],
        vec![C64::ZERO; 2],
    ] {
        assert!(
            model
                .response_with_holes(1000.0, None, &loads, ComposeOptions::default())
                .is_err()
        );
    }
    assert!(
        model
            .response_with_holes_controlled(
                1000.0,
                None,
                &[C64::ZERO],
                ComposeOptions::default(),
                &mut || true
            )
            .is_err()
    );
    for opening in [-0.1, 1.1, f64::NAN] {
        assert!(
            MolecularDuct::new(
                &geometry(HoleState::Vent(opening), 0.01),
                law.clone(),
                None,
                ConeSlices::Auto,
                3
            )
            .is_err()
        );
    }
}
