use super::*;
use fs_duct::compose::ComposeOptions;
use fs_duct::molecular::MolecularWallResolution;
use fs_duct::{ConeSlices, Duct, HoleState, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};
use fs_thermochem::SpeciesId;

fn gas(temperature_k: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature_k, 101325.0).unwrap()
}

fn point(rate_scale: f64, wall: bool, exterior: GasState) -> MolecularDrivingPoint {
    // Authored synthetic data exercise the original input and constitutive
    // owners; these are deliberately not asserted to be measured air constants.
    let data = format!(
        "declared 1 synthetic 2 vibration 1 1300 1 power {} 293.15 0 200 500 \
         vibration 1 1900 1 power {} 293.15 0 200 500",
        101325.0 * 2e-5 * rate_scale,
        101325.0 * 3.2e-4 * rate_scale,
    );
    let state = read_relaxation(&mut data.split_ascii_whitespace())
        .unwrap()
        .evaluate(
            &RelaxingMedium::new(
                gas(293.15),
                vec![(SpeciesId::new("synthetic").unwrap(), 1.0)],
            )
            .unwrap(),
        )
        .unwrap();
    MolecularDrivingPoint::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.001,
                length: 0.15,
            }],
        },
        state,
        wall.then_some(MolecularWallResolution {
            poles_per_field: 16,
            band_rad_s: [1e-4, 30000.0],
            maximum_relative_error: 0.05,
        }),
        ConeSlices::Fixed(8),
        ComposeOptions::default(),
        exterior,
    )
    .unwrap()
}

#[test]
fn molecular_outlet_matches_independent_loaded_cylinder_solution() {
    let model = point(1.0, false, gas(320.0));
    let radius = model.inlet_radius_m();
    let area = core::f64::consts::PI * radius * radius;
    let gas = model.inlet_gas();
    let zc = gas.density * gas.sound_speed / area;
    for omega in [100.0, 1000.0, 6000.0, 24000.0] {
        let k = model.material().state().relaxing_wavenumber(omega).unwrap();
        let phase = k.scale(0.15);
        let cosine = C64::new(
            phase.re.cos() * phase.im.cosh(),
            -phase.re.sin() * phase.im.sinh(),
        );
        let sine = C64::new(
            phase.re.sin() * phase.im.cosh(),
            phase.re.cos() * phase.im.sinh(),
        );
        let wave_z = C64::from_re(omega * gas.density / area) / k;
        let load = fs_duct::compact_radiation_impedance(
            Termination::FlangedOpen,
            model.exterior_gas(),
            radius,
            omega,
        )
        .unwrap();
        // Direct uniform-cylinder solution of [Pin,Uin] = ABCD [Zload,1] Uout
        // and p+ = (Pin + Zc Uin)/2, independent of chain/scattering composition.
        let mixed = wave_z + load.scale(zc) / wave_z;
        let denominator = (load + C64::from_re(zc)) * cosine + C64::new(mixed.im, -mixed.re) * sine;
        let expected = C64::from_re(2.0) / denominator;
        let actual = model.flanged_outlet_flow(omega).unwrap();
        assert!(
            (actual - expected).abs() < 2e-11 * expected.abs(),
            "omega={omega}: {actual:?} versus {expected:?}"
        );
    }
}

#[test]
fn changing_only_reservoir_data_changes_physical_pressure_after_the_frozen_wavefront() {
    let baseline = point(1.0, true, gas(320.0));
    let changed = point(0.1, true, gas(320.0));
    let receiver = CircularOutletReceiver {
        position_m: [0.0, 0.0, 1.0],
        radial_rings: 2,
        angular_points: 8,
        maximum_frequency_hz: 800.0,
    };
    let mut a = baseline.baffled_outlet(8000, 1024, receiver).unwrap();
    let mut b = changed.baffled_outlet(8000, 1024, receiver).unwrap();
    assert_eq!(a.info().gas, gas(320.0));
    assert_eq!(
        a.info().axial_delay_samples,
        (baseline.one_way_wavefront_s() * 8000.0).floor() as usize
    );
    assert_eq!(a.info().axial_delay_samples, b.info().axial_delay_samples);
    let first = a.info().axial_delay_samples + a.info().exterior_delay_samples.0;
    let (mut energy, mut difference) = (0.0, 0.0);
    for n in 0..2000 {
        let input = if n == 0 { 100.0 } else { 0.0 };
        let pa = a.step(input).unwrap();
        let pb = b.step(input).unwrap();
        if n < first {
            assert_eq!(pa, 0.0);
            assert_eq!(pb, 0.0);
        }
        energy += pa * pa;
        difference += (pa - pb).powi(2);
    }
    assert!(
        energy > 1e-12 && difference > 1e-6 * energy,
        "energy={energy}, difference={difference}"
    );
}

fn holed_point(scale: f64, opening: HoleState, wall: bool) -> MolecularDrivingPoint {
    let source = point(scale, false, gas(320.0));
    MolecularDrivingPoint::new(
        &Duct {
            segments: vec![
                Segment::Cylinder {
                    radius: 0.001,
                    length: 0.08,
                },
                Segment::ToneHole {
                    hole_radius: 0.0005,
                    chimney_height: 0.006,
                    bore_radius: 0.001,
                    state: opening,
                },
                Segment::Cylinder {
                    radius: 0.001,
                    length: 0.12,
                },
            ],
        },
        source.material().state().clone(),
        wall.then_some(MolecularWallResolution {
            poles_per_field: 16,
            band_rad_s: [1e-4, 30000.0],
            maximum_relative_error: 0.05,
        }),
        ConeSlices::Fixed(8),
        ComposeOptions::default(),
        gas(320.0),
    )
    .unwrap()
}

// Independent finite-cylinder ABCD solution, without the scattering/branch
// composition owner used in production. The original caloric law supplies k.
fn cylinder(model: &MolecularDrivingPoint, omega: f64, radius: f64, length: f64) -> [C64; 4] {
    let k = model.material().state().relaxing_wavenumber(omega).unwrap();
    let phase = k.scale(length);
    let c = C64::new(
        phase.re.cos() * phase.im.cosh(),
        -phase.re.sin() * phase.im.sinh(),
    );
    let s = C64::new(
        phase.re.sin() * phase.im.cosh(),
        phase.re.cos() * phase.im.sinh(),
    );
    let z =
        C64::from_re(omega * model.inlet_gas().density / (core::f64::consts::PI * radius * radius))
            / k;
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

#[track_caller]
fn near(a: C64, b: C64) {
    assert!(
        (a - b).abs() <= 2e-10 * a.abs().max(b.abs()).max(1e-20),
        "{a:?} != {b:?}"
    );
}

#[test]
fn molecular_openings_match_one_analytic_loaded_network_including_vent_weight() {
    for sigma in [1.0, 0.31] {
        let model = holed_point(1.0, HoleState::Vent(sigma), false);
        let site = model.open_hole_sites()[0];
        let zc = model.inlet_gas().characteristic_impedance
            / (core::f64::consts::PI * 0.001_f64.powi(2));
        for omega in [100.0, 1400.0, 6000.0] {
            let end = model
                .terminal_load(omega, Termination::FlangedOpen)
                .unwrap()
                .unwrap();
            let mouth = fs_duct::compact_radiation_impedance(
                Termination::FlangedOpen,
                model.exterior_gas(),
                site.hole_radius_m,
                omega,
            )
            .unwrap();
            let open = cylinder(
                &model,
                omega,
                site.hole_radius_m,
                site.open_chimney_length_m,
            );
            let closed = cylinder(&model, omega, site.hole_radius_m, site.chimney_height_m);
            let zo = (open[0] * mouth + open[1]) / (open[2] * mouth + open[3]);
            let y = zo.recip().scale(sigma) + (closed[2] / closed[0]).scale(1.0 - sigma);
            let half = C64::new(
                0.0,
                -0.5 * omega
                    * model.inlet_gas().density
                    * fs_phs::side_hole_series_length(site.hole_radius_m, site.bore_radius_m)
                    * sigma
                    / (core::f64::consts::PI * site.bore_radius_m.powi(2)),
            );
            let series = [C64::ONE, half, C64::ZERO, C64::ONE];
            let downstream = multiply(series, cylinder(&model, omega, 0.001, 0.12));
            let junction = multiply(series, [C64::ONE, C64::ZERO, y, C64::ONE]);
            let total = multiply(
                multiply(cylinder(&model, omega, 0.001, 0.08), junction),
                downstream,
            );
            let pressure = total[0] * end + total[1];
            let flow = total[2] * end + total[3];
            let terminal_flow = C64::from_re(2.0) / (pressure + flow.scale(zc));
            let junction_pressure = (downstream[0] * end + downstream[1]) * terminal_flow;
            let mouth_flow = junction_pressure / (open[0] * mouth + open[1]);
            let (actual_terminal, actual_holes) = model.flanged_opening_flows(omega).unwrap();
            assert_eq!(actual_holes.len(), 1);
            near(actual_terminal, terminal_flow);
            near(actual_holes[0], mouth_flow.scale(sigma));
            near(
                model
                    .reflectance(omega, Termination::FlangedOpen, zc)
                    .unwrap(),
                (pressure - flow.scale(zc)) / (pressure + flow.scale(zc)),
            );
        }
    }
}

#[test]
fn every_open_molecular_mouth_radiates_after_its_own_wavefront() {
    let baseline = holed_point(1.0, HoleState::Vent(0.31), true);
    let changed = holed_point(0.1, HoleState::Vent(0.31), true);
    let receiver = CircularOutletReceiver {
        position_m: [0.0, 0.0, 1.0],
        radial_rings: 2,
        angular_points: 8,
        maximum_frequency_hz: 800.0,
    };
    let (mut bell, mut holes) = baseline
        .baffled_outlet_and_open_holes_with_delay_speeds(8000, 1024, receiver, 340.0, &[320.0])
        .unwrap();
    let (mut changed_bell, mut changed_holes) = changed
        .baffled_outlet_and_open_holes_with_delay_speeds(8000, 1024, receiver, 340.0, &[320.0])
        .unwrap();
    assert_eq!(holes.len(), 1);
    let info = holes[0].info();
    assert_eq!(info.radius_m, 0.0005); // No second opening weight.
    assert_eq!(info.gas, gas(320.0));
    assert!((info.receiver.position_m[0] - 0.12).abs() < 1e-14);
    assert_eq!(info.receiver.position_m[1..], [0.0, 1.0]);
    let site = baseline.open_hole_sites()[0];
    assert_eq!(
        info.axial_delay_samples,
        ((site.axial_position_m + site.open_chimney_length_m)
            / baseline.material().state().frozen_sound_speed_m_s()
            * 8000.0)
            .floor() as usize
    );
    let first = info.axial_delay_samples + info.exterior_delay_samples.0;
    let (mut energy, mut difference) = (0.0, 0.0);
    let mut replay = holes[0].clone();
    for n in 0..1800 {
        let input = if n == 0 { 100.0 } else { 0.0 };
        let p = holes[0].step(input).unwrap();
        let q = changed_holes[0].step(input).unwrap();
        assert_eq!(p, replay.step(input).unwrap());
        if n < first {
            assert_eq!(p, 0.0);
        }
        energy += p * p;
        difference += (p - q).powi(2);
        assert!(bell.step(input).unwrap().is_finite());
        assert!(changed_bell.step(input).unwrap().is_finite());
    }
    assert!(
        energy > 1e-12 && difference > 1e-6 * energy,
        "energy={energy}, difference={difference}"
    );
    assert!(
        baseline
            .baffled_outlet_and_open_holes_with_delay_speeds(8000, 1024, receiver, 340.0, &[])
            .is_err()
    );
    let closed = holed_point(1.0, HoleState::Closed, true);
    let (_, holes) = closed
        .baffled_outlet_and_open_holes_delayed_at(8000, 1024, receiver, 340.0)
        .unwrap();
    assert!(holes.is_empty());
}

#[test]
fn closed_terminal_keeps_the_open_chimneys_static_pressure_release() {
    let open = holed_point(1.0, HoleState::Open, true);
    let sealed = holed_point(1.0, HoleState::Closed, true);
    let zc = open.inlet_gas().characteristic_impedance
        / (core::f64::consts::PI * open.inlet_radius_m().powi(2));
    let expected = open.reflectance(1e-3, Termination::Closed, zc).unwrap().re;
    assert!(expected < 0.0);
    assert!(
        sealed
            .reflectance(1e-3, Termination::Closed, zc)
            .unwrap()
            .re
            > 0.99
    );
    let (line, _) = open
        .realize_line(
            Termination::Closed,
            8000,
            zc,
            &crate::driving_point::LineLadder::default(),
        )
        .unwrap();
    let (_, taps) = line.quiescent_fir().unwrap();
    assert!((taps.iter().sum::<f64>() - expected).abs() < 2e-11);
}
