use super::*;
use crate::{HoleState, Segment};
use fs_material::gas::{ConductivityModel, GasSpec};
use fs_math::det;

fn gases(pressure: f64) -> (GasState, GasState) {
    let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, pressure).unwrap();
    let heavy = GasState::try_new(
        &GasSpec {
            molar_mass: 0.044,
            gamma: 1.3,
            sutherland_beta: 1.5e-6,
            sutherland_s: 240.0,
            conductivity: ConductivityModel::Eucken,
        },
        300.0,
        pressure,
    )
    .unwrap();
    (air, heavy)
}

fn duct(state: HoleState) -> Duct {
    Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.002,
                length: 0.09,
            },
            Segment::ToneHole {
                hole_radius: 0.0006,
                chimney_height: 0.003,
                bore_radius: 0.002,
                state,
            },
            Segment::Cylinder {
                radius: 0.002,
                length: 0.15,
            },
        ],
    }
}

fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs().max(1e-30)
}
fn bits(r: DuctResponse) -> [u64; 9] {
    [
        r.omega,
        r.impedance.re,
        r.impedance.im,
        r.min_shear_number,
        r.mouth_ka,
        r.u_mouth_over_u_in.re,
        r.u_mouth_over_u_in.im,
        r.p_mouth_over_u_in.re,
        r.p_mouth_over_u_in.im,
    ]
    .map(f64::to_bits)
}

#[test]
fn g0_identical_media_keep_response_and_every_hole_flow_bitwise() {
    let (gas, _) = gases(100000.0);
    for state in [HoleState::Closed, HoleState::Open, HoleState::Vent(0.35)] {
        for loss in [LossModel::Lossless, LossModel::Bessel] {
            for termination in [Termination::Closed, Termination::FlangedOpen] {
                for frequency in [231.0, 739.0, 1493.0] {
                    let omega = core::f64::consts::TAU * frequency;
                    let original = crate::input_impedance_with_wall_and_holes(
                        &duct(state),
                        &gas,
                        omega,
                        loss,
                        termination,
                        WallModel::Rigid,
                    )
                    .unwrap();
                    let explicit = input_impedance_and_holes(
                        &duct(state),
                        &gas,
                        &gas,
                        omega,
                        loss,
                        termination,
                        WallModel::Rigid,
                    )
                    .unwrap();
                    assert_eq!(bits(original.0), bits(explicit.0));
                    assert_eq!(original.1, explicit.1);
                }
            }
        }
    }
}

// Independent lossless p/U chain under exp(-i omega t), including the finite
// chimney. Radiation stays with its existing owner; this checks how its load
// crosses media, not a second copy of the radiation model.
fn multiply(a: [C64; 4], b: [C64; 4]) -> [C64; 4] {
    [
        a[0] * b[0] + a[1] * b[2],
        a[0] * b[1] + a[1] * b[3],
        a[2] * b[0] + a[3] * b[2],
        a[2] * b[1] + a[3] * b[3],
    ]
}
fn cylinder(gas: &GasState, radius: f64, length: f64, omega: f64) -> [C64; 4] {
    let phase = omega * length / gas.sound_speed;
    let z = gas.characteristic_impedance / (core::f64::consts::PI * radius * radius);
    let c = C64::from_re(det::cos(phase));
    [
        c,
        C64::new(0.0, -z * det::sin(phase)),
        C64::new(0.0, -det::sin(phase) / z),
        c,
    ]
}

#[test]
fn g1_exterior_chimney_matches_independent_finite_line_and_flow_continuity() {
    let (inside, outside) = gases(100000.0);
    for sigma in [0.35, 1.0] {
        for frequency in [231.0, 739.0, 1493.0] {
            let omega = core::f64::consts::TAU * frequency;
            let state = HoleState::Vent(sigma);
            let (actual, holes) = input_impedance_and_holes(
                &duct(state),
                &inside,
                &outside,
                omega,
                LossModel::Lossless,
                Termination::Closed,
                WallModel::Rigid,
            )
            .unwrap();
            let (a, b, h) = (0.002, 0.0006, 0.003);
            let radiation =
                crate::compact_radiation_impedance(Termination::FlangedOpen, &outside, b, omega)
                    .unwrap();
            let chimney = cylinder(&inside, b, h + fs_phs::side_hole_inner_length(b, a), omega);
            let open =
                (chimney[0] * radiation + chimney[1]) / (chimney[2] * radiation + chimney[3]);
            let capped = cylinder(&inside, b, h, omega);
            let closed = capped[0] / capped[2];
            let shunt_y = open.recip().scale(sigma) + closed.recip().scale(1.0 - sigma);
            let half = C64::new(
                0.0,
                -omega * inside.density * fs_phs::side_hole_series_length(b, a) * sigma
                    / (2.0 * core::f64::consts::PI * a * a),
            );
            let series = [C64::ONE, half, C64::ZERO, C64::ONE];
            let junction = multiply(
                multiply(series, [C64::ONE, C64::ZERO, shunt_y, C64::ONE]),
                series,
            );
            let prefix = cylinder(&inside, a, 0.09, omega);
            let whole = multiply(
                multiply(prefix, junction),
                cylinder(&inside, a, 0.15, omega),
            );
            let zin = whole[0] / whole[2];
            assert!(relative(actual.impedance, zin) < 2e-10);
            let before_p = prefix[3] * zin - prefix[1];
            let before_u = -prefix[2] * zin + prefix[0];
            let midpoint_p = before_p - half * before_u;
            let chimney_in = midpoint_p.scale(sigma) / open;
            let chimney_out = chimney_in / (chimney[2] * radiation + chimney[3]);
            assert_eq!(holes.len(), 1);
            assert!(relative(holes[0].u_chimney_in_over_u_in, chimney_in) < 2e-10);
            assert!(relative(holes[0].u_hole_over_u_in, chimney_out) < 2e-10);
            assert!(actual.impedance.re >= 0.0);
            let original = crate::input_impedance(
                &duct(state),
                &inside,
                omega,
                LossModel::Lossless,
                Termination::Closed,
            )
            .unwrap();
            assert!(
                relative(actual.impedance, original.impedance) > 1e-4,
                "a closed bell must not discard radiation from open side holes"
            );
        }
    }
}

#[test]
fn g3_pressure_scaling_and_closed_hole_isolation() {
    let (inside, outside) = gases(100000.0);
    let (dense_inside, dense_outside) = gases(200000.0);
    let omega = core::f64::consts::TAU * 739.0;
    let solve = |a: &GasState, b: &GasState| {
        input_impedance_and_holes(
            &duct(HoleState::Open),
            a,
            b,
            omega,
            LossModel::Lossless,
            Termination::Closed,
            WallModel::Rigid,
        )
        .unwrap()
    };
    let a = solve(&inside, &outside);
    let b = solve(&dense_inside, &dense_outside);
    assert!(relative(b.0.impedance, a.0.impedance.scale(2.0)) < 2e-12);
    assert!(relative(b.1[0].u_hole_over_u_in, a.1[0].u_hole_over_u_in) < 2e-12);
    let closed = |gas| {
        input_impedance(
            &duct(HoleState::Closed),
            &inside,
            gas,
            omega,
            LossModel::Bessel,
            Termination::Closed,
            WallModel::Rigid,
        )
        .unwrap()
    };
    assert_eq!(bits(closed(&inside)), bits(closed(&outside)));
    assert!(
        input_impedance(
            &duct(HoleState::Open),
            &inside,
            &dense_outside,
            omega,
            LossModel::Bessel,
            Termination::Closed,
            WallModel::Rigid
        )
        .is_err()
    );
    let mut bad = outside;
    bad.sound_speed *= 2.0;
    assert!(
        input_impedance(
            &duct(HoleState::Open),
            &inside,
            &bad,
            omega,
            LossModel::Bessel,
            Termination::Closed,
            WallModel::Rigid
        )
        .is_err()
    );
}

#[test]
fn g0_independent_exterior_preserves_finite_reflection_at_a_closed_duct_pole() {
    let (inside, outside) = gases(100000.0);
    let (radius, length) = (0.002, 0.2);
    let geometry = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let omega = core::f64::consts::PI * inside.sound_speed / length;
    let zc = inside.characteristic_impedance / (core::f64::consts::PI * radius * radius);
    let reflection = input_reflection(
        &geometry,
        &inside,
        &outside,
        omega,
        LossModel::Lossless,
        Termination::Closed,
        WallModel::Rigid,
    )
    .unwrap()
    .referenced_to(zc)
    .unwrap();
    assert!((reflection - C64::ONE).abs() < 1e-12);
}
