//! Independent DD analytic ABCD oracles; production never calls this code.
use super::*;
use fs_material::gas::GasSpec;
use fs_math::dd::Dd;
use std::{
    f64::consts::{PI, TAU},
    time::Instant,
};
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn opts(rule: Composition) -> ComposeOptions {
    ComposeOptions {
        rule,
        preserve_well_conditioned_bits: false,
        ..ComposeOptions::default()
    }
}
fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs().max(f64::MIN_POSITIVE)
}
fn emit(
    case: &str,
    fixture: &str,
    rung: usize,
    old: Option<f64>,
    error: f64,
    condition: f64,
    tol: f64,
    start: Instant,
    extra: &str,
) {
    let old = old
        .filter(|v| v.is_finite())
        .map_or_else(|| "null".into(), |v| format!("{v:.17e}"));
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"{case}\",\"fixture\":\"{fixture}\",\"inputs\":{{\"fixture\":\"{fixture}\",\"provenance\":\"analytic SI data and independent double-double ABCD\"}},\"axis\":\"composition\",\"rung\":{rung},\"observed\":{error:.17e},\"err\":{error:.17e},\"oracle\":\"closed form or DD ABCD\",\"abcd_rel_err\":{old},\"stable_rel_err\":{error:.17e},\"cond\":{condition:.17e},\"tol\":{tol:.17e},\"verdict\":\"{}\",\"elapsed_us\":{}{extra}}}",
        if error <= tol { "pass" } else { "fail" },
        start.elapsed().as_micros()
    );
}
#[derive(Clone, Copy)]
struct Dc {
    r: Dd,
    i: Dd,
}
impl Dc {
    const ZERO: Self = Self {
        r: Dd::ZERO,
        i: Dd::ZERO,
    };
    const ONE: Self = Self {
        r: Dd::ONE,
        i: Dd::ZERO,
    };
    fn from(z: C64) -> Self {
        Self {
            r: Dd::from_f64(z.re),
            i: Dd::from_f64(z.im),
        }
    }
    fn real(x: f64) -> Self {
        Self::from(C64::from_re(x))
    }
    fn scale(self, x: f64) -> Self {
        Self {
            r: self.r * Dd::from_f64(x),
            i: self.i * Dd::from_f64(x),
        }
    }
    fn rounded(self) -> C64 {
        C64::new(self.r.to_f64(), self.i.to_f64())
    }
    fn exp(self) -> Self {
        let mut n = 0;
        let mut z = self;
        while z.r.hi.abs() + z.i.hi.abs() > 0.25 {
            z = z.scale(0.5);
            n += 1;
        }
        let mut term = Self::ONE;
        let mut sum = term;
        for j in 1..=100 {
            term = term * z / Self::real(f64::from(j));
            sum = sum + term;
            if term.r.hi.abs() + term.i.hi.abs() < 1e-34 {
                break;
            }
        }
        for _ in 0..n {
            sum = sum * sum;
        }
        sum
    }
}
impl std::ops::Add for Dc {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            r: self.r + b.r,
            i: self.i + b.i,
        }
    }
}
impl std::ops::Sub for Dc {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self {
            r: self.r - b.r,
            i: self.i - b.i,
        }
    }
}
impl std::ops::Mul for Dc {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self {
            r: self.r * b.r - self.i * b.i,
            i: self.r * b.i + self.i * b.r,
        }
    }
}
impl std::ops::Div for Dc {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        let d = b.r * b.r + b.i * b.i;
        Self {
            r: (self.r * b.r + self.i * b.i) / d,
            i: (self.i * b.r - self.r * b.i) / d,
        }
    }
}
fn dd_mul(a: [Dc; 4], b: [Dc; 4]) -> [Dc; 4] {
    [
        a[0] * b[0] + a[1] * b[2],
        a[0] * b[1] + a[1] * b[3],
        a[2] * b[0] + a[3] * b[2],
        a[2] * b[1] + a[3] * b[3],
    ]
}
fn dd_cylinder(wave: SegmentWave, length: f64) -> [Dc; 4] {
    let phase = Dc::from(C64::new(-wave.wavenumber.im, wave.wavenumber.re)).scale(length);
    let forward = phase.exp();
    let backward = (Dc::ZERO - phase).exp();
    let c = (forward + backward).scale(0.5);
    let minus_i_sin = (backward - forward).scale(0.5);
    let z = Dc::from(wave.characteristic_impedance);
    [c, z * minus_i_sin, minus_i_sin / z, c]
}
fn dd_matrix(duct: &Duct, gas: &GasState, w: f64, loss: LossModel) -> Result<[Dc; 4], DuctError> {
    let mut total = [Dc::ONE, Dc::ZERO, Dc::ZERO, Dc::ONE];
    for s in &duct.segments {
        let m = match *s {
            Segment::Cylinder { radius, length } => {
                dd_cylinder(crate::segment_wave(gas, radius, w, loss).unwrap(), length)
            }
            Segment::ToneHole {
                hole_radius,
                chimney_height,
                state: HoleState::Closed,
                ..
            } => {
                let loss = if loss == LossModel::WideTube {
                    LossModel::Bessel
                } else {
                    loss
                };
                let m = dd_cylinder(
                    crate::segment_wave(gas, hole_radius, w, loss).unwrap(),
                    chimney_height,
                );
                [Dc::ONE, Dc::ZERO, m[2] / m[0], Dc::ONE]
            }
            _ => {
                return Err(DuctError::BadParameter {
                    what: "independent DD fixture requires cylinders/closed chimneys",
                });
            }
        };
        total = dd_mul(total, m);
    }
    Ok(total)
}
fn dd_response(m: [Dc; 4], load: Option<C64>) -> (C64, C64, C64) {
    let (p, u) = load.map_or((Dc::ONE, Dc::ZERO), |z| (Dc::from(z), Dc::ONE));
    let d = m[2] * p + m[3] * u;
    (
        ((m[0] * p + m[1] * u) / d).rounded(),
        (u / d).rounded(),
        (p / d).rounded(),
    )
}

#[test]
fn two_medium_open_and_vented_holes_match_independent_junction_flow() {
    let local = gas();
    let exterior = GasState::try_new(
        &fs_material::gas_species::GasSpecies::CarbonDioxide.spec(), 303.15, local.pressure,
    ).unwrap();
    let (a, b, h) = (0.004, 0.0015, 0.006);
    let (left_length, right_length) = (0.21, 0.31);
    for sigma in [1.0, 0.25] {
        let duct = Duct { segments: vec![
            Segment::Cylinder { radius: a, length: left_length },
            Segment::ToneHole { hole_radius: b, chimney_height: h, bore_radius: a,
                state: HoleState::Vent(sigma) },
            Segment::Cylinder { radius: a, length: right_length },
        ] };
        for frequency in [173.0, 731.0, 2031.0] {
            let omega = TAU * frequency;
            let bore = crate::segment_wave(&local, a, omega, LossModel::Lossless).unwrap();
            let chimney = crate::segment_wave(&local, b, omega, LossModel::Lossless).unwrap();
            let load = crate::flanged_pipe_impedance_with_exterior(
                &local, &exterior, a, omega, crate::FLANGED_MODES,
            ).unwrap();
            let hole_load = crate::flanged_pipe_impedance_with_exterior(
                &local, &exterior, b, omega, crate::FLANGED_MODES,
            ).unwrap();
            let (z_open, chimney_flow, _) = dd_response(dd_cylinder(
                chimney, h + crate::side_hole_inner_length(b, a),
            ), Some(hole_load));
            let (z_closed, _, _) = dd_response(dd_cylinder(chimney, h), None);
            let admittance = z_open.recip().scale(sigma) + z_closed.recip().scale(1.0 - sigma);
            let half_series = Dc::from(C64::new(0.0,
                -omega * local.density * crate::side_hole_series_length(b, a) * sigma
                    / (2.0 * PI * a * a)));
            let series = [Dc::ONE, half_series, Dc::ZERO, Dc::ONE];
            let left = dd_mul(dd_cylinder(bore, left_length), series);
            let right = dd_mul(series, dd_cylinder(bore, right_length));
            let shunt = [Dc::ONE, Dc::ZERO, Dc::from(admittance), Dc::ONE];
            let (z_in, mouth_flow, mouth_pressure) =
                dd_response(dd_mul(dd_mul(left, shunt), right), Some(load));
            let middle_pressure = (right[0] * Dc::from(load) + right[1]) * Dc::from(mouth_flow);
            let entering_hole = middle_pressure.rounded() * z_open.recip().scale(sigma);
            let outward_hole = entering_hole * chimney_flow;

            let (response, holes) = crate::exterior::input_impedance_and_holes(
                &duct, &local, &exterior, omega, LossModel::Lossless,
                Termination::FlangedOpen, WallModel::Rigid,
            ).unwrap();
            assert_eq!(holes.len(), 1);
            for (actual, expected) in [
                (response.impedance, z_in),
                (response.u_mouth_over_u_in, mouth_flow),
                (response.p_mouth_over_u_in, mouth_pressure),
                (holes[0].u_chimney_in_over_u_in, entering_hole),
                (holes[0].u_hole_over_u_in, outward_hole),
            ] {
                assert!(relative(actual, expected) < 3.0e-11,
                    "sigma={sigma} frequency={frequency}: {actual:?} vs {expected:?}");
            }
            let reflection = crate::exterior::input_reflection(&duct, &local, &exterior,
                omega, LossModel::Lossless, Termination::FlangedOpen, WallModel::Rigid).unwrap();
            let zc = C64::from_re(reflection.reference_impedance);
            assert!(relative(reflection.reflection, (z_in - zc) / (z_in + zc)) < 3.0e-11);
            assert!(reflection.reflection.abs() <= 1.0 + 1.0e-12);

            let (legacy, legacy_holes) = crate::input_impedance_with_wall_and_holes(
                &duct, &local, omega, LossModel::Bessel, Termination::FlangedOpen, WallModel::Rigid,
            ).unwrap();
            let (same, same_holes) = crate::exterior::input_impedance_and_holes(
                &duct, &local, &local, omega, LossModel::Bessel,
                Termination::FlangedOpen, WallModel::Rigid,
            ).unwrap();
            assert_eq!(same.impedance, legacy.impedance);
            assert_eq!(same.u_mouth_over_u_in, legacy.u_mouth_over_u_in);
            assert_eq!(same_holes, legacy_holes);
        }
    }
}
// Measure today's actual primitive+ABCD operation graph, rather than rounding
// the independent oracle and calling that a measured legacy defect.
fn legacy_matrix(
    duct: &Duct,
    gas: &GasState,
    w: f64,
    loss: LossModel,
) -> Result<[C64; 4], DuctError> {
    let mut m = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
    for s in &duct.segments {
        let block = match *s {
            Segment::Cylinder { radius, .. } => {
                crate::segment_matrix(s, &crate::segment_wave(gas, radius, w, loss)?)?
            }
            Segment::ToneHole {
                hole_radius,
                chimney_height,
                bore_radius,
                state: HoleState::Closed,
            } => crate::tone_hole_t_junction(
                gas,
                hole_radius,
                chimney_height,
                HoleState::Closed,
                w,
                loss,
                bore_radius,
                0.0,
                None,
            )?,
            _ => {
                return Err(DuctError::BadParameter {
                    what: "legacy measurement fixture primitive",
                });
            }
        };
        m = crate::mul2(m, block);
    }
    Ok(m)
}
fn stable(
    duct: &Duct,
    gas: &GasState,
    w: f64,
    loss: LossModel,
    load: Option<C64>,
    rule: Composition,
) -> ChainResult {
    let o = opts(rule);
    let chain = duct_chain(duct, gas, w, loss, None, o, true, || false).unwrap();
    let r = compose_chain(&chain.blocks, load, o).unwrap();
    assert!(
        [r.z_in, r.u_out_over_u_in, r.p_out_over_u_in]
            .iter()
            .all(|&z| finite(z))
    );
    assert!(r.condition_estimate.is_finite());
    r
}
fn old_error(
    duct: &Duct,
    gas: &GasState,
    w: f64,
    loss: LossModel,
    load: Option<C64>,
    oracle: C64,
) -> Option<f64> {
    crate::input_impedance_legacy(
        duct,
        gas,
        w,
        loss,
        load.map_or(LoadSpec::Analytic(Termination::Closed), LoadSpec::Fixed),
        None,
    )
    .ok()
    .map(|r| relative(r.impedance, oracle))
}
fn difficult() -> Vec<(&'static str, Duct, f64)> {
    let gas = gas();
    let mut horn = Vec::new();
    let length = gas.sound_speed / (4.0 * 18000.0);
    // Declared increasing envelope with 25 narrow/wide stop-band cells.
    for i in 0..50 {
        horn.push(Segment::Cylinder {
            radius: (if i % 2 == 0 { 0.001 } else { 0.012 }) * (1.0 + 0.01 * i as f64),
            length,
        });
    }
    let mut holes = Vec::new();
    let bore = 0.004;
    let length = gas.sound_speed / (4.0 * 15000.0);
    for _ in 0..24 {
        holes.push(Segment::Cylinder {
            radius: bore,
            length,
        });
        holes.push(Segment::ToneHole {
            hole_radius: 0.003,
            chimney_height: length * 1.0001,
            bore_radius: bore,
            state: HoleState::Closed,
        });
    }
    holes.push(Segment::Cylinder {
        radius: bore,
        length,
    });
    vec![
        (
            "a-10m-1mm-Bessel-20kHz",
            Duct {
                segments: vec![Segment::Cylinder {
                    radius: 0.001,
                    length: 10.0,
                }],
            },
            20000.0,
        ),
        (
            "b-50-segment-staircase-18kHz",
            Duct { segments: horn },
            18000.0,
        ),
        (
            "c-24-closed-hole-lattice-15kHz",
            Duct { segments: holes },
            15000.0,
        ),
    ]
}
#[test]
fn stable_composition_ill_conditioned_cases() {
    let gas = gas();
    for (rung, (name, duct, hz)) in difficult().into_iter().enumerate() {
        let start = Instant::now();
        let w = TAU * hz;
        let z = crate::segment_wave(
            &gas,
            duct.segments.last().unwrap().outlet_radius(),
            w,
            LossModel::Bessel,
        )
        .unwrap()
        .characteristic_impedance;
        let m = dd_matrix(&duct, &gas, w, LossModel::Bessel).unwrap();
        let want = dd_response(m, Some(z));
        let actual = stable(
            &duct,
            &gas,
            w,
            LossModel::Bessel,
            Some(z),
            Composition::ImpedanceRecursion,
        );
        let error = relative(actual.z_in, want.0);
        let transfer_error =
            relative(actual.u_out_over_u_in, want.1).max(relative(actual.p_out_over_u_in, want.2));
        let raw = legacy_matrix(&duct, &gas, w, LossModel::Bessel);
        let det_error = raw
            .as_ref()
            .ok()
            .map(|m| (m[0] * m[3] - m[1] * m[2] - C64::ONE).abs())
            .filter(|v| v.is_finite());
        let raw_status = match &raw {
            Err(_) => "refused",
            Ok(m) if m.iter().any(|&v| !finite(v)) => "overflow",
            _ => "finite",
        };
        emit(
            "compose-ill-conditioned",
            name,
            rung,
            old_error(&duct, &gas, w, LossModel::Bessel, Some(z), want.0),
            error,
            actual.condition_estimate,
            1e-12,
            start,
            &format!(
                ",\"abcd_determinant_error\":{},\"abcd_status\":\"{raw_status}\",\"transfer_rel_err\":{transfer_error:.17e}",
                det_error.map_or_else(|| "null".into(), |v| format!("{v:.17e}"))
            ),
        );
        assert!(
            error <= 1e-12 && transfer_error <= 1e-12,
            "{name}: Z={error}, transfer={transfer_error}"
        );
        let default =
            crate::input_impedance_load(&duct, &gas, w, LossModel::Bessel, z, None).unwrap();
        assert_eq!(default.composition, CompositionMode::ImpedanceRecursion);
        assert!(relative(default.impedance, want.0) <= 1e-12);
        if rung == 0 {
            let wave = crate::segment_wave(&gas, 0.001, w, LossModel::Bessel).unwrap();
            assert!(wave.wavenumber.im * 10.0 > 10.0);
            assert!(relative(actual.z_in, wave.characteristic_impedance) <= 1e-12);
        }
    }
    for (rung, length) in [10.0, 100.0, 1000.0].into_iter().enumerate() {
        let start = Instant::now();
        let w = TAU * 20000.0;
        let duct = Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.001,
                length,
            }],
        };
        let wave = crate::segment_wave(&gas, 0.001, w, LossModel::Bessel).unwrap();
        let actual = stable(
            &duct,
            &gas,
            w,
            LossModel::Bessel,
            None,
            Composition::ImpedanceRecursion,
        );
        let error = relative(actual.z_in, wave.characteristic_impedance);
        assert!(error <= 1e-12);
        emit(
            "compose-attenuation-length",
            "1mm-Bessel-20kHz",
            rung,
            old_error(
                &duct,
                &gas,
                w,
                LossModel::Bessel,
                None,
                wave.characteristic_impedance,
            ),
            error,
            actual.condition_estimate,
            1e-12,
            start,
            &format!(",\"length_m\":{length}"),
        );
    }
}
fn familiar() -> Vec<Duct> {
    let mut holes = vec![Segment::Cylinder {
        radius: 0.006,
        length: 0.12,
    }];
    for state in [HoleState::Closed, HoleState::Open, HoleState::Vent(0.4)] {
        holes.push(Segment::ToneHole {
            hole_radius: 0.003,
            chimney_height: 0.004,
            bore_radius: 0.006,
            state,
        });
        holes.push(Segment::Cylinder {
            radius: 0.006,
            length: 0.12,
        });
    }
    let profile = std::sync::Arc::new(crate::RadiusProfile::exponential(0.006, 2.0, 0.2).unwrap());
    vec![
        Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.006,
                length: 0.6,
            }],
        },
        Duct {
            segments: vec![Segment::Cone {
                inlet_radius: 0.006,
                outlet_radius: 0.006 * (1.0 + 1e-13),
                length: 0.6,
            }],
        },
        Duct {
            segments: vec![Segment::Cone {
                inlet_radius: 0.004,
                outlet_radius: 0.008,
                length: 0.3,
            }],
        },
        Duct {
            segments: vec![Segment::Cone {
                inlet_radius: 0.008,
                outlet_radius: 0.004,
                length: 0.3,
            }],
        },
        Duct { segments: holes },
        Duct {
            segments: vec![Segment::Profile {
                profile,
                steps: crate::RiccatiSteps(64),
            }],
        },
    ]
}
#[test]
fn stable_composition_matches_abcd_when_well_conditioned() {
    let gas = gas();
    for (index, duct) in familiar().iter().enumerate() {
        for loss in [LossModel::Lossless, LossModel::WideTube, LossModel::Bessel] {
            for hz in [110.0, 440.0, 880.0] {
                let w = TAU * hz;
                let z = C64::new(1.7e6, -2.1e5);
                let start = Instant::now();
                let old =
                    crate::input_impedance_legacy(duct, &gas, w, loss, LoadSpec::Fixed(z), None)
                        .unwrap();
                let a = stable(
                    duct,
                    &gas,
                    w,
                    loss,
                    Some(z),
                    Composition::ImpedanceRecursion,
                );
                let error = relative(a.z_in, old.impedance)
                    .max(relative(a.u_out_over_u_in, old.u_mouth_over_u_in))
                    .max(relative(a.p_out_over_u_in, old.p_mouth_over_u_in));
                emit(
                    "compose-well-conditioned",
                    &format!("existing-fixture-{index}-{loss:?}-{hz}Hz"),
                    index,
                    Some(0.0),
                    error,
                    a.condition_estimate,
                    1e-13,
                    start,
                    "",
                );
                assert!(
                    error <= 1e-13,
                    "fixture={index} loss={loss:?} hz={hz}: {error}"
                );
                let default = crate::input_impedance_load(duct, &gas, w, loss, z, None).unwrap();
                assert_eq!(default.composition, CompositionMode::LegacyAbcd);
                for (a, b) in [
                    (default.impedance, old.impedance),
                    (default.u_mouth_over_u_in, old.u_mouth_over_u_in),
                    (default.p_mouth_over_u_in, old.p_mouth_over_u_in),
                ] {
                    assert_eq!(
                        (a.re.to_bits(), a.im.to_bits()),
                        (b.re.to_bits(), b.im.to_bits())
                    );
                }
            }
        }
    }
}
#[test]
fn stable_composition_transfer_functions() {
    let gas = gas();
    let start = Instant::now();
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for count in [1, 8, 64] {
            let (radius, length, w) = (0.004, 0.73, TAU * 731.0);
            let duct = Duct {
                segments: vec![
                    Segment::Cylinder {
                        radius,
                        length: length / count as f64
                    };
                    count
                ],
            };
            let wave = crate::segment_wave(&gas, radius, w, loss).unwrap();
            let matrix = dd_cylinder(wave, length);
            let reference = gas.density * gas.sound_speed / (PI * radius * radius);
            for load in [
                Some(C64::ZERO),
                Some(C64::new(1.2e6, -0.4e6)),
                None,
                Some(C64::from_re(-reference)),
            ] {
                let want = dd_response(matrix, load);
                let a = stable(&duct, &gas, w, loss, load, Composition::ImpedanceRecursion);
                let uerr = if load.is_none() {
                    assert_eq!(a.u_out_over_u_in, C64::ZERO);
                    0.0
                } else {
                    relative(a.u_out_over_u_in, want.1)
                };
                let perr = if load == Some(C64::ZERO) {
                    assert_eq!(a.p_out_over_u_in, C64::ZERO);
                    0.0
                } else {
                    relative(a.p_out_over_u_in, want.2)
                };
                let error = relative(a.z_in, want.0).max(uerr).max(perr);
                assert!(
                    error <= 1e-12,
                    "{loss:?} count={count} load={load:?}: {error}"
                );
                emit(
                    "compose-cylinder-transfers",
                    "cascaded-exact-cylinder-SI",
                    count,
                    None,
                    error,
                    a.condition_estimate,
                    1e-12,
                    start,
                    "",
                );
            }
        }
    }
}
#[test]
fn scattering_rule_equals_impedance_recursion() {
    let gas = gas();
    let start = Instant::now();
    for (name, duct, hz) in difficult().into_iter().chain(
        familiar()
            .into_iter()
            .map(|d| ("existing-fixture", d, 440.0)),
    ) {
        let w = TAU * hz;
        let load = Some(C64::from_re(1e6));
        let a = stable(
            &duct,
            &gas,
            w,
            LossModel::Bessel,
            load,
            Composition::ImpedanceRecursion,
        );
        let b = stable(
            &duct,
            &gas,
            w,
            LossModel::Bessel,
            load,
            Composition::Scattering,
        );
        let error = relative(a.z_in, b.z_in)
            .max(relative(a.u_out_over_u_in, b.u_out_over_u_in))
            .max(relative(a.p_out_over_u_in, b.p_out_over_u_in));
        assert!(error <= 1e-12, "{name}: {error}");
        emit(
            "compose-two-rules",
            name,
            0,
            None,
            error,
            a.condition_estimate,
            1e-12,
            start,
            "",
        );
    }
}
#[test]
fn condition_estimate_refuses_above_bound() {
    let segment = Segment::Cylinder {
        radius: 1.0,
        length: PI + 1e-6,
    };
    let wave = SegmentWave {
        wavenumber: C64::ONE,
        characteristic_impedance: C64::ONE,
        specific_impedance: C64::from_re(PI),
        shear_number: f64::INFINITY,
    };
    let b = TwoPort::from_wave(&segment, &wave, 1.0).unwrap();
    let o = opts(Composition::ImpedanceRecursion);
    let c = compose_chain(&[b], None, o).unwrap().condition_estimate;
    assert!(c > 1e5);
    assert!(
        compose_chain(
            &[b],
            None,
            ComposeOptions {
                maximum_condition: c,
                ..o
            }
        )
        .is_ok()
    );
    assert!(matches!(
        compose_chain(
            &[b],
            None,
            ComposeOptions {
                maximum_condition: c.next_down(),
                ..o
            }
        ),
        Err(DuctError::Singular)
    ));
    let through = TwoPort::through(1.0).unwrap();
    let blocks = vec![through; MAX_COMPOSITION_BLOCKS + 1];
    assert!(compose_chain(&blocks[..MAX_COMPOSITION_BLOCKS], Some(C64::ONE), o).is_ok());
    assert!(compose_chain(&blocks, Some(C64::ONE), o).is_err());
    assert!(compose_chain(&[b], Some(C64::new(f64::NAN, 0.0)), o).is_err());
    assert!(
        TwoPort::from_wave(
            &segment,
            &SegmentWave {
                wavenumber: C64::new(1.0, -1e-12),
                ..wave
            },
            1.0
        )
        .is_err()
    );
    assert!(compose_chain_controlled(&[b], None, o, &mut || true).is_err());
    let mut calls = 0;
    assert!(
        compose_chain_controlled(&[b, b], None, o, &mut || {
            calls += 1;
            calls > 1
        })
        .is_err()
    );
    emit(
        "compose-condition-bound",
        "near-singular-closed-cylinder",
        0,
        None,
        0.0,
        c,
        0.0,
        Instant::now(),
        "",
    );
}
#[test]
fn stable_composition_bitwise_repeat() {
    let start = Instant::now();
    let gas = gas();
    for (_, duct, hz) in difficult() {
        let a = stable(
            &duct,
            &gas,
            TAU * hz,
            LossModel::Bessel,
            Some(C64::from_re(1e6)),
            Composition::ImpedanceRecursion,
        );
        for _ in 0..3 {
            let b = stable(
                &duct,
                &gas,
                TAU * hz,
                LossModel::Bessel,
                Some(C64::from_re(1e6)),
                Composition::ImpedanceRecursion,
            );
            for (a, b) in [
                (a.z_in, b.z_in),
                (a.u_out_over_u_in, b.u_out_over_u_in),
                (a.p_out_over_u_in, b.p_out_over_u_in),
            ] {
                assert_eq!(
                    (a.re.to_bits(), a.im.to_bits()),
                    (b.re.to_bits(), b.im.to_bits())
                );
            }
            assert_eq!(
                a.condition_estimate.to_bits(),
                b.condition_estimate.to_bits()
            );
        }
    }
    emit(
        "compose-bitwise-repeat",
        "three-ill-conditioned-fixtures",
        3,
        Some(0.0),
        0.0,
        1.0,
        0.0,
        start,
        "",
    );
}

#[test]
fn wall_ports_use_the_same_local_wave_owner() {
    let gas = gas();
    let wall = WallPin {
        surface_density: 1.5,
        stiffness_per_area: 2e5,
        resistance: 100.0,
    };
    let load = C64::new(1.2e6, -0.3e6);
    let options = opts(Composition::ImpedanceRecursion);
    for duct in familiar().into_iter().take(5) {
        let omega = TAU * 220.0;
        let chain = duct_chain(
            &duct,
            &gas,
            omega,
            LossModel::Bessel,
            Some(&wall),
            options,
            true,
            || false,
        )
        .unwrap();
        let actual = compose_chain(&chain.blocks, Some(load), options).unwrap();
        let old = crate::input_impedance_legacy(
            &duct,
            &gas,
            omega,
            LossModel::Bessel,
            LoadSpec::Fixed(load),
            Some(&wall),
        )
        .unwrap();
        assert!(relative(actual.z_in, old.impedance) <= 1e-13);
        assert!(relative(actual.u_out_over_u_in, old.u_mouth_over_u_in) <= 1e-13);
    }
}
#[test]
fn multiport_redheffer_and_modal_decay_share_the_owner() {
    let n = 2;
    let a: Vec<_> = (0..16)
        .map(|i| C64::new(0.02 * (i + 1) as f64, 0.003 * ((i % 3) as f64 - 1.0)))
        .collect();
    let b: Vec<_> = (0..16)
        .map(|i| C64::new(0.015 * ((i * 7) % 17) as f64, -0.002 * (i % 4) as f64))
        .collect();
    let (s, c) = redheffer_star(&a, &b, n, 1e12, || false).unwrap();
    assert!(c.is_finite());
    for col in 0..4 {
        let mut matrix = vec![C64::ZERO; 16];
        let mut rhs = vec![C64::ZERO; 4];
        for i in 0..2 {
            matrix[i * 4 + i] = C64::ONE;
            matrix[(2 + i) * 4 + 2 + i] = C64::ONE;
            for j in 0..2 {
                matrix[i * 4 + 2 + j] = -a[(2 + i) * 4 + 2 + j];
                matrix[(2 + i) * 4 + j] = -b[i * 4 + j];
            }
            if col < 2 {
                rhs[i] = a[(2 + i) * 4 + col];
            } else {
                rhs[2 + i] = b[i * 4 + col];
            }
        }
        let lu = lu_complex(&matrix, 4).unwrap();
        lu.solve(&mut rhs);
        for row in 0..4 {
            let mut want = if row < 2 && col < 2 {
                a[row * 4 + col]
            } else if row >= 2 && col >= 2 {
                b[row * 4 + col]
            } else {
                C64::ZERO
            };
            for j in 0..2 {
                want = want
                    + if row < 2 {
                        a[row * 4 + 2 + j] * rhs[2 + j]
                    } else {
                        b[row * 4 + j] * rhs[j]
                    };
            }
            assert!(relative(s[row * 4 + col], want) < 1e-12);
        }
    }
    let mut r = vec![C64::new(0.4, 0.2); 4];
    let before = r.clone();
    let e = [
        crate::exp_i(C64::new(0.5, 0.1)),
        crate::exp_i(C64::new(0.0, 900.0)),
    ];
    propagate_reflection(&mut r, &e).unwrap();
    for i in 0..2 {
        for j in 0..2 {
            assert_eq!(r[i * 2 + j], e[i] * before[i * 2 + j] * e[j]);
        }
    }
    assert!(propagate_reflection(&mut r, &[C64::from_re(1.1), C64::ONE]).is_err());
    for n in [MAX_COMPOSITION_MODES, MAX_COMPOSITION_MODES + 1] {
        let mut through = vec![C64::ZERO; 4 * n * n];
        for i in 0..n {
            through[i * 2 * n + n + i] = C64::ONE;
            through[(n + i) * 2 * n + i] = C64::ONE;
        }
        let result = redheffer_star(&through, &through, n, 1e12, || false);
        if n == MAX_COMPOSITION_MODES {
            assert_eq!(result.unwrap().0, through);
        } else {
            assert!(result.is_err());
        }
    }
}
#[test]
fn opaque_network_edge_keeps_finite_admittance() {
    let gas = gas();
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.001,
            length: 1000.0,
        }],
    };
    let w = TAU * 20000.0;
    let (y, c, _) = crate::network::network_chain_admittance(
        &duct,
        &gas,
        w,
        LossModel::Bessel,
        crate::network::JunctionCorrections::ChainLocal,
    )
    .unwrap();
    let wave = crate::segment_wave(&gas, 0.001, w, LossModel::Bessel).unwrap();
    assert!(
        relative(y[0], wave.characteristic_impedance.recip()) < 1e-12
            && relative(y[3], wave.characteristic_impedance.recip()) < 1e-12
    );
    assert!(y[1].abs() < 1e-300 && y[2].abs() < 1e-300 && c.is_finite());
}

/// G0/G2, DEF-01: a lossless rigid-capped cylinder at kL = pi has an input
/// impedance pole. Z refuses there (|1-r| ~ 1e-16 exceeds the 1e12 bound),
/// while the input reflection is exact: R = exp(2ikL) = +1. Off the pole,
/// the Z-free reference change agrees with (Z-Zc)/(Z+Zc) from the Z route.
#[test]
fn input_reflection_is_exact_where_lossless_impedance_poles_refuse() {
    let start = Instant::now();
    let g = gas();
    let (radius, length) = (0.005, 0.30);
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let zc = g.characteristic_impedance / (PI * radius * radius);
    let pole = PI * g.sound_speed / length;
    let at = |omega: f64| {
        let z = crate::input_impedance_wall(
            &duct,
            &g,
            omega,
            LossModel::Lossless,
            Termination::Closed,
            None,
        );
        let r = input_reflection_wall(
            &duct,
            &g,
            omega,
            LossModel::Lossless,
            Termination::Closed,
            None,
        );
        (z, r)
    };
    let (z, r) = at(pole);
    assert!(matches!(z, Err(DuctError::Singular)), "{z:?}");
    let r = r.unwrap();
    let unit = r.referenced_to(zc).unwrap();
    let pole_error = (unit - C64::ONE).abs();
    assert!((r.reflection.abs() - 1.0).abs() < 1e-12 && pole_error < 1e-9);
    emit(
        "def01-input-reflection-at-lossless-pole",
        "cylinder-0.30m-closed",
        0,
        None,
        pole_error,
        r.condition_estimate,
        1e-9,
        start,
        "",
    );
    for omega in [pole / 3.0, 0.7 * pole, 1.37 * pole] {
        let (z, r) = at(omega);
        let z = z.unwrap().impedance;
        let from_z = (z - C64::from_re(zc)) / (z + C64::from_re(zc));
        let error = relative(r.unwrap().referenced_to(zc).unwrap(), from_z);
        assert!(error < 1e-12, "omega {omega}: {error}");
    }
    assert!(r.referenced_to(0.0).is_err() && r.referenced_to(f64::NAN).is_err());
    // An active load at the reflection-coordinate pole Z_L = -R refuses here;
    // the Z route keeps its projective contraction for that algebraic probe.
    let through = [TwoPort::through(1.0).unwrap()];
    let active = Some(C64::from_re(-1.0));
    assert!(matches!(
        compose_chain_reflection(&through, active, ComposeOptions::default()),
        Err(DuctError::Singular)
    ));
    assert!(compose_chain(&through, active, ComposeOptions::default()).is_ok());
}

/// G0, DEF-09: fs-phs derives each Helmholtz neck end correction from the
/// constants of this crate's terminations. Unflanged is the identical
/// Levine-Schwinger fit (0.6133 a); flanged matches the multimodal model's
/// low-ka value within its converged spread (0.8216 vs ~0.822 a).
#[test]
fn neck_end_corrections_match_the_duct_terminations() {
    let g = gas();
    let a = 0.01;
    let omega = 1.0e-3 * g.sound_speed / a;
    let area = PI * a * a;
    let delta = |termination| {
        let (z, _) = crate::termination_impedance(termination, &g, a, omega).unwrap();
        -z.unwrap().im * area / (omega * g.density)
    };
    let phs = |flange| fs_phs::neck_end_correction(a, flange).unwrap();
    let unflanged = delta(Termination::UnflangedOpen) / phs(fs_phs::MouthFlange::Unflanged);
    assert!(
        (unflanged - 1.0).abs() < 1e-12,
        "unflanged ratio {unflanged}"
    );
    // The duct's flanged load is the 16-mode flanged pipe, qualified to
    // 1e-3 a against 0.8216 a by `flanged_mode_ladder_converges_to_0p8216`
    // (production 0.8224 a); the neck must sit inside that same bound.
    let flanged = (delta(Termination::FlangedOpen) - phs(fs_phs::MouthFlange::Flanged)) / a;
    assert!(flanged.abs() < 1e-3, "flanged end correction gap {flanged} a");
}

#[test]
fn open_tone_hole_flows_conserve_dc_and_match_prefix_suffix_split() {
    let g = gas();
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.0075,
                length: 0.14,
            },
            Segment::ToneHole {
                hole_radius: 0.0032,
                chimney_height: 0.006,
                bore_radius: 0.0075,
                state: HoleState::Closed,
            },
            Segment::Cylinder {
                radius: 0.0075,
                length: 0.06,
            },
            Segment::ToneHole {
                hole_radius: 0.0036,
                chimney_height: 0.0065,
                bore_radius: 0.0075,
                state: HoleState::Open,
            },
            Segment::Cylinder {
                radius: 0.0075,
                length: 0.05,
            },
            Segment::ToneHole {
                hole_radius: 0.0040,
                chimney_height: 0.007,
                bore_radius: 0.0075,
                state: HoleState::Vent(0.64),
            },
            Segment::Cone {
                inlet_radius: 0.0075,
                outlet_radius: 0.025,
                length: 0.09,
            },
        ],
    };
    let wall = WallModel::from_pin(None);
    // Low-frequency volume-flow conservation: sum(u_hole/u_in) + u_bell/u_in -> 1.
    let (dc_resp, dc_holes) = crate::input_impedance_with_wall_and_holes(
        &duct,
        &g,
        TAU * 0.25,
        LossModel::Bessel,
        Termination::FlangedOpen,
        wall,
    )
    .unwrap();
    assert_eq!(dc_holes.len(), 2);
    assert!((dc_holes[0].axial_position_m - 0.20).abs() < 1e-15);
    assert!((dc_holes[1].axial_position_m - 0.25).abs() < 1e-15);
    assert!((dc_holes[1].effective_radius_m - 0.0040 * 0.8).abs() < 1e-15);
    let total_dc =
        dc_resp.u_mouth_over_u_in + dc_holes[0].u_hole_over_u_in + dc_holes[1].u_hole_over_u_in;
    assert!(
        (total_dc - C64::ONE).abs() < 1e-5,
        "low-frequency volume flow sum {total_dc:?}"
    );
    // Midpoint-pressure formula agrees with explicit chain-splitting flow difference.
    let options = ComposeOptions::default();
    for hz in [110.0, 440.0, 1175.0, 2637.0] {
        let omega = TAU * hz;
        let expected_resp = crate::input_impedance_with_wall(
            &duct,
            &g,
            omega,
            LossModel::Bessel,
            Termination::FlangedOpen,
            wall,
        )
        .unwrap();
        let (resp, holes) = crate::input_impedance_with_wall_and_holes(
            &duct,
            &g,
            omega,
            LossModel::Bessel,
            Termination::FlangedOpen,
            wall,
        )
        .unwrap();
        assert_eq!(
            (resp.impedance.re.to_bits(), resp.impedance.im.to_bits()),
            (
                expected_resp.impedance.re.to_bits(),
                expected_resp.impedance.im.to_bits()
            )
        );
        assert_eq!(
            (
                resp.u_mouth_over_u_in.re.to_bits(),
                resp.u_mouth_over_u_in.im.to_bits()
            ),
            (
                expected_resp.u_mouth_over_u_in.re.to_bits(),
                expected_resp.u_mouth_over_u_in.im.to_bits()
            )
        );
        let (chain, prepared_holes) = super::duct_chain_with_wall_and_holes(
            &duct,
            &g,
            omega,
            LossModel::Bessel,
            wall,
            options,
            true,
            || false,
        )
        .unwrap();
        let (z_load, _) =
            crate::termination_impedance(Termination::FlangedOpen, &g, 0.025, omega).unwrap();
        for (flow, prep) in holes.iter().zip(&prepared_holes) {
            let b = prep.block_index;
            let suffix = compose_chain(&chain.blocks[b..], z_load, options).unwrap();
            let prefix = compose_chain(&chain.blocks[..b], Some(suffix.z_in), options).unwrap();
            let after = compose_chain(&chain.blocks[b + 1..], z_load, options).unwrap();
            let through_hole = compose_chain(&chain.blocks[..=b], Some(after.z_in), options).unwrap();
            let u_shunt = prefix.u_out_over_u_in - through_hole.u_out_over_u_in;
            let seg = &duct.segments[if b == 3 { 3 } else { 5 }];
            let Segment::ToneHole {
                hole_radius,
                chimney_height,
                bore_radius,
                state,
            } = *seg
            else {
                unreachable!()
            };
            let z_h = crate::tone_hole_shunt_with_wall(
                &g,
                hole_radius,
                chimney_height,
                state,
                omega,
                LossModel::Bessel,
                bore_radius,
                wall,
            )
            .unwrap();
            let want = u_shunt
                * (z_h.scale(prep.sigma) / prep.open_chimney_impedance)
                * prep.chimney_u_mouth_over_u_in;
            let err = relative(flow.u_hole_over_u_in, want);
            assert!(err < 1e-12, "hz={hz} block={b}: err={err}");
        }
    }
}
