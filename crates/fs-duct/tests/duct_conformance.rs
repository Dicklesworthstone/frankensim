//! fs-duct CONFORMANCE battery (music bead
//! `frankensim-music-v8-root-3ez8g.13.2`): the reimplementation
//! contract, exercised entirely through the public surface with
//! INDEPENDENT oracles. The inline `#[cfg(test)]` unit modules stay
//! untouched (the conformance surface GROWS; nothing is weakened) —
//! this file is what a consumer or reimplementer reads to learn the
//! claim surface, next to `tests/ernoult_corpus.rs` (the registered
//! corpus fixture; not duplicated here).
//!
//! Cases:
//! - dt-001: the closed lossless cylinder matches the analytic
//!   cot form `Z_in = +i Z0 cot(kL)` (e^{-iωt}) to 1e-10.
//! - dt-002: a vanishing-taper cone degenerates to the cylinder, and
//!   the closed cylinder's first sweep peak sits at the quarter-wave
//!   `c/(4L)`.
//! - dt-003: viscothermal Q follows the sqrt-frequency law across
//!   the first two resonances (an independent scaling oracle for the
//!   wide-tube loss model).
//! - dt-004: refusals by name — the narrow-tube shear-number floor
//!   and non-physical parameters.
//! - dt-005: radiation end corrections are QUANTITATIVE — the
//!   unflanged (0.6133 a) and flanged (0.8216 a) loads move the
//!   quarter-wave peak by the textbook effective-length ratios.
//! - dt-008: the diffusive limit is Poiseuille. An open narrow tube
//!   reads `Z_in(0⁺) → R'L = 8μL/(πa⁴)` (real), with a rigid or a
//!   yielding wall, and a wall-loaded line stays passive.

use fs_duct::{
    Duct, DuctError, LossModel, Segment, Termination, impedance_peaks, impedance_sweep,
    input_impedance,
};
use fs_material::gas::{GasSpec, GasState};

const TAU: f64 = core::f64::consts::TAU;

/// G2: the existing Duct chain embeds smooth DATA between cylinders;
/// an independent closed-form ABCD product carries the same table load.
#[test]
fn dt_007_smooth_profile_chain() {
    use fs_duct::{RadiusProfile, RiccatiSteps, TabulatedLoad, input_impedance_tabulated};
    use fs_math::{c64::C64, det};
    use std::sync::Arc;
    let start=std::time::Instant::now();
    let gas=air();
    let (r0,flare,length)=(0.005,2.0,0.5);
    let rout=r0*det::exp(flare*length/2.0);
    let profile=Arc::new(RadiusProfile::exponential(r0,flare,length).unwrap());
    let load=C64::new(1.5e6,0.0);
    let table=TabulatedLoad::try_new(vec![TAU*10.0,TAU*1500.0],vec![load,load],"analytic resistive mouth",0.0).unwrap();
    let mul=|a:[C64;4],b:[C64;4]| [
        a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3],
        a[2]*b[0]+a[3]*b[2],a[2]*b[1]+a[3]*b[3],
    ];
    let cylinder=|radius:f64,l:f64,w:f64| {
        let theta=w*l/gas.sound_speed;
        let z=gas.density*gas.sound_speed/(core::f64::consts::PI*radius*radius);
        [C64::from_re(det::cos(theta)),C64::new(0.0,-z*det::sin(theta)),
         C64::new(0.0,-det::sin(theta)/z),C64::from_re(det::cos(theta))]
    };
    let exact=|w:f64| {
        // p=exp(-m*x/2)y, y''+beta^2*y=0, and U=S*p'/(i*w*rho).
        let b2=(w/gas.sound_speed)*(w/gas.sound_speed)-flare*flare/4.0;
        let (c,s)=if b2>0.0 {
            let b=det::sqrt(b2); (det::cos(b*length),det::sin(b*length)/b)
        } else if b2<0.0 {
            let b=det::sqrt(-b2); let e=det::exp(b*length);
            ((e+1.0/e)/2.0,(e-1.0/e)/(2.0*b))
        } else { (1.0,length) };
        let area=core::f64::consts::PI*r0*r0;
        let e=det::exp(flare*length/2.0);
        let horn=[C64::from_re(e*(c-flare*s/2.0)),C64::new(0.0,-w*gas.density/area/e*s),
            C64::new(0.0,-e*area*w/(gas.density*gas.sound_speed*gas.sound_speed)*s),
            C64::from_re((c+flare*s/2.0)/e)];
        let m=mul(mul(cylinder(r0,0.05,w),horn),cylinder(rout,0.04,w));
        (m[0]*load+m[1])/(m[2]*load+m[3])
    };
    // Independent log-parabolic peak positions on the declared uniform
    // frequency grid, applied equally to the two computed impedance rows.
    let peaks=|rows:&[C64]|->Vec<f64> {
        rows.windows(3).enumerate().filter_map(|(i,row)| {
            let (a,b,c)=(det::ln(row[0].abs()),det::ln(row[1].abs()),det::ln(row[2].abs()));
            if b>a && b>c {
                let offset=0.5*(a-c)/(a-2.0*b+c);
                Some(40.0+960.0*(i as f64+1.0+offset)/99.0)
            } else {None}
        }).collect()
    };
    let expected:Vec<_>=(0..100).map(|i|exact(TAU*(40.0+960.0*f64::from(i)/99.0))).collect();
    let expected_peaks=peaks(&expected);
    assert!(expected_peaks.len()>=2);
    for steps in [1024,2048,4096,8192] {
        let duct=Duct {segments:vec![Segment::Cylinder {radius:r0,length:0.05},
            Segment::Profile {profile:Arc::clone(&profile),steps:RiccatiSteps(steps)},
            Segment::Cylinder {radius:rout,length:0.04}]};
        let mut worst=0.0_f64;
        let mut observed=Vec::new();
        for (i,&want) in expected.iter().enumerate() {
            let omega=TAU*(40.0+960.0*i as f64/99.0);
            let got=input_impedance_tabulated(&duct,&gas,omega,LossModel::Lossless,&table,None).unwrap();
            let error=(got.impedance-want).abs()/want.abs();
            worst=worst.max(error); observed.push(got.impedance);
        }
        let actual_peaks=peaks(&observed);
        assert_eq!(actual_peaks.len(),expected_peaks.len());
        let peak_error=actual_peaks.iter().zip(&expected_peaks).map(|(a,b)|(a/b-1.0).abs()).fold(0.0,f64::max);
        assert!(worst<1e-9 && peak_error<1e-9,"steps={steps}: Zerr={worst:e}, peakerr={peak_error:e}");
        eprintln!("{{\"suite\":\"fs-duct\",\"case\":\"dt-007\",\"profile\":\"exp\",\"inputs\":\"analytic cylinder .05m, exp S0*exp(2x) .5m, cylinder .04m; resistive TabulatedLoad1.5e6 Pa s/m3; air293.15K\",\"axis\":\"Magnus steps\",\"rung\":{steps},\"steps\":{steps},\"observed\":{worst:e},\"z_rel_err\":{worst:e},\"peak_rel_err\":{peak_error:e},\"oracle\":\"independent closed-form ABCD product\",\"tol\":1e-9,\"verdict\":\"pass\",\"elapsed_us\":{}}}",start.elapsed().as_micros());
    }
}

/// G2: the additive exact-WH termination transports through a lossless
/// cylinder according to the independent closed-form acoustic two-port.
#[test]
fn dt_006_levine_schwinger_thin_wall_termination() {
    use fs_duct::levine_schwinger::{
        LS_MAX_KA, LsQuadrature, levine_schwinger, levine_schwinger_impedance,
    };
    use fs_math::{c64::C64, det};
    let state = air();
    let radius = 0.008;
    let length = 0.30;
    let duct = cylinder(radius, length);
    let q = LsQuadrature { nodes: 256 };
    let zc = state.density * state.sound_speed / (core::f64::consts::PI * radius * radius);
    for i in 0..50 {
        let start = std::time::Instant::now();
        let ka = 0.01 + (3.8 - 0.01) * f64::from(i) / 49.0;
        let omega = ka * state.sound_speed / radius;
        let load = levine_schwinger_impedance(&state, radius, omega, q).expect("exact mouth load");
        let phase = omega / state.sound_speed * length;
        let (c, s) = (det::cos(phase), det::sin(phase));
        // Algebraically the requested tan(kL) transfer, scaled by cos(kL)
        // to avoid a removable tangent pole at a quarter-wave frequency.
        let want = (load.scale(c) - C64::new(0.0, zc * s))
            / (C64::from_re(zc * c) - C64::new(0.0, s) * load).scale(1.0 / zc);
        let got = input_impedance(
            &duct,
            &state,
            omega,
            LossModel::Lossless,
            Termination::UnflangedThinWall { quadrature: q },
        )
        .expect("exact LS TMM termination");
        let error = (got.impedance - want).abs() / want.abs();
        assert!(
            error < 1e-10,
            "ka={ka}: {:?} vs {want:?}, relative {error}",
            got.impedance
        );
        let ls = levine_schwinger(ka, q).unwrap();
        let reference = levine_schwinger(ka, LsQuadrature { nodes: 512 }).unwrap();
        let quadrature_error = (ls.r_abs / reference.r_abs - 1.0)
            .abs()
            .max((ls.l_over_a / reference.l_over_a - 1.0).abs());
        assert!(quadrature_error < 1e-3);
        eprintln!(
            "{{\"suite\":\"fs-duct\",\"case\":\"dt-006\",\"inputs\":{{\"fixture\":\"0.30m x 8mm lossless cylinder\",\"provenance\":\"LS1948 p383/V16; independent analytic two-port transport\",\"ka\":{ka}}},\"ka\":{ka},\"nodes\":256,\"r_abs\":{},\"l_over_a\":{},\"ref_r_abs\":{},\"ref_l_over_a\":{},\"observed\":{{\"input_re\":{},\"input_im\":{}}},\"oracle\":{{\"input_re\":{},\"input_im\":{}}},\"rel_err\":{quadrature_error},\"err\":{error},\"tol\":1e-10,\"reference_tol\":1e-3,\"verdict\":\"pass\",\"elapsed_us\":{}}}",
            ls.r_abs,
            ls.l_over_a,
            reference.r_abs,
            reference.l_over_a,
            got.impedance.re,
            got.impedance.im,
            want.re,
            want.im,
            start.elapsed().as_micros()
        );
    }
    let omega = (*LS_MAX_KA * (1.0 + 1e-9)) * state.sound_speed / radius;
    assert!(matches!(
        input_impedance(
            &duct,
            &state,
            omega,
            LossModel::Lossless,
            Termination::UnflangedThinWall { quadrature: q }
        ),
        Err(DuctError::LevineSchwingerCutoff { .. })
    ));
}

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

fn cylinder(radius: f64, length: f64) -> Duct {
    Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    }
}

/// First sweep peak [Hz] of a duct under the given loss/termination.
fn first_peak_hz(duct: &Duct, state: &GasState, loss: LossModel, term: Termination) -> f64 {
    let sweep =
        impedance_sweep(duct, state, TAU * 60.0, TAU * 2500.0, 24_000, loss, term).expect("sweep");
    let peaks = impedance_peaks(&sweep);
    assert!(!peaks.is_empty(), "no impedance peak in the sweep window");
    sweep[peaks[0]].omega / TAU
}

#[test]
fn dt_001_closed_cylinder_matches_the_cot_form() {
    let state = air();
    let (radius, length) = (0.008, 0.45);
    let duct = cylinder(radius, length);
    let z0 = state.density * state.sound_speed / (core::f64::consts::PI * radius * radius);
    let mut worst = 0.0f64;
    for f_hz in [123.0f64, 411.0, 977.0] {
        let omega = TAU * f_hz;
        let k = omega / state.sound_speed;
        let resp = input_impedance(
            &duct,
            &state,
            omega,
            LossModel::Lossless,
            Termination::Closed,
        )
        .expect("closed cylinder");
        // e^{-iωt}: Z_in = +i Z0 cot(kL) (the low-frequency closed
        // pipe is a compliance, Z = 1/(-iωC) = +i/(ωC) — positive
        // imaginary under this convention; the -i form belongs to
        // e^{+iωt}).
        let cot = (k * length).cos() / (k * length).sin();
        let (want_re, want_im) = (0.0, z0 * cot);
        let dre = (resp.impedance.re - want_re).abs() / z0;
        let dim = (resp.impedance.im - want_im).abs() / z0.max(want_im.abs());
        worst = worst.max(dre.max(dim));
        assert!(
            dre < 1.0e-10 && dim < 1.0e-10,
            "{f_hz} Hz: Z = ({}, {}) vs +i Z0 cot(kL) = (0, {want_im:.6e})",
            resp.impedance.re,
            resp.impedance.im
        );
    }
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"dt-001-cot-form\",\"worst_rel\":{worst:.3e},\
         \"convention\":\"e^{{-i omega t}}\",\"verdict\":\"pass\"}}"
    );
}

#[test]
fn dt_002_cone_degenerates_and_quarter_wave_pins() {
    let state = air();
    let (radius, length) = (0.008, 0.45);
    // Vanishing taper: the truncated cone must land on the cylinder.
    let cone = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: radius,
            outlet_radius: radius * 1.000_001,
            length,
        }],
    };
    let cyl = cylinder(radius, length);
    let mut worst = 0.0f64;
    for f_hz in [200.0f64, 700.0, 1500.0] {
        let omega = TAU * f_hz;
        let zc = input_impedance(
            &cone,
            &state,
            omega,
            LossModel::Lossless,
            Termination::IdealOpen,
        )
        .expect("cone")
        .impedance;
        let zy = input_impedance(
            &cyl,
            &state,
            omega,
            LossModel::Lossless,
            Termination::IdealOpen,
        )
        .expect("cylinder")
        .impedance;
        let num = ((zc.re - zy.re).powi(2) + (zc.im - zy.im).powi(2)).sqrt();
        let den = (zy.re.powi(2) + zy.im.powi(2)).sqrt().max(1.0e-300);
        let rel = num / den;
        worst = worst.max(rel);
        assert!(
            rel < 1.0e-4,
            "cone-cylinder degeneracy at {f_hz} Hz: {rel:.3e}"
        );
    }
    // Quarter-wave pin: ideal-open pipe peaks at c/(4L).
    let f1 = first_peak_hz(&cyl, &state, LossModel::Lossless, Termination::IdealOpen);
    let quarter = state.sound_speed / (4.0 * length);
    let dev = (f1 / quarter - 1.0).abs();
    assert!(dev < 1.0e-3, "first peak {f1:.2} vs c/4L {quarter:.2}");
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"dt-002-cone-quarterwave\",\
         \"cone_degeneracy_worst\":{worst:.3e},\"quarter_wave_dev\":{dev:.3e},\
         \"verdict\":\"pass\"}}"
    );
}

#[test]
fn dt_003_viscothermal_q_follows_the_sqrt_frequency_law() {
    let state = air();
    let (radius, length) = (0.006, 0.60);
    let duct = cylinder(radius, length);
    let sweep = impedance_sweep(
        &duct,
        &state,
        TAU * 80.0,
        TAU * 900.0,
        60_000,
        LossModel::WideTube,
        Termination::IdealOpen,
    )
    .expect("sweep");
    let peaks = impedance_peaks(&sweep);
    assert!(
        peaks.len() >= 2,
        "need two resonances, found {}",
        peaks.len()
    );
    let q_at = |peak_idx: usize| -> (f64, f64) {
        let mag = |i: usize| (sweep[i].impedance.re.powi(2) + sweep[i].impedance.im.powi(2)).sqrt();
        let m0 = mag(peak_idx);
        let target = m0 / 2.0f64.sqrt();
        let mut lo = peak_idx;
        while lo > 0 && mag(lo) > target {
            lo -= 1;
        }
        let mut hi = peak_idx;
        while hi + 1 < sweep.len() && mag(hi) > target {
            hi += 1;
        }
        let f0 = sweep[peak_idx].omega / TAU;
        let bw = (sweep[hi].omega - sweep[lo].omega) / TAU;
        (f0, f0 / bw.max(1.0e-9))
    };
    let (f1, q1) = q_at(peaks[0]);
    let (f2, q2) = q_at(peaks[1]);
    // Wide-tube wall losses: alpha ~ sqrt(f) so Q = k/(2 alpha) ~ sqrt(f).
    let want = (f2 / f1).sqrt();
    let got = q2 / q1;
    let dev = (got / want - 1.0).abs();
    assert!(
        dev < 0.05,
        "Q scaling: Q2/Q1 = {got:.3} vs sqrt(f2/f1) = {want:.3} (f {f1:.0}/{f2:.0}, Q {q1:.0}/{q2:.0})"
    );
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"dt-003-q-scaling\",\"f1\":{f1:.1},\"q1\":{q1:.1},\
         \"f2\":{f2:.1},\"q2\":{q2:.1},\"ratio\":{got:.4},\"sqrt_law\":{want:.4},\
         \"dev\":{dev:.4},\"verdict\":\"pass\"}}"
    );
}

#[test]
fn dt_004_refusals_fire_by_name() {
    let state = air();
    // Narrow-tube shear-number floor: a capillary at low frequency.
    let narrow = cylinder(1.0e-4, 0.1);
    assert!(matches!(
        input_impedance(
            &narrow,
            &state,
            TAU * 20.0,
            LossModel::WideTube,
            Termination::Closed
        ),
        Err(DuctError::TooNarrow { .. })
    ));
    // Non-physical parameters.
    let bad = cylinder(-0.01, 0.1);
    assert!(matches!(
        input_impedance(
            &bad,
            &state,
            TAU * 100.0,
            LossModel::Lossless,
            Termination::Closed
        ),
        Err(DuctError::BadParameter { .. })
    ));
    let empty = Duct { segments: vec![] };
    assert!(
        input_impedance(
            &empty,
            &state,
            TAU * 100.0,
            LossModel::Lossless,
            Termination::Closed
        )
        .is_err()
    );
    println!("{{\"suite\":\"fs-duct\",\"case\":\"dt-004-refusals\",\"verdict\":\"pass\"}}");
}

#[test]
fn dt_005_end_corrections_are_quantitative() {
    let state = air();
    let (radius, length) = (0.010, 0.40);
    let duct = cylinder(radius, length);
    let f_ideal = first_peak_hz(&duct, &state, LossModel::Lossless, Termination::IdealOpen);
    let f_unflanged = first_peak_hz(
        &duct,
        &state,
        LossModel::Lossless,
        Termination::UnflangedOpen,
    );
    let f_flanged = first_peak_hz(&duct, &state, LossModel::Lossless, Termination::FlangedOpen);
    // Effective lengths: L + 0.6133 a and L + 0.8216 a.
    let want_unflanged = f_ideal * length / (0.6133f64.mul_add(radius, length));
    let want_flanged = f_ideal * length / (0.8216f64.mul_add(radius, length));
    let dev_u = (f_unflanged / want_unflanged - 1.0).abs();
    let dev_f = (f_flanged / want_flanged - 1.0).abs();
    assert!(
        dev_u < 0.01,
        "unflanged peak {f_unflanged:.2} vs effective-length prediction {want_unflanged:.2}"
    );
    assert!(
        dev_f < 0.01,
        "flanged peak {f_flanged:.2} vs effective-length prediction {want_flanged:.2}"
    );
    assert!(f_flanged < f_unflanged && f_unflanged < f_ideal);
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"dt-005-end-corrections\",\"f_ideal\":{f_ideal:.2},\
         \"f_unflanged\":{f_unflanged:.2},\"f_flanged\":{f_flanged:.2},\
         \"dev_unflanged\":{dev_u:.4},\"dev_flanged\":{dev_f:.4},\"verdict\":\"pass\"}}"
    );
}

/// G2: below the diffusion corner `ω ≪ 1/(R'C'L²)` an open (`p = 0`)
/// narrow tube is a pure Poiseuille resistor: `Z_in = R'L (1 + O(ω))`
/// with `R' = 8μ/(πa⁴)`, independent of the wall (the shunt enters at
/// `O(Z'Y'L²)`). The basis pair is `e^{-iωt}`, so the telegraph `Zc`
/// must be conjugated into it; reusing the `+iω` value left `k Zc`
/// real and read `Z_in(0⁺) = −i R'L`, and the wall fold flipped the
/// series resistance (`Re Z_in < 0`, `|R| > 1`).
#[test]
fn dt_008_diffusive_limit_is_poiseuille_and_walls_stay_passive() {
    use fs_duct::input_impedance_wall;
    use fs_phs::WallPin;
    let start = std::time::Instant::now();
    let state = air();
    let (radius, length) = (0.001, 0.5);
    let duct = cylinder(radius, length);
    let r_poiseuille =
        8.0 * state.dynamic_viscosity * length / (core::f64::consts::PI * radius * radius * radius * radius);
    let soft = WallPin {
        surface_density: 15.0,
        stiffness_per_area: 3.0e5,
        resistance: 8000.0,
    };
    let omega = 1.0e-3;
    for (label, wall, tol) in [("rigid", None, 1.0e-4), ("soft", Some(&soft), 1.0e-3)] {
        let z = input_impedance_wall(
            &duct,
            &state,
            omega,
            LossModel::Bessel,
            Termination::IdealOpen,
            wall,
        )
        .expect("diffusive Z_in")
        .impedance;
        let re_dev = (z.re / r_poiseuille - 1.0).abs();
        let im_rel = z.im.abs() / r_poiseuille;
        assert!(
            re_dev < tol && im_rel < tol,
            "{label} wall: Z_in(ω=1e-3) = {:.6e}{:+.6e}i, want R'L = {r_poiseuille:.6e} \
             (Re dev {re_dev:.3e}, |Im|/R'L {im_rel:.3e}, tol {tol:e})",
            z.re,
            z.im
        );
    }
    // Passivity over the band where the soft wall's mass-spring corner
    // (√(K/σ)/2π ≈ 22.5 Hz) and the diffusion corner both live.
    let mut min_re = f64::INFINITY;
    for i in 0..=240 {
        let f = 0.05 * fs_math::det::exp(core::f64::consts::LN_10 * f64::from(i) / 40.0);
        for wall in [None, Some(&soft)] {
            let z = input_impedance_wall(
                &duct,
                &state,
                TAU * f,
                LossModel::Bessel,
                Termination::IdealOpen,
                wall,
            )
            .expect("sweep Z_in")
            .impedance;
            assert!(
                z.re >= 0.0,
                "Re Z_in = {:.3e} < 0 at {f:.3} Hz (wall {})",
                z.re,
                wall.is_some()
            );
            min_re = min_re.min(z.re / r_poiseuille);
        }
    }
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"dt-008-diffusive-poiseuille\",\
         \"r_poiseuille\":{r_poiseuille:.6e},\"min_re_over_rp\":{min_re:.3e},\
         \"elapsed_ms\":{},\"verdict\":\"pass\"}}",
        start.elapsed().as_millis()
    );
}
