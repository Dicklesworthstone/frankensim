//! Manufactured resistive valve paths; no measured instrument or microphone.
use super::*;
use fs_duct::{TabulatedLoad, TabulatedMouth};
use fs_material::gas::GasSpec;
use fs_math::c64::C64;

struct Fixture {
    gas: GasState,
    ducts: Vec<Duct>,
    tables: Vec<TabulatedLoad>,
    line: MmLineConfig,
}
fn lip() -> LipIsland {
    LipIsland { mass_kg: 1.8e-4, stiffness_n_m: 350.0, damping_n_s_m: 0.05,
        width_m: 0.012, rest_gap_m: 5e-4, face_area_m2: 2e-5,
        provenance: "manufactured modal mouth-flow test".into() }
}
fn config() -> BrassSequenceConfig {
    BrassSequenceConfig { samples: 256, max_block_len: 256,
        initial_blowing_pressure_pa: 1500.0, initial_lip_tension: 1.0, initial_combo: 0 }
}
fn events() -> Vec<BrassEvent> {
    vec![BrassEvent { sample: 0, control: BrassControl::SetBlowingPressure(400.0) },
        BrassEvent { sample: 31, control: BrassControl::SetValve { combo: 1, fade_samples: 11 } },
        BrassEvent { sample: 31, control: BrassControl::SetBlowingPressure(800.0) },
        // Same-path assignment must not restart the active fade.
        BrassEvent { sample: 33, control: BrassControl::SetValve { combo: 1, fade_samples: 99 } },
        BrassEvent { sample: 127, control: BrassControl::SetValve { combo: 0, fade_samples: 0 } },
        BrassEvent { sample: 192, control: BrassControl::SetValve { combo: 1, fade_samples: 16 } },
        BrassEvent { sample: 255, control: BrassControl::SetBlowingPressure(0.0) }]
}
impl Fixture {
    fn new() -> Self {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
        let line = MmLineConfig { sample_rate_hz: 48000, n_modes: 2,
            extra_slices: 1, loss: LossModel::Lossless };
        let dt = 1.0 / f64::from(line.sample_rate_hz);
        let radius = 0.002;
        let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
        let mut ducts = Vec::new();
        let mut tables = Vec::new();
        for (transit, resistance) in [(8.0, 0.7), (12.0, 1.3)] {
            let length = transit * gas.sound_speed * dt;
            ducts.push(Duct { segments: vec![Segment::Cylinder { radius, length }] });
            tables.push(TabulatedLoad::try_new(
                vec![core::f64::consts::TAU / (65_536.0 * dt), core::f64::consts::PI / dt],
                vec![C64::from_re(resistance * zc); 2], "resistive manufactured control", 0.0,
            ).unwrap().with_mouth(TabulatedMouth { axial_position_m: length,
                radius_m: radius, gas, dc_impedance: resistance * zc }).unwrap());
        }
        Self { gas, ducts, tables, line }
    }
    fn loads(&self) -> Vec<MmLoad<'_>> { self.tables.iter().map(MmLoad::Tabulated).collect() }
    fn build(&self, events: Vec<BrassEvent>) -> BrassMouthSequence {
        BrassMouthSequence::new(&self.ducts, &["short", "long"], &self.gas,
            MmBankLoad::PerCombo(&self.loads()), &self.line, lip(), 1.5e-6,
            config(), events, MouthFlowConfig::default()).unwrap()
    }
    fn source(&self) -> BrassSequence {
        BrassSequence::new(&self.ducts, &["short", "long"], &self.gas,
            MmBankLoad::PerCombo(&self.loads()), &self.line, lip(), 1.5e-6,
            config(), events()).unwrap()
    }
}

#[test]
fn mouth_flows_use_actual_outgoing_history_and_exact_original_valve_ramps() {
    let fixture = Fixture::new();
    let mut sequence = fixture.build(events());
    assert!(fixture.tables.iter().all(|table| table.z_at(1e-3).is_err()),
        "DC must use each separate static load, not table extrapolation");
    let taps: Vec<Vec<Vec<f64>>> = sequence.filters.iter().map(|path|
        path.iter().map(|filter| filter.quiescent_fir().unwrap().1.to_vec()).collect()).collect();
    let zc = sequence.source().voice().port_impedance();
    for (combo, resistance) in [0.7, 1.3].into_iter().enumerate() {
        let dc: f64 = taps[combo][0].iter().sum();
        let expected = 2.0 / ((1.0 + resistance) * zc);
        assert!((dc - expected).abs() < 2e-10 * expected);
        assert!(taps[combo][1].iter().all(|x| *x == 0.0));
        let receipt = &sequence.receipts()[combo];
        assert!(receipt.rungs.len() >= 2);
        let last = receipt.rungs.last().unwrap();
        assert!(last.tail_energy.iter().all(|x| *x < 1e-10));
        assert!(last.tap_change.iter().all(|x| *x < 1e-5));
    }
    let mut pressure = vec![0.0; 256];
    let mut outgoing = pressure.clone();
    let mut flow = vec![0.0; 512];
    sequence.block(&mut pressure, &mut outgoing, &mut flow).unwrap();
    let mut reference_pressure = vec![0.0; 256];
    let mut reference_outgoing = reference_pressure.clone();
    fixture.source().block_outgoing(&mut reference_pressure, &mut reference_outgoing).unwrap();
    assert_eq!(pressure, reference_pressure, "observer must not change mechanics");
    assert_eq!(outgoing, reference_outgoing);
    assert!(flow.iter().any(|q| q.abs() > 1e-8));
    let mut active = 0;
    let mut previous = 0;
    let mut begin = 0;
    let mut fade = 0;
    for sample in 0..256 {
        if let Some((combo, duration)) = match sample {
            31 => Some((1, 11)), 127 => Some((0, 0)), 192 => Some((1, 16)), _ => None,
        } {
            previous = active; active = combo; begin = sample; fade = duration;
        }
        for mode in 0..2 {
            let convolution = |combo: usize| -> f64 {
                taps[combo][mode].iter().take(sample + 1).enumerate()
                    .map(|(lag, &tap)| tap * outgoing[sample - lag]).sum()
            };
            let expected = if sample - begin < fade {
                let weight = (sample - begin + 1) as f64 / fade as f64;
                (1.0 - weight) * convolution(previous) + weight * convolution(active)
            } else { convolution(active) };
            assert!((flow[2*sample+mode] - expected).abs() < 2e-12 * expected.abs().max(1e-12),
                "sample {sample}, mode {mode}");
        }
    }
    assert_eq!(sequence.samples_rendered(), 256);
    assert_eq!(sequence.source().voice().lift_log().len(), 3);
}

#[test]
fn callback_partitions_preserve_pressure_waves_and_every_mouth_mode() {
    let fixture = Fixture::new();
    let render = |block| {
        let mut sequence = fixture.build(events());
        let mut p = vec![0.0; 256];
        let mut a = p.clone();
        let mut u = vec![0.0; 512];
        for ((p, a), u) in p.chunks_mut(block).zip(a.chunks_mut(block)).zip(u.chunks_mut(2*block)) {
            sequence.block(p, a, u).unwrap();
        }
        (p, a, u)
    };
    let baseline = render(256);
    assert_eq!(render(17), baseline);
    assert_eq!(render(1), baseline);
}

#[test]
fn budgets_loads_and_callback_dimensions_refuse_before_advancing() {
    let fixture = Fixture::new();
    let loads = fixture.loads();
    for config in [MouthFlowConfig { max_total_taps: 0, ..MouthFlowConfig::default() },
        MouthFlowConfig { max_taps: 300, ..MouthFlowConfig::default() },
        MouthFlowConfig { tail_tolerance: 0.0, ..MouthFlowConfig::default() }]
    {
        assert!(realize::admit(&fixture.ducts, &["a", "b"], &fixture.gas,
            MmBankLoad::PerCombo(&loads), &fixture.line, config).is_err());
    }
    assert!(realize::admit(&fixture.ducts, &["a", "b"], &fixture.gas,
        MmBankLoad::PerCombo(&loads[..1]), &fixture.line, MouthFlowConfig::default()).is_err());
    assert!(realize::admit(&fixture.ducts, &["a", "b"], &fixture.gas,
        MmBankLoad::Shared(loads[0]), &fixture.line, MouthFlowConfig::default()).is_err(),
        "longer path must not reuse the short path's mouth reference");
    assert!(realize::admit(&fixture.ducts, &["a", "b"], &fixture.gas,
        MmBankLoad::Shared(MmLoad::Matched), &fixture.line, MouthFlowConfig::default()).is_err());
    let mut sequence = fixture.build(events());
    let mut p = [7.0; 2]; let mut a = [9.0; 2]; let mut u = [11.0; 3];
    assert!(sequence.block(&mut p, &mut a, &mut u).is_err());
    assert_eq!(p, [7.0; 2]); assert_eq!(a, [9.0; 2]); assert_eq!(u, [11.0; 3]);
    assert_eq!(sequence.samples_rendered(), 0);
    assert_eq!(sequence.source().samples_rendered(), 0);
    assert!(sequence.block(&mut [], &mut [], &mut []).is_err());
    let mut u = [0.0; 4];
    sequence.block(&mut p, &mut a, &mut u).unwrap();
    let mut expected = fixture.build(events());
    let mut q = [0.0; 2]; let mut b = q; let mut v = [0.0; 4];
    expected.block(&mut q, &mut b, &mut v).unwrap();
    assert_eq!((p, a, u), (q, b, v));
}

#[test]
fn failed_physical_prefix_cannot_continue_with_stale_output_history() {
    let fixture = Fixture::new();
    let mut sequence = fixture.build(vec![BrassEvent { sample: 8,
        control: BrassControl::SetBlowingPressure(f64::MAX) }]);
    assert!(sequence.block(&mut [0.0; 16], &mut [0.0; 16], &mut [0.0; 32]).is_err());
    assert!(matches!(sequence.block(&mut [0.0], &mut [0.0], &mut [0.0; 2]),
        Err(RenderError::Poisoned)));
}

// ---- complete flanged modal bell load (the lip performer's `flanged_modal` mouth) -------------------------------
// Manufactured flare, no instrument claimed. Gauntlet: G0 closed form, G1 plane-wave limit, G2 FIR vs direct solve.
use fs_duct::Termination;
use fs_duct::modal::outlet;

fn air_at(temperature: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 101325.0).unwrap()
}
/// Cylinder of 12 mm radius, then a cone to `mouth_radius`.
fn flare(mouth_radius: f64) -> Duct {
    Duct { segments: vec![
        Segment::Cylinder { radius: 0.012, length: 0.25 },
        Segment::Cone { inlet_radius: 0.012, outlet_radius: mouth_radius, length: 0.15 },
    ] }
}

#[test]
fn flanged_matrix_plane_entry_is_the_closed_form_baffled_piston_at_any_mode_count() {
    // The same quantity the nonlinear bore's causal radiation weights pin from the other side (fs-duct
    // `piston_weights_sum_to_one_and_reproduce_the_struve_impedance`): Z00/(rho c/S) = 1 - J1(2ka)/ka -/+ i H1(2ka)/ka,
    // the imaginary sign being this crate's e^{-i omega t} convention (a mass-like load has Im Z < 0).
    // J1 and H1 come from their own power series here, independent of fs_phs (ka <= 0.4 keeps 2ka <= 0.8).
    let gas = air_at(293.15);
    let radius = 0.03;
    let z0 = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    for ka in [0.02, 0.05, 0.1, 0.2, 0.4] {
        let x: f64 = 2.0 * ka;
        let q = (x / 2.0) * (x / 2.0);
        let (mut j1, mut tj) = (0.0, x / 2.0);
        let (mut h1, mut th) = (0.0, 8.0 * q / (3.0 * core::f64::consts::PI));
        for m in 0..16 {
            j1 += tj;
            h1 += th;
            let m = f64::from(m);
            tj *= -q / ((m + 1.0) * (m + 2.0));
            th *= -q / ((m + 1.5) * (m + 2.5));
        }
        let expected = C64::new(1.0 - j1 / ka, -h1 / ka).scale(z0);
        let omega = ka * gas.sound_speed / radius;
        for modes in [1, 3, 8] {
            let z = fs_duct::flanged::flanged_modal_impedance(&gas, radius, omega, modes).unwrap();
            assert!((z[0] - expected).abs() < 1e-9 * z0, "ka {ka}, modes {modes}: {:?} vs {expected:?}", z[0]);
        }
    }
}

#[test]
fn flanged_modal_bore_reduces_to_todays_flanged_load_in_the_plane_wave_limit() {
    // Stated before the first run. A uniform tube many radii long: every higher mode is evanescent
    // (kappa L = 76 at ka = 0.4), so the complete 8-mode matrix and today's `FlangedOpen` (Zorumski, the 16-mode
    // Schur-reduced plane load) are the same physics up to mode truncation, which only moves the end correction
    // (fs-duct `flanged_mode_ladder_converges_to_0p8216`: 0.8488 a at one mode to 0.8216 a, so at eight modes
    // the residual is below 0.01 a). Phase error of the inlet reflectance is then 2 ka (dl/a) <= 2 ka 1e-2.
    let gas = air_at(293.15);
    let radius = 0.03;
    let duct = Duct { segments: vec![Segment::Cylinder { radius, length: 0.6 }] };
    let z_in = C64::from_re(gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius));
    let reflectance = |z: C64| (z - z_in) / (z + z_in);
    let mut errors = Vec::new();
    for ka in [0.02, 0.05, 0.1, 0.2, 0.4] {
        let omega = ka * gas.sound_speed / radius;
        let full = outlet::mouth_response_flanged(&duct, &gas, &gas, omega, LossModel::Lossless, 8, 1).unwrap();
        let today = outlet::mouth_response(&duct, &gas, omega, LossModel::Lossless,
            Termination::FlangedOpen, 8, 1).unwrap();
        let error = (reflectance(full.input.plane_impedance) - reflectance(today.input.plane_impedance)).abs();
        eprintln!("FLANGED_PLANE_LIMIT ka={ka} error={error:e} bound={:e}", 2.0 * ka * 1e-2);
        assert!(error <= 2.0 * ka * 1e-2, "ka {ka}: {error:e}");
        errors.push(error);
    }
    assert!(errors.windows(2).all(|pair| pair[1] > pair[0]), "the difference must shrink with ka: {errors:?}");
    assert!(errors[0] < 0.2 * errors[4], "{errors:?}");
}

#[test]
fn flanged_mouth_fir_reproduces_the_coupled_mouth_flow_between_dft_bins() {
    // Same property as the analytic outlet FIR test, for the complete radiation matrix and every retained mode.
    let gas = air_at(293.15);
    let duct = flare(0.05);
    let line = MmLineConfig { sample_rate_hz: 48000, n_modes: 3, extra_slices: 1, loss: LossModel::Bessel };
    let reference = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.012 * 0.012);
    // Where the refusal used to fire: admission and the realized path both accept the flanged load now.
    let output = MouthFlowConfig::default();
    realize::admit(core::slice::from_ref(&duct), &["flare"], &gas, MmBankLoad::Shared(MmLoad::FlangedModal(&gas)),
        &line, output).unwrap();
    // Mode zero (total flow) is the observed one: it meets the 1e-10 / 1e-5 gates on this flare, while the higher
    // modes' own time images do not at 65536 taps (measured 5e-7 and 9e-7), and no on-axis receiver sees them.
    let (filters, receipt) = realize::path(&duct, "flare", &gas, MmLoad::FlangedModal(&gas), &line, reference,
        output, 1).unwrap();
    assert_eq!(filters.len(), 1);
    assert!(realize::path(&duct, "flare", &gas, MmLoad::FlangedModal(&gas), &line, reference, output, 0).is_err());
    assert!(realize::path(&duct, "flare", &gas, MmLoad::FlangedModal(&gas), &line, reference, output, 4).is_err());
    assert!(receipt.rungs.len() >= 2);
    // The staircase's actual last station (the flare's mouth side, not the 12 mm inlet), never beyond the cone's outlet.
    assert!(receipt.mouth_radius_m > 0.03 && receipt.mouth_radius_m <= 0.05, "{}", receipt.mouth_radius_m);
    let last = receipt.rungs.last().unwrap();
    assert!(last.tail_energy.iter().all(|x| *x < 1e-10) && last.tap_change.iter().all(|x| *x < 1e-5));
    let dt = 1.0 / f64::from(line.sample_rate_hz);
    let frequencies = [1.0, 7.0, 310.0, 1013.0, 2207.0, 3119.0, 9100.0];
    let direct: Vec<Vec<C64>> = frequencies.iter().map(|f| outlet::mouth_response_flanged(&duct, &gas, &gas,
        core::f64::consts::TAU * f, LossModel::Bessel, 3, 1).unwrap()
        .per_incident_plane_pressure(reference).unwrap()).collect();
    for (mode, filter) in filters.iter().enumerate() {
        let (_, taps) = filter.quiescent_fir().unwrap();
        let scale = direct.iter().map(|row| row[mode].abs()).fold(0.0, f64::max);
        assert!(scale > 0.0, "mode {mode} carries no flow");
        let mut worst = 0.0_f64;
        for (f, row) in frequencies.iter().zip(&direct) {
            let omega = core::f64::consts::TAU * f;
            let (mut re, mut im) = (0.0, 0.0);
            for (n, h) in taps.iter().enumerate() {
                let phase = omega * n as f64 * dt;
                re += h * phase.cos();
                im += h * phase.sin();
            }
            worst = worst.max((C64::new(re, im) - row[mode]).abs());
        }
        eprintln!("FLANGED_MOUTH_FIR mode={mode} taps={} worst={worst:e} scale={scale:e} ratio={:e}",
            taps.len(), worst / scale);
        assert!(worst < 5e-5 * scale, "mode {mode}: {worst:e} against scale {scale:e}");
    }
}

#[test]
fn flanged_admission_ceiling_tracks_the_mouth_radius_and_one_static_pressure() {
    // Planted negative: the exterior ka ceiling (FLANGED_MAX_KA = 24 at Nyquist) belongs to the MOUTH radius.
    // The admitted clock must scale as 1/R and sit at 24 c / (pi R); a code path that read the 12 mm inlet radius
    // instead would admit clocks five times too high and fail both checks.
    let gas = air_at(293.15);
    let admitted = |mouth: f64, rate: u32| {
        let line = MmLineConfig { sample_rate_hz: rate, n_modes: 2, extra_slices: 1, loss: LossModel::Lossless };
        realize::admit(&[flare(mouth)], &["flare"], &gas, MmBankLoad::Shared(MmLoad::FlangedModal(&gas)),
            &line, MouthFlowConfig::default()).is_ok()
    };
    let boundary = |mouth: f64| {
        let (mut lo, mut hi) = (8_000_u32, 192_000_u32);
        assert!(admitted(mouth, lo) && !admitted(mouth, hi));
        while hi - lo > 1 {
            let mid = lo + (hi - lo) / 2;
            if admitted(mouth, mid) { lo = mid } else { hi = mid }
        }
        lo
    };
    let (small, large) = (boundary(0.05), boundary(0.06));
    for (mouth, found) in [(0.05, small), (0.06, large)] {
        let predicted = 24.0 * gas.sound_speed / (core::f64::consts::PI * mouth);
        assert!((f64::from(found) - predicted).abs() <= 1.0, "R {mouth}: {found} vs {predicted}");
    }
    assert!((f64::from(small) / f64::from(large) - 1.2).abs() < 1e-3, "{small} / {large}");
    let other_pressure = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 90_000.0).unwrap();
    let line = MmLineConfig { sample_rate_hz: 32_000, n_modes: 2, extra_slices: 1, loss: LossModel::Lossless };
    assert!(realize::admit(&[flare(0.05)], &["flare"], &gas, MmBankLoad::Shared(MmLoad::FlangedModal(&other_pressure)),
        &line, MouthFlowConfig::default()).is_err(), "interface needs one static pressure");
    // The other refusals of the same function still fire for the flanged load, and for the loads it still refuses.
    let flanged = MmBankLoad::Shared(MmLoad::FlangedModal(&gas));
    assert!(realize::admit(&[flare(0.05)], &["flare"], &gas, flanged, &line, MouthFlowConfig::default()).is_ok());
    for config in [MouthFlowConfig { max_total_taps: 0, ..MouthFlowConfig::default() },
        MouthFlowConfig { max_taps: 300, ..MouthFlowConfig::default() },
        MouthFlowConfig { tail_tolerance: 0.0, ..MouthFlowConfig::default() }] {
        assert!(realize::admit(&[flare(0.05)], &["flare"], &gas, flanged, &line, config).is_err());
    }
    let too_many = MmLineConfig { sample_rate_hz: 32_000, n_modes: 9, extra_slices: 1, loss: LossModel::Lossless };
    assert!(realize::admit(&[flare(0.05)], &["flare"], &gas, flanged, &too_many, MouthFlowConfig::default())
        .is_err(), "more modes than the staircase image keeps");
    assert!(realize::admit(&[flare(0.05), flare(0.05)], &["a"], &gas, flanged, &line,
        MouthFlowConfig::default()).is_err(), "one label per path");
    assert!(realize::admit(&[flare(0.05)], &["flare"], &gas, MmBankLoad::Shared(MmLoad::Matched), &line,
        MouthFlowConfig::default()).is_err(), "a matched continuation is still refused");
    // A bore too short for its clock still refuses before any load is evaluated.
    let stub = Duct { segments: vec![Segment::Cylinder { radius: 0.02, length: 1e-4 }] };
    assert!(realize::admit(&[stub], &["stub"], &gas, flanged, &line, MouthFlowConfig::default()).is_err());
}
