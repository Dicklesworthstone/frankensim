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
