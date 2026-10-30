//! G1/G3: native molecular input reaches physical side-mouth pressure and replay.
//! Reservoir coefficients are explicitly synthetic, not measured N2/O2 data.
use fs_couple::pcm_wav::observation::ObservationPolicy;
use fs_couple::render::schedule::reed::{ReedPerformance, ReedPerformanceError};

const OPENINGS: &str = include_str!("../examples/reed-molecular-openings-synthetic.performance");
const STRAIGHT: &str = include_str!("../examples/reed-molecular-synthetic.performance");

fn input(source: &str) -> String {
    source
        .replace("audio 48000 9606", "audio 8000 1601")
        .replace("diffusive 32 0.0001 160000", "diffusive 16 0.0001 30000")
}

fn performance(source: &str, block: usize) -> ReedPerformance {
    ReedPerformance::from_bytes_with_policy(source.as_bytes(), block, ObservationPolicy::Acceptance)
        .unwrap()
}

fn render(mut source: ReedPerformance, block: usize) -> Vec<f64> {
    let mut pressure = vec![0.0; source.info().samples as usize];
    for samples in pressure.chunks_mut(block) {
        source.pressure_block(samples).unwrap();
    }
    assert_eq!(source.remaining_samples(), 0);
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| *p != 0.0));
    pressure
}

#[test]
fn data_only_open_vent_closed_and_reservoir_changes_reach_physical_pressure() {
    let text = input(OPENINGS);
    let baseline = performance(&text, 37);
    let point = baseline.molecular_driving_point().unwrap();
    assert_eq!(point.open_hole_sites().len(), 2);
    assert_eq!(baseline.molecular_open_hole_observers().count(), 2);
    assert_eq!(point.inlet_gas(), baseline.ambient_gas());
    assert_eq!(
        *point.exterior_gas(),
        baseline.outlet_receiver().unwrap().gas
    );
    let pa = render(baseline, 37);
    let norm: f64 = pa.iter().map(|p| p * p).sum();
    for (changed, mouth_count) in [
        (
            text.replace(
                "hole 0.00015 0.0015 0.001 0.5",
                "hole 0.00015 0.0015 0.001 0",
            ),
            1,
        ),
        (
            text.replace("hole 0.0002 0.001 0.001 1", "hole 0.0002 0.001 0.001 0.25"),
            2,
        ),
        (text.replace("power 2.0265", "power 0.20265"), 2),
    ] {
        let source = performance(&changed, 37);
        assert_eq!(source.molecular_open_hole_observers().count(), mouth_count);
        assert_eq!(
            source
                .molecular_driving_point()
                .unwrap()
                .open_hole_sites()
                .len(),
            mouth_count
        );
        let pb = render(source, 37);
        let delta: f64 = pa.iter().zip(&pb).map(|(a, b)| (a - b).powi(2)).sum();
        assert!(
            delta > 1e-12 * norm,
            "data changes must reach nonzero physical pressure: norm={norm}, delta={delta}"
        );
    }
    assert_eq!(pa, render(performance(&text, 113), 113));
}

#[test]
fn each_mouth_keeps_its_own_receiver_range_and_observation_cannot_change_reed_feedback() {
    let text = input(OPENINGS).replace(
        "baffled_outlet 0 0 1 2 8 800",
        "baffled_outlet 0 0 1 2 8 800\nexterior mixture 320 101325 mole 2 N2 0.7 O2 0.3",
    );
    let filtered_text = text.replace("observer_law none", "observer_law stokes-kirchhoff");
    let mut plain = performance(&text, 37);
    let mut filtered = performance(&filtered_text, 37);
    let point = filtered.molecular_driving_point().unwrap();
    let terminal = filtered.outlet_receiver().unwrap();
    assert_eq!(point.inlet_gas(), filtered.ambient_gas());
    assert_ne!(point.inlet_gas(), point.exterior_gas());
    assert_eq!(point.exterior_gas().temperature, 320.0);
    assert_eq!(*point.exterior_gas(), terminal.gas);
    let mouths: Vec<_> = filtered.molecular_open_hole_observers().collect();
    assert_eq!(mouths.len(), 2);
    for ((receiver, observer), site) in mouths.iter().zip(point.open_hole_sites()) {
        assert_eq!(receiver.gas, terminal.gas);
        assert_eq!(receiver.radius_m, site.hole_radius_m);
        let p = receiver.receiver.position_m;
        let distance = p[0].hypot(p[1]).hypot(p[2]);
        assert_eq!(observer.filter.range_m, distance);
        assert!(observer.filter.range_m > filtered.observer().filter.range_m);
        assert_eq!(observer.filter.sample_rate_hz, 8000);
        assert!(!observer.internal_part_unpropagated);
        assert!(observer.range_spread_m.unwrap() > 0.0);
        assert!(observer.filter.receiver_delay_s > 0.0);
    }
    assert_ne!(
        mouths[0].0.receiver.position_m,
        mouths[1].0.receiver.position_m
    );
    assert_ne!(mouths[0].1.filter.range_m, mouths[1].1.filter.range_m);
    filtered.observer_admission().unwrap();

    let mut plain_pressure = vec![0.0; plain.info().samples as usize];
    let mut filtered_pressure = vec![0.0; plain_pressure.len()];
    for (a, b) in plain_pressure
        .chunks_mut(37)
        .zip(filtered_pressure.chunks_mut(37))
    {
        plain.pressure_block(a).unwrap();
        filtered.pressure_block(b).unwrap();
        let internal = |performance: &ReedPerformance| {
            performance
                .renderer()
                .context()
                .reed_voice(0)
                .unwrap()
                .last_observer_parts()
                .internal_pa
        };
        assert_eq!(internal(&plain).to_bits(), internal(&filtered).to_bits());
    }
    assert!(plain_pressure.iter().all(|x| x.is_finite()));
    assert!(filtered_pressure.iter().all(|x| x.is_finite()));
    assert_ne!(plain_pressure, filtered_pressure);
    assert_eq!(plain.remaining_samples(), 0);
    assert_eq!(filtered.remaining_samples(), 0);
}

#[test]
fn straight_limit_and_explicit_source_refusals_preserve_the_native_contract() {
    let straight = performance(&input(STRAIGHT), 37);
    assert!(
        straight
            .molecular_driving_point()
            .unwrap()
            .open_hole_sites()
            .is_empty()
    );
    assert_eq!(straight.molecular_open_hole_observers().count(), 0);
    let text = input(OPENINGS);
    let regional = text.replacen(
        "cylinder 0.001 0.05",
        "cylinder 0.001 0.05 gas mixture 293.15 101325 mole 2 N2 0.7 O2 0.3",
        1,
    );
    assert!(
        matches!(ReedPerformance::from_bytes(regional.as_bytes(), 37),
        Err(ReedPerformanceError::Input { what, .. }) if what.contains("homogeneous ambient species"))
    );
    for refused in [
        text.replace("segments 5", "segments 4"),
        text.replace("O2 1 vibration", "CO2 1 vibration"),
        text.replace("hole 0.0002 0.001 0.001 1", "hole 0.0002 0.001 0.001 1.1"),
    ] {
        assert!(ReedPerformance::from_bytes(refused.as_bytes(), 37).is_err());
    }
}
