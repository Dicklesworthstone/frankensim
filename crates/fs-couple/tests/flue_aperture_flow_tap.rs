//! The aperture flow tap observes and never acts (bead
//! `frankensim-irp-v5-instrument-realism-pexfp.104`, H4 test 1): the committed
//! published-Boehm flute D4 performance renders bit-identically with and without
//! it, and the tap holds one finite outward volume flow per observed opening per
//! mechanical frame, in the order its header names them.
use fs_couple::flue_loop::network::performance::aperture_flows::APERTURE_FLOWS_SCHEMA;
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_exec::CancelGate;
use std::io::Write;
use std::sync::{Arc, Mutex};

/// 50 ms at the 48 kHz observation clock.
const OUTPUT_SAMPLES: usize = 2400;
const BLOCK: usize = 480;

#[derive(Clone, Default)]
struct Shared(Arc<Mutex<Vec<u8>>>);
impl Write for Shared {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.0.lock().expect("sink").extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

struct Run {
    pressure: Vec<f64>,
    frames: Option<u64>,
    outlets: usize,
    mouth: bool,
    ratio: usize,
}

fn render(tap: Option<Shared>) -> Run {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/instruments/descriptions/flute/flute-d4.performance");
    let bytes = std::fs::read(path).expect("committed flute D4 performance");
    let mut p = FlueFilePerformance::from_bytes(&bytes, 512, &CancelGate::new())
        .expect("flute D4 admits")
        .into_renderer();
    let tapped = tap.is_some();
    if let Some(sink) = tap {
        p = p
            .with_aperture_flow_tap(Box::new(sink))
            .expect("tap attaches");
    }
    let mut pressure = vec![0.0; OUTPUT_SAMPLES];
    for chunk in pressure.chunks_mut(BLOCK) {
        p.block(chunk).expect("flute renders");
    }
    if tapped {
        p.flush_aperture_flow_tap().expect("tap flushes");
    }
    Run {
        pressure,
        frames: p.aperture_flow_frames(),
        outlets: p.outlet_nodes().len(),
        mouth: p.mouth_radiation().is_some(),
        ratio: p.mechanics_per_output(),
    }
}

#[test]
fn aperture_flow_tap_leaves_the_flute_render_bit_identical() {
    let plain = render(None);
    assert_eq!(plain.frames, None);
    assert!(
        plain.pressure.iter().any(|&p| p != 0.0),
        "the window is not silent"
    );

    let sink = Shared::default();
    let tapped = render(Some(sink.clone()));
    let identical = plain
        .pressure
        .iter()
        .zip(&tapped.pressure)
        .all(|(a, b)| a.to_bits() == b.to_bits());
    assert!(identical, "the tap changed the microphone pressure");

    let bytes = sink.0.lock().expect("sink").clone();
    let split = bytes
        .windows(6)
        .position(|w| w == b"\ndata\n")
        .expect("header ends in a data line")
        + 6;
    let header = std::str::from_utf8(&bytes[..split]).expect("text header");
    assert!(header.starts_with(APERTURE_FLOWS_SCHEMA));
    let apertures = header
        .lines()
        .filter(|l| l.starts_with("aperture "))
        .count();
    assert_eq!(
        apertures,
        tapped.outlets + usize::from(tapped.mouth),
        "{header}"
    );
    assert!(tapped.mouth, "the flute file observes its embouchure");

    let frames = (OUTPUT_SAMPLES * tapped.ratio) as u64;
    assert_eq!(tapped.frames, Some(frames));
    let data = &bytes[split..];
    assert_eq!(data.len() as u64, frames * apertures as u64 * 8);
    let flows: Vec<f64> = data
        .chunks_exact(8)
        .map(|b| f64::from_le_bytes(b.try_into().expect("8 bytes")))
        .collect();
    assert!(flows.iter().all(|q| q.is_finite()));
    for a in 0..apertures {
        assert!(
            flows.iter().skip(a).step_by(apertures).any(|&q| q != 0.0),
            "aperture {a} carries flow"
        );
    }
}
