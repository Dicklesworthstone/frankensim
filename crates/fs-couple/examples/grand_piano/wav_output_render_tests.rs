//! Direct WAV publication from the actual hammer/string/board pressure stream.
use super::*;

fn u16_at(bytes: &[u8], offset: usize) -> u16 {
    u16::from_le_bytes(bytes[offset..offset + 2].try_into().unwrap())
}

fn u32_at(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(bytes[offset..offset + 4].try_into().unwrap())
}

#[test]
fn wav_formats_preserve_physical_mono_and_stereo_pressure_and_quiet_samples() {
    // G3: only the storage encoding changes. This uses the same illustrative
    // board/surface fixture as the existing stereo tests, not another engine.
    let course = geometry::demonstration_scale().unwrap()[48];
    let modes = board::demonstration();
    let surface: Vec<_> = [-0.25, 0.25].into_iter().enumerate().map(|(side, x)| {
        board_geometry::SurfaceSample {
            position_m: [x, 0.0, 0.0],
            area_m2: 0.1,
            mode_shape: modes.iter().enumerate().map(|(i, mode)| {
                mode.volume / 0.2 * if (i + side) % 2 == 0 { 1.0 } else { -0.5 }
            }).collect(),
        }
    }).collect();
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)
        .unwrap().as_nanos();
    let directory = std::env::temp_dir().join(format!(
        "frankensim-wav-format-render-{}-{stamp}", std::process::id(),
    ));
    std::fs::create_dir(&directory).unwrap();
    let frames = 1400;
    let rate = 48_000;
    // One declared scale for every encoding and both receivers. The large
    // headroom makes quiet quantization observable without changing mechanics.
    let full_scale_pa = 64.0;
    let mut mono_pressure = Vec::new();
    for channels in [1, 2] {
        let mut reference_csv = None;
        let mut pressure = Vec::new();
        let mut decoded = Vec::new();
        for format in [wav_output::Format::Pcm16, wav_output::Format::Pcm24,
            wav_output::Format::Float32] {
            let stem = format!("{channels}-{}", format.label());
            let wav_path = directory.join(format!("{stem}.wav"));
            let csv_path = directory.join(format!("{stem}.csv"));
            let options = Options {
                render: Some(wav_path.to_str().unwrap().to_owned()),
                note: Some(69), velocity: Some(2.0),
                duration: frames as f64 / f64::from(rate),
                sample_rate: rate, modes: 8, wav_format: format,
                pcm_full_scale_pa: full_scale_pa,
                microphone: Some([-0.4, 0.1, 0.8]),
                microphone_right: (channels == 2).then_some([0.65, 0.15, 1.2]),
                receiver_pressure_csv: Some(csv_path.to_str().unwrap().to_owned()),
                ..Options::default()
            };
            render(wav_path.to_str().unwrap(), vec![course], &modes,
                Some(&surface), &options).unwrap();
            let csv = std::fs::read_to_string(&csv_path).unwrap();
            if let Some(reference) = &reference_csv {
                assert_eq!(&csv, reference, "encoding changed the physical pressure trace");
            } else {
                let mut rows = csv.lines();
                assert_eq!(rows.next().unwrap(), if channels == 1 {
                    "sample,time_s,pressure_left_pa"
                } else { "sample,time_s,pressure_left_pa,pressure_right_pa" });
                for (sample, row) in rows.enumerate() {
                    let columns: Vec<_> = row.split(',').collect();
                    assert_eq!(columns.len(), channels + 2);
                    assert_eq!(columns[0].parse::<usize>().unwrap(), sample);
                    assert_eq!(columns[1].parse::<f64>().unwrap(), sample as f64 / f64::from(rate));
                    pressure.extend(columns[2..].iter().map(|p| p.parse::<f64>().unwrap()));
                }
                assert_eq!(pressure.len(), frames * channels);
                assert!(pressure.iter().all(|p| p.is_finite() && p.abs() < full_scale_pa));
                reference_csv = Some(csv);
            }

            let wav = std::fs::read(&wav_path).unwrap();
            let (data_offset, bytes_per_sample, quantizer) = match format {
                wav_output::Format::Pcm16 => (44, 2, 32_767.0),
                wav_output::Format::Pcm24 => (68, 3, 8_388_607.0),
                wav_output::Format::Float32 => (80, 4, 1.0),
            };
            assert_eq!(&wav[..4], b"RIFF");
            assert_eq!(&wav[8..16], b"WAVEfmt ");
            assert_eq!(u32_at(&wav, 4) as usize + 8, wav.len());
            assert_eq!(u16_at(&wav, 22) as usize, channels);
            assert_eq!(u32_at(&wav, 24), rate);
            assert_eq!(u16_at(&wav, 32) as usize, channels * bytes_per_sample);
            assert_eq!(u16_at(&wav, 34) as usize, 8 * bytes_per_sample);
            if format != wav_output::Format::Pcm16 {
                assert_eq!(u16_at(&wav, 20), 0xfffe);
                assert_eq!(u32_at(&wav, 44), if format == wav_output::Format::Float32 { 3 } else { 1 });
            }
            assert_eq!(&wav[data_offset - 8..data_offset - 4], b"data");
            assert_eq!(u32_at(&wav, data_offset - 4) as usize, pressure.len() * bytes_per_sample);
            assert_eq!(wav.len(), data_offset + pressure.len() * bytes_per_sample);
            let mut samples = Vec::with_capacity(pressure.len());
            for (bytes, p) in wav[data_offset..].chunks_exact(bytes_per_sample).zip(&pressure) {
                let actual = match format {
                    wav_output::Format::Pcm16 => f64::from(i16::from_le_bytes(bytes.try_into().unwrap())) / quantizer,
                    wav_output::Format::Pcm24 => f64::from(i32::from_le_bytes([
                        bytes[0], bytes[1], bytes[2], if bytes[2] & 0x80 != 0 { 0xff } else { 0 },
                    ])) / quantizer,
                    wav_output::Format::Float32 => f64::from(f32::from_le_bytes(bytes.try_into().unwrap())),
                };
                let expected = p / full_scale_pa;
                if format == wav_output::Format::Float32 {
                    assert_eq!(actual, f64::from(expected as f32));
                } else {
                    assert!((actual - expected).abs() <= 0.5 / quantizer + 1e-15);
                }
                samples.push(actual);
            }
            decoded.push(samples);
        }
        for channel in 0..channels {
            assert!((channel..pressure.len()).step_by(channels).any(|i| {
                pressure[i] != 0.0 && decoded[0][i] == 0.0
                    && decoded[1][i] != 0.0 && decoded[2][i] != 0.0
            }), "channel {channel} needs a physical quiet sample lost by PCM16 but retained by both finer encodings");
        }
        if channels == 1 {
            mono_pressure = pressure;
        } else {
            assert!(pressure.chunks_exact(2).any(|frame| frame[0] != frame[1]));
            assert_eq!(pressure.chunks_exact(2).map(|frame| frame[0]).collect::<Vec<_>>(), mono_pressure);
        }
    }
}
