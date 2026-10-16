//! Actual supplied geometry/card preparation and the main piano render path.
//! The soft panel and sealed box are authored integration fixtures.
use super::*;
use std::{fmt::Write as _, io::Write as _, path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering}};

fn fresh_path(extension: &str) -> PathBuf {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)
        .unwrap().as_nanos();
    std::env::temp_dir().join(format!("fs-piano-cavity-{}-{stamp}-{}.{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed), extension))
}

fn input(text: &str, extension: &str) -> PathBuf {
    let path = fresh_path(extension);
    std::fs::OpenOptions::new().create_new(true).write(true).open(&path).unwrap()
        .write_all(text.as_bytes()).unwrap();
    path
}

fn panel() -> String {
    let mut text = String::from("frankensim-board-geometry-si-v1\nsource,estimated,three-DOF cavity render panel\nsupport,clamped\npretension,0\ndamping,0.01\n");
    for (i, p) in [[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0],[0.43,0.54]].iter().enumerate() {
        writeln!(text, "node,{i},{},{}", p[0], p[1]).unwrap();
    }
    for (i, triangle) in [[0,1,4],[1,2,4],[2,3,4],[3,0,4]].iter().enumerate() {
        writeln!(text, "triangle,{i},{},{},{},0.08,450,1e6,8e4,0.3,6e4,0.27",
            triangle[0], triangle[1], triangle[2]).unwrap();
    }
    text.push_str("fixed,0\nfixed,1\nfixed,2\nfixed,3\nbridge,69,0,0,0,1\n");
    text
}

fn pressure(path: &Path, frames: usize, rate: u32) -> Vec<[f64; 2]> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some("sample,time_s,pressure_left_pa,pressure_right_pa"));
    let values: Vec<_> = lines.enumerate().map(|(sample, line)| {
        let columns: Vec<f64> = line.split(',').map(|v| v.parse().unwrap()).collect();
        assert_eq!(columns.len(), 4);
        assert!(columns.iter().all(|v| v.is_finite()));
        assert_eq!(columns[0], sample as f64);
        assert_eq!(columns[1], sample as f64 / f64::from(rate));
        [columns[2], columns[3]]
    }).collect();
    assert_eq!(values.len(), frames);
    values
}

#[test]
fn supplied_cavity_changes_rendered_stereo_pressure_for_full_p1_and_reduced_cubic_boards() {
    let board_text = panel();
    let board_path = input(&board_text, "fsb");
    let course = geometry::demonstration_scale().unwrap()[48];
    let scale_text = geometry::write_scale(&[course]);
    let scale_path = input(&scale_text, "csv");
    let cavity_path = input("frankensim-piano-cavity-si-v1\nsource,estimated,authored sealed box\ninterface-origin-m,0,0,0\ndimensions-m,1,1,0.1\nmodes,2\ndamping-ratio,0.02\ngas,dry-air-ussa1976,293.15,101325\n", "fspc");
    for reduced_cubic in [false, true] {
        let loaded_wav = fresh_path("wav");
        let bare_wav = fresh_path("wav");
        let loaded_csv = fresh_path("csv");
        let bare_csv = fresh_path("csv");
        let mut args: Vec<String> = ["--scale", scale_path.to_str().unwrap(),
            "--board-geometry", board_path.to_str().unwrap(), "--render", loaded_wav.to_str().unwrap(),
            "--cavity", cavity_path.to_str().unwrap(), "--note", "69", "--velocity", "0.5",
            "--duration", "0.025", "--modes", "8", "--board-band-hz", "100",
            "--equilibrate-board-mass", "--microphone", "0.3,0.7,1.2",
            "--microphone-right", "0.8,0.25,1.0", "--receiver-pressure-csv", loaded_csv.to_str().unwrap()]
            .iter().map(|s| s.to_string()).collect();
        if reduced_cubic {
            args.extend(["--edge-cubic-board-mass", "--board-reduction", "2,1,10,30"]
                .iter().map(|s| s.to_string()));
        }
        let mut options = Options::parse(&args).unwrap();
        let scale = selected_scale(Some(&scale_text), &options).unwrap();
        let board = prepare_geometric_board_with_reduction(&board_text, &[69], options.board_band_hz,
            true, false, reduced_cubic, 0, true, options.board_reduction.as_ref()).unwrap();
        assert_eq!(board.motion.as_ref().unwrap().is_edge_cubic(), reduced_cubic);
        if reduced_cubic {
            assert!(board.reduction.as_ref().unwrap().source_modes > board.modes.len());
            assert!(board.physical_damping.is_some());
        }
        let enclosed = cavity::Specification::load(cavity_path.to_str().unwrap()).unwrap()
            .project(&board).unwrap();
        render_with_cavity(loaded_wav.to_str().unwrap(), scale.clone(), &board.modes,
            Some(&board.surface), &options, None, None, board.physical_damping.as_deref(),
            Some(&enclosed)).unwrap();
        options.cavity = None;
        options.render = Some(bare_wav.to_str().unwrap().to_owned());
        options.receiver_pressure_csv = Some(bare_csv.to_str().unwrap().to_owned());
        render_with_cavity(bare_wav.to_str().unwrap(), scale, &board.modes,
            Some(&board.surface), &options, None, None, board.physical_damping.as_deref(), None).unwrap();
        let frames = (options.duration * f64::from(options.sample_rate)).round() as usize;
        let loaded = pressure(&loaded_csv, frames, options.sample_rate);
        let bare = pressure(&bare_csv, frames, options.sample_rate);
        assert!(loaded.iter().flatten().any(|p| p.abs() > 1e-9));
        assert!(loaded.iter().any(|p| (p[0] - p[1]).abs() > 1e-12));
        let difference = loaded.iter().flatten().zip(bare.iter().flatten())
            .map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        assert!(difference > 1e-10, "cavity was lost before pressure rendering: reduced_cubic={reduced_cubic}");
        for path in [&loaded_wav, &bare_wav] {
            let wav = std::fs::read(path).unwrap();
            assert_eq!(&wav[..4], b"RIFF");
            assert_eq!(&wav[8..12], b"WAVE");
            assert_eq!(u16::from_le_bytes([wav[22], wav[23]]), 2);
            assert_eq!(wav.len(), 44 + frames * 4);
        }
    }
}

#[test]
fn supplied_cavity_cannot_disappear_at_the_render_front_door() {
    for args in [vec!["--cavity", "box.fspc"],
        vec!["--render", "p.wav", "--cavity", "box.fspc"],
        vec!["--render", "p.wav", "--board", "modal.csv", "--cavity", "box.fspc"],
        vec!["--render", "p.wav", "--board-geometry", "panel.fsb", "--cavity", ""],
        vec!["--render", "p.wav", "--board-geometry", "panel.fsb", "--cavity", "box.fspc", "--cavity", "second.fspc"]] {
        assert!(Options::parse(&args.iter().map(|s| s.to_string()).collect::<Vec<_>>()).is_err());
    }
    let mut options = Options::default();
    options.cavity = Some("supplied.fspc".into());
    let path = fresh_path("wav");
    let course = geometry::demonstration_scale().unwrap()[48];
    let error = render_with_cavity(path.to_str().unwrap(), vec![course], &board::demonstration(),
        None, &options, None, None, None, None).unwrap_err();
    assert!(error.contains("cavity was not admitted and projected"));
    assert!(!path.exists());
}
