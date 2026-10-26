//! The one-way skin receiver on the original closed shell. The shell's own outward triangles are its skin, each vertex
//! anchored at its own node, so every anchored point is the node itself and the receiver replays the accepted body modes.
use super::shell_observer::shell_pack;
use super::*;

const SHELL: &str = include_str!("../closed-shell.body");
const BARE: &str = include_str!("../closed-shell-cavity.material-strings");
const RECORD: &str = "one_way_observer front 0.2 0.2 1.0 0.01 120 400 256 4 1024 2048 4096 8192 1e-10 1e-4 65536 65536 triangle 1048576 80 1e-9 1000000000000 4 0.1 0.1\n";

fn own_skin() -> String {
    let nodes: Vec<Vec<&str>> = SHELL
        .lines()
        .filter(|l| l.starts_with("node "))
        .map(|l| l.split_ascii_whitespace().skip(1).take(3).collect())
        .collect();
    let triangles: Vec<Vec<usize>> = SHELL
        .lines()
        .filter(|l| l.starts_with("triangle "))
        .map(|l| {
            l.split_ascii_whitespace()
                .skip(1)
                .take(3)
                .map(|v| v.parse().unwrap())
                .collect()
        })
        .collect();
    let mut skin = format!("acoustic-skin-v1\nskin-vertices {}\n", nodes.len());
    for (node, p) in nodes.iter().enumerate() {
        let (triangle, corner) = triangles
            .iter()
            .enumerate()
            .find_map(|(t, tri)| tri.iter().position(|&v| v == node).map(|c| (t, c)))
            .unwrap();
        let mut b = [0; 3];
        b[corner] = 1;
        skin += &format!(
            "skin-vertex {} {} {} 1\nskin-anchor 1 0 {triangle} {} {} {}\n",
            p[0], p[1], p[2], b[0], b[1], b[2]
        );
    }
    skin += &format!("skin-triangles {}\n", triangles.len());
    for t in &triangles {
        skin += &format!("skin-triangle {} {} {}\n", t[0], t[1], t[2]);
    }
    skin + "skin-apertures 0\n"
}

fn source(owner: &str) -> String {
    BARE.replace("audio 48000 1024", "audio 8000 1024")
        .replace("cavity none\n", &format!("{owner}{RECORD}cavity none\n"))
}

fn profile(skin: &str) -> body::Profile {
    body::Profile::parse(format!("{}\n{skin}", SHELL.trim_end()).as_bytes()).unwrap()
}

fn skinned(text: &str, cx: &Cx<'_>, gate: &CancelGate) -> Ready {
    input::Spec::parse(text.as_bytes())
        .unwrap()
        .compile_with_body(cx, gate, Some(profile(&own_skin())), |name| {
            Ok(if name == "closed-shell.fsmatpk" {
                shell_pack(name)
            } else {
                pack(name, 7800.0, 200_000.0, None)
            })
        })
        .unwrap()
}

#[test]
fn one_way_listener_replays_the_accepted_body_modes_through_the_shell_s_own_skin() {
    scoped(|cx, gate| {
        let mut ready = skinned(&source(""), cx, gate);
        assert!(ready.performance.material_shell_observer().is_none());
        let report = ready.performance.one_way_shell_observer().unwrap().clone();
        assert_eq!(report.gas, ready.spec.gas);
        assert_eq!(report.body_modes, ready.spec.body.modes);
        assert_eq!((report.surface_panels, report.apertures), (80, 0));
        assert_eq!(report.largest_anchor_distance_m, 0.0);
        assert_eq!(report.modes.len(), report.body_modes);
        assert!(report.boundary_terms <= report.options.maximum_boundary_terms);
        assert!(report.replay_terms <= report.options.maximum_replay_terms);
        assert_eq!(report.transfer().filter_latency_samples, 256);
        assert!(
            report
                .transfer()
                .geometric_delays_s
                .iter()
                .all(|delay| *delay > 0.0)
        );
        for mode in &report.modes {
            assert!(
                mode.radiation_zeta > 0.0 && mode.added_mass_shift.is_finite(),
                "{mode:?}"
            );
        }
        let mut reference = report.transfer().observer().unwrap();
        let body = ready.spec.strings.len();
        let (mut expected, mut actual) = (Vec::new(), Vec::new());
        for _ in 0..ready.spec.samples {
            let mut output = [0.0];
            ready.performance.block(&mut output).unwrap();
            let velocity = ready.performance.bank().system().components()[body]
                .states()
                .iter()
                .map(|state| state.velocity_m_sqrt_kg_per_s)
                .collect::<Vec<_>>();
            expected.push(reference.step(&velocity).unwrap()[0]);
            actual.push(output[0]);
        }
        assert_eq!(
            actual, expected,
            "the receiver must consume the accepted body modes once per sample"
        );
        assert!(actual.iter().any(|pressure| pressure.abs() > 1e-14));
        // Prepared-linear keeps the same modes, so the same skin gives the same receipt.
        let prepared = skinned(
            &source("prepared-linear 512 4294967296 1000000000000000\n"),
            cx,
            gate,
        );
        assert!(
            prepared
                .performance
                .bank()
                .source()
                .bank()
                .prepared_linear_budget()
                .is_some()
        );
        assert_eq!(
            prepared.performance.one_way_shell_observer().unwrap().modes,
            report.modes
        );
    });
}

#[test]
fn one_way_source_refuses_loaded_listeners_unmatched_skins_and_bad_records_before_material_io() {
    let one_way = input::Spec::parse(source("").as_bytes()).unwrap();
    let options = one_way.one_way_observer.as_ref().unwrap();
    assert_eq!(
        (
            options.receiver.id.as_str(),
            options.workers,
            options.maximum_skin_panels
        ),
        ("front", 4, 80)
    );
    assert_eq!(
        options.panel_rule,
        fs_bem::helmholtz::PanelRule::Triangle {
            max_pairs: 1_048_576
        }
    );
    assert!(
        input::Spec::parse(BARE.as_bytes())
            .unwrap()
            .one_way_observer
            .is_none()
    );
    for loaded in [
        include_str!("../closed-shell-pluck.material-strings"),
        include_str!("../closed-shell-observed.material-strings"),
    ] {
        let text = loaded.replace("cavity none\n", &format!("{RECORD}cavity none\n"));
        assert!(
            input::Spec::parse(text.as_bytes()).is_err(),
            "a loaded listener and a one-way one together"
        );
    }
    for (from, to) in [
        ("triangle 1048576", "square 1048576"),
        ("triangle 1048576", "centroid"),
        ("triangle 1048576", "triangle 0"),
        (" 4 0.1 0.1\n", " 0 0.1 0.1\n"),
        ("0.1 0.1\n", "0.1 NaN\n"),
        ("0.1 0.1\n", "0.1 0.1 7\n"),
        (" 80 1e-9", " 9000 1e-9"),
        ("1e-10 1e-4", "1e-9 1e-4"),
        ("0.01 120 400 256", "0.01 120 400 0"),
        ("front 0.2 0.2 1.0", "front NaN 0.2 1.0"),
    ] {
        assert!(
            input::Spec::parse(source("").replace(from, to).as_bytes()).is_err(),
            "{to}"
        );
    }
    let skin = own_skin();
    for broken in [
        skin.replace("skin-apertures 0\n", "skin-apertures 0\nskin-vertices 1\n"),
        skin.replacen("skin-triangle 0 ", "skin-triangle 42 ", 1),
        skin.replacen("skin-vertex ", "skin-vertex NaN ", 1),
    ] {
        assert!(
            body::Profile::parse(format!("{}\n{broken}", SHELL.trim_end()).as_bytes()).is_err()
        );
    }
    scoped(|cx, gate| {
        let cases = [
            (source(""), None),
            (source(""), Some(profile(""))),
            (
                source("").replace(" 80 1e-9", " 79 1e-9"),
                Some(profile(&skin)),
            ),
            (
                source("").replace("cavity none", "cavity neck 0.01 0.01 0.01 0 1"),
                Some(profile(&skin)),
            ),
        ];
        for (text, body) in cases {
            let mut reads = 0;
            assert!(
                input::Spec::parse(text.as_bytes())
                    .unwrap()
                    .compile_with_body(cx, gate, body, |_| {
                        reads += 1;
                        Err("material reads are forbidden on an unmatched one-way skin".into())
                    })
                    .is_err()
            );
            assert_eq!(reads, 0);
        }
    });
}

#[test]
fn a_kept_setup_renders_the_computed_one_bit_for_bit_and_a_damaged_entry_refuses() {
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir =
        std::env::temp_dir().join(format!("fs-one-way-setups-{}-{stamp}", std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    scoped(|cx, gate| {
        let kept = |text: &str| {
            let mut body = profile(&own_skin());
            body.keep_observer_setups(&dir).unwrap();
            input::Spec::parse(text.as_bytes())
                .unwrap()
                .compile_with_body(cx, gate, Some(body), |name| {
                    Ok(if name == "closed-shell.fsmatpk" {
                        shell_pack(name)
                    } else {
                        pack(name, 7800.0, 200_000.0, None)
                    })
                })
        };
        let render = |ready: &mut Ready| {
            (0..ready.spec.samples)
                .map(|_| {
                    let mut output = [0.0];
                    ready.performance.block(&mut output).unwrap();
                    output[0].to_bits()
                })
                .collect::<Vec<_>>()
        };
        let mut plain = skinned(&source(""), cx, gate);
        let mut computed = kept(&source("")).unwrap();
        let mut loaded = kept(&source("")).unwrap();
        let reports = [&plain, &computed, &loaded]
            .map(|ready| ready.performance.one_way_shell_observer().unwrap().clone());
        assert_eq!(
            reports.each_ref().map(|r| r.setup_loaded),
            [false, false, true]
        );
        let kernels = reports.each_ref().map(|r| {
            r.transfer()
                .kernels
                .iter()
                .map(|k| {
                    let taps: Vec<u64> = k.taps.iter().map(|t| t.to_bits()).collect();
                    (k.input, k.output, k.integer_delay_samples, taps)
                })
                .collect::<Vec<_>>()
        });
        for (report, bits) in reports.iter().zip(&kernels).skip(1) {
            assert_eq!(report.setup_key, reports[0].setup_key);
            assert_eq!(report.modes, reports[0].modes);
            assert_eq!(*bits, kernels[0]);
        }
        let entry = computed
            .body_profile
            .as_ref()
            .unwrap()
            .observer_store_entry(&reports[0].setup_key)
            .unwrap();
        assert_eq!(
            std::fs::read_dir(&dir).unwrap().count(),
            1,
            "one setup, one entry"
        );
        let expected = render(&mut plain);
        assert!(expected.iter().any(|&p| f64::from_bits(p).abs() > 1e-14));
        assert_eq!(
            render(&mut computed),
            expected,
            "keeping a setup changes no sample"
        );
        assert_eq!(
            render(&mut loaded),
            expected,
            "a loaded setup renders the computed one bit for bit"
        );
        // Another receiver is another setup, under its own key.
        let moved = kept(&source("").replace("front 0.2 0.2 1.0", "front 0.2 0.3 1.0")).unwrap();
        let other = moved.performance.one_way_shell_observer().unwrap();
        assert!(!other.setup_loaded && other.setup_key != reports[0].setup_key);
        let second = moved
            .body_profile
            .as_ref()
            .unwrap()
            .observer_store_entry(&other.setup_key)
            .unwrap();
        // One flipped byte refuses, and nothing is recomputed behind it.
        let mut bytes = std::fs::read(&entry).unwrap();
        let middle = bytes.len() / 2;
        bytes[middle] ^= 1;
        std::fs::write(&entry, &bytes).unwrap();
        assert!(kept(&source("")).is_err(), "a damaged kept setup");
        assert_eq!(
            std::fs::read(&entry).unwrap(),
            bytes,
            "a damaged entry is left as it is"
        );
        for path in [&entry, &second] {
            std::fs::remove_file(path).unwrap();
        }
        std::fs::remove_dir(&dir).unwrap();
    });
}
