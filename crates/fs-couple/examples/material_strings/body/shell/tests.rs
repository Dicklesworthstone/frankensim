//! G3/G4: complete source data -> curved shell/string/gas mechanics -> pressure.
use super::*;
use super::super::shell::{MAX_MODES, MAX_NODES, MAX_SHELL_PROJECTION_WORK, MAX_TRIANGLES};

fn shell_text() -> &'static str { include_str!("../../arched.body") }
fn shell_profile() -> Profile { Profile::parse(shell_text().as_bytes()).unwrap() }
fn shell_spec() -> input::Spec {
    let bytes = include_str!("../../cylinder-flow.material-strings")
        .replace("48000 4096 0.01 128", "48000 96 0.01 128")
        .replace("cylinder-gas-loading 0.001 30 12", "cylinder-gas-loading 0.001 30 4");
    input::Spec::parse(bytes.as_bytes()).unwrap()
}
fn ready(cx: &Cx<'_>, gate: &CancelGate, loss: f64, shear: f64) -> Ready {
    shell_spec().compile_with_body(cx, gate, Some(shell_profile()),
        |name| Ok(if name == "body.fsmatpk" { pack_with_shear(name, loss, shear, None) }
            else { pack(name, 1.0, None) })).unwrap()
}

#[test]
fn curved_material_body_keeps_layers_planes_torsion_cavity_and_gas_in_actual_pressure() {
    scoped(|cx, gate| {
        let mut a = ready(cx, gate, 1.0, 1.0);
        let source = a.performance.bank().source().bank();
        let shell = source.material_shell().unwrap();
        let image = source.material_shell_image().unwrap();
        assert_eq!(a.body_profile.as_ref().unwrap().bytes, shell_text().as_bytes());
        assert_eq!(shell.specimen.mesh().nodes.len(), 25);
        assert_eq!(shell.specimen.mesh().tris.len(), 32);
        assert_eq!(shell.specimen.mesh().nodes[12], [0.2, 0.15, 0.006]);
        assert_eq!(shell.specimen.regions()[0].input().axis1_global, [1.0, 0.0, 0.0]);
        assert!(shell.specimen.mass_kg() > 0.0);
        assert_eq!(shell.force_ports.len(), 1);
        assert_eq!(shell.force_ports[0].len(), 150);
        assert_eq!(shell.force_ports[0][6*12+2], 1.0);
        assert_eq!(image.physical_modes().len(), 2);
        assert!(image.physical_modes().iter().all(|v| v.len() == 150));
        assert!(image.volume_shapes_m2_per_sqrt_kg().iter().any(|v| *v != 0.0));
        assert_eq!(image.damping().mode_count(), 2);
        let MaterialStringSource::WithPolarizations(planes) = &source.strings()[0].source else {
            panic!("two physical string planes disappeared");
        };
        assert!(matches!(planes.source(), MaterialStringSource::WithTorsion(_)));
        assert_eq!(source.cylinder_gas_loading().unwrap().strings()[0].flow().branches().len(), 4);
        assert_eq!(source.cylinder_gas_loading().unwrap().strings()[0].flow().gas(), source.gas());
        assert_eq!(source.material_cavity().unwrap().cavity.gas(), source.gas());
        assert!(a.performance.preload_energy_j().unwrap() > 0.0);

        let samples = pressure(&mut a, 7);
        assert!(samples.iter().any(|p| p.abs() > 1e-14));
        assert!(a.performance.bank().system().components()[1].states().iter()
            .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0));
        let frame = a.performance.bank().system().last_frame().unwrap();
        assert_eq!(frame.external_work_j, 0.0);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);

        let mut b = ready(cx, gate, 1.0, 1.0);
        let encoded = fs_couple::pcm_wav::encode_pcm16_wav(&samples, b.spec.rate, b.spec.full_scale).unwrap().0;
        let (wav, report) = b.write(std::io::Cursor::new(Vec::<u8>::new()), gate).unwrap();
        assert_eq!(wav.into_inner(), encoded);
        assert_eq!(report.samples, 96);

        // Change only sourced transverse shear, then only sourced viscosity.
        // The original material string, its gas and its installed load persist.
        let mut softer = ready(cx, gate, 1.0, 0.1);
        let mut lossy = ready(cx, gate, 4.0, 1.0);
        assert_eq!(a.performance.bank().cards()[0].fundamental_hz(),
            softer.performance.bank().cards()[0].fundamental_hz());
        for name in ["wire.fsmatpk", "sleeve.fsmatpk"] {
            assert_eq!(a.catalog[name].pack.content_hash(), softer.catalog[name].pack.content_hash());
            assert_eq!(a.catalog[name].pack.content_hash(), lossy.catalog[name].pack.content_hash());
        }
        let modes = |r: &Ready| r.performance.bank().system().components()[1].modes().to_vec();
        assert_ne!(modes(&a), modes(&softer));
        assert_eq!(modes(&a), modes(&lossy));
        assert_ne!(a.performance.bank().source().bank().material_shell_image().unwrap()
            .damping().matrix_per_s(), lossy.performance.bank().source().bank()
            .material_shell_image().unwrap().damping().matrix_per_s());
        assert!(samples.iter().zip(pressure(&mut softer, 31)).any(|(a, b)| *a != b));
        assert!(samples.iter().zip(pressure(&mut lossy, 31)).any(|(a, b)| *a != b));
    });
}

#[test]
fn shell_geometry_and_resource_records_refuse_without_borrowing_missing_material_data() {
    let over_nodes = format!("nodes {}", MAX_NODES + 1);
    let over_triangles = format!("triangles {}", MAX_TRIANGLES + 1);
    let over_work = format!("projection-work {}", MAX_SHELL_PROJECTION_WORK + 1);
    for (from, to) in [
        ("nodes 25", over_nodes.as_str()),
        ("triangles 32", over_triangles.as_str()),
        ("mass consistent", "mass automatic"),
        ("director-min-dot 0.9", "director-min-dot 0"),
        ("drilling 0.001", "drilling 0"),
        ("projection-work 67108864", over_work.as_str()),
        ("orthotropic12 1 0.17", "orthotropic12 0 0.17"),
        ("0.8333333333333334 1 0 0", "0.8333333333333334 2 0 0"),
        ("region body @body", "region body ../body.fsmatpk"),
        ("triangle 0 1 6 body", "triangle 0 1 6 unknown"),
        ("support 0 111111", "support 0 11111"),
        ("bridge 12 0 0 1", "bridge 0 0 0 1"),
        ("bridge 12 0 0 1", "bridge 12 0 0 2"),
    ] {
        assert!(Profile::parse(shell_text().replace(from, to).as_bytes()).is_err(), "accepted {to}");
    }
    assert!(Profile::parse(format!("{}extra\n", shell_text()).as_bytes()).is_err());
    scoped(|cx, gate| {
        // Topology, actual physical directors and bridge membership are checked
        // before the first source pack is requested by the executable.
        for text in [
            shell_text().replace("triangle 0 1 6 body", "triangle 0 0 6 body"),
            shell_text().replace("triangle 0 1 6 body", "triangle 0 6 1 body"),
            shell_text().replace("bridges 1\nbridge 12 0 0 1", "bridges 2\nbridge 12 0 0 1\nbridge 11 0 0 1"),
        ] {
            let mut loads = 0;
            assert!(shell_spec().compile_with_body(cx, gate, Some(Profile::parse(text.as_bytes()).unwrap()),
                |name| { loads += 1; Ok(pack(name, 1.0, None)) }).is_err());
            assert_eq!(loads, 0);
        }
        for key in ["shear_modulus_13", "transverse_shear_viscosity_23"] {
            assert!(shell_spec().compile_with_body(cx, gate, Some(shell_profile()),
                |name| Ok(pack(name, 1.0, Some(key)))).is_err(), "borrowed missing {key}");
        }
        let too_small = shell_text().replace("projection-work 67108864", "projection-work 1");
        assert!(shell_spec().compile_with_body(cx, gate,
            Some(Profile::parse(too_small.as_bytes()).unwrap()),
            |name| Ok(pack(name, 1.0, None))).is_err());
    });
}

#[test]
fn shell_rotational_chord_uses_actual_3d_nodes_and_common_force_direction() {
    let mut spec = shell_spec();
    spec.strings.push(spec.strings[0].clone());
    let text = shell_text().replace("bridges 1\nbridge 12 0 0 1",
        "bridges 2\nbridge 12 0 0 1\nbridge 11 0 0 1");
    let profile = Profile::parse(text.as_bytes()).unwrap();
    let lever = profile.bridge_separation_m(&spec, 0, 1).unwrap();
    assert!((lever - 0.1).abs() < 1e-15);
    assert!(lever < fs_math::det::hypot(0.1, 0.0015));
    // Both points at equal height with a force direction parallel to their
    // baseline measure extension, not a work-conjugate angular coordinate.
    let parallel = text.replace("bridge 12 0 0 1", "bridge 13 1 0 0")
        .replace("bridge 11 0 0 1", "bridge 11 1 0 0");
    let parallel = Profile::parse(parallel.as_bytes()).unwrap();
    assert!(parallel.bridge_separation_m(&spec, 0, 1).is_err());
    let nearly_parallel = text.replace("bridge 12 0 0 1", "bridge 13 1 0 1e-16")
        .replace("bridge 11 0 0 1", "bridge 11 1 0 1e-16");
    let nearly_parallel = Profile::parse(nearly_parallel.as_bytes()).unwrap();
    assert!(nearly_parallel.bridge_separation_m(&spec, 0, 1).is_err());
    let different = Profile::parse(text.replace("bridge 11 0 0 1", "bridge 11 1 0 0").as_bytes()).unwrap();
    assert!(different.bridge_separation_m(&spec, 0, 1).is_err());
    let coincident = Profile::parse(text.replace("bridge 11", "bridge 12").as_bytes()).unwrap();
    assert!(coincident.bridge_separation_m(&spec, 0, 1).is_err());
}

#[path = "interface_tests.rs"]
mod interfaces;

#[path = "pressure_tests.rs"]
mod pressure_surfaces;

/// A flat `nx * ny` grid body in the layout of `arched.body`: `nodes` rows, `2 (nx-1)(ny-1) + extra`
/// triangle rows (the extras repeat triangle 0), one support and one bridge. Only counts and indices are
/// meaningful: this feeds the parser, not the geometry validator.
fn grid_body(nx: usize, ny: usize, declared_nodes: usize, extra_triangles: usize) -> String {
    let mut text: String = shell_text().lines().take(8).map(|l| format!("{l}\n")).collect();
    text.push_str(&format!("nodes {declared_nodes}\n"));
    for j in 0..ny {
        for i in 0..nx {
            text.push_str(&format!("node {} {} 0 0 0 1\n",
                i as f64 / (nx - 1) as f64, j as f64 / (ny - 1) as f64));
        }
    }
    let mut triangles = Vec::new();
    for j in 0..ny - 1 {
        for i in 0..nx - 1 {
            let n = j * nx + i;
            triangles.push([n, n + 1, n + nx + 1]);
            triangles.push([n, n + nx + 1, n + nx]);
        }
    }
    for _ in 0..extra_triangles { triangles.push(triangles[0]); }
    text.push_str(&format!("triangles {}\n", triangles.len()));
    for t in &triangles { text.push_str(&format!("triangle {} {} {} body\n", t[0], t[1], t[2])); }
    text.push_str(&format!("supports 1\nsupport 0 111111\nbridges 1\nbridge {} 0 0 1\n", nx * ny - 1));
    text
}

#[test]
fn shell_profile_admits_exactly_the_node_ceiling_and_refuses_one_more() {
    assert_eq!(MAX_NODES % 32, 0, "the grid below is 32 rows high");
    let nx = MAX_NODES / 32;
    assert!(Profile::parse(grid_body(nx, 32, MAX_NODES, 0).as_bytes()).is_ok());
    // One more node is refused when it is only declared and when it is really present.
    assert!(Profile::parse(grid_body(nx, 32, MAX_NODES + 1, 0).as_bytes()).is_err());
    assert_eq!((MAX_NODES + 1) % 3, 0, "MAX_NODES + 1 must factor as 3 * m for the real over-cap grid");
    assert!(Profile::parse(grid_body(3, (MAX_NODES + 1) / 3, MAX_NODES + 1, 0).as_bytes()).is_err());
}

#[test]
fn shell_profile_admits_exactly_the_triangle_ceiling_and_refuses_one_more() {
    let nx = MAX_NODES / 32;
    let base = 2 * (nx - 1) * 31;
    assert!(base <= MAX_TRIANGLES, "the node-cap grid must leave room for the triangle cap");
    assert!(Profile::parse(grid_body(nx, 32, MAX_NODES, MAX_TRIANGLES - base).as_bytes()).is_ok());
    assert!(Profile::parse(grid_body(nx, 32, MAX_NODES, MAX_TRIANGLES - base + 1).as_bytes()).is_err());
}

#[test]
fn shell_profile_admits_exactly_the_mode_ceiling_and_refuses_one_more() {
    let profile = shell_profile();
    let mut spec = shell_spec();
    spec.body.modes = MAX_MODES;
    assert!(profile.preflight(&spec).is_ok());
    spec.body.modes = MAX_MODES + 1;
    assert!(profile.preflight(&spec).is_err());
    spec.body.modes = 0;
    assert!(profile.preflight(&spec).is_err());
    // A flat plate still admits no more than its own 16.
    let plate = Profile::parse(include_str!("../../regional.body").as_bytes()).unwrap();
    spec.body.modes = input::MAX_PLATE_BODY_MODES;
    assert!(plate.preflight(&spec).is_ok());
    spec.body.modes = input::MAX_PLATE_BODY_MODES + 1;
    assert!(plate.preflight(&spec).is_err());
}

#[test]
fn shell_profile_admits_exactly_the_byte_ceiling_and_refuses_one_more() {
    let text = shell_text();
    // Pad one numeric token so the whole profile is exactly the ceiling, then one byte over.
    let padded = |extra: usize| text.replace("drilling 0.001", &format!("drilling 0.001{}", "0".repeat(extra)));
    let to_cap = super::super::MAX_SHELL_BYTES - text.len();
    assert_eq!(padded(to_cap).len(), super::super::MAX_SHELL_BYTES);
    assert!(Profile::parse(padded(to_cap).as_bytes()).is_ok());
    assert!(Profile::parse(padded(to_cap + 1).as_bytes()).is_err());
}

#[test]
fn shell_profile_admits_exactly_the_projection_work_ceiling_and_refuses_one_more() {
    let with = |work: usize| shell_text().replace("projection-work 67108864", &format!("projection-work {work}"));
    assert!(Profile::parse(with(MAX_SHELL_PROJECTION_WORK).as_bytes()).is_ok());
    assert!(Profile::parse(with(MAX_SHELL_PROJECTION_WORK + 1).as_bytes()).is_err());
    // A cap-sized component's own projection work fits the ceiling.
    let cost = fs_plate::shell::mitc::viscosity::projection_terms(MAX_NODES, MAX_TRIANGLES, MAX_MODES).unwrap();
    assert!(cost < MAX_SHELL_PROJECTION_WORK, "{cost} terms");
}

/// A flat clamped `nx * ny` plate in the layout of `arched.body`, with one bridge at its centre node.
/// Unlike `grid_body` this is a complete, geometrically valid shell: every boundary node is clamped.
fn clamped_plate_body(nx: usize, ny: usize, window: &str, work: usize) -> String {
    let mut text: String = shell_text().lines().take(8).map(|l| format!("{l}\n")).collect();
    text = text.replace("window 1 100000000", window)
        .replace("projection-work 67108864", &format!("projection-work {work}"));
    text.push_str(&format!("nodes {}\n", nx * ny));
    for j in 0..ny {
        for i in 0..nx {
            text.push_str(&format!("node {} {} 0 0 0 1\n",
                i as f64 / (nx - 1) as f64, j as f64 / (ny - 1) as f64));
        }
    }
    let mut triangles = Vec::new();
    for j in 0..ny - 1 {
        for i in 0..nx - 1 {
            let n = j * nx + i;
            triangles.push([n, n + 1, n + nx + 1]);
            triangles.push([n, n + nx + 1, n + nx]);
        }
    }
    text.push_str(&format!("triangles {}\n", triangles.len()));
    for t in &triangles { text.push_str(&format!("triangle {} {} {} body\n", t[0], t[1], t[2])); }
    let clamped: Vec<usize> = (0..nx * ny)
        .filter(|n| { let (i, j) = (n % nx, n / nx); i == 0 || j == 0 || i == nx - 1 || j == ny - 1 }).collect();
    text.push_str(&format!("supports {}\n", clamped.len()));
    for n in &clamped { text.push_str(&format!("support {n} 111111\n")); }
    text.push_str(&format!("bridges 1\nbridge {} 0 0 1\n", (ny / 2) * nx + nx / 2));
    text
}

/// Rendered samples of the end-to-end ceiling test. A COST choice, not a limit: nothing caps the sample count. Measured on ts1
/// (release, load average about 60, 98-coordinate prepared network, 2026-10-23, bead .358): construction takes 16 s and rendering
/// about 1.0 ms per sample; before the pre-summed loss (A3) and the one-integral construction they were 115 s and 19.5 ms. 96 is
/// the length the curved-body tests in this file already render, and the length the route's before/after pressure-bit checks
/// compare. Every dimension the ceilings bound (nodes, triangles, modes, gas memory, projection work) is at full value.
const CEILING_SAMPLES: u64 = 96;

/// The ceiling fixture: the clamped 15 x 11 plate with `modes` retained modes under the cylinder-flow performance (its 12 gas-memory
/// intervals) for `samples` samples, with the prepared-linear record that any network above 42 modal coordinates needs.
fn ceiling_inputs(modes: usize, samples: u64) -> (Profile, input::Spec) {
    let text = clamped_plate_body(15, 11, "window 1 400000000", 2_000_000_000);
    let profile = Profile::parse(text.as_bytes()).unwrap();
    // 70 body modes, the cavity, two string planes, twist and the gas memory (2 planes x 12 intervals) make 98: declare room for 128.
    let performance = include_str!("../../cylinder-flow.material-strings")
        .replace("48000 4096 0.01 128", &format!("48000 {samples} 0.01 128"))
        .replace("cylinder-gas-loading 0.001 30 12",
            "cylinder-gas-loading 0.001 30 12\nprepared-linear 128 1073741824 10000000000000");
    let mut spec = input::Spec::parse(performance.as_bytes()).unwrap();
    spec.body.modes = modes;
    (profile, spec)
}

#[test]
fn a_mesh_and_mode_count_beyond_every_former_ceiling_run_through_the_whole_pipeline() {
    // 15 x 11 = 165 nodes (former ceiling 128), 280 triangles (256), 70 retained modes (16), and a projection
    // allowance of 2e9 terms (268,435,456). The window holds well over 70 bending modes of the fixture plate.
    let (nx, ny) = (15, 11);
    let (profile, spec) = ceiling_inputs(70, CEILING_SAMPLES);
    scoped(|cx, gate| {
        let constructing = std::time::Instant::now();
        let mut r = spec.compile_with_body(cx, gate, Some(profile), |name| Ok(if name == "body.fsmatpk" {
            pack_with_shear(name, 1.0, 1.0, None) } else { pack(name, 1.0, None) })).unwrap();
        // Timing lines go to stderr (visible with --nocapture): this test is slow by the cost measured above.
        eprintln!("ceiling e2e: {}", cost_context(&r.spec));
        eprintln!("ceiling e2e: constructed in {:.1?}", constructing.elapsed());
        let source = r.performance.bank().source().bank();
        let image = source.material_shell_image().unwrap();
        assert_eq!(source.material_shell().unwrap().specimen.mesh().nodes.len(), nx * ny);
        assert_eq!(source.material_shell().unwrap().specimen.mesh().tris.len(), 2 * (nx - 1) * (ny - 1));
        assert_eq!(image.physical_modes().len(), 70);
        assert!(image.physical_modes().iter().all(|v| v.len() == 6 * nx * ny));
        assert_eq!(image.damping().mode_count(), 70);
        let rendering = std::time::Instant::now();
        let samples = pressure(&mut r, 7);
        eprintln!("ceiling e2e: {} samples in {:.1?}", samples.len(), rendering.elapsed());
        assert!(samples.iter().all(|p| p.is_finite()));
        assert!(samples.iter().any(|p| p.abs() > 1e-14), "the body radiated nothing");
        eprintln!("ceiling e2e: pressure bits fnv1a {:#018x}", pressure_bits(&samples));
    });
}

/// FNV-1a over the exact bits of every pressure sample: an identity of the rendered sound, not a tolerance. It is printed,
/// not pinned: an exact restructuring of the route is checked by comparing it before and after the change.
fn pressure_bits(samples: &[f64]) -> u64 {
    samples.iter().flat_map(|p| p.to_bits().to_le_bytes())
        .fold(0xcbf2_9ce4_8422_2325_u64, |hash, byte| (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3))
}

/// What a wall time was measured on: host, logical CPUs, the load average when the line is printed, the render's thread,
/// and the complete network the prepared route propagates (coordinates, connections).
fn cost_context(spec: &input::Spec) -> String {
    let host = std::fs::read_to_string("/etc/hostname").map_or_else(|_| "unknown".into(), |s| s.trim().to_string());
    let cpus = std::thread::available_parallelism().map_or(0, |n| n.get());
    let load = std::fs::read_to_string("/proc/loadavg")
        .map_or_else(|_| "unavailable".into(), |s| s.split_whitespace().take(3).collect::<Vec<_>>().join(" "));
    let (coordinates, links, _) = spec.network_dimensions(true).unwrap();
    format!("host {host}, {cpus} logical CPUs, load average {load}, render on the test thread, \
        {coordinates} coordinates, {links} connections")
}

/// The ceiling fixture rendered for half a second (24000 samples at 48 kHz) through the prepared exact route: the wall time
/// per audio second with what it was measured on, and the pressure bits of the whole half second.
#[test]
#[ignore = "half-second ceiling render: run with --ignored --nocapture"]
fn the_half_second_ceiling_render_reports_its_cost_and_pressure_bits() {
    let (profile, spec) = ceiling_inputs(70, 24_000);
    scoped(|cx, gate| {
        eprintln!("half-second ceiling: {}", cost_context(&spec));
        let constructing = std::time::Instant::now();
        let mut r = spec.compile_with_body(cx, gate, Some(profile), |name| Ok(if name == "body.fsmatpk" {
            pack_with_shear(name, 1.0, 1.0, None) } else { pack(name, 1.0, None) })).unwrap();
        let built = constructing.elapsed();
        let rendering = std::time::Instant::now();
        let samples = pressure(&mut r, 7);
        let rendered = rendering.elapsed();
        assert_eq!(samples.len(), 24_000);
        assert!(samples.iter().all(|p| p.is_finite()));
        let audio_s = samples.len() as f64 / 48_000.0;
        eprintln!("half-second ceiling: constructed in {built:.1?}, {} samples rendered in {rendered:.1?} \
            ({:.1} s render per audio second, {:.1} s with construction)", samples.len(),
            rendered.as_secs_f64() / audio_s, (built + rendered).as_secs_f64() / audio_s);
        eprintln!("half-second ceiling: {}", cost_context(&r.spec));
        eprintln!("half-second ceiling: pressure bits fnv1a {:#018x}", pressure_bits(&samples));
    });
}

/// Cost probe, not a correctness check (`-- --ignored --nocapture`): prints how the prepared route's construction time and time
/// per sample grow with the retained body modes of the same plate, so the growth is measured rather than assumed.
#[test]
#[ignore = "timing probe: run with --ignored --nocapture"]
fn prepared_body_cost_against_retained_modes() {
    for modes in [10, 20, 40, 70] {
        let (profile, spec) = ceiling_inputs(modes, 28);
        scoped(|cx, gate| {
            let constructing = std::time::Instant::now();
            let mut r = spec.compile_with_body(cx, gate, Some(profile), |name| Ok(if name == "body.fsmatpk" {
                pack_with_shear(name, 1.0, 1.0, None) } else { pack(name, 1.0, None) })).unwrap();
            let built = constructing.elapsed();
            let rendering = std::time::Instant::now();
            let samples = pressure(&mut r, 7);
            eprintln!("cost probe: {modes} body modes: constructed in {built:.1?}, {:.1?} per sample",
                rendering.elapsed() / samples.len() as u32);
        });
    }
}
