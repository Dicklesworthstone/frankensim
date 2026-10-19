//! G0/G1/G2/G5 declared-band bake qualification. Native execution is required.
use fs_bem::helmholtz::iterative::IterativeOptions;
use fs_bem::helmholtz::{MAX_DENSE_PANELS, Medium};
use fs_bem::panel3d::SpherePanels;
use fs_bem::radiation_bake::declared::{
    BakeDeclaration, BakeError, BakeMesh, BakeStatus, ModeContent, PanelLadder, PanelSource,
    ReferencePlane, bake_declared, bake_declared_routed,
};
use fs_blake3::{ContentHash, hash_domain};
use fs_math::c64::C64;

struct Sphere {
    radius: f64,
}
impl PanelSource for Sphere {
    fn geometry_hash(&self) -> ContentHash {
        hash_domain(
            "analytic-sphere-radius-metres",
            &self.radius.to_bits().to_le_bytes(),
        )
    }
    fn reference_plane(&self) -> ReferencePlane {
        ReferencePlane::WholeSurface
    }
    fn surface_area_m2(&self) -> f64 {
        4.0 * std::f64::consts::PI * self.radius * self.radius
    }
    fn panel_count(&self, resolution: usize) -> Result<usize, BakeError> {
        4usize
            .checked_pow(
                u32::try_from(resolution).map_err(|_| BakeError::Invalid("sphere resolution"))?,
            )
            .and_then(|n| n.checked_mul(20))
            .ok_or(BakeError::Invalid("sphere panel count overflow"))
    }
    fn mesh(&self, resolution: usize) -> Result<BakeMesh, BakeError> {
        let surface = SpherePanels::icosphere(self.radius, resolution as u32)
            .map_err(|_| BakeError::Invalid("sphere mesh refused"))?;
        let driven = vec![true; surface.areas().len()];
        Ok(BakeMesh { surface, driven })
    }
}

fn declare(source: &dyn PanelSource, kas: &[f64], tol: f64) -> BakeDeclaration {
    let medium = Medium::air();
    let frequencies: Vec<_> = kas
        .iter()
        .map(|ka| ka * medium.sound_speed / std::f64::consts::TAU)
        .collect();
    BakeDeclaration::new(
        source.geometry_hash(),
        (frequencies[0], frequencies[frequencies.len() - 1]),
        frequencies,
        tol,
        source.reference_plane(),
        ModeContent::Plane,
        2,
        medium,
    )
    .unwrap()
}

fn print_rows(bake: &fs_bem::radiation_bake::declared::DeclaredBake) {
    let host = std::env::var("HOSTNAME")
        .ok()
        .filter(|s| !s.trim().is_empty())
        .unwrap_or_else(|| {
            let output = std::process::Command::new("uname")
                .arg("-n")
                .output()
                .expect("native runner hostname");
            assert!(output.status.success());
            String::from_utf8(output.stdout)
                .expect("hostname UTF-8")
                .trim()
                .to_string()
        });
    assert!(!host.is_empty());
    let command = std::env::args().collect::<Vec<_>>().join(" ");
    for row in bake.receipt_lines(&host, &command) {
        println!("{row}");
    }
}

/// G1: real solver calls, no supplied/fabricated operator rows. This broad declared
/// software-test tolerance exercises completion; it does NOT qualify 1e-3 physics.
#[test]
fn bake_band_complete_ladder() {
    let source = Sphere { radius: 1.0 };
    let declaration = declare(&source, &[0.05, 0.1], 0.25);
    let ladder = PanelLadder::new(&source, vec![1, 2, 3]).unwrap();
    let bake = bake_declared(&declaration, &ladder, || false).unwrap();
    print_rows(&bake);
    assert_eq!(bake.status(), BakeStatus::Complete);
    assert_eq!(bake.rows().len(), declaration.frequencies_hz().len());
    for &frequency in declaration.frequencies_hz() {
        let rows: Vec<_> = bake
            .refinement()
            .iter()
            .filter(|r| r.frequency_hz == frequency)
            .collect();
        assert!(rows.len() >= 3);
        assert!(rows.windows(2).all(|r| r[0].panels < r[1].panels));
        for row in &rows[rows.len() - 2..] {
            assert!(row.z_rel_change.unwrap() <= declaration.operator_tolerance());
            assert!(row.sh_rel_change.unwrap() <= declaration.operator_tolerance());
        }
    }
}

/// G2 strict qualification. If the current dense owner caps before meeting this
/// test, that is a genuine blocked qualification, never a relaxed tolerance.
#[test]
fn bake_pulsating_sphere_band() {
    let source = Sphere { radius: 1.0 };
    let declaration = declare(&source, &[0.05, 0.5, 1.0], 1e-3);
    let ladder = PanelLadder::new(&source, vec![1, 2, 3, 4, 5]).unwrap();
    let bake = bake_declared(&declaration, &ladder, || false).unwrap();
    print_rows(&bake);
    assert_eq!(
        bake.status(),
        BakeStatus::Complete,
        "strict sphere qualification is blocked by the dense route"
    );
    for row in bake.rows() {
        let ka = row.k * source.radius;
        let oracle = C64::new(ka * ka, -ka).scale(
            declaration.medium().density * declaration.medium().sound_speed / (1.0 + ka * ka),
        );
        let error = (row.z_specific - oracle).abs() / oracle.abs();
        println!(
            "{{\"case\":\"bake_pulsating_sphere_band\",\"inputs\":\"analytic pulsating sphere a=1m; e^-iwt\",\"axis\":\"ka\",\"rung\":{ka:e},\"observed_re\":{:e},\"observed_im\":{:e},\"err\":{error:e},\"tol\":1e-3,\"verdict\":\"{}\"}}",
            row.z_specific.re,
            row.z_specific.im,
            if error <= 1e-3 { "pass" } else { "fail" }
        );
        assert!(error <= 1e-3, "ka={ka}, analytic relative error {error}");
    }
}

/// G2 (M2.5-W5 route): the same pulsating-sphere oracle on a Burton-Miller
/// band (mouth ka >= 0.5) whose ladder needs the 20480-panel rung, past the
/// dense cap. The iterative owner solves that rung. Plain-CBIE rows (mouth
/// ka < 0.5, the low end of bake_pulsating_sphere_band) have no iterative
/// quadrature and refuse there.
///
/// The declared operator tolerance is 5e-3, and the closed-form check below
/// stays at 1e-3. At 1e-3 the three-rung rule cannot complete on this ladder.
/// trj b181 measured successive changes at ka = 0.5 of 1.1e-2, 2.5e-3 and
/// 4.1e-4 (SH 1.2e-2, 3.1e-3, 7.5e-4) over 320, 1280, 5120 and 20480 panels.
/// Completing at 1e-3 needs icosphere-6 (81920 panels) at an FMM order whose
/// M2L work exceeds MAX_TRANSLATION_WORK. That remains the strict
/// qualification's blocker (bake_pulsating_sphere_band).
#[test]
fn bake_routes_burton_miller_rows_past_the_dense_cap() {
    let source = Sphere { radius: 1.0 };
    let declaration = declare(&source, &[0.5, 1.0], 5e-3);
    let ladder = PanelLadder::new(&source, vec![2, 3, 4, 5]).unwrap();
    let iterative = IterativeOptions {
        tol: 1e-9,
        fmm_order: 6,
        max_iters: 900,
        ..IterativeOptions::default()
    };
    let bake = bake_declared_routed(&declaration, &ladder, Some(&iterative), || false).unwrap();
    print_rows(&bake);
    assert_eq!(bake.status(), BakeStatus::Complete);
    assert!(
        bake.refinement()
            .iter()
            .any(|r| r.route == "iterative" && r.panels > MAX_DENSE_PANELS),
        "the qualification must exercise the iterative route"
    );
    for row in bake.rows() {
        let ka = row.k * source.radius;
        let oracle = C64::new(ka * ka, -ka).scale(
            declaration.medium().density * declaration.medium().sound_speed / (1.0 + ka * ka),
        );
        let error = (row.z_specific - oracle).abs() / oracle.abs();
        println!(
            "{{\"case\":\"bake_routes_burton_miller_rows_past_the_dense_cap\",\"inputs\":\"analytic pulsating sphere a=1m; e^-iwt\",\"axis\":\"ka\",\"rung\":{ka:e},\"err\":{error:e},\"tol\":1e-3}}"
        );
        assert!(error <= 1e-3, "ka={ka}, analytic relative error {error}");
    }
}

#[test]
fn bake_dense_cap_blocks_and_routes() {
    let source = Sphere { radius: 1.0 };
    let declaration = declare(&source, &[0.05, 0.1], 1e-3);
    // These are real source sphere sizes; preflight must block before allocating
    // the first 20480-panel surface or its >6 GiB dense matrix.
    let ladder = PanelLadder::new(&source, vec![5, 6, 7]).unwrap();
    let bake = bake_declared(&declaration, &ladder, || false).unwrap();
    print_rows(&bake);
    assert_eq!(
        bake.status(),
        BakeStatus::Blocked {
            at_hz: declaration.band_hz().0,
            needed_panels: 20480,
            route: "M2.5-W5"
        }
    );
    assert!(20480 > MAX_DENSE_PANELS);
    assert!(bake.rows().is_empty());
    assert!(bake.load_rows(declaration.hash()).is_err());
}

#[test]
fn bake_declaration_hash_binds_artifact() {
    let source = Sphere { radius: 1.0 };
    let original = declare(&source, &[0.05, 0.1], 1e-3);
    let changed = declare(&source, &[0.05, 0.2], 1e-3);
    let ladder = PanelLadder::new(&source, vec![5, 6, 7]).unwrap();
    let bake = bake_declared(&original, &ladder, || false).unwrap();
    assert_eq!(bake.declaration().hash(), original.hash());
    assert_ne!(original.hash(), changed.hash());
    assert_eq!(
        bake.require_declaration(changed.hash()),
        Err(BakeError::DeclarationMismatch)
    );
    assert_eq!(bake.require_declaration(original.hash()), Ok(()));
}

#[test]
fn higher_modes_are_an_explicit_dependency_and_cancellation_refuses() {
    let source = Sphere { radius: 1.0 };
    let plain = declare(&source, &[0.05, 0.1], 1e-3);
    let modal = BakeDeclaration::new(
        source.geometry_hash(),
        plain.band_hz(),
        plain.frequencies_hz().to_vec(),
        1e-3,
        ReferencePlane::WholeSurface,
        ModeContent::Azimuthal {
            modes: vec![(0, 0), (1, 0)],
        },
        2,
        plain.medium(),
    )
    .unwrap();
    let ladder = PanelLadder::new(&source, vec![1, 2, 3]).unwrap();
    assert!(matches!(
        bake_declared(&modal, &ladder, || false),
        Err(BakeError::Dependency { route: "M2.1-W3c" })
    ));
    assert!(matches!(
        bake_declared(&plain, &ladder, || true),
        Err(BakeError::Cancelled)
    ));
}

#[test]
fn bake_bitwise_repeat() {
    let source = Sphere { radius: 1.0 };
    let declaration = declare(&source, &[0.05, 0.1], 0.25);
    let ladder = PanelLadder::new(&source, vec![1, 2, 3]).unwrap();
    let a = bake_declared(&declaration, &ladder, || false).unwrap();
    let b = bake_declared(&declaration, &ladder, || false).unwrap();
    assert_eq!(a.status(), b.status());
    assert_eq!(a.rows().len(), 2);
    for (a, b) in a.rows().iter().zip(b.rows()) {
        assert_eq!(a.z_acoustic.re.to_bits(), b.z_acoustic.re.to_bits());
        assert_eq!(a.z_acoustic.im.to_bits(), b.z_acoustic.im.to_bits());
        for (x, y) in a
            .directivity
            .coefficients
            .iter()
            .zip(&b.directivity.coefficients)
        {
            assert_eq!(x.re.to_bits(), y.re.to_bits());
            assert_eq!(x.im.to_bits(), y.im.to_bits());
        }
    }
}

/// Authored closed-cylinder piston fixture, not an open finite-wall pipe. Each
/// profile/azimuth mesh is rebuilt from the same exact cylinder parameters.
struct CappedCylinder {
    length: f64,
    radius: f64,
}
impl PanelSource for CappedCylinder {
    fn geometry_hash(&self) -> ContentHash {
        let mut bytes = self.length.to_bits().to_le_bytes().to_vec();
        bytes.extend_from_slice(&self.radius.to_bits().to_le_bytes());
        hash_domain("authored-closed-cylinder-piston-SI", &bytes)
    }
    fn reference_plane(&self) -> ReferencePlane {
        ReferencePlane::DuctMouth {
            origin_m: [self.length, 0.0, 0.0],
            outward_normal: [1.0, 0.0, 0.0],
        }
    }
    fn surface_area_m2(&self) -> f64 {
        std::f64::consts::TAU * self.radius * (self.length + self.radius)
    }
    fn panel_count(&self, resolution: usize) -> Result<usize, BakeError> {
        if resolution > 12 {
            return Err(BakeError::Invalid("bounded cylinder resolution"));
        }
        let stations = 1usize << resolution;
        let around = 8 * stations;
        let rings = (self.radius / (self.length / stations as f64))
            .ceil()
            .max(1.0) as usize;
        Ok(2 * stations * around + 2 * around * (2 * rings - 1))
    }
    fn mesh(&self, resolution: usize) -> Result<BakeMesh, BakeError> {
        let stations = 1usize << resolution;
        let profile: Vec<_> = (0..=stations)
            .map(|i| (self.length * i as f64 / stations as f64, self.radius))
            .collect();
        let bell = fs_bem::radiation_bake::lathe_profile(&profile, 8 * stations)?;
        let driven = bell.driven_mask();
        let surface = SpherePanels::from_triangles(bell.triangles)
            .map_err(|_| BakeError::Invalid("cylinder panel geometry"))?;
        Ok(BakeMesh { surface, driven })
    }
}

#[test]
fn partial_band_artifact_refuses_outside() {
    let source = CappedCylinder {
        length: 0.2,
        radius: 0.01,
    };
    let frequencies = vec![10.0, 20.0, 1_000_000.0];
    let declaration = BakeDeclaration::new(
        source.geometry_hash(),
        (10.0, 1_000_000.0),
        frequencies,
        0.25,
        source.reference_plane(),
        ModeContent::Plane,
        2,
        Medium::air(),
    )
    .unwrap();
    let ladder = PanelLadder::new(&source, vec![1, 2, 3, 4, 5]).unwrap();
    let bake = bake_declared(&declaration, &ladder, || false).unwrap();
    print_rows(&bake);
    assert!(matches!(
        bake.status(),
        BakeStatus::Blocked {
            at_hz: 1_000_000.0,
            route: "M2.5-W5",
            ..
        }
    ));
    assert!(bake.load_rows(declaration.hash()).is_err());
    let partial = bake.partial_band(declaration.hash()).unwrap();
    let (lo, hi) = partial.band_refusal();
    assert!((lo - 10.0).abs() < 1e-12 && (hi - 20.0).abs() < 1e-12);
    let (omegas, z) = partial.load_rows().unwrap();
    let table =
        fs_duct::TabulatedLoad::try_new(omegas, z, &declaration.hash().to_hex(), 0.25).unwrap();
    let (a, b) = table.support();
    assert!(table.z_at(a).is_ok() && table.z_at(b).is_ok());
    assert!(matches!(
        table.z_at(a.next_down()),
        Err(fs_duct::DuctError::OutOfTable { .. })
    ));
    assert!(matches!(
        table.z_at(b.next_up()),
        Err(fs_duct::DuctError::OutOfTable { .. })
    ));
    assert!(matches!(
        table.z_at(std::f64::consts::TAU * 1_000_000.0),
        Err(fs_duct::DuctError::OutOfTable { .. })
    ));
}

/// Closed solid finite-wall tube: outer wall, inward bore wall, annular lip,
/// backed piston at the throat. The bore is open to the exterior at x=length.
struct FiniteWallPipe {
    radius: f64,
    length: f64,
    thickness: f64,
}
impl PanelSource for FiniteWallPipe {
    fn geometry_hash(&self) -> ContentHash {
        let mut bytes = Vec::new();
        for x in [self.radius, self.length, self.thickness] {
            bytes.extend_from_slice(&x.to_bits().to_le_bytes());
        }
        hash_domain("analytic-finite-wall-open-pipe-backed-piston-SI", &bytes)
    }
    fn reference_plane(&self) -> ReferencePlane {
        ReferencePlane::DuctMouth {
            origin_m: [0.0; 3],
            outward_normal: [1.0, 0.0, 0.0],
        }
    }
    fn surface_area_m2(&self) -> f64 {
        let a = self.radius;
        let b = a + self.thickness;
        let t = self.thickness;
        std::f64::consts::TAU * (b * (self.length + t) + a * self.length)
            + std::f64::consts::PI * (b * b + a * a + b * b - a * a)
    }
    fn panel_count(&self, resolution: usize) -> Result<usize, BakeError> {
        if resolution > 10 {
            return Err(BakeError::Invalid("finite-wall pipe resolution budget"));
        }
        let radial = 1usize << resolution;
        let around = 8 * radial;
        let axial = (self.length / self.radius).ceil() as usize * radial;
        let outer_radial =
            ((self.radius + self.thickness) / self.radius * radial as f64).ceil() as usize;
        Ok(4 * axial * around + 2 * around * (radial + outer_radial))
    }
    fn mesh(&self, resolution: usize) -> Result<BakeMesh, BakeError> {
        use fs_math::det;
        let radial = 1usize << resolution;
        let around = 8 * radial;
        let axial = (self.length / self.radius).ceil() as usize * radial;
        let a = self.radius;
        let b = a + self.thickness;
        let t = self.thickness;
        let outer_radial = (b / a * radial as f64).ceil() as usize;
        let ring = |x: f64, r: f64| {
            (0..around)
                .map(|j| {
                    let theta = std::f64::consts::TAU * j as f64 / around as f64;
                    [x, r * det::cos(theta), r * det::sin(theta)]
                })
                .collect::<Vec<_>>()
        };
        let mut triangles = Vec::new();
        for (outer, start, radius) in [(true, -t, b), (false, 0.0, a)] {
            for segment in 0..axial {
                let x0 = if segment == 0 {
                    start
                } else {
                    start + (self.length - start) * segment as f64 / axial as f64
                };
                let x1 = if segment + 1 == axial {
                    self.length
                } else {
                    start + (self.length - start) * (segment + 1) as f64 / axial as f64
                };
                let lo = ring(x0, radius);
                let hi = ring(x1, radius);
                for j in 0..around {
                    let next = (j + 1) % around;
                    let mut faces = [[lo[j], lo[next], hi[j]], [lo[next], hi[next], hi[j]]];
                    if !outer {
                        for face in &mut faces {
                            face.swap(1, 2);
                        }
                    }
                    triangles.extend(faces);
                }
            }
        }
        let inner = ring(self.length, a);
        let outer = ring(self.length, b);
        for j in 0..around {
            let next = (j + 1) % around;
            triangles.extend([
                [inner[j], outer[j], inner[next]],
                [inner[next], outer[j], outer[next]],
            ]);
        }
        // Rigid back first, driven throat piston last, both at wall-scale h.
        let mut mouth_start = 0;
        for (x, radius, rings, forward) in [(-t, b, outer_radial, false), (0.0, a, radial, true)] {
            if forward {
                mouth_start = triangles.len();
            }
            let first = ring(x, radius / rings as f64);
            for j in 0..around {
                let mut face = [[x, 0.0, 0.0], first[j], first[(j + 1) % around]];
                if !forward {
                    face.swap(1, 2);
                }
                triangles.push(face);
            }
            for shell in 1..rings {
                let lo = ring(x, radius * shell as f64 / rings as f64);
                let hi = ring(
                    x,
                    if shell + 1 == rings {
                        radius
                    } else {
                        radius * (shell + 1) as f64 / rings as f64
                    },
                );
                for j in 0..around {
                    let next = (j + 1) % around;
                    let mut faces = [[lo[j], hi[j], lo[next]], [lo[next], hi[j], hi[next]]];
                    if !forward {
                        for face in &mut faces {
                            face.swap(1, 2);
                        }
                    }
                    triangles.extend(faces);
                }
            }
        }
        let driven = (0..triangles.len()).map(|i| i >= mouth_start).collect();
        let surface = SpherePanels::from_triangles(triangles)
            .map_err(|_| BakeError::Invalid("finite-wall pipe surface"))?;
        Ok(BakeMesh { surface, driven })
    }
}

/// G2: finite wall thickness is an independent physical axis. The original LS
/// mouth impedance is transported to the ACTUAL driven plane by fs-duct's
/// lossless cylinder owner, never by an invented radiation end correction.
/// A dense-cap refusal is an honestly failed qualification until W5 is available.
#[test]
fn bake_thin_pipe_vs_levine_schwinger() {
    use fs_duct::levine_schwinger::{LsQuadrature, levine_schwinger_impedance};
    use fs_material::gas::{GasSpec, GasState};
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let medium = Medium {
        density: gas.density,
        sound_speed: gas.sound_speed,
    };
    let (radius, length) = (0.007, 0.224); // Explicit L/a=32; finite-length error remains part of the oracle check.
    let kas = [0.05, 0.5, 1.0, 2.0, 3.8];
    let frequencies: Vec<_> = kas
        .iter()
        .map(|ka| ka * gas.sound_speed / (std::f64::consts::TAU * radius))
        .collect();
    let duct = fs_duct::Duct {
        segments: vec![fs_duct::Segment::Cylinder { radius, length }],
    };
    let mut errors = Vec::new();
    for thickness_ratio in [0.08, 0.04, 0.02] {
        let source = FiniteWallPipe {
            radius,
            length,
            thickness: thickness_ratio * radius,
        };
        let declaration = BakeDeclaration::new(
            source.geometry_hash(),
            (frequencies[0], frequencies[4]),
            frequencies.clone(),
            1e-3,
            source.reference_plane(),
            ModeContent::Plane,
            8,
            medium,
        )
        .unwrap();
        let ladder = PanelLadder::new(&source, vec![0, 1, 2, 3]).unwrap();
        let bake = bake_declared(&declaration, &ladder, || false).unwrap();
        print_rows(&bake);
        assert_eq!(
            bake.status(),
            BakeStatus::Complete,
            "finite-wall LS qualification needs the W5 route; t/a={thickness_ratio}"
        );
        let mut worst = 0.0_f64;
        for row in bake.rows() {
            let mouth =
                levine_schwinger_impedance(&gas, radius, row.omega, LsQuadrature { nodes: 1024 })
                    .unwrap();
            let oracle = fs_duct::input_impedance_load(
                &duct,
                &gas,
                row.omega,
                fs_duct::LossModel::Lossless,
                mouth,
                None,
            )
            .unwrap()
            .impedance;
            let error = (row.z_acoustic - oracle).abs() / oracle.abs();
            worst = worst.max(error);
            println!(
                "{{\"case\":\"bake_thin_pipe_vs_levine_schwinger\",\"inputs\":\"Levine-Schwinger PhysRev.73.383 full WH integral; finite-wall pipe; lossless reference-plane transport\",\"axis\":\"t_over_a\",\"rung\":{thickness_ratio:e},\"ka\":{:e},\"observed_re\":{:e},\"observed_im\":{:e},\"err\":{error:e},\"tol\":1e-3,\"verdict\":\"{}\"}}",
                row.mouth_ka,
                row.z_acoustic.re,
                row.z_acoustic.im,
                if error <= 1e-3 { "pass" } else { "fail" }
            );
        }
        errors.push(worst);
    }
    assert!(
        errors[2] <= 1e-3 && errors[2] < errors[1] && errors[1] < errors[0],
        "thin-wall ladder {errors:?}"
    );
}
