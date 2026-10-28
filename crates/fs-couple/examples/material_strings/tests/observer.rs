//! The original material shell, load, accepted velocities and common FIR reach playback.
use super::*;

const OBSERVED: &str = include_str!("../closed-shell-observed.material-strings");

pub(super) fn shell_pack(name: &str) -> NormalizedPack {
    let provenance = Provenance {
        source: "synthetic closed-shell observer integration".into(),
        license: "CC0-1.0".into(),
        artifact: Some(fs_blake3::ContentHash([7; 32])),
    };
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen: name.into(),
        method: "manufactured source coefficients".into(),
        artifact: fs_blake3::ContentHash([7; 32]),
        caveats: "not measured instrument data".into(),
        provenance: provenance.clone(),
    }).unwrap();
    let mut entries = vec![
        ("density", Density::DIMS, 700.0),
        ("young_modulus", Pressure::DIMS, 100_000.0),
        ("poisson_ratio", Dims::NONE, 0.3),
        ("shear_modulus_13", Pressure::DIMS, 100_000.0 / 2.6),
        ("shear_modulus_23", Pressure::DIMS, 100_000.0 / 2.6),
    ];
    for key in [
        "plane_stress_viscosity_11", "plane_stress_viscosity_22",
        "plane_stress_viscosity_66", "plane_stress_viscosity_12",
        "plane_stress_viscosity_16", "plane_stress_viscosity_26",
        "transverse_shear_viscosity_13", "transverse_shear_viscosity_23",
        "transverse_shear_viscosity_13_23",
    ] {
        entries.push((key, DynViscosity::DIMS, 0.0));
    }
    for (key, dims, value) in entries {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained()
                .with("T", 280.0, 320.0).with("omega", 0.0, 200_000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation.clone()],
            provenance: provenance.clone(),
        }).unwrap();
    }
    let pack = NormalizedPack::new(
        name, "synthetic-closed-shell-v1", fs_blake3::ContentHash([7; 32]),
        "synthetic redistribution permitted", claims, vec![], vec![],
    ).unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}

fn observed_source() -> input::Spec {
    input::Spec::parse(
        OBSERVED.replace("audio 8000 4096", "audio 8000 1024")
            .replace("shell_angular_refinement 0","shell_angular_refinement 1").as_bytes(),
    ).unwrap()
}

fn shell_ready(source: input::Spec, cx: &Cx<'_>, gate: &CancelGate) -> Ready {
    source.compile_with_body(
        cx, gate,
        Some(body::Profile::parse(include_bytes!("../closed-shell.body")).unwrap()),
        |name| Ok(if name == "closed-shell.fsmatpk" {
            shell_pack(name)
        } else {
            pack(name, 7800.0, 200_000.0, None)
        }),
    ).unwrap()
}

fn observed_ready(cx: &Cx<'_>, gate: &CancelGate) -> Ready {
    shell_ready(observed_source(), cx, gate)
}

#[test]
fn shell_listener_replays_the_original_body_velocity_once_per_accepted_sample() {
    scoped(|cx, gate| {
        let mut ready = observed_ready(cx, gate);
        let report = ready.performance.material_shell_observer().unwrap();
        let binding = ready.performance.bank().source().bank().body_radiation().unwrap();
        assert_eq!(report.basis_identity, binding.basis_identity());
        assert_eq!(&report.gas, binding.gas());
        assert_eq!(report.gas, ready.spec.gas);
        assert_eq!(report.body_modes, ready.spec.body.modes);
        assert_eq!(report.surface_panels, 80);
        let load=binding.material_shell_bem().unwrap();
        assert_eq!(load.reciprocal.angular_refinement,1);
        assert_eq!(load.reciprocal.unknowns,122);
        assert_eq!(report.quadrature.relative_tolerance,1e-7);
        assert!(report.sweep_terms<=report.options.maximum_sweep_terms);
        assert_eq!(report.transfer().filter_latency_samples, 256);
        assert!(report.transfer().geometric_delays_s.iter().all(|delay| *delay > 0.0));
        let mut reference = report.transfer().observer().unwrap();
        let body = ready.spec.strings.len();
        let mut expected = Vec::new();
        let mut actual = Vec::new();
        for _ in 0..ready.spec.samples {
            let mut output = [0.0];
            ready.performance.block(&mut output).unwrap();
            let velocity = ready.performance.bank().system().components()[body].states()
                .iter().map(|state| state.velocity_m_sqrt_kg_per_s).collect::<Vec<_>>();
            expected.push(reference.step(&velocity).unwrap()[0]);
            actual.push(output[0]);
        }
        assert_eq!(actual, expected, "listener must consume the original accepted body modes");
        assert!(actual.iter().any(|pressure| pressure.abs() > 1e-14));
        let mut partitioned = observed_ready(cx, gate);
        assert_eq!(pressure(&mut partitioned, 37), actual);
        let mut source = observed_source();
        let options = source.shell_observer.take().unwrap();
        let mut advanced = shell_ready(source, cx, gate);
        assert!(advanced.performance.material_shell_observer().is_none());
        advanced.performance.block(&mut [0.0]).unwrap();
        assert!(
            advanced.performance.with_material_shell_observer(options, gate).is_err(),
            "new zero observer history must not bind to an advanced unobserved performance"
        );
    });
}

#[test]
fn listener_source_refuses_unmatched_bands_cavities_and_unbounded_fft_before_material_io() {
    for (from, to) in [
        ("body-radiation bem 0.5 400 65 80 12000000000 36\n", ""),
        ("shell_angular_refinement 0","shell_angular_refinement 3"),
        ("shell_quadrature 1e-7 14 2000000","shell_quadrature 1e-7 17 2000000"),
        ("shell_quadrature 1e-7 14 2000000","shell_quadrature 1e-7 14 79"),
        ("0.01 120 400 256", "0.01 120 401 256"),
        ("body-radiation bem 0.5 400", "body-radiation bem 1 400"),
        ("cavity none", "cavity neck 0.01 0.01 0.01 0 1"),
        ("4 1024 2048 4096 8192", "4 1024 2048 4096 10000"),
        ("4 1024 2048 4096 8192", "3 1024 2048 4096"),
        ("4 1024 2048 4096 8192", "999999999 1024"),
        ("1e-10 1e-4", "1e-3 1e-4"),
        ("0.01 120 400 256", "0.01 120 400 0"),
        ("front 0.2 0.2 1.0", "front NaN 0.2 1.0"),
    ] {
        assert!(input::Spec::parse(OBSERVED.replace(from, to).as_bytes()).is_err(), "{to}");
    }
    scoped(|cx, gate| {
        let mut source = observed_source();
        source.body_radiation = None;
        let mut reads = 0;
        assert!(source.compile::<NormalizedPack>(cx, gate, |_| {
            reads += 1;
            Err("material reads are forbidden on unmatched observer input".into())
        }).is_err());
        assert_eq!(reads, 0);
    });
    let ordinary = input::Spec::parse(include_bytes!("../closed-shell-pluck.material-strings")).unwrap();
    assert!(ordinary.shell_observer.is_none());
}
