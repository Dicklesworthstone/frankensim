//! Synthetic regression specimens, not measurements of cane or an instrument.
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::{TetElasticMaterial, TetMaterialField};
use fs_solid::reduce::material_field;
use fs_solid::reduce::{OrificeSpec, ValveCardRequest};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 77,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

struct Mesh {
    nodes: Vec<[f64; 3]>,
    tets: Vec<[usize; 4]>,
    fixed: Vec<usize>,
    tip: Vec<usize>,
}

impl Mesh {
    fn new() -> Self {
        let (nx, ny, nz) = (2usize, 2usize, 2usize);
        let node = |i, j, k| (i * (ny + 1) + j) * (nz + 1) + k;
        let mut nodes = Vec::new();
        for i in 0..=nx {
            for j in 0..=ny {
                for k in 0..=nz {
                    nodes.push([
                        0.010 * i as f64 / nx as f64,
                        0.012 * j as f64 / ny as f64,
                        0.004 * k as f64 / nz as f64,
                    ]);
                }
            }
        }
        let mut tets = Vec::new();
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let v = [
                        node(i, j, k),
                        node(i + 1, j, k),
                        node(i, j + 1, k),
                        node(i + 1, j + 1, k),
                        node(i, j, k + 1),
                        node(i + 1, j, k + 1),
                        node(i, j + 1, k + 1),
                        node(i + 1, j + 1, k + 1),
                    ];
                    for t in [
                        [0, 1, 3, 7],
                        [0, 3, 2, 7],
                        [0, 2, 6, 7],
                        [0, 6, 4, 7],
                        [0, 4, 5, 7],
                        [0, 5, 1, 7],
                    ] {
                        tets.push(t.map(|index| v[index]));
                    }
                }
            }
        }
        let mut fixed = Vec::new();
        let mut tip = Vec::new();
        for (n, p) in nodes.iter().enumerate() {
            if p[0] == 0.0 {
                fixed.extend([3 * n, 3 * n + 1, 3 * n + 2]);
            }
            if p[2] == 0.0 && p[0] >= 0.005 {
                tip.push(n);
            }
        }
        Self {
            nodes,
            tets,
            fixed,
            tip,
        }
    }

    fn request<'a>(
        &'a self,
        material: TetMaterialField<'a>,
    ) -> ValveCardRequest<'a, TetMaterialField<'a>> {
        ValveCardRequest {
            element_order: fs_solid::reduce::ValveElementOrder::Tet4,
            nodes_m: &self.nodes,
            tetrahedra: &self.tets,
            fixed_dofs: &self.fixed,
            material,
            orifice: OrificeSpec {
                face_nodes: &self.tip,
                opening_axis: [0.0, 0.0, 1.0],
                width_axis: [0.0, 1.0, 0.0],
                opposing_plane_offset_m: -0.0004,
            },
            window_hz: (1.0, 50_000.0),
            retained_compliance_floor: 0.5,
            loss_factor: 0.1,
            source_id: "synthetic/material-field-regression/v1",
        }
    }
}

fn material(e: f64, rho: f64, identity: u8) -> TetElasticMaterial {
    TetElasticMaterial::try_new(rho, e, 0.30, [identity; 32]).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e}"
    );
}

use fs_couple::material_assembly::valve::compile_pressure_loaded_reed_assembly;
use fs_couple::material_assembly::{
    AcousticMaterialBindings, MaterialSource, StringMaterialBinding, StringPrestressBinding,
};
use fs_couple::string_specimen::{StringGeometryConstraint, StringPrestress};
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::{
    AcousticAssembly, AmbientGas, BeatingReed, CylinderSegment, Listener, PrestressedString,
    ViscothermalDuct, WaveguideEnd,
};
use fs_solid::reduce::{PressureLoadedValve, SingleModeValveAdmission};

fn admission() -> SingleModeValveAdmission {
    // Acceptance for this SYNTHETIC coarse slab, not measured cane accuracy.
    SingleModeValveAdmission {
        minimum_opening_compliance_fraction: 0.90,
        minimum_pressure_compliance_fraction: 0.90,
        maximum_static_opening_relative_error: 0.10,
    }
}
fn valve(cx: &Cx<'_>, e: f64) -> PressureLoadedValve {
    let mesh = Mesh::new();
    let mat = material(e, 1050.0, 1);
    let mut load = vec![[0.0; 3]; mesh.nodes.len()];
    for &n in &mesh.tip {
        load[n][2] = 1e-4 / mesh.tip.len() as f64;
    }
    material_field::reduce_pressure_loaded_valve(
        &mesh.request(TetMaterialField::Uniform(&mat)),
        &load,
        cx,
    )
    .unwrap()
}
fn template() -> AcousticAssembly {
    AcousticAssembly {
        ambient: AmbientGas {
            temperature_k: 293.15,
            pressure_pa: 101_325.0,
            relative_humidity: 0.0,
        },
        string: None,
        pluck: None,
        bow: None,
        blow: None,
        duct: Some(ViscothermalDuct {
            segments: vec![CylinderSegment::cylinder(0.007, 0.50)],
            tone_holes: vec![],
            termination: WaveguideEnd::UnflangedOpen,
            wall: None,
        }),
        // These intentionally unrelated mechanics must not leak into the
        // material-bound specimen. Only the two playing controls survive.
        reed: Some(BeatingReed {
            rest_opening_m: 0.02,
            width_m: 0.03,
            closing_pressure_pa: 1.0,
            blowing_pressure_pa: 100.0,
            attack_s: 0.008,
            mass_kg: 1.0,
            stiffness_n_m: 1.0,
            damping_ratio: 1.0,
        }),
        soundboard: None,
        body_modes: vec![],
        plate: None,
        cavity: None,
        obstacles: vec![],
        contact_texture: None,
        listener: Listener { distance_m: 1.0 },
        sample_rate_hz: 48_000,
        duration_s: 4096.0 / 48_000.0,
    }
}

#[test]
fn compiled_material_reed_reaches_the_physical_exterior_microphone() {
    use fs_couple::acoustic_realize::{RealizeOptions, duct_observer::BaffledDuctObserver};
    use fs_couple::air_path::GasAbsorptionModel;
    use fs_couple::pcm_wav::baffled::CircularOutletReceiver;
    with_cx(|cx| {
        let mut input = template();
        input.duct.as_mut().unwrap().termination = WaveguideEnd::FlangedOpen;
        input.duration_s = 1024.0 / f64::from(input.sample_rate_hz);
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.15, 101_325.0).unwrap();
        let exterior = GasState::try_new(&GasSpec::dry_air_ussa1976(), 283.15, 101_325.0).unwrap();
        let options = RealizeOptions {
            gas: Some(gas), observer_absorption: Some(GasAbsorptionModel::Lossless),
            baffled_duct: Some(BaffledDuctObserver {
                segment_gases: None, exterior_gas: Some(exterior), material_wall: None,
                receiver: CircularOutletReceiver {
                    position_m: [0.0, 0.0, 0.5], radial_rings: 4, angular_points: 16,
                    maximum_frequency_hz: 1000.0,
                },
            }), ..RealizeOptions::default()
        };
        let render = |young| {
            let compiled = compile_pressure_loaded_reed_assembly(&input, valve(cx, young),
                admission(), &AcousticMaterialBindings::default()).unwrap();
            let resolved = compiled.realize_with_options(&options).unwrap();
            assert_eq!(resolved.gas, exterior);
            assert_eq!(resolved.duct_observer.as_ref().unwrap().segment_gases, vec![gas]);
            assert!(resolved.pressure_pa.iter().all(|p| p.is_finite()));
            assert!(resolved.pressure_pa.iter().any(|p| *p != 0.0));
            resolved
        };
        // The same actual gas and microphone retain sensitivity to the original
        // solid-reduced valve's mass/stiffness, without a second synthesis loop.
        assert_ne!(render(1e6).pressure_pa, render(4e6).pressure_pa);
    });
}

#[test]
fn material_swap_replaces_all_derived_reed_coefficients_and_keeps_the_template() {
    with_cx(|cx| {
        let template = template();
        let snapshot = template.clone();
        let a = compile_pressure_loaded_reed_assembly(
            &template,
            valve(cx, 1e6),
            admission(),
            &AcousticMaterialBindings::default(),
        )
        .unwrap();
        let b = compile_pressure_loaded_reed_assembly(
            &template,
            valve(cx, 4e6),
            admission(),
            &AcousticMaterialBindings::default(),
        )
        .unwrap();
        assert_eq!(template, snapshot);
        let ar = a.compiled().assembly().reed.unwrap();
        let br = b.compiled().assembly().reed.unwrap();
        let source = a.valve().card();
        assert_eq!(ar.rest_opening_m, source.rest_gap_m);
        assert_eq!(ar.width_m, source.width_m);
        assert_eq!(ar.mass_kg, source.modes[0].mass_kg);
        assert_eq!(ar.stiffness_n_m, source.modes[0].stiffness_n_m);
        assert_eq!(ar.damping_ratio, 0.5 * source.loss_factor);
        assert_eq!(
            ar.blowing_pressure_pa,
            snapshot.reed.unwrap().blowing_pressure_pa
        );
        assert_eq!(ar.attack_s, snapshot.reed.unwrap().attack_s);
        close(
            ar.stiffness_n_m * ar.rest_opening_m / ar.closing_pressure_pa,
            a.valve().effective_areas_m2()[0],
            1e-12,
        );
        close(br.stiffness_n_m, 4.0 * ar.stiffness_n_m, 1e-6);
        close(br.closing_pressure_pa, 4.0 * ar.closing_pressure_pa, 1e-6);
        close(br.mass_kg, ar.mass_kg, 1e-6);
        assert_ne!(a.valve().identity(), b.valve().identity());
        assert_eq!(a.compiled().assembly().duct, template.duct);
    });
}

#[test]
fn checked_reed_composes_with_independent_source_bound_string() {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_material::state_point::MaterialPropertySelection;
    use fs_qty::{Density, Dims, Pressure};
    with_cx(|cx| {
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            ("density", Density::DIMS, 1000.0),
            ("young_modulus", Pressure::DIMS, 1e8),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    provenance: Provenance {
                        source: "synthetic integration fixture".into(),
                        license: "CC0-1.0".into(),
                        artifact: None,
                    },
                    observations: vec![],
                })
                .unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic-string".into(),
                phase: "solid".into(),
                process: "test".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let point = QueryPoint::new().with("T", 293.15).unwrap();
        let bindings = AcousticMaterialBindings {
            string: Some(StringMaterialBinding {
                source: MaterialSource {
                    card: &card,
                    point: &point,
                    selection: MaterialPropertySelection::SingleClaimOnly,
                },
                geometry: StringGeometryConstraint::FixedRadius(0.0005),
                prestress: StringPrestressBinding::Prescribed(StringPrestress::FixedTension(30.0)),
            }),
            plate: None,
        };
        let mut template = template();
        template.string = Some(PrestressedString {
            length_m: 0.3,
            tension_n: 1.0,
            lin_density_kg_m: 1.0,
            axial_stiffness_n: 1.0,
            width_m: 0.01,
            n_modes: 4,
            damping_ratio: 0.0,
            rayleigh: None,
            bending_stiffness_n_m2: 0.0,
            kelvin_voigt_bending: None,
            relaxation_bending: None,
            polarization_detune: 0.0,
            moving_end: false,
            tension_loss_factor: None,
        });
        let snapshot = template.clone();
        let result = compile_pressure_loaded_reed_assembly(
            &template,
            valve(cx, 1e6),
            admission(),
            &bindings,
        )
        .unwrap();
        let string = result.compiled().assembly().string.as_ref().unwrap();
        assert!(result.compiled().string().is_some());
        close(
            string.lin_density_kg_m,
            1000.0 * core::f64::consts::PI * 0.0005_f64.powi(2),
            1e-12,
        );
        assert_eq!(string.tension_n, 30.0);
        assert_ne!(result.compiled().assembly().reed, snapshot.reed);
        assert_eq!(template, snapshot);
    });
}

#[test]
fn same_compiled_reed_uses_the_explicit_gas_renderer_without_a_second_loop() {
    use fs_couple::acoustic_realize::{AcousticRealizeError, BoreRealization, RealizeOptions};
    use fs_couple::air_path::GasAbsorptionModel;
    with_cx(|cx| {
        let mut input = template();
        input.ambient.temperature_k = -1.0; // invalid, and ignored ONLY by explicit-gas realization
        // An unflanged 7 mm mouth at 48 kHz has ka = 3.08 at Nyquist, past the
        // unflanged model's ka <= 1 (MAX_RADIATION_KA), so realization refuses.
        // This test is about the renderer, not the mouth: the flanged mouth's
        // Rayleigh piston covers the whole band.
        input.duct.as_mut().unwrap().termination = WaveguideEnd::FlangedOpen;
        let result = compile_pressure_loaded_reed_assembly(
            &input,
            valve(cx, 1e6),
            admission(),
            &AcousticMaterialBindings::default(),
        )
        .unwrap();
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 80_000.0).unwrap();
        let actual = result.realize_with_gas(&gas).unwrap();
        let direct = result.compiled().realize_with_gas(&gas).unwrap();
        assert_eq!(actual, direct);
        assert_eq!(actual.gas, gas);
        assert!(!actual.pressure_pa.is_empty());
        assert!(actual.pressure_pa.iter().all(|x| x.is_finite()));
        assert!(actual.pressure_pa.iter().any(|x| *x != 0.0));
        let denser = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
        let changed = result.realize_with_gas(&denser).unwrap();
        assert_eq!(changed.gas, denser);
        assert_ne!(actual.pressure_pa, changed.pressure_pa);
        assert!(changed.pressure_pa.iter().all(|x| x.is_finite()));
        assert!(result.realize().is_err());
        let mut invalid = gas;
        invalid.density = 0.0;
        assert!(result.realize_with_gas(&invalid).is_err());
        let (dry_air, air) =
            fs_material::iso9613::TerrestrialAir::moist_air(293.15, 80_000.0, 0.0).unwrap();
        assert_eq!(dry_air, gas);
        let options = RealizeOptions {
            gas: Some(gas),
            observer_absorption: Some(GasAbsorptionModel::Iso9613Air {
                relative_humidity: 0.5,
                air,
            }),
            ..RealizeOptions::default()
        };
        let selected = result.realize_with_options(&options).unwrap();
        assert_eq!(
            selected,
            result.compiled().realize_with_options(&options).unwrap()
        );
        assert_eq!(selected.gas, gas);
        assert_ne!(selected.pressure_pa, actual.pressure_pa);
        // A material-derived valve must retain the caller's numerical request:
        // the underresolved ladder refuses rather than silently becoming the
        // default characteristic line when compiled solid sources are present.
        // One cell is refused earlier as an invalid description ("at least 2
        // cells"). Two cells pass that check, and their 181.8-cent dispersion
        // is far over the 1-cent ceiling.
        let ladder = result.realize_with_options(&RealizeOptions {
            bore: BoreRealization::lumped_ladder(2),
            ..options
        });
        assert!(
            matches!(ladder, Err(AcousticRealizeError::LumpedLadderBand { .. })),
            "{:?}",
            ladder.as_ref().map(|_| "rendered")
        );
    });
}

#[test]
fn missing_components_bad_controls_and_overstrict_reductions_do_not_publish() {
    with_cx(|cx| {
        let base = template();
        let v = valve(cx, 1e6);
        let bindings = AcousticMaterialBindings::default();
        let mut candidate = base.clone();
        candidate.reed = None;
        assert!(
            compile_pressure_loaded_reed_assembly(&candidate, v.clone(), admission(), &bindings)
                .is_err()
        );
        candidate = base.clone();
        candidate.duct = None;
        assert!(
            compile_pressure_loaded_reed_assembly(&candidate, v.clone(), admission(), &bindings)
                .is_err()
        );
        candidate = base.clone();
        candidate.reed.as_mut().unwrap().blowing_pressure_pa = -1.0;
        assert!(
            compile_pressure_loaded_reed_assembly(&candidate, v.clone(), admission(), &bindings)
                .is_err()
        );
        let strict = SingleModeValveAdmission {
            minimum_pressure_compliance_fraction: 0.999999,
            ..admission()
        };
        assert!(
            compile_pressure_loaded_reed_assembly(&base, v.clone(), strict, &bindings).is_err()
        );
        let original = base.clone();
        compile_pressure_loaded_reed_assembly(&base, v, admission(), &bindings).unwrap();
        assert_eq!(base, original);
    });
}
