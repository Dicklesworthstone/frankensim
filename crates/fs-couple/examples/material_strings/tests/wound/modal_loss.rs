//! Original core loss must survive a sliding winding, shared body, and gas inertia.
//! All coefficients here are manufactured; these are constitutive integration
//! checks, not an acoustic qualification of a real bowed instrument.
use super::*;

const INTRINSIC: &str = "intrinsic_bending_loss";
const EXPANSION: &str = "linear_thermal_expansion_coefficient";
const CAPACITY: &str = "specific_heat_capacity";
const CONDUCTIVITY: &str = "thermal_conductivity";

#[derive(Clone, Copy)]
struct Sources {
    tangent: f64,
    loss_band: f64,
    wrap_density: f64,
    wrap_band: f64,
    omit: Option<&'static str>,
}
impl Default for Sources {
    fn default() -> Self {
        Self {
            tangent: 0.003,
            loss_band: 100_000.0,
            wrap_density: 8900.0,
            wrap_band: 200_000.0,
            omit: None,
        }
    }
}

fn selected(selection: &str) -> input::Spec {
    input::Spec::parse(
        String::from_utf8(wound_spec().bytes)
            .unwrap()
            .replace("cavity neck", &format!("{selection}\ncavity neck"))
            .as_bytes(),
    )
    .unwrap()
}

fn source_pack(name: &str, source: Sources) -> NormalizedPack {
    if name == "body.fsmatpk" {
        return pack(name, 7800.0, 200_000.0, None);
    }
    let wrap = name == "wrap.fsmatpk";
    let band = if wrap { source.wrap_band } else { 200_000.0 };
    let base = pack(
        name,
        if wrap { source.wrap_density } else { 7800.0 },
        band,
        Some(ETA),
    );
    let kelvin = fs_qty::QuantitySpec::semantic(fs_qty::semantic::SemanticType::new(
        fs_qty::semantic::QuantityKind::AbsoluteTemperature,
        fs_qty::semantic::ValueForm::Static,
    ));
    let validity = |maximum| {
        ValidityDomain::unconstrained()
            .with_quantity("T", kelvin, 280.0, 320.0)
            .with("omega", 0.0, maximum)
    };
    // Core and wrap declare the SAME typed physical environment. No adapter
    // may invent temperature semantics to admit a thermal law or join sources.
    let mut claims = ClaimSet::new();
    for id in base.claims().observation_ids() {
        assert_eq!(
            claims
                .register_observation(base.claims().observation(id).unwrap().clone())
                .unwrap(),
            id
        );
    }
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        claim.validity = validity(band);
        claims.insert_claim(claim).unwrap();
    }
    if !wrap {
        for (key, dims, value, maximum) in [
            (INTRINSIC, Dims::NONE, source.tangent, source.loss_band),
            (
                EXPANSION,
                fs_material::state_point::INVERSE_TEMPERATURE_DIMS,
                12e-6,
                200_000.0,
            ),
            (
                CAPACITY,
                fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS,
                450.0,
                150_000.0,
            ),
            (
                CONDUCTIVITY,
                fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS,
                50.0,
                200_000.0,
            ),
        ] {
            if source.omit == Some(key) {
                continue;
            }
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(key, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: validity(maximum),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: vec![],
                    provenance: Provenance {
                        source: "manufactured wound-core modal-loss integration".into(),
                        license: "CC0-1.0".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    let pack = NormalizedPack::new(
        base.pack_id(),
        "synthetic-wound-core-modal-loss-v1",
        base.source_artifact(),
        base.redistribution_terms(),
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}

fn build(cx: &Cx<'_>, gate: &CancelGate, spec: input::Spec, source: Sources) -> Result<Ready> {
    spec.compile(cx, gate, |name| Ok(source_pack(name, source)))
}

#[test]
fn core_loss_angle_uses_helix_mass_and_changes_shared_pressure_without_retuning() {
    scoped(|cx, gate| {
        let selection = format!("bending-loss loss-angle {INTRINSIC}");
        let prepare = |tangent| {
            build(
                cx,
                gate,
                selected(&selection),
                Sources {
                    tangent,
                    ..Sources::default()
                },
            )
            .unwrap()
        };
        let mut a = prepare(0.002);
        let mut b = prepare(0.02);
        let section = specimen(&a).section();
        let loss = section.core().modal_loss().unwrap();
        assert_eq!(loss.loss_angle().unwrap().loss_property(), INTRINSIC);
        assert!(loss.thermoelastic().is_none());
        assert_eq!(section.omega_band_rad_s(), (0.0, 100_000.0));
        assert!(loss.elastic().layers()[0].material.property(ETA).is_none());
        assert!(section.winding_material().property(INTRINSIC).is_none());
        assert!(specimen(&a).string().kelvin_voigt_bending.is_none());
        assert!(specimen(&a).string().relaxation_bending.is_none());
        assert_eq!(
            section.winding_material().identity(),
            specimen(&b).section().winding_material().identity()
        );
        let (left, right) = (
            &a.performance.bank().cards()[0],
            &b.performance.bank().cards()[0],
        );
        assert_eq!(left.viscous_bending_n_m2_s, 0.0);
        assert_eq!(right.viscous_bending_n_m2_s, 0.0);
        assert_eq!(left.tension_n, 20.0);
        assert_eq!(right.tension_n, left.tension_n);
        assert_eq!(left.linear_density_kg_m, section.line_density_kg_m());
        assert!(left.linear_density_kg_m > loss.elastic().line_density_kg_m());
        assert_eq!(left.linear_density_kg_m, right.linear_density_kg_m);
        assert_eq!(left.bending_stiffness_n_m2, right.bending_stiffness_n_m2);
        for mode in 1..=left.mode_count {
            let omega = left.mode_omega_rad_s(mode);
            assert_eq!(omega, right.mode_omega_rad_s(mode));
            let wave = mode as f64 * std::f64::consts::PI / left.length_m;
            // Independent per-cycle bending energy: the winding supplies no
            // EI but DOES move, so core-only modal mass would be incorrect.
            let expected = (0.02 - 0.002) * left.bending_stiffness_n_m2 * wave.powi(4)
                / (2.0 * left.linear_density_kg_m * omega * omega);
            assert!(
                ((right.zetas[mode - 1] - left.zetas[mode - 1]) / expected - 1.0).abs() < 1e-10
            );
        }
        assert_eq!(a.spec.bytes, b.spec.bytes);
        assert_eq!(
            a.catalog["body.fsmatpk"].pack.content_hash(),
            b.catalog["body.fsmatpk"].pack.content_hash()
        );
        let p = pressure(&mut a, 37);
        let q = pressure(&mut b, 128);
        assert!(p.iter().any(|p| p.abs() > 1e-12));
        assert!(p.iter().zip(&q).any(|(p, q)| (p - q).abs() > 1e-12));
        assert_eq!(p, pressure(&mut prepare(0.002), 1));
        for ready in [&a, &b] {
            let frame = ready.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    });
}

#[test]
fn core_thermal_diffusion_and_excess_use_original_core_radius_and_loaded_total_mass() {
    scoped(|cx, gate| {
        for wrap_density in [8900.0, 17_800.0] {
            let prepare = |selection| {
                let mut spec = selected(selection);
                spec.cylinder_gas_loading = Some(CylinderFlowResolution {
                    minimum_wavenumber_radius: 0.001,
                    maximum_wavenumber_radius: 30.0,
                    intervals: 2,
                });
                build(
                    cx,
                    gate,
                    spec,
                    Sources {
                        wrap_density,
                        ..Sources::default()
                    },
                )
                .unwrap()
            };
            let a = prepare("bending-loss thermoelastic");
            let b = prepare(&format!("bending-loss thermoelastic excess {INTRINSIC}"));
            let section = specimen(&b).section();
            let combined = section.core().modal_loss().unwrap();
            let thermal = combined.thermoelastic().unwrap();
            let loss = combined.loss_angle().unwrap();
            assert_eq!(thermal.elastic(), loss.elastic());
            assert_eq!(thermal.elastic().outer_radius_m(), 0.0002);
            assert!(section.outer_radius_m() > thermal.elastic().outer_radius_m());
            assert_eq!(thermal.law().t0, 293.15);
            let expected_time = 7800.0 * 450.0 * 0.0002_f64.powi(2) / 50.0;
            assert!((thermal.base_time_s() / expected_time - 1.0).abs() < 2e-14);
            assert!(section.winding_material().property(CONDUCTIVITY).is_none());
            assert_eq!(section.omega_band_rad_s(), (0.0, 100_000.0));
            assert_eq!(specimen(&a).section().omega_band_rad_s(), (0.0, 150_000.0));
            assert!(combined.hysteretic_stiffness_n_m2(120_000.0).is_err());
            let loading = b
                .performance
                .bank()
                .source()
                .bank()
                .cylinder_gas_loading()
                .unwrap();
            assert_eq!(
                loading.strings()[0].flow().radius_m(),
                section.outer_radius_m()
            );
            assert_eq!(
                loading.strings()[0].solid_card().linear_density_kg_m,
                section.line_density_kg_m()
            );
            let (plain, total) = (
                &a.performance.bank().cards()[0],
                &b.performance.bank().cards()[0],
            );
            assert!(total.linear_density_kg_m > section.line_density_kg_m());
            assert_eq!(plain.linear_density_kg_m, total.linear_density_kg_m);
            assert_eq!(total.viscous_bending_n_m2_s, 0.0);
            for mode in 1..=total.mode_count {
                let omega = total.mode_omega_rad_s(mode);
                assert_eq!(omega, plain.mode_omega_rad_s(mode));
                let wave = mode as f64 * std::f64::consts::PI / total.length_m;
                let expected = 0.5 * loss.hysteretic_stiffness_n_m2() * wave.powi(4)
                    / total.linear_density_kg_m
                    / omega
                    / omega;
                assert!(
                    ((total.zetas[mode - 1] - plain.zetas[mode - 1]) / expected - 1.0).abs()
                        < 1e-10
                );
                // The cylinder owner replaces sampled gas drag, then rebuilds
                // the SAME thermal law at the loaded frequency. A retained dry
                // thermal rate or omitted core law cannot satisfy this balance.
                let thermal_stiffness =
                    thermal.elastic().bending_stiffness_n_m2() * thermal.loss_tangent(omega);
                let thermal_zeta = 0.5 * thermal_stiffness * wave.powi(4)
                    / total.linear_density_kg_m
                    / omega
                    / omega;
                let gas_zeta = 0.5
                    * loading.strings()[0]
                        .flow()
                        .low_tail_resistance_per_length_n_s_m2()
                    / total.linear_density_kg_m
                    / omega;
                assert!((plain.zetas[mode - 1] / (thermal_zeta + gas_zeta) - 1.0).abs() < 1e-10);
            }
        }
    });
}

#[test]
fn missing_core_coefficients_original_wrap_band_and_duplicate_bending_owner_refuse() {
    scoped(|cx, gate| {
        let angle = format!("bending-loss loss-angle {INTRINSIC}");
        for (selection, missing) in [
            (angle.as_str(), INTRINSIC),
            ("bending-loss thermoelastic", EXPANSION),
            ("bending-loss thermoelastic", CAPACITY),
            ("bending-loss thermoelastic", CONDUCTIVITY),
        ] {
            assert!(
                build(
                    cx,
                    gate,
                    selected(selection),
                    Sources {
                        omit: Some(missing),
                        ..Sources::default()
                    }
                )
                .is_err()
            );
        }
        for source in [
            Sources {
                tangent: -0.001,
                ..Sources::default()
            },
            Sources {
                loss_band: 1000.0,
                ..Sources::default()
            },
            Sources {
                wrap_band: 1000.0,
                ..Sources::default()
            },
        ] {
            assert!(build(cx, gate, selected(&angle), source).is_err());
        }
        let text = String::from_utf8(selected(&angle).bytes).unwrap();
        let duplicate = text.replace(
            "cavity neck",
            concat!(
                "bending-memory 4 1\n",
                "string-memory 0 1 relaxing_bending_modulus bending_relaxation_time\n",
                "cavity neck"
            ),
        );
        assert!(input::Spec::parse(duplicate.as_bytes()).is_err());
    });
}
