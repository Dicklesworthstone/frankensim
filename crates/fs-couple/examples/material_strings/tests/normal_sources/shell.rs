//! Synthetic source/body/contact integration; no measured fingerboard claim.
use super::*;

const RECEIVER: &str = "shell-receiver board 0 11 0 0 1 0.02 0 0";

fn stop_text() -> String {
    PROFILE
        .replace(
            "material-stop finger",
            &format!("{RECEIVER}\nmaterial-stop finger"),
        )
        .replace(
            "positive fixed -1e-9 receiver.fsmatpk",
            "positive shell:board -1e-9 @receiver",
        )
}

fn stops(text: &str, gate: &CancelGate, receiving: &str) -> stops::Profile {
    let mut profile = stops::Profile::parse(text.as_bytes()).unwrap();
    profile
        .load_sources(gate, |path| {
            Ok(match path {
                "normal.fsintpk" => interface("wire.fsmatpk", receiving).to_bytes(),
                "normal.fsmodelpk" => models().to_bytes(),
                _ => panic!("unexpected contact source {path}"),
            })
        })
        .unwrap();
    profile
}

fn shell_text() -> String {
    include_str!("../../arched.body").replace("orthotropic12", "isotropic")
}

fn shell_bulk(name: &str, young: f64) -> NormalizedPack {
    let is_body = name.contains("body");
    let mut claims = ClaimSet::new();
    for (key, value, dims) in [
        (
            "density",
            if is_body { 600.0 } else { 7800.0 },
            Density::DIMS,
        ),
        (
            "young_modulus",
            if is_body { young } else { 200e9 },
            Pressure::DIMS,
        ),
        ("poisson_ratio", 0.3, Dims::NONE),
        ("yield_stress", 1e10, Pressure::DIMS),
        (ETA, 10.0, DynViscosity::DIMS),
        (SHEAR, 80e9, Pressure::DIMS),
        (SHEAR_ETA, 1000.0, DynViscosity::DIMS),
        ("young_modulus_1", 8e9, Pressure::DIMS),
        ("young_modulus_2", 2e9, Pressure::DIMS),
        ("poisson_ratio_12", 0.2, Dims::NONE),
        ("shear_modulus_12", 0.6e9, Pressure::DIMS),
        ("shear_modulus_13", 0.5e9, Pressure::DIMS),
        ("shear_modulus_23", 0.25e9, Pressure::DIMS),
    ] {
        claims.insert_claim(claim(key, value, dims)).unwrap();
    }
    for (key, value) in fs_material::plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES
        .into_iter()
        .zip([3e5, 6e4, 4e4, 1e4, 4e3, -3e3])
    {
        claims
            .insert_claim(claim(key, value, DynViscosity::DIMS))
            .unwrap();
    }
    for (key, value) in fs_couple::shear_shell::TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES
        .into_iter()
        .zip([9000.0, 7000.0, 0.0])
    {
        claims
            .insert_claim(claim(key, value, DynViscosity::DIMS))
            .unwrap();
    }
    pack_claims(name, claims)
}

fn ready(cx: &Cx<'_>, gate: &CancelGate, young: f64) -> Ready {
    resting_spec()
        .compile_with_profiles(
            cx,
            gate,
            Some(body::Profile::parse(shell_text().as_bytes()).unwrap()),
            Some(stops(&stop_text(), gate, "body.fsmatpk")),
            |name| Ok(shell_bulk(name, young)),
        )
        .unwrap()
}

fn shell_pair_bulk(name:&str,model:PlateMaterialModel)->NormalizedPack {
    let base=shell_bulk(name,6e9);
    if !name.contains("body") { return base; }
    let mut claims=ClaimSet::new();
    for (_,original) in base.claims().claims_ordered() {
        if (original.key.name=="poisson_ratio" && model!=PlateMaterialModel::Isotropic)
            || (original.key.name=="young_modulus" && model==PlateMaterialModel::IsotropicBulkShear) {
            continue;
        }
        let property=match original.key.name.as_str() {
            "poisson_ratio"=>claim("poisson_ratio",0.25,Dims::NONE),
            "shear_modulus"=>claim("shear_modulus",2.4e9,Pressure::DIMS),
            _=>original.clone(),
        };
        claims.insert_claim(property).unwrap();
    }
    if model==PlateMaterialModel::IsotropicBulkShear {
        claims.insert_claim(claim("bulk_modulus",4e9,Pressure::DIMS)).unwrap();
    }
    pack_claims(name,claims)
}

#[test]
fn independent_shell_elastic_pairs_reach_the_same_source_bound_contact_and_body_reaction() {
    scoped(|cx,gate| {
        let mut stiffnesses=Vec::new();let mut forces=Vec::new();
        for (token,model) in [("isotropic",PlateMaterialModel::Isotropic),
            ("isotropic-young-shear",PlateMaterialModel::IsotropicYoungShear),
            ("isotropic-bulk-shear",PlateMaterialModel::IsotropicBulkShear)] {
            let body=include_str!("../../arched.body").replace("orthotropic12",token);
            let mut instrument=resting_spec().compile_with_profiles(cx,gate,
                Some(body::Profile::parse(body.as_bytes()).unwrap()),
                Some(stops(&stop_text(),gate,"body.fsmatpk")),
                |name|Ok(shell_pair_bulk(name,model))).unwrap();
            let bank=instrument.performance.bank();
            let original=&bank.source().bank().material_shell().unwrap()
                .specimen.regions()[0].input().region.material;
            let contact=bank.system().obstacle_source(0).unwrap();
            let solid=contact.model().normal_state().elastic().surface_b_state();
            assert_eq!(solid.resolved().properties().len(),original.properties().len()+1);
            assert!(original.properties().iter().all(|p|solid.resolved().properties().contains(p)));
            assert_eq!(solid.yield_stress_pa(),1e10);
            if model!=PlateMaterialModel::Isotropic {
                assert!(solid.resolved().property("poisson_ratio").is_none());
                assert_eq!(solid.elastic_state().derivation().unwrap().parameters(),
                    &model.isotropic_parameters().unwrap().unwrap());
            }
            if model==PlateMaterialModel::IsotropicBulkShear {
                assert!(solid.resolved().property("young_modulus").is_none());
            }
            stiffnesses.push(contact.stiffness_n_per_m_3_2());
            let row=bank.source().stops()[0].body.as_ref().unwrap().shapes.clone();
            let mut samples=vec![0.0;instrument.spec.samples as usize];
            instrument.performance.block(&mut samples[..1]).unwrap();
            let system=instrument.performance.bank().system();
            forces.push(system.last_obstacle_frames().unwrap()[0].normal_force_n);
            assert!(forces.last().unwrap()> &0.0);
            let velocity=row.iter().zip(system.components()[1].states())
                .map(|(shape,state)|shape*state.velocity_m_sqrt_kg_per_s).sum::<f64>();
            assert!(velocity>0.0,"the contact must react onto the original shell");
            let frame=system.last_frame().unwrap();
            assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
            for block in samples[1..].chunks_mut(37) {instrument.performance.block(block).unwrap();}
            assert!(samples.iter().any(|value|*value!=0.0));
        }
        for i in 1..3 {
            assert!((stiffnesses[i]/stiffnesses[0]-1.0).abs()<2e-13);
            assert!((forces[i]/forces[0]-1.0).abs()<1e-10);
        }
    });
}

#[test]
fn material_receiver_uses_original_shell_point_and_reacts_in_the_same_step() {
    scoped(|cx, gate| {
        let mut instrument = ready(cx, gate, 8e9);
        assert_eq!(
            instrument.catalog.len(),
            2,
            "no separate receiver bulk was loaded"
        );
        let bank = instrument.performance.bank();
        let body = bank.source().stops()[0].body.as_ref().unwrap();
        assert_eq!(body.component, 0);
        let modes = bank
            .source()
            .bank()
            .material_shell_image()
            .unwrap()
            .physical_modes();
        // n=+z, r=+0.02x gives r cross n=-0.02y, including shell rotation.
        let expected: Vec<_> = modes
            .iter()
            .map(|mode| mode[6 * 11 + 2] - 0.02 * mode[6 * 11 + 4])
            .collect();
        for (actual, expected) in body.shapes.iter().zip(expected) {
            assert!((actual - expected).abs() <= 1e-14 * expected.abs().max(1e-300));
        }
        let receiver = bank
            .system()
            .obstacle_source(0)
            .unwrap()
            .model()
            .normal_state()
            .elastic()
            .surface_b_state()
            .resolved();
        assert_eq!(
            receiver.card_identity(),
            instrument.catalog["body.fsmatpk"].card.content_hash()
        );
        let row = body.shapes.clone();
        let mut samples = vec![0.0; instrument.spec.samples as usize];
        instrument.performance.block(&mut samples[..1]).unwrap();
        let system = instrument.performance.bank().system();
        assert!(system.last_obstacle_frames().unwrap()[0].normal_force_n > 0.0);
        let body_velocity = row
            .iter()
            .zip(system.components()[1].states())
            .map(|(shape, state)| shape * state.velocity_m_sqrt_kg_per_s)
            .sum::<f64>();
        assert!(
            body_velocity > 0.0,
            "the same accepted contact must push the receiver along its positive normal"
        );
        let frame = system.last_frame().unwrap();
        assert_eq!(frame.external_work_j, 0.0);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        for block in samples[1..].chunks_mut(37) {
            instrument.performance.block(block).unwrap();
        }
        assert!(samples.iter().any(|value| *value != 0.0));
        let mut partitioned = ready(cx, gate, 8e9);
        assert_eq!(samples, pressure(&mut partitioned, 128));
        let mut replacement = ready(cx, gate, 4e9);
        assert_eq!(instrument.spec.bytes, replacement.spec.bytes);
        assert_eq!(
            instrument.stop_profile.as_ref().unwrap().bytes,
            replacement.stop_profile.as_ref().unwrap().bytes
        );
        assert_eq!(
            instrument.catalog["wire.fsmatpk"].card,
            replacement.catalog["wire.fsmatpk"].card
        );
        let old = instrument
            .performance
            .bank()
            .system()
            .obstacle_source(0)
            .unwrap();
        let new = replacement
            .performance
            .bank()
            .system()
            .obstacle_source(0)
            .unwrap();
        assert!(new.stiffness_n_per_m_3_2() < old.stiffness_n_per_m_3_2());
        assert_ne!(samples, pressure(&mut replacement, 37));
    });
}

#[test]
fn assembled_receiver_keeps_original_part_material_and_common_mode_basis() {
    scoped(|cx, gate| {
        let body = include_str!("../../assembled.body")
            .replace("orthotropic12", "isotropic")
            .replace("region right @body", "region right right-body.fsmatpk");
        let stop = stop_text().replace("board 0 11", "board 1 4");
        let build = |receiving| {
            resting_spec().compile_with_profiles(
                cx,
                gate,
                Some(body::Profile::parse(body.as_bytes()).unwrap()),
                Some(stops(&stop, gate, receiving)),
                |name| Ok(shell_bulk(name, 8e9)),
            )
        };
        let mut instrument = build("right-body.fsmatpk").unwrap();
        let source = instrument.performance.bank().source();
        let row = &source.stops()[0].body.as_ref().unwrap().shapes;
        let modes = &source
            .bank()
            .material_shell_assembly_image()
            .unwrap()
            .physical_modes()[1];
        assert_eq!(
            row.len(),
            instrument.performance.bank().system().components()[1]
                .modes()
                .len()
        );
        for (actual, mode) in row.iter().zip(modes) {
            let expected = mode[6 * 4 + 2] - 0.02 * mode[6 * 4 + 4];
            assert!((actual - expected).abs() <= 1e-14 * expected.abs().max(1e-300));
        }
        let contact = instrument
            .performance
            .bank()
            .system()
            .obstacle_source(0)
            .unwrap();
        assert_eq!(
            contact
                .model()
                .normal_state()
                .elastic()
                .surface_b_state()
                .resolved()
                .card_identity(),
            instrument.catalog["right-body.fsmatpk"].card.content_hash()
        );
        assert!(
            pressure(&mut instrument, 31)
                .iter()
                .any(|value| *value != 0.0)
        );
        assert!(
            build("body.fsmatpk").is_err(),
            "a tied part cannot borrow another part's bulk identity"
        );
    });
}

#[test]
fn geometry_material_seams_and_actual_contact_depth_refuse_instead_of_becoming_fixed() {
    scoped(|cx, gate| {
        let build = |body: Option<body::Profile>, text: &str, spec: input::Spec| {
            spec.compile_with_profiles(
                cx,
                gate,
                body,
                Some(stops(text, gate, "body.fsmatpk")),
                |name| Ok(shell_bulk(name, 8e9)),
            )
        };
        assert!(build(None, &stop_text(), resting_spec()).is_err());
        for invalid in [
            stop_text().replace("board 0 11", "board 1 11"),
            stop_text().replace("board 0 11", "board 0 9999"),
        ] {
            assert!(
                build(
                    Some(body::Profile::parse(shell_text().as_bytes()).unwrap()),
                    &invalid,
                    resting_spec()
                )
                .is_err()
            );
        }
        let mut thin = resting_spec();
        thin.body.thickness = 0.00002;
        let error = build(
            Some(body::Profile::parse(shell_text().as_bytes()).unwrap()),
            &stop_text(),
            thin,
        )
        .err()
        .unwrap();
        assert!(
            error
                .to_string()
                .contains("receiver's actual section thickness"),
            "{error}"
        );
        let mixed = shell_text().replace("regions 1", "regions 2")
            .replace("nodes 25", "region patch patch-body.fsmatpk isotropic 1 0.17 0.8333333333333334 1 0 0\nnodes 25")
            .replace("triangle 6 12 11 body", "triangle 6 12 11 patch");
        let error = build(
            Some(body::Profile::parse(mixed.as_bytes()).unwrap()),
            &stop_text(),
            resting_spec(),
        )
        .err()
        .unwrap();
        assert!(
            error
                .to_string()
                .contains("crosses different resolved material states"),
            "{error}"
        );
        let error = build(
            Some(body::Profile::parse(include_bytes!("../../arched.body")).unwrap()),
            &stop_text(),
            resting_spec(),
        )
        .err()
        .unwrap();
        assert!(
            error.to_string().contains("actual isotropic elastic state"),
            "{error}"
        );
    });
}

#[test]
fn shell_receiver_grammar_requires_explicit_geometry_and_actual_bulk_binding() {
    for invalid in [
        stop_text().replace("0 11 0 0 1", "0 11 0 0 2"),
        stop_text().replace("0.02 0 0\n", "NaN 0 0\n"),
        stop_text().replace("shell:board", "shell:missing"),
        stop_text().replace("@receiver", "receiver.fsmatpk"),
        PROFILE.replace("receiver.fsmatpk", "@receiver"),
        stop_text().replace(RECEIVER, &format!("{RECEIVER}\n{RECEIVER}")),
    ] {
        assert!(
            stops::Profile::parse(invalid.as_bytes()).is_err(),
            "accepted {invalid}"
        );
    }
}
