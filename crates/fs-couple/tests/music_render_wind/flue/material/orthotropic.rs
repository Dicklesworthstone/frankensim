//! Principal-axis source changes through the original jet, gas and microphone.
use super::*;
use fs_couple::bernoulli_aperture::wall::material::WallMaterialModel;

fn material(omit: Option<&str>) -> NormalizedMaterialCardPack {
    let base = pack(1000.0, 1e6, None, 100_000.0);
    let mut claims = base.card().claims().clone();
    let template = claims.claims_for("young_modulus")[0].1.clone();
    for (name, dims, value) in [("young_modulus_1", Pressure::DIMS, 1e6),
        ("young_modulus_2", Pressure::DIMS, 4e6), ("poisson_ratio_12", Dims::NONE, 0.1),
        ("wall_viscosity_time", Dims([0,0,1,0,0,0]), 2e-6)] {
        if omit == Some(name) { continue; }
        let mut claim = template.clone();
        claim.key = PropertyKey::new(name, dims);
        claim.value = PropertyValue::Scalar { value, dims };
        claims.insert_claim(claim).unwrap();
    }
    let normalized = NormalizedPack::new("synthetic-directional-wall", "regression-v1",
        ContentHash([11;32]), "synthetic redistribution permitted", claims, vec![], vec![]).unwrap();
    NormalizedMaterialCardPack::new(base.card().id().clone(), normalized).unwrap()
}

fn input(axes: &str) -> String {
    source().replace("wall 1 isotropic", "wall 1 orthotropic")
        .replace("stress-free kelvin-voigt", &format!("strain-fixed {axes} proportional-kelvin-voigt wall_viscosity_time"))
}

fn normal_material(cross: f64, omit: Option<&str>) -> NormalizedMaterialCardPack {
    let base = material(None);
    let mut claims = base.card().claims().clone();
    let template = claims.claims_for("young_modulus_1")[0].1.clone();
    for (name, value) in [("wall_eta_hh", 400.0), ("wall_eta_ha", cross), ("wall_eta_aa", 90.0)] {
        if omit == Some(name) { continue; }
        let mut claim = template.clone();
        claim.key = PropertyKey::new(name, DynViscosity::DIMS);
        claim.value = PropertyValue::Scalar { value, dims:DynViscosity::DIMS };
        claims.insert_claim(claim).unwrap();
    }
    let normalized = NormalizedPack::new("synthetic-normal-viscous-wall", "regression-v1",
        base.source_artifact(), "synthetic redistribution permitted", claims, vec![], vec![]).unwrap();
    NormalizedMaterialCardPack::new(base.card().id().clone(), normalized).unwrap()
}

fn normal_input() -> String {
    input("1 2").replace("strain-fixed 1 2 proportional-kelvin-voigt wall_viscosity_time",
        "stress-free 1 2 normal-kelvin-voigt 1 0.8 wall_eta_hh wall_eta_ha wall_eta_aa")
}

#[test]
fn normal_viscosity_source_memory_reaches_jet_feedback_microphone_and_native_pcm() {
    use fs_couple::bernoulli_aperture::wall::material::WallMaterialLoss;
    let text = normal_input();
    let mut traces = Vec::new();
    let mut patches = Vec::new();
    let dir = directory();
    let source = dir.join("normal-wall.performance");
    std::fs::write(&source, &text).unwrap();
    for (cross, name) in [(-60.0, "negative"), (60.0, "positive")] {
        let material = normal_material(cross, None);
        let file = loaded(&text, &material, 37);
        let wall = &file.renderer().voice().network().unwrap().material_walls()[0];
        assert_eq!(wall.source().model, WallMaterialModel::OrthotropicNormal { axes:[0,1],
            properties:["wall_eta_hh", "wall_eta_ha", "wall_eta_aa"].map(str::to_owned),
            max_dt_over_tau:1.0, max_angular_step:0.8 });
        let WallMaterialLoss::Maxwell(spectrum) = &wall.source().loss else { panic!("lost exact axial memory") };
        assert_eq!(spectrum.material().terms.len(), 1);
        assert_eq!(wall.load().compliance_terms().len(), 1);
        assert!(wall.patch().wall.resistance > 0.0);
        patches.push(wall.patch());
        // The original network's actual memory stores energy while the full
        // pressure/flow work identity is checked independently by trace().
        let mut renderer = file.into_renderer();
        let mut maximum_memory = 0.0_f64;
        let mut pressure = Vec::new();
        while renderer.remaining_samples() != 0 {
            let mut sample = [0.0]; renderer.block(&mut sample).unwrap();
            pressure.push(sample[0]);
            let energy = renderer.voice().network().unwrap().wall_relaxation_energy_j(1).unwrap();
            assert!(energy.is_finite() && energy >= 0.0);
            maximum_memory = maximum_memory.max(energy);
        }
        assert!(maximum_memory > 0.0);
        let reference = trace(loaded(&text, &material, 37));
        assert_eq!(pressure, reference.pressure);
        let mut chunked = loaded(&text, &material, 512).into_renderer();
        let mut actual = vec![0.0; pressure.len()];
        for chunk in actual.chunks_mut(512) { chunked.block(chunk).unwrap(); }
        assert_eq!(actual, pressure);
        std::fs::write(dir.join("wall.fsmcdpk"), material.to_bytes()).unwrap();
        let output = dir.join(format!("normal-{name}.wav"));
        success(run(&source, &output, 37, false));
        assert_eq!(std::fs::read(&output).unwrap(), encode_pcm16_wav(&pressure, 48_000, 1.0).unwrap().0);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["\"loss\":\"maxwell\"", "\"maxwell_terms\":[{", "\"normal_viscosity\":{",
            "\"property\":\"wall_eta_hh\"", "\"property\":\"wall_eta_ha\"", "\"property\":\"wall_eta_aa\""] {
            assert!(metadata.contains(token), "{token}: {metadata}");
        }
        assert!(!metadata.contains("C_viscous=tau*C_elastic"));
        traces.push(reference);
    }
    assert_eq!(patches[0], patches[1], "only axial memory changes for this cross-viscosity sign swap");
    assert_ne!(traces[0].incoming, traces[1].incoming);
    assert_ne!(traces[0].pressure, traces[1].pressure);
}

#[test]
fn normal_viscosity_ingress_refuses_missing_nonpassive_or_unresolved_source_laws() {
    let text = normal_input();
    let gate = CancelGate::new_clock_free();
    for material in [normal_material(-60.0, Some("wall_eta_hh")),
        normal_material(-60.0, Some("wall_eta_ha")), normal_material(-60.0, Some("wall_eta_aa")),
        normal_material(500.0, None)] {
        assert!(FlueFilePerformance::from_bytes_with_sources(text.as_bytes(), 37, &gate,
            None, |_| Ok(material.clone())).is_err());
    }
    for changed in [text.replace("normal-kelvin-voigt 1 0.8", "normal-kelvin-voigt 0.01 0.8"),
        text.replace("normal-kelvin-voigt 1 0.8", "normal-kelvin-voigt 1 0.0001"),
        text.replace("wall_eta_ha", "young_modulus_1")] {
        assert!(FlueFilePerformance::from_bytes_with_sources(changed.as_bytes(), 37, &gate,
            None, |_| Ok(normal_material(-60.0, None))).is_err());
    }
}

#[test]
fn directional_wall_data_reaches_native_feedback_microphone_and_pcm() {
    let pack = material(None);
    let mut traces = Vec::new();
    let dir = directory();
    std::fs::write(dir.join("wall.fsmcdpk"), pack.to_bytes()).unwrap();
    for (axes, directions, e_hoop, name) in [("1 2", [0,1], 1e6, "hoop-1"), ("2 1", [1,0], 4e6, "hoop-2")] {
        let text = input(axes);
        let file = loaded(&text, &pack, 37);
        let network = file.renderer().voice().network().unwrap();
        let wall = &network.material_walls()[0];
        assert_eq!(wall.source().model, WallMaterialModel::Orthotropic { axes:directions,
            proportional_viscosity_time_property:Some("wall_viscosity_time".into()) });
        assert!(wall.source().material.property("young_modulus").is_none(),
            "unused isotropic claims must not become the directional law");
        assert_eq!(file.wall_records()[0].pack(), &pack);
        let expected = e_hoop / 0.96 * 0.0002 / 0.006f64.powi(2);
        assert!((wall.patch().wall.stiffness_per_area/expected-1.0).abs()<1e-14);
        assert!((wall.patch().wall.resistance/(expected*2e-6)-1.0).abs()<1e-14);
        let reference = trace(file); // checks original network work and live wall storage
        for block in [1,37,512] {
            let mut renderer = loaded(&text, &pack, block).into_renderer();
            let mut actual = vec![0.0; reference.pressure.len()];
            for chunk in actual.chunks_mut(block) { renderer.block(chunk).unwrap(); }
            assert_eq!(actual, reference.pressure);
        }
        let source = dir.join(format!("{name}.performance"));
        let output = dir.join(format!("{name}.wav"));
        std::fs::write(&source, &text).unwrap();
        success(run(&source, &output, 37, false));
        assert_eq!(std::fs::read(&output).unwrap(),
            encode_pcm16_wav(&reference.pressure, 48_000, 1.0).unwrap().0);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["orthotropic-principal-local-ring-membrane".to_string(),
            format!("\"material_axes\":{{\"hoop\":{},\"axial\":{}",directions[0]+1,directions[1]+1),
            "\"time_property\":\"wall_viscosity_time\"".to_string(), pack.content_hash().to_hex()] {
            assert!(metadata.contains(&token), "{token}: {metadata}");
        }
        traces.push(reference);
    }
    assert_ne!(traces[0].incoming, traces[1].incoming);
    assert_ne!(traces[0].pressure, traces[1].pressure);
}

#[test]
fn directional_ingress_refuses_missing_data_and_unsupported_constitutive_images() {
    let text = input("1 2");
    let gate = CancelGate::new_clock_free();
    for missing in ["young_modulus_1", "young_modulus_2", "poisson_ratio_12", "wall_viscosity_time"] {
        assert!(FlueFilePerformance::from_bytes_with_sources(text.as_bytes(), 37, &gate, None,
            |_| Ok(material(Some(missing)))).is_err());
    }
    for changed in [input("1 1"), input("0 2"), input("1 4"),
        text.replace("proportional-kelvin-voigt wall_viscosity_time", "kelvin-voigt"),
        text.replace("proportional-kelvin-voigt wall_viscosity_time", "maxwell 1 0.1 0.2 boundary_modulus boundary_time")] {
        assert!(FlueFilePerformance::from_bytes_with_sources(changed.as_bytes(), 37, &gate, None,
            |_| Ok(material(None))).is_err());
    }
}
