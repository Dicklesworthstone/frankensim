//! Native combined-physics regressions: the original import/flow/thermal path,
//! not a replacement solver. Keep the selected vertex for all physical twins.
use super::*;
use fs_cli::cards::{CardPackKind, CardPackSet, RawCardPack};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack,
    NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel};
use fs_qty::QtyAny;
use fs_project::{ConductionRadiation, RadiatingSurface, spec::dims};

fn fixture(nonlinear: bool, radiation: bool) -> Fixture {
    let mut f = Fixture::new();
    let source = fs_blake3::hash_bytes(b"synthetic combined-thermal fixture; k(T) and gray emissivity");
    let provenance = || Provenance { source: "synthetic combined-physics test, not experimental data".into(),
        license: "CC-BY-4.0; redistribution permitted".into(), artifact: Some(source) };
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen: "combined-thermal synthetic solid".into(), method: "declared k(T) and constant emissivity".into(),
        artifact: source, caveats: "numerical regression only, not a validated material".into(),
        provenance: provenance(),
    }).unwrap();
    let kd = fs_conduction::CONDUCTIVITY_DIMS;
    let conductivity = if nonlinear {
        PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: dims::TEMPERATURE,
            knots: vec![(250.0,6.0),(450.0,16.0)], dims: kd }
    } else { PropertyValue::Scalar { value: 10.0, dims: kd } };
    for (property, dims, value, interpolation) in [
        ("thermal-conductivity", kd, conductivity,
            if nonlinear { InterpolationPolicy::LinearInside } else { InterpolationPolicy::ConstantWithinValidity }),
        (fs_conduction::SURFACE_EMISSIVITY_PROPERTY, fs_qty::Dims::NONE,
            PropertyValue::Scalar { value: 0.8, dims: fs_qty::Dims::NONE }, InterpolationPolicy::ConstantWithinValidity),
    ] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(property,dims), value,
            validity: fs_evidence::ValidityDomain::unconstrained().with("T",250.0,450.0),
            uncertainty: UncertaintyModel::Unstated, interpolation,
            observations: vec![observation], provenance: provenance(),
        }).unwrap();
    }
    let pack = NormalizedMaterialCardPack::new(MaterialStateId {
        chemistry: "combined-thermal".into(), phase: "solid".into(), process: "synthetic".into(), revision: 0,
    }, NormalizedPack::new("combined-thermal", "synthetic-v1", source,
        "CC-BY-4.0; redistribution permitted", claims, Vec::new(), Vec::new()).unwrap()).unwrap();
    let bytes = pack.to_bytes();
    let cards = CardPackSet::admit(vec![RawCardPack { kind: CardPackKind::Material,
        source: "combined-thermal-fixture".into(), bytes: bytes.clone(), expect: None }]).unwrap();
    let card = &cards.materials()[0];
    // This is the private fixture's newly created input file, not the tracked
    // reference pack. Its new exact identity travels with both physics laws.
    std::fs::write(f.sources.join("aa6061.fsmcdpk"),bytes).unwrap();
    let binding = &mut f.project.materials.as_mut().unwrap()[0];
    binding.card = card.card().to_hex(); binding.state = card.identity().to_string();
    binding.claim = None; binding.temp_lo.value = 250.0; binding.temp_hi.value = 450.0;
    binding.source = "synthetic combined-thermal regression".into(); binding.conductivity_tolerance = None;
    f.project.power.as_mut().unwrap()[0].watts.value = 300.0;
    f.project.solver.as_mut().unwrap().tolerance_rel = 1e-8;
    if radiation {
        f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().radiation = Some(ConductionRadiation {
            surfaces: vec![RadiatingSurface { name: "gray-wall".into(), target: "air".into(),
                card: card.card().to_hex(), claim: None,
                query_temperature: QtyAny::new(300.0,dims::TEMPERATURE),
                reservoir_temperature: QtyAny::new(285.0,dims::TEMPERATURE) }],
            max_iterations: 128, temperature_tolerance: QtyAny::new(1e-11,dims::TEMPERATURE),
            heat_tolerance: QtyAny::new(1e-8,dims::POWER), relaxation: 0.5,
        });
    }
    f
}

#[test]
fn coupled_thermal_native_adjoint_matches_retained_vertex_control_resolves() {
    for (nonlinear,radiation) in [(true,false),(false,true),(true,true)] {
        let f = fixture(nonlinear,radiation);
        let (baseline,field) = f.solve(&f.project,0);
        let mut request = f.project.clone();
        request.outputs.get_or_insert_with(Vec::new).push(fs_project::spec::OutputRequest {
            name: "temperature-max-adjoint".into(), kind: "report".into(), region: None,
        });
        let (receipt,with_goal) = f.solve(&request,1);
        assert_eq!(field.get("temperature"),with_goal.get("temperature"));
        for key in ["energy","conjugate","radiation"] { assert_eq!(baseline.get(key),receipt.get(key),"{key}"); }
        let adjoint = receipt.get("nominal_adjoint").unwrap();
        assert_eq!(adjoint.str_field("mode"),Some(if radiation { "radiation-full-air-feedback" }
            else { "nonlinear-solid-full-air-feedback" }));
        assert_eq!(adjoint.str_field("authority"),Some("Estimated"));
        assert!(adjoint.f64_field("true_relative_residual").unwrap() < 1e-10);
        let vertex = adjoint.f64_field("selected_vertex").unwrap() as usize;
        let rows = adjoint.get("parameters").unwrap().as_array().unwrap();
        assert_eq!(rows.iter().filter(|r| r.str_field("target")==Some("air-inlet-temperature")).count(),1);
        let mut controls = vec![("power",0.25),("air-inlet-temperature",0.01)];
        if radiation { controls.push(("radiation-reservoir-temperature",0.01)); }
        for (i,(target,step)) in controls.into_iter().enumerate() {
            let g = rows.iter().find(|r| r.str_field("target")==Some(target)).unwrap().f64_field("derivative").unwrap();
            let mut values = Vec::new();
            for (side,sign) in [-1.0,1.0].into_iter().enumerate() {
                let mut shifted = f.project.clone();
                match target {
                    "power" => shifted.power.as_mut().unwrap()[0].watts.value += sign*step,
                    "air-inlet-temperature" => inlet(&mut shifted,297.0+sign*step),
                    _ => shifted.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                        .radiation.as_mut().unwrap().surfaces[0].reservoir_temperature.value += sign*step,
                }
                let (_,field) = f.solve(&shifted,2+2*i+side);
                values.push(field.get("temperature").unwrap().as_array().unwrap()[vertex].as_f64().unwrap());
            }
            let expected = (values[1]-values[0])/(2.0*step);
            assert!((g-expected).abs() < 2e-3*expected.abs().max(1e-4),
                "nonlinear={nonlinear} radiation={radiation} {target}: adjoint {g:e}, physical difference {expected:e}");
        }
        if radiation {
            let rad = receipt.get("radiation").unwrap();
            let air = receipt.get("conjugate").unwrap();
            assert!(rad.f64_field("radiative_out_w").unwrap().abs() > 0.01);
            assert!((air.f64_field("air_total_w").unwrap()-rad.f64_field("convective_out_w").unwrap()).abs() < 1e-5,
                "only convection, never radiation, may heat the air");
        }
        assert!(adjoint.get("unsupported").unwrap().as_array().unwrap().iter()
            .any(|r| r.str_field("target")==Some("fan-speed-ratio")),
            "a fixed-flow thermal derivative must not claim a hydraulic control");
    }
}
