//! Axial source spectra must reach the original physical pressure trajectory.
use super::*;

const MEMORY: &str = "longitudinal-memory 8 young_modulus 1 axial_delta axial_time";
fn base() -> String {
    text().replace("prepared-linear 64 67108864 200000000000\n", "")
}
fn selected() -> String { base().replace(RECORD, &format!("{RECORD}\n{MEMORY}")) }
fn source_pack(name: &str, delta: f64, tau: f64, high: f64) -> NormalizedPack {
    let base = pack(name, 7800.0, 200_000.0, None);
    if name != "wire.fsmatpk" && name != "sleeve.fsmatpk" { return base; }
    let mut claims = base.claims().clone();
    for (name, value, dims) in [("axial_delta", delta, Pressure::DIMS),
        ("axial_time", tau, Dims([0,0,1,0,0,0]))] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,high),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
            provenance: Provenance { source: "manufactured axial spectrum, not instrument calibration".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    NormalizedPack::new(base.pack_id(), "synthetic-axial-memory-v1", base.source_artifact(),
        base.redistribution_terms(), claims, vec![], vec![]).unwrap()
}
fn build(cx: &Cx<'_>, gate: &CancelGate, input: &str, delta: f64, sleeve_tau: f64, high: f64)
    -> Result<Ready> {
    input::Spec::parse(input.as_bytes())?.compile(cx, gate, |name| {
        Ok(source_pack(name, delta, if name == "sleeve.fsmatpk" { sleeve_tau } else { 0.001 }, high))
    })
}

#[test]
fn g1_axial_layer_memory_matches_the_independent_complex_modulus_continuum() {
    scoped(|cx, gate| {
        let ready = build(cx, gate, &selected(), 20e9, 0.004, 200_000.0).unwrap();
        let source = ready.performance.bank().source().bank();
        let axial = source.longitudinal_waves().unwrap();
        let memory = axial.material_relaxation().unwrap();
        assert_eq!(memory.channels().len(), 4);
        assert_eq!(memory.sections()[0].spectra()[0].terms, vec![(20e9,0.001)]);
        assert_eq!(memory.sections()[0].spectra()[1].terms, vec![(20e9,0.004)]);
        let law = &axial.strings()[0];
        let section = memory.sections()[0].elastic();
        // Independent area integral for a uniform axial strain-gradient mode:
        // K*(omega) = k^2/mu * sum_l A_l E_l*(omega). No bending I, gas mass,
        // fitted loss factor, or material-memory response routine is used.
        for omega in [10.0_f64,250.0,1000.0,20_000.0] {
            for (mode, &n) in law.model().axial_modes().iter().enumerate() {
                let k = n as f64 * std::f64::consts::PI / 0.5;
                let scale = k*k/law.solid_linear_density_kg_m();
                let mut real = section.axial_stiffness_n()*scale;
                let mut imaginary = 0.0;
                for (geometry, tau) in section.layer_integrals().iter().zip([0.001,0.004]) {
                    let x = omega*tau;
                    real += 20e9*geometry.area_m2()*scale*x*x/(1.0+x*x);
                    imaginary += 20e9*geometry.area_m2()*scale*x/(1.0+x*x);
                }
                let mut actual_real = law.model().axial_omegas()[mode].powi(2);
                let mut actual_imaginary = 0.0;
                for channel in memory.channels().iter().filter(|c|c.mode==mode) {
                    let x = omega*channel.relaxation_time_s;
                    let stiffness = channel.coupling_per_s.powi(2);
                    actual_real += stiffness*x*x/(1.0+x*x);
                    actual_imaginary += stiffness*x/(1.0+x*x);
                }
                assert!((actual_real/real-1.0).abs()<2e-13);
                assert!((actual_imaginary/imaginary-1.0).abs()<2e-13);
            }
        }
        assert!(source.cylinder_gas_loading().is_some());
        assert!(law.transverse_linear_density_kg_m() > law.solid_linear_density_kg_m());
        assert!(memory.channels().iter().all(|channel|
            channel.axial_component==law.axial_component()));
    });
}

#[test]
fn g3_axial_material_substitution_changes_shared_pressure_and_preserves_callback_replay() {
    scoped(|cx, gate| {
        let input = selected();
        let mut slow = build(cx,gate,&input,20e9,0.004,200_000.0).unwrap();
        let mut fast = build(cx,gate,&input,20e9,0.0002,200_000.0).unwrap();
        assert_eq!(slow.spec.bytes,fast.spec.bytes);
        assert_eq!(slow.catalog["wire.fsmatpk"].pack.content_hash(),fast.catalog["wire.fsmatpk"].pack.content_hash());
        assert_ne!(slow.catalog["sleeve.fsmatpk"].pack.content_hash(),fast.catalog["sleeve.fsmatpk"].pack.content_hash());
        assert_eq!(slow.performance.bank().system().mode_count(),fast.performance.bank().system().mode_count());
        let expected = pressure(&mut slow,37);
        let changed = pressure(&mut fast,128);
        assert!(expected.iter().any(|p|p.abs()>1e-15));
        assert!(expected.iter().zip(&changed).any(|(a,b)|(a-b).abs()>1e-20));
        assert!(slow.performance.bank().system().material_memory().iter().any(|z|z.abs()>1e-20));
        assert_eq!(expected,pressure(&mut build(cx,gate,&input,20e9,0.004,200_000.0).unwrap(),1));
        let frame=slow.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
        let (wav,summary)=build(cx,gate,&input,20e9,0.004,200_000.0).unwrap()
            .write(Cursor::new(Vec::new()),gate).unwrap();
        assert_eq!(wav.into_inner(),fs_couple::pcm_wav::encode_pcm16_wav(&expected,48000,0.01).unwrap().0);
        assert_eq!(summary.samples,129);
    });
}

#[test]
fn g0_axial_memory_requires_original_waves_valid_sources_and_complete_history_budget() {
    let input=selected();
    for declaration in ["longitudinal-memory 0 young_modulus 1 axial_delta axial_time",
        "longitudinal-memory 3 young_modulus 1 axial_delta axial_time",
        "longitudinal-memory 8 young_modulus 0",
        "longitudinal-memory 8 young_modulus 1 same same",
        "longitudinal-memory 8 young_modulus 2 axial_delta axial_time axial_delta axial_time"] {
        assert!(input::Spec::parse(input.replace(MEMORY,declaration).as_bytes()).is_err());
    }
    assert!(input::Spec::parse(input.replace(&format!("{RECORD}\n"),"").as_bytes()).is_err());
    assert!(input::Spec::parse(input.replace(RECORD,
        &format!("prepared-linear 64 67108864 200000000000\n{RECORD}")).as_bytes()).is_err());
    scoped(|cx,gate| {
        for (delta,tau,high) in [(-1.0,0.004,200_000.0),(20e9,0.0,200_000.0),
            (20e9,1e-8,200_000.0),(20e9,0.004,1000.0)] {
            assert!(build(cx,gate,&input,delta,tau,high).is_err());
        }
        assert!(build(cx,gate,&input.replace("axial_time","absent_time"),20e9,0.004,200_000.0).is_err());
        let mut no_memory=build(cx,gate,&base(),0.0,0.004,200_000.0).unwrap();
        let mut zero=build(cx,gate,&input,0.0,0.004,200_000.0).unwrap();
        assert_eq!(pressure(&mut no_memory,37),pressure(&mut zero,1));
        let mut programmatic=input::Spec::parse(input.as_bytes()).unwrap();
        programmatic.longitudinal_memory.as_mut().unwrap().max_states=3;
        let mut reads=0;
        assert!(programmatic.compile(cx,gate,|name| { reads+=1;
            Ok(source_pack(name,20e9,0.004,200_000.0)) }).is_err());
        assert_eq!(reads,0);
    });
}
