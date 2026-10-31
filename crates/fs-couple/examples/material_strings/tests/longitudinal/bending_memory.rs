//! Original performance data selects the shared hereditary string mechanics.
use super::*;
const MEMORY: &str = "bending-memory 16 1\nstring-memory 0 1 bending_delta bending_time\n";
fn base() -> String {
    text().replace("prepared-linear 64 67108864 200000000000\n", "")
}
fn selected() -> String { base().replace("cavity ",&format!("{MEMORY}cavity ")) }
fn spectral_pack(name: &str, delta: f64, tau: f64) -> NormalizedPack {
    layer_spectral_pack(name, &[("bending_delta",delta,"bending_time",tau)])
}
fn layer_spectral_pack(name: &str, branches: &[(&str, f64, &str, f64)]) -> NormalizedPack {
    let mut claims=pack(name,7800.0,200_000.0,None).claims().clone();
    let artifact=fs_blake3::ContentHash([37;32]);
    let provenance=Provenance {source:"manufactured bending memory, not measured wire data".into(),
        license:"CC0-1.0".into(),artifact:Some(artifact)};
    let observation=claims.register_observation(ObservationDataset {specimen:name.into(),
        method:"synthetic spectrum".into(),artifact,caveats:"not measured".into(),
        provenance:provenance.clone()}).unwrap();
    for (name,value,dims) in branches.iter().flat_map(|&(modulus,delta,time,tau)| [
        (modulus,delta,Pressure::DIMS), (time,tau,Dims([0,0,1,0,0,0]))
    ]) {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,200_000.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![observation],provenance:provenance.clone()}).unwrap();
    }
    let pack=NormalizedPack::new(name,"synthetic-string-memory-v1",artifact,
        "synthetic redistribution permitted",claims,vec![],vec![]).unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}
fn build(cx:&Cx<'_>,gate:&CancelGate,text:&str,delta:f64,tau:f64) -> Result<Ready> {
    input::Spec::parse(text.as_bytes())?.compile(cx,gate,|name|Ok(match name {
        "wire.fsmatpk"=>spectral_pack(name,delta,0.002),
        "sleeve.fsmatpk"=>spectral_pack(name,delta,tau),
        _=>pack(name,7800.0,200_000.0,None), // No spectrum demanded from the body.
    }))
}
#[test]
fn file_selected_memory_keeps_gas_axial_sources_and_callback_independent_pressure() {
    scoped(|cx,gate| {
        let text=selected();
        let mut a=build(cx,gate,&text,40e9,0.004).unwrap();
        let mut b=build(cx,gate,&text,40e9,0.012).unwrap();
        let source=a.performance.bank().source().bank();
        assert!(source.cylinder_gas_loading().is_some());
        assert!(source.longitudinal_waves().is_some());
        assert!(source.material_cavity().is_some());
        let memory=source.string_bending_relaxation().unwrap();
        assert_eq!(memory.channels().len(),4);
        assert_eq!(memory.spectra()[0][0].terms[0].1,0.002);
        assert_eq!(memory.spectra()[0][1].terms[0].1,0.004);
        assert_eq!(a.spec.bytes,b.spec.bytes);
        assert_eq!(a.catalog["wire.fsmatpk"].pack.content_hash(),b.catalog["wire.fsmatpk"].pack.content_hash());
        assert_ne!(a.catalog["sleeve.fsmatpk"].pack.content_hash(),b.catalog["sleeve.fsmatpk"].pack.content_hash());
        let x=pressure(&mut a,37); let y=pressure(&mut b,128);
        assert!(x.iter().any(|p|p.abs()>1e-15));
        assert!(x.iter().zip(&y).any(|(x,y)|(x-y).abs()>1e-18));
        assert_eq!(y,pressure(&mut build(cx,gate,&text,40e9,0.012).unwrap(),1));
        let frame=b.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
        let (wav,summary)=build(cx,gate,&text,40e9,0.004).unwrap()
            .write(Cursor::new(Vec::new()),gate).unwrap();
        assert_eq!(wav.into_inner(),fs_couple::pcm_wav::encode_pcm16_wav(&x,48000,0.01).unwrap().0);
        assert_eq!(summary.samples,129);
    });
}
#[test]
fn file_zero_spectrum_preserves_the_original_elastic_performance() {
    scoped(|cx,gate| {
        let mut original=build(cx,gate,&base(),0.0,0.004).unwrap();
        let mut zero=build(cx,gate,&selected(),0.0,0.004).unwrap();
        assert!(original.performance.bank().source().bank().string_bending_relaxation().is_none());
        assert!(zero.performance.bank().source().bank().string_bending_relaxation().is_some());
        assert_eq!(pressure(&mut original,37),pressure(&mut zero,128));
    });
}
#[test]
fn malformed_memory_and_incompatible_backend_refuse_without_source_io() {
    let text=selected();
    for record in ["bending-memory 0 1\nstring-memory 0 1 bending_delta bending_time\n",
        "bending-memory 16 0\n", "bending-memory 16 1\nstring-memory 1 1 bending_delta bending_time\n",
        "bending-memory 16 1\nstring-memory 0 0\n",
        "bending-memory 16 1\nstring-memory 0 1 same same\n",
        "bending-memory 16 1\nstring-memory 0 2 bending_delta bending_time bending_delta bending_time\n"] {
        assert!(input::Spec::parse(text.replace(MEMORY,record).as_bytes()).is_err());
    }
    assert!(input::Spec::parse(text.replace(RECORD,
        &format!("prepared-linear 64 67108864 200000000000\n{RECORD}")).as_bytes()).is_ok());
    scoped(|cx,gate| {
        let missing=text.replace("bending_delta bending_time","bending_delta absent_time");
        assert!(build(cx,gate,&missing,40e9,0.004).is_err());
        assert!(build(cx,gate,&text,-1.0,0.004).is_err());
        assert!(build(cx,gate,&text,40e9,0.0).is_err());
        // Programmatic input still obeys the complete prepared mode allowance.
        let mut spec=input::Spec::parse(text.as_bytes()).unwrap();
        spec.prepared_linear=Some(PreparedModalBudget {max_modes:1,..MEMORY_PREPARATION});
        let mut reads=0;
        assert!(spec.compile(cx,gate,|name| {reads+=1;Ok(pack(name,7800.0,200_000.0,None))}).is_err());
        assert_eq!(reads,0);
    });
}

const LAYER_MEMORY: &str = "bending-memory 16 1\nlayered-string-memory 0 2\nlayer-memory young_modulus 0\nlayer-memory young_modulus 2 coating_delta_1 coating_tau_1 coating_delta_2 coating_tau_2\n";
fn layered_build(cx: &Cx<'_>, gate: &CancelGate, text: &str, tau: f64) -> Result<Ready> {
    input::Spec::parse(text.as_bytes())?.compile(cx,gate,|name| Ok(match name {
        "sleeve.fsmatpk" => layer_spectral_pack(name,&[
            ("coating_delta_1",12e9,"coating_tau_1",0.002),
            ("coating_delta_2",28e9,"coating_tau_2",tau),
        ]),
        _ => pack(name,7800.0,200_000.0,None),
    }))
}
#[test]
fn data_selected_independent_layer_spectra_reach_pressure_without_dummy_core_data() {
    scoped(|cx,gate| {
        let text=base().replace("cavity ",&format!("{LAYER_MEMORY}cavity "));
        let mut a=layered_build(cx,gate,&text,0.004).unwrap();
        let mut b=layered_build(cx,gate,&text,0.012).unwrap();
        let memory=a.performance.bank().source().bank().string_bending_relaxation().unwrap();
        assert!(memory.spectra()[0][0].terms.is_empty());
        assert_eq!(memory.spectra()[0][1].terms.len(),2);
        assert_eq!(memory.channels().len(),4);
        assert_eq!(a.spec.bytes,b.spec.bytes);
        assert_eq!(a.catalog["wire.fsmatpk"].pack.content_hash(),b.catalog["wire.fsmatpk"].pack.content_hash());
        assert_ne!(a.catalog["sleeve.fsmatpk"].pack.content_hash(),b.catalog["sleeve.fsmatpk"].pack.content_hash());
        let x=pressure(&mut a,37); let y=pressure(&mut b,128);
        assert!(x.iter().any(|p|p.abs()>1e-15));
        assert!(x.iter().zip(&y).any(|(x,y)|(x-y).abs()>1e-18));
        assert_eq!(y,pressure(&mut layered_build(cx,gate,&text,0.012).unwrap(),1));
        let frame=b.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
        let (wav,summary)=layered_build(cx,gate,&text,0.004).unwrap()
            .write(Cursor::new(Vec::new()),gate).unwrap();
        assert_eq!(wav.into_inner(),fs_couple::pcm_wav::encode_pcm16_wav(&x,48000,0.01).unwrap().0);
        assert_eq!(summary.samples,129);
    });
}
#[test]
fn explicit_uniform_layers_preserve_existing_file_arithmetic() {
    scoped(|cx,gate| {
        let uniform=selected();
        let explicit=uniform.replace(MEMORY,"bending-memory 16 1\nlayered-string-memory 0 2\nlayer-memory young_modulus 1 bending_delta bending_time\nlayer-memory young_modulus 1 bending_delta bending_time\n");
        let mut a=build(cx,gate,&uniform,40e9,0.004).unwrap();
        let mut b=build(cx,gate,&explicit,40e9,0.004).unwrap();
        assert_eq!(pressure(&mut a,37),pressure(&mut b,128));
    });
}
#[test]
fn layer_selector_mismatch_refuses_before_source_io_and_missing_data_does_not_fall_back() {
    let text=base().replace("cavity ",&format!("{LAYER_MEMORY}cavity "));
    for record in [
        "bending-memory 16 1\nlayered-string-memory 0 1\nlayer-memory young_modulus 0\n",
        "bending-memory 16 1\nlayered-string-memory 0 0\n",
        "bending-memory 16 1\nlayered-string-memory 0 2\nlayer-memory young_modulus 0\nlayer-memory young_modulus 1 coating_delta_1 young_modulus\n",
    ] {
        assert!(input::Spec::parse(text.replace(LAYER_MEMORY,record).as_bytes()).is_err());
    }
    scoped(|cx,gate| {
        let missing=text.replace("coating_tau_2","absent_time");
        assert!(layered_build(cx,gate,&missing,0.004).is_err());
        assert!(layered_build(cx,gate,&text,-0.001).is_err());
    });
}
