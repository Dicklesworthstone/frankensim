//! G1/G3/G4: original film work, finite heat and pressure share one accepted clock.
//! Synthetic liquid curves test causality; the sphere transfer remains the
//! independent observation fixture, not a measured lubricated-body radiator.
use super::*;
use crate::render::plate::impact::material::heating::{
    CaloricReservoir, SqueezeFilmHeatRoute,
};
use crate::render::plate::impact::material::heating::liquid_film::LiquidFilmBinding;
use crate::render::plate::impact::squeeze::{FilmCell, FilmChannel, FilmLimits, GapPort};
use crate::render::plate::impact::squeeze::material::MaterialSqueezeFilm;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EnthalpyPhaseKnot, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

fn card(hot_viscosity: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, value, dims, quantity) in [
        ("density", 1000.0, Dims([-3,1,0,0,0,0]),
            QuantitySpec::dimensional(Dims([-3,1,0,0,0,0]))),
        ("specific_heat_capacity", 0.2, Dims([2,0,-2,-1,0,0]),
            QuantitySpec::dimensional(Dims([2,0,-2,-1,0,0]))),
        ("dynamic_viscosity", 0.01, Dims([-1,1,-1,0,0,0]),
            QuantitySpec::dimensional(Dims([-1,1,-1,0,0,0]))),
        ("thermal_conductivity", 0.6, Dims([1,1,-3,-1,0,0]),
            QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static))),
    ] {
        let (value, interpolation) = if name == "dynamic_viscosity" {
            (PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, value), (300.001, hot_viscosity)], dims },
                InterpolationPolicy::LinearInside)
        } else {
            (PropertyValue::Scalar { value, dims }, InterpolationPolicy::ConstantWithinValidity)
        };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(name, quantity), value,
            validity: ValidityDomain::unconstrained().with("T",300.0,300.001)
                .with("p",100000.0,100000.0),
            uncertainty: UncertaintyModel::Unstated, interpolation, observations: vec![],
            provenance: Provenance { source: "synthetic liquid for thermal/acoustic integration".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic liquid".into(),
        phase: "liquid".into(), process: "numerical fixture".into(), revision: 0 },
        claims, vec![]).unwrap()
}

fn law(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
    LiquidTemperatureLaw::try_new(card,
        &QueryPoint::new().with("T",300.0).unwrap().with("p",100000.0).unwrap(),
        "T", LiquidViscosityBasis::Dynamic, MaterialPropertySelection::SingleClaimOnly).unwrap()
}

fn chart(card: &MaterialCard, upper_h: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(card.content_hash(), SolidLiquidPhase::Liquid,
        vec![
            EnthalpyPhaseKnot { specific_enthalpy_j_kg: 0.0, temperature_k: 300.0,
                liquid_mass_fraction: 1.0, bulk_density_kg_m3: 1000.0 },
            EnthalpyPhaseKnot { specific_enthalpy_j_kg: upper_h, temperature_k: 300.0+upper_h/0.2,
                liquid_mass_fraction: 1.0, bulk_density_kg_m3: 1000.0 },
        ]).unwrap()
}

fn mechanical(source: &ResolvedLiquidState) -> super::super::super::PreparedImpactSystem {
    let (driver, _) = ImpactBody::free_mass(0.2,0.0,0.001).unwrap();
    let (receiver, _) = ImpactBody::free_mass(0.4,0.0,0.0).unwrap();
    let system = ImpactSystem::new(vec![driver,receiver],vec![],vec![],vec![],ImpactConfig {
        dt_s: 1.0/f64::from(RATE), max_steps: 1028, maximum_energy_j: 1.0,
        energy_absolute_tolerance_j: 1e-14, energy_relative_tolerance: 1e-7,
        maximum_generalized_force: 1.0,
    }).unwrap();
    let gap = GapPort { reference_m: 2e-4, closure: vec![1.0/0.2_f64.sqrt(),-1.0/0.4_f64.sqrt()] };
    let binding = MaterialSqueezeFilm::new(source,
        vec![FilmCell { area_m2: 1e-4, gap: gap.clone() }],
        vec![FilmChannel { from: 0, to: None, width_m: 0.01, length_m: 0.02, gap }], 2,
        FilmLimits { minimum_cell_gap_m: 1e-5, maximum_gap_m: 5e-4, maximum_pressure_pa: 10000.0 }).unwrap();
    binding.attach(system).unwrap().0.prepare_analytic().unwrap()
}

fn heated<'m>(source: &ResolvedLiquidState, chart: &'m EquilibriumEnthalpyPhaseCurve)
    -> HeatedImpactSystem<'m>
{
    let reservoir = CaloricReservoir::new("liquid".into(),chart,
        chart.state_at_specific_enthalpy(0.0).unwrap(),0.001,1e-6,0.0006).unwrap();
    HeatedImpactSystem::new_with_films(mechanical(source),vec![reservoir],vec![],
        vec![SqueezeFilmHeatRoute { film: 0, to_a: "liquid".into(), to_b: "liquid".into(),
            fraction_to_a: 1.0 }],1e-13).unwrap()
}

fn coupled<'m>(card: &'m MaterialCard, chart: &'m EquilibriumEnthalpyPhaseCurve)
    -> LiquidFilmImpactSystem<'m>
{
    let law = law(card);
    LiquidFilmImpactSystem::new(heated(law.reference_state(),chart),vec![LiquidFilmBinding {
        film: 0, reservoir: "liquid".into(), law, maximum_relative_viscosity_change: 0.1,
    }]).unwrap()
}

fn observer<M: ImpactSource>(source: M, radiation: &ImpactRadiation) -> ImpactPressureRenderer<'_, M> {
    ImpactPressureRenderer::new(source,radiation,vec![VelocityProjection {
        input_id: "surface_velocity_m_s".into(), weights: vec![0.0,1.0/0.4_f64.sqrt()],
    }],vec![0.0;2],listener(1.0),RATE,97).unwrap()
}

#[test]
fn finite_heat_observes_the_original_pressure_and_retains_complete_thermal_frames() {
    let card = card(0.01);
    let law = law(&card);
    let chart = chart(&card,0.0002);
    let radiation = artifact();
    let mut raw = observer(mechanical(law.reference_state()),&radiation);
    let mut thermal = observer(heated(law.reference_state(),&chart),&radiation);
    let initial_energy = thermal.mechanics().mechanics().stored_energy_j();
    let mut peak = 0.0_f64;
    for _ in 0..32 {
        let (mut a,mut b) = ([0.0;16],[0.0;16]);
        raw.block(&mut a).unwrap();
        thermal.block(&mut b).unwrap();
        assert_eq!(a,b,"heat reception cannot re-integrate or change this fixed-viscosity trajectory");
        assert_eq!(raw.mechanics().state(),thermal.mechanics().state());
        let frame = thermal.last_mechanical_frame().unwrap();
        assert_eq!(frame.mechanical.sample,thermal.samples_rendered());
        assert!((frame.mechanical.time_s-thermal.samples_rendered() as f64/f64::from(RATE)).abs()<1e-15);
        assert!(frame.combined_energy_residual_j.abs()<1e-12);
        peak = peak.max(b.iter().fold(0.0_f64,|a,b|a.max(b.abs())));
    }
    let store = &thermal.mechanics().reservoirs()[0];
    assert!(peak>1e-6 && store.received_heat_j()>1e-8);
    assert!(store.state().temperature_k()>300.0);
    assert!((thermal.mechanics().mechanics().stored_energy_j()+store.storage_change_j()
        -initial_energy).abs()<1e-12);
}

#[test]
fn sourced_viscosity_feedback_changes_the_same_motion_and_pressure() {
    let card = card(0.001);
    let law = law(&card);
    let chart = chart(&card,0.0002);
    let radiation = artifact();
    let mut live = observer(coupled(&card,&chart),&radiation);
    let mut frozen = observer(heated(law.reference_state(),&chart),&radiation);
    let mut difference = 0.0_f64;
    for _ in 0..512 {
        let (mut a,mut b) = ([0.0],[0.0]);
        live.block(&mut a).unwrap();
        frozen.block(&mut b).unwrap();
        difference = difference.max((a[0]-b[0]).abs());
        let source = live.mechanics();
        let temperature = source.system().reservoirs()[0].state().temperature_k();
        let expected = law.at_temperature(Temperature::new(temperature)).unwrap();
        assert_eq!(source.sources()[0],expected);
        assert_eq!(source.system().mechanics().squeeze_film(0).unwrap().viscosity_pa_s(),
            expected.dynamic_viscosity_pa_s());
        let frame = live.last_mechanical_frame().unwrap();
        assert_eq!(frame.thermal.impact.mechanical.sample,live.samples_rendered());
        assert!(frame.thermal.impact.combined_energy_residual_j.abs()<1e-12);
    }
    assert!(live.mechanics().sources()[0].dynamic_viscosity_pa_s()<0.009);
    assert!(difference>1e-6,"changing material viscosity must reach the acoustic observer");
}

#[test]
fn heated_pcm_partition_and_cancel_resume_preserve_all_accepted_owners() {
    use crate::pcm_wav::observation::DecimatedRenderer;
    use crate::pcm_wav::stream::{Pcm16WavStream,render_pressure_pcm16,ScheduledWavProgress};
    use std::io::Cursor;
    let card = card(0.001);
    let chart = chart(&card,0.0002);
    let radiation = artifact();
    let make = || DecimatedRenderer::new(observer(coupled(&card,&chart),&radiation),RATE,RATE/4,43).unwrap();
    let mut baseline = make();
    let mut expected = vec![0.0;257];
    for block in expected.chunks_mut(17) { baseline.block(block).unwrap(); }
    let (expected_wav,_) = crate::pcm_wav::encode_pcm16_wav(&expected,RATE/4,0.02).unwrap();
    for block in [7,43] {
        let mut actual = make();
        let mut stream = Pcm16WavStream::new(Cursor::new(Vec::<u8>::new()),RATE/4,0.02,43).unwrap();
        let mut scratch = vec![0.0;block];
        let gate = CancelGate::new_clock_free();
        render_pressure_pcm16(&mut actual,&mut stream,&gate,&mut scratch,31).unwrap();
        let motion = actual.source().mechanics().state().to_vec();
        let stores = actual.source().mechanics().system().reservoirs().to_vec();
        let sources = actual.source().mechanics().sources().to_vec();
        let stopped = CancelGate::new_clock_free(); stopped.request();
        assert_eq!(render_pressure_pcm16(&mut actual,&mut stream,&stopped,&mut scratch,226).unwrap(),
            ScheduledWavProgress::Cancelled { samples: 0 });
        assert_eq!(actual.source().mechanics().state(),motion);
        assert_eq!(actual.source().mechanics().system().reservoirs(),stores);
        assert_eq!(actual.source().mechanics().sources(),sources);
        assert_eq!(stream.samples_written(),31);
        render_pressure_pcm16(&mut actual,&mut stream,&gate,&mut scratch,226).unwrap();
        let (bytes,summary) = stream.finish().unwrap();
        assert_eq!(bytes.into_inner(),expected_wav);
        assert_eq!(summary.samples,257);
        assert_eq!(actual.source().mechanics().samples_rendered(),1028);
        assert_eq!(actual.source().mechanics().state(),baseline.source().mechanics().state());
        assert_eq!(actual.source().mechanics().system().reservoirs(),
            baseline.source().mechanics().system().reservoirs());
        assert_eq!(actual.source().mechanics().sources(),baseline.source().mechanics().sources());
    }
}

#[test]
fn caloric_refusal_cannot_publish_a_mechanical_or_acoustic_sample() {
    let card = card(0.001);
    let chart = chart(&card,1e-10);
    let radiation = artifact();
    let mut source = coupled(&card,&chart);
    let motion = source.state().to_vec();
    let stores = source.system().reservoirs().to_vec();
    let sources = source.sources().to_vec();
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(matches!(source.advance(&[0.0;2],&cancelled),Err(ImpactError::Cancelled)));
    assert_eq!(source.state(),motion);
    assert_eq!(source.system().reservoirs(),stores);
    assert_eq!(source.sources(),sources);
    let mut renderer = observer(source,&radiation);
    let mut output = [42.0];
    assert!(renderer.block(&mut output).is_err());
    assert_eq!(renderer.samples_rendered(),0);
    assert_eq!(renderer.mechanics().samples_rendered(),0);
    assert_eq!(renderer.mechanics().state(),motion);
    assert_eq!(renderer.mechanics().system().reservoirs(),stores);
    assert_eq!(renderer.mechanics().sources(),sources);
    assert!(renderer.last_mechanical_frame().is_none());
    assert_eq!(output,[42.0]);
    assert!(matches!(renderer.block(&mut output),Err(RenderError::Poisoned)));
}
