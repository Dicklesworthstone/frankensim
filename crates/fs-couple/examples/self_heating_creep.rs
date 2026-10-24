//! A loaded uniform specimen heats itself; its temperature changes its creep.
//! Existing Maxwell, caloric and convection/radiation owners do all physics.
//! Synthetic demonstration data, NOT a calibrated polymer or metal specimen.
#[path = "self_heating_creep/source.rs"]
mod source;
use fs_blake3::DomainHasher;
use fs_conduction::lumped::{LumpedThermalEnvironment, LumpedThermalTransport};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase};
use fs_material::thermovisco::{Limits, ThermalKnot, ThermoMaxwell};
use fs_material::thermovisco::thermal::{CaloricPoint, FeedbackConfig, MechanicalControl};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const LENGTH_M: f64 = 0.1;
const AREA_M2: f64 = 1e-5;
const VOLUME_M3: f64 = LENGTH_M * AREA_M2;
const DENSITY_KG_M3: f64 = 1000.0;
fn surface_area() -> f64 { 4.0 * AREA_M2.sqrt() * LENGTH_M + 2.0 * AREA_M2 }

#[derive(Clone, Copy)]
struct Request { cp: f64, htc: f64, stress: f64, duration: f64, steps: usize }
impl Request {
    fn validate(self) -> Result<Self> {
        if !self.cp.is_finite() || self.cp <= 0.0 || !self.htc.is_finite() || self.htc < 0.0
            || !self.stress.is_finite() || !self.duration.is_finite() || self.duration <= 0.0
            || self.duration > 1000.0 || !(2..=10_000).contains(&self.steps) || self.steps % 2 != 0 {
            return Err("positive Cp/duration, nonnegative film coefficient, finite stress and 2..10000 even steps required; duration <=1000 s".into());
        }
        Ok(self)
    }
}
fn config() -> FeedbackConfig {
    FeedbackConfig { max_iterations: 100, temperature_tolerance_k: 1e-8,
        energy_tolerance_j: 1e-8, relaxation: 1.0 }
}
fn material(cp: f64) -> Result<(ThermoMaxwell, EquilibriumEnthalpyPhaseCurve)> {
    let law = ThermoMaxwell::new(1e8, vec![(3e8, 1.0)], vec![
        ThermalKnot { temperature_k: 280.0, modulus_scale: 1.0, relaxation_scale: 1.0, free_strain: 0.0 },
        ThermalKnot { temperature_k: 300.0, modulus_scale: 1.0, relaxation_scale: 1.0, free_strain: 0.0 },
        ThermalKnot { temperature_k: 320.0, modulus_scale: 1.0, relaxation_scale: 0.05, free_strain: 0.0 },
    ], Limits { max_strain: 0.04, max_strain_rate_s: 10.0 })?;
    // Bind the synthetic caloric input, not invented measured card authority.
    let mut identity = DomainHasher::new("org.frankensim.example.self-heating-creep.v1");
    identity.update(&cp.to_bits().to_le_bytes());
    let knots: Vec<_> = [280.0, 300.0, 320.0].into_iter().map(|temperature_k| HeatCapacityKnot {
        temperature_k, specific_heat_capacity_j_kg_k: cp, bulk_density_kg_m3: DENSITY_KG_M3,
    }).collect();
    let curve = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(identity.finalize(),
        SolidLiquidPhase::Solid, -20.0 * cp, &knots, 1e-7)?;
    Ok((law, curve))
}
fn environment(htc: f64) -> Result<LumpedThermalEnvironment> {
    environment_with_transport(htc, LumpedThermalTransport::try_declared(10.0, 0.0)?)
}
fn environment_with_transport(htc: f64, transport: LumpedThermalTransport) -> Result<LumpedThermalEnvironment> {
    Ok(LumpedThermalEnvironment { temperature_k: 300.0, radiation_temperature_k: 300.0,
        convection_w_per_m2_k: htc, transport,
        // Radiation is explicitly disabled in both input modes. Conductivity
        // still drives the original isothermal-body admission, not an invented gradient.
        enthalpy_tolerance_j_kg: 1e-9, maximum_thermal_residual_j: 1e-10 })
}
fn context<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 0, kernel_id: 0x534843, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        f(&gate, &cx)
    })
}

#[derive(Debug)]
struct Record {
    time: f64, temperature: f64, strain: f64, stress: f64, target_stress: f64,
    physical_heat: f64, numerical_loss: f64, mechanical_work: f64, boundary_heat: f64,
    elastic_storage: f64, caloric_storage: f64, combined_defect: f64, proposals: usize,
    reversible_heat: f64, deposited_heat: f64, mechanical_internal_change: f64,
    mechanical_entropy: f64,
}
fn simulate(request: Request, gate: &CancelGate, cx: &Cx<'_>) -> Result<Vec<Record>> {
    let request = request.validate()?;
    let (law, curve) = material(request.cp)?;
    let environment = environment(request.htc)?;
    let initial = curve.state_at_specific_enthalpy(0.0)?;
    let state = CaloricPoint::try_new(law.birth(initial.temperature_k(), 0.0)?, &curve,
        0.0, DENSITY_KG_M3 * VOLUME_M3, VOLUME_M3, surface_area())?;
    simulate_state(request, state, &environment, gate, cx)
}
fn simulate_source(request: Request, source: &source::Source, gate: &CancelGate, cx: &Cx<'_>) -> Result<Vec<Record>> {
    let state = source.mechanical.caloric_point(&source.caloric,300.0,0.0,VOLUME_M3,surface_area())?;
    let environment = environment_with_transport(request.htc,source.transport.clone())?;
    simulate_state(request.validate()?,state,&environment,gate,cx)
}
fn simulate_state(request: Request, mut state: CaloricPoint<'_>, environment: &LumpedThermalEnvironment,
    gate: &CancelGate, cx: &Cx<'_>) -> Result<Vec<Record>> {
    let initial = *state.caloric_state();
    let initial_mechanical_internal = state.mechanical_internal_energy_j();
    let mut records = Vec::with_capacity(request.steps);
    let (mut work, mut external, mut numerical) = (0.0, 0.0, 0.0);
    for index in 1..=request.steps {
        // The discontinuity is exactly an interval boundary; the remaining
        // half removes the load but does not reset strain, heat or branch memory.
        let target_stress = if index <= request.steps/2 { request.stress } else { 0.0 };
        let target_time = request.duration * (index as f64 / request.steps as f64);
        let step = state.advance_controlled(
            MechanicalControl::Stress { stress_pa: target_stress, tolerance_pa: 1e-6 },
            target_time-state.time_s(), config(), |input| environment.advance(cx, input),
            || gate.is_requested())?;
        work += step.external_mechanical_work_j;
        external += step.transport.external_heat_j;
        numerical += VOLUME_M3 * step.mechanical_energy.numerical_dissipation;
        let next = step.next;
        let elastic = VOLUME_M3 * next.mechanical().stored_energy_j_m3();
        let mechanical_internal_change = next.mechanical_internal_energy_j()
            - initial_mechanical_internal;
        let caloric = next.mass_kg() * (next.caloric_state().specific_enthalpy_j_kg()
            - initial.specific_enthalpy_j_kg());
        // Temperature-dependent elastic storage is Helmholtz free energy.
        // The original material owner also retains its mechanical entropy;
        // the first law uses U = psi + T*s, evaluated from the actual endpoint.
        let defect = mechanical_internal_change + caloric + numerical - work - external;
        if ![work, external, numerical, defect].iter().all(|v| v.is_finite())
            || defect.abs() > index as f64 * config().energy_tolerance_j {
            return Err("cumulative thermo-mechanical energy gate failed".into());
        }
        records.push(Record { time: next.time_s(), temperature: next.caloric_state().temperature_k(),
            strain: next.mechanical().total_strain(), stress: next.mechanical().stress_pa(), target_stress,
            physical_heat: VOLUME_M3 * step.mechanical_energy.dissipated,
            numerical_loss: VOLUME_M3 * step.mechanical_energy.numerical_dissipation,
            mechanical_work: step.external_mechanical_work_j, boundary_heat: step.transport.external_heat_j,
            elastic_storage: elastic, caloric_storage: caloric, combined_defect: defect,
            proposals: step.iterations, reversible_heat: step.reversible_heat_j,
            deposited_heat: step.applied_heat_j, mechanical_internal_change,
            mechanical_entropy: next.mechanical_entropy_j_k() });
        state = next;
    }
    Ok(records)
}
fn run(args: &[String]) -> Result<()> {
    let mode=args.first().map(String::as_str).unwrap_or("");
    let offset=match (mode,args.len()) {
        ("--demo" | "--table",6) => 0,
        ("--pack",7) => 1,
        _ => return Err("usage: self_heating_creep --demo CP|--table FILE HTC STRESS DURATION EVEN_STEPS; or --pack FILE BRANCH_COUNT HTC STRESS DURATION EVEN_STEPS".into()),
    };
    let (records,description)=context(|gate,cx| -> Result<_> {
        let supplied=match mode {
            "--table" => Some(source::Source::table(&args[1],gate)?),
            "--pack" => Some(source::Source::pack(&args[1],args[2].parse()?,gate)?),
            _ => None,
        };
        // In file mode this validator value is the actual initial source Cp;
        // the march still uses the complete caloric curve, not this one scalar.
        let cp=match &supplied { Some(s) => s.caloric.secant_heat_capacity(300.0,300.0)?, None => args[1].parse()? };
        let request=Request {cp,htc:args[2+offset].parse()?,stress:args[3+offset].parse()?,
            duration:args[4+offset].parse()?,steps:args[5+offset].parse()?}.validate()?;
        let records=match &supplied {Some(s)=>simulate_source(request,s,gate,cx)?,None=>simulate(request,gate,cx)?};
        let description=match supplied {
            Some(s)=>format!("material_card={}, branches={}, source_states={}, caloric_curve={}, thermal_transport={}; temperature_span=[{},{}] K; initial_Cp={} J/(kg K); no material qualification inferred",
                s.card.content_hash().to_hex(),s.mechanical.law().branch_count(),s.mechanical.samples().len(),
                s.caloric.curve().identity().to_hex(),s.transport.identity().to_hex(),
                s.caloric.curve().knots()[0].temperature_k,s.caloric.curve().knots().last().unwrap().temperature_k,cp),
            None=>format!("SYNTHETIC data; rho={DENSITY_KG_M3} kg/m3; E_inf=1e8 Pa, branch=(3e8 Pa,1 s); shift=1 at 300 K,0.05 at 320 K; Cp={cp} J/(kg K); k=10 W/(m K)"),
        };
        Ok((records,format!("{description}; L={LENGTH_M} m, A={AREA_M2} m2, surface={} m2; film={} W/(m2 K), ambient=300 K, radiation disabled; load removed at {} s; fixed thermal geometry, small-strain uniform specimen, no yield or spatial/dynamic claim",
            surface_area(),request.htc,request.duration/2.0)))
    })?;
    eprintln!("{description}");
    println!("time_s,temperature_K,total_strain,extension_m,stress_Pa,target_stress_Pa,physical_heat_J,numerical_loss_J,mechanical_work_J,boundary_heat_J,elastic_free_energy_J,caloric_storage_change_J,cumulative_defect_J,thermal_proposals,reversible_heat_J,deposited_heat_J,mechanical_internal_energy_change_J,mechanical_entropy_J_per_K");
    for r in records {
        println!("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}", r.time, r.temperature, r.strain,
            LENGTH_M*r.strain, r.stress, r.target_stress, r.physical_heat, r.numerical_loss,
            r.mechanical_work, r.boundary_heat, r.elastic_storage, r.caloric_storage,
            r.combined_defect, r.proposals, r.reversible_heat, r.deposited_heat,
            r.mechanical_internal_change, r.mechanical_entropy);
    }
    Ok(())
}
fn main() -> Result<()> { run(&std::env::args().skip(1).collect::<Vec<_>>()) }

#[cfg(test)]
mod tests {
    use super::*;
    fn request(cp: f64, htc: f64) -> Request {
        Request { cp, htc, stress: 2e6, duration: 4.8, steps: 48 }
    }
    #[test]
    fn actual_thermal_owner_changes_creep_when_only_caloric_data_changes() {
        context(|gate, cx| {
            let a = simulate(request(5.0, 2.0), gate, cx).unwrap();
            let b = simulate(request(10.0, 2.0), gate, cx).unwrap();
            assert!(a[23].temperature > b[23].temperature + 0.2);
            assert!(a[23].strain > b[23].strain + 1e-5);
            for history in [&a, &b] {
                for (i, row) in history.iter().enumerate() {
                    assert!(row.combined_defect.abs() < (i+1) as f64 * config().energy_tolerance_j);
                    assert!((row.stress-row.target_stress).abs() <= 1e-6);
                    assert!(row.physical_heat >= 0.0 && row.numerical_loss >= 0.0);
                }
                assert!(history[24].strain > 0.0 && history[24].strain < history[23].strain);
                assert!(history[24].temperature > 300.0 && history[24].physical_heat > 0.0);
                assert!(history[24..].iter().all(|row| row.mechanical_work == 0.0));
            }
        });
    }
    #[test]
    fn changing_only_boundary_cooling_changes_the_later_material_response() {
        context(|gate, cx| {
            let insulated = simulate(request(5.0, 0.0), gate, cx).unwrap();
            let cooled = simulate(request(5.0, 20.0), gate, cx).unwrap();
            assert!(insulated[23].temperature > cooled[23].temperature + 0.5);
            assert!(insulated[23].strain > cooled[23].strain + 1e-5);
            assert!(cooled[23].boundary_heat < 0.0);
            assert_eq!(insulated[23].boundary_heat, 0.0);
        });
    }
    #[test]
    fn real_transport_domain_and_cancellation_refusals_keep_both_histories() {
        context(|gate, cx| {
            let (law, curve) = material(5.0).unwrap();
            let state = CaloricPoint::try_new(law.birth(300.0, 0.0).unwrap(), &curve, 0.0,
                DENSITY_KG_M3*VOLUME_M3, VOLUME_M3, surface_area()).unwrap();
            let old = state.clone();
            let invalid = environment(1e8).unwrap(); // Actual Biot refusal, no mock.
            let control = MechanicalControl::Stress { stress_pa: 2e6, tolerance_pa: 1e-6 };
            assert!(matches!(state.advance_controlled(control, 0.1, config(),
                |input| invalid.advance(cx,input), || gate.is_requested()),
                Err(fs_material::thermovisco::thermal::FeedbackError::Transport(_))));
            let valid = environment(2.0).unwrap();
            assert!(matches!(state.advance_controlled(control, 0.1, config(), |input| {
                let result = valid.advance(cx, input);
                gate.request();
                result
            }, || gate.is_requested()), Err(fs_material::thermovisco::thermal::FeedbackError::Cancelled)));
            assert_eq!(state, old);
        });
        assert!(Request { steps: 3, ..request(5.0, 2.0) }.validate().is_err());
    }

    #[test]
    fn material_file_values_reproduce_demo_and_data_only_changes_reach_creep() {
        context(|gate,cx| {
            let low=source::Source::table_text(include_str!("self_heating_creep/cp5.txt"),gate).unwrap();
            let high=source::Source::table_text(include_str!("self_heating_creep/cp10.txt"),gate).unwrap();
            let a=simulate_source(request(5.0,2.0),&low,gate,cx).unwrap();
            let b=simulate_source(request(10.0,2.0),&high,gate,cx).unwrap();
            let demo=simulate(request(5.0,2.0),gate,cx).unwrap();
            for (actual,expected) in a.iter().zip(demo) {
                assert!((actual.temperature-expected.temperature).abs()<1e-6);
                assert!((actual.strain-expected.strain).abs()<1e-9);
                assert!(actual.combined_defect.abs()<1e-7);
            }
            assert!(a[23].temperature>b[23].temperature+0.2);
            assert!(a[23].strain>b[23].strain+1e-5);
            assert_eq!(low.transport.material_card_identity(),Some(low.card.content_hash()));
            assert_eq!(low.caloric.curve().material_card_identity(),low.card.content_hash());
            assert!(!low.transport.conductivity_receipts().is_empty());
        });
    }
    #[test]
    fn affine_material_file_changes_creep_and_closes_internal_energy_with_real_transport() {
        context(|gate,cx| {
            let fixed=source::Source::table_text(include_str!("self_heating_creep/cp5.txt"),gate).unwrap();
            let affine=source::Source::table_text(include_str!("self_heating_creep/cp5_affine.txt"),gate).unwrap();
            let fixed_rows=simulate_source(request(5.0,2.0),&fixed,gate,cx).unwrap();
            let rows=simulate_source(request(5.0,2.0),&affine,gate,cx).unwrap();
            // Both files have the same initial stiffness, heat capacity, density,
            // relaxation law, geometry, load and actual convection owner.
            // Only S(T)=2.5-0.005*T differs from the constant reference file.
            assert_eq!(fixed.mechanical.law().birth(300.0,0.0).unwrap()
                .linearize(300.0,0.1).unwrap().tangent_pa,
                affine.mechanical.law().birth(300.0,0.0).unwrap()
                .linearize(300.0,0.1).unwrap().tangent_pa);
            assert!(rows[0].physical_heat>0.0 && rows[0].reversible_heat<0.0);
            assert!(rows[0].deposited_heat<0.0 && rows[0].temperature<300.0);
            assert!(rows[24].reversible_heat>0.0);
            assert!((rows[23].temperature-fixed_rows[23].temperature).abs()>0.1);
            assert!((rows[23].strain-fixed_rows[23].strain).abs()>1e-5);
            for (index,row) in rows.iter().enumerate() {
                let scale=2.5-0.005*row.temperature;
                // Independently, psi=S*Xi, s=-S'*Xi and U=(S-T*S')*Xi.
                let unscaled=row.elastic_storage/scale;
                assert!((row.mechanical_internal_change-2.5*unscaled).abs()<1e-12);
                assert!((row.mechanical_entropy-0.005*unscaled).abs()<1e-12);
                assert!((row.physical_heat+row.reversible_heat-row.deposited_heat).abs()
                    <=config().energy_tolerance_j);
                assert!(row.physical_heat>=0.0 && row.numerical_loss>=0.0);
                assert!((row.stress-row.target_stress).abs()<=1e-6);
                assert!(row.combined_defect.abs()<(index+1) as f64*config().energy_tolerance_j);
            }
            // Using free energy in the first law would miss this finite term.
            assert!((rows[23].mechanical_internal_change-rows[23].elastic_storage).abs()>1e-3);
        });
    }
    #[test]
    fn sourced_density_changes_reference_heat_storage_without_changing_geometry() {
        context(|gate,cx| {
            let text=include_str!("self_heating_creep/cp5.txt");
            let a=source::Source::table_text(text,gate).unwrap();
            let b=source::Source::table_text(&text.replace("1000 10","2000 10"),gate).unwrap();
            let sa=a.mechanical.caloric_point(&a.caloric,300.0,0.0,VOLUME_M3,surface_area()).unwrap();
            let sb=b.mechanical.caloric_point(&b.caloric,300.0,0.0,VOLUME_M3,surface_area()).unwrap();
            assert_eq!(sb.mass_kg(),2.0*sa.mass_kg());
            let ra=simulate_source(request(5.0,2.0),&a,gate,cx).unwrap();
            let rb=simulate_source(request(5.0,2.0),&b,gate,cx).unwrap();
            assert!(ra[23].temperature>rb[23].temperature+0.2);
            assert!(ra[23].strain>rb[23].strain);
        });
    }
    #[test]
    fn source_file_refusals_and_cooling_admission_remain_real() {
        context(|gate,cx| {
            let text=include_str!("self_heating_creep/cp5.txt");
            for bad in [text.replace("equilibrium 100000000", ""),
                text.replace("state 300", "state 280"),text.replace("1000 10", "1000 NaN"),
                text.replace("state 320 1 0.05 0", "state 320 0.5 0.05 0")] {
                assert!(source::Source::table_text(&bad,gate).is_err());
            }
            let low_k=source::Source::table_text(&text.replace("1000 10","1000 0.000001"),gate).unwrap();
            assert!(simulate_source(request(5.0,2.0),&low_k,gate,cx).is_err());
            gate.request(); assert!(source::Source::table_text(text,gate).is_err());
        });
        assert!(run(&["--demo".into()]).is_err());
        assert!(source::Source::pack("/definitely-absent-material.fsmat",1,&CancelGate::new()).is_err());
    }
}
