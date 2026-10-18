//! Short source/reservoir pulses must survive latent and radiative stepping.
mod support;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryMethod, ThermalHistoryProblem, ThermalHistoryRun};
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryBudget, EnthalpyHistoryMaterials, EnthalpyHistoryRun};
use fs_conduction::bc::history::march::enthalpy::radiation::RadiativeEnthalpyHistoryRun;
use fs_conduction::bc::history::march::radiation::{RadiationRegionHistory, RadiativeThermalHistoryRun};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::ReferenceEnthalpyMaterial;
use fs_conduction::{AmbientRadiationConfig, AmbientRadiationPatch, ConductionMesh, ConductionProblem,
    ConductivityModel, ConductivityTable, EMISSIVITY_DIMS, LinearConfig, ScalarField, SurfaceEmissivity,
    SURFACE_EMISSIVITY_PROPERTY, TEMPERATURE_DIMS, VOLUMETRIC_SOURCE_DIMS, ThermalBc,
    ThermalBoundaryBuilder, STEFAN_BOLTZMANN_W_M2_K4};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::Cx;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, SelectionPolicy, UncertaintyModel};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_qty::Dims;
use fs_rep_mesh::TetComplex;
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{Interp, TimeSignal};
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(4, vec![[0,1,2,3]]),
        vec![[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.]]).unwrap()
}
fn history<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(mesh)
        .remainder("surface", ThermalBc::adiabatic()).unwrap(), &Environment::earth_lab(), vec![], 0.).unwrap()
}
fn held(times: Vec<f64>, values: Vec<f64>, dims: Dims) -> TimeSignal {
    TimeSignal::Table { times, values, dims, interp: Interp::Hold }
}
fn curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([73; 32]),
        [(-100.,300.,0.), (0.,310.,0.), (1000.,310.,1.), (1100.,320.,1.)].into_iter()
            .map(|(h,t,f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg: h, temperature_k: t,
                liquid_mass_fraction: f, bulk_density_kg_m3: 1. }).collect()).unwrap()
}
fn config() -> EnthalpyStepConfig {
    EnthalpyStepConfig { newton: NewtonKrylovConfig { absolute_tolerance: 1e-10,
        relative_tolerance: 1e-12, linear_restart: 4, max_linear_cycles: 16,
        forcing_minimum: 1e-12, forcing_maximum: 1e-3, ..Default::default() },
        max_newton_iterations: 64, energy_tolerance_j: 1e-7 }
}
fn problem<'h, 'm>(boundary: &'h ThermalBoundaryHistory<'m>, material: &'h ConductivityModel,
    source: HeatSourceHistory) -> ThermalHistoryProblem<'h, 'm>
{
    ThermalHistoryProblem { boundary, material, source, element_materials: None, interfaces: None }
}
fn phase_run<'h, 'm, 'c>(cx: &Cx<'_>, problem: ThermalHistoryProblem<'h, 'm>,
    curve: &'c EquilibriumEnthalpyPhaseCurve, density: f64, max_steps: usize) -> EnthalpyHistoryRun<'h, 'm, 'c>
{
    EnthalpyHistoryRun::new(cx, problem, EnthalpyHistoryMaterials {
        table: &[ReferenceEnthalpyMaterial { curve, reference_density_kg_m3: density }],
        of_element: &[0], phase_conductivity: None }, &[200.; 4],
        EnthalpyHistoryBudget { spatial: EnthalpyBudget { max_vertices: 4, max_elements: 1 }, max_steps }).unwrap()
}
fn zero() -> HeatSourceHistory { HeatSourceHistory::Fixed(ScalarField::Uniform(0.)) }
fn patch(final_temperature: f64) -> RadiationRegionHistory {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim { key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Scalar { value: 0.8, dims: EMISSIVITY_DIMS },
        validity: ValidityDomain::unconstrained().with("T", 250., 500.), uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
        provenance: Provenance { source: "synthetic pulse fixture".into(), license: "internal-test-use".into(), artifact: None } }).unwrap();
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "gray body".into(), phase: "solid".into(),
        process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
    RadiationRegionHistory { region: "surface".into(),
        emissivity: SurfaceEmissivity::from_card("surface", &card, 310., SelectionPolicy::SingleClaimOnly).unwrap(),
        ambient_temperature: held(vec![0.,0.125,0.25], vec![310.,350.,final_temperature], TEMPERATURE_DIMS) }
}

#[test]
fn g1_short_heating_and_cooling_pulses_change_phase_not_plateau_temperature() {
    let mesh = mesh(); let material = ConductivityModel::isotropic_declared(20.).unwrap(); let curve = curve();
    with_cx(|cx| {
        let history = history(cx, &mesh);
        let source = HeatSourceHistory::Uniform(held(vec![0.,0.125,0.25,0.5,0.625],
            vec![0.,120.,0.,-120.,0.], VOLUMETRIC_SOURCE_DIMS));
        let mut driver = phase_run(cx, problem(&history, &material, source), &curve, 1., 5);
        for (end, expected_h) in [(0.125,200.),(0.25,215.),(0.5,215.),(0.625,200.),(1.,200.)] {
            let step = driver.advance_event_step(cx, 1., 1., config()).unwrap();
            assert_eq!(step.end_time_s, end);
            assert!(step.report.specific_enthalpy_j_kg.iter().all(|&h| (h-expected_h).abs()<1e-8));
            assert!(step.report.temperature.iter().all(|&t| t==310.));
            assert!(step.report.liquid_mass_fraction.iter().all(|&f| (f-expected_h/1000.).abs()<1e-10));
            assert!(step.report.energy_residual_j.abs()<=config().energy_tolerance_j);
        }
    });
}

#[test]
fn g1_radiative_temperature_events_match_manual_pre_switch_reservoirs() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(20.).unwrap();
    let nonlinear = ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(250.,15.),(500.,40.)]).unwrap());
    with_cx(|cx| {
        let history = history(cx, &mesh);
        let config = StepConfig { linear: LinearConfig { tolerance: 1e-12, ..Default::default() }, energy_tolerance_j: 1e-7 };
        for (material, policy) in [(&material, None), (&nonlinear, Some(NonlinearStepConfig::default()))] {
            let patch = patch(310.);
            let capacity = VolumetricHeatCapacity::declared(30.).unwrap();
            let method = policy.map_or(ThermalHistoryMethod::Linear, ThermalHistoryMethod::Nonlinear);
            let base = ThermalHistoryRun::new_with_method(cx, problem(&history, material, zero()), capacity, &[310.;4], 3, method).unwrap();
            let mut driver = RadiativeThermalHistoryRun::new(cx, base, vec![patch.clone()], AmbientRadiationConfig::default()).unwrap();
            let manual = BackwardEuler::uniform(cx, &mesh, capacity).unwrap();
            let mut old = vec![310.;4]; let mut start = 0.;
            for (end, ambient) in [(0.125,310.),(0.25,350.),(1.,310.)] {
                let boundary = history.at(cx, end).unwrap(); let source = ScalarField::Uniform(0.);
                let sampled = AmbientRadiationPatch::new("surface", patch.emissivity.clone(), ambient).unwrap();
                let expected = manual.advance_prescribed_with_ambient_radiation(cx, ConductionProblem {
                    mesh: &mesh, boundary: &boundary, material, source: &source, element_materials: None },
                    None, &old, end-start, config, policy, &[sampled], AmbientRadiationConfig::default()).unwrap();
                let step = driver.advance_event_step(cx, 1., 1., config).unwrap();
                assert_eq!(step.end_time_s, end);
                assert_eq!(step.report.conduction.temperature, expected.conduction.temperature);
                assert_eq!(step.report.radiation, expected.radiation);
                assert_eq!(step.report.convective_out_w, 0.);
                assert!(step.report.convective_robin_fluxes.is_empty());
                old = expected.conduction.temperature; start=end;
            }
        }
    });
}

#[test]
fn g1_radiative_pulse_deposits_analytic_latent_energy_without_temperature_rise() {
    let mesh = mesh(); let material = ConductivityModel::isotropic_declared(20.).unwrap(); let curve = curve();
    with_cx(|cx| {
        let history = history(cx, &mesh);
        let base = phase_run(cx, problem(&history, &material, zero()), &curve, 10., 3);
        let mut driver = RadiativeEnthalpyHistoryRun::new(cx, base, vec![patch(310.)], AmbientRadiationConfig::default()).unwrap();
        let mut stored = 0.;
        for (end, ambient) in [(0.125,310.),(0.25,350.),(1.,310.)] {
            let step = driver.advance_event_step(cx, 1., 1., config()).unwrap();
            assert_eq!(step.end_time_s, end);
            assert_eq!(step.report.radiation.patches[0].patch.ambient_temperature_k(), ambient);
            assert!(step.report.conduction.temperature.iter().all(|&t| t==310.));
            assert!(step.report.physical_energy_residual_j.abs()<=config().energy_tolerance_j);
            stored += step.report.conduction.stored_energy_change_j;
        }
        let area = (3.0+3.0_f64.sqrt())/2.0;
        let expected = 0.125*0.8*STEFAN_BOLTZMANN_W_M2_K4*area*(350.0_f64.powi(4)-310.0_f64.powi(4));
        assert!((stored-expected).abs()<1e-7);
        let state = driver.conduction();
        let direct: f64 = state.specific_enthalpy_j_kg().iter().zip(state.storage().reference_nodal_masses_kg())
            .map(|(&h,&m)|m*(h-200.)).sum();
        assert!((direct-expected).abs()<1e-7);
        assert!(state.liquid_mass_fraction().iter().all(|&f| f>0.2 && f<0.4));
    });
}

#[test]
fn g0_invalid_future_reservoir_keeps_completed_pulse_and_budget_resume() {
    let mesh=mesh(); let material=ConductivityModel::isotropic_declared(20.).unwrap(); let curve=curve();
    with_cx(|cx| {
        let history=history(cx,&mesh);
        let base=phase_run(cx,problem(&history,&material,zero()),&curve,10.,1);
        let mut driver=RadiativeEnthalpyHistoryRun::new(cx,base,vec![patch(-1.)],AmbientRadiationConfig::default()).unwrap();
        with_cancelled_cx(|cancel| assert!(driver.advance_event_step(cancel,1.,1.,config()).is_err()));
        assert_eq!(driver.conduction().time_s(),0.);
        driver.advance_event_step(cx,1.,1.,config()).unwrap();
        assert!(driver.advance_event_step(cx,1.,1.,config()).is_err());
        driver.extend_step_budget(3).unwrap();
        // The valid hot interval is accepted even though the next interval is bad.
        driver.advance_event_step(cx,1.,1.,config()).unwrap();
        let accepted=driver.conduction().specific_enthalpy_j_kg().to_vec();
        assert!(driver.advance_event_step(cx,1.,1.,config()).is_err());
        assert_eq!(driver.conduction().specific_enthalpy_j_kg(),accepted.as_slice());
        assert_eq!(driver.conduction().time_s(),0.25);
        assert_eq!(driver.conduction().accepted_steps(),2);
    });
}
