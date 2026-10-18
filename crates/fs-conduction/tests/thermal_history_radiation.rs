//! Scheduled radiative reservoirs, physical heat splits and transactional resume.
mod support;

use fs_conduction::bc::history::{ThermalBoundaryHistory, ThermalRegionHistory};
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryMethod, ThermalHistoryProblem, ThermalHistoryRun};
use fs_conduction::bc::history::march::radiation::{RadiationRegionHistory, RadiativeThermalHistoryRun};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
use fs_conduction::{AmbientRadiationConfig, ConductionError, ConductionMesh, ConductionProblem,
    ConductivityModel, ConductivityTable, EMISSIVITY_DIMS, HEAT_FLUX_DIMS, LinearConfig,
    ScalarField, SurfaceEmissivity, SURFACE_EMISSIVITY_PROPERTY, TEMPERATURE_DIMS,
    ThermalBc, ThermalBoundaryBuilder};
use fs_evidence::ValidityDomain;
use fs_exec::Cx;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, SelectionPolicy, UncertaintyModel};
use fs_qty::QtyAny;
use fs_rep_mesh::TetComplex;
use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition, Physics};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{Interp, TimeSignal};
use support::{with_cancelled_cx, with_cx};

fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]).unwrap()
}
fn controls() -> StepConfig {
    StepConfig { linear: LinearConfig { tolerance: 1e-12, ..LinearConfig::default() },
        energy_tolerance_j: 1e-7 }
}
fn emissivity() -> SurfaceEmissivity {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Scalar { value: 0.8, dims: EMISSIVITY_DIMS },
        validity: ValidityDomain::unconstrained().with("T", 250.0, 500.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: Provenance { source: "radiative history fixture".into(),
            license: "internal-test-use".into(), artifact: None },
    }).unwrap();
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "gray-surface".into(),
        phase: "solid".into(), process: "test".into(), revision: 0 }, claims, Vec::new()).unwrap();
    SurfaceEmissivity::from_card("surface", &card, 350.0, SelectionPolicy::SingleClaimOnly).unwrap()
}
fn reservoir(region: &str, signal: TimeSignal) -> RadiationRegionHistory {
    RadiationRegionHistory { region: region.into(), emissivity: emissivity(), ambient_temperature: signal }
}
fn constant(value: f64) -> TimeSignal { TimeSignal::Constant(QtyAny::new(value, TEMPERATURE_DIMS)) }
fn driven<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh) -> ThermalBoundaryHistory<'m> {
    let builder = ThermalBoundaryBuilder::new(mesh)
        .region("driven", |f| (f.centroid.iter().sum::<f64>()-1.0).abs()<1e-12,
            ThermalBc::dirichlet(300.0).unwrap()).unwrap()
        .remainder("surface", ThermalBc::robin(5.0, 300.0).unwrap()).unwrap();
    let row = BoundaryCondition { region: "driven".into(), physics: Physics::Thermal,
        kind: BcKind::Dirichlet, compatibility: None, frame: 0,
        value: Some(BcValue::Signal(TimeSignal::Ramp { t_start: 0.0, t_end: 1.0,
            from: QtyAny::new(300.0, TEMPERATURE_DIMS), to: QtyAny::new(360.0, TEMPERATURE_DIMS) })) };
    ThermalBoundaryHistory::new(cx, builder, &Environment::earth_lab(),
        vec![ThermalRegionHistory { row, reference_temperature: None }], 0.0).unwrap()
}
fn run<'h, 'm>(cx: &Cx<'_>, history: &'h ThermalBoundaryHistory<'m>, material: &'h ConductivityModel,
    initial: &[f64], max_steps: usize, method: ThermalHistoryMethod) -> ThermalHistoryRun<'h, 'm>
{
    ThermalHistoryRun::new_with_method(cx, ThermalHistoryProblem { boundary: history, material,
        element_materials: None, interfaces: None, source: HeatSourceHistory::Fixed(ScalarField::Uniform(0.0)) },
        VolumetricHeatCapacity::declared(15.0).unwrap(), initial, max_steps, method).unwrap()
}

#[test]
fn g1_radiative_histories_match_direct_endpoints_and_keep_material_and_heat_reports() {
    let mesh = mesh();
    let linear = ConductivityModel::isotropic_declared(20.0).unwrap();
    let curved = ConductivityModel::isotropic(ConductivityTable::declared_curve(
        vec![(250.0, 15.0), (500.0, 40.0)]).unwrap());
    with_cx(|cx| {
        let history = driven(cx, &mesh);
        for (material, policy) in [(&linear, None), (&curved, Some(NonlinearStepConfig::default()))] {
            let method = policy.map_or(ThermalHistoryMethod::Linear, ThermalHistoryMethod::Nonlinear);
            let patch = reservoir("surface", TimeSignal::Table {
                times: vec![0.0, 0.5, 0.75], values: vec![300.0, 450.0, 300.0],
                dims: TEMPERATURE_DIMS, interp: Interp::Hold });
            let mut driver = RadiativeThermalHistoryRun::new(cx,
                run(cx, &history, material, &[300.0; 4], 3, method),
                vec![patch.clone()], AmbientRadiationConfig::default()).unwrap();
            let manual = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(15.0).unwrap()).unwrap();
            let mut old = vec![300.0; 4];
            let mut previous = 0.0;
            for (time, ambient) in [(0.25, 300.0), (0.5, 450.0), (0.75, 300.0)] {
                let boundary = history.at(cx, time).unwrap();
                let source = ScalarField::Uniform(0.0);
                let expected = manual.advance_prescribed_with_ambient_radiation(cx, ConductionProblem {
                    mesh: &mesh, boundary: &boundary, material, source: &source, element_materials: None },
                    None, &old, time-previous, controls(), policy, &[patch.at(time).unwrap()],
                    AmbientRadiationConfig::default()).unwrap();
                let actual = driver.advance_to(cx, time, controls()).unwrap();
                assert_eq!(actual.conduction.temperature, expected.conduction.temperature);
                assert_eq!(actual.radiation, expected.radiation);
                assert_eq!(actual.physical_dirichlet_in_w, expected.physical_dirichlet_in_w);
                assert_eq!(actual.physical_energy_residual_j, expected.physical_energy_residual_j);
                assert!(actual.physical_energy_residual_j.abs() <= controls().energy_tolerance_j);
                assert!(actual.physical_residual_norm_j <= actual.physical_residual_tolerance_j);
                assert_eq!(actual.nonlinear.is_some(), policy.is_some());
                assert_eq!(actual.radiation.patches[0].patch.ambient_temperature_k(), ambient);
                assert_eq!(actual.radiation.patches[0].patch.emissivity(), &patch.emissivity);
                assert_eq!(actual.convective_robin_fluxes[0].mean_reference_temperature_k, 300.0);
                assert_eq!(actual.radiation.nonlinear_radiation_out_w < 0.0, ambient == 450.0);
                old = expected.conduction.temperature;
                previous = time;
            }
            assert_eq!(driver.conduction().temperature(), old.as_slice());
            assert_eq!(driver.conduction().accepted_steps(), 3);
            assert_eq!(driver.patches(), &[patch]);
        }
    });
}

#[test]
fn g1_vacuum_facing_history_cools_without_inventing_convection() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(20.0).unwrap();
    with_cx(|cx| {
        let history = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .remainder("surface", ThermalBc::adiabatic()).unwrap(), &Environment::earth_lab(),
            Vec::new(), 0.0).unwrap();
        let mut driver = RadiativeThermalHistoryRun::new(cx,
            run(cx, &history, &material, &[350.0; 4], 2, ThermalHistoryMethod::Linear),
            vec![reservoir("surface", constant(300.0))], AmbientRadiationConfig::default()).unwrap();
        for time in [0.25, 0.5] {
            let step = driver.advance_to(cx, time, controls()).unwrap();
            assert_eq!(step.convective_out_w, 0.0);
            assert!(step.convective_robin_fluxes.is_empty());
            assert!(step.radiation.nonlinear_radiation_out_w > 0.0);
            assert!(step.conduction.stored_energy_change_j < 0.0);
            assert!(step.physical_energy_residual_j.abs() <= controls().energy_tolerance_j);
        }
    });
}

#[test]
fn g0_bad_reservoir_cancel_and_budget_refusals_preserve_radiative_resume() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(20.0).unwrap();
    let history = with_cx(|cx| driven(cx, &mesh));
    let patch = reservoir("surface", TimeSignal::Table {
        times: vec![0.0, 1.0], values: vec![300.0, -1.0], dims: TEMPERATURE_DIMS, interp: Interp::Hold });
    let mut driver = with_cx(|cx| RadiativeThermalHistoryRun::new(cx,
        run(cx, &history, &material, &[300.0; 4], 1, ThermalHistoryMethod::Linear),
        vec![patch.clone()], AmbientRadiationConfig::default()).unwrap());
    with_cancelled_cx(|cx| assert!(matches!(driver.advance_to(cx, 0.25, controls()),
        Err(ConductionError::Cancelled { .. }))));
    assert_eq!(driver.conduction().time_s(), 0.0);
    with_cx(|cx| {
        driver.advance_to(cx, 0.25, controls()).unwrap();
        let accepted = driver.conduction().temperature().to_vec();
        assert!(driver.advance_to(cx, 0.5, controls()).is_err());
        driver.extend_step_budget(2).unwrap();
        assert!(driver.advance_to(cx, 1.0, controls()).is_err());
        assert_eq!(driver.conduction().temperature(), accepted.as_slice());
        assert_eq!(driver.conduction().time_s(), 0.25);
        assert_eq!(driver.conduction().accepted_steps(), 1);
        driver.advance_to(cx, 0.5, controls()).unwrap();
        let mut direct = RadiativeThermalHistoryRun::new(cx,
            run(cx, &history, &material, &[300.0; 4], 2, ThermalHistoryMethod::Linear),
            vec![patch], AmbientRadiationConfig::default()).unwrap();
        direct.advance_to(cx, 0.25, controls()).unwrap();
        direct.advance_to(cx, 0.5, controls()).unwrap();
        assert_eq!(driver.conduction().temperature(), direct.conduction().temperature());
    });
}

#[test]
fn g0_binding_and_outer_solve_refusals_do_not_enable_an_incomplete_radiation_model() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(20.0).unwrap();
    with_cx(|cx| {
        let history = driven(cx, &mesh);
        let base = || run(cx, &history, &material, &[300.0; 4], 2, ThermalHistoryMethod::Linear);
        let good = reservoir("surface", constant(300.0));
        let bad_units = reservoir("surface", TimeSignal::Constant(QtyAny::new(300.0, HEAT_FLUX_DIMS)));
        for patches in [Vec::new(), vec![reservoir("absent", constant(300.0))],
            vec![reservoir("driven", constant(300.0))], vec![good.clone(), good.clone()], vec![bad_units]] {
            assert!(RadiativeThermalHistoryRun::new(cx, base(), patches, AmbientRadiationConfig::default()).is_err());
        }
        let mut advanced = base();
        advanced.advance_to(cx, 0.25, controls()).unwrap();
        assert!(RadiativeThermalHistoryRun::new(cx, advanced, vec![good.clone()], AmbientRadiationConfig::default()).is_err());
        let mut limited = RadiativeThermalHistoryRun::new(cx, base(), vec![good],
            AmbientRadiationConfig { max_iterations: 1, ..AmbientRadiationConfig::default() }).unwrap();
        assert!(limited.advance_to(cx, 0.25, controls()).is_err());
        assert_eq!(limited.conduction().temperature(), &[300.0; 4]);
        assert_eq!(limited.conduction().time_s(), 0.0);
        assert_eq!(limited.conduction().accepted_steps(), 0);
    });
}
