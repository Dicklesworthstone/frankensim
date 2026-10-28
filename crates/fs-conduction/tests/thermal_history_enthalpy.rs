//! Physical latent storage, material ownership and transactional schedule execution.
#![allow(clippy::float_cmp)] // Exact published-state, plateau and resume comparisons.
mod support;

#[path = "thermal_history_enthalpy/radiation.rs"]
mod radiation;

use fs_blake3::ContentHash;
use fs_conduction::bc::history::{ThermalBoundaryHistory, ThermalRegionHistory};
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{
    EnthalpyHistoryBudget, EnthalpyHistoryMaterials, EnthalpyHistoryRun,
};
use fs_conduction::material::LiquidMassFractionConductivity;
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyError, EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, HeterogeneousEnthalpyError, ReferenceEnthalpyMaterial,
};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    ConductivityTable, ScalarField, ThermalBc, ThermalBoundaryBuilder,
    HEAT_FLUX_DIMS, HTC_DIMS, TEMPERATURE_DIMS, VOLUMETRIC_SOURCE_DIMS};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_qty::{Dims, QtyAny};
use fs_rep_mesh::TetComplex;
use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition, Physics};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{Interp, TimeSignal};
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

fn mesh(pair: bool) -> ConductionMesh {
    let mut positions = vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
    let mut tets = vec![[0,1,2,3]];
    if pair {
        positions.extend([[3.,0.,0.],[4.,0.,0.],[3.,1.,0.],[3.,0.,1.]]);
        tets.push([4,5,6,7]);
    }
    ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap()
}
fn curve(latent: f64, id: u8) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([id;32]),
        [(-200.,280.,0.),(100.,310.,0.),(100.+latent,310.,1.),(800.+latent,380.,1.)]
            .into_iter().map(|(h,t,f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg:h,
                temperature_k:t, liquid_mass_fraction:f, bulk_density_kg_m3:9. }).collect()).unwrap()
}
fn config() -> EnthalpyStepConfig {
    EnthalpyStepConfig { newton: NewtonKrylovConfig { absolute_tolerance:1e-11,
        relative_tolerance:1e-12, linear_restart:8, max_linear_cycles:16,
        forcing_minimum:1e-12, forcing_maximum:1e-3, ..Default::default() },
        max_newton_iterations:64, energy_tolerance_j:1e-8 }
}
fn budget(mesh: &ConductionMesh, max_steps: usize) -> EnthalpyHistoryBudget {
    EnthalpyHistoryBudget { spatial: EnthalpyBudget {
        max_vertices:mesh.vertex_count(), max_elements:mesh.element_count() }, max_steps }
}
fn problem<'h,'m>(boundary: &'h ThermalBoundaryHistory<'m>, material: &'h ConductivityModel,
    source: HeatSourceHistory) -> ThermalHistoryProblem<'h,'m>
{
    ThermalHistoryProblem { boundary, material, source, element_materials:None, interfaces:None }
}
fn ramp(from: f64, to: f64, dims: Dims) -> TimeSignal {
    TimeSignal::Ramp { t_start:0., t_end:2., from:QtyAny::new(from,dims), to:QtyAny::new(to,dims) }
}
fn natural<'m>(cx: &fs_exec::Cx<'_>, mesh: &'m ConductionMesh) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(mesh).adiabatic_remainder(),
        &Environment::earth_lab(), Vec::new(), 0.).unwrap()
}

#[test]
fn g1_scheduled_heat_cycle_retains_distinct_material_latent_storage_and_mass() {
    let mesh = mesh(true);
    let curves = [curve(100.,41), curve(400.,42)];
    let table = [ReferenceEnthalpyMaterial { curve:&curves[0], reference_density_kg_m3:1. },
        ReferenceEnthalpyMaterial { curve:&curves[1], reference_density_kg_m3:2. }];
    let conductivity = ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        let boundary = natural(cx,&mesh);
        let source = HeatSourceHistory::Uniform(TimeSignal::Table {
            times:vec![0.,1.5], values:vec![150.,-150.], dims:VOLUMETRIC_SOURCE_DIMS, interp:Interp::Hold });
        let mut run = EnthalpyHistoryRun::new(cx, problem(&boundary,&conductivity,source),
            EnthalpyHistoryMaterials { table:&table, of_element:&[0,1], phase_conductivity:None },
            &[80.;8], budget(&mesh,10)).unwrap();
        let (mut expected_h, mut total) = ([80.,80.],0.);
        for ordinal in 1..=10 {
            let time = f64::from(ordinal)*0.25;
            let q = if time<1.5 {150.} else {-150.};
            let step = run.advance_to(cx,time,config()).unwrap();
            for (body, (rho, latent)) in [(1.,100.),(2.,400.)].into_iter().enumerate() {
                expected_h[body] += q*0.25/rho;
                let h = expected_h[body];
                let (t,f) = if h<100. {(300.+h/10.,0.)} else if h<100.+latent {
                    (310.,(h-100.)/latent)
                } else {(310.+(h-100.-latent)/10.,1.)};
                for node in body*4..body*4+4 {
                    assert!((step.specific_enthalpy_j_kg[node]-h).abs()<1e-8);
                    assert!((step.temperature[node]-t).abs()<1e-9);
                    assert!((step.liquid_mass_fraction[node]-f).abs()<1e-9);
                    assert!((run.storage().reference_nodal_masses_kg()[node]-rho/24.).abs()<1e-15);
                }
            }
            if ordinal==5 {
                assert_eq!(run.liquid_mass_fraction()[0],1.);
                assert!(run.liquid_mass_fraction()[4]>0. && run.liquid_mass_fraction()[4]<1.);
                assert!(run.temperature()[0]>310.);
                assert_eq!(run.temperature()[4],310.);
            }
            assert!((step.source_w-q/3.).abs()<1e-11);
            assert!((step.stored_energy_change_j-0.25*q/3.).abs()<1e-8);
            assert!(step.energy_residual_j.abs()<=config().energy_tolerance_j);
            total+=step.stored_energy_change_j;
        }
        assert!(run.specific_enthalpy_j_kg().iter().all(|h| (h-80.).abs()<1e-8));
        assert!(total.abs()<1e-8);
        assert_eq!(run.storage().element_material_ids(), &[0,1]);
        assert_eq!(run.accepted_steps(),10);
    });
}

#[test]
fn g3_nonlinear_phase_conductivity_and_robin_histories_match_original_endpoint_solver() {
    let mesh = mesh(false);
    let curve = curve(100.,43);
    let table = [ReferenceEnthalpyMaterial { curve:&curve, reference_density_kg_m3:1. }];
    let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(280.,1.6),(380.,3.6)]).unwrap());
    let phase = LiquidMassFractionConductivity::declared(1.,3.).unwrap();
    with_cx(|cx| {
        let row = BoundaryCondition { region:"air".into(), physics:Physics::Thermal,
            kind:BcKind::Robin, value:Some(BcValue::Signal(ramp(10.,30.,HTC_DIMS))), compatibility:None, frame:0 };
        let history = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .remainder("air",ThermalBc::robin(10.,290.).unwrap()).unwrap(), &Environment::earth_lab(),
            vec![ThermalRegionHistory { row, reference_temperature:Some(ramp(290.,310.,TEMPERATURE_DIMS)) }],0.).unwrap();
        let footprint = vec![40.,60.,100.,80.];
        let mut old = vec![80.,135.,155.,240.];
        let mut run = EnthalpyHistoryRun::new(cx, problem(&history,&material,
            HeatSourceHistory::Scaled { field:ScalarField::Nodal(footprint.clone()), scale:ramp(1.,2.,Dims([0;6])) }),
            EnthalpyHistoryMaterials { table:&table, of_element:&[0], phase_conductivity:Some(vec![phase]) },
            &old,budget(&mesh,3)).unwrap();
        let direct = HeterogeneousEnthalpyBackwardEuler::new(cx,&mesh,&table,&[0],budget(&mesh,3).spatial).unwrap()
            .with_phase_conductivity(cx,vec![phase]).unwrap();
        let mut previous=0.;
        for time in [0.0625,0.125,0.25] {
            let boundary=ThermalBoundaryBuilder::new(&mesh)
                .remainder("air",ThermalBc::robin(10.+10.*time,290.+10.*time).unwrap()).unwrap().finish().unwrap();
            let source=ScalarField::Nodal(footprint.iter().map(|q| q*(1.+0.5*time)).collect());
            let expected=direct.advance(cx,ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,
                element_materials:None,source:&source},None,&old,time-previous,config()).unwrap();
            let actual=run.advance_to(cx,time,config()).unwrap();
            assert_eq!(actual.specific_enthalpy_j_kg,expected.specific_enthalpy_j_kg);
            assert_eq!(actual.temperature,expected.temperature);
            assert_eq!(actual.liquid_mass_fraction,expected.liquid_mass_fraction);
            assert_eq!(actual.robin_out_w,expected.robin_out_w);
            assert_eq!(actual.energy_residual_j,expected.energy_residual_j);
            old=expected.specific_enthalpy_j_kg;
            previous=time;
        }
    });
}

#[test]
fn g0_cancel_budget_bad_endpoint_and_material_refusals_preserve_all_accepted_fields() {
    let mesh=mesh(false); let curve=curve(100.,44);
    let material=ConductivityModel::isotropic_declared(2.).unwrap();
    let table=[ReferenceEnthalpyMaterial {curve:&curve,reference_density_kg_m3:1.}];
    let history=with_cx(|cx| natural(cx,&mesh));
    let create=|cx: &fs_exec::Cx<'_>, max_steps| EnthalpyHistoryRun::new(cx,
        problem(&history,&material,HeatSourceHistory::Fixed(ScalarField::Uniform(150.))),
        EnthalpyHistoryMaterials {table:&table,of_element:&[0],phase_conductivity:None},
        &[80.;4],budget(&mesh,max_steps)).unwrap();
    let mut run=with_cx(|cx| create(cx,1));
    with_cancelled_cx(|cx| assert!(matches!(run.advance_to(cx,0.25,config()),
        Err(HeterogeneousEnthalpyError::Enthalpy(EnthalpyError::Conduction(ConductionError::Cancelled {..}))))));
    assert_eq!(run.accepted_steps(),0);
    with_cx(|cx| {
        run.advance_to(cx,0.25,config()).unwrap();
        let old_h=run.specific_enthalpy_j_kg().to_vec();
        let old_t=run.temperature().to_vec(); let old_f=run.liquid_mass_fraction().to_vec();
        assert!(run.advance_to(cx,0.5,config()).is_err());
        run.extend_step_budget(2).unwrap();
        for time in [f64::NAN,f64::INFINITY,0.25,-1.,100.] {
            assert!(run.advance_to(cx,time,config()).is_err());
            assert_eq!(run.specific_enthalpy_j_kg(),old_h.as_slice());
            assert_eq!(run.temperature(),old_t.as_slice());
            assert_eq!(run.liquid_mass_fraction(),old_f.as_slice());
            assert_eq!(run.time_s(),0.25); assert_eq!(run.accepted_steps(),1);
        }
        run.advance_to(cx,0.5,config()).unwrap();
        let mut direct=create(cx,2);
        direct.advance_to(cx,0.25,config()).unwrap(); direct.advance_to(cx,0.5,config()).unwrap();
        assert_eq!(run.specific_enthalpy_j_kg(),direct.specific_enthalpy_j_kg());
        assert_eq!(run.liquid_mass_fraction(),direct.liquid_mass_fraction());
    });
}

#[test]
fn g0_ambiguous_temperature_constraints_and_invalid_initial_enthalpy_are_refused() {
    let mesh=mesh(false); let curve=curve(100.,45);
    let material=ConductivityModel::isotropic_declared(2.).unwrap();
    let table=[ReferenceEnthalpyMaterial {curve:&curve,reference_density_kg_m3:1.}];
    with_cx(|cx| {
        let fixed=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&mesh)
            .remainder("fixed",ThermalBc::dirichlet(310.).unwrap()).unwrap(),
            &Environment::earth_lab(),Vec::new(),0.).unwrap();
        assert!(matches!(EnthalpyHistoryRun::new(cx,problem(&fixed,&material,
            HeatSourceHistory::Fixed(ScalarField::Uniform(0.))),
            EnthalpyHistoryMaterials {table:&table,of_element:&[0],phase_conductivity:None},
            &[150.;4],budget(&mesh,1)), Err(HeterogeneousEnthalpyError::Enthalpy(EnthalpyError::AmbiguousPrescribedTemperature { .. }))));
        let boundary=natural(cx,&mesh);
        for bad in [f64::NAN,f64::INFINITY,-201.,1001.] {
            assert!(EnthalpyHistoryRun::new(cx,problem(&boundary,&material,
                HeatSourceHistory::Fixed(ScalarField::Uniform(0.))),
                EnthalpyHistoryMaterials {table:&table,of_element:&[0],phase_conductivity:None},
                &[bad;4],budget(&mesh,1)).is_err());
        }
        assert!(EnthalpyHistoryRun::new(cx,problem(&boundary,&material,
            HeatSourceHistory::Uniform(TimeSignal::Constant(QtyAny::new(150.,HEAT_FLUX_DIMS)))),
            EnthalpyHistoryMaterials {table:&table,of_element:&[0],phase_conductivity:None},
            &[80.;4],budget(&mesh,1)).is_err());
    });
}

#[test]
fn g1_prescribed_history_preserves_initial_phase_and_charges_boundary_storage_once() {
    let mesh = mesh(false);
    let chart = curve(100., 46);
    let table = [ReferenceEnthalpyMaterial {
        curve: &chart, reference_density_kg_m3: 2.,
    }];
    let material = ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        let boundary = ThermalBoundaryHistory::new(
            cx,
            ThermalBoundaryBuilder::new(&mesh)
                .remainder("bath", ThermalBc::dirichlet(330.).unwrap()).unwrap(),
            &Environment::earth_lab(), Vec::new(), 0.,
        ).unwrap();
        let initial_h = [150.; 4];
        let mut run = EnthalpyHistoryRun::new(
            cx,
            problem(&boundary, &material,
                HeatSourceHistory::Fixed(ScalarField::Uniform(30.))),
            EnthalpyHistoryMaterials {
                table: &table, of_element: &[0], phase_conductivity: None,
            },
            &initial_h, budget(&mesh, 2),
        ).unwrap();
        // Admission must preserve the latent initial state, even though the
        // prescribed bath has a different, uniquely invertible temperature.
        assert_eq!(run.specific_enthalpy_j_kg(), initial_h.as_slice());
        assert_eq!(run.temperature(), &[310.; 4]);
        assert_eq!(run.liquid_mass_fraction(), &[0.5; 4]);
        assert_eq!(run.time_s(), 0.);
        assert_eq!(run.accepted_steps(), 0);

        // This reference tetrahedron has volume 1/6 and mass 2/6.
        // The original chart gives h(330 K) = 400 J/kg. Generation is 30/6 W.
        let expected_storage_j = (2. / 6.) * (400. - 150.);
        let first = run.advance_to(cx, 0.25, config()).unwrap();
        assert!(first.specific_enthalpy_j_kg.iter().all(|h| (h - 400.).abs() < 1e-9));
        assert!(first.temperature.iter().all(|t| (t - 330.).abs() < 1e-10));
        assert_eq!(first.liquid_mass_fraction, vec![1.; 4]);
        assert!((first.source_w - 5.).abs() < 1e-12);
        assert!((first.stored_energy_change_j - expected_storage_j).abs() < 1e-9);
        assert!((first.dirichlet_in_w - (expected_storage_j / 0.25 - 5.)).abs() < 1e-9);
        assert!(first.energy_residual_j.abs() <= config().energy_tolerance_j);
        assert_eq!(run.specific_enthalpy_j_kg(), first.specific_enthalpy_j_kg.as_slice());
        assert_eq!(run.temperature(), first.temperature.as_slice());
        assert_eq!(run.liquid_mass_fraction(), first.liquid_mass_fraction.as_slice());

        // A later, longer interval has no repeated storage jump. The bath
        // removes all ongoing source heat while every constrained state stays put.
        let second = run.advance_to(cx, 0.75, config()).unwrap();
        assert_eq!(second.specific_enthalpy_j_kg, first.specific_enthalpy_j_kg);
        assert_eq!(second.temperature, first.temperature);
        assert_eq!(second.liquid_mass_fraction, first.liquid_mass_fraction);
        assert!(second.stored_energy_change_j.abs() < 1e-12);
        assert!((second.source_w - 5.).abs() < 1e-12);
        assert!((second.dirichlet_in_w + 5.).abs() < 1e-9);
        assert!(second.energy_residual_j.abs() <= config().energy_tolerance_j);
        assert!((first.stored_energy_change_j + second.stored_energy_change_j
            - expected_storage_j).abs() < 1e-9);
        assert_eq!(run.time_s(), 0.75);
        assert_eq!(run.accepted_steps(), 2);
    });
}

#[test]
fn g1_prescribed_history_held_events_use_material_charts_and_rollback_ambiguous_bath() {
    let mesh = mesh(true);
    let charts = [curve(100., 47), curve(400., 48)];
    let table = [
        ReferenceEnthalpyMaterial { curve: &charts[0], reference_density_kg_m3: 1. },
        ReferenceEnthalpyMaterial { curve: &charts[1], reference_density_kg_m3: 2. },
    ];
    let material = ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        let row = BoundaryCondition {
            region: "bath".into(), physics: Physics::Thermal, kind: BcKind::Dirichlet,
            value: Some(BcValue::Signal(TimeSignal::Table {
                times: vec![0., 0.25, 0.5, 0.75],
                values: vec![300., 340., 290., 310.],
                dims: TEMPERATURE_DIMS, interp: Interp::Hold,
            })),
            compatibility: None, frame: 0,
        };
        let boundary = ThermalBoundaryHistory::new(
            cx,
            ThermalBoundaryBuilder::new(&mesh)
                .remainder("bath", ThermalBc::dirichlet(300.).unwrap()).unwrap(),
            &Environment::earth_lab(),
            vec![ThermalRegionHistory { row, reference_temperature: None }], 0.,
        ).unwrap();
        let mut run = EnthalpyHistoryRun::new(
            cx,
            problem(&boundary, &material,
                HeatSourceHistory::Fixed(ScalarField::Uniform(0.))),
            EnthalpyHistoryMaterials {
                table: &table, of_element: &[0, 1], phase_conductivity: None,
            },
            &[80.; 8], budget(&mesh, 4),
        ).unwrap();
        assert_eq!(run.specific_enthalpy_j_kg(), &[80.; 8]);
        assert!(run.temperature().iter().all(|t| (t - 308.).abs() < 1e-10));
        assert_eq!(run.liquid_mass_fraction(), &[0.; 8]);

        let mut previous_time = 0.;
        let mut previous_h = [80.; 2];
        let mut stored_total_j = 0.;
        let mut input_total_j = 0.;
        // Every call requests time 1 with a cap larger than the entire window.
        // Exact event endpoints still use the PRE-jump bath temperature.
        for (ordinal, (end, expected_h, temperature, fraction)) in [
            (0.25, [0., 0.], 300., 0.),
            (0.5, [500., 800.], 340., 1.),
            (0.75, [-100., -100.], 290., 0.),
        ].into_iter().enumerate() {
            let step = run.advance_event_step(cx, 1., 2., config()).unwrap();
            assert_eq!(step.start_time_s, previous_time);
            assert_eq!(step.end_time_s, end);
            for body in 0..2 {
                for vertex in body * 4..body * 4 + 4 {
                    assert!((step.report.specific_enthalpy_j_kg[vertex] - expected_h[body]).abs() < 1e-9);
                    assert!((step.report.temperature[vertex] - temperature).abs() < 1e-10);
                    assert_eq!(step.report.liquid_mass_fraction[vertex], fraction);
                }
            }
            // Each tetrahedron has volume 1/6; its own frozen density and
            // original chart determine storage, including unequal latent heat.
            let expected_storage_j = (expected_h[0] - previous_h[0]) / 6.
                + (expected_h[1] - previous_h[1]) / 3.;
            let dt = end - previous_time;
            let net_w = step.report.source_w + step.report.dirichlet_in_w
                - step.report.neumann_out_w - step.report.robin_out_w;
            assert!(step.report.source_w.abs() < 1e-12);
            assert!((step.report.stored_energy_change_j - expected_storage_j).abs() < 1e-9);
            assert!((step.report.dirichlet_in_w - expected_storage_j / dt).abs() < 1e-8);
            assert!((expected_storage_j - dt * net_w).abs() < 1e-8);
            assert!(step.report.energy_residual_j.abs() <= config().energy_tolerance_j);
            assert_eq!(run.specific_enthalpy_j_kg(), step.report.specific_enthalpy_j_kg.as_slice());
            assert_eq!(run.time_s(), end);
            assert_eq!(run.accepted_steps(), ordinal + 1);
            stored_total_j += step.report.stored_energy_change_j;
            input_total_j += dt * net_w;
            previous_time = end;
            previous_h = expected_h;
        }
        assert!((stored_total_j + 90.).abs() < 1e-8);
        assert!((input_total_j + 90.).abs() < 1e-8);

        // The 310 K plateau takes effect only AFTER the 0.75 event. It has
        // multiple admissible h values, so the next endpoint must not publish.
        let accepted_h = run.specific_enthalpy_j_kg().to_vec();
        let accepted_t = run.temperature().to_vec();
        let accepted_f = run.liquid_mass_fraction().to_vec();
        assert!(matches!(
            run.advance_event_step(cx, 1., 2., config()),
            Err(HeterogeneousEnthalpyError::Enthalpy(
                EnthalpyError::AmbiguousPrescribedTemperature { .. }
            ))
        ));
        assert_eq!(run.specific_enthalpy_j_kg(), accepted_h.as_slice());
        assert_eq!(run.temperature(), accepted_t.as_slice());
        assert_eq!(run.liquid_mass_fraction(), accepted_f.as_slice());
        assert_eq!(run.time_s(), 0.75);
        assert_eq!(run.accepted_steps(), 3);
    });
}
