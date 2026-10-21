//! G1/G3: finite thermal stores power the original reciprocal field actuator.
//! These source-authored coefficients are manufactured, not measured alloys.
use super::*;
#[path = "../electrical_contacts/source.rs"]
mod field_source;
use fs_couple::circuit_enthalpy::thermoelectric::motion;
use fs_couple::material_circuit::motion::MaterialFieldMotion;
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_ga::Vec3;
use fs_material::field_storage::{MAGNETIC_PERMEABILITY_DIMS,
    MAGNETIC_PERMEABILITY_PROPERTY, ResolvedInductor};
use fs_mbd_material::{RegionBinding, VolumeMoments, resolve_body};

fn actuator(material: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'_>,
    linkage: f64, velocity: f64, limits: [f64; 2]) -> MaterialFieldMotion {
    let field = field_source::scalar_card(MAGNETIC_PERMEABILITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS, 1.0);
    let storage = MaterialStorage { branch: 4,
        element: MaterialStorageElement::Inductor(ResolvedInductor::resolve(&field,
            &QueryPoint::new(), MaterialPropertySelection::SingleClaimOnly, 1.0, 1.0, 1).unwrap()) };
    let circuit = MaterialCircuit::new_with_storage(CircuitGraph { node_count: 5, branches: vec![
        (1, 2, Branch::VoltageSource { port: 0 }), (2, 3, Branch::Resistor { ohms: f64::NAN }),
        (3, 4, Branch::VoltageSource { port: 1 }), (4, 0, Branch::Resistor { ohms: f64::NAN }),
        (1, 0, Branch::Inductor { henries: f64::NAN }),
    ], transformers: vec![] }, material.resolve(thermal, || false).unwrap().resistors,
        vec![storage], vec![linkage], vec![], &[0.0; 2]).unwrap();
    let density = field_source::scalar_card("density", Density::DIMS, 2.0);
    let body = resolve_body(&[RegionBinding { key: "guided pole", material: &density,
        point: &QueryPoint::new(), selection: MaterialPropertySelection::SingleClaimOnly,
        geometry: VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap() }], 1, || false).unwrap();
    MaterialFieldMotion::new(circuit, 4, body, velocity, limits).unwrap()
}

fn stored(motion: &MaterialFieldMotion) -> f64 {
    // Independent uniform-gap energy: mu=N=A=1, mass=2 kg in this fixture.
    let linkage = motion.circuit().fluxes_wb()[0];
    0.5*linkage*linkage*motion.gap_m() + motion.velocity_m_s().powi(2)
}

#[test]
fn g1_finite_thermal_gradient_powers_motion_and_affine_kelvin_heat_closes() {
    let run = |seebeck: f64| {
        let (mut cards, mut curves) = fixture(seebeck, -0.01, 3.0, 1.0);
        let slope = 1e-4;
        (cards[2], curves[2]) = source_value(PropertyValue::Curve {
            abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, seebeck-150.0*slope), (2000.0, seebeck+1650.0*slope)],
            dims: SEEBECK_DIMS,
        }, InterpolationPolicy::LinearInside, 2.0, 1.0);
        let context = point(350.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key);
        let resistors = resistors(); let links = links();
        let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
        let mut thermal = initial(&curves, [400.0, 300.0, 350.0, 350.0]);
        let mut moving = actuator(&material, &thermal, 0.0, 0.0, [0.5, 1.5]);
        let initial_energy = thermal.energy_j();
        let mut currents = Vec::new();
        for index in 1..=8 {
            let target = f64::from(index)*0.1;
            let dt = target-thermal.time_s();
            let next = motion::step_candidate(&material, &thermal, target, &moving,
                &[0.0; 2], 0.0, &legs(), &BASES, config(), &[0.0; 4], || false).unwrap();
            let report = next.report();
            let sample = &report.coupled.source.legs[0];
            let current = sample.transfer.current_a;
            let a = next.thermal().temperature_k(CellId(0)).unwrap();
            let b = next.thermal().temperature_k(CellId(1)).unwrap();
            let middle = next.thermal().temperature_k(CellId(2)).unwrap();
            // Independent affine-S integrals on the actual two thermal edges.
            let left = seebeck+slope*(0.5*(a+middle)-350.0);
            let right = seebeck+slope*(0.5*(middle+b)-350.0);
            let voltage = (seebeck+slope*(0.5*(a+b)-350.0))*(b-a);
            near(sample.transfer.voltage_drop_v, voltage);
            near(sample.transfer.terminal_heat_w[0], -left*a*current);
            near(sample.transfer.terminal_heat_w[1], right*b*current);
            near(sample.transfer.bulk_heat_w, middle*(left-right)*current);
            near(sample.transfer.terminal_heat_w[0]/a + sample.transfer.bulk_heat_w/middle
                + sample.transfer.terminal_heat_w[1]/b, 0.0);
            near(report.coupled.source.legs[1].transfer.current_a, current);
            near(report.coupled.deposited_joule_j, dt*current*current*(2.0+3.0));
            // KVL with the actual linkage derivative, not an authored RL current.
            let second_voltage = -0.01*(a-b);
            near((next.circuit().fluxes_wb()[0]-moving.circuit().fluxes_wb()[0])/dt,
                5.0*current+voltage+second_voltage);
            near(report.coupled.source.reversible_electrical_work_j,
                -dt*current*(voltage+second_voltage));
            near(report.coupled.source.reversible_heat_j,
                -report.coupled.source.reversible_electrical_work_j);
            near(report.external_electrical_work_j, 0.0);
            near(report.external_mechanical_work_j, 0.0);
            near(report.external_heat_j, 0.0);
            near(report.total_energy_residual_j, 0.0);
            assert!((next.thermal().energy_j()+stored(next.motion())-initial_energy).abs()<3e-5);
            assert!(report.coupled.electrical.field_work_on_body_j > 0.0);
            assert_eq!(sample.material.material()[0].card_identity(), cards[2].content_hash());
            currents.push(current);
            let (next_thermal, next_motion, _) = next.into_parts();
            thermal = next_thermal; moving = next_motion;
        }
        assert!(moving.gap_m()<1.0 && moving.velocity_m_s()<0.0);
        assert!(stored(&moving)>1e-3 && thermal.energy_j()<initial_energy);
        assert!(currents.iter().all(|i| *i>0.0));
        (moving.gap_m(), moving.circuit().fluxes_wb()[0], thermal.temperature_k(CellId(0)).unwrap())
    };
    let original = run(0.02);
    let changed = run(0.04);
    assert!(changed.0<original.0 && changed.1<original.1);
    assert!((changed.2-original.2).abs()>1e-4);
}

#[test]
fn g3_zero_seebeck_is_the_exact_existing_electrothermal_motion_path() {
    let (cards, curves) = fixture(0.0, 0.0, 3.0, 1.0);
    let context = point(350.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let resistors = resistors(); let links = links();
    let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let thermal = initial(&curves, [400.0, 300.0, 350.0, 350.0]);
    let moving = actuator(&material, &thermal, 0.8, 0.02, [0.5, 1.5]);
    let original = material.field_motion_candidate(&thermal, 0.1, &moving, &[0.0; 2], 0.2,
        &BASES, config(), &[0.0; 4], || false).unwrap();
    let joined = motion::step_candidate(&material, &thermal, 0.1, &moving, &[0.0; 2], 0.2,
        &legs(), &BASES, config(), &[0.0; 4], || false).unwrap();
    assert_eq!(joined.thermal().cells(), original.thermal().cells());
    assert_eq!(joined.motion(), original.motion());
    near(joined.report().coupled.source.reversible_heat_j, 0.0);
    near(joined.report().coupled.source.reversible_electrical_work_j, 0.0);
    near(joined.report().external_mechanical_work_j, 0.2*(joined.motion().gap_m()-moving.gap_m()));
    near(joined.thermal().energy_j()-thermal.energy_j()+stored(joined.motion())-stored(&moving),
        joined.report().external_mechanical_work_j);
}

#[test]
fn g3_seebeck_polarity_reverses_current_but_preserves_heat_and_attraction() {
    let run = |sign| {
        let (cards, curves) = fixture(sign*0.02, -sign*0.01, 3.0, 1.0);
        let context = point(350.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key);
        let resistors = resistors(); let links = links();
        let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
        let thermal = initial(&curves, [400.0, 300.0, 350.0, 350.0]);
        let moving = actuator(&material, &thermal, 0.0, 0.0, [0.5, 1.5]);
        let next = motion::step_candidate(&material, &thermal, 0.2, &moving, &[0.0; 2], 0.0,
            &legs(), &BASES, config(), &[0.0; 4], || false).unwrap();
        let sample = &next.report().coupled.source.legs[0];
        (sample.transfer.current_a, next.circuit().fluxes_wb()[0], next.motion().gap_m(),
            next.motion().velocity_m_s(), sample.transfer.terminal_heat_w,
            [0,1,2,3].map(|i| next.thermal().temperature_k(CellId(i)).unwrap()))
    };
    let a = run(1.0); let b = run(-1.0);
    assert!(a.0>0.0 && b.0<0.0);
    near(a.0, -b.0); near(a.1, -b.1); near(a.2, b.2); near(a.3, b.3);
    for (a,b) in a.4.into_iter().chain(a.5).zip(b.4.into_iter().chain(b.5)) { near(a,b); }
}

#[test]
fn g4_last_checkpoint_budget_and_gap_refusal_preserve_all_original_owners() {
    let (cards, curves) = fixture(0.02, -0.01, 3.0, 1.0);
    let context = point(350.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let resistors = resistors(); let links = links();
    let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let thermal = initial(&curves, [400.0, 300.0, 350.0, 350.0]);
    let moving = actuator(&material, &thermal, 0.4, 0.0, [0.5, 1.5]);
    let original_cells = thermal.cells().clone(); let original_motion = moving.clone();
    let mut count = 0;
    let first = motion::step_candidate(&material, &thermal, 0.1, &moving, &[0.0; 2], 0.0,
        &legs(), &BASES, config(), &[0.0; 4], || { count+=1; false }).unwrap();
    let mut actual = 0;
    assert!(matches!(motion::step_candidate(&material, &thermal, 0.1, &moving, &[0.0; 2], 0.0,
        &legs(), &BASES, config(), &[0.0; 4], || { actual+=1; actual==count }),
        Err(CircuitEnthalpyError::Cancelled)));
    let mut small = config(); small.maximum_circuit_trials = 1;
    assert!(motion::step_candidate(&material, &thermal, 0.1, &moving, &[0.0; 2], 0.0,
        &legs(), &BASES, small, &[0.0; 4], || false).is_err());
    let restricted = actuator(&material, &thermal, 0.4, -0.1, [0.999, 1.001]);
    let saved = restricted.clone();
    assert!(motion::step_candidate(&material, &thermal, 0.1, &restricted, &[0.0; 2], 0.0,
        &legs(), &BASES, config(), &[0.0; 4], || false).is_err());
    assert!(motion::step_candidate(&material, &thermal, 0.1, &moving, &[1.0, 0.0], 0.0,
        &legs(), &BASES, config(), &[0.0; 4], || false).is_err());
    assert_eq!(thermal.cells(), &original_cells); assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(moving, original_motion); assert_eq!(restricted, saved);
    let replay = motion::step_candidate(&material, &thermal, 0.1, &moving, &[0.0; 2], 0.0,
        &legs(), &BASES, config(), &[0.0; 4], || false).unwrap();
    assert_eq!(first.thermal().cells(), replay.thermal().cells());
    assert_eq!(first.motion(), replay.motion());
}
