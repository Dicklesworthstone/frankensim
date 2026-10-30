//! Independent coupled field/contact/conduction/radiosity references.
use super::*;
#[path = "../../../fs-conduction/tests/radiation_material/common.rs"]
mod radiation_geometry;
use fs_conduction::transient::deposition::{contact::*, radiation::*};
use fs_conduction::{SurfaceEmissivity, AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP};
use fs_couple::circuit_enthalpy::contacts::heat::{MaterialHeatNetwork, PreparedHeatNetwork};
use fs_exec::{CancelGate, Cx};

fn with_heat<R>(magnetic: bool, alpha: f64, emissivity: f64, resistance: f64, paths: bool,
    run: impl FnOnce(&Cx<'_>, &CancelGate, PreparedHeatNetwork<'_, '_>,
        &PreparedElectricalContactNetwork<'_, '_>, DepositionSimulation<'_>,
        MaterialCircuit, FieldThermalBinding) -> R,
) -> R {
    fixture_with_emissivity(magnetic, alpha, Some(2.0), 0.3, 2.0, emissivity,
        |electrical, material, thermal, circuit, field, network| {
            let cells = material.cells;
            let mut claims = ClaimSet::new();
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(RP, RD),
                value: PropertyValue::Scalar { value: resistance, dims: RD },
                validity: ValidityDomain::unconstrained().with("Ta", 250.0, 500.0)
                    .with("Tb", 250.0, 500.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
                provenance: Provenance { source: "manufactured thermal contact".into(),
                    license: "test-only".into(), artifact: None },
            }).unwrap();
            let interface = InterfaceSystemCard::assemble(
                SurfaceSpec { material: cells[0].point.material.id().clone(), texture_frame: "A".into() },
                SurfaceSpec { material: cells[1].point.material.id().clone(), texture_frame: "B".into() },
                SystemContext { medium: "test".into(), third_body: None,
                    environment: "test".into(), history: "fixed".into() }, claims, vec![],
            ).unwrap();
            let context = QueryPoint::new().with("Ta", 300.0).unwrap().with("Tb", 300.0).unwrap();
            let contacts = [MaterialContactLink {
                link: 3,
                point: ThermalInterfacePoint { interface: &interface,
                    materials: [cells[0].point.material, cells[1].point.material], context: &context,
                    temperature_axes: ["Ta", "Tb"], liquid_fraction_axes: [None, None],
                    mean_temperature_axis: None },
                sides: [0, 1].map(|i| ContactCell {
                    cell: cells[i].cell, material_id: cells[i].material_id, curve: cells[i].curve,
                }),
                selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            }];
            let surfaces = [0, 1].map(|i| RadiatingCell {
                surface: i, cell: cells[i].cell, material_id: cells[i].material_id,
                curve: cells[i].curve, material: cells[i].point,
                selection: MaterialPropertySelection::SingleClaimOnly,
            });
            let emissivities = [0, 1].map(|i| SurfaceEmissivity::from_thermal_state(
                "test", cells[i].point,
                cells[i].curve.state_at_specific_enthalpy(thermal.cells()[&cells[i].cell].specific_h_j_kg).unwrap(),
                MaterialPropertySelection::SingleClaimOnly,
            ).unwrap().0);
            let enclosure = radiation_geometry::enclosure(emissivities);
            let radiation = MaterialRadiation { enclosure: &enclosure, cells: &surfaces,
                power_tolerance_w: 1e-8, radiosity_tolerance_w_m2: 1e-8 };
            let assembly = MaterialHeatNetwork { electrical: network,
                thermal_contacts: if paths { &contacts } else { &[] },
                radiation: paths.then_some(&radiation) };
            let gate = CancelGate::new_clock_free();
            let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
            pool.scope(|arena| {
                let cx = Cx::new(&gate, arena,
                    fs_exec::StreamKey { seed: 73, kernel_id: 51, tile: 0, iteration: 0 },
                    fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
                let prepared = assembly.prepare(&cx, &thermal).unwrap();
                run(&cx, &gate, prepared, electrical, thermal, circuit, field)
            })
        })
}

// Independent SI equations for two unit-area gray parallel plates. The native
// enclosure uses explicit limiting view factors; no finite-cube claim is made.
const SIGMA: f64 = 5.670374419e-8;
fn radiative_power(ta: f64, tb: f64, epsilon: f64) -> f64 {
    SIGMA * (ta - tb) * (ta + tb) * (ta * ta + tb * tb)
        / (1.0 / epsilon + 1.0 / 0.6 - 1.0)
}
fn bisect(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) > 0.0 { hi = mid; } else { lo = mid; }
    }
    0.5 * (lo + hi)
}
fn heat_reference(magnetic: bool, alpha: f64, epsilon: f64, rt: f64,
    old: [f64; 3], dt: f64,
) -> [f64; 3] {
    let g = 1.0 / (1.0 / 2.0 + rt);
    let state = |ta: f64| {
        let k = 1.0 + alpha * (ta - 300.0);
        let rc = 2.0 * (1.0 + BETA * (ta - 300.0));
        let r = 1.0 + rc;
        let ratio = if magnetic { dt * r / (2.0 * k) } else { dt / (2.0 * k * r) };
        let z = old[2] * (1.0 - ratio) / (1.0 + ratio);
        let i = 0.5 * (old[2] + z) / k / if magnetic { 1.0 } else { r };
        let contact_heat = dt * rc * i * i;
        let tb = bisect(250.0, 500.0, |tb|
            CP[1] * (tb - old[1]) - 0.7 * contact_heat
                - dt * (g * (ta - tb) + radiative_power(ta, tb, epsilon)));
        [ta, tb, z]
    };
    let ta = bisect(250.0, 450.0, |ta| {
        let y = state(ta);
        CP[0] * (ta - old[0]) + CP[1] * (y[1] - old[1])
            + field_u(alpha, ta, y[2]) - field_u(alpha, old[0], old[2])
    });
    state(ta)
}

#[test]
fn complete_heat_transfer_and_both_field_kinds_match_independent_coupled_equations() {
    for magnetic in [false, true] {
        for alpha in [-0.001, 0.001] {
            for (epsilon, rt) in [(0.2, 0.0), (0.8, 3.0)] {
                with_heat(magnetic, alpha, epsilon, rt, true,
                    |cx, _, prepared, _, thermal, circuit, field| {
                        let expected = heat_reference(magnetic, alpha, epsilon, rt,
                            [300.0, 300.0, 10.0], DT);
                        let next = prepared.field_storage_candidate(cx, &thermal, DT,
                            &circuit, &[], &[field], controls(), &[0.0; 2]).unwrap();
                        near(next.thermal().temperature_k(A).unwrap(), expected[0]);
                        near(next.thermal().temperature_k(B).unwrap(), expected[1]);
                        near(coordinate(next.circuit(), magnetic), expected[2]);
                        let report = next.report();
                        let heat = &report.coupled.source.network;
                        near(heat.thermal_contacts[0].conductance_w_k, 1.0 / (0.5 + rt));
                        let radiation = heat.radiation.as_ref().unwrap();
                        near(radiation.sources_w()[&A], -radiative_power(expected[0], expected[1], epsilon));
                        near(radiation.sources_w()[&B], radiative_power(expected[0], expected[1], epsilon));
                        near(next.thermal().energy_j() + report.coupled.source.field_internal_energy_j[1],
                            thermal.energy_j() + field_u(alpha, 300.0, 10.0));
                        near(report.external_heat_j, 0.0);
                        near(report.total_energy_residual_j, 0.0);
                        assert_eq!(report.coupled.electrical.contacts.len(), 1);
                        assert!(report.coupled.source.reversible_heat_j.abs() > 0.01);
                    });
            }
        }
    }
}

#[test]
fn surface_and_interface_data_change_current_and_finite_field_history() {
    let mut endpoints = Vec::new();
    for (epsilon, rt) in [(0.2, 0.3), (0.8, 0.3), (0.2, 3.0)] {
        with_heat(false, 0.001, epsilon, rt, true,
            |cx, _, prepared, _, mut thermal, mut circuit, field| {
                let fields = [field];
                let mut expected = [300.0, 300.0, 10.0];
                let total = thermal.energy_j() + field_u(0.001, 300.0, 10.0);
                for n in 1..=10 {
                    let target = n as f64 * 0.1;
                    expected = heat_reference(false, 0.001, epsilon, rt, expected, target - thermal.time_s());
                    let next = prepared.field_storage_candidate(cx, &thermal, target,
                        &circuit, &[], &fields, controls(), &[0.0; 2]).unwrap();
                    near(next.thermal().temperature_k(A).unwrap(), expected[0]);
                    near(next.thermal().temperature_k(B).unwrap(), expected[1]);
                    near(next.circuit().charges_c()[0], expected[2]);
                    near(next.thermal().energy_j() + next.report().coupled.source.field_internal_energy_j[1], total);
                    let (t, c, _) = next.into_parts(); thermal = t; circuit = c;
                }
                endpoints.push(expected);
            });
    }
    for changed in &endpoints[1..] {
        for (a, b) in endpoints[0].iter().zip(changed) { assert!((a - b).abs() > 1e-6); }
    }
}

#[test]
fn no_heat_paths_preserve_contact_field_arithmetic_and_failed_steps_preserve_history() {
    with_heat(false, 0.001, 0.3, 0.2, false,
        |cx, gate, prepared, electrical, mut thermal, mut circuit, field| {
            let fields = [field];
            let first = prepared.field_storage_candidate(cx, &thermal, 0.1,
                &circuit, &[], &fields, controls(), &[0.0; 2]).unwrap();
            let original = electrical.field_storage_candidate(&thermal, 0.1,
                &circuit, &[], &fields, controls(), &[0.0; 2], || false).unwrap();
            assert_eq!(first.circuit(), original.circuit());
            assert_eq!(first.thermal().cells(), original.thermal().cells());
            assert_eq!(first.report().coupled.source.fields, original.report().coupled.source.fields);
            let (t, c, _) = first.into_parts(); thermal = t; circuit = c;
            let saved = thermal.cells().clone(); let saved_circuit = circuit.clone();
            let mut exhausted = controls(); exhausted.maximum_circuit_trials = 1;
            assert!(matches!(prepared.field_storage_candidate(cx, &thermal, 0.2,
                &circuit, &[], &fields, exhausted, &[0.0; 2]), Err(CircuitEnthalpyError::CircuitBudget)));
            let retry = prepared.field_storage_candidate(cx, &thermal, 0.2,
                &circuit, &[], &fields, controls(), &[0.0; 2]).unwrap();
            let replay = prepared.field_storage_candidate(cx, &thermal, 0.2,
                &circuit, &[], &fields, controls(), &[0.0; 2]).unwrap();
            assert_eq!(retry.circuit(), replay.circuit());
            assert_eq!(retry.thermal().cells(), replay.thermal().cells());
            gate.request();
            assert!(matches!(prepared.field_storage_candidate(cx, &thermal, 0.2,
                &circuit, &[], &fields, controls(), &[0.0; 2]), Err(CircuitEnthalpyError::Cancelled)));
            assert_eq!(thermal.cells(), &saved); assert_eq!(circuit, saved_circuit);
            assert_eq!(thermal.time_s(), 0.1);
        });
}
