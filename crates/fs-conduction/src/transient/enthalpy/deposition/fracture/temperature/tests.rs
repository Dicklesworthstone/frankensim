use super::*;
use crate::transient::deposition::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};
use crate::transient::enthalpy::deposition::fracture::ContactPolicy;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::source::*;
use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
use fs_qty::{Dimensionless, Pressure, Temperature};

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=600_000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0 + h / 2000.0)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.0005
    }
}
static HEAT: Sensible = Sensible;

fn card(label: &str, mixed: bool, viscous: bool) -> InterfaceSystemCard {
    let stiffness = interface_stiffness_per_area_dims();
    let mut properties = if mixed {
        vec![
            (COHESIVE_NORMAL_STIFFNESS_PROPERTY, stiffness, [1e9, 4e8]),
            (
                COHESIVE_TANGENTIAL_STIFFNESS_PROPERTY,
                stiffness,
                [5e8, 2e8],
            ),
            (
                COHESIVE_NORMAL_STRENGTH_PROPERTY,
                Pressure::DIMS,
                [2e5, 2e4],
            ),
            (
                COHESIVE_TANGENTIAL_STRENGTH_PROPERTY,
                Pressure::DIMS,
                [3e5, 3e4],
            ),
            (
                COHESIVE_MODE_I_FRACTURE_ENERGY_PROPERTY,
                FRACTURE_ENERGY_PER_AREA_DIMS,
                [100.0, 1.0],
            ),
            (
                COHESIVE_MODE_II_FRACTURE_ENERGY_PROPERTY,
                FRACTURE_ENERGY_PER_AREA_DIMS,
                [240.0, 2.4],
            ),
            (COHESIVE_BK_EXPONENT_PROPERTY, Dimensionless::DIMS, [2.0; 2]),
        ]
    } else {
        vec![
            (COHESIVE_STIFFNESS_PROPERTY, stiffness, [1e9, 4e8]),
            (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, [2e5, 2e4]),
            (
                COHESIVE_FRACTURE_ENERGY_PROPERTY,
                FRACTURE_ENERGY_PER_AREA_DIMS,
                [100.0, 1.0],
            ),
        ]
    };
    properties.push((COHESIVE_COMPRESSION_STIFFNESS_PROPERTY, stiffness, [2e9; 2]));
    if viscous {
        properties.extend([
            (
                COHESIVE_NORMAL_VISCOSITY_PROPERTY,
                COHESIVE_VISCOSITY_DIMS,
                [1e8, 2e7],
            ),
            (
                COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY,
                COHESIVE_VISCOSITY_DIMS,
                [2e8, 4e7],
            ),
        ]);
    }
    let mut claims = ClaimSet::new();
    for (name, dims, values) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (400.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: label.into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "manufactured fixture; not calibrated".into(),
            revision: 1,
        },
        texture_frame: "fixed patch axes".into(),
    };
    InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "declared bonded interface".into(),
            third_body: None,
            environment: "synthetic thermal-adhesion test".into(),
            history: "initial bond".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn source(card: &InterfaceSystemCard, mixed: bool, viscous: bool) -> ResolvedBilinearCohesive {
    let query = QueryPoint::new().with("T", 300.0).unwrap();
    let select = MaterialPropertySelection::SingleClaimOnly;
    if viscous {
        resolve_viscous_bilinear_cohesive(
            card,
            &query,
            select,
            if mixed {
                CohesiveEnvelope::FrozenOnsetMixedMode
            } else {
                CohesiveEnvelope::CommonMode
            },
            0.01,
            0.1,
        )
        .unwrap()
    } else if mixed {
        resolve_mixed_mode_bilinear_cohesive(card, &query, select, 0.01).unwrap()
    } else {
        resolve_bilinear_cohesive(card, &query, select, 0.01).unwrap()
    }
}
fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        source_w,
        reservoir: None,
    }
}
fn thermal(temperatures: [f64; 2], heat_cycle: bool) -> DepositionSimulation<'static> {
    let mut schedule = Vec::new();
    for (index, temperature) in temperatures.into_iter().enumerate() {
        let id = index as u64 + 1;
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: id,
                material_id: MaterialId(1),
                mass_kg: 0.001,
                specific_h_j_kg: 2000.0 * (temperature - 300.0),
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, id as f64 * 0.001],
                    end_m: [0.01, 0.0, id as f64 * 0.001],
                    volume_m3: 1e-6,
                    surface_area_m2: 0.001,
                },
                boundary: boundary(if heat_cycle { 2.0 } else { 0.0 }),
            }),
        });
    }
    schedule.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(link()),
    });
    if heat_cycle {
        for id in [CellId(1), CellId(2)] {
            schedule.push(TimedEvent {
                time_s: 100.0,
                event: Event::Boundary(id, boundary(-2.0)),
            });
        }
    }
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), &HEAT)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 5,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}
fn link() -> ThermalLink {
    ThermalLink {
        id: 7,
        a: CellId(1),
        b: CellId(2),
        area_m2: 0.0002,
        conductance_w_k: 0.2,
    }
}
fn thermal_config() -> StepConfig {
    StepConfig {
        max_dt_s: 1.0,
        energy_tolerance_j: 1e-8,
    }
}
fn config() -> TemperatureConfig {
    TemperatureConfig {
        local_energy_tolerance_j_m2: 1e-8,
        total_energy_tolerance_j: 1e-10,
    }
}
fn state(
    source: &ResolvedBilinearCohesive,
    thermal: DepositionSimulation<'static>,
    gap: [f64; 3],
    open_gap: f64,
) -> FractureCooling<'static, '_> {
    let point = source
        .law()
        .bonded_point()
        .trial_over_time(gap, 1.0, 1e-8)
        .unwrap()
        .equilibrium
        .point;
    FractureCooling::new(
        thermal,
        vec![Interface {
            reference: link(),
            point,
            closed: true,
            policy: ContactPolicy {
                open_gap_m: open_gap,
                reclose_gap_m: 0.0,
                reclosed_conductance_w_k: 0.05,
                maximum_slip_m: 0.005,
            },
        }],
        1,
        || Ok(()),
    )
    .unwrap()
}
fn update<'a, 'b>(
    previous: &'a ResolvedBilinearCohesive,
    successor: &'b ResolvedBilinearCohesive,
    weight: f64,
) -> [TemperatureUpdate<'a, 'b>; 1] {
    [TemperatureUpdate {
        interface: 7,
        previous,
        successor,
        temperature_weight_a: weight,
    }]
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-10, "{a:e} != {b:e}");
}

// G1/G3: actual heat weakens the retained bond, including mixed and viscous
// sources. A zero-duration coefficient change does no dashpot work.
#[test]
fn actual_heat_debonds_and_cooling_preserves_failure_and_energy_history() {
    for (mixed, viscous) in [(false, false), (true, false), (false, true), (true, true)] {
        let card = card("manufactured thermal cycle", mixed, viscous);
        let cold = source(&card, mixed, viscous);
        let initial = state(&cold, thermal([300.0; 2], true), [0.0002, 0.0, 0.0], 0.0001);
        let warmed = initial
            .advance_to(100.0, thermal_config(), || Ok(()))
            .unwrap()
            .0;
        let hot = cold
            .at_temperature(&card, "T", warmed.interface_temperature_k(7, 0.5).unwrap())
            .unwrap();
        assert_eq!(hot.resolved().state_coordinate("T"), Some(400.0));
        // A timed segment cannot be silently reinterpreted under new laws.
        assert!(
            warmed
                .transfer_temperature(100.0, &update(&cold, &hot, 0.5), "T", config(), || Ok(()))
                .is_err()
        );
        let point = warmed.interfaces()[&7].point;
        let held = warmed
            .synchronize_over_time(100.0, 100.0, &[(7, point)], 1e-8, || Ok(()))
            .unwrap()
            .0;
        let (failed, r) = held
            .transfer_temperature(100.0, &update(&cold, &hot, 0.5), "T", config(), || Ok(()))
            .unwrap();
        assert!(failed.interfaces()[&7].point.failed());
        assert!(!failed.interfaces()[&7].closed);
        assert!(!failed.thermal().links().contains_key(&7));
        assert_eq!(r.contact.opened, 1);
        assert_eq!(r.contact.viscous_dissipation_j, 0.0);
        close(r.parameter_work_j, -12.0 * link().area_m2);
        close(r.contact.fracture_dissipation_j, 8.0 * link().area_m2);
        close(r.stored_change_j, -20.0 * link().area_m2);
        close(r.defect_j, 0.0);
        assert_eq!(failed.thermal().cells(), held.thermal().cells());
        assert_eq!(failed.last_synchronized_time_s(), 100.0);
        assert!(!held.interfaces()[&7].point.failed());

        let cooled = failed
            .advance_to(200.0, thermal_config(), || Ok(()))
            .unwrap()
            .0;
        let point = cooled.interfaces()[&7].point;
        let cooled = cooled
            .synchronize_over_time(200.0, 100.0, &[(7, point)], 1e-8, || Ok(()))
            .unwrap()
            .0;
        let again = hot
            .at_temperature(&card, "T", cooled.interface_temperature_k(7, 0.5).unwrap())
            .unwrap();
        assert_eq!(again.resolved().state_coordinate("T"), Some(300.0));
        let (cooled, r) = cooled
            .transfer_temperature(200.0, &update(&hot, &again, 0.5), "T", config(), || Ok(()))
            .unwrap();
        assert!(cooled.interfaces()[&7].point.failed());
        assert!(!cooled.interfaces()[&7].closed);
        close(cooled.interfaces()[&7].point.dissipated_j_m2(), 8.0);
        close(r.contact.fracture_dissipation_j, 0.0);
        close(r.parameter_work_j, 0.0);
        assert_eq!(cooled.interfaces()[&7].reference, link());
    }
}

// G1: a still-closed failed interface retains its explicit lower conductance;
// an opened interface transfers no heat. Neither changes original area.
#[test]
fn failed_closed_patch_uses_declared_conductance_and_open_patch_stops_heat() {
    let card = card("manufactured failed patch", false, false);
    let cold = source(&card, false, false);
    let hot = cold.at_temperature(&card, "T", 400.0).unwrap();
    for (open_gap, expected_hot, expected_cold, opened) in
        [(0.0003, 397.5, 302.5, 0), (0.0001, 400.0, 300.0, 1)]
    {
        let original = state(
            &cold,
            thermal([400.0, 300.0], false),
            [0.0002, 0.0, 0.0],
            open_gap,
        );
        let (next, r) = original
            .transfer_temperature(0.0, &update(&cold, &hot, 1.0), "T", config(), || Ok(()))
            .unwrap();
        assert_eq!(r.contact.opened, opened);
        assert!(next.interfaces()[&7].point.failed());
        if opened == 0 {
            assert_eq!(r.contact.changed_conductances, 1);
            assert_eq!(next.thermal().links()[&7].conductance_w_k, 0.05);
            assert_eq!(next.thermal().links()[&7].area_m2, link().area_m2);
        }
        let cooled = next.advance_to(1.0, thermal_config(), || Ok(())).unwrap().0;
        // C1=C2=2 J/K; a closed patch transfers .05*(400-300)=5 J.
        close(
            cooled.thermal().temperature_k(CellId(1)).unwrap(),
            expected_hot,
        );
        close(
            cooled.thermal().temperature_k(CellId(2)).unwrap(),
            expected_cold,
        );
    }
}

#[test]
fn weighted_temperature_preserves_partial_mixed_mode_history() {
    let card = card("manufactured partial mixed bond", true, false);
    let cold = source(&card, true, false);
    let original = state(
        &cold,
        thermal([400.0, 300.0], false),
        [0.0004, 0.0002, -0.0001],
        0.0005,
    );
    assert!(original.interfaces()[&7].point.damage() > 0.0);
    let t = original.interface_temperature_k(7, 0.25).unwrap();
    assert_eq!(t, 325.0);
    let changed = cold.at_temperature(&card, "T", t).unwrap();
    let expected = original.interfaces()[&7]
        .point
        .transfer_temperature(changed.law(), 1e-8)
        .unwrap();
    let (next, r) = original
        .transfer_temperature(
            0.0,
            &update(&cold, &changed, 0.25),
            "T",
            config(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(next.interfaces()[&7].point, expected.point);
    close(
        r.parameter_work_j,
        link().area_m2 * expected.parameter_work_j_m2,
    );
    close(
        r.contact.fracture_dissipation_j,
        link().area_m2 * expected.damage_dissipated_j_m2,
    );
    assert_eq!(
        next.interfaces()[&7]
            .point
            .mixed_mode_history()
            .unwrap()
            .onset_direction,
        original.interfaces()[&7]
            .point
            .mixed_mode_history()
            .unwrap()
            .onset_direction
    );
    let (_, repeat) = next
        .transfer_temperature(
            0.0,
            &update(&changed, &changed, 0.25),
            "T",
            config(),
            || Ok(()),
        )
        .unwrap();
    close(repeat.contact.fracture_dissipation_j, 0.0);
    close(repeat.parameter_work_j, 0.0);
}

#[test]
fn source_temperature_clock_and_cancellation_refusals_preserve_the_input() {
    let card = card("original manufactured source", false, false);
    let cold = source(&card, false, false);
    let hot = cold.at_temperature(&card, "T", 400.0).unwrap();
    let original = state(
        &cold,
        thermal([400.0, 300.0], false),
        [0.0002, 0.0, 0.0],
        0.0001,
    );
    let other_card = self::card("other provenance, same coefficients", false, false);
    let other = source(&other_card, false, false)
        .at_temperature(&other_card, "T", 400.0)
        .unwrap();
    assert!(matches!(
        original.transfer_temperature(0.0, &update(&cold, &other, 1.0), "T", config(), || Ok(())),
        Err(Error::Source { .. })
    ));
    assert!(matches!(
        original.transfer_temperature(0.0, &update(&cold, &hot, 0.5), "T", config(), || Ok(())),
        Err(Error::Temperature {
            expected_k: 350.0,
            ..
        })
    ));
    assert!(
        original
            .transfer_temperature(0.0, &update(&hot, &hot, 1.0), "T", config(), || Ok(()))
            .is_err()
    );
    assert!(
        original
            .transfer_temperature(0.0, &[], "T", config(), || Ok(()))
            .is_err()
    );
    assert!(
        original
            .transfer_temperature(1.0, &update(&cold, &hot, 1.0), "T", config(), || Ok(()))
            .is_err()
    );
    for weight in [f64::NAN, -0.1, 1.1] {
        assert!(
            original
                .transfer_temperature(0.0, &update(&cold, &hot, weight), "T", config(), || Ok(()))
                .is_err()
        );
    }
    let mut polls = 0;
    original
        .transfer_temperature(0.0, &update(&cold, &hot, 1.0), "T", config(), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    for at in 1..=polls {
        let mut count = 0;
        let result =
            original.transfer_temperature(0.0, &update(&cold, &hot, 1.0), "T", config(), || {
                count += 1;
                if count == at {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            });
        assert!(matches!(
            result,
            Err(Error::Contact(super::super::Error::Thermal(
                DepositionError::Cancelled
            )))
        ));
        assert!(!original.interfaces()[&7].point.failed());
        assert_eq!(original.thermal().links()[&7], link());
        assert_eq!(original.thermal().temperature_k(CellId(1)).unwrap(), 400.0);
        assert_eq!(original.last_synchronized_time_s(), 0.0);
    }
}
