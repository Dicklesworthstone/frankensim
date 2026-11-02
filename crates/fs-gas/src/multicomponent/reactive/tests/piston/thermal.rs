use super::molecular::molecular;
use super::*;
use crate::multicomponent::dissipation::{
    ConstantTransport, DissipationControl, DissipationOperator,
};
use crate::multicomponent::wall_exchange::{
    MaterialWallCaloric, WallContact, WallControl, WallDrag, WallExchangeError, WallThermal,
};
use fs_thermal::enthalpy::{
    Boundary as ThermalBoundary, CellGeometry, CellId, DepositedCell, DepositionError, MaterialId,
    ThermalLaw, ThermalStorageBasis,
};

const ADIABATIC: WallContact = WallContact {
    drag: WallDrag::None,
    thermal: WallThermal::Adiabatic,
};

fn finite(temperature_k: f64, capacity: f64, conductance: f64) -> WallContact {
    WallContact {
        drag: WallDrag::None,
        thermal: WallThermal::Lumped {
            temperature_k,
            heat_capacity_j_per_k: capacity,
            conductance_w_per_k: conductance,
        },
    }
}

fn step(
    flow: &mut ReactingEuler1d,
    molecular: &mut DissipationOperator<ConstantTransport>,
    thermal: &mut PistonThermalEnds<'_>,
    piston: &mut Piston,
    cancelled: impl FnMut() -> bool,
) -> Result<crate::dynamics::piston::PistonStep<ThermalReactiveStepReport>, ReactiveError> {
    flow.step_with_piston_dissipation_and_heat(
        molecular,
        thermal,
        piston,
        Time::new(0.005),
        ReactiveStepControl::default(),
        DissipationControl::default(),
        WallControl::default(),
        load(),
        piston_controls(),
        cancelled,
    )
}

fn physical_change(gas: &MixtureGas, before: Totals, after: Totals) -> f64 {
    sum(
        core::iter::once(after.energy - before.energy).chain((0..gas.species_count()).map(|k| {
            (after.species_masses[k] - before.species_masses[k]) * gas.energy_references()[k]
        })),
    )
}

fn wall_temperature(contact: WallContact) -> f64 {
    match contact.thermal {
        WallThermal::Lumped { temperature_k, .. } => temperature_k,
        _ => panic!("expected finite temperature storage"),
    }
}

#[test]
fn g1_finite_end_heat_changes_actual_pressure_and_piston_motion_with_closed_energy() {
    let mechanism = mechanism(0.0, 0.0, 0.0);
    let mut insulated = closed(&mechanism, 1, chemistry());
    let mut insulated_operator = molecular(insulated.flow());
    let mut insulated_piston = Piston::new(0.2, 0.0).unwrap();
    insulated
        .step_with_piston_and_dissipation(
            &mut insulated_operator,
            &mut insulated_piston,
            Time::new(0.005),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            load(),
            piston_controls(),
            || false,
        )
        .unwrap();
    close(insulated_piston.velocity_m_s(), 0.0, 1e-12);
    for initial_wall_t in [3.0, 0.4] {
        let mut flow = closed(&mechanism, 1, chemistry());
        let mut operator = molecular(flow.flow());
        let mut walls = PistonThermalEnds::new(
            flow.flow(),
            [ADIABATIC, finite(initial_wall_t, 0.4, 1.0)],
            [None, None],
            || false,
        )
        .unwrap();
        let mut piston = Piston::new(0.2, 0.0).unwrap();
        let before = flow.flow().totals();
        let mut exterior_work = 0.0;
        for index in 0..4 {
            let report = step(&mut flow, &mut operator, &mut walls, &mut piston, || false).unwrap();
            if index == 0 {
                // Independent two-body backward-Euler exchange. The fixture's
                // gas mass is 0.1 kg and cv=2.5 J/(kg K), so Cg=0.25 J/K.
                let hdt = 0.5 * report.gas.reactive.elapsed.value();
                let expected =
                    hdt * (initial_wall_t - 1.0) / (1.0 + hdt * (1.0 / 0.25 + 1.0 / 0.4));
                close(report.gas.heat.first.ends[1].heat_to_gas_j, expected, 2e-10);
            }
            assert_eq!(report.gas.heat.external_heat_to_gas_j, 0.0);
            for half in [report.gas.heat.first, report.gas.heat.second] {
                assert_eq!(
                    half.elapsed.value(),
                    0.5 * report.gas.reactive.elapsed.value()
                );
                for end in half.ends {
                    assert_eq!(end.wall_impulse_ns, 0.0);
                    assert_eq!(end.drag_heating_j, 0.0);
                }
            }
            close(report.gas.heat.energy_residual_j, 0.0, 4e-10);
            close(report.gas.reactive.absolute_energy_residual_j, 0.0, 4e-11);
            close(
                report.total_energy_residual_j,
                0.0,
                piston_controls().energy_tolerance_j,
            );
            close(
                report.mechanical_energy_residual_j,
                0.0,
                piston_controls().energy_tolerance_j,
            );
            close(
                report.combined_momentum_residual_n_s,
                0.0,
                piston_controls().momentum_tolerance_n_s,
            );
            assert_eq!(walls.last_report(), Some(report.gas.heat));
            assert_eq!(walls.mesh(), flow.flow().mesh());
            assert_eq!(operator.mesh(), flow.flow().mesh());
            exterior_work += report.exterior_pressure_work_j;
        }
        let gas = flow
            .flow()
            .gas()
            .primitive(flow.flow().cell(0).unwrap(), || false)
            .unwrap();
        let sign = (initial_wall_t - 1.0).signum();
        assert!(sign * piston.velocity_m_s() > 0.0);
        assert!(sign * (gas.pressure - 1.0) > 0.0);
        assert!(sign * (initial_wall_t - wall_temperature(walls.contacts()[1])) > 0.0);
        let wall_change = 0.4 * (wall_temperature(walls.contacts()[1]) - initial_wall_t);
        let kinetic = 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2);
        close(
            physical_change(flow.flow().gas(), before, flow.flow().totals())
                + wall_change
                + kinetic
                + exterior_work,
            0.0,
            4.0 * piston_controls().energy_tolerance_j,
        );
    }
}

// Manufactured original thermal source u(T)=offset+a*T²/2. It has a genuinely
// varying capacity, retained reference geometry and an arbitrary energy datum.
struct QuadraticStore {
    a: f64,
    offset: f64,
    basis: ThermalStorageBasis,
}
impl ThermalLaw for QuadraticStore {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.basis
    }
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        if !energy.is_finite() || energy <= self.offset {
            return Err(DepositionError::MaterialDomain);
        }
        Ok((2.0 * (energy - self.offset) / self.a).sqrt())
    }
    fn maximum_temperature_slope(&self) -> f64 {
        100.0
    }
    fn validate_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        if cell.mass_kg != 0.4 || cell.geometry.volume_m3 != 0.1 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }
}

fn caloric(offset: f64, basis: ThermalStorageBasis, lower: f64) -> MaterialWallCaloric<'static> {
    MaterialWallCaloric::from_thermal_law(
        QuadraticStore {
            a: 0.4,
            offset,
            basis,
        },
        DepositedCell {
            id: CellId(17),
            segment_id: 3,
            body_id: 8,
            material_id: MaterialId(29),
            mass_kg: 0.4,
            specific_h_j_kg: offset + 0.2 * 9.0,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [1.0, 0.0, 0.0],
                volume_m3: 0.1,
                surface_area_m2: 1.0,
            },
            boundary: ThermalBoundary {
                ambient_k: 1.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        },
        [offset + lower, offset + 20.0],
    )
    .unwrap()
}

fn source_fixture(
    reference: f64,
    wall_offset: f64,
) -> (
    ReactingEuler1d,
    DissipationOperator<ConstantTransport>,
    PistonThermalEnds<'static>,
    Piston,
) {
    let flow = closed(&mechanism(8.0, 30.0, reference), 3, chemistry());
    let operator = molecular(flow.flow());
    let source = caloric(wall_offset, ThermalStorageBasis::InternalEnergy, 0.01);
    let contact = source.initial_contact(0.6, WallDrag::None).unwrap();
    let walls = PistonThermalEnds::new(
        flow.flow(),
        [
            WallContact {
                drag: WallDrag::None,
                thermal: WallThermal::Reservoir {
                    temperature_k: 4.0,
                    conductance_w_per_k: 0.3,
                },
            },
            contact,
        ],
        [None, Some(source)],
        || false,
    )
    .unwrap();
    (flow, operator, walls, Piston::new(0.2, 0.0).unwrap())
}

fn wall_energy(contact: WallContact) -> (f64, f64) {
    match contact.thermal {
        WallThermal::Caloric {
            mass_kg,
            specific_energy_j_kg,
            ..
        } => (mass_kg, specific_energy_j_kg),
        _ => panic!("expected original energy coordinate"),
    }
}

#[test]
fn g1_original_wall_energy_gauge_and_reservoir_heat_preserve_chemical_mechanical_balance() {
    let (mut a, mut da, mut wa, mut pa) = source_fixture(0.0, 0.0);
    let (mut b, mut db, mut wb, mut pb) = source_fixture(1000.0, 2000.0);
    let before = a.flow().totals();
    let (wall_mass, wall_u0) = wall_energy(wa.contacts()[1]);
    let mut reservoir_heat = 0.0;
    let mut exterior_work = 0.0;
    for _ in 0..3 {
        let left = step(&mut a, &mut da, &mut wa, &mut pa, || false).unwrap();
        let right = step(&mut b, &mut db, &mut wb, &mut pb, || false).unwrap();
        let heat = left.gas.heat;
        close(
            heat.external_heat_to_gas_j,
            heat.first.ends[0].heat_to_gas_j + heat.second.ends[0].heat_to_gas_j,
            2e-14,
        );
        assert_eq!(heat.first.ends[1].external_heat_to_gas_j, 0.0);
        assert_eq!(heat.second.ends[1].external_heat_to_gas_j, 0.0);
        assert!(heat.external_heat_to_gas_j > 0.0);
        close(
            left.total_energy_residual_j,
            0.0,
            piston_controls().energy_tolerance_j,
        );
        close(
            left.gas.gas_energy_change_j,
            physical_change(a.flow().gas(), left.gas.before, left.gas.after),
            2e-14,
        );
        close(
            left.gas.reactive.elapsed.value(),
            right.gas.reactive.elapsed.value(),
            1e-15,
        );
        close(pa.velocity_m_s(), pb.velocity_m_s(), 3e-9);
        close(
            heat.wall_energy_change_j,
            right.gas.heat.wall_energy_change_j,
            3e-9,
        );
        for i in 0..a.flow().cell_count() {
            let qa = a
                .flow()
                .gas()
                .primitive(a.flow().cell(i).unwrap(), || false)
                .unwrap();
            let qb = b
                .flow()
                .gas()
                .primitive(b.flow().cell(i).unwrap(), || false)
                .unwrap();
            close(qa.pressure, qb.pressure, 3e-9);
            close(qa.temperature, qb.temperature, 3e-9);
        }
        reservoir_heat += heat.external_heat_to_gas_j;
        exterior_work += left.exterior_pressure_work_j;
    }
    // Independent original endpoint u(T), not a heat-capacity tangent or the
    // sum of reported heat sources. Chemical formation energy is counted once.
    let (_, wall_u1) = wall_energy(wa.contacts()[1]);
    let wall_temperature = wa.heat_storage()[1]
        .as_ref()
        .unwrap()
        .temperature_k(wa.contacts()[1].thermal)
        .unwrap();
    close(wall_u1, 0.2 * wall_temperature.powi(2), 2e-13);
    close(
        physical_change(a.flow().gas(), before, a.flow().totals())
            + wall_mass * (wall_u1 - wall_u0)
            + 0.5 * pa.mass_kg() * pa.velocity_m_s().powi(2)
            + exterior_work
            - reservoir_heat,
        0.0,
        3.0 * piston_controls().energy_tolerance_j,
    );
}

#[test]
fn g4_final_cancellation_preserves_moved_thermal_stores_geometry_and_repeatability() {
    let (mut a, mut da, mut wa, mut pa) = source_fixture(0.0, 0.0);
    let (mut b, mut db, mut wb, mut pb) = source_fixture(0.0, 0.0);
    step(&mut a, &mut da, &mut wa, &mut pa, || false).unwrap();
    step(&mut b, &mut db, &mut wb, &mut pb, || false).unwrap();
    let cells = a.accepted.states.clone();
    let mesh = a.flow().mesh().clone();
    let time = a.flow().time();
    let piston = pa;
    let walls = *wa.contacts();
    let heat = wa.last_report();
    let mut polls = 0;
    let expected = step(&mut b, &mut db, &mut wb, &mut pb, || {
        polls += 1;
        false
    })
    .unwrap();
    assert!(polls > 50);
    let mut count = 0;
    let error = step(&mut a, &mut da, &mut wa, &mut pa, || {
        count += 1;
        count == polls
    })
    .unwrap_err();
    assert_eq!(error, ReactiveError::Cancelled);
    assert_eq!(count, polls);
    assert_eq!(a.accepted.states, cells);
    assert_eq!(a.flow().mesh(), &mesh);
    assert_eq!(a.flow().time(), time);
    assert_eq!(da.mesh(), &mesh);
    assert_eq!(wa.mesh(), &mesh);
    assert_eq!(wa.contacts(), &walls);
    assert_eq!(wa.last_report(), heat);
    assert_eq!(pa, piston);
    let resumed = step(&mut a, &mut da, &mut wa, &mut pa, || false).unwrap();
    assert_eq!(resumed, expected);
    assert_eq!(a.accepted.states, b.accepted.states);
    assert_eq!(pa, pb);
    assert_eq!(wa.contacts(), wb.contacts());
    assert_eq!(wa.mesh(), wb.mesh());
    assert_eq!(da.mesh(), db.mesh());

    // An original storage interval that cannot surrender any heat must refuse
    // the entire later step, retaining accepted gas, mechanical and wall state.
    let source = caloric(0.0, ThermalStorageBasis::InternalEnergy, 0.2 * 9.0);
    let contact = source.initial_contact(0.6, WallDrag::None).unwrap();
    let mut limited =
        PistonThermalEnds::new(a.flow(), [ADIABATIC, contact], [None, Some(source)], || {
            false
        })
        .unwrap();
    let cells = a.accepted.states.clone();
    let mesh = a.flow().mesh().clone();
    let time = a.flow().time();
    let piston = pa;
    let refusal = a
        .step_with_piston_dissipation_and_heat(
            &mut da,
            &mut limited,
            &mut pa,
            Time::new(0.005),
            ReactiveStepControl {
                maximum_coupling_retries: 0,
                ..ReactiveStepControl::default()
            },
            DissipationControl::default(),
            WallControl::default(),
            load(),
            piston_controls(),
            || false,
        )
        .unwrap_err();
    assert!(matches!(refusal, ReactiveError::RetryLimit { .. }));
    assert_eq!(a.accepted.states, cells);
    assert_eq!(a.flow().time(), time);
    assert_eq!(a.flow().mesh(), &mesh);
    assert_eq!(da.mesh(), &mesh);
    assert_eq!(limited.mesh(), &mesh);
    assert_eq!(limited.contacts(), &[ADIABATIC, contact]);
    assert_eq!(limited.last_report(), None);
    assert_eq!(pa, piston);
}

#[test]
fn g3_thermal_piston_refuses_unrepresented_enthalpy_work_and_duplicate_drag() {
    let flow = closed(&mechanism(0.0, 0.0, 0.0), 1, chemistry());
    let source = caloric(0.0, ThermalStorageBasis::Enthalpy, 0.01);
    let contact = source.initial_contact(0.5, WallDrag::None).unwrap();
    assert!(matches!(
        PistonThermalEnds::new(
            flow.flow(),
            [ADIABATIC, contact],
            [None, Some(source)],
            || false,
        ),
        Err(WallExchangeError::Invalid(_))
    ));
    let dragged = WallContact {
        drag: WallDrag::Linear { rate_per_s: 1.0 },
        ..finite(3.0, 0.4, 0.5)
    };
    assert!(matches!(
        PistonThermalEnds::new(flow.flow(), [ADIABATIC, dragged], [None, None], || false,),
        Err(WallExchangeError::Invalid(_))
    ));
    assert!(matches!(
        PistonThermalEnds::new(flow.flow(), [ADIABATIC, contact], [None, None], || false,),
        Err(WallExchangeError::Invalid(_))
    ));
}
