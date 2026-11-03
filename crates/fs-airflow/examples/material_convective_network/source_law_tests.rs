//! Source-pinned property queries survive the original reaction/Joule/motion solve.
use super::*;
use fs_airflow::graph::thermal::medium::convection::finite_gas::{
    film::{FilmGasHeatLaw, FilmGasSurface},
    sphere::ForcedSphereGasSurface,
};
use fs_material::fluid_transport::FluidTransport;

#[test]
fn source_laws_follow_reaction_composition_without_extra_heat_or_electrical_storage() {
    fixture(|binding, thermal, circuit| {
        scoped(|cx, _| {
            let electrical = ElectricalContactNetwork {
                bulk: binding,
                contacts: &[],
                electrical_bases: &BASES,
            };
            let assembly = MaterialExchangeNetwork {
                electrical: &electrical,
                thermal_contacts: &[],
                radiation: None,
            };
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let reaction = chemical(50.0);
            let mut gas = reacting_gases(&reaction);
            // Pure synthetic B has Pr > 0.71. The original A/B half-mixture does
            // not, so it must not be used to fabricate a successful sphere query.
            let inert_bank = gas_source::model(0.0);
            gas.insert(
                20,
                GasPiston::new(
                    ClosedGasParcel::new(
                        inert_bank.clone(),
                        gas_source::point(400.0, 100000.0),
                        0.001,
                    )
                    .unwrap(),
                    1.0,
                    0.01,
                    0.0,
                )
                .unwrap(),
            );
            let saved = gas.clone();
            let initial_energy = storage(&thermal, &circuit, &gas);
            let base = contacts(1.0);
            let patches = [
                ExchangeContact {
                    contact: base[0].contact,
                    law: GasHeatLaw::FilmConvection(FilmGasSurface {
                        model: reaction.catalog(),
                        law: FilmGasHeatLaw::Forced(base[0].surface),
                    }),
                },
                ExchangeContact {
                    contact: base[1].contact,
                    law: GasHeatLaw::FilmConvection(FilmGasSurface {
                        model: reaction.catalog(),
                        law: FilmGasHeatLaw::Natural(NaturalGasSurface {
                            geometry: fs_convection::NaturalConvectionGeometry::HorizontalCylinder,
                            area_m2: base[1].surface.area_m2,
                            characteristic_length: Length::new(0.1),
                            gravity_m_s2: 9.81,
                        }),
                    }),
                },
                ExchangeContact {
                    contact: base[2].contact,
                    law: GasHeatLaw::ForcedSphereConvection(ForcedSphereGasSurface {
                        model: &inert_bank,
                        area_m2: base[2].surface.area_m2,
                        diameter: Length::new(0.2),
                        speed_m_s: 1.0,
                    }),
                },
            ];
            let solve = chemistry(&gas);
            let expected = reaction
                .parcel_candidate(gas[&10].gas(), DT, solve, || false)
                .unwrap();
            let next = exchange::step_candidate(
                &prepared,
                cx,
                &thermal,
                DT,
                &circuit,
                &[],
                &chambers(&gas, DT),
                &[ChamberReaction {
                    chamber: 10,
                    reaction: &reaction,
                    controls: solve,
                }],
                &patches,
                controls(),
                &[0.0; 2],
            )
            .unwrap();
            assert_eq!(&next.reactions()[0].reaction, expected.step());
            let mut accepted = BTreeMap::new();
            let mut outside = 0.0;
            for chamber in next.chambers() {
                for patch in &chamber.surfaces {
                    let cell = next
                        .report()
                        .source
                        .material
                        .material
                        .bulk
                        .cells
                        .iter()
                        .find(|cell| cell.cell == patch.source.0.cell)
                        .unwrap();
                    let (sample, property) = match &patch.source.1 {
                        GasHeatSample::FilmConvection(sample) => {
                            let property = sample.property_gas.as_ref().unwrap();
                            assert!(
                                (property.state().temperature
                                    - 0.5
                                        * (sample.interface.body.temperature_k()
                                            + sample.interface.gas.state().temperature))
                                    .abs()
                                    < 1e-12
                            );
                            assert_eq!(
                                sample.interface.evaluation.as_ref().unwrap().transport,
                                property.transport_state().unwrap()
                            );
                            assert_eq!(sample.interface.gas.components().len(), 2);
                            (&sample.interface, property)
                        }
                        GasHeatSample::ForcedSphereConvection(sample) => {
                            let property = sample.surface_gas.as_ref().unwrap();
                            assert_eq!(property.state().temperature, cell.phase.temperature_k());
                            assert_eq!(
                                sample.interface.evaluation.as_ref().unwrap().transport,
                                sample.interface.gas.transport_state().unwrap()
                            );
                            (&sample.interface, property)
                        }
                        _ => panic!("source-pinned law was replaced"),
                    };
                    assert_eq!(sample.gas, chamber.boundary.piston.midpoint);
                    assert_eq!(sample.body, cell.phase);
                    assert_eq!(property.state().pressure, sample.gas.state().pressure);
                    assert_eq!(
                        property.thermodynamics().receipt().mole_fraction_bits(),
                        sample.gas.thermodynamics().receipt().mole_fraction_bits()
                    );
                    assert_eq!(property.source_packs(), sample.gas.source_packs());
                    let heat = DT
                        * patch.conductance_w_k
                        * (sample.body.temperature_k() - sample.gas.state().temperature);
                    assert!((patch.heat_j - heat).abs() <= 1e-9);
                }
                outside += chamber.boundary.piston.exterior_pressure_work_j;
                accepted.insert(chamber.id, chamber.piston.clone());
            }
            assert_ne!(next.circuit().charges_c(), circuit.charges_c());
            assert_ne!(accepted[&10].velocity_m_s(), gas[&10].velocity_m_s());
            assert!(
                (storage(next.thermal(), next.circuit(), &accepted) + outside - initial_energy)
                    .abs()
                    <= 5e-6
            );
            assert!(
                next.step().total_energy_residual_j.abs() <= controls().total_energy_tolerance_j
            );
            assert_eq!(gas, saved);
            assert_eq!(thermal.time_s(), 0.0);
        })
    });
}
