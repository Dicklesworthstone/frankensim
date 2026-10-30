use super::*;
use crate::compartment::pyrolysis::{VolatileMixture, VolatileTransferControl};
use crate::compartment::radiation::spectral::{SpectralBand, SpectralOptics, SpectralSolve};
use crate::compartment::radiation::{
    AbsorptionCoefficient, LayerAbsorption, OpticalPath, OpticalSegment, RadiationSolve,
};
use crate::compartment::wall::{
    HeatTransferCoefficient, SpecificEnthalpy, ThermalConductivity, WallGeometry,
};
use crate::compartment::{Ambient, NetworkControl, NetworkError, Room, ZoneNetwork};
use fs_conduction::radiation::{EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, ViewFactorMatrix};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Acceleration, Area, Density, Length, Pressure};
use fs_thermal::enthalpy::pyrolysis::slab::SlabGeometry;
use fs_thermal::enthalpy::pyrolysis::{
    Arrhenius, CondensedCaloric, PyrolysisMaterial, PyrolysisSimulation, VolatileEnthalpy,
};
use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event, MaterialId,
    StepConfig, ThermalLaw, ThermalLink, TimedEvent, series_conductance_w_k,
};

#[path = "../../../../examples/compartment_fire/air.rs"]
mod air;

const AREA: f64 = 0.01;
const LENGTH: f64 = 0.02;
const DENSITY: f64 = 500.0;
const CONDUCTIVITY: f64 = 1.0;

struct Caloric {
    offset: f64,
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + (h - self.offset) / 1000.0;
        if (200.0..=2200.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.001
    }
}
impl CondensedCaloric for Caloric {
    fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, DepositionError> {
        if (200.0..=2200.0).contains(&t) {
            Ok(self.offset + (t - 300.0) * 1000.0)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
}
static LAW: Caloric = Caloric { offset: 0.0 };

fn thermal_boundary(power: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: power,
    }
}
fn slab<'m>(
    law: &'m Caloric,
    volatile: &'m dyn VolatileEnthalpy,
    rate: f64,
    t: f64,
    cells: usize,
) -> RecedingSlab<'m> {
    let dx = LENGTH / cells as f64;
    let mut events = Vec::new();
    for i in 0..cells {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: 1,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: DENSITY * AREA * dx,
                specific_h_j_kg: law.specific_enthalpy_j_kg(t).unwrap(),
                geometry: CellGeometry {
                    start_m: [i as f64 * dx, 0.0, 0.0],
                    end_m: [(i + 1) as f64 * dx, 0.0, 0.0],
                    volume_m3: AREA * dx,
                    surface_area_m2: 2.0 * AREA,
                },
                boundary: thermal_boundary(0.0),
            }),
        });
    }
    for i in 1..cells {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: i as u64,
                a: CellId(i as u64 - 1),
                b: CellId(i as u64),
                area_m2: AREA,
                conductance_w_k: series_conductance_w_k(AREA, [0.5 * dx; 2], [CONDUCTIVITY; 2])
                    .unwrap(),
            }),
        });
    }
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: cells,
            max_links: cells,
            max_events: 2 * cells,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let original = PyrolysisSimulation::new(
        thermal,
        &[PyrolysisMaterial {
            material_id: MaterialId(1),
            virgin: law,
            char_caloric: None,
            volatile,
            char_mass_yield: 0.0,
            kinetics: Arrhenius {
                prefactor_s_inv: rate,
                activation_energy_j_mol: 0.0,
            },
            temperature_range_k: [200.0, 1000.0],
            source: "manufactured conservation/reference fixture",
        }],
        || Ok(()),
    )
    .unwrap();
    RecedingSlab::new(
        original,
        SlabGeometry {
            area_m2: AREA,
            density_kg_m3: DENSITY,
            conductivity_w_m_k: CONDUCTIVITY,
            source: "authored homogeneous fixture",
        },
        || Ok(()),
    )
    .unwrap()
}
fn wall(t: f64, cells: usize) -> Wall1d<'static> {
    Wall1d::new(
        &LAW,
        WallGeometry {
            area: Area::new(AREA),
            thickness: Length::new(LENGTH),
            density: Density::new(DENSITY),
            conductivity: ThermalConductivity::new(CONDUCTIVITY),
            cells,
        },
        SpecificEnthalpy::new((t - 300.0) * 1000.0),
        || false,
    )
    .unwrap()
}
fn boundary(temperatures: [Temperature; 2], h: f64, external: f64) -> WallBoundary {
    WallBoundary {
        gas: temperatures.map(|temperature| GasContact {
            area: Area::new(AREA / 2.0),
            temperature,
            htc: HeatTransferCoefficient::new(h),
        }),
        ambient_temperature: Temperature::new(300.0),
        ambient_htc: HeatTransferCoefficient::new(0.0),
        external_radiant_power: Power::new(external),
    }
}
fn slab_boundary(boundary: WallBoundary) -> SlabBoundary {
    SlabBoundary {
        gas: boundary.gas,
        external_radiant_power: boundary.external_radiant_power,
    }
}
fn controls(dt: f64) -> SlabSurfaceSolve {
    SlabSurfaceSolve {
        pyrolysis: PyrolysisControls {
            maximum_dt_s: dt,
            ..Default::default()
        },
        ..Default::default()
    }
}
fn model() -> ParticipatingEnclosure {
    ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(AREA).unwrap(),
        vec![OpticalPath {
            surfaces: [0, 1],
            exchange_area: Area::new(AREA),
            segments: vec![
                OpticalSegment {
                    layer: Layer::Lower,
                    length: Length::new(0.5),
                },
                OpticalSegment {
                    layer: Layer::Upper,
                    length: Length::new(0.5),
                },
            ],
        }],
        "authored reciprocal normal chord",
        1e-12,
        || false,
    )
    .unwrap()
}
fn absorption() -> [LayerAbsorption; 2] {
    [0.2, 0.6]
        .map(|k| LayerAbsorption::new(AbsorptionCoefficient::new(k), "authored optics").unwrap())
}
fn emissivity(t: Temperature) -> SurfaceEmissivity {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: 0.8,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 2200.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "authored gray analytic surface".into(),
                license: "internal-test".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "analytic".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    SurfaceEmissivity::from_card("test", &card, t.value(), SelectionPolicy::SingleClaimOnly)
        .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn slab_result<'a, 'm>(c: &'a EnclosureSurfaceCandidate<'m, '_>) -> &'a SlabSurfaceCandidate<'m> {
    match &c.surfaces[0] {
        ThermalSurfaceCandidate::Slab(s) => s,
        _ => panic!("slab expected"),
    }
}
fn wall_result<'a, 'm>(c: &'a EnclosureSurfaceCandidate<'m, '_>) -> &'a WallCandidate<'m> {
    match &c.surfaces[1] {
        ThermalSurfaceCandidate::Wall(w) => w,
        _ => panic!("wall expected"),
    }
}

#[test]
fn inert_slab_matches_original_wall_and_heterogeneous_wall_only_path() {
    let volatile = |t| LAW.specific_enthalpy_j_kg(t);
    let solid = slab(&LAW, &volatile, 0.0, 700.0, 3);
    let walls = [wall(700.0, 3), wall(300.0, 3)];
    let model = model();
    let gas = [Temperature::new(500.0), Temperature::new(900.0)];
    let boundaries = [boundary(gas, 8.0, 10.0), boundary(gas, 8.0, 0.0)];
    let kappa = absorption();
    let optics = WallOptics::Gray {
        absorption: &kappa,
        solve: RadiationSolve::default(),
    };
    let original = model
        .candidate_walls(
            &walls,
            Time::new(0.1),
            &boundaries,
            gas,
            optics,
            |_, t| Ok(emissivity(t)),
            None,
            EnclosureWallSolve::default(),
            &mut model.wall_workspace(&walls).unwrap(),
            || false,
        )
        .unwrap();
    let inert = [
        ThermalSurface::Wall {
            wall: &walls[0],
            boundary: boundaries[0],
        },
        ThermalSurface::Wall {
            wall: &walls[1],
            boundary: boundaries[1],
        },
    ];
    let generic = model
        .candidate_surfaces(
            &inert,
            Time::new(0.1),
            gas,
            optics,
            |_, t| Ok(emissivity(t)),
            None,
            EnclosureWallSolve::default(),
            &mut model.surface_workspace(&inert).unwrap(),
            || false,
        )
        .unwrap();
    assert_eq!(original.report, generic.report);
    let surfaces = [
        ThermalSurface::Slab {
            slab: &solid,
            boundary: slab_boundary(boundaries[0]),
            solve: controls(0.1),
        },
        inert[1],
    ];
    let mixed = model
        .candidate_surfaces(
            &surfaces,
            Time::new(0.1),
            gas,
            optics,
            |_, t| Ok(emissivity(t)),
            None,
            EnclosureWallSolve::default(),
            &mut model.surface_workspace(&surfaces).unwrap(),
            || false,
        )
        .unwrap();
    for (i, cell) in slab_result(&mixed)
        .slab
        .state()
        .cells()
        .values()
        .enumerate()
    {
        near(
            slab_result(&mixed)
                .slab
                .state()
                .temperature_k(cell.id)
                .unwrap(),
            original.walls[0].wall.cell_temperature(i).unwrap().value(),
            2e-7,
        );
    }
    near(
        slab_result(&mixed).surface_temperature.value(),
        original.walls[0].report.surface_temperatures[0].value(),
        2e-7,
    );
    near(
        wall_result(&mixed).wall.simulation().energy_j(),
        original.walls[1].wall.simulation().energy_j(),
        2e-6,
    );
    assert_eq!(slab_result(&mixed).report.pyrolysis.volatile_mass_kg, 0.0);
    assert_eq!(mixed.report.exported_volatile_enthalpy.value(), 0.0);
}

#[test]
fn heated_reactive_surface_closes_original_mass_enthalpy_and_endpoint_flux() {
    let volatile = |t| Ok(LAW.specific_enthalpy_j_kg(t)? + 200_000.0);
    let solid = slab(&LAW, &volatile, 0.05, 600.0, 4);
    let inert = wall(400.0, 4);
    let model = model();
    let gas = [Temperature::new(700.0), Temperature::new(900.0)];
    let b = boundary(gas, 10.0, 500.0);
    let surfaces = [
        ThermalSurface::Slab {
            slab: &solid,
            boundary: slab_boundary(b),
            solve: controls(0.1),
        },
        ThermalSurface::Wall {
            wall: &inert,
            boundary: boundary(gas, 10.0, 0.0),
        },
    ];
    let kappa = absorption();
    let result = model
        .candidate_surfaces(
            &surfaces,
            Time::new(0.1),
            gas,
            WallOptics::Gray {
                absorption: &kappa,
                solve: RadiationSolve::default(),
            },
            |_, t| Ok(emissivity(t)),
            None,
            EnclosureWallSolve::default(),
            &mut model.surface_workspace(&surfaces).unwrap(),
            || false,
        )
        .unwrap();
    let slab = slab_result(&result);
    assert!(slab.report.pyrolysis.volatile_mass_kg > 0.0);
    assert!(slab.slab.recession_m() > 0.0);
    near(
        slab.slab.mass_kg() + slab.report.pyrolysis.volatile_mass_kg,
        solid.mass_kg(),
        1e-14,
    );
    let incoming = result.input_powers[0].value() * 0.1;
    near(
        slab.slab.energy_j() - solid.energy_j() + slab.report.pyrolysis.volatile_enthalpy_j,
        incoming,
        1e-6,
    );
    assert_eq!(
        result.report.exported_volatile_enthalpy.value(),
        slab.report.pyrolysis.volatile_enthalpy_j
    );
    assert_eq!(
        slab.coupled_radiant_energy.value(),
        0.1 * result.report.radiation.surface_heat_power()[0].value()
    );
    near(
        result.report.wall_storage_change.value()
            + result.report.exported_volatile_enthalpy.value()
            + result
                .report
                .gas_heat_energy
                .iter()
                .map(|e| e.value())
                .sum::<f64>()
            - result.report.external_radiant_energy.value(),
        result.report.combined_energy_residual_j,
        1e-10,
    );
    assert!(
        result
            .report
            .surface_boundary_power_residual_w
            .iter()
            .all(|v| v.abs() <= 1e-4)
    );
    assert!(result.report.maximum_temperature_residual_k <= 1e-7);
    assert_eq!(solid.time_s(), 0.0);
}

fn network_run(maximum_dt: f64, force_retry: bool) -> (Vec<f64>, EnclosureWallReport<'static>) {
    let (gas, air) = air::source().unwrap();
    let fractions = vec![1.0, 0.0];
    let volatile = VolatileMixture::new(&gas, &fractions, || false).unwrap();
    let law = Caloric {
        offset: volatile.specific_enthalpy_j_kg(300.0).unwrap() - 200_000.0,
    };
    let solid = slab(&law, &volatile, 0.001, 600.0, 3);
    let wall = wall(300.0, 3);
    let room = TwoZone::new(
        gas,
        Area::new(2.0),
        Length::new(2.0),
        Length::new(1.0),
        Pressure::new(101325.0),
        [Temperature::new(400.0), Temperature::new(800.0)],
        [air.clone(), air.clone()],
        || false,
    )
    .unwrap();
    let initial_energy =
        room.physical_internal_energy().value() + solid.energy_j() + wall.simulation().energy_j();
    let initial_mass =
        room.layer(Layer::Lower).mass_kg() + room.layer(Layer::Upper).mass_kg() + solid.mass_kg();
    let mut network = ZoneNetwork::new(
        vec![Room {
            gas: room,
            floor_elevation: Length::new(0.0),
        }],
        Ambient {
            pressure: Pressure::new(101325.0),
            reference_height: Length::new(0.0),
            temperature: Temperature::new(300.0),
            mass_fractions: air,
        },
        vec![],
        vec![],
        Acceleration::new(9.81),
        || false,
    )
    .unwrap();
    let model = model();
    let gas_t = [Temperature::new(400.0), Temperature::new(800.0)];
    let surfaces = [
        ThermalSurface::Slab {
            slab: &solid,
            boundary: slab_boundary(boundary(gas_t, 8.0, 25.0)),
            solve: controls(maximum_dt),
        },
        ThermalSurface::Wall {
            wall: &wall,
            boundary: boundary(gas_t, 8.0, 0.0),
        },
    ];
    let kappa = absorption();
    let mut work = model.surface_workspace(&surfaces).unwrap();
    let (step, result) = network
        .step_with_sources(
            Time::new(maximum_dt),
            NetworkControl::default(),
            |rooms, dt, rates, cancel| {
                let result = model.candidate_room_surfaces(
                    &rooms[0].gas,
                    &surfaces,
                    dt,
                    WallOptics::Gray {
                        absorption: &kappa,
                        solve: RadiationSolve::default(),
                    },
                    |_, t| Ok(emissivity(t)),
                    None,
                    EnclosureWallSolve::default(),
                    &mut work,
                    &mut *cancel,
                )?;
                result.report.add_gas_heat(&mut rates[0])?;
                for release in slab_result(&result).report.pyrolysis.releases.values() {
                    volatile
                        .release_rates(
                            &rooms[0].gas,
                            *release,
                            dt,
                            VolatileTransferControl::default(),
                            &mut *cancel,
                        )?
                        .add_to(Layer::Upper, &mut rates[0])?;
                }
                if force_retry && dt.value() > 0.01 {
                    return Err(NetworkError::Invalid(
                        "authored complete coupled-surface retry",
                    ));
                }
                Ok(result)
            },
            || false,
        )
        .unwrap();
    if force_retry {
        assert_eq!(step.halvings, 1);
    }
    let room = &network.rooms()[0].gas;
    let slab = slab_result(&result);
    let wall = wall_result(&result);
    near(
        room.physical_internal_energy().value()
            + slab.slab.energy_j()
            + wall.wall.simulation().energy_j()
            - initial_energy,
        result.report.external_radiant_energy.value(),
        2e-5,
    );
    near(
        room.layer(Layer::Lower).mass_kg()
            + room.layer(Layer::Upper).mass_kg()
            + slab.slab.mass_kg(),
        initial_mass,
        1e-12,
    );
    assert_eq!(room.time().value(), slab.slab.time_s());
    assert_eq!(room.time(), wall.wall.time());
    let mut state = vec![
        room.pressure().value(),
        room.interface_height().value(),
        slab.slab.energy_j(),
        slab.slab.mass_kg(),
        slab.slab.recession_m(),
        wall.wall.simulation().energy_j(),
    ];
    for layer in [Layer::Lower, Layer::Upper] {
        state.push(room.layer(layer).temperature().value());
        state.extend_from_slice(room.layer(layer).species_masses());
    }
    (state, result.report)
}

#[test]
fn actual_nasa_volatile_transfer_and_network_retry_close_all_three_owners() {
    let (retry, retry_report) = network_run(0.02, true);
    let (direct, direct_report) = network_run(0.01, false);
    assert_eq!(
        retry.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        direct.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
    );
    assert_eq!(retry_report, direct_report);
    assert!(retry_report.exported_volatile_enthalpy.value() > 0.0);
}

#[test]
fn spectral_reactive_coupling_and_refusals_reuse_original_owners() {
    let volatile = |t| Ok(LAW.specific_enthalpy_j_kg(t)? + 200_000.0);
    let mut solid = slab(&LAW, &volatile, 0.01, 600.0, 2);
    let inert = wall(300.0, 2);
    let model = model();
    let gas = [Temperature::new(600.0), Temperature::new(800.0)];
    let b = boundary(gas, 8.0, 20.0);
    let kappa = absorption();
    let optics = SpectralOptics::new(
        vec![SpectralBand {
            wavelengths: [Length::new(1e-7), Length::new(1e-3)],
            absorption: kappa,
        }],
        [None, None],
        "authored constant gray emissivity over wavelength",
    )
    .unwrap();
    let controls = EnclosureWallSolve::default();
    let surfaces = [
        ThermalSurface::Slab {
            slab: &solid,
            boundary: slab_boundary(b),
            solve: self::controls(0.01),
        },
        ThermalSurface::Wall {
            wall: &inert,
            boundary: boundary(gas, 8.0, 0.0),
        },
    ];
    let mut work = model.surface_workspace(&surfaces).unwrap();
    let mut polls = 0;
    assert!(matches!(
        model.candidate_surfaces(
            &surfaces,
            Time::new(0.01),
            gas,
            WallOptics::Spectral {
                optics: &optics,
                solve: SpectralSolve::default()
            },
            |_, t| Ok(emissivity(t)),
            None,
            controls,
            &mut work,
            || {
                polls += 1;
                polls >= 40
            }
        ),
        Err(EnclosureWallError::Cancelled)
    ));
    assert!(
        polls >= 40,
        "cancel after entering an original thermal candidate"
    );
    let result = model
        .candidate_surfaces(
            &surfaces,
            Time::new(0.01),
            gas,
            WallOptics::Spectral {
                optics: &optics,
                solve: SpectralSolve::default(),
            },
            |_, t| Ok(emissivity(t)),
            None,
            controls,
            &mut work,
            || false,
        )
        .unwrap();
    assert!(result.report.work.wavelength_evaluations > 0);
    assert!(result.report.combined_energy_residual_j.abs() <= 1e-4);
    assert!(matches!(
        model.candidate_surfaces(
            &surfaces,
            Time::new(0.01),
            gas,
            WallOptics::Spectral {
                optics: &optics,
                solve: SpectralSolve::default()
            },
            |_, t| Ok(emissivity(t)),
            None,
            EnclosureWallSolve {
                maximum_wall_evaluations: 1,
                ..controls
            },
            &mut work,
            || false
        ),
        Err(EnclosureWallError::Budget { .. })
    ));
    solid
        .set_boundary(solid.back_cell_id(), thermal_boundary(1.0), || Ok(()))
        .unwrap();
    let invalid = [
        ThermalSurface::Slab {
            slab: &solid,
            boundary: slab_boundary(b),
            solve: self::controls(0.01),
        },
        ThermalSurface::Wall {
            wall: &inert,
            boundary: boundary(gas, 8.0, 0.0),
        },
    ];
    assert!(matches!(
        model.surface_workspace(&invalid),
        Err(EnclosureWallError::Invalid(
            "radiating slab requires insulated back"
        ))
    ));
}
