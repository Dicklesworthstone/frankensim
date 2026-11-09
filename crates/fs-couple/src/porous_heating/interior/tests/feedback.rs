//! Original material/FE/thermal owners, with an independently orchestrated
//! held-source reference. Authored storage constants are numerical fixtures.

use super::*;
use crate::cavity_modal::interior::heterogeneous::porous::{
    AcousticVolumeRegion, PorousRegionsInterior,
};
use crate::porous_heating::PorousHeatingError;
use crate::porous_heating::feedback::{IsobaricPorousGas, PorousTemperatureFluid};
use crate::porous_heating::interior::feedback::*;
use fs_feec::acoustics::helmholtz::heterogeneous::time::mechanical::MechanicalBoundary;
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_phs::{WallPin, acoustic_outlet::AcousticWallRelaxation};
use std::cell::Cell;
use std::rc::Rc;

mod molecular;
mod sources;
mod walls;

static SPEC: GasSpec = GasSpec::dry_air_ussa1976();
const SCALE: f64 = 0.03;
const FREQUENCY: f64 = 800.0;
const WALL: CellId = CellId(55);
const WALL_AREA: f64 = 0.5 * SCALE * SCALE;

#[derive(Clone)]
struct Fluid {
    source: IsobaricPorousGas<'static>,
    calls: Rc<Cell<usize>>,
}
impl PorousTemperatureFluid for Fluid {
    type State = GasState;
    fn evaluate(&self, temperature: f64) -> Result<Self::State, PorousHeatingError> {
        self.calls.set(self.calls.get() + 1);
        PorousTemperatureFluid::evaluate(&self.source, temperature)
    }
    fn temperature_k(&self, state: &Self::State) -> f64 {
        state.temperature
    }
    fn volume_response(
        &self,
        state: &Self::State,
        material: Option<&PorousFluidLaw>,
        omega: f64,
    ) -> Result<
        crate::cavity_modal::interior::heterogeneous::porous::AcousticVolumeRegionState,
        fs_material::MaterialError,
    > {
        self.source.volume_response(state, material, omega)
    }
}
fn fluid(support: [f64; 2]) -> Fluid {
    Fluid {
        source: IsobaricPorousGas::try_new(&SPEC, 101325.0, support).unwrap(),
        calls: Rc::new(Cell::new(0)),
    }
}
fn fixture() -> (TetComplex, Vec<[f64; 3]>) {
    let (mesh, mut points) = mesh();
    for point in &mut points {
        *point = point.map(|x| x * SCALE);
    }
    (mesh, points)
}
fn laws() -> [PorousFluidLaw; 2] {
    core::array::from_fn(|index| {
        let mut card = PorousCard::melamine_uf_nguyen2024();
        card.name = "authored full-volume thermal feedback specimen";
        card.provenance = "numerical source composition, not specimen calibration";
        card.flow_resistivity *= [1.0, 1.7][index];
        PorousFluidLaw::try_new(card, 1.81e-5).unwrap()
    })
}
fn mapping<'t, 'm>(
    law: &'t dyn ThermalLaw,
    porous: Option<&'m [PorousFluidLaw; 2]>,
) -> [TemperatureInteriorCell<'t, 'm>; 2] {
    core::array::from_fn(|index| TemperatureInteriorCell {
        destination: InteriorHeatCell {
            cell: IDS[index],
            caloric: law,
        },
        porous: porous.map(|laws| &laws[index]),
    })
}
fn wall_routes(law: &dyn ThermalLaw) -> [InteriorWallHeatCell<'_>; 1] {
    [InteriorWallHeatCell {
        cell: WALL,
        caloric: law,
        acoustic_area_m2: WALL_AREA,
    }]
}
fn problem(with_wall: bool) -> InteriorTemperatureAcoustics {
    InteriorTemperatureAcoustics {
        frequency_hz: FREQUENCY,
        inward_velocity: vec![([0, 2, 3], [C64::new(1.5, -0.4); 3])],
        boundary_drives: Vec::new(),
        volume_drives: Vec::new(),
        traction_drives: Vec::new(),
        admittances: Vec::new(),
        mechanical_walls: if with_wall {
            vec![MechanicalBoundary {
                faces: vec![[0, 1, 3]],
                wall: WallPin {
                    surface_density: 0.6,
                    stiffness_per_area: 2e8,
                    resistance: 500.0,
                },
                relaxation: vec![AcousticWallRelaxation {
                    stiffness_per_area: 5e7,
                    relaxation_time_s: 0.002,
                }],
            }]
        } else {
            Vec::new()
        },
        pressure_release_faces: Vec::new(),
        controls: EquivalentFluidControls::default(),
    }
}
fn events(h: [f64; 2], authored: [f64; 2], with_wall: bool) -> Vec<TimedEvent> {
    let mut events = schedule(authored);
    for (index, event) in events.iter_mut().enumerate() {
        let Event::Birth(cell) = &mut event.event else {
            panic!("birth")
        };
        cell.specific_h_j_kg = h[index];
        cell.geometry.start_m = cell.geometry.start_m.map(|x| x * SCALE);
        cell.geometry.end_m = cell.geometry.end_m.map(|x| x * SCALE);
        cell.geometry.volume_m3 *= SCALE.powi(3);
        cell.geometry.surface_area_m2 *= SCALE * SCALE;
        cell.mass_kg = 30.0 * cell.geometry.volume_m3;
    }
    if with_wall {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: WALL,
                segment_id: 210,
                body_id: 9,
                material_id: MaterialId(6),
                mass_kg: 0.6 * WALL_AREA,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.0, 0.001, 0.0],
                    volume_m3: WALL_AREA * 0.001,
                    surface_area_m2: 2.0 * WALL_AREA + 0.0001,
                },
                boundary: boundary(0.05),
            }),
        });
    }
    events
}

fn original_field(
    geometry: &HelmholtzFe<'_>,
    thermal: &DepositionSimulation<'_>,
    laws: &[PorousFluidLaw; 2],
    problem: &InteriorTemperatureAcoustics,
) -> fs_feec::acoustics::helmholtz::heterogeneous::EquivalentFluidPressureField {
    let states = IDS
        .map(|id| GasState::try_new(&SPEC, thermal.temperature_k(id).unwrap(), 101325.0).unwrap());
    let sources = core::array::from_fn::<_, 2, _>(|index| AcousticVolumeRegion::Porous {
        law: &laws[index],
        fluid: &states[index],
    });
    let original = PorousRegionsInterior::new_with_coupled_mechanical_walls(
        geometry,
        problem.frequency_hz,
        &sources,
        &[0, 1],
        &problem.admittances,
        &problem.mechanical_walls,
        &problem.pressure_release_faces,
        problem.controls,
        &mut || false,
    )
    .unwrap();
    original
        .operator()
        .solve(&problem.inward_velocity, &mut || false)
        .unwrap()
}

#[test]
fn distinct_native_temperatures_refresh_original_materials_and_complete_wall_field() {
    let (mesh, positions) = fixture();
    let geometry = geometry(&mesh, &positions);
    let law = caloric();
    let porous = laws();
    let thermal = initial(events([0.0, 40.0 * CP], [0.12, -0.01], true), &law);
    let source = fluid([280.0, 500.0]);
    let counter = source.calls.clone();
    let problem = problem(true);
    let heating = TemperatureInteriorHeating::try_new(
        &geometry,
        &thermal,
        &mapping(&law, Some(&porous)),
        &wall_routes(&law),
        source,
        problem.clone(),
        3,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(counter.get(), 2);
    assert!(Rc::ptr_eq(&heating.source().calls, &counter));
    for index in 0..2 {
        let gas =
            GasState::try_new(&SPEC, thermal.temperature_k(IDS[index]).unwrap(), 101325.0).unwrap();
        assert_eq!(heating.fluid_states()[index], gas);
        assert!(core::ptr::eq(
            heating.cells()[index].porous.unwrap(),
            &porous[index]
        ));
        let expected = porous[index]
            .effective_fluid(&gas, core::f64::consts::TAU * FREQUENCY)
            .unwrap();
        assert_eq!(
            heating.material_responses()[index]
                .coefficients()
                .dynamic_bulk_modulus_pa,
            expected.dynamic_bulk_modulus_pa()
        );
        assert_eq!(
            heating.material_responses()[index]
                .coefficients()
                .dynamic_density_kg_m3,
            expected.dynamic_density_kg_m3()
        );
    }
    let reference = original_field(&geometry, &thermal, &porous, &problem);
    assert_eq!(heating.field().pressure_pa, reference.pressure_pa);
    assert_eq!(
        heating.field().walls[0].normal_velocity_m_s,
        reference.walls[0].normal_velocity_m_s
    );
    assert_eq!(heating.volume_heating().len(), 2);
    assert_eq!(heating.wall_heating().len(), 1);
    assert!(heating.field().power_balance_residual_w().abs() < 1e-9);
}

#[test]
fn accepted_start_heat_and_endpoint_field_share_native_events_and_split_history() {
    let (mesh, positions) = fixture();
    let geometry = geometry(&mesh, &positions);
    let law = caloric();
    let porous = laws();
    let mut schedule = events([0.0, 40.0 * CP], [0.12, -0.01], true);
    schedule.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 3,
            a: IDS[0],
            b: WALL,
            area_m2: WALL_AREA,
            conductance_w_k: 0.01,
        }),
    });
    schedule.push(TimedEvent {
        time_s: 0.125,
        event: Event::Boundary(WALL, boundary(0.08)),
    });
    let thermal = initial(schedule, &law);
    let original_cells = thermal.cells().clone();
    let problem = problem(true);
    let heating = TemperatureInteriorHeating::try_new(
        &geometry,
        &thermal,
        &mapping(&law, Some(&porous)),
        &wall_routes(&law),
        fluid([280.0, 500.0]),
        problem.clone(),
        3,
        || Ok(()),
    )
    .unwrap();
    let result = heating
        .advance_to(&thermal, 0.25, config(0.0625), || Ok(()))
        .unwrap();
    assert_eq!(result.work.field_solves, 4);
    assert_eq!(result.work.fluid_evaluations, 8);
    assert_eq!(result.work.porous_evaluations, 8);
    assert_eq!(result.step.heat_report.substeps, 4);
    assert_ne!(
        result.heating.field().pressure_pa,
        heating.field().pressure_pa
    );

    // Direct reference through the original source adapter and frozen native
    // heating owner. Each solve occurs BEFORE its one physical thermal interval.
    let mut reference = thermal.clone();
    let mut expected = [0.0; 4];
    for index in 1..=4 {
        let states = IDS.map(|id| {
            GasState::try_new(&SPEC, reference.temperature_k(id).unwrap(), 101325.0).unwrap()
        });
        let sources = core::array::from_fn::<_, 2, _>(|i| AcousticVolumeRegion::Porous {
            law: &porous[i],
            fluid: &states[i],
        });
        let original = PorousRegionsInterior::new_with_coupled_mechanical_walls(
            &geometry,
            FREQUENCY,
            &sources,
            &[0, 1],
            &[],
            &problem.mechanical_walls,
            &[],
            problem.controls,
            &mut || false,
        )
        .unwrap();
        let field = original
            .operator()
            .solve(&problem.inward_velocity, &mut || false)
            .unwrap();
        for (sum, watts) in expected.iter_mut().zip([
            field.density_dissipation_w,
            field.compressional_dissipation_w,
            field.wall_accounting.resistance_dissipation_w,
            field.wall_accounting.relaxation_dissipation_w,
        ]) {
            *sum += watts * 0.0625;
        }
        let held = InteriorAcousticHeating::try_new_with_mechanical_walls(
            original.operator(),
            &field.pressure_pa,
            &field.walls,
            &reference,
            &routes(&law),
            &wall_routes(&law),
            3,
            || Ok(()),
        )
        .unwrap();
        reference = held
            .advance_to(&reference, index as f64 * 0.0625, config(0.0625), || Ok(()))
            .unwrap()
            .thermal;
    }
    for id in IDS.into_iter().chain([WALL]) {
        close(
            result.step.thermal.cells()[&id].specific_h_j_kg,
            reference.cells()[&id].specific_h_j_kg,
            2e-10,
        );
        assert_eq!(
            result.step.thermal.cells()[&id].geometry,
            original_cells[&id].geometry
        );
        assert_eq!(
            result.step.thermal.cells()[&id].mass_kg,
            original_cells[&id].mass_kg
        );
    }
    let actual = result.step.acoustic;
    for (actual, expected) in [
        actual.density_dissipation_j,
        actual.compressional_dissipation_j,
        actual.wall_resistance_dissipation_j,
        actual.wall_relaxation_dissipation_j,
    ]
    .into_iter()
    .zip(expected)
    {
        close(actual, expected, 2e-12);
    }
    close(actual.deposited_j, expected.iter().sum(), 2e-12);
    close(
        result.step.heat_report.external_heat_j,
        actual.deposited_j + (0.12 - 0.01) * 0.25 + 0.05 * 0.125 + 0.08 * 0.125,
        2e-11,
    );
    close(result.step.heat_report.residual_j, 0.0, 2e-11);
    let final_reference = original_field(&geometry, &reference, &porous, &problem);
    for (actual, expected) in result
        .heating
        .field()
        .pressure_pa
        .iter()
        .zip(final_reference.pressure_pa)
    {
        assert!((*actual - expected).abs() < 1e-8);
    }
    for (index, state) in result.heating.fluid_states().iter().enumerate() {
        assert_eq!(
            state.temperature,
            result.step.thermal.temperature_k(IDS[index]).unwrap()
        );
    }
    let first = heating
        .advance_to(&thermal, 0.125, config(0.0625), || Ok(()))
        .unwrap();
    let second = first
        .heating
        .advance_to(&first.step.thermal, 0.25, config(0.0625), || Ok(()))
        .unwrap();
    assert_eq!(second.step.thermal.cells(), result.step.thermal.cells());
    assert_eq!(
        second.heating.field().pressure_pa,
        result.heating.field().pressure_pa
    );
    assert_eq!(
        first.work.field_solves + second.work.field_solves,
        result.work.field_solves
    );
    close(
        first.step.acoustic.deposited_j + second.step.acoustic.deposited_j,
        actual.deposited_j,
        2e-12,
    );
    assert_eq!(thermal.cells(), &original_cells);
}

struct Plateau;
impl ThermalLaw for Plateau {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && (0.0..=1e5).contains(&h) {
            Ok(300.0)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.0
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([77; 32])
    }
}

#[test]
fn unchanged_temperature_reuses_factor_and_one_changed_cell_refreshes_only_its_source() {
    let (mesh, positions) = fixture();
    let geometry = geometry(&mesh, &positions);
    let plateau = Plateau;
    let porous = laws();
    let thermal = initial(events([0.0; 2], [0.12, 0.1], false), &plateau);
    let source = fluid([280.0, 500.0]);
    let calls = source.calls.clone();
    let heating = TemperatureInteriorHeating::try_new(
        &geometry,
        &thermal,
        &mapping(&plateau, Some(&porous)),
        &[],
        source,
        problem(false),
        2,
        || Ok(()),
    )
    .unwrap();
    let result = heating
        .advance_to(&thermal, 0.25, config(0.0625), || Ok(()))
        .unwrap();
    assert_eq!(result.work, InteriorFeedbackWork::default());
    assert_eq!(calls.get(), 2);
    assert!(core::ptr::eq(heating.prepared(), result.heating.prepared()));
    assert!(result.step.acoustic.deposited_j > 0.0);

    // A lossless free-fluid background changes its EOS coherently when only
    // one externally heated native cell changes temperature. No viscosity is
    // converted into an invented bulk absorption rate.
    let law = caloric();
    let thermal = initial(events([0.0, 40.0 * CP], [0.12, 0.0], false), &law);
    let source = fluid([280.0, 500.0]);
    let calls = source.calls.clone();
    let heating = TemperatureInteriorHeating::try_new(
        &geometry,
        &thermal,
        &mapping(&law, None),
        &[],
        source,
        problem(false),
        2,
        || Ok(()),
    )
    .unwrap();
    let result = heating
        .advance_to(&thermal, 0.25, config(0.0625), || Ok(()))
        .unwrap();
    assert_eq!(result.work.fluid_evaluations, 4);
    assert_eq!(result.work.porous_evaluations, 0);
    assert_eq!(result.work.field_solves, 4);
    assert_eq!(calls.get(), 6);
    assert_eq!(result.heating.fluid_states()[1], heating.fluid_states()[1]);
    assert_eq!(result.step.acoustic.deposited_j, 0.0);
    assert_ne!(
        result.heating.field().pressure_pa,
        heating.field().pressure_pa
    );
}

#[test]
fn native_admission_and_endpoint_source_refusal_publish_no_partial_feedback() {
    let (mesh, positions) = fixture();
    let geometry = geometry(&mesh, &positions);
    let law = caloric();
    let porous = laws();
    let thermal = initial(events([0.0; 2], [2.0, 1.0], false), &law);
    let source = fluid([299.0, 300.1]);
    let calls = source.calls.clone();
    let mut bad_problem = problem(false);
    bad_problem.controls.max_volume_points = 1;
    assert!(
        TemperatureInteriorHeating::try_new(
            &geometry,
            &thermal,
            &mapping(&law, Some(&porous)),
            &[],
            source.clone(),
            bad_problem,
            2,
            || Ok(()),
        )
        .is_err()
    );
    assert_eq!(calls.get(), 0);
    let heating = TemperatureInteriorHeating::try_new(
        &geometry,
        &thermal,
        &mapping(&law, Some(&porous)),
        &[],
        source,
        problem(false),
        2,
        || Ok(()),
    )
    .unwrap();
    let original = thermal.cells().clone();
    let pressure = heating.field().pressure_pa.clone();
    assert!(matches!(
        heating.advance_to(&thermal, 0.01, config(0.01), || Ok(())),
        Err(InteriorFeedbackError::Fluid { element: 0, .. })
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells(), &original);
    assert_eq!(heating.field().pressure_pa, pressure);
    assert_eq!(heating.fluid_states()[0].temperature, 300.0);

    let mut changed = events([0.0; 2], [2.0, 1.0], false);
    let Event::Birth(cell) = &mut changed[0].event else {
        panic!("birth")
    };
    cell.mass_kg *= 1.01;
    let changed = initial(changed, &law);
    let before = calls.get();
    assert!(matches!(
        heating.advance_to(&changed, 0.001, config(0.001), || Ok(())),
        Err(InteriorFeedbackError::Heating(
            InteriorHeatingError::Invalid(_)
        ))
    ));
    assert_eq!(calls.get(), before);
}

#[test]
fn cancellation_during_endpoint_rebuild_leaves_original_field_and_thermal_state_reusable() {
    let (mesh, positions) = fixture();
    let geometry = geometry(&mesh, &positions);
    let law = caloric();
    let porous = laws();
    let thermal = initial(events([0.0, 40.0 * CP], [0.12, -0.01], true), &law);
    let heating = TemperatureInteriorHeating::try_new(
        &geometry,
        &thermal,
        &mapping(&law, Some(&porous)),
        &wall_routes(&law),
        fluid([280.0, 500.0]),
        problem(true),
        3,
        || Ok(()),
    )
    .unwrap();
    let original = thermal.cells().clone();
    let pressure = heating.field().pressure_pa.clone();
    let mut polls = 0;
    let accepted = heating
        .advance_to(&thermal, 0.125, config(0.03125), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        heating.advance_to(&thermal, 0.125, config(0.03125), || {
            calls += 1;
            if calls >= polls * 3 / 4 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(InteriorFeedbackError::Heating(
            InteriorHeatingError::Thermal(DepositionError::Cancelled)
        ))
    ));
    assert_eq!(thermal.cells(), &original);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(heating.field().pressure_pa, pressure);
    let retry = heating
        .advance_to(&thermal, 0.125, config(0.03125), || Ok(()))
        .unwrap();
    assert_eq!(retry.step.thermal.cells(), accepted.step.thermal.cells());
    assert_eq!(
        retry.heating.field().pressure_pa,
        accepted.heating.field().pressure_pa
    );
    assert_eq!(retry.step.acoustic, accepted.step.acoustic);
    assert_eq!(retry.work, accepted.work);
}

#[test]
fn thermal_timestep_refinement_resolves_original_material_and_field_feedback() {
    let (mesh, positions) = fixture();
    let geometry = geometry(&mesh, &positions);
    let law = caloric();
    let porous = laws();
    let thermal = initial(events([0.0, 20.0 * CP], [2.0, 1.0], false), &law);
    let heating = TemperatureInteriorHeating::try_new(
        &geometry,
        &thermal,
        &mapping(&law, Some(&porous)),
        &[],
        fluid([280.0, 500.0]),
        problem(false),
        2,
        || Ok(()),
    )
    .unwrap();
    // Direct reference repeatedly builds the original source/FE owner and
    // applies its physical watts through the native frozen heat interface.
    let dt = 1.0 / 512.0;
    let mut reference = thermal.clone();
    for index in 1..=128 {
        let states = IDS.map(|id| {
            GasState::try_new(&SPEC, reference.temperature_k(id).unwrap(), 101325.0).unwrap()
        });
        let sources = core::array::from_fn::<_, 2, _>(|i| AcousticVolumeRegion::Porous {
            law: &porous[i],
            fluid: &states[i],
        });
        let original = PorousRegionsInterior::new(
            &geometry,
            FREQUENCY,
            &sources,
            &[0, 1],
            &[],
            &[],
            EquivalentFluidControls::default(),
            &mut || false,
        )
        .unwrap();
        let field = original
            .operator()
            .solve(&problem(false).inward_velocity, &mut || false)
            .unwrap();
        let held = InteriorAcousticHeating::try_new(
            original.operator(),
            &field.pressure_pa,
            &reference,
            &routes(&law),
            2,
            || Ok(()),
        )
        .unwrap();
        reference = held
            .advance_to(&reference, index as f64 * dt, config(dt), || Ok(()))
            .unwrap()
            .thermal;
    }
    let mut errors = Vec::new();
    for dt in [0.0625, 0.03125, 0.015625] {
        let result = heating
            .advance_to(&thermal, 0.25, config(dt), || Ok(()))
            .unwrap();
        let error = IDS
            .into_iter()
            .map(|id| {
                (result.step.thermal.temperature_k(id).unwrap()
                    - reference.temperature_k(id).unwrap())
                .abs()
            })
            .fold(0.0, f64::max);
        errors.push(error);
        assert!(result.work.field_solves > 0);
        close(result.step.heat_report.residual_j, 0.0, 1e-9);
    }
    assert!(errors[0] > 1e-9, "feedback axis inactive: {errors:?}");
    assert!(
        errors[1] < 0.7 * errors[0] && errors[2] < 0.7 * errors[1],
        "{errors:?}"
    );
}
