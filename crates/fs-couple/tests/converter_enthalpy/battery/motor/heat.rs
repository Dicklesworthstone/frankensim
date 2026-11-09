//! Finite coolant calorimetry and source transport on the original drive solve.
//! Material laws are manufactured fixtures, not measured motor/coolant data.
use super::*;
use fs_convection::{CorrelationId, CorrelationInputs, ThermalConductivity};
use fs_couple::circuit_enthalpy::{
    CircuitThermalCoefficients,
    contacts::{
        ContactThermalCoefficients,
        heat::{
            HeatAssembly, HeatNetworkError, HeatNetworkSample, MaterialHeatNetwork,
            PreparedHeatNetwork,
        },
        material::ElectricalContactNetwork,
    },
    converter::battery::motor::heat::{BatteryMotorHeatError, BatteryMotorHeatNetwork},
    material::MaterialNetwork,
};
use fs_exec::{CancelGate, Cx};
use fs_qty::Length;
use std::collections::BTreeMap;

const COOLANT: CellId = CellId(151);
const FLUID_CP: f64 = 4000.0;
const FLUID_MASS: f64 = 0.001;
const PATCH_AREA: f64 = 2e-5;
const LENGTH: f64 = 0.1;

fn coolant() -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = material(1.0, FLUID_CP, 0.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "density" {
            claim.value = PropertyValue::Scalar {
                value: 1000.0,
                dims: Density::DIMS,
            };
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic coolant".into(),
            phase: "liquid".into(),
            process: "authored integration fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(2000.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        0.0,
        1e-8,
    )
    .unwrap();
    (card, curve)
}

fn with_coolant<R>(
    cp: [f64; 10],
    alpha: f64,
    run: impl FnOnce(&Cx<'_>, &CancelGate, &MaterialHeatNetwork<'_>, DepositionSimulation<'_>) -> R,
) -> R {
    let fixture = DriveFixture::new(cp, alpha);
    let (fluid_card, fluid_curve) = coolant();
    let context = point(300.0);
    let fluid_context = point(290.0);
    let mut cells = fixture.cells(&context).into_iter().collect::<Vec<_>>();
    cells.push(MaterialCell {
        cell: COOLANT,
        material_id: MaterialId(10),
        curve: fluid_curve.curve(),
        point: ThermalMaterialPoint {
            material: &fluid_card,
            context: &fluid_context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    });
    let original = fixture.body();
    let mut events = original
        .cells()
        .values()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell.clone()),
        })
        .collect::<Vec<_>>();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: COOLANT,
            segment_id: 10,
            body_id: 10,
            material_id: MaterialId(10),
            mass_kg: FLUID_MASS,
            specific_h_j_kg: 90.0 * FLUID_CP,
            geometry: CellGeometry {
                start_m: [20.0, 0.0, 0.0],
                end_m: [20.01, 0.0, 0.0],
                volume_m3: 1e-6,
                surface_area_m2: 6e-4,
            },
            boundary: Boundary {
                ambient_k: 290.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        }),
    });
    for (i, id) in DRIVE_IDS.into_iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: i as u64,
                a: id,
                b: COOLANT,
                area_m2: PATCH_AREA,
                conductance_w_k: 0.0,
            }),
        });
    }
    let mut laws = fixture
        .curves
        .iter()
        .enumerate()
        .map(|(i, curve)| (MaterialId(i as u64), curve.curve() as &dyn ThermalLaw))
        .collect::<Vec<_>>();
    laws.push((MaterialId(10), fluid_curve.curve() as &dyn ThermalLaw));
    let thermal = DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: 11,
            max_links: 10,
            max_events: 32,
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
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 19,
                kernel_id: 73,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        run(&cx, &gate, &heat, thermal)
    })
}

#[derive(Debug)]
enum CoolingError {
    Heat(HeatNetworkError),
    Correlation(fs_convection::CorrelationError),
    FilmDomain,
    Refused,
}
impl From<HeatNetworkError> for CoolingError {
    fn from(error: HeatNetworkError) -> Self {
        Self::Heat(error)
    }
}
struct CoolingSample {
    heat: HeatNetworkSample,
    temperature_k: [f64; 11],
    conductance_w_k: [f64; 10],
    into_coolant_w: [f64; 10],
    reservoir_heat_w: f64,
    resolution: usize,
}
impl AsRef<HeatNetworkSample> for CoolingSample {
    fn as_ref(&self) -> &HeatNetworkSample {
        &self.heat
    }
}

// This exercises the same admitted link/provider seam as downstream
// fs-airflow::SourceConvectionHeatAssembly without introducing a dependency
// cycle. The original fs-convection owner evaluates every heat coefficient;
// the native thermal cells own BOTH inventories and all time integration.
struct Cooling<'p, 'n, 'a> {
    prepared: &'p PreparedHeatNetwork<'n, 'a>,
    speed_m_s: f64,
    property_slope: f64,
    reservoir_g: f64,
    calls: Cell<usize>,
    refuse_at: Cell<Option<usize>>,
    cancel_at: Cell<Option<usize>>,
    gate: &'p CancelGate,
}
impl<'p, 'n, 'a> Cooling<'p, 'n, 'a> {
    fn new(prepared: &'p PreparedHeatNetwork<'n, 'a>, gate: &'p CancelGate) -> Self {
        Self {
            prepared,
            gate,
            speed_m_s: 0.1,
            property_slope: 0.0,
            reservoir_g: 0.0,
            calls: Cell::new(0),
            refuse_at: Cell::new(None),
            cancel_at: Cell::new(None),
        }
    }
}
impl HeatAssembly for Cooling<'_, '_, '_> {
    type Sample = CoolingSample;
    type Error = CoolingError;
    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self.prepared
    }
    fn resolve_heat(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<CoolingSample>, CircuitEnthalpyError<CoolingError>> {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let resolution = self.calls.get() + 1;
        self.calls.set(resolution);
        if self.refuse_at.get() == Some(resolution) {
            return Err(CircuitEnthalpyError::Source(CoolingError::Refused));
        }
        let mut temperature_k = [0.0; 11];
        for (i, id) in DRIVE_IDS.into_iter().chain([COOLANT]).enumerate() {
            temperature_k[i] = thermal.temperature_k(id)?;
        }
        let mut conductance_w_k = [0.0; 10];
        for i in 0..10 {
            let film = f64::midpoint(temperature_k[i], temperature_k[10]);
            if !(280.0..=320.0).contains(&film) {
                return Err(CircuitEnthalpyError::Source(CoolingError::FilmDomain));
            }
            let k = 0.5 + self.property_slope * 0.002 * (film - 300.0);
            let mu = 0.001 * (1.0 - self.property_slope * 0.01 * (film - 300.0));
            let groups =
                CorrelationInputs::forced(1000.0 * self.speed_m_s * LENGTH / mu, mu * FLUID_CP / k);
            let nu = fs_convection::evaluate(CorrelationId::FlatPlateLaminarAverage, groups)
                .map_err(|error| CircuitEnthalpyError::Source(CoolingError::Correlation(error)))?;
            conductance_w_k[i] = PATCH_AREA
                * nu.heat_transfer_coefficient(ThermalConductivity::new(k), Length::new(LENGTH))
                    .map_err(|error| {
                        CircuitEnthalpyError::Source(CoolingError::Correlation(error))
                    })?
                    .value
                    .value();
        }
        let links = conductance_w_k
            .into_iter()
            .enumerate()
            .map(|(i, g)| (i as u64, g))
            .collect();
        let mut resolved = self
            .prepared
            .resolve_with_link_conductances(cx, thermal, &links)
            .map_err(|error| match error {
                HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
                other => CircuitEnthalpyError::Source(CoolingError::Heat(other)),
            })?;
        let reservoir_heat_w = self.reservoir_g * (285.0 - temperature_k[0]);
        *resolved.bulk.sources_w.entry(DRIVE_IDS[0]).or_default() += reservoir_heat_w;
        if self.cancel_at.get() == Some(resolution) {
            self.gate.request();
        }
        Ok(ContactThermalCoefficients {
            contacts: resolved.contacts,
            bulk: CircuitThermalCoefficients {
                resistors: resolved.bulk.resistors,
                sources_w: resolved.bulk.sources_w,
                conductances_w_k: resolved.bulk.conductances_w_k,
                source: CoolingSample {
                    heat: resolved.bulk.source,
                    temperature_k,
                    conductance_w_k,
                    into_coolant_w: core::array::from_fn(|i| {
                        conductance_w_k[i] * (temperature_k[i] - temperature_k[10])
                    }),
                    reservoir_heat_w,
                    resolution,
                },
            },
        })
    }
}
fn extra_links() -> BTreeMap<u64, f64> {
    (0..10).map(|i| (i, 0.0)).collect()
}
fn cooling_controls() -> CircuitEnthalpySolve {
    let mut controls = drive_controls();
    controls.maximum_cells = 11;
    controls.newton.linear_restart = 11;
    controls
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-8, "{a:.15e} != {b:.15e}");
}

#[test]
fn finite_coolant_matches_eliminated_caloric_equations_and_retains_continuing_drive() {
    let cp = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
    with_coolant(cp, 0.0, |cx, gate, heat, mut thermal| {
        let prepared = heat
            .prepare_with_link_conductances(cx, &thermal, &extra_links())
            .unwrap();
        let provider = Cooling::new(&prepared, gate);
        let key = key();
        let channels = drive_paths(&key);
        let network = BatteryMotorHeatNetwork {
            assembly: &provider,
            channels: &channels,
            bearing_cell: DRIVE_IDS[9],
        };
        let (native, before, torques, load, heat_w) = steady_drive(1.0, 0.1);
        let mut drive = network.bind(cx, &thermal, &native, before).unwrap();
        let bare = BatteryMotorThermalNetwork {
            cells: heat.electrical.bulk.cells,
            links: &[],
            channels: &channels,
            bearing_cell: DRIVE_IDS[9],
        };
        assert!(
            bare.bind(&thermal, &native, before, || false).is_err(),
            "ordinary Fourier bindings must not silently insulate coolant links"
        );
        let dt = 0.02;
        // Independent Nu=0.664 sqrt(Re) cbrt(Pr), with held synthetic properties.
        let g = PATCH_AREA * 0.664 * 100.0 * 2.0 * 0.5 / LENGTH;
        let fluid_capacity = FLUID_MASS * FLUID_CP;
        let mut expected = [300.0; 11];
        expected[10] = 290.0;
        for step in 1..=2 {
            let old = expected;
            let w = cp.map(|capacity| dt * g * capacity / (capacity + dt * g));
            expected[10] = (fluid_capacity * old[10]
                + (0..10)
                    .map(|i| w[i] * old[i] + dt * g * dt * heat_w[i] / (cp[i] + dt * g))
                    .sum::<f64>())
                / (fluid_capacity + w.iter().sum::<f64>());
            for i in 0..10 {
                expected[i] =
                    (cp[i] * old[i] + dt * heat_w[i] + dt * g * expected[10]) / (cp[i] + dt * g);
            }
            let next = network
                .step_candidate(
                    cx,
                    &thermal,
                    dt * step as f64,
                    &drive,
                    native.electrical().modulation(),
                    &torques,
                    load,
                    electrical_controls(),
                    cooling_controls(),
                    &[0.0; 11],
                )
                .unwrap();
            let report = next.report();
            for (i, id) in DRIVE_IDS.into_iter().chain([COOLANT]).enumerate() {
                close(next.thermal().temperature_k(id).unwrap(), expected[i]);
                close(report.source.network.temperature_k[i], expected[i]);
            }
            close(
                fluid_capacity * (expected[10] - old[10]),
                dt * report.source.network.into_coolant_w.iter().sum::<f64>(),
            );
            close(
                next.thermal().energy_j() - thermal.energy_j(),
                report.deposited_joule_j,
            );
            close(report.deposited_joule_j, dt * heat_w.iter().sum::<f64>());
            close(report.non_electrical_heat_j, 0.0);
            close(report.total_energy_residual_j, 0.0);
            for value in report.source.network.conductance_w_k {
                close(value, g);
            }
            let (next_thermal, next_drive, _) = next.into_parts();
            thermal = next_thermal;
            drive = next_drive;
        }
        assert!(thermal.temperature_k(COOLANT).unwrap() > 290.0);
        assert!(drive.state()[0] < before[0]);
        assert_eq!(drive.time_s(), thermal.time_s());
    });
}

#[test]
fn nonlinear_cooling_changes_all_resistances_and_motor_current_with_external_heat_separate() {
    let (_, _, _, _, heat_w) = steady_drive(1.0, 0.1);
    let cp = heat_w.map(|power| 0.02 * power);
    let mut endpoints = Vec::new();
    for speed in [0.025, 1.0] {
        endpoints.push(with_coolant(cp, 0.0039, |cx, gate, heat, thermal| {
            let prepared = heat
                .prepare_with_link_conductances(cx, &thermal, &extra_links())
                .unwrap();
            let mut provider = Cooling::new(&prepared, gate);
            provider.speed_m_s = speed;
            provider.property_slope = 1.0;
            provider.reservoir_g = 0.002;
            let key = key();
            let channels = drive_paths(&key);
            let network = BatteryMotorHeatNetwork {
                assembly: &provider,
                channels: &channels,
                bearing_cell: DRIVE_IDS[9],
            };
            let (native, state, torques, load, _) = steady_drive(1.0, 0.1);
            let drive = network.bind(cx, &thermal, &native, state).unwrap();
            let dt = 0.05;
            let next = network
                .step_candidate(
                    cx,
                    &thermal,
                    dt,
                    &drive,
                    native.electrical().modulation(),
                    &torques,
                    load,
                    electrical_controls(),
                    cooling_controls(),
                    &[0.0; 11],
                )
                .unwrap();
            let report = next.report();
            let sample = &report.source.network;
            for i in 0..10 {
                let t = sample.temperature_k[i];
                let film = 0.5 * (t + sample.temperature_k[10]);
                let k = 0.5 + 0.002 * (film - 300.0);
                let mu = 0.001 * (1.0 - 0.01 * (film - 300.0));
                let g = PATCH_AREA
                    * 0.664
                    * (1000.0 * speed * LENGTH / mu).sqrt()
                    * (mu * FLUID_CP / k).cbrt()
                    * k
                    / LENGTH;
                close(sample.conductance_w_k[i], g);
                let external = if i == 0 {
                    dt * sample.reservoir_heat_w
                } else {
                    0.0
                };
                close(
                    cp[i] * (t - 300.0),
                    report.electrical.heat_joules()[i] - dt * g * (t - sample.temperature_k[10])
                        + external,
                );
                if i < 9 {
                    close(
                        next.drive().channels()[i].resistance_ohm(),
                        DRIVE_R[i] * (1.0 + 0.0039 * (t - 300.0)),
                    );
                }
            }
            close(
                FLUID_MASS * FLUID_CP * (sample.temperature_k[10] - 290.0),
                dt * sample.into_coolant_w.iter().sum::<f64>(),
            );
            assert!(sample.reservoir_heat_w < 0.0);
            close(report.non_electrical_heat_j, dt * sample.reservoir_heat_w);
            close(
                next.thermal().energy_j() - thermal.energy_j(),
                report.deposited_joule_j + report.non_electrical_heat_j,
            );
            close(
                report.deposited_joule_j,
                report.electrical.heat_joules().iter().sum(),
            );
            close(report.total_energy_residual_j, 0.0);
            (sample.temperature_k[8], next.drive().state()[6])
        }));
    }
    assert!(
        endpoints[1].0 < endpoints[0].0,
        "stronger physical cooling lowers winding temperature"
    );
    assert!(
        (endpoints[1].1 - endpoints[0].1).abs() > 1e-7,
        "cooling must reach actual native motor current/flux"
    );
}

#[test]
fn empty_extra_heat_paths_recover_the_original_nine_state_entry_point() {
    let fixture = DriveFixture::new([1.0; 10], 0.0039);
    let context = point(300.0);
    let cells = fixture.cells(&context);
    let thermal = fixture.body();
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 19,
                kernel_id: 73,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let prepared = heat.prepare(&cx, &thermal).unwrap();
        let key = key();
        let channels = drive_paths(&key);
        let network = BatteryMotorHeatNetwork {
            assembly: &prepared,
            channels: &channels,
            bearing_cell: DRIVE_IDS[9],
        };
        let bare = BatteryMotorThermalNetwork {
            cells: &cells,
            links: &[],
            channels: &channels,
            bearing_cell: DRIVE_IDS[9],
        };
        let (native, state, torques, load, _) = steady_drive(1.0, 0.1);
        let drive = network.bind(&cx, &thermal, &native, state).unwrap();
        let complete = network
            .step_candidate(
                &cx,
                &thermal,
                0.02,
                &drive,
                native.electrical().modulation(),
                &torques,
                load,
                electrical_controls(),
                drive_controls(),
                &[0.0; 10],
            )
            .unwrap();
        let original = bare
            .step_candidate(
                &thermal,
                0.02,
                &drive,
                native.electrical().modulation(),
                &torques,
                load,
                electrical_controls(),
                drive_controls(),
                &[0.0; 10],
                || false,
            )
            .unwrap();
        assert_eq!(complete.thermal().cells(), original.thermal().cells());
        assert_eq!(complete.drive().state(), original.drive().state());
        assert_eq!(
            complete.report().electrical.native,
            original.report().electrical.native
        );
        assert_eq!(
            complete.report().electrical.heat_joules(),
            original.report().electrical.heat_joules()
        );
    });
}

#[test]
fn late_cooling_source_refusal_and_cancellation_publish_no_drive_or_coolant_history() {
    with_coolant([1.0; 10], 0.0039, |cx, gate, heat, thermal| {
        let prepared = heat
            .prepare_with_link_conductances(cx, &thermal, &extra_links())
            .unwrap();
        let provider = Cooling::new(&prepared, gate);
        let key = key();
        let channels = drive_paths(&key);
        let network = BatteryMotorHeatNetwork {
            assembly: &provider,
            channels: &channels,
            bearing_cell: DRIVE_IDS[9],
        };
        let (native, state, torques, load, _) = steady_drive(1.0, 0.1);
        let drive = network.bind(cx, &thermal, &native, state).unwrap();
        let before = thermal.cells().clone();
        let attempt = || {
            network.step_candidate(
                cx,
                &thermal,
                0.02,
                &drive,
                native.electrical().modulation(),
                &torques,
                load,
                electrical_controls(),
                cooling_controls(),
                &[0.0; 11],
            )
        };
        provider.calls.set(0);
        let last = attempt().unwrap().report().source.network.resolution;
        assert!(last > 3);
        provider.calls.set(0);
        provider.refuse_at.set(Some(last));
        assert!(matches!(
            attempt(),
            Err(CircuitEnthalpyError::Source(BatteryMotorHeatError::Heat(
                CoolingError::Refused
            )))
        ));
        assert_eq!(thermal.cells(), &before);
        assert_eq!(drive.state(), &state);
        provider.calls.set(0);
        provider.refuse_at.set(None);
        assert_eq!(attempt().unwrap().report().source.network.resolution, last);
        provider.calls.set(0);
        provider.cancel_at.set(Some(last));
        assert!(matches!(attempt(), Err(CircuitEnthalpyError::Cancelled)));
        assert_eq!(thermal.cells(), &before);
        assert_eq!(drive.state(), &state);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(drive.time_s(), 0.0);
    });
}

#[test]
fn switched_drive_heats_finite_coolant_using_the_same_complete_heat_owner() {
    use fs_phs::circuit::converter::switching::battery_motor::{
        BatteryMotorSwitchCommand, PreparedBatteryMotorSwitchCycle,
    };
    let cp = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
    with_coolant(cp, 0.0, |cx, gate, heat, thermal| {
        let prepared = heat
            .prepare_with_link_conductances(cx, &thermal, &extra_links())
            .unwrap();
        let provider = Cooling::new(&prepared, gate);
        let key = key();
        let channels = drive_paths(&key);
        let network = BatteryMotorHeatNetwork {
            assembly: &provider,
            channels: &channels,
            bearing_cell: DRIVE_IDS[9],
        };
        let (native, before, torques, load, _) = steady_drive(1.0, 0.1);
        let drive = network.bind(cx, &thermal, &native, before).unwrap();
        let dt = 0.02;
        let intervals = centered_pwm(dt, native.electrical().modulation()).unwrap();
        let commands = vec![
            BatteryMotorSwitchCommand {
                shaft_torque_nm: torques,
                auxiliary_load: load,
            };
            intervals.len()
        ];
        let cycle = SwitchCycleControls {
            absolute_energy_tolerance_j: 1e-9,
            relative_energy_tolerance: 1e-11,
            ..SwitchCycleControls::default()
        };
        let mut original = PreparedBatteryMotorSwitchCycle::new_controlled(
            native,
            &intervals,
            4,
            electrical_controls(),
            cycle,
            || false,
        )
        .unwrap();
        let mut state = [0.0; 9];
        let native_step = original
            .step_interval_commands_controlled(&before, &commands, &mut state, || false)
            .unwrap();
        let q = native_step.heat_joules();
        let mut controls = cooling_controls();
        controls.maximum_circuit_trials = 20_000;
        let next = network
            .step_switch_candidate(
                cx,
                &thermal,
                dt,
                &drive,
                &intervals,
                &commands,
                4,
                electrical_controls(),
                cycle,
                controls,
                &[0.0; 11],
            )
            .unwrap();
        assert_eq!(*next.drive().state(), state);
        assert_eq!(next.report().electrical.heat_joules(), q);
        // Independently eliminate the ten endpoint body temperatures from
        // Cp*dT = Q - dt*g*(T-Tfluid); coolant is a finite internal inventory.
        let g = PATCH_AREA * 0.664 * 100.0 * 2.0 * 0.5 / LENGTH;
        let fluid_capacity = FLUID_MASS * FLUID_CP;
        let w = cp.map(|capacity| dt * g * capacity / (capacity + dt * g));
        let fluid = (fluid_capacity * 290.0
            + (0..10)
                .map(|i| 300.0 * w[i] + dt * g * q[i] / (cp[i] + dt * g))
                .sum::<f64>())
            / (fluid_capacity + w.iter().sum::<f64>());
        close(next.thermal().temperature_k(COOLANT).unwrap(), fluid);
        for i in 0..10 {
            let t = (cp[i] * 300.0 + q[i] + dt * g * fluid) / (cp[i] + dt * g);
            close(next.thermal().temperature_k(DRIVE_IDS[i]).unwrap(), t);
        }
        close(
            next.thermal().energy_j() - thermal.energy_j(),
            q.iter().sum(),
        );
        close(next.report().deposited_joule_j, q.iter().sum());
        close(next.report().non_electrical_heat_j, 0.0);
        close(next.report().total_energy_residual_j, 0.0);
        assert!(fluid > 290.0);
    });
}
