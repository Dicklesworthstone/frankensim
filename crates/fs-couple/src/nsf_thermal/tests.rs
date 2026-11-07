use super::*;

mod interval;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{NsfLaw, NsfSpec, WallThermal};
use fs_material::gas::GasSpec;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, Event, MaterialId, StepConfig,
    ThermalLink, TimedEvent,
};

const DX: f64 = 1e-4;
const CP: f64 = 1000.0;
const FACE_AREA: f64 = DX * DX;
const PLANE_AREA: f64 = 4.0 * FACE_AREA;

struct Law {
    datum: f64,
}
impl ThermalLaw for Law {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let temperature = 300.0 + (h - self.datum) / CP;
        if !(273.15..=333.15).contains(&temperature) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(temperature)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / CP
    }
}
struct K(f64);
impl ConductivityLaw for K {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(self.0)
    }
}

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 0x4e53_4654,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn fluid() -> NsfGrid3<D3q27> {
    NsfGrid3::new(
        NsfSpec {
            dims: [2, 1, 2],
            periodic: [true, false, true],
            gas: GasSpec::dry_air_ussa1976(),
            temperature_k: 300.0,
            pressure_pa: 100_000.0,
            dx_m: DX,
            temperature_window_k: [273.15, 333.15],
            maximum_mach: 0.2,
            maximum_temperature_ratio: 0.2,
            bulk_to_shear_ratio: 0.73,
            bulk_ratio_source: "Cramer 2012 N2 ratio; finite wall exchange fixture",
            maximum_cells: 4,
            wall: WallThermal::Adiabatic,
            workers: 1,
            transport_points: 17,
        },
        |_| false,
        |_| (1.0, [0.0; 3], 300.0),
    )
    .unwrap()
}

fn solid<'a>(
    law: &'a Law,
    values: &[f64],
    mass: f64,
    ambient: bool,
    linked: bool,
) -> DepositionSimulation<'a> {
    let mut events: Vec<_> = values
        .iter()
        .enumerate()
        .map(|(id, &temperature)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id as u64),
                segment_id: id as u64,
                body_id: 1,
                material_id: MaterialId(7),
                mass_kg: mass,
                specific_h_j_kg: law.datum + CP * (temperature - 300.0),
                geometry: CellGeometry {
                    start_m: [id as f64 * DX, 0.0, 0.0],
                    end_m: [(id + 1) as f64 * DX, 0.0, 0.0],
                    volume_m3: 2.0 * PLANE_AREA * 1e-5,
                    surface_area_m2: 2.0 * PLANE_AREA,
                },
                boundary: Boundary {
                    ambient_k: 320.0,
                    htc_w_m2_k: if ambient { 25.0 } else { 0.0 },
                    reservoir: None,
                    source_w: if ambient { 3e-5 } else { 0.0 },
                },
            }),
        })
        .collect();
    if linked {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 19,
                a: CellId(0),
                b: CellId(1),
                area_m2: 0.1 * PLANE_AREA,
                conductance_w_k: 1e-6,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(7), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 2,
            max_events: 8,
            max_substeps: 10_000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn contacts<'a>(
    law: &'a Law,
    layer: NsfSolidLayer<'a>,
    separate: bool,
) -> Vec<NsfSolidContact<'a>> {
    (0..2)
        .flat_map(|x| {
            (0..2).flat_map(move |z| {
                [false, true].map(move |high| NsfSolidContact {
                    fluid_cell: [x, 0, z],
                    axis: 1,
                    high,
                    solid_cell: CellId(u64::from(separate && high)),
                    caloric: law,
                    layer,
                })
            })
        })
        .collect()
}

fn controls() -> NsfThermalControls {
    NsfThermalControls {
        maximum_fluid_cells: 4,
        maximum_solid_cells: 2,
        maximum_contacts: 8,
        maximum_endpoint_evaluations: 256,
        maximum_newton_iterations: 16,
        energy_scale_j: 1e-9,
        finite_difference_step: 1e-4,
        cell_tolerance_j: 1e-16,
        thermal_tolerance_j: 1e-16,
        fluid_tolerance_j: 1e-17,
        interface_tolerance_j: 1e-20,
        total_tolerance_j: 2e-16,
        solver: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
    }
}

fn unforced(
    _: &Cx<'_>,
    _: &DepositionSimulation<'_>,
) -> Result<BoundaryStepCoefficients<()>, DepositionError> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    }
    .into())
}

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

// Independent scalar series-resistance and two-finite-capacity equations.
// The only reused numerical data are the original material's k(T) values.
fn physical_k(law: &NsfLaw, temperature: f64) -> f64 {
    law.kappa_lattice(law.theta_of_kelvin(temperature))
        * law.gas_constant
        * law.rho_ref_kg_m3
        * law.dx_m.powi(2)
        / law.dt_s
}
fn scalar_flux(law: &NsfLaw, fluid_k: f64, solid_k: f64, resistance: f64) -> f64 {
    if resistance == 0.0 {
        return physical_k(law, 0.5 * (fluid_k + solid_k)) * (solid_k - fluid_k) / (0.5 * DX);
    }
    let (mut lo, mut hi) = (fluid_k.min(solid_k), fluid_k.max(solid_k));
    for _ in 0..80 {
        let wall = f64::midpoint(lo, hi);
        let into_fluid = physical_k(law, 0.5 * (fluid_k + wall)) * (wall - fluid_k) / (0.5 * DX);
        let from_solid = (solid_k - wall) / resistance;
        if into_fluid > from_solid {
            hi = wall;
        } else {
            lo = wall;
        }
    }
    (solid_k - f64::midpoint(lo, hi)) / resistance
}

#[test]
fn finite_wall_matches_independent_two_capacity_surface_balance() {
    with_cx(|cx, _| {
        let law = Law { datum: 1e8 };
        let conductivity = K(0.1);
        let mut gas = fluid();
        let gas_capacity = gas.total_mass_kg() * gas.law().gas_constant / (gas.law().gamma - 1.0);
        let wall_capacity = 2.0 * gas_capacity;
        let mass = wall_capacity / CP;
        let mut wall = solid(&law, &[310.0], mass, false, false);
        let initial_h = wall.cells()[&CellId(0)].specific_h_j_kg;
        let mapping = contacts(
            &law,
            NsfSolidLayer::Conducting {
                half_distance_m: 1e-5,
                contact_resistance_m2_k_w: 2e-5,
                conductivity: &conductivity,
            },
            false,
        );
        let resistance = 1e-5 / 0.1 + 2e-5;
        let (mut low, mut high) = (0.0, 10.0 / (1.0 / gas_capacity + 1.0 / wall_capacity));
        for _ in 0..90 {
            let heat = f64::midpoint(low, high);
            let t_wall = 310.0 - heat / wall_capacity;
            let t_gas = 300.0 + heat / gas_capacity;
            let transferred = 0.5
                * gas.law().dt_s
                * (2.0 * PLANE_AREA)
                * (scalar_flux(gas.law(), 300.0, t_wall, resistance)
                    + scalar_flux(gas.law(), t_gas, t_wall, resistance));
            if heat > transferred {
                high = heat;
            } else {
                low = heat;
            }
        }
        let expected_heat = f64::midpoint(low, high);
        let report = advance_nsf_solid(
            cx,
            &mut gas,
            &mut wall,
            &mapping,
            controls(),
            &[0.0],
            unforced,
        )
        .unwrap();
        near(report.coupled_heat_into_fluid_j, expected_heat, 1e-18, 2e-6);
        near(
            wall.temperature_k(CellId(0)).unwrap(),
            310.0 - expected_heat / wall_capacity,
            2e-7,
            0.0,
        );
        near(
            gas.cell([0, 0, 0]).2,
            300.0 + expected_heat / gas_capacity,
            2e-7,
            0.0,
        );
        near(
            mass * (wall.cells()[&CellId(0)].specific_h_j_kg - initial_h),
            -expected_heat,
            2e-16,
            2e-6,
        );
        assert!(
            report
                .exchanges
                .iter()
                .all(|entry| entry.face.heat_into_fluid_j > 0.0)
        );
        assert_eq!(report.exchanges.len(), 8);
        assert_eq!(report.solid_external_heat_j, 0.0);
        assert!(report.total_energy_residual_j.abs() < controls().total_tolerance_j);
        assert_eq!(gas.steps(), 1);
        assert_eq!(wall.time_s(), gas.law().dt_s);
        // The next call keeps the original histories, law and common clock.
        advance_nsf_solid(
            cx,
            &mut gas,
            &mut wall,
            &mapping,
            controls(),
            &[0.0],
            unforced,
        )
        .unwrap();
        assert_eq!(gas.steps(), 2);
        assert_eq!(wall.time_s(), 2.0 * gas.law().dt_s);
    });
}

#[test]
fn opposing_wall_heat_keeps_spatial_solid_conduction_and_remaining_ambient_area() {
    with_cx(|cx, _| {
        let law = Law { datum: 0.0 };
        let mut gas = fluid();
        let mass = 1e-11;
        let mut wall = solid(&law, &[310.0, 290.0], mass, true, true);
        let old_h = [
            wall.cells()[&CellId(0)].specific_h_j_kg,
            wall.cells()[&CellId(1)].specific_h_j_kg,
        ];
        let mapping = contacts(&law, NsfSolidLayer::Isothermal, true);
        let report = advance_nsf_solid(
            cx,
            &mut gas,
            &mut wall,
            &mapping,
            controls(),
            &[0.0; 2],
            |_, _| {
                Ok::<_, DepositionError>(
                    StepCoefficients {
                        sources_w: [(CellId(1), 5e-5)].into(),
                        conductances_w_k: BTreeMap::new(),
                        payload: 42,
                    }
                    .into(),
                )
            },
        )
        .unwrap();
        let dt = gas.law().dt_s;
        let t = [
            wall.temperature_k(CellId(0)).unwrap(),
            wall.temperature_k(CellId(1)).unwrap(),
        ];
        let mut q = [0.0; 2];
        for entry in &report.exchanges {
            q[entry.solid_cell.0 as usize] += entry.face.heat_into_fluid_j;
        }
        assert!(q[0] > 0.0 && q[1] < 0.0);
        let conduction = dt * 1e-6 * (t[0] - t[1]);
        let remaining = 2.0 * PLANE_AREA - 0.1 * PLANE_AREA - PLANE_AREA;
        let external = [
            dt * (remaining * 25.0 * (320.0 - t[0]) + 3e-5),
            dt * (remaining * 25.0 * (320.0 - t[1]) + 3e-5 + 5e-5),
        ];
        for i in 0..2 {
            let stored = mass * (wall.cells()[&CellId(i as u64)].specific_h_j_kg - old_h[i]);
            let expected = external[i] - q[i] + if i == 0 { -conduction } else { conduction };
            near(stored, expected, 1e-16, 0.0);
        }
        near(
            report.solid_external_heat_j,
            external.iter().sum(),
            1e-24,
            1e-12,
        );
        near(
            report.thermal.stored_energy_change_j + report.coupled_heat_into_fluid_j,
            external.iter().sum(),
            2e-16,
            0.0,
        );
        assert_eq!(report.heat, 42);
    });
}

#[test]
fn source_insulation_and_finite_resistance_change_the_original_heat() {
    with_cx(|cx, _| {
        let law = Law { datum: 0.0 };
        let mut heat = Vec::new();
        for k in [0.0, 0.01, 10.0] {
            let conductivity = K(k);
            let mut gas = fluid();
            let mut wall = solid(&law, &[310.0], 1e-11, false, false);
            let mapping = contacts(
                &law,
                NsfSolidLayer::Conducting {
                    half_distance_m: 1e-5,
                    contact_resistance_m2_k_w: 0.0,
                    conductivity: &conductivity,
                },
                false,
            );
            let report = advance_nsf_solid(
                cx,
                &mut gas,
                &mut wall,
                &mapping,
                controls(),
                &[0.0],
                unforced,
            )
            .unwrap();
            heat.push(report.coupled_heat_into_fluid_j);
            if k == 0.0 {
                assert_eq!(wall.temperature_k(CellId(0)).unwrap(), 310.0);
                assert!(
                    report
                        .exchanges
                        .iter()
                        .all(|entry| entry.face.heat_into_fluid_j == 0.0)
                );
            }
        }
        assert_eq!(heat[0], 0.0);
        assert!(heat[1] > 0.0 && heat[2] > 1.2 * heat[1], "heat={heat:?}");
    });
}

#[test]
fn invalid_bindings_and_exhausted_work_publish_neither_owner() {
    with_cx(|cx, _| {
        let law = Law { datum: 0.0 };
        let other_law = Law { datum: 0.0 };
        let mut gas = fluid();
        let mut wall = solid(&law, &[310.0], 1e-11, false, false);
        let mapping = contacts(&law, NsfSolidLayer::Isothermal, false);
        let before = (gas.state_hash(), wall.cells().clone(), wall.time_s());
        let mut duplicate = mapping.clone();
        duplicate[1] = duplicate[0];
        assert!(
            advance_nsf_solid(
                cx,
                &mut gas,
                &mut wall,
                &duplicate,
                controls(),
                &[0.0],
                unforced
            )
            .is_err()
        );
        let mut wrong = mapping.clone();
        wrong[0].caloric = &other_law;
        assert!(
            advance_nsf_solid(
                cx,
                &mut gas,
                &mut wall,
                &wrong,
                controls(),
                &[0.0],
                unforced
            )
            .is_err()
        );
        let mut budget = controls();
        budget.maximum_endpoint_evaluations = 1;
        assert!(matches!(
            advance_nsf_solid(cx, &mut gas, &mut wall, &mapping, budget, &[0.0], unforced),
            Err(NsfThermalError::WorkBudget { performed: 1, .. })
        ));
        assert_eq!(
            before,
            (gas.state_hash(), wall.cells().clone(), wall.time_s())
        );
        // A numerical source failure after an actual native candidate also rolls back.
        let calls = Cell::new(0);
        let error = advance_nsf_solid(
            cx,
            &mut gas,
            &mut wall,
            &mapping,
            controls(),
            &[0.0],
            |cx, candidate| {
                calls.set(calls.get() + 1);
                if calls.get() == 2 {
                    Err(DepositionError::MaterialDomain)
                } else {
                    unforced(cx, candidate)
                }
            },
        )
        .unwrap_err();
        assert!(matches!(
            error,
            NsfThermalError::Resolve(DepositionError::MaterialDomain)
        ));
        assert_eq!(calls.get(), 2);
        assert_eq!(
            before,
            (gas.state_hash(), wall.cells().clone(), wall.time_s())
        );
        advance_nsf_solid(
            cx,
            &mut gas,
            &mut wall,
            &mapping,
            controls(),
            &[0.0],
            unforced,
        )
        .unwrap();
    });
}

#[test]
fn overlapping_ambient_owner_and_late_cancellation_are_transactional() {
    with_cx(|cx, gate| {
        let law = Law { datum: 0.0 };
        let mut gas = fluid();
        let mut wall = solid(&law, &[310.0], 1e-11, false, false);
        let mapping = contacts(&law, NsfSolidLayer::Isothermal, false);
        let before = (gas.state_hash(), wall.cells().clone(), wall.time_s());
        let error = advance_nsf_solid(
            cx,
            &mut gas,
            &mut wall,
            &mapping,
            controls(),
            &[0.0],
            |cx, candidate| {
                let mut heat = unforced(cx, candidate)?;
                heat.boundaries.ambient_w_k.insert(CellId(0), 0.0);
                Ok::<_, DepositionError>(heat)
            },
        )
        .unwrap_err();
        assert!(matches!(error, NsfThermalError::Invalid(_)));
        let calls = Cell::new(0);
        let error = advance_nsf_solid(
            cx,
            &mut gas,
            &mut wall,
            &mapping,
            controls(),
            &[0.0],
            |cx, candidate| {
                calls.set(calls.get() + 1);
                if calls.get() == 2 {
                    gate.request();
                }
                unforced(cx, candidate)
            },
        )
        .unwrap_err();
        assert!(matches!(error, NsfThermalError::Cancelled));
        assert_eq!(
            before,
            (gas.state_hash(), wall.cells().clone(), wall.time_s())
        );
    });
}

#[test]
fn dynamic_channel_transports_heat_with_finite_solid_storage_and_one_clock() {
    use fs_lbm::nsf3::NsfChannelDynamicThermalSgsConfig;
    use fs_lbm::sgs3::SgsUse;
    use fs_lbm::sgs3::dynamic::channel::{ChannelDynamicScalarConfig, ChannelDynamicSgsConfig};

    fn noise(index: usize, salt: u64) -> f64 {
        let mut value = (index as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ salt;
        value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        value ^= value >> 31;
        ((value >> 11) as f64 / (1_u64 << 53) as f64) * 2.0 - 1.0
    }
    let dynamic = NsfChannelDynamicThermalSgsConfig {
        momentum: ChannelDynamicSgsConfig {
            maximum_coefficient_squared: 16.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        heat: ChannelDynamicScalarConfig {
            maximum_coefficient: 16.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    };
    let make = || {
        NsfGrid3::<D3q27>::new(
            NsfSpec {
                dims: [5, 4, 5],
                periodic: [true, true, false],
                gas: GasSpec::dry_air_ussa1976(),
                temperature_k: 300.0,
                pressure_pa: 100_000.0,
                dx_m: DX,
                temperature_window_k: [273.15, 333.15],
                maximum_mach: 0.2,
                maximum_temperature_ratio: 0.2,
                bulk_to_shear_ratio: 0.73,
                bulk_ratio_source: "Cramer N2 ratio; finite solid/dynamic channel fixture",
                maximum_cells: 100,
                wall: WallThermal::Adiabatic,
                workers: 1,
                transport_points: 17,
            },
            |_| false,
            |p| {
                let index = (p[2] * 4 + p[1]) * 5 + p[0];
                let height = p[2] as f64 + 0.5;
                let envelope = 4.0 * height * (5.0 - height) / 25.0;
                (
                    1.0,
                    std::array::from_fn(|axis| {
                        0.002 * envelope * noise(index, 31 * axis as u64 + 7)
                    }),
                    300.0 + 0.1 * (1.0 - height / 5.0 + 0.05 * envelope * noise(index, 101)),
                )
            },
        )
        .unwrap()
        .with_channel_dynamic_thermal_sgs(dynamic, 20_000)
        .unwrap()
    };
    with_cx(|cx, _| {
        let law = Law { datum: 0.0 };
        let mass = 1e-11;
        let mut gas = make();
        let mut direct = make();
        assert!((0..5).any(|z| (0..4).any(|y| (0..5).any(|x| {
            gas.sgs_transport([x, y, z])
                .unwrap()
                .unwrap()
                .turbulent_conductivity_lu
                > 0.0
        }))));
        let mut wall = solid(&law, &[310.0, 290.0], mass, false, false);
        let mapping =
            [([1, 1, 0], false), ([2, 2, 4], true)].map(|(fluid_cell, high)| NsfSolidContact {
                fluid_cell,
                axis: 2,
                high,
                solid_cell: CellId(u64::from(high)),
                caloric: &law,
                layer: NsfSolidLayer::Isothermal,
            });
        let mut budget = controls();
        budget.maximum_fluid_cells = 100;
        budget.maximum_contacts = 2;
        for _ in 0..2 {
            let before = gas.total_energy_j();
            let h = [
                wall.cells()[&CellId(0)].specific_h_j_kg,
                wall.cells()[&CellId(1)].specific_h_j_kg,
            ];
            let report = advance_nsf_solid(
                cx, &mut gas, &mut wall, &mapping, budget, &[0.0; 2], unforced,
            )
            .unwrap();
            let mut heat = [0.0; 2];
            for exchange in &report.exchanges {
                heat[exchange.solid_cell.0 as usize] += exchange.face.heat_into_fluid_j;
            }
            assert!(heat[0] > 0.0 && heat[1] < 0.0);
            let mut stored = 0.0;
            for (index, old) in h.into_iter().enumerate() {
                let change = mass * (wall.cells()[&CellId(index as u64)].specific_h_j_kg - old);
                near(change, -heat[index], budget.cell_tolerance_j, 0.0);
                stored += change;
            }
            near(
                gas.total_energy_j() - before + stored,
                0.0,
                budget.total_tolerance_j,
                0.0,
            );
            near(
                report.coupled_heat_into_fluid_j,
                report.fluid.wall_heat_j,
                1e-22,
                1e-12,
            );
            assert_eq!(report.solid_external_heat_j, 0.0);
            assert_eq!(gas.channel_dynamic_thermal_sgs_config(), Some(dynamic));
            assert_eq!(gas.dynamic_sgs_workspace_bytes(), Some(20_000));
            assert_eq!(wall.time_s(), gas.steps() as f64 * gas.law().dt_s);

            // The original dynamic fluid, stepped independently under the
            // actually accepted solid endpoint, must give the identical state.
            direct
                .set_thermal_walls(
                    mapping
                        .iter()
                        .map(|contact| NsfThermalWallFace {
                            cell: contact.fluid_cell,
                            axis: contact.axis,
                            high: contact.high,
                            thermal: NsfWallThermal::Isothermal(
                                wall.temperature_k(contact.solid_cell).unwrap(),
                            ),
                        })
                        .collect(),
                )
                .unwrap();
            assert_eq!(direct.step().unwrap(), report.fluid);
            assert_eq!(gas.state_hash(), direct.state_hash());
        }
        let accepted = (gas.state_hash(), wall.cells().clone(), wall.time_s());
        budget.maximum_endpoint_evaluations = 1;
        assert!(matches!(
            advance_nsf_solid(
                cx, &mut gas, &mut wall, &mapping, budget, &[0.0; 2], unforced
            ),
            Err(NsfThermalError::WorkBudget { performed: 1, .. })
        ));
        assert_eq!(
            (gas.state_hash(), wall.cells().clone(), wall.time_s()),
            accepted
        );
        assert_eq!(gas.channel_dynamic_thermal_sgs_config(), Some(dynamic));
    });
}
