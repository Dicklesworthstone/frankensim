//! Native storage, independent acoustic flux, coherent phase and spatial refinement.
//! Authored numerical material fixtures are not measured heating validation.
use core::cell::Cell;
use fs_couple::porous_heating::PorousHeatCell;
use fs_couple::porous_heating::coherent::*;
use fs_material::MaterialError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::porous::layers::{
    AngularSpectrumBudget, IncidentPlaneWave, Layer, LayerAbsorptionBudget,
    PreparedLayerAcousticSpectrum, TransmittingLayers,
};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_math::c64::C64;
use fs_thermal::enthalpy::*;

const LENGTH: f64 = 0.06;
const OMEGA: f64 = core::f64::consts::TAU * 1000.0;
const CP: f64 = 200.0;
const DENSITY: f64 = 100.0;
const MID: [AcousticHeatNode; 1] = [AcousticHeatNode {
    coordinate: 0.5,
    weight: 1.0,
}];

struct Caloric(u8);
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h / CP;
        if t.is_finite() && t > 0.0 {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / CP
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([self.0; 32])
    }
}
static CALORIC: Caloric = Caloric(31);
static OTHER_CALORIC: Caloric = Caloric(32);

struct Source {
    state: ThermoviscousFluidState,
    samples: Cell<usize>,
}
impl FluidAcoustics for Source {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.samples.set(self.samples.get() + 1);
        Ok(self.state.acoustic())
    }
}
impl FluidTransport for Source {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.state.transport())
    }
}
impl ThermoviscousFluid for Source {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        self.samples.set(self.samples.get() + 1);
        Ok(self.state)
    }
}
fn source() -> Source {
    Source {
        state: ThermoviscousFluidState::from_heat_capacity_ratio(
            FluidAcousticState::try_new(1.2, 340.0).unwrap(),
            FluidTransportState::try_new(1.2, 1005.0, 0.0257, 1.81e-5).unwrap(),
            1.4,
        )
        .unwrap(),
        samples: Cell::new(0),
    }
}
fn material() -> PorousFluidLaw {
    PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap()
}
fn column() -> CoherentThermalColumn {
    CoherentThermalColumn {
        origin_m: [0.2, -0.3, 0.1],
        transverse_axes: [[0.0, 0.0, 1.0], [0.8, -0.6, 0.0]],
        lower_m: [-0.01, -0.02],
        extent_m: [0.02, 0.04],
    }
}
fn point(c: CoherentThermalColumn, z: f64) -> [f64; 3] {
    let x = c.lower_m[0] + 0.5 * c.extent_m[0];
    let y = c.lower_m[1] + 0.5 * c.extent_m[1];
    // Independent explicit fixture frame: normal is [0.6,0.8,0].
    [
        c.origin_m[0] + 0.8 * y + 0.6 * z,
        c.origin_m[1] - 0.6 * y + 0.8 * z,
        c.origin_m[2] + x,
    ]
}
fn id(i: usize) -> CellId {
    CellId(900 - 7 * i as u64)
}
fn partition(n: usize) -> Vec<PorousHeatCell> {
    (0..n)
        .map(|i| PorousHeatCell {
            cell: id(i),
            end_depth_m: if i + 1 == n {
                LENGTH
            } else {
                LENGTH * (i + 1) as f64 / n as f64
            },
        })
        .collect()
}
fn config() -> StepConfig {
    StepConfig {
        max_dt_s: 0.07,
        energy_tolerance_j: 1e-10,
    }
}
fn initial(
    c: CoherentThermalColumn,
    routes: &[PorousHeatCell],
    volume_scale: f64,
    mut extra: Vec<TimedEvent>,
    law: &'static dyn ThermalLaw,
) -> DepositionSimulation<'static> {
    let area = c.extent_m[0] * c.extent_m[1];
    let mut start = 0.0;
    let mut events = Vec::new();
    for (i, route) in routes.iter().enumerate() {
        let width = route.end_depth_m - start;
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: route.cell,
                segment_id: i as u64,
                body_id: 70,
                material_id: MaterialId(5),
                mass_kg: DENSITY * area * width,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: point(c, start),
                    end_m: point(c, route.end_depth_m),
                    volume_m3: volume_scale * area * width,
                    surface_area_m2: 2.0 * area + 2.0 * (c.extent_m[0] + c.extent_m[1]) * width,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
        start = route.end_depth_m;
    }
    events.append(&mut extra);
    DepositionSimulation::new(
        events,
        &[(MaterialId(5), law)],
        DepositionBudget {
            max_cells: 32,
            max_links: 32,
            max_events: 64,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}
fn budget() -> CoherentHeatingBudget {
    CoherentHeatingBudget {
        max_cells: 32,
        max_component_evaluations: 1_000_000,
    }
}
fn midpoint(n: usize) -> Vec<AcousticHeatNode> {
    (0..n)
        .map(|i| AcousticHeatNode {
            coordinate: (i as f64 + 0.5) / n as f64,
            weight: 1.0 / n as f64,
        })
        .collect()
}
fn simpson(n: usize) -> Vec<AcousticHeatNode> {
    assert_eq!(n % 2, 0);
    (0..=n)
        .map(|i| AcousticHeatNode {
            coordinate: i as f64 / n as f64,
            weight: if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            } / (3.0 * n as f64),
        })
        .collect()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:.3e}"
    );
}
fn total(heat: &CoherentLayerHeating<'_, '_, '_, '_>) -> f64 {
    heat.cells()
        .iter()
        .map(|cell| cell.dissipation.absorbed_power_w)
        .sum()
}

#[test]
fn original_volume_heat_matches_independent_face_work_and_native_storage() {
    let fluid = source();
    let law = material();
    let stack = TransmittingLayers::new(
        &fluid,
        &fluid,
        OMEGA,
        &[Layer::Porous {
            law: &law,
            thickness_m: LENGTH,
        }],
    )
    .unwrap();
    let source_queries = fluid.samples.get();
    let wave = IncidentPlaneWave {
        direction: [0.36, 0.48, 0.8],
        pressure_pa: C64::new(30.0, -7.0),
    };
    let beam = stack
        .prepare_angular_spectrum(
            &[wave],
            AngularSpectrumBudget { max_components: 1 },
            &mut || false,
        )
        .unwrap();
    let field = beam.prepare_layer_field(0, &mut || false).unwrap();
    let c = column();
    let area = c.extent_m[0] * c.extent_m[1];
    let incident = area * 0.5 * wave.pressure_pa.norm_sq() * 0.8 / (1.2 * 340.0);
    let absorption = stack
        .layer_absorption(
            0.8,
            incident,
            LayerAbsorptionBudget { max_layers: 1 },
            || false,
        )
        .unwrap();
    let depth = simpson(128);
    for n in [1, 3, 11] {
        let routes = partition(n);
        let thermal = initial(c, &routes, 1.0, vec![], &CALORIC);
        let heat = CoherentLayerHeating::try_new(
            &field,
            &thermal,
            c,
            &routes,
            CoherentHeatQuadrature {
                axes: [&MID, &MID, &depth],
            },
            budget(),
            || Ok(()),
        )
        .unwrap();
        assert!(core::ptr::eq(heat.field(), &field));
        close(total(&heat), absorption.absorbed_power_w, 2e-9 * incident);
        for cell in heat.cells() {
            // For one plane wave the transverse intensity is spatially constant.
            // Actual acoustic normal work, independently of local heat quadrature.
            let front = field
                .field([0.0, 0.0, cell.start_depth_m], &mut || false)
                .unwrap();
            let rear = field
                .field([0.0, 0.0, cell.end_depth_m], &mut || false)
                .unwrap();
            let work = area * (front.active_intensity_w_m2[2] - rear.active_intensity_w_m2[2]);
            close(cell.dissipation.absorbed_power_w, work, 2e-9 * incident);
        }
        let step = heat.advance_to(&thermal, 1.3, config(), || Ok(())).unwrap();
        let expected = 1.3 * total(&heat);
        close(step.acoustic.deposited_j, expected, 3e-15 * expected);
        close(
            step.acoustic.viscous_j + step.acoustic.thermal_j,
            expected,
            3e-15 * expected,
        );
        close(step.heat_report.external_heat_j, expected, 3e-15 * expected);
        close(step.thermal.energy_j(), expected, 3e-15 * expected);
        for (bound, (cell, joules)) in heat.cells().iter().zip(&step.cell_energy_j) {
            assert_eq!(bound.cell, *cell);
            let original = &thermal.cells()[cell];
            close(
                *joules,
                1.3 * bound.dissipation.absorbed_power_w,
                3e-15 * expected,
            );
            close(
                step.thermal.temperature_k(*cell).unwrap(),
                300.0 + joules / (CP * original.mass_kg),
                6e-13,
            );
            assert_eq!(step.thermal.cells()[cell].geometry, original.geometry);
            assert_eq!(step.thermal.cells()[cell].boundary.source_w, 0.0);
        }
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.energy_j(), 0.0);
    }
    assert_eq!(
        fluid.samples.get(),
        source_queries,
        "warm acoustic/thermal work must retain frozen sources"
    );
}

#[test]
fn native_volume_and_footprint_scale_watts_without_changing_uniform_temperature_rise() {
    let fluid = source();
    let law = material();
    let stack = TransmittingLayers::new(
        &fluid,
        &fluid,
        OMEGA,
        &[Layer::Porous {
            law: &law,
            thickness_m: LENGTH,
        }],
    )
    .unwrap();
    let beam = stack
        .prepare_angular_spectrum(
            &[IncidentPlaneWave {
                direction: [0.0, 0.0, 1.0],
                pressure_pa: C64::from_re(30.0),
            }],
            AngularSpectrumBudget { max_components: 1 },
            &mut || false,
        )
        .unwrap();
    let field = beam.prepare_layer_field(0, &mut || false).unwrap();
    let depth = simpson(128);
    let routes = partition(1);
    let mut previous = None;
    for factor in [1.0, 3.0] {
        let mut c = column();
        c.extent_m[0] *= factor;
        let thermal = initial(c, &routes, 1.0, vec![], &CALORIC);
        let heat = CoherentLayerHeating::try_new(
            &field,
            &thermal,
            c,
            &routes,
            CoherentHeatQuadrature {
                axes: [&MID, &MID, &depth],
            },
            budget(),
            || Ok(()),
        )
        .unwrap();
        let step = heat.advance_to(&thermal, 1.0, config(), || Ok(())).unwrap();
        let watts = total(&heat);
        let temperature = step.thermal.temperature_k(id(0)).unwrap();
        if let Some((first_watts, first_temperature)) = previous {
            close(watts, 3.0 * first_watts, 2e-14 * watts);
            close(temperature, first_temperature, 6e-13);
        }
        previous = Some((watts, temperature));
    }
}

// Independent closed transverse integral for symmetric crossing waves. The
// original single-wave fields supply depth dependence; analytical cos²/sin²
// integrals retain transverse velocity rather than treating pressure nodes as
// total silence. Simpson depth integration is refined far beyond production.
fn crossing_reference(
    single: &PreparedLayerAcousticSpectrum<'_, '_, '_>,
    c: CoherentThermalColumn,
    phase: f64,
) -> [f64; 2] {
    let response = single
        .spectrum()
        .stack()
        .prepared_layers()
        .material_response(0)
        .unwrap();
    let weights = [
        -0.5 * OMEGA * response.dynamic_density_kg_m3().im,
        -0.5 * OMEGA * response.dynamic_bulk_modulus_pa().recip().im,
    ];
    let kx = OMEGA * 0.6 / 340.0;
    let a = c.lower_m[0];
    let b = a + c.extent_m[0];
    let cosine =
        0.5 + ((2.0 * kx * b + phase).sin() - (2.0 * kx * a + phase).sin()) / (4.0 * kx * (b - a));
    let mut result = [0.0; 2];
    for node in simpson(1024) {
        let sample = single
            .field([0.0, 0.0, LENGTH * node.coordinate], &mut || false)
            .unwrap();
        result[0] += 4.0
            * node.weight
            * weights[0]
            * (cosine * sample.particle_velocity_m_s[2].norm_sq()
                + (1.0 - cosine) * sample.particle_velocity_m_s[0].norm_sq());
        result[1] += 4.0 * node.weight * weights[1] * cosine * sample.pressure_pa.norm_sq();
    }
    result.map(|q| q * c.extent_m[0] * c.extent_m[1] * LENGTH)
}

#[test]
fn coherent_cross_angle_phase_and_positive_spatial_refinement_match_analytic_window() {
    let fluid = source();
    let law = material();
    let stack = TransmittingLayers::new(
        &fluid,
        &fluid,
        OMEGA,
        &[Layer::Porous {
            law: &law,
            thickness_m: LENGTH,
        }],
    )
    .unwrap();
    let first = IncidentPlaneWave {
        direction: [0.6, 0.0, 0.8],
        pressure_pa: C64::from_re(30.0),
    };
    let one = stack
        .prepare_angular_spectrum(
            &[first],
            AngularSpectrumBudget { max_components: 1 },
            &mut || false,
        )
        .unwrap();
    let single = one.prepare_layer_field(0, &mut || false).unwrap();
    let mut c = column();
    c.lower_m[0] = -0.02;
    c.extent_m[0] = 0.15;
    let routes = partition(1);
    let depth = simpson(128);
    let thermal = initial(c, &routes, 1.0, vec![], &CALORIC);
    let mut phase_powers = Vec::new();
    for phase in [0.0, 1.3, core::f64::consts::PI] {
        let beam = stack
            .prepare_angular_spectrum(
                &[
                    first,
                    IncidentPlaneWave {
                        direction: [-0.6, 0.0, 0.8],
                        pressure_pa: C64::new(phase.cos(), phase.sin()).scale(30.0),
                    },
                ],
                AngularSpectrumBudget { max_components: 2 },
                &mut || false,
            )
            .unwrap();
        let field = beam.prepare_layer_field(0, &mut || false).unwrap();
        let expected = crossing_reference(&single, c, phase);
        let mut prior_error = f64::INFINITY;
        for resolution in [2, 8, 32, 128] {
            let lateral = midpoint(resolution);
            let heat = CoherentLayerHeating::try_new(
                &field,
                &thermal,
                c,
                &routes,
                CoherentHeatQuadrature {
                    axes: [&lateral, &MID, &depth],
                },
                budget(),
                || Ok(()),
            )
            .unwrap();
            let actual = heat.cells()[0].dissipation;
            assert!(actual.viscous_power_w >= 0.0 && actual.thermal_power_w >= 0.0);
            let error = (actual.viscous_power_w - expected[0]).abs()
                + (actual.thermal_power_w - expected[1]).abs();
            assert!(
                error < prior_error / 8.0,
                "phase {phase}, resolution {resolution}, {error} >= {prior_error}/8"
            );
            prior_error = error;
            if resolution == 128 {
                close(
                    actual.viscous_power_w,
                    expected[0],
                    3e-5 * (expected[0] + expected[1]),
                );
                close(
                    actual.thermal_power_w,
                    expected[1],
                    3e-5 * (expected[0] + expected[1]),
                );
                phase_powers.push(actual.absorbed_power_w);
            }
        }
    }
    assert!((phase_powers[0] - phase_powers[2]).abs() > 0.1 * phase_powers[0]);
}

#[test]
fn coherent_duplicate_components_cancel_before_native_heat_is_credited() {
    let fluid = source();
    let law = material();
    let stack = TransmittingLayers::new(
        &fluid,
        &fluid,
        OMEGA,
        &[Layer::Porous {
            law: &law,
            thickness_m: LENGTH,
        }],
    )
    .unwrap();
    let wave = IncidentPlaneWave {
        direction: [0.6, 0.0, 0.8],
        pressure_pa: C64::new(30.0, 7.0),
    };
    let beam = stack
        .prepare_angular_spectrum(
            &[
                wave,
                IncidentPlaneWave {
                    pressure_pa: -wave.pressure_pa,
                    ..wave
                },
            ],
            AngularSpectrumBudget { max_components: 2 },
            &mut || false,
        )
        .unwrap();
    let field = beam.prepare_layer_field(0, &mut || false).unwrap();
    let routes = partition(3);
    let c = column();
    let nodes = midpoint(4);
    let thermal = initial(c, &routes, 1.0, vec![], &CALORIC);
    let heat = CoherentLayerHeating::try_new(
        &field,
        &thermal,
        c,
        &routes,
        CoherentHeatQuadrature {
            axes: [&nodes, &nodes, &nodes],
        },
        budget(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(total(&heat), 0.0);
    let step = heat.advance_to(&thermal, 1.0, config(), || Ok(())).unwrap();
    assert_eq!(step.acoustic.deposited_j, 0.0);
    assert_eq!(step.acoustic.viscous_j, 0.0);
    assert_eq!(step.acoustic.thermal_j, 0.0);
    assert_eq!(step.thermal.energy_j(), 0.0);
}

#[test]
fn native_events_contacts_and_authored_heat_keep_single_credit_and_cached_source() {
    let fluid = source();
    let law = material();
    let stack = TransmittingLayers::new(
        &fluid,
        &fluid,
        OMEGA,
        &[Layer::Porous {
            law: &law,
            thickness_m: LENGTH,
        }],
    )
    .unwrap();
    let beam = stack
        .prepare_angular_spectrum(
            &[IncidentPlaneWave {
                direction: [0.0, 0.0, 1.0],
                pressure_pa: C64::from_re(30.0),
            }],
            AngularSpectrumBudget { max_components: 1 },
            &mut || false,
        )
        .unwrap();
    let field = beam.prepare_layer_field(0, &mut || false).unwrap();
    let c = column();
    let routes = partition(2);
    let nodes = simpson(64);
    let boundary = |source_w| Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    };
    let events = vec![
        TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 8,
                a: id(0),
                b: id(1),
                area_m2: c.extent_m[0] * c.extent_m[1],
                conductance_w_k: 0.3,
            }),
        },
        TimedEvent {
            time_s: 0.0,
            event: Event::Boundary(id(0), boundary(2.0)),
        },
        TimedEvent {
            time_s: 0.4,
            event: Event::Boundary(id(0), boundary(5.0)),
        },
        TimedEvent {
            time_s: 1.0,
            event: Event::Boundary(id(1), boundary(7.0)),
        },
    ];
    let thermal = initial(c, &routes, 1.0, events, &CALORIC);
    let heat = CoherentLayerHeating::try_new(
        &field,
        &thermal,
        c,
        &routes,
        CoherentHeatQuadrature {
            axes: [&MID, &MID, &nodes],
        },
        budget(),
        || Ok(()),
    )
    .unwrap();
    let queries = fluid.samples.get();
    let first = heat.advance_to(&thermal, 0.4, config(), || Ok(())).unwrap();
    let second = heat
        .advance_to(&first.thermal, 1.0, config(), || Ok(()))
        .unwrap();
    let once = heat.advance_to(&thermal, 1.0, config(), || Ok(())).unwrap();
    let expected = total(&heat) + 2.0 * 0.4 + 5.0 * 0.6;
    close(once.heat_report.external_heat_j, expected, 1e-13);
    close(once.thermal.energy_j(), expected, 1e-13);
    close(second.thermal.energy_j(), once.thermal.energy_j(), 1e-13);
    close(
        first.acoustic.deposited_j + second.acoustic.deposited_j,
        once.acoustic.deposited_j,
        1e-16,
    );
    assert_eq!(once.thermal.cells()[&id(0)].boundary.source_w, 5.0);
    assert_eq!(once.thermal.cells()[&id(1)].boundary.source_w, 7.0);
    assert_eq!(once.thermal.links(), thermal.links());
    assert!(once.thermal.temperature_k(id(1)).unwrap() > 300.0);
    assert_eq!(fluid.samples.get(), queries);
}

#[test]
fn work_geometry_weights_binding_and_cancellation_refuse_transactionally() {
    let fluid = source();
    let law = material();
    let stack = TransmittingLayers::new(
        &fluid,
        &fluid,
        OMEGA,
        &[Layer::Porous {
            law: &law,
            thickness_m: LENGTH,
        }],
    )
    .unwrap();
    let beam = stack
        .prepare_angular_spectrum(
            &[IncidentPlaneWave {
                direction: [0.0, 0.0, 1.0],
                pressure_pa: C64::from_re(30.0),
            }],
            AngularSpectrumBudget { max_components: 1 },
            &mut || false,
        )
        .unwrap();
    let field = beam.prepare_layer_field(0, &mut || false).unwrap();
    let c = column();
    let routes = partition(2);
    let nodes = midpoint(8);
    let thermal = initial(c, &routes, 1.0, vec![], &CALORIC);
    let rule = CoherentHeatQuadrature {
        axes: [&nodes, &nodes, &nodes],
    };
    let mut polls = 0;
    assert!(
        CoherentLayerHeating::try_new(
            &field,
            &thermal,
            c,
            &routes,
            rule,
            CoherentHeatingBudget {
                max_component_evaluations: 100,
                ..budget()
            },
            || {
                polls += 1;
                Ok(())
            }
        )
        .is_err()
    );
    assert_eq!(
        polls, 1,
        "aggregate work admission precedes quadrature/field work"
    );
    let bad_weight = [AcousticHeatNode {
        coordinate: 0.5,
        weight: 0.9,
    }];
    assert!(
        CoherentLayerHeating::try_new(
            &field,
            &thermal,
            c,
            &routes,
            CoherentHeatQuadrature {
                axes: [&MID, &MID, &bad_weight]
            },
            budget(),
            || Ok(())
        )
        .is_err()
    );
    let wrong_volume = initial(c, &routes, 1.1, vec![], &CALORIC);
    assert!(
        CoherentLayerHeating::try_new(&field, &wrong_volume, c, &routes, rule, budget(), || Ok(()))
            .is_err()
    );
    let mut wrong_frame = c;
    wrong_frame.transverse_axes[1] = wrong_frame.transverse_axes[0];
    assert!(
        CoherentLayerHeating::try_new(
            &field,
            &thermal,
            wrong_frame,
            &routes,
            rule,
            budget(),
            || Ok(())
        )
        .is_err()
    );
    for stop in [2, 60, 1000] {
        let mut polls = 0;
        let result =
            CoherentLayerHeating::try_new(&field, &thermal, c, &routes, rule, budget(), || {
                polls += 1;
                if polls == stop {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            });
        assert!(matches!(
            result,
            Err(CoherentHeatingError::Thermal(DepositionError::Cancelled))
        ));
    }
    let heat =
        CoherentLayerHeating::try_new(&field, &thermal, c, &routes, rule, budget(), || Ok(()))
            .unwrap();
    let other = initial(c, &routes, 1.0, vec![], &OTHER_CALORIC);
    assert!(heat.advance_to(&other, 1.0, config(), || Ok(())).is_err());
    for stop in [2, 20, 100] {
        let mut polls = 0;
        let result = heat.advance_to(&thermal, 1.0, config(), || {
            polls += 1;
            if polls == stop {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            result,
            Err(CoherentHeatingError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.energy_j(), 0.0);
    }
    assert!(heat.advance_to(&thermal, 1.0, config(), || Ok(())).is_ok());
}

#[test]
fn lossless_layer_has_exactly_zero_deposited_heat() {
    let fluid = source();
    let stack = TransmittingLayers::new(
        &fluid,
        &fluid,
        OMEGA,
        &[Layer::FluidGap {
            thickness_m: LENGTH,
        }],
    )
    .unwrap();
    let beam = stack
        .prepare_angular_spectrum(
            &[IncidentPlaneWave {
                direction: [0.6, 0.0, 0.8],
                pressure_pa: C64::from_re(30.0),
            }],
            AngularSpectrumBudget { max_components: 1 },
            &mut || false,
        )
        .unwrap();
    let field = beam.prepare_layer_field(0, &mut || false).unwrap();
    let c = column();
    let routes = partition(1);
    let thermal = initial(c, &routes, 1.0, vec![], &CALORIC);
    let heat = CoherentLayerHeating::try_new(
        &field,
        &thermal,
        c,
        &routes,
        CoherentHeatQuadrature {
            axes: [&MID, &MID, &MID],
        },
        budget(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(total(&heat), 0.0);
    let step = heat.advance_to(&thermal, 1.0, config(), || Ok(())).unwrap();
    assert_eq!(step.heat_report.external_heat_j, 0.0);
    assert_eq!(step.acoustic.deposited_j, 0.0);
}
