//! Complete sourced thermal transport retains moving nonlinear diagnostics.
use super::*;
use fs_couple::circuit_enthalpy::{
    CircuitThermalCoefficients,
    conduction::{
        heat::ConductionHeatError,
        spatial::response::step_spatial_moving_conduction_response_candidate,
    },
    contacts::ContactThermalCoefficients,
};
use fs_em::{
    conduction::response::terminals::{
        ConductionTerminalControl, ConductionTerminalDrive, ConductionTerminalLimits,
    },
    spatial::conduction::{motion::SpatialConductionMotion, response::SpatialConductionResponse},
};
use fs_material::conductor::response::ScalarConductorResponse;

struct TaggedSample {
    heat: HeatNetworkSample,
    marker: u64,
}
impl AsRef<HeatNetworkSample> for TaggedSample {
    fn as_ref(&self) -> &HeatNetworkSample {
        &self.heat
    }
}
#[derive(Debug)]
struct ProviderError(HeatNetworkError);
impl From<HeatNetworkError> for ProviderError {
    fn from(error: HeatNetworkError) -> Self {
        Self(error)
    }
}
struct TaggedHeat<'p, 'n, 'a> {
    prepared: &'p PreparedHeatNetwork<'n, 'a>,
    refuse_hot: bool,
}
impl HeatAssembly for TaggedHeat<'_, '_, '_> {
    type Sample = TaggedSample;
    type Error = ProviderError;

    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self.prepared
    }

    fn resolve_heat(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        if self.refuse_hot && thermal.temperature_k(CellId(30)).unwrap() > 300.000_01 {
            return Err(CircuitEnthalpyError::Source(ProviderError(
                HeatNetworkError::Invalid("authored hot-source refusal"),
            )));
        }
        let coefficients = self
            .prepared
            .resolve_heat(cx, thermal)
            .map_err(|error| error.map_source(ProviderError))?;
        let bulk = coefficients.bulk;
        Ok(ContactThermalCoefficients {
            contacts: coefficients.contacts,
            bulk: CircuitThermalCoefficients {
                resistors: bulk.resistors,
                sources_w: bulk.sources_w,
                conductances_w_k: bulk.conductances_w_k,
                source: TaggedSample {
                    heat: bulk.source,
                    marker: 73,
                },
            },
        })
    }
}

// Independent C=[1,2], J=2, dt=.5 energy elimination with original pressure
// resistance and closed two-surface gray radiation; only theta_1 is bisected.
fn endpoint_with_radiation(old: [f64; 2], emissivity: f64, pressure: f64) -> [f64; 2] {
    let resistance = 2.0 - 1.8 * (pressure - 100000.0) / 100000.0;
    let residual = |first: f64| {
        let second = (old[0] + 2.0 * old[1] + 8.0 - 0.996 * first) / 2.0;
        let (t1, t2) = (300.0 + first, 300.0 + second);
        let (e1, e2) = (
            emissivity * (1.0 + 0.0001 * first),
            0.6 * (1.0 + 0.0001 * second),
        );
        let transfer = (first - second) / (1.0 + resistance)
            + 5.670_374_419e-8 * (t1.powi(4) - t2.powi(4)) / (1.0 / e1 + 1.0 / e2 - 1.0);
        (
            second,
            first - old[0] - 0.5 * (4.0 * (1.0 + 0.002 * first) - transfer),
        )
    };
    let (mut low, mut high) = (0.0, 32.0);
    assert!(residual(low).1 < 0.0 && residual(high).1 > 0.0);
    for _ in 0..100 {
        let mid = 0.5 * (low + high);
        if residual(mid).1 < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    let first = 0.5 * (low + high);
    [first, residual(first).0]
}

#[test]
fn complete_assembly_retains_typed_samples_sources_and_signed_moving_work() {
    for (emissivity, pressure, speed) in [
        (0.2, 100000.0, 2.0),
        (0.8, 100000.0, -10.0),
        (0.2, 200000.0, 0.0),
    ] {
        with_heat(
            emissivity,
            pressure,
            2000.0,
            |cx, gate, heat, mut thermal, field, regions| {
                let sources = field
                    .conductors()
                    .iter()
                    .enumerate()
                    .map(|(tet, conductor)| {
                        ScalarConductorResponse::new(
                            heat.electrical.bulk.cells[tet / 6].point.material,
                            conductor.clone(),
                            &[],
                        )
                        .unwrap()
                    })
                    .collect();
                let mut problem = SpatialConductionResponse::new(
                    field.positions_m().to_vec(),
                    field.mesh().tets.clone(),
                    sources,
                    limits(),
                    &mut || false,
                )
                .unwrap();
                let prepared = heat.prepare(cx, &thermal).unwrap();
                let provider = TaggedHeat {
                    prepared: &prepared,
                    refuse_hot: false,
                };
                let [positive, negative] = electrodes(&field, 2.0);
                let terminals = [negative.as_slice(), positive.as_slice()];
                let drive = ConductionTerminalDrive {
                    terminals: &terminals,
                    control: ConductionTerminalControl::Currents(&[-2.0, 2.0]),
                    limits: ConductionTerminalLimits {
                        maximum_terminals: 2,
                        maximum_workspace_values: 20_000,
                    },
                };
                let velocities = [[0.0, speed, 0.0]; 12];
                let magnetic = [[0.0, 0.0, 1.0]; 12];
                let motion = SpatialConductionMotion {
                    velocity_m_s: &velocities,
                    flux_density_t: &magnetic,
                };
                let mut old = [0.0; 2];
                for interval in 1..=2 {
                    let target = interval as f64 * 0.5;
                    let before = thermal.energy_j();
                    let next = step_spatial_moving_conduction_response_candidate(
                        &provider,
                        cx,
                        &thermal,
                        target,
                        &problem,
                        drive,
                        motion,
                        &[],
                        regions,
                        1e-12,
                        limits(),
                        field_controls(),
                        crate::response::controls(),
                        config(),
                        &[0.0; 2],
                    )
                    .unwrap();
                    let theta = endpoint_with_radiation(old, emissivity, pressure);
                    let report = next.report();
                    let step = &report.electrical;
                    let native = &step.field.field;
                    for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
                        near(
                            next.thermal().temperature_k(cell).unwrap(),
                            300.0 + theta[i],
                        );
                    }
                    near(native.current_density_a_m2()[0][0], -2.0);
                    near(
                        native.conductor_field_v_m()[0][0],
                        -2.0 * (1.0 + 0.002 * theta[0]),
                    );
                    near(native.lorentz_force_density_n_m3()[0][1], 2.0);
                    let heat_j = 2.0 * (4.0 + 0.002 * theta[0]);
                    near(step.heat_j, heat_j);
                    near(step.mechanical_work_j, 2.0 * speed);
                    near(step.supplied_work_j, heat_j + 2.0 * speed);
                    near(step.net_supplied_work_j, heat_j);
                    near(next.thermal().energy_j() - before, heat_j);
                    near(
                        0.5 * step.field.terminals.as_ref().unwrap().supplied_power_w,
                        step.supplied_work_j,
                    );
                    assert_eq!(report.source.marker, 73);
                    let sample = report.source.as_ref();
                    near(
                        sample.thermal_contacts[0].interface.value_m2_k_per_w(),
                        2.0 - 1.8 * (pressure - 100000.0) / 100000.0,
                    );
                    assert!(sample.radiation.is_some());
                    near(
                        report.total_energy_residual_j + 0.5 * sample.radiation_net_power_w(),
                        0.0,
                    );
                    if speed < 0.0 {
                        assert!(step.supplied_work_j < 0.0);
                        assert!(step.mechanical_work_j < 0.0);
                        assert!(step.heat_j > 0.0);
                    }
                    if interval == 1 {
                        let refused = TaggedHeat {
                            prepared: &prepared,
                            refuse_hot: true,
                        };
                        assert!(matches!(
                            step_spatial_moving_conduction_response_candidate(
                                &refused,
                                cx,
                                &thermal,
                                target,
                                &problem,
                                drive,
                                motion,
                                &[],
                                regions,
                                1e-12,
                                limits(),
                                field_controls(),
                                crate::response::controls(),
                                config(),
                                &[0.0; 2],
                            ),
                            Err(CircuitEnthalpyError::Source(ConductionHeatError::Heat(
                                ProviderError(HeatNetworkError::Invalid(
                                    "authored hot-source refusal"
                                ))
                            )))
                        ));
                        assert_eq!(thermal.energy_j(), before);
                    }
                    (thermal, problem, _) = next.into_parts();
                    old = theta;
                }
                let before = thermal.energy_j();
                gate.request();
                assert!(matches!(
                    step_spatial_moving_conduction_response_candidate(
                        &provider,
                        cx,
                        &thermal,
                        1.5,
                        &problem,
                        drive,
                        motion,
                        &[],
                        regions,
                        1e-12,
                        limits(),
                        field_controls(),
                        crate::response::controls(),
                        config(),
                        &[0.0; 2],
                    ),
                    Err(CircuitEnthalpyError::Cancelled)
                ));
                assert_eq!(thermal.energy_j(), before);
                near(thermal.time_s(), 1.0);
            },
        );
    }
}
