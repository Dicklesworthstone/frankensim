//! Independently reconstructed physical-volume diffusion on a tapered channel.
use super::*;
use crate::multicomponent::{Cell, Conserved};
use fs_qty::Pressure;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;

fn tapered(gas: MixtureGas, states: Vec<Conserved>) -> MixtureEuler1d {
    // Piecewise-linear area, integrated separately in each nonuniform cell.
    let mesh = Mesh1d::new(vec![0.0, 0.3, 1.0], vec![1.0, 2.0, 3.0], vec![0.45, 1.75]).unwrap();
    MixtureEuler1d::new_on_mesh(
        gas,
        mesh,
        states,
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        || false,
    )
    .unwrap()
}

fn advance<L: TransportLaw>(
    stabilized: bool,
    op: &mut DissipationOperator<L>,
    flow: &mut MixtureEuler1d,
    duration: f64,
    controls: DissipationControl,
    baths: [Option<Cell<'_>>; 2],
    cancelled: impl FnMut() -> bool,
) -> Result<DissipationReport, DissipationError> {
    if stabilized {
        op.advance_stabilized_with_reservoirs_for(
            flow,
            Time::new(duration),
            StabilizedDissipationControl {
                transport: controls,
                ..Default::default()
            },
            baths,
            cancelled,
        )
        .map(|report| report.transport)
    } else {
        op.advance_with_reservoirs_for(flow, Time::new(duration), controls, baths, cancelled)
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

fn closed_relaxation(stabilized: bool, thermal: bool, step: f64) -> f64 {
    let g = gas(if thermal { 1 } else { 2 }, 0.0, false);
    let initial = if thermal { [2.3, 1.8] } else { [0.8, 0.2] };
    let states = initial
        .map(|value| {
            if thermal {
                g.conserved(&[1.0], 0.0, value).unwrap()
            } else {
                g.conserved(&[value, 1.0 - value], 0.0, 2.0).unwrap()
            }
        })
        .to_vec();
    let mut flow = tapered(g, states);
    let material = if thermal {
        law(&flow, 0.0, 0.0, 0.05, 0.0)
    } else {
        law(&flow, 0.0, 0.0, 0.0, 0.03)
    };
    let mut op = DissipationOperator::new_diffusion_only(&flow, material, || false).unwrap();
    let report = advance(
        stabilized,
        &mut op,
        &mut flow,
        1.0,
        DissipationControl {
            maximum_step_s: step,
            ..Default::default()
        },
        [None; 2],
        || false,
    )
    .unwrap();
    // Independent two-volume ODE: V0*dq0/dt=G*(q1-q0),
    // V1*dq1/dt=-V0*dq0/dt; G=interface_area*diffusivity/center_distance.
    let diffusivity = if thermal { 0.05 / 2.5 } else { 0.03 };
    let volumes = [0.45, 1.75];
    let mean = (volumes[0] * initial[0] + volumes[1] * initial[1]) / 2.2;
    let rate: f64 = 2.0 * diffusivity / 0.5 * (1.0 / volumes[0] + 1.0 / volumes[1]);
    let expected = initial.map(|q| mean + (q - mean) * (-rate).exp());
    let mut error: f64 = 0.0;
    for (i, expected) in expected.into_iter().enumerate() {
        let cell = flow.cell(i).unwrap();
        let primitive = flow.gas().primitive(cell, || false).unwrap();
        let value = if thermal {
            primitive.temperature
        } else {
            cell.partial_densities[0]
        };
        error = error.max((value - expected).abs());
        near(cell.density(), 1.0, 3e-12);
        assert_eq!(cell.momentum, 0.0);
        if !thermal {
            near(primitive.temperature, 2.0, 3e-11);
        }
    }
    for k in 0..report.before.species_count + 2 {
        near(report.conservation_residual.get(k), 0.0, 2e-11);
        assert_eq!(report.boundary_change.get(k), 0.0);
    }
    assert_eq!(flow.time().value(), 0.0);
    error
}

#[test]
fn g1_variable_area_heat_and_species_follow_volume_weighted_relaxation() {
    for stabilized in [false, true] {
        for thermal in [false, true] {
            let coarse = closed_relaxation(stabilized, thermal, 0.1);
            let fine = closed_relaxation(stabilized, thermal, 0.05);
            assert!(coarse < 5e-5, "coarse={coarse:e}");
            assert!(fine < 0.28 * coarse, "coarse={coarse:e}, fine={fine:e}");
        }
    }
}

#[test]
fn g1_unequal_end_areas_close_each_bath_transfer_at_steady_heat_flow() {
    for stabilized in [false, true] {
        let g = gas(1, 0.0, false);
        let conductivity = 0.05;
        // Three thermal resistances in series, using actual half-cell lengths.
        let resistances = [
            0.15 / conductivity,
            0.5 / (2.0 * conductivity),
            0.35 / (3.0 * conductivity),
        ];
        let heat_rate = 2.0 / resistances.iter().sum::<f64>();
        let temperatures = [
            3.0 - heat_rate * resistances[0],
            1.0 + heat_rate * resistances[2],
        ];
        let left = g.conserved(&[1.0], 0.0, 3.0).unwrap();
        let right = g.conserved(&[1.0], 0.0, 1.0).unwrap();
        let states = temperatures
            .map(|t| g.conserved(&[1.0], 0.0, t).unwrap())
            .to_vec();
        let mut flow = tapered(g, states);
        let material = law(&flow, 0.0, 0.0, conductivity, 0.0);
        let mut op = DissipationOperator::new_diffusion_only(&flow, material, || false).unwrap();
        let report = advance(
            stabilized,
            &mut op,
            &mut flow,
            2.0,
            DissipationControl {
                maximum_step_s: 0.2,
                ..Default::default()
            },
            [Some(left.as_cell()), Some(right.as_cell())],
            || false,
        )
        .unwrap();
        near(report.left_boundary_change.energy, 2.0 * heat_rate, 3e-11);
        near(report.right_boundary_change.energy, -2.0 * heat_rate, 3e-11);
        near(report.boundary_change.energy, 0.0, 3e-11);
        near(report.absolute_energy_residual_j, 0.0, 3e-11);
        for (i, expected) in temperatures.into_iter().enumerate() {
            near(
                flow.gas()
                    .primitive(flow.cell(i).unwrap(), || false)
                    .unwrap()
                    .temperature,
                expected,
                3e-11,
            );
        }
    }
}

#[test]
fn g3_uniform_geometry_retains_original_transport_bit_for_bit() {
    for stabilized in [false, true] {
        let make = || {
            field(gas(2, 0.0, false), 8, Boundaries::Periodic, |x| {
                let a = 0.5 + 0.1 * (core::f64::consts::TAU * x).sin();
                (
                    vec![a, 1.0 - a],
                    0.1,
                    2.0 + 0.1 * (core::f64::consts::TAU * x).cos(),
                )
            })
        };
        let mut original = make();
        let mut scalar = make();
        let coupled = |flow: &MixtureEuler1d| {
            let heats = flow
                .gas()
                .source()
                .models()
                .iter()
                .zip([2.0, -1.0])
                .map(|(model, h)| {
                    (
                        model.species().clone(),
                        fs_thermochem::MassSpecificEnergyQuantityV1::new(h),
                    )
                })
                .collect();
            SoretDufourTransport::new(
                law(flow, 0.0, 0.0, 0.05, 0.03),
                "authored-cross-transport",
                1,
                TransportDomain {
                    minimum_temperature_k: 0.001,
                    maximum_temperature_k: 50.0,
                    minimum_pressure_pa: 1e-8,
                    maximum_pressure_pa: 1e8,
                },
                heats,
            )
            .unwrap()
        };
        let mut a = DissipationOperator::new(&original, coupled(&original), || false).unwrap();
        let material = coupled(&scalar);
        let mut b = DissipationOperator::new_diffusion_only(&scalar, material, || false).unwrap();
        let controls = DissipationControl::default();
        let expected = advance(
            stabilized,
            &mut a,
            &mut original,
            0.1,
            controls,
            [None; 2],
            || false,
        )
        .unwrap();
        let actual = advance(
            stabilized,
            &mut b,
            &mut scalar,
            0.1,
            controls,
            [None; 2],
            || false,
        )
        .unwrap();
        assert_eq!(scalar.states, original.states);
        assert_eq!(actual, expected);
    }
}

struct CompositionViscosity {
    inner: ConstantTransport,
    bulk: bool,
    admitted: usize,
    refused: usize,
}
impl TransportLaw for CompositionViscosity {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.inner.source()
    }
    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        mole_fractions: &[f64],
        binary: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        let mut value =
            self.inner
                .evaluate(temperature, pressure, mole_fractions, binary, cancelled)?;
        if mole_fractions[0] > 0.5 {
            self.refused += 1;
            if self.bulk {
                value.bulk_viscosity = DynViscosity::new(0.01);
            } else {
                value.dynamic_viscosity = DynViscosity::new(0.01);
            }
        } else {
            self.admitted += 1;
        }
        Ok(value)
    }
}

#[test]
fn g4_viscosity_appearing_in_a_later_stage_refuses_without_changing_inventory() {
    for stabilized in [false, true] {
        for bulk in [false, true] {
            let g = gas(2, 0.0, false);
            let bath = g.conserved(&[0.8, 0.2], 0.0, 2.0).unwrap();
            let states = vec![g.conserved(&[0.2, 0.8], 0.0, 2.0).unwrap(); 2];
            let mut flow = tapered(g, states);
            let before = flow.states.clone();
            let material = CompositionViscosity {
                inner: law(&flow, 0.0, 0.0, 0.0, 0.03),
                bulk,
                admitted: 0,
                refused: 0,
            };
            let mut op =
                DissipationOperator::new_diffusion_only(&flow, material, || false).unwrap();
            let error = advance(
                stabilized,
                &mut op,
                &mut flow,
                0.1,
                DissipationControl {
                    maximum_halvings: 0,
                    ..Default::default()
                },
                [Some(bath.as_cell()), None],
                || false,
            )
            .unwrap_err();
            assert!(
                matches!(error, DissipationError::StepFailure { .. }),
                "{error:?}"
            );
            assert!(op.law().admitted > 0 && op.law().refused > 0);
            assert_eq!(flow.states, before);
            assert_eq!(flow.time().value(), 0.0);
        }
    }
}

#[test]
fn g4_variable_area_cancellation_is_atomic_and_workspace_is_reusable() {
    for stabilized in [false, true] {
        let make = || {
            let g = gas(1, 0.0, false);
            let states = [2.3, 1.8]
                .map(|t| g.conserved(&[1.0], 0.0, t).unwrap())
                .to_vec();
            tapered(g, states)
        };
        let build = |flow: &MixtureEuler1d| {
            DissipationOperator::new_diffusion_only(flow, law(flow, 0.0, 0.0, 0.05, 0.0), || false)
                .unwrap()
        };
        let controls = DissipationControl {
            maximum_step_s: 0.05,
            ..Default::default()
        };
        let mut expected = make();
        let mut reference = build(&expected);
        let mut polls = 0;
        advance(
            stabilized,
            &mut reference,
            &mut expected,
            0.2,
            controls,
            [None; 2],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
        let mut flow = make();
        let before = flow.states.clone();
        let mut op = build(&flow);
        let mut seen = 0;
        let error = advance(
            stabilized,
            &mut op,
            &mut flow,
            0.2,
            controls,
            [None; 2],
            || {
                seen += 1;
                seen == polls
            },
        )
        .unwrap_err();
        assert_eq!(error, DissipationError::Cancelled);
        assert_eq!(flow.states, before);
        advance(
            stabilized,
            &mut op,
            &mut flow,
            0.2,
            controls,
            [None; 2],
            || false,
        )
        .unwrap();
        assert_eq!(flow.states, expected.states);
        assert!(matches!(
            op.admit_motion(&flow, &mut || false),
            Err(DissipationError::Invalid(_))
        ));
    }
}
