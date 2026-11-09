use super::*;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
use fs_material::fluid_transport::FluidTransportState;
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_math::c64::C64;
use fs_phs::linear_exact::DenseExponentialBudget;

const RADIUS: f64 = 0.02;
const LENGTH: f64 = 0.03;
const LEFT: PortId = PortId(17);
const RIGHT: PortId = PortId(19);
fn geometry() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: RADIUS,
            length: LENGTH,
        }],
    }
}
fn fluid(gamma: f64) -> ThermoviscousFluidState {
    ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::from_adiabatic_bulk_modulus(1.204, 1.4 * 101_325.0).unwrap(),
        FluidTransportState::try_new(1.204, 1005.0, 1.81e-5 * 1005.0 / 0.71, 1.81e-5).unwrap(),
        gamma,
    )
    .unwrap()
}
fn model(state: &ThermoviscousFluidState) -> PorousDiffusive {
    let law = PorousFluidLaw::try_new(PorousCard::melamine_uf_nguyen2024(), 1.81e-5).unwrap();
    PorousDiffusive::new(
        &law,
        state,
        8,
        [
            100.0 * core::f64::consts::TAU,
            4000.0 * core::f64::consts::TAU,
        ],
        5e-5,
    )
    .unwrap()
}
fn duct(
    gamma: f64,
    inlet: DgBoundary,
    outlet: DgBoundary,
) -> DgDuct<DgPorousMaterial<ThermoviscousFluidState>> {
    let state = fluid(gamma);
    DgDuct::from_porous_material(
        &geometry(),
        &model(&state),
        &state,
        DgMesh {
            elements: 1,
            order: 1,
        },
        inlet,
        outlet,
    )
    .unwrap()
}
fn budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 64 << 20,
        max_multiply_terms: 40_000_000_000,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-30),
        "{a:e} != {b:e}"
    );
}

#[test]
fn source_bound_porous_storage_and_darcy_heat_use_full_area_once() {
    let dg = duct(1.4, DgBoundary::IdealOpen, DgBoundary::IdealOpen);
    assert_eq!(dg.porous_material().source(), &fluid(1.4));
    let model = dg.porous_diffusive_model().unwrap();
    assert_eq!(dg.viscothermal_band_rad_s(), Some(model.band_rad_s()));
    assert!(
        dg.viscothermal_model().is_none(),
        "no added rigid circular-wall law"
    );
    let flow = 1e-6;
    let mut state = dg.project_initial_state(|_| 0.0, |_| flow);
    let m = dg.field_dofs;
    for (arm, term) in model.series_terms().iter().enumerate() {
        for i in 0..m {
            state.tilde[2 * m + arm * m + i] =
                det::sqrt(term.gain_per_s / term.rate_per_s) * state.tilde[m + i];
        }
    }
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    let sigma = model.fluid().transport().dynamic_viscosity_pa_s()
        / model.material().viscous_permeability_m2();
    let power = dg.porous_dissipated_power_w(&state).unwrap();
    close(power.viscous_w, sigma * LENGTH * flow * flow / area, 2e-12);
    assert_eq!(power.thermal_w, 0.0);
    close(
        dg.wall_dissipated_power_w(&state).unwrap(),
        power.viscous_w,
        1e-14,
    );
    assert!(
        dg.continuous_operator()
            .admit_frequency(50.0 * core::f64::consts::TAU)
            .is_err()
    );
    assert!(
        DgDuct::from_porous_material(
            &geometry(),
            model,
            &fluid(1.0),
            DgMesh {
                elements: 1,
                order: 1
            },
            DgBoundary::Closed,
            DgBoundary::Closed
        )
        .is_err()
    );
}

#[test]
fn both_porous_heat_channels_close_native_energy_and_split_intervals() {
    let dg = duct(1.4, DgBoundary::Closed, DgBoundary::Closed);
    let initial = dg.project_initial_state(
        |x| 1.0 + 0.2 * det::cos(core::f64::consts::PI * x / LENGTH),
        |_| 1e-6,
    );
    let step = dg.propagator(2e-6, 1e-11).unwrap();
    let heat = dg
        .prepare_porous_heat(&step, budget(), &mut || false)
        .unwrap();
    let observed = heat
        .integrated_heat_j(&initial, &[], &mut || false)
        .unwrap();
    assert!(observed.viscous_j > 0.0 && observed.thermal_j > 0.0);
    let mut end = initial.clone();
    step.step(&mut end, &[]).unwrap();
    close(
        observed.total_j(),
        initial.stored_energy_j() - end.stored_energy_j(),
        2e-8,
    );

    let half = dg.propagator(step.dt() / 2.0, 1e-11).unwrap();
    let half_heat = dg
        .prepare_porous_heat(&half, budget(), &mut || false)
        .unwrap();
    let first = half_heat
        .integrated_heat_j(&initial, &[], &mut || false)
        .unwrap();
    let mut middle = initial.clone();
    half.step(&mut middle, &[]).unwrap();
    let second = half_heat
        .integrated_heat_j(&middle, &[], &mut || false)
        .unwrap();
    close(first.viscous_j + second.viscous_j, observed.viscous_j, 3e-9);
    close(first.thermal_j + second.thermal_j, observed.thermal_j, 3e-9);
    assert_eq!(
        initial.step_count(),
        0,
        "observation never advances mechanics"
    );
}

#[test]
fn forced_native_power_matches_separate_heat_without_warm_allocation() {
    let dg = duct(1.4, DgBoundary::Port(LEFT), DgBoundary::Closed);
    let initial = dg.project_initial_state(|_| 0.0, |_| 0.0);
    let step = dg.propagator(2e-6, 1e-11).unwrap();
    let heat = dg
        .prepare_porous_heat(&step, budget(), &mut || false)
        .unwrap();
    let input = [PortInput {
        id: LEFT,
        effort: 1e-6,
    }];
    let mut scratch = vec![0.0; heat.scratch_len()];
    let observed = heat
        .integrated_heat_into_j(&initial, &input, &mut scratch, &mut || false)
        .unwrap();
    assert_eq!(
        observed,
        heat.integrated_heat_j(&initial, &input, &mut || false)
            .unwrap()
    );
    let mut state = initial.clone();
    let record = step.step(&mut state, &input).unwrap();
    assert!(observed.viscous_j > 0.0 && observed.thermal_j > 0.0);
    close(
        observed.total_j(),
        record.ledger.port_work_step_j - state.stored_energy_j(),
        2e-8,
    );
    assert_eq!(initial.step_count(), 0);
}

#[test]
fn g1_element_heat_partitions_original_sources_and_preserves_physical_volume() {
    let fluid = fluid(1.4);
    let dg = DgDuct::from_porous_material(
        &geometry(),
        &model(&fluid),
        &fluid,
        DgMesh {
            elements: 2,
            order: 1,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap();
    let step = dg.propagator(2e-6, 1e-11).unwrap();
    let whole = dg
        .prepare_porous_heat(&step, budget(), &mut || false)
        .unwrap();
    let cells: Vec<_> = (0..2)
        .map(|e| {
            dg.prepare_porous_element_heat(&step, e, budget(), &mut || false)
                .unwrap()
        })
        .collect();
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    for (index, cell) in cells.iter().enumerate() {
        assert_eq!(cell.element(), index);
        assert_eq!(
            cell.interval_m(),
            [
                index as f64 * LENGTH / 2.0,
                (index + 1) as f64 * LENGTH / 2.0
            ]
        );
        close(cell.area_m2(), area, 1e-15);
        assert_eq!(cell.heat().operator_identity(), whole.operator_identity());
    }
    // Uniform compression relaxes identically in two equal physical volumes.
    // This checks the native element mass map, not arbitrary fractions of total
    // heat: the same caches below also see strongly nonuniform illumination.
    let uniform = dg.project_initial_state(|_| 1.0, |_| 0.0);
    let left = cells[0]
        .heat()
        .integrated_heat_j(&uniform, &[], &mut || false)
        .unwrap();
    let right = cells[1]
        .heat()
        .integrated_heat_j(&uniform, &[], &mut || false)
        .unwrap();
    assert!(left.thermal_j > 0.0);
    close(left.thermal_j, right.thermal_j, 2e-8);
    let total = whole
        .integrated_heat_j(&uniform, &[], &mut || false)
        .unwrap();
    close(left.thermal_j + right.thermal_j, total.thermal_j, 2e-8);

    let local = dg.project_initial_state(|x| if x < LENGTH / 2.0 { 1.0 } else { 0.0 }, |_| 0.0);
    let left = cells[0]
        .heat()
        .integrated_heat_j(&local, &[], &mut || false)
        .unwrap();
    let right = cells[1]
        .heat()
        .integrated_heat_j(&local, &[], &mut || false)
        .unwrap();
    let total = whole.integrated_heat_j(&local, &[], &mut || false).unwrap();
    assert!(
        left.thermal_j > 100.0 * right.thermal_j,
        "{left:?} {right:?}"
    );
    close(left.viscous_j + right.viscous_j, total.viscous_j, 2e-8);
    close(left.thermal_j + right.thermal_j, total.thermal_j, 2e-8);
    let mut end = local.clone();
    step.step(&mut end, &[]).unwrap();
    close(
        total.total_j(),
        local.stored_energy_j() - end.stored_energy_j(),
        3e-8,
    );
    assert_eq!(local.step_count(), 0);
}

#[test]
fn g3_element_heat_keeps_interelement_transport_and_held_port_work() {
    let fluid = fluid(1.0);
    let dg = DgDuct::from_porous_material(
        &geometry(),
        &model(&fluid),
        &fluid,
        DgMesh {
            elements: 2,
            order: 1,
        },
        DgBoundary::Port(LEFT),
        DgBoundary::Closed,
    )
    .unwrap();
    let initial = dg.zero_state();
    let step = dg.propagator(2e-5, 1e-11).unwrap();
    let half = dg.propagator(step.dt() / 2.0, 1e-11).unwrap();
    let input = [PortInput {
        id: LEFT,
        effort: 1e-6,
    }];
    let mut middle = initial.clone();
    half.step(&mut middle, &input).unwrap();
    let mut sum = 0.0;
    for element in 0..2 {
        let full = dg
            .prepare_porous_element_heat(&step, element, budget(), &mut || false)
            .unwrap();
        let partial = dg
            .prepare_porous_element_heat(&half, element, budget(), &mut || false)
            .unwrap();
        let total = full
            .heat()
            .integrated_heat_j(&initial, &input, &mut || false)
            .unwrap();
        let first = partial
            .heat()
            .integrated_heat_j(&initial, &input, &mut || false)
            .unwrap();
        let second = partial
            .heat()
            .integrated_heat_j(&middle, &input, &mut || false)
            .unwrap();
        assert!(total.viscous_j > 0.0, "drive must reach element {element}");
        assert_eq!(total.thermal_j, 0.0);
        close(first.viscous_j + second.viscous_j, total.viscous_j, 2e-8);
        sum += total.viscous_j;
    }
    let mut end = initial.clone();
    let record = step.step(&mut end, &input).unwrap();
    close(
        sum,
        record.ledger.port_work_step_j - end.stored_energy_j(),
        3e-8,
    );
}

#[test]
fn g4_element_heat_refuses_foreign_operator_bounds_budget_and_cancellation() {
    let dg = duct(1.0, DgBoundary::Closed, DgBoundary::Closed);
    let state = dg.project_initial_state(|_| 1.0, |_| 1e-6);
    let original = state.clone();
    let step = dg.propagator(2e-6, 1e-11).unwrap();
    assert!(
        dg.prepare_porous_element_heat(&step, 1, budget(), &mut || false)
            .is_err()
    );
    assert!(
        dg.prepare_porous_element_heat(&step, 0, budget(), &mut || true)
            .is_err()
    );
    assert!(
        dg.prepare_porous_element_heat(
            &step,
            0,
            DenseExponentialBudget {
                max_dense_bytes: 0,
                ..budget()
            },
            &mut || false,
        )
        .is_err()
    );
    let other = duct(1.4, DgBoundary::Closed, DgBoundary::Closed);
    let other_step = other.propagator(2e-6, 1e-11).unwrap();
    assert!(
        dg.prepare_porous_element_heat(&other_step, 0, budget(), &mut || false)
            .is_err()
    );
    let heat = dg
        .prepare_porous_element_heat(&step, 0, budget(), &mut || false)
        .unwrap();
    assert!(
        heat.heat()
            .integrated_heat_j(&other.zero_state(), &[], &mut || false)
            .is_err()
    );
    assert_eq!(state, original);
}

#[test]
fn zero_thermal_coupling_and_refusals_preserve_original_state() {
    let dg = duct(1.0, DgBoundary::Closed, DgBoundary::Closed);
    let state = dg.project_initial_state(|_| 1.0, |_| 1e-6);
    let original = state.clone();
    let step = dg.propagator(2e-6, 1e-11).unwrap();
    let heat = dg
        .prepare_porous_heat(&step, budget(), &mut || false)
        .unwrap();
    let result = heat.integrated_heat_j(&state, &[], &mut || false).unwrap();
    assert!(result.viscous_j > 0.0);
    assert_eq!(result.thermal_j, 0.0);
    assert!(heat.integrated_heat_j(&state, &[], &mut || true).is_err());
    assert!(
        heat.integrated_heat_j(
            &state,
            &[PortInput {
                id: LEFT,
                effort: 1.0
            }],
            &mut || false
        )
        .is_err()
    );
    let other =
        duct(1.4, DgBoundary::Closed, DgBoundary::Closed).project_initial_state(|_| 1.0, |_| 0.0);
    assert!(heat.integrated_heat_j(&other, &[], &mut || false).is_err());
    assert!(
        dg.prepare_porous_heat(
            &step,
            DenseExponentialBudget {
                max_dimension: 1,
                ..budget()
            },
            &mut || false
        )
        .is_err()
    );
    assert_eq!(state, original);
}

#[test]
fn complete_porous_dg_scattering_converges_to_original_coherent_layer() {
    let fluid = fluid(1.4);
    let model = model(&fluid);
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    let z = fluid.acoustic().characteristic_impedance_pa_s_m() / area;
    let incoming = [C64::new(1.0, 0.3), C64::new(-0.4, 0.2)];
    for hz in [300.0, 1500.0, 3500.0] {
        let omega = core::f64::consts::TAU * hz;
        let reference =
            crate::porous::PorousLayer::try_new(model.material(), &fluid, omega, LENGTH, area, z)
                .unwrap()
                .power_for_incident_pressure(incoming)
                .unwrap();
        let mut errors = Vec::new();
        for elements in [1, 2, 4] {
            let dg = DgDuct::from_porous_material(
                &geometry(),
                &model,
                &fluid,
                DgMesh { elements, order: 3 },
                DgBoundary::Port(LEFT),
                DgBoundary::Port(RIGHT),
            )
            .unwrap();
            let factor = dg
                .scattering_factor(
                    omega,
                    &[(LEFT, z), (RIGHT, z)],
                    100_000_000,
                    1e-10,
                    &mut || false,
                )
                .unwrap();
            let result = factor
                .response(&[(LEFT, incoming[0]), (RIGHT, incoming[1])], &mut || false)
                .unwrap();
            let error = result
                .waves
                .iter()
                .zip(reference.outgoing_pressure_pa)
                .map(|(wave, exact)| (wave.reflected_pressure_pa - exact).abs())
                .fold(0.0_f64, f64::max);
            errors.push(error);
            assert!(result.physical.wall_dissipation_w > 0.0);
            assert!(result.physical.power_residual_w.abs() < 1e-9 * reference.incident_power_w);
            if elements == 4 {
                assert!(error < 2e-4, "{hz} Hz: {errors:?}");
                close(
                    result.physical.wall_dissipation_w,
                    reference.absorbed_power_w,
                    5e-4,
                );
            }
        }
        assert!(errors[2] < errors[0].max(5e-5), "{hz} Hz: {errors:?}");
    }
}
