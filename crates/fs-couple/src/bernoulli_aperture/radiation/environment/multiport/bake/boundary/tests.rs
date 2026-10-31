use super::*;
use core::f64::consts::TAU;
use fs_bem::helmholtz::{IterativeError, IterativeOptions};
use fs_bem::near_field::Options as ReceiverQuadrature;
use fs_bem::observation_transfers::{
    IterativeObservationOptions, LocalBoundary, LocalWall, LowerBandWindow, ReceiverSpec,
    TransferError,
};
use fs_bem::panel3d::SpherePanels;
use fs_material::fluid_acoustics::{FluidAcousticState, thermoviscous::ThermoviscousFluidState};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_phs::WallPin;
use fs_vfit::causality::CausalityOptions;
use std::sync::OnceLock;

pub(crate) fn gas() -> GasState {
    GasState::try_new_moist_air(293.15, 101325.0, 0.0).unwrap()
}
pub(crate) fn options() -> BaffledMimoBakeOptions {
    BaffledMimoBakeOptions {
        frequency_samples: 65,
        orders: [4, 8, 12],
        fit_iterations: 20,
        maximum_relative_error: 0.02,
        max_fit_memory_bytes: 64 * 1024 * 1024,
        causality: CausalityOptions {
            subtractions: 2,
            tail_bound_ratio: 100.0,
            relative_data_tolerance: 1e-8,
        },
    }
}
fn limits() -> MimoRadiationLimits {
    MimoRadiationLimits {
        max_states: 64,
        max_memory_bytes: 1 << 20,
        sampled_passivity_points: 257,
    }
}
pub(crate) fn sampling() -> BoundaryMimoSamplingLimits {
    BoundaryMimoSamplingLimits {
        maximum_retained_boundary_entries: 4096,
        maximum_port_matrix_entries: 4,
        maximum_source_solves: 2 * (3 * options().frequency_samples - 1),
        maximum_projection_panel_products: 1_000_000,
        maximum_reciprocity_defect: 1e-10,
    }
}
pub(crate) fn observer_options() -> IterativeObservationOptions {
    IterativeObservationOptions {
        highest_omega: TAU * 2000.0,
        minimum_clearance_m: 1e-4,
        receiver_quadrature: ReceiverQuadrature::default(),
        lower_window: None,
        solver: IterativeOptions {
            tol: 1e-10,
            restart: 80,
            max_iters: 600,
            fmm_order: 7,
            leaf_capacity: 256,
            ..IterativeOptions::default()
        },
    }
}

// Reflection-related finite sources retain identical physical work scales.
// The third, passive body changes both pressure reactions. These are synthetic
// boundary profiles, not an asserted real instrument aperture geometry.
pub(crate) fn geometry() -> (
    SpherePanels,
    Vec<SourceId>,
    Vec<Vec<f64>>,
    Vec<ReceiverSpec>,
) {
    let radiator = SpherePanels::icosphere(0.006, 0).unwrap();
    let obstacle = SpherePanels::icosphere(0.005, 0).unwrap();
    let mut triangles = Vec::new();
    for (body, offset) in [
        (&radiator, [-0.015, 0.0, 0.0]),
        (&radiator, [0.015, 0.0, 0.0]),
        (&obstacle, [0.0, 0.015, 0.0]),
    ] {
        triangles.extend(body.triangles().unwrap().iter().map(|triangle| {
            triangle.map(|point| std::array::from_fn(|axis| point[axis] + offset[axis]))
        }));
    }
    let surface = SpherePanels::from_triangles(triangles).unwrap();
    let sites = radiator.areas().len();
    let mut traces = vec![vec![0.0; surface.areas().len()]; 2];
    for (source, trace) in traces.iter_mut().enumerate() {
        let area: f64 = surface.areas()[source * sites..(source + 1) * sites]
            .iter()
            .sum();
        trace[source * sites..(source + 1) * sites].fill(1.0 / area);
    }
    (
        surface,
        vec![
            SourceId::ApertureVolumeFlow("left".into()),
            SourceId::ApertureVolumeFlow("right".into()),
        ],
        traces,
        vec![ReceiverSpec {
            id: "microphone".into(),
            position_m: [0.1, 0.07, 0.05],
        }],
    )
}
fn boundary(surface: &SpherePanels, wall: LocalWall) -> LocalBoundary {
    let mut patches = vec![0; surface.areas().len()];
    patches[2 * surface.areas().len() / 3..].fill(1);
    LocalBoundary::try_new_iterative(Medium::from(&gas()), vec![LocalWall::Rigid, wall], patches)
        .unwrap()
}
pub(crate) fn mechanical_boundary(surface: &SpherePanels) -> LocalBoundary {
    boundary(
        surface,
        LocalWall::Mechanical(WallPin {
            surface_density: 0.1,
            stiffness_per_area: 4e5,
            resistance: 200.0,
        }),
    )
}

/// Shared by the actual DG integration tests: cache only the expensive cold
/// physical bake, then clone its quiescent original runtime for each test.
pub(crate) fn baked_fixture() -> BoundaryMimoLoad {
    static SOURCE: OnceLock<BoundaryMimoLoad> = OnceLock::new();
    SOURCE
        .get_or_init(|| {
            let (surface, ids, traces, receivers) = geometry();
            let gate = CancelGate::new_clock_free();
            let owner = IterativeBoundaryObservation::for_surface(
                &surface,
                Medium::from(&gas()),
                &ids,
                &traces,
                &receivers,
                mechanical_boundary(&surface),
                observer_options(),
                &gate,
            )
            .unwrap();
            bake_boundary_mimo(
                &owner,
                gas(),
                4000,
                [20.0, 2000.0],
                options(),
                limits(),
                sampling(),
                &gate,
            )
            .unwrap()
        })
        .clone()
}

fn sampler<'a, 's>(owner: &'a IterativeBoundaryObservation<'s>) -> BoundarySampler<'a, 's> {
    BoundarySampler {
        observer: owner,
        sampling: sampling(),
        products_per_frequency: 4 * owner.surface_panel_count() as u64,
        statistics: RefCell::new(SamplingState {
            frequencies: 0,
            solves: 0,
            products: 0,
            residual: 0.0,
            reciprocity: 0.0,
            margin: f64::INFINITY,
            hash: DomainHasher::new("test-material-boundary-samples"),
        }),
    }
}
fn relative(actual: &[C64], expected: &[C64]) -> f64 {
    let error: f64 = actual
        .iter()
        .zip(expected)
        .map(|(&a, &b)| (a - b).norm_sq())
        .sum();
    let norm: f64 = expected.iter().map(|z| z.norm_sq()).sum();
    fs_math::det::sqrt(error / norm)
}

#[test]
fn physical_boundary_samples_keep_full_mutual_work_and_one_conjugation() {
    let (surface, ids, traces, receivers) = geometry();
    let gate = CancelGate::new_clock_free();
    let owner = IterativeBoundaryObservation::for_surface(
        &surface,
        Medium::from(&gas()),
        &ids,
        &traces,
        &receivers,
        mechanical_boundary(&surface),
        observer_options(),
        &gate,
    )
    .unwrap();
    let sampling_owner = sampler(&owner);
    let hz = [100.0, 600.0, 1800.0];
    let rows = sampling_owner.rows(&hz, &gate).unwrap();
    for (&frequency, row) in hz.iter().zip(&rows) {
        let response = owner.solve_frequency(TAU * frequency, 4096, &gate).unwrap();
        let physical = response
            .velocity_port_impedance(4, 1_000_000, &gate)
            .unwrap();
        assert_eq!(
            row,
            &physical
                .impedance()
                .iter()
                .map(|z| z.conj())
                .collect::<Vec<_>>()
        );
        assert!(row[1].abs() > 1e-3 * row[0].abs());
        assert!(row[0].im > 0.0);
        let rates = [C64::new(2e-7, -1e-7), C64::new(-3e-7, 4e-7)];
        let fit_rates = rates.map(|v| v.conj());
        let mut work = 0.0;
        for i in 0..2 {
            for j in 0..2 {
                work += 0.5 * (fit_rates[i].conj() * row[2 * i + j] * fit_rates[j]).re;
            }
        }
        let expected = physical.input_power_w(&rates).unwrap();
        assert!((work - expected).abs() <= 1e-12 * expected.abs().max(1e-300));
    }
    let report = sampling_owner.into_report();
    assert_eq!(report.frequency_samples, 3);
    assert_eq!(report.source_solves, 6);
    assert_eq!(
        report.projection_panel_products,
        12 * surface.areas().len() as u64
    );
}

#[test]
fn actual_material_bake_retains_fit_ladder_and_original_causal_runtime() {
    let source = baked_fixture();
    let record = source.record();
    assert_eq!(record.gas(), &gas());
    assert_eq!(record.time_step_s(), 1.0 / 4000.0);
    assert_eq!(
        record
            .report()
            .rungs
            .iter()
            .map(|r| r.order)
            .collect::<Vec<_>>(),
        options().orders
    );
    assert_eq!(
        record.sampling_report().frequency_samples,
        3 * options().frequency_samples - 1
    );
    assert_eq!(record.source_sample_band_hz(), [2000.0 / 129.0, 2000.0]);
    assert!(record.report().continuous.passive && record.report().sampled.passive);
    assert_eq!(
        source.load().reports().1.class,
        fs_vfit::passivity::CertificateClass::GridOnly
    );
    let (surface, ids, traces, receivers) = geometry();
    let gate = CancelGate::new_clock_free();
    let owner = IterativeBoundaryObservation::for_surface(
        &surface,
        Medium::from(&gas()),
        &ids,
        &traces,
        &receivers,
        mechanical_boundary(&surface),
        observer_options(),
        &gate,
    )
    .unwrap();
    assert_eq!(record.operator_hash(), owner.operator_hash());
    assert_eq!(
        record.aperture_support_areas_m2(),
        owner.aperture_support_areas_m2(&gate).unwrap()
    );
    // Fresh physical frequencies belong to neither identification nor its
    // interleaved validation grids. This is discrete-operator fit fidelity.
    let hz = [87.0, 337.0, 1173.0, 1969.0];
    let data = sampler(&owner).rows(&hz, &gate).unwrap();
    for (&frequency, row) in hz.iter().zip(data) {
        assert!(
            relative(&source.load().source_model().eval_iw(TAU * frequency), &row)
                <= options().maximum_relative_error
        );
    }
    let mut original = source
        .load()
        .source_model()
        .zoh(record.time_step_s())
        .unwrap();
    let mut actual = source.load().clone();
    for input in [[2e-7, -3e-7], [-1e-7, 4e-7], [0.0, 0.0]] {
        let expected = original.step_mean(&input);
        let mut output = [0.0; 2];
        actual.step_mean_into(&input, &mut output).unwrap();
        assert_eq!(output.as_slice(), expected.as_slice());
        assert_eq!(actual.pole_states(), original.pole_states());
    }
}

#[test]
fn physical_sampling_budgets_and_gas_precede_failed_source_solves() {
    let (surface, ids, traces, receivers) = geometry();
    let gate = CancelGate::new_clock_free();
    let mut config = observer_options();
    config.solver.max_iters = 0;
    let owner = IterativeBoundaryObservation::for_surface(
        &surface,
        Medium::from(&gas()),
        &ids,
        &traces,
        &receivers,
        mechanical_boundary(&surface),
        config,
        &gate,
    )
    .unwrap();
    let build = |budgets| {
        bake_boundary_mimo(
            &owner,
            gas(),
            4000,
            [20.0, 2000.0],
            options(),
            limits(),
            budgets,
            &gate,
        )
    };
    let mut short = sampling();
    short.maximum_retained_boundary_entries = 4 * surface.areas().len() - 1;
    assert!(matches!(
        build(short),
        Err(BoundaryMimoBakeError::Work {
            kind: "retained boundary entries",
            ..
        })
    ));
    short = sampling();
    short.maximum_source_solves -= 1;
    assert!(matches!(
        build(short),
        Err(BoundaryMimoBakeError::Work {
            kind: "complete source solves",
            ..
        })
    ));
    short = sampling();
    short.maximum_projection_panel_products =
        (3 * options().frequency_samples - 1) as u64 * 4 * surface.areas().len() as u64 - 1;
    assert!(matches!(
        build(short),
        Err(BoundaryMimoBakeError::Work {
            kind: "complete source-pair panel products",
            ..
        })
    ));
    let other_gas = GasState::try_new_moist_air(293.15, 81060.0, 0.0).unwrap();
    assert!(matches!(
        bake_boundary_mimo(
            &owner,
            other_gas,
            4000,
            [20.0, 2000.0],
            options(),
            limits(),
            sampling(),
            &gate
        ),
        Err(BoundaryMimoBakeError::Invalid(_))
    ));
    let Err(BoundaryMimoBakeError::Sample {
        frequency_hz,
        source: IterativeObservationError::Solver(error),
    }) = build(sampling())
    else {
        panic!("expected the original failed first source solve");
    };
    assert_eq!(frequency_hz, 20.0);
    assert!(matches!(*error, IterativeError::NotConverged { .. }));
    gate.request();
    assert!(matches!(
        build(sampling()),
        Err(BoundaryMimoBakeError::Cancelled)
    ));
}

#[test]
fn observation_windows_cannot_supply_missing_causal_load_support() {
    let (surface, ids, traces, receivers) = geometry();
    let gate = CancelGate::new_clock_free();
    let air = gas();
    let fluid = ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::try_new(air.density, air.sound_speed).unwrap(),
        FluidTransportState::try_new(
            air.density,
            air.specific_heat_cp,
            air.thermal_conductivity,
            air.dynamic_viscosity,
        )
        .unwrap(),
        air.gamma,
    )
    .unwrap();
    let law =
        PorousFluidLaw::try_new(PorousCard::basotect_foam02(), air.dynamic_viscosity).unwrap();
    let material = boundary(
        &surface,
        LocalWall::PorousRigid {
            law,
            fluid,
            thickness_m: 0.03,
        },
    );
    let mut config = observer_options();
    config.lower_window = Some(LowerBandWindow::try_new(150.0, 400.0).unwrap());
    let owner = IterativeBoundaryObservation::for_surface(
        &surface,
        Medium::from(&air),
        &ids,
        &traces,
        &receivers,
        material,
        config,
        &gate,
    )
    .unwrap();
    assert!(matches!(
        bake_boundary_mimo(
            &owner,
            air,
            4000,
            [400.0, 2000.0],
            options(),
            limits(),
            sampling(),
            &gate
        ),
        Err(BoundaryMimoBakeError::Source(
            IterativeObservationError::Transfer(TransferError::Invalid(_))
        ))
    ));

    let mechanical = boundary(
        &surface,
        LocalWall::Mechanical(WallPin {
            surface_density: 0.1,
            stiffness_per_area: 4e5,
            resistance: 0.0,
        }),
    );
    let mut mechanical_options = observer_options();
    mechanical_options.lower_window = Some(LowerBandWindow::try_new(400.0, 600.0).unwrap());
    let owner = IterativeBoundaryObservation::for_surface(
        &surface,
        Medium::from(&air),
        &ids,
        &traces,
        &receivers,
        mechanical,
        mechanical_options,
        &gate,
    )
    .unwrap();
    // Its native zero impedance lies strictly between grid samples. Entire
    // interval admission refuses before the fitter can miss that singularity.
    assert!(matches!(
        bake_boundary_mimo(
            &owner,
            air,
            4000,
            [20.0, 2000.0],
            options(),
            limits(),
            sampling(),
            &gate
        ),
        Err(BoundaryMimoBakeError::Source(_))
    ));
}
