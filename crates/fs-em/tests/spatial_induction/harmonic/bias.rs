//! G1/G3: a resolved DC magnetic field restores the missing carrier force.
use super::*;
use fs_em::spatial::vector_potential::{
    SpatialVectorField,
    induction::{
        harmonic::SpatialHarmonicInductionResponse,
        mechanical::{SpatialLorentzLoadError, SpatialLorentzLoadLimits},
    },
};
use fs_ga::Se3;

const LOADS: SpatialLorentzLoadLimits = SpatialLorentzLoadLimits {
    max_vertices: 8,
    max_cells: 6,
};
const ALL: [usize; 6] = [0, 1, 2, 3, 4, 5];

fn virtual_velocity(p: [f64; 3]) -> [f64; 3] {
    [
        0.3 + 2.0 * p[0] - 0.2 * p[2],
        -0.4 + 0.7 * p[1],
        0.8 - 0.5 * p[0] + 1.3 * p[2],
    ]
}

// Reconstruct real J and total real B BEFORE multiplication. Independent
// degree-two tetrahedral quadrature then integrates force, torque and virtual
// work; neither Fourier-force formulas nor production P1 mass weights appear.
fn instantaneous_reference(
    responses: &[&SpatialHarmonicInductionResponse<'_, '_>],
    bias: &SpatialVectorField<'_>,
    selected: &[usize],
    time: f64,
) -> (Vec<[f64; 3]>, [f64; 6], f64) {
    let problem = responses[0].problem();
    let mut nodal = vec![[0.0; 3]; problem.mesh().vertex_count];
    let mut wrench = [0.0; 6];
    let mut power = 0.0;
    let high = (5.0 + 3.0 * 5.0_f64.sqrt()) / 20.0;
    let low = (5.0 - 5.0_f64.sqrt()) / 20.0;
    for &cell in selected {
        let original = responses[0]
            .cells()
            .iter()
            .find(|r| r.cell == cell)
            .unwrap();
        let tet = problem.mesh().tets[cell].map(|v| v as usize);
        let mut b = bias.flux_density_t()[cell];
        let mut current = [[0.0; 3]; 4];
        for response in responses {
            let (sine, cosine) = (response.omega_rad_s() * time).sin_cos();
            let real = |z: C64| z.re * cosine + z.im * sine;
            let record = response.cells().iter().find(|r| r.cell == cell).unwrap();
            for axis in 0..3 {
                b[axis] += real(response.flux_density_t()[cell][axis]);
                for vertex in 0..4 {
                    current[vertex][axis] += real(record.current_density_a_m2[vertex][axis]);
                }
            }
        }
        for point in 0..4 {
            let shape: [f64; 4] = std::array::from_fn(|i| if i == point { high } else { low });
            let j = std::array::from_fn(|axis| (0..4).map(|i| shape[i] * current[i][axis]).sum());
            let position = std::array::from_fn(|axis| {
                (0..4)
                    .map(|i| shape[i] * problem.positions_m()[tet[i]][axis])
                    .sum()
            });
            let force = cross(j, b);
            let torque = cross(position, force);
            let velocity = virtual_velocity(position);
            let weight = original.volume_m3 / 4.0;
            for axis in 0..3 {
                wrench[axis] += weight * torque[axis];
                wrench[axis + 3] += weight * force[axis];
                power += weight * force[axis] * velocity[axis];
                for vertex in 0..4 {
                    nodal[tet[vertex]][axis] += weight * shape[vertex] * force[axis];
                }
            }
        }
    }
    (nodal, wrench, power)
}

#[test]
fn static_bias_restores_carrier_force_torque_and_nonrigid_virtual_work() {
    let (problem, sources) = cube(1, None);
    let initial = problem
        .solve(
            &[[0.0; 3]; 6],
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let response = state
        .solve_harmonic(
            &[[C64::new(3.2, -2.4); 3]; 6],
            &[[0.0; 3]; 8],
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let virtual_velocity: Vec<_> = problem
        .positions_m()
        .iter()
        .copied()
        .map(virtual_velocity)
        .collect();
    // Boundary-supported original static fields in two distinct directions.
    // These conducting cells have no impressed DC current.
    for bias_field in [[0.23, -0.17, 0.61], [-0.3, 0.4, 0.2]] {
        let bias = problem
            .solve(
                &[[0.0; 3]; 6],
                &trace(&problem, bias_field, false),
                static_controls(),
                &mut || false,
            )
            .unwrap();
        for selected in [&ALL[..], &[5, 0, 3][..]] {
            let ac = response
                .lorentz_loads(selected, LOADS, &mut || false)
                .unwrap();
            let dc = response
                .static_bias_lorentz_loads(&bias, selected, LOADS, &mut || false)
                .unwrap();
            near(dc.omega_rad_s(), 3.0, 1e-15);
            let coefficients = [
                ac.mean(),
                ac.second_harmonic_real(),
                ac.second_harmonic_imaginary(),
                dc.real(),
                dc.imaginary(),
            ];
            assert!(
                coefficients[3..]
                    .iter()
                    .flat_map(|c| c.nodal_forces_n())
                    .flatten()
                    .any(|v| v.abs() > 1e-3)
            );
            let wrenches =
                coefficients.map(|c| c.wrench(&Se3::identity(), [0.0; 3]).unwrap().to_array());
            let powers = coefficients.map(|c| c.power_w(&virtual_velocity, &mut || false).unwrap());
            for phase in [0.0_f64, 0.31, 1.2, 2.8, 4.7] {
                let (sine, cosine) = phase.sin_cos();
                let (sine2, cosine2) = (2.0_f64 * phase).sin_cos();
                let weights = [1.0, cosine2, sine2, cosine, sine];
                let (expected_nodes, expected_wrench, expected_power) =
                    instantaneous_reference(&[&response], &bias, selected, phase / 3.0);
                for (node, expected) in expected_nodes.iter().enumerate() {
                    for axis in 0..3 {
                        let actual = (0..5)
                            .map(|i| weights[i] * coefficients[i].nodal_forces_n()[node][axis])
                            .sum();
                        near(actual, expected[axis], 2e-10);
                    }
                }
                for axis in 0..6 {
                    near(
                        (0..5).map(|i| weights[i] * wrenches[i][axis]).sum(),
                        expected_wrench[axis],
                        2e-10,
                    );
                }
                near(
                    (0..5).map(|i| weights[i] * powers[i]).sum(),
                    expected_power,
                    2e-10,
                );
            }
        }
    }
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.potential_wb());
}

#[test]
fn static_bias_preserves_gauge_sign_source_support_and_atomic_refusal() {
    let (problem, sources) = cube(1, None);
    let initial = problem
        .solve(
            &[[0.0; 3]; 6],
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let current = [[C64::new(4.0, 1.0); 3]; 6];
    let response = state
        .solve_harmonic(
            &current,
            &[[0.0; 3]; 8],
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let solve_bias = |b, gauge| {
        problem
            .solve(
                &[[0.0; 3]; 6],
                &trace(&problem, b, gauge),
                static_controls(),
                &mut || false,
            )
            .unwrap()
    };
    let bias = solve_bias([0.2, -0.1, 0.8], false);
    let reversed = solve_bias([-0.2, 0.1, -0.8], false);
    let gauged = solve_bias([0.2, -0.1, 0.8], true);
    let mut calls = 0;
    let reference = response
        .static_bias_lorentz_loads(&bias, &ALL, LOADS, &mut || {
            calls += 1;
            false
        })
        .unwrap();
    for (field, scale) in [(&reversed, -1.0), (&gauged, 1.0), (&initial, 0.0)] {
        let actual = response
            .static_bias_lorentz_loads(field, &[5, 4, 3, 2, 1, 0], LOADS, &mut || false)
            .unwrap();
        for (a, b) in [
            (actual.real(), reference.real()),
            (actual.imaginary(), reference.imaginary()),
        ] {
            for (a, b) in a.nodal_forces_n().iter().zip(b.nodal_forces_n()) {
                for axis in 0..3 {
                    near(a[axis], scale * b[axis], 2e-10);
                }
            }
        }
    }
    let (other, _) = cube(1, None);
    let foreign = other
        .solve(
            &[[0.0; 3]; 6],
            &trace(&other, [0.2, -0.1, 0.8], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let mut partial_trace = trace(&problem, [0.0; 3], false);
    partial_trace.potential.truncate(1);
    let different_support = problem
        .solve(
            &[[0.0; 3]; 6],
            &partial_trace,
            static_controls(),
            &mut || false,
        )
        .unwrap();
    assert_ne!(
        different_support.prescribed_edges(),
        initial.prescribed_edges()
    );
    let dc_driven = problem
        .solve(
            &[[1.7; 3]; 6],
            &trace(&problem, [0.2, -0.1, 0.8], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    for field in [&foreign, &different_support, &dc_driven] {
        assert!(matches!(
            response.static_bias_lorentz_loads(field, &ALL, LOADS, &mut || false),
            Err(SpatialLorentzLoadError::Invalid(_))
        ));
    }
    let moving = state
        .solve_harmonic(
            &current,
            &[[0.1, 0.0, 0.0]; 8],
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    assert!(matches!(
        moving.static_bias_lorentz_loads(&bias, &ALL, LOADS, &mut || false),
        Err(SpatialLorentzLoadError::Invalid(_))
    ));
    for selected in [&[][..], &[0, 0][..], &[6][..]] {
        assert!(matches!(
            response.static_bias_lorentz_loads(&bias, selected, LOADS, &mut || false),
            Err(SpatialLorentzLoadError::Invalid(_))
        ));
    }
    for limit in [
        SpatialLorentzLoadLimits {
            max_vertices: 7,
            ..LOADS
        },
        SpatialLorentzLoadLimits {
            max_cells: 5,
            ..LOADS
        },
    ] {
        assert!(matches!(
            response.static_bias_lorentz_loads(&bias, &ALL, limit, &mut || false),
            Err(SpatialLorentzLoadError::Invalid(_))
        ));
    }
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        assert!(matches!(
            response.static_bias_lorentz_loads(&bias, &ALL, LOADS, &mut || {
                count += 1;
                count == stop
            }),
            Err(SpatialLorentzLoadError::Cancelled)
        ));
    }
    let replay = response
        .static_bias_lorentz_loads(&bias, &ALL, LOADS, &mut || false)
        .unwrap();
    assert_eq!(
        reference.real().nodal_forces_n(),
        replay.real().nodal_forces_n()
    );
    assert_eq!(
        reference.imaginary().nodal_forces_n(),
        replay.imaginary().nodal_forces_n()
    );
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.potential_wb());
}

#[cfg(feature = "harmonic-exposure")]
mod exposure {
    use super::*;
    use fs_em::spatial::vector_potential::induction::harmonic::coherent_exposure::{
        SpatialHarmonicInductionExposureError as Error,
        SpatialHarmonicInductionExposureLimits as Limits,
        SpatialHarmonicInductionExposureResult as Exposure,
    };
    use fs_fft::harmonic_window::{HarmonicProduct, HarmonicWindow, HarmonicWindowProfile};

    fn limits() -> Limits {
        Limits {
            maximum_lines: 8,
            maximum_pair_products: 1_000_000,
            maximum_workspace_bytes: 1_000_000,
        }
    }

    #[test]
    fn biased_multitone_impulse_matches_real_space_time_fields_and_interval_additivity() {
        let (problem, sources) = cube(1, None);
        let initial = problem
            .solve(
                &[[0.0; 3]; 6],
                &trace(&problem, [0.0; 3], false),
                static_controls(),
                &mut || false,
            )
            .unwrap();
        let state =
            SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
        let bias = problem
            .solve(
                &[[0.0; 3]; 6],
                &trace(&problem, [0.23, -0.17, 0.61], false),
                static_controls(),
                &mut || false,
            )
            .unwrap();
        let a = state
            .solve_harmonic(
                &[[C64::new(4.0, -1.2); 3]; 6],
                &[[0.0; 3]; 8],
                3.0,
                harmonic_controls(),
                &mut || false,
            )
            .unwrap();
        let b = state
            .solve_harmonic(
                &[[C64::new(2.5, 1.7); 3]; 6],
                &[[0.0; 3]; 8],
                4.0,
                harmonic_controls(),
                &mut || false,
            )
            .unwrap();
        let bank = state
            .prepare_harmonic_exposure(&[&b, &a], limits(), &mut || false)
            .unwrap();
        let window = HarmonicWindow::new(0.13, 0.47).unwrap();
        let observed = bank
            .observe_with_static_bias(&bias, window, &mut || false)
            .unwrap();
        let ac = bank.observe(window, &mut || false).unwrap();
        assert_eq!(observed.ohmic_heat_j.to_bits(), ac.ohmic_heat_j.to_bits());
        assert_eq!(observed.mechanical_work_j, 0.0);
        for (a, b) in observed.cells().iter().zip(ac.cells()) {
            assert_eq!(a.ohmic_heat_j.to_bits(), b.ohmic_heat_j.to_bits());
            assert!(std::ptr::eq(a.source, b.source));
        }
        let loads = observed.lorentz_loads(&ALL, LOADS, &mut || false).unwrap();
        let ac_loads = ac.lorentz_loads(&ALL, LOADS, &mut || false).unwrap();
        assert!(
            loads
                .nodal_forces_n()
                .iter()
                .zip(ac_loads.nodal_forces_n())
                .any(|(a, b)| (0..3).any(|axis| (a[axis] - b[axis]).abs() > 1e-3))
        );
        let mut nodes = [[0.0; 3]; 8];
        let mut wrench = [0.0; 6];
        let mut power = 0.0;
        // Composite Simpson, independent of the analytic temporal kernels.
        let steps = 1024;
        for step in 0..=steps {
            let time = window.start_s() + window.duration_s() * step as f64 / steps as f64;
            let weight = if step == 0 || step == steps {
                1.0
            } else if step % 2 == 0 {
                2.0
            } else {
                4.0
            } / (3 * steps) as f64;
            let (n, w, p) = instantaneous_reference(&[&a, &b], &bias, &ALL, time);
            for node in 0..8 {
                for axis in 0..3 {
                    nodes[node][axis] += weight * n[node][axis];
                }
            }
            for axis in 0..6 {
                wrench[axis] += weight * w[axis];
            }
            power += weight * p;
        }
        for (actual, expected) in loads.nodal_forces_n().iter().zip(nodes) {
            for axis in 0..3 {
                near(
                    window.duration_s() * actual[axis],
                    window.duration_s() * expected[axis],
                    2e-10,
                );
            }
        }
        let actual_wrench = loads.wrench(&Se3::identity(), [0.0; 3]).unwrap().to_array();
        for axis in 0..6 {
            near(actual_wrench[axis], wrench[axis], 2e-10);
        }
        let velocity: Vec<_> = problem
            .positions_m()
            .iter()
            .copied()
            .map(virtual_velocity)
            .collect();
        near(
            loads.power_w(&velocity, &mut || false).unwrap(),
            power,
            2e-10,
        );
        let left = bank
            .observe_with_static_bias(&bias, HarmonicWindow::new(0.13, 0.2).unwrap(), &mut || {
                false
            })
            .unwrap();
        let right = bank
            .observe_with_static_bias(&bias, HarmonicWindow::new(0.33, 0.27).unwrap(), &mut || {
                false
            })
            .unwrap();
        let left = left.lorentz_loads(&ALL, LOADS, &mut || false).unwrap();
        let right = right.lorentz_loads(&ALL, LOADS, &mut || false).unwrap();
        for node in 0..8 {
            for axis in 0..3 {
                near(
                    0.47 * loads.nodal_forces_n()[node][axis],
                    0.2 * left.nodal_forces_n()[node][axis]
                        + 0.27 * right.nodal_forces_n()[node][axis],
                    2e-12,
                );
            }
        }
        // One full carrier period removes the DC/AC mean, while a short
        // interval approaches the instantaneous real-field force.
        let single = state
            .prepare_harmonic_exposure(&[&a], limits(), &mut || false)
            .unwrap();
        let period = HarmonicWindow::new(0.17, std::f64::consts::TAU / 3.0).unwrap();
        let biased = single
            .observe_with_static_bias(&bias, period, &mut || false)
            .unwrap();
        let unbiased = single.observe(period, &mut || false).unwrap();
        for (a, b) in biased.cells().iter().zip(unbiased.cells()) {
            for vertex in 0..4 {
                for axis in 0..3 {
                    near(
                        a.mean_lorentz_force_density_n_m3[vertex][axis],
                        b.mean_lorentz_force_density_n_m3[vertex][axis],
                        2e-10,
                    );
                }
            }
        }
        let tiny = HarmonicWindow::new(0.23, 1e-10).unwrap();
        let observed = bank
            .observe_with_static_bias(&bias, tiny, &mut || false)
            .unwrap();
        let actual = observed.lorentz_loads(&ALL, LOADS, &mut || false).unwrap();
        let (expected, _, _) = instantaneous_reference(&[&a, &b], &bias, &ALL, 0.23 + 0.5e-10);
        for (a, b) in actual.nodal_forces_n().iter().zip(expected) {
            for axis in 0..3 {
                near(a[axis], b[axis], 2e-10);
            }
        }
        // Bias force scales linearly with a duplicated coherent carrier;
        // original AC self-force scales quadratically.
        let single_biased = single
            .observe_with_static_bias(&bias, window, &mut || false)
            .unwrap();
        let single_ac = single.observe(window, &mut || false).unwrap();
        let doubled = state
            .prepare_harmonic_exposure(&[&a, &a], limits(), &mut || false)
            .unwrap();
        let doubled_biased = doubled
            .observe_with_static_bias(&bias, window, &mut || false)
            .unwrap();
        let doubled_ac = doubled.observe(window, &mut || false).unwrap();
        for cell in 0..6 {
            for vertex in 0..4 {
                for axis in 0..3 {
                    let value = |r: &Exposure<'_, '_>| {
                        r.cells()[cell].mean_lorentz_force_density_n_m3[vertex][axis]
                    };
                    near(
                        value(&doubled_biased) - value(&doubled_ac),
                        2.0 * (value(&single_biased) - value(&single_ac)),
                        2e-12,
                    );
                }
            }
        }
    }

    #[test]
    fn biased_observation_retains_authored_work_memory_source_and_cancellation_bounds() {
        let (problem, sources) = cube(1, None);
        let initial = problem
            .solve(
                &[[0.0; 3]; 6],
                &trace(&problem, [0.0; 3], false),
                static_controls(),
                &mut || false,
            )
            .unwrap();
        let state =
            SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
        let bias = problem
            .solve(
                &[[0.0; 3]; 6],
                &trace(&problem, [0.2, -0.1, 0.8], false),
                static_controls(),
                &mut || false,
            )
            .unwrap();
        let current = [[C64::new(4.0, 1.0); 3]; 6];
        let a = state
            .solve_harmonic(
                &current,
                &[[0.0; 3]; 8],
                3.0,
                harmonic_controls(),
                &mut || false,
            )
            .unwrap();
        let window = HarmonicWindow::new(0.2, 0.7).unwrap();
        // Two authored lines coalesce, but cannot evade either work admission.
        let exact_work = 39 * 6 * 4 + 12 * 6 * 2;
        let exact = state
            .prepare_harmonic_exposure(
                &[&a, &a],
                Limits {
                    maximum_pair_products: exact_work,
                    ..limits()
                },
                &mut || false,
            )
            .unwrap();
        assert_eq!(exact.distinct_frequencies(), 1);
        assert!(
            exact
                .observe_with_static_bias(&bias, window, &mut || false)
                .is_ok()
        );
        let short = state
            .prepare_harmonic_exposure(
                &[&a, &a],
                Limits {
                    maximum_pair_products: exact_work - 1,
                    ..limits()
                },
                &mut || false,
            )
            .unwrap();
        assert!(short.observe(window, &mut || false).is_ok());
        assert!(matches!(
            short.observe_with_static_bias(&bias, window, &mut || false),
            Err(Error::Invalid(_))
        ));
        // Discover the native AC payload admission through the public API;
        // the bias must additionally admit every authored moment kernel.
        let (mut low, mut high) = (0, limits().maximum_workspace_bytes);
        while low < high {
            let middle = low + (high - low) / 2;
            let fits = state
                .prepare_harmonic_exposure(
                    &[&a, &a],
                    Limits {
                        maximum_workspace_bytes: middle,
                        ..limits()
                    },
                    &mut || false,
                )
                .is_ok();
            if fits {
                high = middle;
            } else {
                low = middle + 1;
            }
        }
        let exact_bytes = low + 2 * core::mem::size_of::<HarmonicProduct>();
        for (bytes, succeeds) in [(exact_bytes - 1, false), (exact_bytes, true)] {
            let bank = state
                .prepare_harmonic_exposure(
                    &[&a, &a],
                    Limits {
                        maximum_workspace_bytes: bytes,
                        ..limits()
                    },
                    &mut || false,
                )
                .unwrap();
            assert!(bank.observe(window, &mut || false).is_ok());
            assert_eq!(
                bank.observe_with_static_bias(&bias, window, &mut || false)
                    .is_ok(),
                succeeds
            );
        }
        let (foreign_problem, _) = cube(1, None);
        let foreign = foreign_problem
            .solve(
                &[[0.0; 3]; 6],
                &trace(&foreign_problem, [0.2, -0.1, 0.8], false),
                static_controls(),
                &mut || false,
            )
            .unwrap();
        assert!(
            exact
                .observe_with_static_bias(&foreign, window, &mut || false)
                .is_err()
        );
        let mut partial_trace = trace(&problem, [0.0; 3], false);
        partial_trace.potential.truncate(1);
        let different_support = problem
            .solve(
                &[[0.0; 3]; 6],
                &partial_trace,
                static_controls(),
                &mut || false,
            )
            .unwrap();
        assert!(
            exact
                .observe_with_static_bias(&different_support, window, &mut || false)
                .is_err()
        );
        let dc_driven = problem
            .solve(
                &[[1.7; 3]; 6],
                &trace(&problem, [0.2, -0.1, 0.8], false),
                static_controls(),
                &mut || false,
            )
            .unwrap();
        assert!(
            exact
                .observe_with_static_bias(&dc_driven, window, &mut || false)
                .is_err()
        );
        let moving = state
            .solve_harmonic(
                &current,
                &[[0.1, 0.0, 0.0]; 8],
                3.0,
                harmonic_controls(),
                &mut || false,
            )
            .unwrap();
        let moving = state
            .prepare_harmonic_exposure(&[&moving], limits(), &mut || false)
            .unwrap();
        assert!(
            moving
                .observe_with_static_bias(&bias, window, &mut || false)
                .is_err()
        );
        let shaped =
            HarmonicWindow::with_profile(0.2, 0.7, HarmonicWindowProfile::Triangular).unwrap();
        assert!(
            exact
                .observe_with_static_bias(&bias, shaped, &mut || false)
                .is_err()
        );
        let mut calls = 0;
        let reference = exact
            .observe_with_static_bias(&bias, window, &mut || {
                calls += 1;
                false
            })
            .unwrap();
        for stop in [1, calls / 2, calls] {
            let mut count = 0;
            assert!(matches!(
                exact.observe_with_static_bias(&bias, window, &mut || {
                    count += 1;
                    count == stop
                }),
                Err(Error::Cancelled)
            ));
        }
        let replay = exact
            .observe_with_static_bias(&bias, window, &mut || false)
            .unwrap();
        for (a, b) in reference.cells().iter().zip(replay.cells()) {
            assert_eq!(
                a.mean_lorentz_force_density_n_m3,
                b.mean_lorentz_force_density_n_m3
            );
        }
        assert_eq!(state.time_s(), 0.0);
        assert_eq!(state.potential_wb(), initial.potential_wb());
    }
}
