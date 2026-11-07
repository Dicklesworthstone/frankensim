//! G1/G3/G4: current-area pressure, retained thermal history and original inertia.
use super::*;
use fs_solid::hyper3::plastic::FinitePlasticTetError;
use fs_solid::hyper3::plastic::porous::thermal::FiniteThermalPorousTetError;
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};

fn at_pressure(pressure: f64, temperature: f64, old: Hydro) -> Hydro {
    // Cauchy pressure supplies P=-p*s², not a constant nominal stress.
    let residual = |log_stretch| {
        let state = Hydro::at(log_stretch, temperature, Some(old));
        state.piola + pressure * state.stretch.powi(2)
    };
    let (mut lo, mut hi) = (-0.35, 0.2);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Hydro::at(0.5 * (lo + hi), temperature, Some(old))
}

#[test]
fn pressure_heating_and_release_follow_original_thermal_hydrostatics() {
    with_cx(|cx, _| {
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let mut solid = initial(cx);
        let mut old = Hydro::at(0.0, 300.0, None);
        let target = Hydro::at(-0.025, 300.0, Some(old));
        let pressure = -target.piola / target.stretch.powi(2);
        let supports: Vec<_> = (0..12)
            .filter(|i| ![3, 7, 11].contains(i))
            .map(|i| (i, 0.0))
            .collect();
        let mut hot_alpha = 0.0;
        for (index, (temperature, pressure)) in [
            (300.0, pressure),
            (450.0, pressure),
            (450.0, 0.0),
            (300.0, 0.0),
        ]
        .into_iter()
        .enumerate()
        {
            let expected = at_pressure(pressure, temperature, old);
            let step = solid
                .advance_with_pressure(
                    &[temperature],
                    &supports,
                    &[],
                    [0.0; 3],
                    0.2,
                    &surface,
                    pressure,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            let account = step.pressure.unwrap();
            let work = -pressure * V * (expected.stretch.powi(3) - old.stretch.powi(3));
            let endpoint =
                -3.0 * pressure * V * expected.stretch.powi(2) * (expected.stretch - old.stretch);
            let [_, reversible, defect, coefficients] = expected.energy(old);
            close(step.solid.element_states()[0].alpha, expected.alpha, 2e-9);
            close(
                step.solid.element_jacobians()[0],
                expected.stretch.powi(3),
                2e-9,
            );
            close(
                step.solid.element_plastic_jacobians()[0],
                expected.log_jp.exp(),
                2e-9,
            );
            for i in [0, 4, 8] {
                close(
                    step.solid.element_first_piola_pa()[0][i],
                    expected.piola,
                    2e-6,
                );
            }
            close(
                step.energy.internal_energy_change_j,
                V * (expected.internal - old.internal),
                2e-6,
            );
            close(step.energy.plastic_dissipation_j, V * expected.heat, 2e-7);
            close(step.energy.reversible_heat_j, V * reversible, 2e-6);
            close(
                step.energy.caloric_heat_j,
                V * (expected.heat + reversible),
                2e-6,
            );
            close(step.energy.endpoint_work_defect_j, V * defect, 2e-7);
            close(
                step.energy.coefficient_storage_change_j,
                V * coefficients,
                2e-7,
            );
            close(account.work_j, work, 2e-7);
            close(account.endpoint_work_j, endpoint, 2e-7);
            close(account.endpoint_work_correction_j, work - endpoint, 2e-7);
            close(
                account.external_work_defect_j,
                V * defect + work - endpoint,
                2e-7,
            );
            close(step.energy.equilibrium_work_residual_j, 0.0, 1e-6);
            close(step.energy.balance_error_j, 0.0, 1e-6);
            close(
                account.work_j
                    + step.energy.dead_load_endpoint_work_j
                    + step.energy.support_endpoint_work_j,
                step.energy.internal_energy_change_j
                    + step.energy.caloric_heat_j
                    + account.external_work_defect_j,
                1e-6,
            );
            assert_eq!(step.solid.families(), solid.families());
            assert_eq!(step.solid.temperatures_k(), &[temperature]);
            assert!(step.inertia.is_none());
            if index == 1 {
                // Thermal expansion and pressure-sensitive plastic volume are
                // different coordinates: heating increases retained compaction.
                assert!(expected.alpha > old.alpha);
                assert!(expected.log_jp < old.log_jp);
                hot_alpha = expected.alpha;
            }
            if index >= 2 {
                close(step.energy.plastic_dissipation_j, 0.0, 1e-10);
                close(step.solid.element_states()[0].alpha, hot_alpha, 2e-9);
            }
            old = expected;
            solid = step.solid;
        }
        assert!(solid.element_jacobians()[0] < 0.9);
        close(
            solid.element_jacobians()[0],
            solid.element_plastic_jacobians()[0],
            2e-9,
        );
    });
}

// Centered regular tetrahedron: grad(N_a)=X_a/4, V0=8/3 and M*X=rho*V0*X/20.
const REGULAR: [[f64; 3]; 4] = [
    [1.0, 1.0, 1.0],
    [1.0, -1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [-1.0, 1.0, -1.0],
];
const REGULAR_VOLUME: f64 = 8.0 / 3.0;

#[test]
fn current_area_pressure_drives_free_thermal_compaction_with_consistent_inertia() {
    with_cx(|cx, _| {
        let rho = 2.0;
        let dt = 0.2;
        let stretch = 0.94_f64;
        let old = Hydro::at(0.0, 300.0, None);
        let expected = Hydro::at(stretch.ln(), 350.0, Some(old));
        // P+p*s²+rho*(s-1)/(5*dt²)=0 follows from the original consistent mass.
        let pressure = -expected.piola / stretch.powi(2)
            + rho * (1.0 - stretch) / (5.0 * dt * dt * stretch.powi(2));
        let solid = FiniteThermalPorousTetSolid::try_new(
            &REGULAR,
            &TETS,
            vec![family()],
            &[300.0],
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap()
        .with_inertia(&[rho], &[[0.0; 3]; 4], cx)
        .unwrap();
        let surface =
            ClosedTetPressureSurface::new(&REGULAR, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let step = solid
            .advance_with_pressure(
                &[350.0],
                &[],
                &[],
                [0.0; 3],
                dt,
                &surface,
                pressure,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        assert!(step.equilibrium.reactions_n.is_empty());
        assert!(step.equilibrium.steps[0].linear_iterations > 0);
        for (node, x) in REGULAR.iter().enumerate() {
            for (axis, coordinate) in x.iter().enumerate() {
                close(
                    step.solid.displacement_m()[node][axis],
                    (stretch - 1.0) * coordinate,
                    3e-8,
                );
                close(
                    step.solid.velocity_mps().unwrap()[node][axis],
                    (stretch - 1.0) * coordinate / dt,
                    2e-7,
                );
            }
        }
        let kinetic = 3.0 * rho * REGULAR_VOLUME * (stretch - 1.0).powi(2) / (10.0 * dt * dt);
        let inertia = step.inertia.as_ref().unwrap();
        close(inertia.kinetic_change_j, kinetic, 2e-7);
        close(inertia.numerical_loss_j, kinetic, 2e-7);
        close(inertia.inertial_work_j, 2.0 * kinetic, 3e-7);
        close(
            step.solid.reference_mass_kg().unwrap(),
            rho * REGULAR_VOLUME,
            1e-12,
        );
        close(
            step.solid.current_density_kg_m3(&[rho], cx).unwrap()[0],
            rho / stretch.powi(3),
            2e-7,
        );
        close(step.solid.element_states()[0].alpha, expected.alpha, 3e-8);
        close(
            step.solid.element_plastic_jacobians()[0],
            expected.log_jp.exp(),
            3e-8,
        );
        let [_, reversible, defect, _] = expected.energy(old);
        close(
            step.energy.caloric_heat_j,
            REGULAR_VOLUME * (expected.heat + reversible),
            5e-6,
        );
        close(
            step.energy.endpoint_work_defect_j,
            REGULAR_VOLUME * defect,
            2e-6,
        );
        let account = step.pressure.unwrap();
        close(
            account.work_j,
            -pressure * REGULAR_VOLUME * (stretch.powi(3) - 1.0),
            3e-6,
        );
        close(
            account.endpoint_work_j,
            -3.0 * pressure * REGULAR_VOLUME * stretch.powi(2) * (stretch - 1.0),
            3e-6,
        );
        close(
            account.work_j,
            step.energy.internal_energy_change_j
                + step.energy.caloric_heat_j
                + account.external_work_defect_j
                + inertia.kinetic_change_j
                + inertia.numerical_loss_j,
            2e-6,
        );
        close(step.energy.equilibrium_work_residual_j, 0.0, 1e-6);
        close(step.energy.balance_error_j, 0.0, 1e-6);
        assert_eq!(solid.displacement_m(), &[[0.0; 3]; 4]);
        assert_eq!(solid.velocity_mps().unwrap(), &[[0.0; 3]; 4]);
        assert_eq!(solid.temperatures_k(), &[300.0]);
    });
}

#[test]
fn pressure_patch_keeps_graded_support_work_heat_and_atomic_admission_separate() {
    with_cx(|cx, gate| {
        let rho = [2.0, 4.0, 6.0, 8.0];
        let solid = initial(cx)
            .with_nodal_inertia(&rho, &[[0.0; 3]; 4], cx)
            .unwrap();
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let stretch = [0.98, 0.99, 0.97];
        let translation = [0.015, -0.02, 0.01];
        let target: Vec<_> = NODES
            .iter()
            .enumerate()
            .flat_map(|(node, x)| {
                (0..3).map(move |axis| {
                    (
                        3 * node + axis,
                        (stretch[axis] - 1.0) * x[axis] + translation[axis],
                    )
                })
            })
            .collect();
        let p = 300.0;
        let mut faces = vec![0.0; surface.faces().len()];
        let patch = surface
            .faces()
            .iter()
            .position(|face| !face.contains(&0))
            .unwrap();
        faces[patch] = p;
        let field = TetPressure::PerFace(&faces);
        let forces = [(3, -1.0)];
        let dt = 0.2;
        let baseline = solid
            .advance(
                &[330.0],
                &target,
                &forces,
                [0.0; 3],
                dt,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let step = solid
            .advance_with_pressure(
                &[330.0],
                &target,
                &forces,
                [0.0; 3],
                dt,
                &surface,
                field,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        // A fixed deformation/temperature has identical constitutive heat,
        // regardless of how external load and support reactions share its work.
        assert_eq!(step.element_energy, baseline.element_energy);
        assert_eq!(step.energy.caloric_heat_j, baseline.energy.caloric_heat_j);
        assert_eq!(
            step.energy.endpoint_work_defect_j,
            baseline.energy.endpoint_work_defect_j
        );
        assert_eq!(step.inertia, baseline.inertia);
        assert_eq!(step.solid.element_states(), baseline.solid.element_states());
        assert_eq!(
            step.solid.reference_nodal_density_kg_m3(),
            Some(rho.as_slice())
        );
        let delta = stretch.map(|s| s - 1.0);
        let mut area_end = [0.0; 3];
        let mut area_integral = [0.0; 3];
        for axis in 0..3 {
            let j = (axis + 1) % 3;
            let k = (axis + 2) % 3;
            area_end[axis] = 0.5 * stretch[j] * stretch[k];
            area_integral[axis] =
                0.5 * (1.0 + 0.5 * (delta[j] + delta[k]) + delta[j] * delta[k] / 3.0);
        }
        // Integrate the original triangle's quadratic current area along the
        // declared straight nodal path, including its nonzero translation work.
        let endpoint = -p / 3.0
            * (0..3)
                .map(|i| area_end[i] * (delta[i] + 3.0 * translation[i]))
                .sum::<f64>();
        let exact = -p / 3.0
            * (0..3)
                .map(|i| area_integral[i] * (delta[i] + 3.0 * translation[i]))
                .sum::<f64>();
        let account = step.pressure.unwrap();
        close(account.endpoint_work_j, endpoint, 2e-10);
        close(account.work_j, exact, 2e-10);
        close(account.endpoint_work_correction_j, exact - endpoint, 2e-10);
        assert!(account.endpoint_work_correction_j.abs() > 1e-3);
        for &(dof, reaction) in &step.equilibrium.reactions_n {
            let base = baseline
                .equilibrium
                .reactions_n
                .iter()
                .find(|&&(i, _)| i == dof)
                .unwrap()
                .1;
            let pressure_force = if dof / 3 == 0 {
                0.0
            } else {
                -p * area_end[dof % 3] / 3.0
            };
            close(reaction, base - pressure_force, 2e-7);
        }
        let inertia = step.inertia.as_ref().unwrap();
        close(
            step.energy.dead_load_endpoint_work_j + step.energy.support_endpoint_work_j + exact,
            step.energy.internal_energy_change_j
                + step.energy.caloric_heat_j
                + account.external_work_defect_j
                + inertia.kinetic_change_j
                + inertia.numerical_loss_j,
            1e-6,
        );
        close(step.energy.balance_error_j, 0.0, 1e-6);
        // Each correction fits separately; their sum does not. Opposite signed
        // contributions cannot hide the need for a smaller physical increment.
        let a = step.energy.endpoint_work_defect_j.abs();
        let b = account.endpoint_work_correction_j.abs();
        assert!(a > 0.0 && b > 0.0);
        let mut strict = controls();
        strict.material.maximum_endpoint_work_defect_j = a.max(b) + 0.5 * a.min(b);
        assert!(matches!(
            solid.advance_with_pressure(
                &[330.0],
                &target,
                &forces,
                [0.0; 3],
                dt,
                &surface,
                field,
                settings(),
                strict,
                cx,
            ),
            Err(FiniteThermalPorousTetError::Solid(
                FinitePlasticTetError::EnergyGate {
                    quantity: "thermal porous pressure/constitutive endpoint integration defect",
                    ..
                }
            ))
        ));
        strict = controls();
        strict.maximum_numerical_kinetic_loss_j = 0.0;
        assert!(
            solid
                .advance_with_pressure(
                    &[330.0],
                    &target,
                    &forces,
                    [0.0; 3],
                    dt,
                    &surface,
                    field,
                    settings(),
                    strict,
                    cx,
                )
                .is_err()
        );
        let retry = solid
            .advance_with_pressure(
                &[330.0],
                &target,
                &forces,
                [0.0; 3],
                dt,
                &surface,
                field,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(retry.equilibrium, step.equilibrium);
        assert_eq!(retry.energy, step.energy);
        assert_eq!(retry.pressure, step.pressure);
        assert!(
            solid
                .advance_with_pressure(
                    &[551.0],
                    &target,
                    &forces,
                    [0.0; 3],
                    dt,
                    &surface,
                    field,
                    settings(),
                    controls(),
                    cx,
                )
                .is_err()
        );
        let shifted = NODES.map(|x| [x[0] + 1.0, x[1], x[2]]);
        let wrong_surface =
            ClosedTetPressureSurface::new(&shifted, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        assert!(
            solid
                .advance_with_pressure(
                    &[330.0],
                    &target,
                    &forces,
                    [0.0; 3],
                    dt,
                    &wrong_surface,
                    field,
                    settings(),
                    controls(),
                    cx,
                )
                .is_err()
        );
        assert!(
            solid
                .advance_with_pressure(
                    &[330.0],
                    &target,
                    &forces,
                    [0.0; 3],
                    dt,
                    &surface,
                    TetPressure::PerFace(&[p]),
                    settings(),
                    controls(),
                    cx,
                )
                .is_err()
        );
        let zero = solid
            .advance_with_pressure(
                &[330.0],
                &target,
                &forces,
                [0.0; 3],
                dt,
                &surface,
                0.0,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(zero.solid.element_states(), baseline.solid.element_states());
        assert_eq!(zero.energy.caloric_heat_j, baseline.energy.caloric_heat_j);
        assert_eq!(zero.pressure.unwrap().work_j, 0.0);
        assert!(baseline.pressure.is_none());
        gate.request();
        assert!(
            solid
                .advance_with_pressure(
                    &[330.0],
                    &target,
                    &forces,
                    [0.0; 3],
                    dt,
                    &surface,
                    field,
                    settings(),
                    controls(),
                    cx,
                )
                .is_err()
        );
        assert_eq!(solid.displacement_m(), &[[0.0; 3]; 4]);
        assert_eq!(solid.velocity_mps().unwrap(), &[[0.0; 3]; 4]);
        assert_eq!(solid.temperatures_k(), &[300.0]);
        assert_eq!(solid.element_states()[0].alpha, 0.0);
        assert_eq!(solid.reference_nodal_density_kg_m3(), Some(rho.as_slice()));
    });
}
