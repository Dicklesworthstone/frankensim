use super::super::vorticity::VorticityBudget3d;
use super::*;

fn evaluate(flow: &MixtureEuler3d) -> Vec<VorticityBudget3d> {
    let bytes = flow.required_vorticity_workspace_bytes().unwrap();
    let mut work = flow.vorticity_workspace(bytes).unwrap();
    flow.vorticity(&mut work, || false).unwrap().to_vec()
}

#[test]
fn g0_isothermal_composition_gradients_generate_all_baroclinic_components() {
    let grid = CartesianGrid::new(4, 5, 3, 0.3, 0.2, 0.4).unwrap();
    let boundaries = Boundaries3d {
        x: OPEN,
        y: OPEN,
        z: OPEN,
    };
    let rotation = [0.2, -0.1, 0.3];
    let expansion = 0.07;
    let linear = 0.4;
    // Each cell's temperature comes from a bounded energy inversion that stops once the
    // energy residual is within rho*abs + rel*|e| (fs-gas eos.rs decode_cell, with the
    // CaloricDomain::new defaults that `gas` uses), so |dT| <= abs/cv_min + rel*T. cv_min is
    // the smaller specific heat at T=2: species B, R_B*(3.5 + linear*T - 1) = 0.5*(3.5 + 0.8 - 1).
    // p = rho*(0.5 + 0.5*Y)*T with 0.5 + 0.5*Y <= 1, so |dp| <= rho_max*|dT|. The pressure
    // gradient on an open face is one-sided, |d(grad p)| <= 2*|dp|/h_min, and the torque is
    // grad(rho) x grad(p)/rho^2 with |grad rho| = 0.2, so |d torque| <= 0.2*|d(grad p)|/rho_min^2.
    let caloric = fs_gas::multicomponent::CaloricDomain::new(
        fs_qty::Temperature::new(0.001),
        fs_qty::Temperature::new(50.0),
    );
    let cv_min = 0.5 * (3.5 + linear * 2.0 - 1.0);
    let temperature_error =
        caloric.absolute_energy_tolerance / cv_min + caloric.relative_energy_tolerance * 2.0;
    let extent = [4.0 * 0.3, 5.0 * 0.2, 3.0 * 0.4];
    let rho_min = 1.2;
    let rho_max = 1.2 + 0.2 * extent.into_iter().fold(0.0_f64, f64::max);
    let h_min = grid.spacing.into_iter().fold(f64::INFINITY, f64::min);
    let torque_error = 0.2 * (2.0 * rho_max * temperature_error / h_min) / (rho_min * rho_min);
    // The two present species have R_A=1 and R_B=1/2. Holding T=2 while
    // varying total density and composition along independent axes makes
    // p = T*rho*(R_B + (R_A-R_B)*Y_A). Its cross gradient has a known sign.
    for [density_axis, fraction_axis, torque_axis] in [[0, 1, 2], [1, 2, 0], [2, 0, 1]] {
        for boost in [[0.0; 3], [0.7, -0.4, 0.2]] {
            let flow = field(standard(linear), grid, boundaries, |x| {
                let rho = 1.2 + 0.2 * x[density_axis];
                let fraction = 0.35 + 0.12 * x[fraction_axis];
                let swirl = [
                    rotation[1] * x[2] - rotation[2] * x[1],
                    rotation[2] * x[0] - rotation[0] * x[2],
                    rotation[0] * x[1] - rotation[1] * x[0],
                ];
                (
                    [rho * fraction, rho * (1.0 - fraction), 0.0],
                    std::array::from_fn(|d| boost[d] + expansion * x[d] + swirl[d]),
                    2.0 * rho * (0.5 + 0.5 * fraction),
                )
            });
            for (i, sample) in evaluate(&flow).iter().enumerate() {
                let x = grid.cell_center(i).unwrap();
                let rho = 1.2 + 0.2 * x[density_axis];
                let temperature = flow
                    .cell(i)
                    .unwrap()
                    .primitive(flow.gas(), || false)
                    .unwrap()
                    .temperature;
                assert!(
                    (temperature - 2.0).abs() <= temperature_error,
                    "{temperature:e} is further from 2 than the inversion tolerance {temperature_error:e}"
                );
                close(sample.divergence, 3.0 * expansion, 2e-12);
                for (d, omega) in rotation.map(|value| 2.0 * value).into_iter().enumerate() {
                    let torque = if d == torque_axis { 0.024 / rho } else { 0.0 };
                    close(sample.vorticity[d], omega, 2e-12);
                    close(sample.stretching[d], expansion * omega, 2e-12);
                    close(sample.dilation[d], -3.0 * expansion * omega, 2e-12);
                    assert!(
                        (sample.baroclinic[d] - torque).abs() <= torque_error,
                        "{:e} is further from {torque:e} than the propagated inversion tolerance {torque_error:e}",
                        sample.baroclinic[d]
                    );
                    // The stretching and dilation terms keep their own 2e-12 closeness; the rest is the torque's.
                    let expected_rate = -2.0 * expansion * omega + torque;
                    assert!(
                        (sample.inviscid_material_rate()[d] - expected_rate).abs()
                            <= torque_error + 2e-11 * (1.0 + expected_rate.abs()),
                        "{:e} is further from {expected_rate:e} than the propagated inversion tolerance",
                        sample.inviscid_material_rate()[d]
                    );
                }
            }
        }
    }
}

#[test]
fn g0_composition_gradient_adds_no_extra_torque_to_a_barotropic_field() {
    let grid = CartesianGrid::new(6, 5, 1, 0.2, 0.3, 0.7).unwrap();
    let flow = field(
        standard(0.6),
        grid,
        Boundaries3d {
            x: OPEN,
            y: OPEN,
            z: AxisBoundaries::Periodic,
        },
        |x| {
            let rho = 1.0 + 0.2 * x[0];
            let fraction = 0.3 + 0.1 * x[1];
            (
                [rho * fraction, rho * (1.0 - fraction), 0.0],
                [0.3, -0.2, 0.1],
                2.0 * rho,
            )
        },
    );
    // Composition, temperature and heat capacity vary, but p=2*rho remains
    // barotropic. Species gradients must not be counted twice in momentum.
    for sample in evaluate(&flow) {
        for value in sample
            .vorticity
            .into_iter()
            .chain(sample.inviscid_material_rate())
        {
            close(value, 0.0, 3e-11);
        }
    }
}

#[test]
fn g1_native_mixture_step_resolves_composition_driven_baroclinic_onset() {
    let mut errors = Vec::new();
    for n in [12, 24, 48] {
        let grid = CartesianGrid::new(n, n, 1, 2.0 / n as f64, 1.3 / n as f64, 0.7).unwrap();
        let kx = std::f64::consts::PI;
        let ky = 2.0 * std::f64::consts::PI / 1.3;
        let mut flow = field(standard(0.3), grid, PERIODIC, |x| {
            let rho = 1.0 + 0.12 * (kx * x[0]).sin();
            let fraction = 0.4 + 0.15 * (ky * x[1]).sin();
            (
                [rho * fraction, rho * (1.0 - fraction), 0.0],
                [0.0; 3],
                2.0 * rho * (0.5 + 0.5 * fraction),
            )
        });
        let mut work = flow
            .vorticity_workspace(flow.required_vorticity_workspace_bytes().unwrap())
            .unwrap();
        let before = flow.vorticity(&mut work, || false).unwrap().to_vec();
        // With zero initial velocity the material derivative is the time
        // derivative. Only the original conservative species Euler step acts.
        let dt = 2e-7;
        let step = flow
            .step(Time::new(dt), StepControl::default(), || false)
            .unwrap();
        close(step.elapsed.value(), dt, 1e-15);
        assert!(!step.first_order_retry);
        let after = flow.vorticity(&mut work, || false).unwrap();
        let mut error = 0.0;
        let mut norm = 0.0;
        let mut correlation = 0.0;
        for (initial, final_state) in before.iter().zip(after) {
            let observed = (final_state.vorticity[2] - initial.vorticity[2]) / dt;
            let expected = initial.baroclinic[2];
            error += (observed - expected).abs();
            norm += expected.abs();
            correlation += observed * expected;
        }
        assert!(
            correlation > 0.0,
            "right-handed composition-driven production"
        );
        errors.push(error / norm);
    }
    // MC-limited fluxes lose pointwise order at smooth extrema; integrated
    // refinement tests the actual native update without asserting unlimited
    // reconstruction there.
    assert!(
        errors[1] < 0.5 * errors[0] && errors[2] < 0.5 * errors[1],
        "{errors:?}"
    );
    assert!(
        errors[2] < 0.03,
        "native mixture onset versus budget: {errors:?}"
    );
}

#[test]
fn g4_mixture_vorticity_cap_cancellation_and_grid_refusal_preserve_state() {
    let grid = CartesianGrid::new(4, 3, 1, 0.2, 0.3, 0.7).unwrap();
    let flow = field(standard(0.4), grid, PERIODIC, |_| {
        ([0.4, 0.6, 0.0], [0.2, -0.1, 0.3], 2.0)
    });
    let bytes = flow.required_vorticity_workspace_bytes().unwrap();
    assert!(flow.vorticity_workspace(bytes - 1).is_err());
    let mut work = flow.vorticity_workspace(bytes).unwrap();
    let capacity = work.workspace_bytes();
    let mut polls = 0;
    let expected = flow
        .vorticity(&mut work, || {
            polls += 1;
            false
        })
        .unwrap()
        .to_vec();
    let accepted = flow.states().to_vec();
    let time = flow.time();
    let report = flow.last_step().copied();
    for cancel_at in [1, polls / 2, polls] {
        let mut seen = 0;
        assert_eq!(
            flow.vorticity(&mut work, || {
                seen += 1;
                seen == cancel_at
            }),
            Err(Error::Cancelled)
        );
        assert_eq!(flow.states(), accepted);
        assert_eq!(flow.time(), time);
        assert_eq!(flow.last_step().copied(), report);
        assert_eq!(flow.vorticity(&mut work, || false).unwrap(), expected);
        assert_eq!(work.workspace_bytes(), capacity);
    }
    let other = field(
        standard(0.4),
        CartesianGrid::new(3, 4, 1, 0.2, 0.3, 0.7).unwrap(),
        PERIODIC,
        |_| ([0.4, 0.6, 0.0], [0.0; 3], 2.0),
    );
    assert!(other.vorticity(&mut work, || false).is_err());
    let unresolved = field(
        standard(0.4),
        CartesianGrid::new(2, 3, 1, 0.2, 0.3, 0.7).unwrap(),
        PERIODIC,
        |_| ([0.4, 0.6, 0.0], [0.0; 3], 2.0),
    );
    assert!(unresolved.vorticity_workspace(usize::MAX).is_err());
}
