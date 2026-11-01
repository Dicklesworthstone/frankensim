//! Source stiffness changes participate in finite heat and native equilibrium.
//! Synthetic material cards and analytical implementation oracles; no measured
//! material validation or exact finite-step entropy conservation is claimed.

mod laminate_thermal_common;
use laminate_thermal_common::*;

use fs_couple::laminate_thermal::LayeredThermalShell;
use fs_couple::laminate_thermal::reciprocal::{
    LaminateThermoelasticControls, LaminateThermoelasticError, ReciprocalLaminateState,
    ReciprocalLayeredThermalShell,
};
use fs_material::laminate::ReferenceSurface;
use fs_material::thermal::transport::stack::transient::implicit::{
    LayeredThermalSolveControls, LayeredThermalWorkspace,
};
use fs_plate::shell::mitc::strain::integration_rule;
use fs_plate::shell::mitc::{ConstraintFrame, ShellConstraint};

const T0: f64 = 300.0;
const ALPHA: f64 = 1e-5;
const SLOPE: f64 = -0.002 * E;

fn controls() -> LaminateThermoelasticControls {
    LaminateThermoelasticControls {
        maximum_iterations: 64,
        heat_relaxation: 1.0,
        cell_heat_tolerance_j: 1e-8,
        energy_tolerance_j: 1e-7,
    }
}

fn field_internal(state: &ReciprocalLaminateState) -> f64 {
    state
        .cell_potentials()
        .iter()
        .map(|p| p.internal_energy_j)
        .sum::<f64>()
        + state.transverse_shear_internal_energy_j()
        + state.drilling_gauge_energy_j()
}

fn scalar_root(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..96 {
        let mid = lo + 0.5 * (hi - lo);
        if residual(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    lo + 0.5 * (hi - lo)
}

fn axial_load(stress: f64) -> [f64; 24] {
    let mut load = [0.0; 24];
    let force = stress * H * 0.6 / 2.0;
    for (node, sign) in [(0, -1.0), (1, 1.0), (2, -1.0), (3, 1.0)] {
        load[6 * node] = sign * force;
    }
    load
}

#[test]
fn g1_restrained_heating_uses_the_affine_modulus_cubic_and_updated_prestress() {
    let cards = [card_with_coefficients(
        "restrained softening sheet",
        SLOPE,
        ALPHA,
        0.0,
        CP,
    )];
    let body = prepare_body(&cards, &[1], 1.0);
    let heat_laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&heat_laws, ids()).unwrap();
    let heater = 6000.0;
    let thermal = initial_with_sources(&bound, &[T0], &[heater]);
    let mut mesh = discretization();
    for node in 1..4 {
        mesh.constraints.push(ShellConstraint {
            node,
            mask: [true, true, false, false, false, true],
            frame: ConstraintFrame::Global,
        });
    }
    let shell = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        mesh,
        limits(),
        || false,
    )
    .unwrap();
    let mut coupled = ReciprocalLayeredThermalShell::try_new(
        &shell,
        &thermal,
        &[0.0; 24],
        &[[280.0, 360.0]],
        0.02,
        controls(),
        || false,
    )
    .unwrap();
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    let step = coupled
        .advance_to(
            &thermal,
            1.0,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || false,
        )
        .unwrap();

    let volume = AREA * H;
    let capacity = volume * RHO * CP;
    // With both normal strains fixed and nu=0, Psi=V E(T) alpha² x².
    // Differentiate that scalar polynomial at fixed total strain, then U=Psi+TS.
    // This gives a cubic balance, independent of the source AD and shell solver.
    let internal = |t: f64| {
        let x = t - T0;
        let young = E + SLOPE * x;
        volume * ALPHA * ALPHA * (young * x * x - t * (SLOPE * x * x + 2.0 * young * x))
    };
    let expected = scalar_root(T0, 320.0, |t| capacity * (t - T0) + internal(t) - heater);
    let actual = coupled.state().shell().temperatures_k()[0];
    close(actual, expected, 3e-8);
    close(
        coupled.state().cell_potentials()[0].internal_energy_j,
        internal(actual),
        1e-7,
    );
    close(step.report.external_mechanical_work_j, 0.0, 1e-12);
    close(step.report.external_heat_j, heater, 1e-7);
    close(step.report.combined_energy_residual_j, 0.0, 1e-7);
    close(step.thermal.mass_kg(), RHO * volume, 1e-12);
    assert!(actual > T0 + heater / capacity);

    // The accepted thermal state must reach the native membrane pencil with
    // current E(T), not the originally factored modulus. A unit corner w on
    // this 1 x .6 patch has integral |grad(w)|² = .3*(1+1/.6²).
    let prestress = coupled
        .state()
        .shell()
        .membrane_prestress(
            fs_plate::shell::mitc::thermal::prestress::MembranePrestressLimits {
                maximum_displacement_gradient: 0.05,
                maximum_rotation_rad: 0.05,
            },
            || false,
        )
        .unwrap();
    let mut corner = [0.0; 24];
    corner[20] = 1.0;
    let direction = prestress.project_load(&corner).unwrap();
    let geometric_work: f64 = direction
        .iter()
        .enumerate()
        .map(|(row, &u)| {
            let (columns, values) = prestress.geometric_stiffness().row(row);
            u * columns
                .iter()
                .zip(values)
                .map(|(&col, &value)| value * direction[col])
                .sum::<f64>()
        })
        .sum();
    let current_young = E + SLOPE * (actual - T0);
    close(
        geometric_work,
        -current_young * H * ALPHA * (actual - T0) * 0.3 * (1.0 + 1.0 / 0.6_f64.powi(2)),
        1e-7,
    );
    assert!(geometric_work < 0.0);
}

#[test]
fn g1_force_controlled_heating_resolves_softening_and_its_entropy_work() {
    // Zero expansion isolates the physical entropy due to changing modulus.
    let cards = [card_with_coefficients(
        "loaded softening sheet",
        SLOPE,
        0.0,
        0.0,
        CP,
    )];
    let body = prepare_body(&cards, &[1], 1.0);
    let heat_laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&heat_laws, ids()).unwrap();
    let heater = 12000.0;
    let thermal = initial_with_sources(&bound, &[T0], &[heater]);
    let shell = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        discretization(),
        limits(),
        || false,
    )
    .unwrap();
    let stress = 1e6;
    let load = axial_load(stress);
    let mut coupled = ReciprocalLayeredThermalShell::try_new(
        &shell,
        &thermal,
        &load,
        &[[280.0, 360.0]],
        0.02,
        controls(),
        || false,
    )
    .unwrap();
    let initial_internal = field_internal(coupled.state());
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    let step = coupled
        .advance_to(
            &thermal,
            1.0,
            &load,
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || false,
        )
        .unwrap();
    let volume = AREA * H;
    let capacity = volume * RHO * CP;
    // For imposed sigma and alpha=nu=0: e=sigma/E, Psi=V sigma²/(2E),
    // S=-V E' sigma²/(2E²); constant-traction work is V sigma² delta(1/E).
    let internal = |t: f64| {
        let young = E + SLOPE * (t - T0);
        0.5 * volume * stress * stress * (1.0 / young - t * SLOPE / (young * young))
    };
    let work = |t: f64| volume * stress * stress * (1.0 / (E + SLOPE * (t - T0)) - 1.0 / E);
    let expected = scalar_root(T0, 320.0, |t| {
        capacity * (t - T0) + internal(t) - internal(T0) - work(t) - heater
    });
    let actual = coupled.state().shell().temperatures_k()[0];
    let current_young = E + SLOPE * (actual - T0);
    close(actual, expected, 3e-8);
    close(
        coupled.state().shell().response().displacement_6n[6],
        stress / current_young,
        2e-11,
    );
    assert!(coupled.state().shell().response().displacement_6n[6] > stress / E);
    assert!(actual < T0 + heater / capacity);
    close(
        coupled.state().cell_potentials()[0].internal_energy_j,
        internal(actual),
        1e-7,
    );
    close(step.report.external_mechanical_work_j, work(actual), 1e-7);
    close(
        step.report.cells[0].stiffness_temperature_work_j,
        0.5 * volume * (stress / E) * (current_young - E) * (stress / current_young),
        1e-7,
    );
    assert!(step.report.cells[0].stiffness_temperature_work_j < 0.0);
    close(step.report.drilling_gauge_parameter_work_j, 0.0, 1e-7);
    close(
        step.thermal.energy_j() + field_internal(coupled.state())
            - thermal.energy_j()
            - initial_internal,
        heater + work(actual),
        1e-7,
    );
    close(step.report.combined_energy_residual_j, 0.0, 1e-7);
}

// Independent scalar sectional oracle: bending stress produces a shear stress
// resultant integral through the thickness; equating that shear energy with
// Q²/(2H_s) gives H_s=D²/integral(r²/G dz). Three-point Gauss is exact for the
// degree-four r² in each layer. This does not reuse the native polynomial or AD.
fn shear_modulus(temperatures: &[f64; 2], slopes: &[f64; 2]) -> f64 {
    let young = core::array::from_fn::<_, 2, _>(|i| E + slopes[i] * (temperatures[i] - T0));
    let h = H / 2.0;
    let middle = [-H / 4.0, H / 4.0];
    let neutral = (young[0] * middle[0] + young[1] * middle[1]) / (young[0] + young[1]);
    let bending: f64 = (0..2)
        .map(|i| young[i] * (h.powi(3) / 12.0 + h * (middle[i] - neutral).powi(2)))
        .sum();
    let gauss = [
        (-(3.0_f64 / 5.0).sqrt(), 5.0 / 9.0),
        (0.0, 8.0 / 9.0),
        ((3.0_f64 / 5.0).sqrt(), 5.0 / 9.0),
    ];
    let mut running = 0.0;
    let mut energy = 0.0;
    for i in 0..2 {
        let bottom = middle[i] - h / 2.0;
        for (point, weight) in gauss {
            let z = middle[i] + h * point / 2.0;
            let resultant =
                running + 0.5 * young[i] * ((z - neutral).powi(2) - (bottom - neutral).powi(2));
            energy += 0.5 * h * weight * resultant * resultant / (0.5 * young[i]);
        }
        running += young[i] * h * (middle[i] - neutral);
    }
    bending * bending / energy
}

fn shear_gradient(temperatures: &[f64; 2], slopes: &[f64; 2], cell: usize) -> f64 {
    // Test-only fourth-order difference, with sources evaluated by the scalar
    // mechanics oracle above; production uses exact source derivatives.
    let step = 0.02;
    let at = |offset: f64| {
        let mut t = *temperatures;
        t[cell] += offset * step;
        shear_modulus(&t, slopes)
    };
    (at(-2.0) - 8.0 * at(-1.0) + 8.0 * at(1.0) - at(2.0)) / (12.0 * step)
}

fn shear_product(a: &ReciprocalLaminateState, b: &ReciprocalLaminateState) -> f64 {
    a.shell()
        .strain_fields()
        .iter()
        .zip(b.shell().strain_fields())
        .map(|(a, b)| {
            integration_rule()
                .iter()
                .map(|&(lambda, weight)| {
                    let a = a.at_barycentric(lambda).unwrap().transverse_shear;
                    let b = b.at_barycentric(lambda).unwrap().transverse_shear;
                    0.3 * weight * (a[0] * b[0] + a[1] * b[1])
                })
                .sum::<f64>()
        })
        .sum()
}

#[test]
fn g1_nonuniform_softening_owns_nonlocal_shear_entropy_and_preserves_late_cancellation() {
    let slopes = [SLOPE, -0.005 * E];
    let cards = [
        card_with_coefficients("softening lower", slopes[0], ALPHA, 0.0, CP),
        card_with_coefficients("softening upper", slopes[1], 2.0 * ALPHA, 0.0, CP),
    ];
    let body = prepare_body(&cards, &[1, 1], 1.0);
    let heat_laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&heat_laws, ids()).unwrap();
    let temperatures = [310.0, 340.0];
    let thermal = initial_with_sources(&bound, &temperatures, &[1000.0, 1000.0]);
    let original_cells = thermal.cells().clone();
    let shell = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        discretization(),
        limits(),
        || false,
    )
    .unwrap();
    let mut load = [0.0; 24];
    load[20] = -1.0;
    let ranges = [[300.0, 350.0]; 2];
    // The hottest declared free strain is 2e-5*50 = .001; the 1 N
    // transverse load adds small bending strains. Retain a .005 whole-field
    // envelope, which the native Bernstein/layer-face gate must still verify,
    // instead of paying for an unused 2% strain box in the capacity bound.
    let fresh = || {
        ReciprocalLayeredThermalShell::try_new(
            &shell,
            &thermal,
            &load,
            &ranges,
            0.005,
            controls(),
            || false,
        )
        .unwrap()
    };
    let initial_owner = fresh();
    let initial_field = field_internal(initial_owner.state());
    let old_g2 = shear_product(initial_owner.state(), initial_owner.state());
    let old_entropy = core::array::from_fn::<_, 2, _>(|i| {
        -0.5 * old_g2 * shear_gradient(&temperatures, &slopes, i)
    });
    let mut coupled = fresh();
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    let mut checkpoints = 0;
    let step = coupled
        .advance_to(
            &thermal,
            0.5,
            &load,
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    let current: [f64; 2] = coupled.state().shell().temperatures_k().try_into().unwrap();
    let g2 = shear_product(coupled.state(), coupled.state());
    let cross = shear_product(initial_owner.state(), coupled.state());
    let entropy =
        core::array::from_fn::<_, 2, _>(|i| -0.5 * g2 * shear_gradient(&current, &slopes, i));
    let free = 0.5 * g2 * shear_modulus(&current, &slopes);
    close(coupled.state().transverse_shear_energy_j(), free, 1e-10);
    assert!(free > 1e-12);
    close(
        coupled.state().transverse_shear_internal_energy_j(),
        free + current.iter().zip(entropy).map(|(t, s)| t * s).sum::<f64>(),
        1e-10,
    );
    for i in 0..2 {
        close(
            coupled.state().transverse_shear_entropy_j_k()[i],
            entropy[i],
            2e-7 * entropy[i].abs() + 1e-17,
        );
        let mut forward = temperatures;
        forward[i] = current[i];
        let mut reverse = current;
        reverse[i] = temperatures[i];
        let delta_h = 0.5
            * (shear_modulus(&forward, &slopes) - shear_modulus(&temperatures, &slopes)
                + shear_modulus(&current, &slopes)
                - shear_modulus(&reverse, &slopes));
        let parameter_work = 0.5 * cross * delta_h;
        let temperature_entropy = current[i] * entropy[i] - temperatures[i] * old_entropy[i];
        close(
            step.report.cells[i].shear_stiffness_temperature_work_j,
            parameter_work,
            1e-10,
        );
        close(
            step.report.cells[i].shear_temperature_entropy_change_j,
            temperature_entropy,
            1e-10,
        );
        close(
            step.report.cells[i].shear_reversible_heat_j,
            -parameter_work - temperature_entropy,
            1e-10,
        );
    }
    // Temperature differentiation acts on the entire heterogeneous section,
    // and is observably different from independently assigning 5/6 G_i h_i.
    let local_guess = -0.5 * g2 * (5.0 / 6.0) * 0.5 * slopes[0] * H / 2.0;
    assert!((entropy[0] - local_guess).abs() > 1e-4 * entropy[0].abs());
    close(
        step.thermal.energy_j() + field_internal(coupled.state())
            - thermal.energy_j()
            - initial_field,
        step.report.external_heat_j
            + step.report.external_mechanical_work_j
            + step.report.drilling_gauge_parameter_work_j,
        1e-7,
    );
    close(step.report.combined_energy_residual_j, 0.0, 1e-7);
    assert_eq!(thermal.cells(), &original_cells);

    // Cancel after the same proposed heat, changed native factorization and
    // final accounting, immediately before publication of either owner.
    let mut cancelled = fresh();
    let original_u = cancelled.state().shell().response().displacement_6n.clone();
    let mut replay = 0;
    assert!(matches!(
        cancelled.advance_to(
            &thermal,
            0.5,
            &load,
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || {
                replay += 1;
                replay == checkpoints
            }
        ),
        Err(LaminateThermoelasticError::Cancelled)
    ));
    assert_eq!(cancelled.time_s(), 0.0);
    assert_eq!(
        cancelled.state().shell().response().displacement_6n,
        original_u
    );
    assert_eq!(thermal.cells(), &original_cells);
}

#[test]
fn g0_source_capacity_and_static_dynamic_domains_refuse_without_committing_heat() {
    for excessive_expansion in [false, true] {
        let cards = [card_with_coefficients(
            "softening admission",
            SLOPE,
            if excessive_expansion { 1e-3 } else { ALPHA },
            0.0,
            CP,
        )];
        let body = prepare_body(&cards, &[1], 1.0);
        let heat_laws = body.thermal_laws(false).unwrap();
        let bound = body.bind(&heat_laws, ids()).unwrap();
        let thermal = initial_with_sources(&bound, &[T0], &[6000.0]);
        let original_cells = thermal.cells().clone();
        let shell = LayeredThermalShell::try_new(
            &bound,
            mechanical(&cards),
            ReferenceSurface::Midsurface,
            discretization(),
            limits(),
            || false,
        )
        .unwrap();
        let range = if excessive_expansion {
            [280.0, 360.0]
        } else {
            [299.0, 301.0]
        };
        let result = ReciprocalLayeredThermalShell::try_new(
            &shell,
            &thermal,
            &[0.0; 24],
            &[range],
            0.02,
            controls(),
            || false,
        );
        if excessive_expansion {
            assert!(
                matches!(result, Err(LaminateThermoelasticError::Invalid(message)) if message.contains("capacity"))
            );
        } else {
            let mut coupled = result.unwrap();
            let old_u = coupled.state().shell().response().displacement_6n.clone();
            let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
            assert!(
                coupled
                    .advance_to(
                        &thermal,
                        1.0,
                        &[0.0; 24],
                        conduction(),
                        LayeredThermalSolveControls::default(),
                        controls(),
                        &mut workspace,
                        || false,
                    )
                    .is_err()
            );
            assert_eq!(coupled.time_s(), 0.0);
            assert_eq!(coupled.state().shell().response().displacement_6n, old_u);
            let mut dynamic = coupled
                .begin_transient(&[0.0; 24], &[[0.0; 2]; 2], 1e-7, || false)
                .unwrap();
            assert!(
                dynamic
                    .advance_to(
                        &thermal,
                        1.0,
                        &[0.0; 24],
                        conduction(),
                        LayeredThermalSolveControls::default(),
                        controls(),
                        &mut workspace,
                        || false,
                    )
                    .is_err()
            );
            assert_eq!(dynamic.time_s(), 0.0);
            assert_eq!(dynamic.state().response().displacement_6n, old_u);
            assert_eq!(dynamic.state().temperatures_k(), &[T0]);
        }
        assert_eq!(thermal.cells(), &original_cells);
        assert_eq!(thermal.time_s(), 0.0);
    }
}
