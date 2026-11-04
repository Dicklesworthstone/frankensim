use super::*;
use crate::{SpeciesId, solution::BinaryRegularSolutionV1};
use fs_qty::{AmountConcentration, Temperature};

fn model(chi: f64, c: f64, kappa: f64) -> BinarySquareGradientV1 {
    BinarySquareGradientV1::new(
        BinaryRegularSolutionV1::new(
            [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
            chi,
            0.0,
            Temperature::new(250.0),
            Temperature::new(800.0),
        )
        .unwrap(),
        Temperature::new(300.0),
        AmountConcentration::new(c),
        kappa,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "{actual:e} != {expected:e} within {relative:e}"
    );
}
fn interface(model: &BinarySquareGradientV1) -> BinaryPlanarInterfaceV1 {
    model
        .equilibrium_interface(BinaryPlanarInterfaceOptionsV1::default(), || false)
        .unwrap()
        .unwrap()
}

#[test]
fn g1_planar_interface_matches_independent_regular_solution_quadrature() {
    let owner = model(3.0, 1000.0, 1e-6);
    let result = interface(&owner);
    // Independent 75-digit decimal log/exp, binodal bisection and 96-point
    // Gauss-Legendre integration of sqrt(2*kappa*c*RT*(f(x)-f(x_alpha))).
    close(result.surface_tension_n_m(), 0.455_281_477_458_744_45, 2e-8);
    close(
        result.central_slope_width_m(),
        1.132_563_344_694_688_8e-6,
        2e-12,
    );
    close(
        result.center_excess_gibbs_j_per_mol(),
        287.333_283_298_325_74,
        2e-12,
    );
    assert_eq!(result.model(), &owner);
    assert_eq!(result.model().temperature(), Temperature::new(300.0));
    let coexistence = owner
        .thermodynamics()
        .coexistence(owner.temperature(), result.options().coexistence, || false)
        .unwrap()
        .unwrap();
    assert_eq!(result.coexistence(), &coexistence);
    assert!(result.estimated_absolute_error_n_m() > 0.0);
    assert!(
        result.estimated_absolute_error_n_m()
            <= result.options().relative_tolerance * result.surface_tension_n_m()
    );
    assert!(
        (5..=result.options().maximum_integrand_evaluations)
            .contains(&result.integrand_evaluations())
    );
}

#[test]
fn g1_interface_retains_near_critical_three_halves_and_inverse_half_limits() {
    // f(1/2+s)=constant-epsilon*s^2+(4/3)*s^4+O(s^6),
    // s_binodal^2=3*epsilon/8. With kappa|grad x|^2/2 this gives
    // sigma=(1/2)*sqrt(kappa*c*RT)*epsilon^(3/2),
    // gap/max_slope=2*sqrt(kappa/(c*RT*epsilon)).
    for offset in [1e-4, 1e-6] {
        let owner = model(2.0 + offset, 1000.0, 1e-6);
        let epsilon = owner.chi() - 2.0;
        let result = interface(&owner);
        let sigma = 0.5
            * (owner.gradient_coefficient_j_m()
                * owner.molar_concentration_mol_m3()
                * owner.rt_j_per_mol())
            .sqrt()
            * epsilon.powf(1.5);
        let width = 2.0
            * (owner.gradient_coefficient_j_m()
                / (owner.molar_concentration_mol_m3() * owner.rt_j_per_mol() * epsilon))
                .sqrt();
        close(result.surface_tension_n_m(), sigma, 3.0 * epsilon + 1e-7);
        close(result.central_slope_width_m(), width, 3.0 * epsilon + 1e-7);
        close(
            result.center_excess_gibbs_j_per_mol(),
            owner.rt_j_per_mol() * 3.0 * epsilon * epsilon / 16.0,
            3.0 * epsilon + 1e-7,
        );
    }
}

#[test]
fn g1_interface_kappa_concentration_and_component_exchange_scaling() {
    let base = model(3.0, 1000.0, 1e-6);
    let original = interface(&base);
    let scaled = interface(&model(3.0, 4000.0, 9e-6));
    close(
        scaled.surface_tension_n_m(),
        6.0 * original.surface_tension_n_m(),
        2e-14,
    );
    close(
        scaled.central_slope_width_m(),
        1.5 * original.central_slope_width_m(),
        2e-14,
    );
    let swapped = BinarySquareGradientV1::new(
        BinaryRegularSolutionV1::new(
            [SpeciesId::new("B").unwrap(), SpeciesId::new("A").unwrap()],
            3.0,
            0.0,
            Temperature::new(250.0),
            Temperature::new(800.0),
        )
        .unwrap(),
        base.temperature(),
        AmountConcentration::new(1000.0),
        1e-6,
    )
    .unwrap();
    let reversed = interface(&swapped);
    assert_eq!(
        reversed.surface_tension_n_m(),
        original.surface_tension_n_m()
    );
    assert_eq!(
        reversed.central_slope_width_m(),
        original.central_slope_width_m()
    );
    assert_eq!(reversed.model(), &swapped);
    assert_ne!(reversed.model(), original.model());
}

#[test]
fn g0_planar_interface_admits_no_artificial_miscible_or_critical_surface() {
    for chi in [-3.0, 0.0, 2.0] {
        assert!(
            model(chi, 1000.0, 0.0)
                .equilibrium_interface(BinaryPlanarInterfaceOptionsV1::default(), || false)
                .unwrap()
                .is_none()
        );
    }
    let owner = model(3.0, 1000.0, 1e-6);
    for controls in [
        BinaryPlanarInterfaceOptionsV1 {
            relative_tolerance: f64::NAN,
            ..Default::default()
        },
        BinaryPlanarInterfaceOptionsV1 {
            absolute_tolerance_n_m: -1.0,
            ..Default::default()
        },
        BinaryPlanarInterfaceOptionsV1 {
            maximum_integrand_evaluations: 4,
            ..Default::default()
        },
        BinaryPlanarInterfaceOptionsV1 {
            coexistence: BinaryCoexistenceOptionsV1 {
                maximum_iterations: 0,
                ..Default::default()
            },
            ..Default::default()
        },
    ] {
        assert!(matches!(
            owner.equilibrium_interface(controls, || false),
            Err(BinarySolutionErrorV1::Invalid(_))
        ));
    }
    // The original coexistence owner refuses when a rich endpoint rounds pure;
    // the interface must propagate that refusal, never invent a clipped phase.
    assert!(matches!(
        model(40.0, 1000.0, 1e-6).equilibrium_interface(Default::default(), || false),
        Err(BinarySolutionErrorV1::Unrepresentable(_))
    ));
}

#[test]
fn g0_interface_work_limit_and_mid_quadrature_cancellation_are_explicit() {
    let owner = model(3.0, 1000.0, 1e-6);
    assert_eq!(
        owner.equilibrium_interface(
            BinaryPlanarInterfaceOptionsV1 {
                relative_tolerance: 1e-12,
                maximum_integrand_evaluations: 5,
                ..Default::default()
            },
            || false,
        ),
        Err(BinarySolutionErrorV1::IterationLimit)
    );
    let mut root_checks = 0;
    owner
        .thermodynamics()
        .coexistence(owner.temperature(), Default::default(), || {
            root_checks += 1;
            false
        })
        .unwrap();
    let mut checks = 0;
    assert_eq!(
        owner.equilibrium_interface(Default::default(), || {
            checks += 1;
            checks > root_checks + 3
        }),
        Err(BinarySolutionErrorV1::Cancelled)
    );
    assert_eq!(checks, root_checks + 4);
}
