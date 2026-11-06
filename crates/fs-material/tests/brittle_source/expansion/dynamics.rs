//! Same-source density and work-conjugate thermal force assembly.
//! These manufactured cases verify the admitted equations, not real material
//! data, a thermal time-convergence order or a full solid trajectory.

use super::super::{card, near, properties, scalar, settings};
use super::{curve_card, query};
use fs_evidence::ValidityDomain;
use fs_matdb::{InterpolationPolicy, MaterialCard, PropertyClaim, PropertyValue};
use fs_material::SmallStrainLaw;
use fs_material::brittle::source::BrittleSourceError;
use fs_material::brittle::source::expansion::dynamics::resolve_thermal_brittle_dynamic_material;
use fs_material::brittle::source::expansion::{ResolvedThermalBrittle, resolve_thermal_brittle};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection,
};
use fs_material::tensor::contract;
use fs_qty::{Density, Pressure, QuantitySpec, Temperature};

fn dynamic_claims(density_300: f64, density_400: f64) -> Vec<PropertyClaim> {
    let mut claims = properties();
    for (name, dims, values) in [
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            [1e-5, 3e-5],
        ),
        (DENSITY_PROPERTY, Density::DIMS, [density_300, density_400]),
    ] {
        let mut property = scalar(name, dims, values[0]);
        property.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, values[0]), (400.0, values[1])],
            dims,
        };
        property.validity = ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            300.0,
            400.0,
        );
        property.interpolation = InterpolationPolicy::LinearInside;
        claims.push(property);
    }
    claims
}

fn dynamic(
    source: &MaterialCard,
    reference: f64,
    current: f64,
) -> Result<ResolvedThermalBrittle, BrittleSourceError> {
    resolve_thermal_brittle_dynamic_material(
        source,
        &query(reference),
        &query(current),
        MaterialPropertySelection::SingleClaimOnly,
        settings(0.01),
    )
}

#[test]
fn density_expansion_and_fracture_keep_one_source_and_one_expansion_reference() {
    let source = card(dynamic_claims(1000.0, 800.0));
    let start = dynamic(&source, 300.0, 300.0).unwrap();
    let hot = start.at_temperature(&source, 400.0).unwrap();
    let current = dynamic(&source, 300.0, 400.0).unwrap();

    near(start.band().reference_density_kg_m3().unwrap(), 1000.0);
    near(hot.band().reference_density_kg_m3().unwrap(), 800.0);
    near(current.band().reference_density_kg_m3().unwrap(), 800.0);
    // Integral of instantaneous alpha, not alpha(current) times delta T.
    near(hot.free_linear_strain(), 0.002);
    assert_eq!(hot, current);
    assert_eq!(start.band().material().properties().len(), 9);
    assert_eq!(hot.band().material().card_identity(), source.content_hash());
    assert_eq!(
        start.expansion().reference().resolved().identity(),
        hot.expansion().reference().resolved().identity(),
    );
    for (old, new) in start
        .band()
        .material()
        .properties()
        .iter()
        .zip(hot.band().material().properties())
    {
        assert_eq!(old.requirement(), new.requirement());
        assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
    }
    assert_eq!(
        start
            .band()
            .material()
            .property(DENSITY_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .selected,
        hot.band()
            .material()
            .property(DENSITY_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .selected,
    );
    assert!(start.at_temperature(&source, 400.1).is_err());
    // Resolution reports a changed density value; it does not supply a mass
    // update. The motion consumer retains its originally bound inventory.
    near(start.band().reference_density_kg_m3().unwrap(), 1000.0);
}

#[test]
fn total_strain_mean_force_has_native_work_and_the_correct_displacement_derivative() {
    let source = card(dynamic_claims(1000.0, 800.0));
    let start = dynamic(&source, 300.0, 300.0).unwrap();
    let next = start.at_temperature(&source, 301.0).unwrap();
    let old_total = [0.00025, 0.0001, -0.00005, 0.00003, 0.00002, 0.00001];
    let loaded = start
        .trial_total(&old_total, &start.band().law().initial_state())
        .unwrap();
    assert_eq!(loaded.state.damage(), 0.0);
    let prepared = start
        .transfer_coefficients(&loaded.state, &next, 1e-6)
        .unwrap();
    near(prepared.parameter_work_j_m3, 0.0);
    let total = [0.00035, 0.00005, 0.00002, 0.00005, 0.00001, 0.00003];
    let mechanical = next.mechanical_strain(&total).unwrap();
    let native = next
        .band()
        .law()
        .trial_with_mean_stress_tangent(&mechanical, &prepared.state)
        .unwrap();
    let force = next.mean_small_strain_law();
    let (path, tangent) = force.trial(&total, &prepared.state).unwrap();
    assert_eq!(path, native.0);
    assert_eq!(tangent, native.1);
    assert_eq!(force.stress(&total, &prepared.state), path.mean_stress_pa);
    assert_eq!(
        force.update_state(&total, &prepared.state),
        path.endpoint.state
    );
    assert_eq!(path.endpoint.state.damage(), 0.0);
    assert!((path.mean_stress_pa[0] - path.endpoint.stress_pa[0]).abs() > 1.0);
    for component in 0..6 {
        near(
            path.mean_stress_pa[component],
            0.5 * (loaded.stress_pa[component] + path.endpoint.stress_pa[component]),
        );
    }

    let change = core::array::from_fn(|i| total[i] - old_total[i]);
    let thermal_change = next.free_linear_strain() - start.free_linear_strain();
    let mechanical_work = contract(&path.mean_stress_pa, &change);
    let eigenstrain_work = -thermal_change
        * (path.mean_stress_pa[0] + path.mean_stress_pa[1] + path.mean_stress_pa[2]);
    near(mechanical_work + eigenstrain_work, path.path_work_j_m3);
    near(path.path_work_j_m3, path.endpoint.energy.work);
    near(
        path.endpoint.energy.work,
        path.endpoint.energy.stored_change + path.endpoint.energy.dissipated,
    );
    assert!(eigenstrain_work < 0.0);
    assert_eq!(&mechanical[3..], &total[3..]);

    // The entire path is in smooth intact tension. The total-strain derivative
    // is half the elastic tensor, including tensor-shear column weights.
    let tangent = tangent.expect("smooth intact path has its classical mean derivative");
    assert_eq!(force.tangent(&total, &prepared.state), tangent);
    for column in 0..6 {
        let mut plus = total;
        let mut minus = total;
        let h = 1e-8;
        plus[column] += h;
        minus[column] -= h;
        let a = force.stress(&plus, &prepared.state);
        let b = force.stress(&minus, &prepared.state);
        for row in 0..6 {
            near(
                tangent[row][column],
                0.5 * path.endpoint.branch_tangent_pa[row][column],
            );
            let derivative = (a[row] - b[row]) / (2.0 * h);
            assert!((derivative - tangent[row][column]).abs() < 100.0);
        }
    }
    assert_eq!(loaded.state.strain(), &old_total);
}

#[test]
fn matched_thermal_motion_is_stress_free_and_restraint_retains_irreversible_failure() {
    let source = card(dynamic_claims(1000.0, 800.0));
    let start = dynamic(&source, 400.0, 400.0).unwrap();
    let cold = start.at_temperature(&source, 300.0).unwrap();
    let virgin = start.band().law().initial_state();
    let prepared = start.transfer_coefficients(&virgin, &cold, 1e-6).unwrap();
    let force = cold.mean_small_strain_law();
    let free = cold.free_linear_strain();
    near(free, -0.002);
    let free_total = [free, free, free, 0.0, 0.0, 0.0];
    let (expanded, _) = force.trial(&free_total, &prepared.state).unwrap();
    assert_eq!(expanded.mean_stress_pa, [0.0; 6]);
    assert_eq!(expanded.endpoint.stress_pa, [0.0; 6]);
    assert_eq!(expanded.endpoint.state.strain(), &[0.0; 6]);
    assert_eq!(expanded.endpoint.state.damage(), 0.0);
    near(expanded.path_work_j_m3, 0.0);

    let (restrained, _) = force.trial(&[0.0; 6], &prepared.state).unwrap();
    assert!(restrained.endpoint.state.damage() > 0.0);
    assert!(restrained.endpoint.energy.dissipated > 0.0);
    let work = -free
        * (restrained.mean_stress_pa[0]
            + restrained.mean_stress_pa[1]
            + restrained.mean_stress_pa[2]);
    near(work, restrained.path_work_j_m3);
    near(
        work,
        restrained.endpoint.state.stored_energy_j_m3()
            + restrained.endpoint.state.dissipated_j_m3(),
    );
    let saved_damage = restrained.endpoint.state.damage();
    let saved_loss = restrained.endpoint.state.dissipated_j_m3();
    let reheated = cold
        .transfer_coefficients(&restrained.endpoint.state, &start, 1e-6)
        .unwrap();
    let (unloaded, _) = start
        .mean_small_strain_law()
        .trial(&[0.0; 6], &reheated.state)
        .unwrap();
    near(unloaded.endpoint.state.damage(), saved_damage);
    near(unloaded.endpoint.state.dissipated_j_m3(), saved_loss);
    near(unloaded.endpoint.energy.dissipated, 0.0);
    near(
        unloaded.path_work_j_m3,
        -restrained.endpoint.state.stored_energy_j_m3(),
    );
    assert_eq!(unloaded.endpoint.stress_pa, [0.0; 6]);
    assert_eq!(virgin.damage(), 0.0);
    assert_eq!(virgin.dissipated_j_m3(), 0.0);
    assert_eq!(prepared.state.damage(), 0.0);
}

#[test]
fn density_admission_and_mean_force_refusals_preserve_native_boundaries() {
    let missing_density = curve_card();
    assert!(dynamic(&missing_density, 300.0, 400.0).is_err());
    // Requiring density for motion does not change the static thermal API.
    assert!(
        resolve_thermal_brittle(
            &missing_density,
            &query(300.0),
            &query(400.0),
            MaterialPropertySelection::SingleClaimOnly,
            settings(0.01),
        )
        .is_ok()
    );
    for (reference, current) in [(0.0, 800.0), (-1.0, 800.0), (1000.0, 0.0)] {
        let source = card(dynamic_claims(reference, current));
        assert!(dynamic(&source, 300.0, 400.0).is_err());
    }
    let mut wrong_units = dynamic_claims(1000.0, 800.0);
    *wrong_units.last_mut().unwrap() = scalar(DENSITY_PROPERTY, Pressure::DIMS, 1000.0);
    assert!(dynamic(&card(wrong_units), 300.0, 400.0).is_err());

    let source = card(dynamic_claims(1000.0, 800.0));
    let start = dynamic(&source, 300.0, 300.0).unwrap();
    let hot = start.at_temperature(&source, 400.0).unwrap();
    let foreign_owner = hot.clone();
    let foreign = foreign_owner.band().law().initial_state();
    let force = hot.mean_small_strain_law();
    let free = hot.free_linear_strain();
    let total = [free, free, free, 0.0, 0.0, 0.0];
    assert!(force.trial(&total, &foreign).is_err());
    assert!(force.stress(&total, &foreign).iter().all(|x| x.is_nan()));
    assert!(
        force
            .tangent(&total, &foreign)
            .iter()
            .flatten()
            .all(|x| x.is_nan())
    );
    assert!(!force.update_state(&total, &foreign).is_valid());
    let virgin = force.initial_state();
    for invalid in [
        [0.03, 0.03, 0.03, 0.0, 0.0, 0.0], // total norm alone exceeds bound
        [-0.028, -0.028, -0.028, 0.0, 0.0, 0.0], // mechanical norm alone exceeds bound
        [f64::NAN; 6],
    ] {
        assert!(force.trial(&invalid, &virgin).is_err());
        assert!(!force.update_state(&invalid, &virgin).is_valid());
    }
    let zero_force = start.mean_small_strain_law();
    let zero_state = zero_force.initial_state();
    let (held, missing_tangent) = zero_force.trial(&[0.0; 6], &zero_state).unwrap();
    assert_eq!(held.mean_stress_pa, [0.0; 6]);
    assert!(missing_tangent.is_none());
    assert!(
        zero_force
            .tangent(&[0.0; 6], &zero_state)
            .iter()
            .flatten()
            .all(|x| x.is_nan())
    );
    assert_eq!(zero_force.update_state(&[0.0; 6], &zero_state), zero_state);
    assert_eq!(virgin.damage(), 0.0);
}
