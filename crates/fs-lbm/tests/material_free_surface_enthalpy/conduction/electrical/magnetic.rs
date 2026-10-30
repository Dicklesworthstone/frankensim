//! Independent Ohm/motional-voltage/force references on actual VOF populations.

use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::electrical::FreeSurfaceMagneticInput;

mod isolated;

fn field(bz: f64, maximum_change: f64) -> FreeSurfaceMagneticInput {
    FreeSurfaceMagneticInput {
        magnetic_flux_density_t: vec![[0.0, 0.0, bz]; 3],
        maximum_velocity_change_m_s: maximum_change,
    }
}

fn velocity(liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>) -> [f64; 3] {
    let native = liquid.flow().lattice().sample_at(CORE[0], CORE[1], CORE[2]).unwrap();
    // dx=dt=0.001 in this independently scaled fixture.
    native.velocity
}

#[test]
fn g1_motional_emf_brakes_native_liquid_and_heats_it_without_electrode_power() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = electrical_card(basis, 1.0, false, false, 500.0);
    with_initial_flow(&card, basis, false, Some([0.0, 0.02, 0.0]), [0.0; 3], |liquid, gate| {
        let mut liquid = liquid.with_electrical_core(current(0.0), gate).unwrap()
            .with_imposed_magnetic_field(field(10.0, 0.03)).unwrap();
        let initial_velocity = velocity(&liquid);
        let report = liquid.step_with_conduction(&[], &[], controls(), gate).unwrap();
        let electric = report.electrical.as_ref().unwrap();
        let magnetic = electric.last_applied.magnetic.as_ref().unwrap();
        // Bulk length=.001 m, area=1e-6 m², rho_e=.001 ohm m: R=1 ohm.
        // u x B adds +0.0002 V along +x. The high-x grounded electrode therefore
        // delivers -0.0002 A, and I*L*B gives a -2e-6 N mechanical reaction.
        let emf = 0.02 * 10.0 * 0.001;
        let current = -emf;
        let force = current * 0.001 * 10.0;
        let heat = current * current * 0.001;
        near(electric.last_applied.solution.electrodes[1].current_a, current, 2e-11);
        near(magnetic.lorentz_force_density_n_m3[1][1], force / 1e-9, 2e-11);
        near(magnetic.mechanical_power_w, force * 0.02, 2e-11);
        near(electric.joule_heat_j, heat, 2e-11);
        near(electric.sampled_lorentz_work_j, -heat, 2e-11);
        near(electric.lorentz_impulse_n_s[1][1], force * 0.001, 2e-11);
        assert_eq!(electric.electrode_energy_j, 0.0);
        near(report.flow.external_heat_j, heat, 2e-11);
        assert!(electric.energy_defect_j.abs() < 1e-18);
        assert!(magnetic.transfer_defect_w.abs() < 1e-14);
        assert_eq!(magnetic.velocity_m_s[1], initial_velocity);
        let accepted_velocity = velocity(&liquid);
        assert!(accepted_velocity[1] < initial_velocity[1]);
        assert!(report.flow.magnetic_motion.unwrap().core_kinetic_energy_change_j.is_finite());
        assert_eq!(liquid.flow().lattice().force_density_at(CORE[0], CORE[1], CORE[2]).unwrap(), [0.0; 3]);
        // No stored half-force is sampled on the next current solve. Its actual
        // initial physical velocity is the previous accepted population state.
        let next = liquid.step_with_conduction(&[], &[], controls(), gate).unwrap();
        let next_magnetic = next.electrical.as_ref().unwrap().last_applied.magnetic.as_ref().unwrap();
        assert_eq!(next_magnetic.velocity_m_s[1], accepted_velocity);
        assert!(next.electrical.as_ref().unwrap().sampled_lorentz_work_j < 0.0);
    });
}

#[test]
fn g1_driven_lorentz_force_reverses_and_retains_original_gravity_without_extra_heat() {
    let basis = ElectricalPropertyBasis::Conductivity;
    let card = electrical_card(basis, 1.0, false, false, 500.0);
    let run = |bz, enabled| with_initial_flow(&card, basis, false,
        Some([0.0; 3]), [0.2, 0.0, 0.0], |liquid, gate| {
            let mut liquid = liquid.with_electrical_core(current(0.002), gate).unwrap();
            if enabled { liquid = liquid.with_imposed_magnetic_field(field(bz, 0.05)).unwrap(); }
            let report = liquid.step_with_conduction(&[], &[], controls(), gate).unwrap();
            let electric = report.electrical.as_ref().unwrap();
            // Initial physical velocity is zero despite the existing gravity
            // source frame, so motional EMF and sampled mechanical work are zero.
            near(electric.joule_heat_j, 4e-9, 2e-11);
            near(electric.electrode_energy_j, 4e-9, 2e-11);
            assert!(electric.sampled_lorentz_work_j.abs() < 1e-18);
            if enabled {
                near(electric.lorentz_impulse_n_s[1][1], 0.002 * 0.001 * bz * 0.001, 2e-11);
                let magnetic = electric.last_applied.magnetic.as_ref().unwrap();
                assert!(magnetic.velocity_m_s[1].iter().all(|u| u.abs() < 1e-15));
            }
            let force = liquid.flow().lattice().force_density_at(CORE[0], CORE[1], CORE[2]).unwrap();
            near(force[0], 0.0002, 2e-13);
            assert_eq!(force[1], 0.0);
            assert_eq!(force[2], 0.0);
            (velocity(&liquid), liquid.flow().lattice().grid().state_bits(), electric.joule_heat_j)
        });
    let unbound = run(0.0, false);
    let zero = run(0.0, true);
    assert_eq!(unbound, zero, "zero imposed B must preserve the original DC/native path");
    let positive = run(10.0, true);
    let negative = run(-10.0, true);
    assert!(positive.0[1] > zero.0[1]);
    assert!(negative.0[1] < zero.0[1]);
    assert_ne!(positive.1, negative.1);
    near(positive.2, negative.2, 2e-11);
}

#[test]
fn g3_conductivity_data_changes_native_magnetic_impulse_and_joule_heat_together() {
    let run = |scale| {
        let basis = ElectricalPropertyBasis::Conductivity;
        let card = electrical_card(basis, scale, false, false, 500.0);
        with_unbound(&card, basis, false, |liquid, gate| {
            let mut liquid = liquid.with_electrical_core(current(0.002), gate).unwrap()
                .with_imposed_magnetic_field(field(10.0, 0.1)).unwrap();
            let report = liquid.step_with_conduction(&[], &[], controls(), gate).unwrap();
            let electric = report.electrical.as_ref().unwrap();
            near(electric.joule_heat_j, 4e-9 * scale, 2e-11);
            near(electric.lorentz_impulse_n_s[1][1], 2e-8 * scale, 2e-11);
            assert_eq!(electric.last_applied.materials[1].source_basis(), basis);
            assert_eq!(electric.last_applied.materials[1].material().properties().len(), 6);
            (velocity(&liquid), liquid.flow().lattice().grid().state_bits())
        })
    };
    let original = run(1.0);
    let conducting = run(2.0);
    assert!(conducting.0[1] > original.0[1]);
    assert_ne!(conducting.1, original.1);
}

#[test]
fn g4_late_magnetic_velocity_refusal_rolls_back_populations_wall_and_gas() {
    let basis = ElectricalPropertyBasis::Conductivity;
    let card = electrical_card(basis, 1.0, false, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    with_unbound(&card, basis, false, |liquid, gate| {
        let mut liquid = liquid.with_time_origin(7.0).unwrap()
            .with_electrical_core(current(0.002), gate).unwrap()
            .with_imposed_magnetic_field(field(10.0, 1e-12)).unwrap();
        let mut atmosphere = gas(100.0);
        let mut bodies = [body::wall(&solid, &chart, &context, 450.0)];
        let bits = liquid.flow().lattice().grid().state_bits();
        let energy = liquid.total_enthalpy_j().unwrap();
        let material = liquid.material_at(CORE).unwrap().clone();
        let gas_before = atmosphere.state().clone();
        let wall_before = bodies[0].specific_enthalpy_j_kg();
        let wall_source = bodies[0].material().clone();
        for _ in 0..2 {
            let error = liquid.step_with_finite_gas_environment(
                &[], &[], controls(), &mut bodies, &[body::contact()], &[],
                &mut atmosphere, &[gas_face()], gas_controls(), false, gate,
            ).unwrap_err();
            assert!(matches!(error, FreeSurfaceEnthalpyError::Invalid(
                "magnetic frozen-velocity split exceeded its velocity-change bound")));
            assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), energy);
            assert_eq!(liquid.material_at(CORE).unwrap(), &material);
            assert_eq!(atmosphere.state(), &gas_before);
            assert_eq!(bodies[0].specific_enthalpy_j_kg(), wall_before);
            assert_eq!(bodies[0].material(), &wall_source);
            assert_eq!(liquid.accepted_steps(), 0);
            assert_eq!(liquid.physical_time_s(), 7.0);
        }
        // Loosening the declared split gate makes this exact initial state
        // retryable; no current or impulse history survived the failed candidate.
        let mut liquid = liquid.with_imposed_magnetic_field(field(10.0, 0.1)).unwrap();
        let accepted = liquid.step_with_finite_gas_environment(
            &[], &[], controls(), &mut bodies, &[body::contact()], &[],
            &mut atmosphere, &[gas_face()], gas_controls(), false, gate,
        ).unwrap().conduction;
        assert!(accepted.flow.magnetic_motion.unwrap().maximum_velocity_change_m_s > 1e-12);
        assert!(accepted.electrical.unwrap().joule_heat_j > 0.0);
        assert_eq!(liquid.accepted_steps(), 1);
        assert_eq!(liquid.physical_time_s(), 7.0 + liquid.flow().scaling().timestep_s());
        assert!(liquid.clone().with_time_origin(0.0).is_err());
    });
}
