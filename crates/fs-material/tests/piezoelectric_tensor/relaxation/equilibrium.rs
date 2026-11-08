//! Independent static Debye polarization and source-domain checks.
use super::*;
use fs_material::piezoelectric::tensor::relaxation::transient::MemoryResponse;

#[test]
fn g1_g3_equilibrium_poles_add_static_energy_without_maxwell_stress_and_corotate() {
    let direction = [0.6, 0.8, 0.0];
    let (card, mechanical, dielectric) = fixture(diagonal_branch(), direction, |_| {});
    let field = [1200.0, -700.0, 300.0];
    let along = 160.0;
    let q = rotation([1.0 / 3.0_f64.sqrt(); 3], 0.73);
    for count in [0, 1, 2] {
        let law = RelaxingTensorPiezoelectric::resolve(
            base_law(&card),
            &card,
            vec![mechanical.clone()],
            vec![dielectric.clone(); count],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let response = law.equilibrium_response(field).unwrap();
        let polarization = count as f64 * PERMITTIVITY_INCREMENT_F_M * along;
        let energy = 0.5 * polarization * along;
        assert_eq!(response.additional_stress_mandel_pa, [0.0; 6]);
        assert_eq!(response.mechanical_memory_energy_j_m3, 0.0);
        for i in 0..3 {
            near(
                response.polarization_c_m2[i],
                direction[i] * polarization,
                1e-20,
                3e-13,
            );
        }
        near(response.dielectric_memory_energy_j_m3, energy, 1e-18, 3e-13);
        assert_eq!(
            response.memory_energy_j_m3,
            response.dielectric_memory_energy_j_m3
        );
        let reversed = law.equilibrium_response(field.map(|e| -e)).unwrap();
        assert_eq!(reversed.memory_energy_j_m3, response.memory_energy_j_m3);
        for i in 0..3 {
            near(
                reversed.polarization_c_m2[i],
                -response.polarization_c_m2[i],
                1e-20,
                3e-13,
            );
        }
        assert_eq!(
            law.equilibrium_response([0.0; 3]).unwrap(),
            MemoryResponse::default()
        );
        assert_eq!(
            law.equilibrium_response([0.0, 0.0, 420.0]).unwrap(),
            MemoryResponse::default()
        );
        let moved = law.rotated(ContentHash([83; 32]), q).unwrap();
        let mapped = moved.equilibrium_response(apply(&q, &field)).unwrap();
        for (actual, expected) in mapped
            .polarization_c_m2
            .iter()
            .zip(apply(&q, &response.polarization_c_m2))
        {
            near(*actual, expected, 1e-20, 3e-12);
        }
        near(
            mapped.memory_energy_j_m3,
            response.memory_energy_j_m3,
            1e-18,
            3e-12,
        );
        // omega=0 in the harmonic evaluator denotes a sinusoidal limit: its
        // mean pole energy is HALF the actual static energy at the same field.
        let harmonic = law
            .harmonic_response(0.0, [C64::from_re(2e-5); 6], field.map(C64::from_re))
            .unwrap();
        near(
            response.memory_energy_j_m3,
            2.0 * harmonic.average_dielectric_memory_energy_j_m3,
            1e-18,
            3e-13,
        );
        assert_eq!(harmonic.average_mechanical_memory_energy_j_m3, 0.0);
    }
}

#[test]
fn g0_equilibrium_requires_dc_source_support_but_not_an_unrestricted_transient_band() {
    for lower in [0.0, 10.0] {
        let (card, mechanical, dielectric) = fixture(diagonal_branch(), [0.0, 0.0, 1.0], |claim| {
            if claim.key.name() == "piezo_polarization_tau" {
                claim.validity = claim.validity.clone().with("omega", lower, 1000.0);
            }
        });
        let law = bound(&card, &mechanical, &dielectric);
        assert_eq!(law.frequency_band_rad_s(), [lower, 1000.0]);
        assert!(law.frequency_coefficients(100.0).is_ok());
        assert!(law.zero_memory().is_err());
        assert_eq!(
            law.equilibrium_response([0.0, 0.0, 2000.0]).is_ok(),
            lower == 0.0
        );
    }
}

#[test]
fn g1_equilibrium_scales_before_squaring_and_refuses_unrepresentable_fields() {
    let (card, mechanical, dielectric) = fixture(diagonal_branch(), [0.0, 0.0, 1.0], |_| {});
    let law = bound(&card, &mechanical, &dielectric);
    let field = 1e155;
    let response = law.equilibrium_response([0.0, 0.0, field]).unwrap();
    let polarization = PERMITTIVITY_INCREMENT_F_M * field;
    near(response.polarization_c_m2[2], polarization, 0.0, 3e-13);
    near(
        response.memory_energy_j_m3,
        (0.5 * polarization) * field,
        0.0,
        3e-13,
    );
    for field in [
        [f64::NAN, 0.0, 0.0],
        [0.0, f64::INFINITY, 0.0],
        [0.0, 0.0, 1e308],
    ] {
        assert!(law.equilibrium_response(field).is_err());
    }
}
