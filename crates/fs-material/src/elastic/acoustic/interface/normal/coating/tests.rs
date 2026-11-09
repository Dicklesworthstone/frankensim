use super::*;
use crate::elastic::acoustic::interface::layer::supported::IsotropicElasticCoating;
use crate::elastic::acoustic::interface::surface::ElasticSurfaceImpedance;

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

// Independently authored normal Christoffel eigensystems embedded in an SPD
// full Mandel material. Remaining xx,yy,xy moduli are independent positive data.
fn tensor_material(
    rho: f64,
    speeds: [f64; 3],
    modes: [[f64; 3]; 3],
    remaining: [f64; 3],
    identity: u8,
    ceiling: f64,
) -> ElasticAcousticState {
    let gamma: [[f64; 3]; 3] = core::array::from_fn(|i| core::array::from_fn(|j| {
        (0..3).map(|m| rho * speeds[m].powi(2) * modes[m][i] * modes[m][j]).sum()
    }));
    let mut c = [[0.0; 6]; 6];
    c[0][0] = remaining[0];
    c[1][1] = remaining[1];
    c[3][3] = remaining[2];
    let indices = [5, 4, 2];
    let scales = [core::f64::consts::SQRT_2, core::f64::consts::SQRT_2, 1.0];
    for i in 0..3 { for j in i..3 {
        let value = scales[i] * scales[j] * gamma[i][j];
        c[indices[i]][indices[j]] = value;
        c[indices[j]][indices[i]] = value;
    } }
    ElasticAcousticState::try_new_mandel(c, rho, ContentHash([identity; 32]), ceiling).unwrap()
}

fn incident_material() -> ElasticAcousticState {
    tensor_material(
        2400.0, [1500.0, 2600.0, 4300.0],
        [[4.0 / 5.0, 3.0 / 5.0, 0.0],
         [-36.0 / 65.0, 48.0 / 65.0, 5.0 / 13.0],
         [3.0 / 13.0, -4.0 / 13.0, 12.0 / 13.0]],
        [20e9, 18e9, 12e9], 81, 0.02,
    )
}

fn coating_material(ceiling: f64) -> ElasticAcousticState {
    tensor_material(
        1800.0, [1100.0, 2100.0, 3400.0],
        [[240.0 / 377.0, 252.0 / 377.0, -5.0 / 13.0],
         [21.0 / 29.0, -20.0 / 29.0, 0.0],
         [100.0 / 377.0, 105.0 / 377.0, 12.0 / 13.0]],
        [16e9, 14e9, 9e9], 82, ceiling,
    )
}

fn isotropic(rho: f64, cp: f64, cs: f64, identity: u8, ceiling: f64) -> ElasticAcousticState {
    ElasticAcousticState::from_bulk_and_shear(
        rho, rho * (cp * cp - (4.0 / 3.0) * cs * cs),
        rho * cs * cs, ContentHash([identity; 32]), ceiling,
    ).unwrap()
}

// R_y(cos=15/17,sin=8/17) diag(2e6+4e6i,3e6-1e6i,5e6+2e6i) R_y^T.
fn support() -> ElasticSurfaceCondition {
    let mut z = [[C64::ZERO; 3]; 3];
    z[0][0] = C64::new(2_664_359.861_591_695_5, 3_557_093.425_605_536_3);
    z[1][1] = C64::new(3e6, -1e6);
    z[2][2] = C64::new(4_335_640.138_408_304, 2_442_906.574_394_463_7);
    z[0][2] = C64::new(1_245_674.740_484_429, -830_449.826_989_619_4);
    z[2][0] = z[0][2];
    ElasticSurfaceCondition::TensorImpedance(ElasticSurfaceImpedanceTensor::new(z).unwrap())
}

fn incoming() -> [C64; 3] {
    [C64::from_re(0.001), C64::new(-0.0002, 0.0001), C64::new(0.0003, -0.0004)]
}

fn close(a: C64, b: C64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:?} != {b:?}");
}

fn scalar_close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:?} != {b:?}");
}

fn strain_norm(field: NormalElasticCoatingField) -> f64 {
    field.strain_mandel.iter().fold(0.0, |n, e| det::hypot(n, e.abs()))
}

// Independent direct Cartesian continuum transfer calculation, using
// [C, iY sin; iD sin, C] from rear to front followed by a 3x3 solve.
// This reference never uses the production augmented wave-boundary matrix.
#[test]
fn full_tensor_coating_matches_independent_cartesian_wave_and_interior_oracle() {
    let incident_state = incident_material();
    let coating_state = coating_material(0.02);
    let incident = NormalElasticCharacteristic::new(&incident_state, IDENTITY).unwrap();
    let coating = NormalElasticCharacteristic::new(&coating_state, IDENTITY).unwrap();
    let owner = NormalElasticCoating::new(&incident, &coating, 0.43, support()).unwrap();
    let prepared = owner.prepare(2000.0 / core::f64::consts::TAU).unwrap();
    let response = prepared.drive_cartesian(incoming()).unwrap();
    let expected_reflected = [
        C64::new(0.0003366177109050661, -0.0001014314239409185),
        C64::new(-0.0002912990675367962, 0.00004852294981711704),
        C64::new(0.00011238885506990482, -0.0002607171153104454),
    ];
    let expected_velocity = [
        C64::new(0.0009950075348079768, -0.00005844025006261862),
        C64::new(-0.0005424828164728077, -0.000317229437654918),
        C64::new(0.0005506830519272021, -0.00033126713749752676),
    ];
    let expected_traction = [
        C64::new(-3269.806797788582, -2513.663153192891),
        C64::new(1944.6778870733413, 409.2054964919463),
        C64::new(-4387.742268698551, 990.0892283745795),
    ];
    for i in 0..3 {
        close(response.reflection.reflected_velocity_m_s[i], expected_reflected[i], 3e-13);
        close(response.support_velocity_m_s[i], expected_velocity[i], 3e-13);
        close(response.support_traction_pa[i], expected_traction[i], 3e-6);
    }
    scalar_close(response.reflection.total_incoming_power_w_m2, 4.428288284023669, 4e-9);
    scalar_close(response.reflection.total_reflected_power_w_m2, 0.9104937875650404, 4e-9);
    scalar_close(response.reflection.absorbed_normal_power_w_m2, 3.51779449645863, 4e-9);
    let field = response.field_at(0.17).unwrap();
    let expected_field_velocity = [
        C64::new(0.001311081639853255, -0.00012298634521126282),
        C64::new(-0.00035986115828431963, -0.00002848906185708692),
        C64::new(0.0003871617827958179, -0.0005126670396884388),
    ];
    let expected_strain = [
        C64::ZERO, C64::ZERO,
        C64::new(-3.3674234305758e-7, 3.005893912087449e-7),
        C64::ZERO,
        C64::new(3.925004400704331e-7, -2.4658620087771044e-7),
        C64::new(-8.863510266263234e-8, -3.608306079676012e-7),
    ];
    for i in 0..3 { close(field.velocity_m_s[i], expected_field_velocity[i], 3e-13); }
    for i in 0..6 { close(field.strain_mandel[i], expected_strain[i], 3e-15); }
    scalar_close(field.normal_flux_w_m2, 3.51779449645863, 4e-9);
    for z in [0.0, 0.1, 0.17, 0.3, 0.43] {
        assert!(strain_norm(response.field_at(z).unwrap()) <= response.coating_strain_norm_bound * (1.0 + 1e-12));
    }
    assert_eq!(prepared.material_state_identities(), [
        incident_state.material_state_identity(), coating_state.material_state_identity(),
    ]);
    assert!(core::ptr::eq(prepared.incident_characteristic().state(), &incident_state));
    assert!(core::ptr::eq(prepared.coating_characteristic().state(), &coating_state));
}

#[test]
fn isotropic_limit_recovers_quarter_wave_transformer_and_original_resolved_coating() {
    let solid = isotropic(2500.0, 4000.0, 2000.0, 83, 0.02);
    let film = isotropic(1800.0, 3000.0, 1500.0, 84, 0.02);
    let incident = NormalElasticCharacteristic::new(&solid, IDENTITY).unwrap();
    let coating = NormalElasticCharacteristic::new(&film, IDENTITY).unwrap();
    let frequency = 1700.0 / core::f64::consts::TAU;
    let h = core::f64::consts::PI * 3000.0 / (2.0 * 1700.0);
    let backing = ElasticSurfaceCondition::Impedance(ElasticSurfaceImpedance::new(
        C64::from_re(2_916_000.0), C64::from_re(1e6),
    ).unwrap());
    let prepared = NormalElasticCoating::new(&incident, &coating, h, backing).unwrap()
        .prepare(frequency).unwrap();
    let response = prepared.drive([C64::ZERO, C64::ZERO, C64::from_re(0.001)]).unwrap();
    close(response.reflection.reflected_velocity_m_s[2], C64::ZERO, 2e-14);
    close(response.support_velocity_m_s[2], C64::new(0.0, 0.001851851851851852), 2e-14);
    scalar_close(response.reflection.absorbed_normal_power_w_m2, 5.0, 2e-10);
    // A nontrivial coherent tensor support compares the physical Cartesian
    // response; the legacy P reflected polarization has the opposite sign.
    let prepared = NormalElasticCoating::new(&incident, &coating, 0.37, support()).unwrap()
        .prepare(frequency).unwrap();
    let native = prepared.drive_cartesian(incoming()).unwrap();
    let old_prepared = IsotropicElasticCoating::new(&solid, &film, 0.37, support()).unwrap()
        .prepare_incidence(ElasticMode::P, 1.0, frequency).unwrap();
    let a = incoming();
    let old = old_prepared.drive([a[2], a[0], a[1]]).unwrap();
    for i in 0..3 {
        close(native.reflection.boundary_velocity_m_s[i], old.reflection.boundary_velocity_m_s[i], 3e-13);
        close(native.reflection.boundary_traction_pa[i], old.reflection.boundary_traction_pa[i], 3e-6);
        close(native.support_velocity_m_s[i], old.support_velocity_m_s[i], 3e-13);
    }
    let new_field = native.field_at(0.15).unwrap();
    let old_field = old.field_at(0.15).unwrap();
    for i in 0..3 {
        close(new_field.velocity_m_s[i], old_field.velocity_m_s[i], 3e-13);
        close(new_field.traction_pa[i], old_field.traction_pa[i], 3e-6);
    }
    for i in 0..6 { close(new_field.strain_mandel[i], old_field.strain_mandel[i], 3e-15); }
    scalar_close(native.reflection.absorbed_normal_power_w_m2, old.reflection.absorbed_normal_power_w_m2, 4e-9);
}

#[test]
fn independent_material_frame_rotations_preserve_surface_response_and_rotate_coating_strain() {
    let solid = incident_material();
    let film = coating_material(0.02);
    let qi = [[0.6, 0.0, 0.8], [0.0, 1.0, 0.0], [-0.8, 0.0, 0.6]];
    let qc = [[1.0, 0.0, 0.0], [0.0, 0.8, -0.6], [0.0, 0.6, 0.8]];
    let rotated_solid = ElasticAcousticState::try_new_mandel(
        crate::tensor::rotate_mandel(*solid.stiffness_mandel_pa(), qi),
        solid.density_kg_m3(), ContentHash([85; 32]), 0.02,
    ).unwrap();
    let rotated_film = ElasticAcousticState::try_new_mandel(
        crate::tensor::rotate_mandel(*film.stiffness_mandel_pa(), qc),
        film.density_kg_m3(), ContentHash([86; 32]), 0.02,
    ).unwrap();
    let incident = NormalElasticCharacteristic::new(&solid, IDENTITY).unwrap();
    let coating = NormalElasticCharacteristic::new(&film, IDENTITY).unwrap();
    let rotated_incident = NormalElasticCharacteristic::new(&rotated_solid, qi).unwrap();
    let rotated_coating = NormalElasticCharacteristic::new(&rotated_film, qc).unwrap();
    let frequency = 2000.0 / core::f64::consts::TAU;
    let prepared = NormalElasticCoating::new(&incident, &coating, 0.43, support()).unwrap()
        .prepare(frequency).unwrap();
    let rotated_prepared = NormalElasticCoating::new(&rotated_incident, &rotated_coating, 0.43, support()).unwrap()
        .prepare(frequency).unwrap();
    let original = prepared.drive_cartesian(incoming()).unwrap();
    let rotated = rotated_prepared.drive_cartesian(incoming()).unwrap();
    for i in 0..3 {
        close(original.support_velocity_m_s[i], rotated.support_velocity_m_s[i], 4e-13);
        close(original.support_traction_pa[i], rotated.support_traction_pa[i], 4e-6);
        close(original.reflection.boundary_velocity_m_s[i], rotated.reflection.boundary_velocity_m_s[i], 4e-13);
    }
    scalar_close(original.reflection.absorbed_normal_power_w_m2, rotated.reflection.absorbed_normal_power_w_m2, 5e-9);
    let first = original.field_at(0.17).unwrap();
    let second = rotated.field_at(0.17).unwrap();
    let transform = crate::tensor::mandel_rotation(qc);
    for i in 0..6 {
        let expected = (0..6).fold(C64::ZERO, |s, j| s + first.strain_mandel[j].scale(transform[i][j]));
        close(second.strain_mandel[i], expected, 3e-15);
    }
    scalar_close(original.coating_strain_norm_bound, rotated.coating_strain_norm_bound, 3e-15);
}

#[test]
fn absent_coating_delegates_exactly_without_inventing_an_interior_strain_domain() {
    let solid = incident_material();
    let film = coating_material(1e-12);
    let incident = NormalElasticCharacteristic::new(&solid, IDENTITY).unwrap();
    let coating = NormalElasticCharacteristic::new(&film, IDENTITY).unwrap();
    let a = incident.surface_velocity_to_modal(incoming());
    for backing in [ElasticSurfaceCondition::TractionFree, ElasticSurfaceCondition::Clamped, support()] {
        let direct = incident.prepare_surface(backing).unwrap();
        let expected = direct.drive(a).unwrap();
        let prepared = NormalElasticCoating::new(&incident, &coating, 0.0, backing).unwrap()
            .prepare(400.0).unwrap();
        let response = prepared.drive(a).unwrap();
        assert_eq!(response.reflection.reflected_velocity_m_s, expected.reflected_velocity_m_s);
        assert_eq!(response.reflection.boundary_velocity_m_s, expected.boundary_velocity_m_s);
        assert_eq!(response.reflection.boundary_traction_pa, expected.boundary_traction_pa);
        assert_eq!(response.reflection.absorbed_normal_power_w_m2, expected.absorbed_normal_power_w_m2);
        assert_eq!(response.reflection.strain_norm_bound, expected.strain_norm_bound);
        assert_eq!(response.support_velocity_m_s, expected.boundary_velocity_m_s);
        assert_eq!(response.support_traction_pa, expected.boundary_traction_pa);
        assert!(response.coating_faces.is_none());
        assert_eq!(response.coating_strain_norm_bound, 0.0);
        assert!(response.field_at(0.0).is_err());
        for m in 0..3 { assert_eq!(prepared.velocity_coefficients(m), direct.velocity_coefficients(m)); }
    }
}

#[test]
fn original_stiff_backing_keeps_tiny_rear_motion_and_positive_work() {
    let solid = isotropic(2500.0, 4000.0, 2000.0, 87, 0.02);
    let characteristic = NormalElasticCharacteristic::new(&solid, IDENTITY).unwrap();
    let resistance = 1e28;
    let backing = ElasticSurfaceCondition::Impedance(ElasticSurfaceImpedance::new(
        C64::from_re(resistance), C64::from_re(resistance),
    ).unwrap());
    let prepared = NormalElasticCoating::new(&characteristic, &characteristic, 0.4, backing).unwrap()
        .prepare(2000.0 / core::f64::consts::TAU).unwrap();
    let response = prepared.drive([C64::ZERO, C64::ZERO, C64::from_re(0.001)]).unwrap();
    let z = characteristic.mode_impedances_pa_s_m()[2];
    let amplitude = 0.001 * (2.0 * (z / resistance) / (1.0 + z / resistance));
    let (sine, cosine) = det::sin_cos(0.2);
    let expected = C64::new(cosine, sine).scale(amplitude);
    close(response.support_velocity_m_s[2], expected, expected.abs() * 3e-13);
    assert!(response.reflection.absorbed_normal_power_w_m2 > 0.0);
    let loss = 0.5 * resistance * amplitude * amplitude;
    scalar_close(response.reflection.absorbed_normal_power_w_m2, loss, loss * 5e-13);
    close(response.support_traction_pa[2], expected.scale(-resistance), 3e-8);
    assert!(loss < 1e-16 * response.reflection.total_incoming_power_w_m2);
}

#[test]
fn complete_layer_strain_envelope_original_incident_domain_and_bounded_queries_are_enforced() {
    let solid = isotropic(2500.0, 4000.0, 2000.0, 88, 0.02);
    let film = isotropic(1800.0, 3000.0, 1500.0, 89, 0.02);
    let limited_film = isotropic(1800.0, 3000.0, 1500.0, 90, 1e-8);
    let limited_solid = isotropic(2500.0, 4000.0, 2000.0, 91, 1e-8);
    let incident = NormalElasticCharacteristic::new(&solid, IDENTITY).unwrap();
    let coating = NormalElasticCharacteristic::new(&film, IDENTITY).unwrap();
    let limited_coating = NormalElasticCharacteristic::new(&limited_film, IDENTITY).unwrap();
    let limited_incident = NormalElasticCharacteristic::new(&limited_solid, IDENTITY).unwrap();
    let h = core::f64::consts::PI * 3000.0 / 1700.0;
    let frequency = 1700.0 / core::f64::consts::TAU;
    let free = ElasticSurfaceCondition::TractionFree;
    let owner = NormalElasticCoating::new(&incident, &coating, h, free).unwrap();
    let prepared = owner.prepare(frequency).unwrap();
    let input = [C64::ZERO, C64::ZERO, C64::from_re(0.001)];
    let response = prepared.drive(input).unwrap();
    let faces = response.coating_faces.unwrap();
    assert!(strain_norm(faces[0]) < 1e-8 && strain_norm(faces[1]) < 1e-8);
    assert!(strain_norm(response.field_at(0.5 * h).unwrap()) > 1e-7);
    let limited = NormalElasticCoating::new(&incident, &limited_coating, h, free).unwrap()
        .prepare(frequency).unwrap();
    assert!(limited.drive(input).is_err());
    let limited = NormalElasticCoating::new(&limited_incident, &coating, h, free).unwrap()
        .prepare(frequency).unwrap();
    assert!(limited.drive(input).is_err());
    for depth in [-0.1, h + 0.1, f64::NAN] { assert!(response.field_at(depth).is_err()); }
    for thickness in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(NormalElasticCoating::new(&incident, &coating, thickness, free).is_err());
    }
    for f in [0.0, -1.0, f64::INFINITY, f64::NAN] { assert!(owner.prepare(f).is_err()); }
    assert!(owner.prepare_cancellable(frequency, &mut || true).is_err());
    let mut polls = 0;
    assert!(owner.prepare_cancellable(frequency, &mut || { polls += 1; polls == 12 }).is_err());
    assert!(polls >= 12);
    assert!(prepared.drive([C64::from_re(f64::NAN), C64::ZERO, C64::ZERO]).is_err());
    assert!(prepared.drive([C64::ZERO; 3]).is_ok());
}
