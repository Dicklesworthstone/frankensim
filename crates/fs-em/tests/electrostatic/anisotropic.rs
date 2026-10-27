//! G1/G3: directional source data reaches native fields, electrodes and storage.
//! Synthetic absolute permittivities, not experimentally qualified specimens.
use super::*;
use fs_material::{
    conductor::planar::PlanarElectricalFrame,
    field_storage::planar::{
        PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES, PlanarDielectricLaw,
    },
};
use fs_em::electrostatic::network::CapacitanceMatrixControls;

pub(super) fn material(principal: [f64; 2]) -> ResolvedMaterialStatePoint {
    let base = source(
        PropertyValue::Scalar { value: principal[0], dims: ELECTRIC_PERMITTIVITY_DIMS },
        ValidityDomain::unconstrained(),
    );
    let mut claims = base.claims().clone();
    for (name, value) in PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES.into_iter().zip(principal) {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, ELECTRIC_PERMITTIVITY_DIMS),
            value: PropertyValue::Scalar { value, dims: ELECTRIC_PERMITTIVITY_DIMS },
            validity: ValidityDomain::unconstrained(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored principal dielectric fixture".into(),
                license: "test-only".into(), artifact: None,
            },
        }).unwrap();
    }
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    let requirements: Vec<_> = [ELECTRIC_PERMITTIVITY_PROPERTY]
        .into_iter().chain(PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES)
        .map(|name| ScalarPropertyRequirement::try_new(
            name, ELECTRIC_PERMITTIVITY_DIMS, ScalarAdmissibility::StrictlyPositive,
        ).unwrap()).collect();
    resolve_material_state_point(
        &card, &QueryPoint::new(), &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}
pub(super) fn law(angle: f64) -> PlanarDielectricLaw {
    // Distinct fixed fixture frame identities, independent of specimen data.
    let a = material([17., 19.]);
    let b = material([23., 29.]);
    PlanarDielectricLaw::Orthotropic(
        PlanarElectricalFrame::new(a.card_identity(), b.card_identity(), angle).unwrap(),
    )
}
fn mesh(angle: f64) -> TriComplex2 {
    let (s, c) = angle.sin_cos();
    TriComplex2::from_indexed_triangles(
        "oriented dielectric specimen",
        [[0.,0.],[1.,0.],[1.,1.],[0.,1.],[0.5,0.5]].into_iter()
            .map(|[x,y]| [c*x-s*y,s*x+c*y]).collect(),
        vec![[0,1,4],[1,2,4],[2,3,4],[3,0,4]],
        Metric2::planar(0.2).unwrap(),
    ).unwrap()
}
fn extraction() -> CapacitanceMatrixControls {
    CapacitanceMatrixControls {
        field: controls(), maximum_terminals: 2, maximum_field_solves: 2,
        maximum_field_values: 128, coefficient_tolerance_f: 1e-15,
    }
}

#[test]
fn oblique_tensor_field_retains_cross_displacement_and_full_stored_energy() {
    let material = material([2e-11, 6e-11]);
    let law = law(core::f64::consts::FRAC_PI_4);
    let problem = PlanarElectrostatics::new_with_laws(
        mesh(0.), vec![material.clone();4], vec![law;4], limits(), &mut || false,
    ).unwrap();
    let field = problem.solve(
        &[0.;4], &[(0,0.),(1,3.),(2,1.),(3,-2.)], controls(), &mut || false,
    ).unwrap();
    // An affine potential V=3x-2y solves the interior equation exactly.
    // Q diag(2,6) Q^T = [[4,-2],[-2,4]] in units of 1e-11 F/m.
    // Thus E=(-3,2), D=(-16,14)e-11 and U=(E.D/2)*0.2.
    near(field.voltage_v()[4], 0.5);
    for cell in 0..4 {
        near(field.electric_field_v_m()[cell][0], -3.);
        near(field.electric_field_v_m()[cell][1], 2.);
        near(field.displacement_c_m2()[cell][0], -16e-11);
        near(field.displacement_c_m2()[cell][1], 14e-11);
        near(field.cell_energy_j()[cell], 19e-12);
        assert_eq!(&problem.materials()[cell], &material);
    }
    near(field.energy_j(), 76e-12);
    assert!(problem.permittivity_f_m().is_none());
    assert_eq!(problem.permittivity_laws(), &[law;4]);

    // Equality of principal values reproduces scalar electrostatics while
    // retaining its explicit directional declaration and original source bundle.
    let equal = self::material([3e-11;2]);
    let scalar = PlanarElectrostatics::new(
        mesh(0.), vec![equal.clone();4], limits(), &mut || false,
    ).unwrap();
    let tensor = PlanarElectrostatics::new_with_laws(
        mesh(0.), vec![equal;4], vec![law;4], limits(), &mut || false,
    ).unwrap();
    let a = scalar.solve(&[0.;4], &[(0,0.),(1,3.),(2,1.),(3,-2.)], controls(), &mut || false).unwrap();
    let b = tensor.solve(&[0.;4], &[(0,0.),(1,3.),(2,1.),(3,-2.)], controls(), &mut || false).unwrap();
    near(a.energy_j(), b.energy_j());
    assert_eq!(scalar.permittivity_f_m().unwrap(), &[3e-11;4]);
    assert!(tensor.permittivity_f_m().is_none());
    for cell in 0..4 {
        for axis in 0..2 { near(a.displacement_c_m2()[cell][axis], b.displacement_c_m2()[cell][axis]); }
    }
}

#[test]
fn dielectric_orientation_changes_native_capacitance_and_circuit_decay() {
    let material = material([2e-6, 8e-6]);
    let mut remaining = Vec::new();
    for specimen_angle in [0., 0.37] {
        for transverse in [false, true] {
            let angle = specimen_angle + if transverse { core::f64::consts::FRAC_PI_2 } else { 0. };
            let problem = PlanarElectrostatics::new_with_laws(
                mesh(specimen_angle), vec![material.clone();4], vec![law(angle);4],
                limits(), &mut || false,
            ).unwrap();
            let capacitance = 0.2 * if transverse { 8e-6 } else { 2e-6 };
            let matrix = problem.capacitance_matrix(
                &[&[0,3], &[1,2]], extraction(), &mut || false,
            ).unwrap();
            for (value, expected) in matrix.maxwell_matrix_f().iter().zip([
                capacitance, -capacitance, -capacitance, capacitance,
            ]) { near(*value, expected); }
            let reduced = matrix.circuit_at_voltages(
                &[0,1], &[0.,3.], 1e-14, 1e-14, &mut || false,
            ).unwrap();
            let mut branches = reduced.branches().to_vec();
            branches.push((1,0,Branch::Resistor { ohms: 1000. }));
            let dae = assemble_circuit(&CircuitGraph {
                node_count:2, branches, transformers: vec![],
            }).unwrap();
            let mut initial = vec![0.;dae.system.state_dim()];
            for (&slot,&charge) in dae.charge_index.iter().zip(reduced.initial_charges_c()) {
                initial[slot]=charge;
            }
            let initial = dae.consistent_initial_state(&initial,&[]).unwrap();
            let dt = 1e-4;
            let (next, defect) = dae.step_audited(&initial,&[],dt).unwrap();
            // Independent two-terminal midpoint RC reference, C=epsilon_parallel*A/L.
            let ratio = (2.*1000.*capacitance-dt)/(2.*1000.*capacitance+dt);
            let voltage = 3.*ratio;
            let reference_energy = 0.5*capacitance*voltage*voltage;
            let reference_heat = dt*((3.+voltage)/2.).powi(2)/1000.;
            near(dae.system.hamiltonian(&next.x),reference_energy);
            near(next.dissipated,reference_heat);
            let physical = matrix.response(&[0.,voltage], &mut || false).unwrap();
            near(physical.energy_j,reference_energy);
            near(physical.terminal_charge_c[1],capacitance*voltage);
            near(physical.energy_j+next.dissipated,0.5*capacitance*9.);
            assert!(defect<1e-14);
            remaining.push(voltage);
        }
    }
    assert!(remaining[1]>remaining[0]);
    near(remaining[0],remaining[2]);
    near(remaining[1],remaining[3]);
}
