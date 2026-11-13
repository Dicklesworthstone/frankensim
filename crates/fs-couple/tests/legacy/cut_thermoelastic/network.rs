//! Same material conduction owner changes the actual reciprocal trajectory.
use super::*;
use fs_couple::circuit_enthalpy::contacts::heat::MaterialHeatNetwork;
use fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork;
use fs_couple::circuit_enthalpy::material::{ConductionLink, MaterialCell, MaterialNetwork};

#[test]
fn prepared_source_conduction_changes_the_same_heated_equilibrium_and_preserves_closed_work() {
    with_cx(|cx,_| {
        let card=card();let law=law(&card);let curve=caloric(&card);let context=point(300.0);
        let cold=source(&law,300.0);
        let mut model=build(false,&cold);let mut bare_model=build(false,&cold);
        let mut thermal=timeline(&model,curve.curve(),true);
        let mut bare=thermal.clone();let loads=forces(&model);
        let mut solid=CutThermoelasticSolid3::try_new(model.source(),&thermal,
            cells(&model,&card,&context,curve.curve(),&law),&loads,
            settings(model.operator().n()),config(),cx).unwrap();
        let mut bare_solid=CutThermoelasticSolid3::try_new(bare_model.source(),&bare,
            cells(&bare_model,&card,&context,curve.curve(),&law),&loads,
            settings(bare_model.operator().n()),config(),cx).unwrap();
        let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
        let material_cells:Vec<_>=(0..2).map(|i|MaterialCell {
            cell:CellId(i+10),material_id:MaterialId(0),curve:curve.curve(),
            point:binding(&card,&context),
            conductivity:Some((&key,MaterialPropertySelection::SingleClaimOnly)),
        }).collect();
        let links=[ConductionLink {link:7,distance_a_m:0.25,distance_b_m:0.25}];
        let bulk=MaterialNetwork {cells:&material_cells,resistors:&[],links:&links};
        let electrical=ElectricalContactNetwork {bulk:&bulk,contacts:&[],electrical_bases:&[]};
        let network=MaterialHeatNetwork {electrical:&electrical,thermal_contacts:&[],radiation:None};
        let prepared=network.prepare(cx,&thermal).unwrap();
        let sample=prepared.resolve(cx,&thermal).unwrap();
        near(sample.bulk.conductances_w_k[&7],2.0*0.73/0.5,1e-13);
        let original_energy=thermal.energy_j();let original_mass=thermal.mass_kg();
        let authored:f64=thermal.cells().values().map(|c|c.boundary.source_w).sum();
        let mut observer=Observe::default();let mut bare_observer=Observe::default();
        bare_solid.advance(bare_model.mutable(),&mut bare,&mut bare_observer,1.0,&loads,
            step(),config(),cx).unwrap();
        let report=solid.advance_heat_network(model.mutable(),&mut thermal,&prepared,&mut observer,
            1.0,&loads,step(),config(),cx).unwrap();
        assert!(thermal.temperature_k(CellId(10)).unwrap()>bare.temperature_k(CellId(10)).unwrap()+1e-4);
        assert!(thermal.temperature_k(CellId(11)).unwrap()<bare.temperature_k(CellId(11)).unwrap()-1e-4);
        assert!(solid.equilibrium().solution.coefficients().iter()
            .zip(bare_solid.equilibrium().solution.coefficients()).any(|(a,b)|(a-b).abs()>1e-7));
        near(report.external_heat_j,authored,1e-9);
        near(thermal.energy_j()-original_energy+report.internal_energy_change_j
            -report.external_work_j+report.numerical_work_j,authored,2e-9);
        assert_eq!(report.closed_radiation_net_j,0.0);
        assert_eq!(thermal.mass_kg(),original_mass);
        assert_eq!(observer.steps,report.thermal_audit.substeps);
        assert_eq!(thermal.links()[&7].conductance_w_k,0.0,
            "source conductance is a native step override, not rewritten authored history");
        assert_eq!(solid.time_s(),thermal.time_s());
        // A failed late observer through the same prepared network must leave
        // thermal, native material, compatible state and independent observer paired.
        let before=thermal.cells().clone();let before_action=action(&model);
        let before_u=solid.equilibrium().solution.coefficients().to_vec();
        let mut failing=Observe {fail_after_s:Some(0.95),..Default::default()};
        let untouched=failing.clone();
        assert!(solid.advance_heat_network(model.mutable(),&mut thermal,&prepared,&mut failing,
            2.0,&loads,step(),config(),cx).is_err());
        assert_eq!(thermal.cells(),&before);assert_eq!(solid.time_s(),1.0);
        assert_eq!(action(&model),before_action);
        assert_eq!(solid.equilibrium().solution.coefficients(),before_u);
        assert_eq!(failing,untouched);
        native!(&model,validate_material_mass,&solid.equilibrium().mass).unwrap();
    });
}
