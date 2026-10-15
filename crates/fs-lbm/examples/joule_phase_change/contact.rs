//! Cold synthetic interface data through the canonical source/card evaluator.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,InterfaceSystemCard,InterpolationPolicy,MaterialStateId,
    PropertyClaim,PropertyKey,PropertyValue,Provenance,QueryPoint,SurfaceSpec,SystemContext,UncertaintyModel};
use fs_material::electrical_interface::{ELECTRICAL_CONTACT_RESISTANCE_PROPERTY as KEY,
    ELECTRICAL_CONTACT_RESISTANCE_DIMS as DIMS,resolve_electrical_interface_state};
use fs_material::state_point::MaterialPropertySelection;
use fs_lbm::conjugate::enthalpy::transport::{EnthalpyTransportSample,interfaces::ElectricalInterfaceResponse};
use fs_lbm::conjugate::material_coupling::interfaces::MaterialContact;
use fs_qty::{Dims,Temperature};

pub(super) struct ContactTable {
    pub heat_fraction:f64,
    axis:&'static str,
    card:InterfaceSystemCard,
    surfaces:[SurfaceSpec;2],
    material_identity:fs_blake3::ContentHash,
}
impl ContactTable {
    // The receiving assignment is owned by the example's physical setup, not
    // copied from whichever interface card happens to be queried later.
    pub(super) fn parse(text:&str,table:&Table)->Result<Self,Box<dyn Error>> {
        if text.len()>65_536 {return Err("contact table exceeds 64 KiB".into());}
        let rows:Vec<_>=text.lines().map(|s|s.split('#').next().unwrap_or("").trim())
            .filter(|s|!s.is_empty()).collect();
        if !(3..=257).contains(&rows.len()) {return Err("contact needs an axis/share header and 2..256 rows".into());}
        let header:Vec<_>=rows[0].split_whitespace().collect();
        if header.len()!=2 {return Err("expected temperature_b_K|liquid_fraction_b heat_fraction_to_first".into());}
        let (axis,axis_dims)=match header[0] {
            "temperature_b_K"=>("T_b",Temperature::DIMS),
            "liquid_fraction_b"=>("fraction_b",Dims::NONE),
            _=>return Err("contact axis must explicitly select second-surface temperature or phase fraction".into()),
        };
        let heat_fraction:f64=header[1].parse()?;
        if !heat_fraction.is_finite() || !(0.0..=1.0).contains(&heat_fraction) {return Err("contact heat fraction must be in [0,1]".into());}
        let mut knots=Vec::new();
        for row in &rows[1..] {
            let values:Vec<f64>=row.split_whitespace().map(str::parse).collect::<Result<_,_>>()?;
            if values.len()!=2 || !values.iter().all(|x|x.is_finite()) || values[1]<0.0
                || (axis=="T_b" && values[0]<=0.0)
                || (axis=="fraction_b" && !(0.0..=1.0).contains(&values[0]))
                || knots.last().is_some_and(|&(x,_)|values[0]<=x) {
                return Err("ordered finite contact coordinates and nonnegative area-specific resistance required".into());
            }
            knots.push((values[0],values[1]));
        }
        let material_identity=table.curve.material_card_identity();
        let surface=SurfaceSpec {material:MaterialStateId {
            chemistry:"synthetic table material".into(),phase:"equilibrium solid-liquid fixture".into(),
            process:format!("declared input {:?}",material_identity),revision:0},texture_frame:"declared planar cell face".into()};
        let surfaces=[surface.clone(),surface];
        let temperature_bounds=(table.curve.knots()[0].temperature_k,table.curve.knots().last().unwrap().temperature_k);
        let mut validity=ValidityDomain::unconstrained().with("T_a",temperature_bounds.0,temperature_bounds.1)
            .with("T_b",temperature_bounds.0,temperature_bounds.1);
        // For a temperature-axis source, use the intersection, not a wider
        // overwritten support interval. The original table query also refuses
        // extrapolation independently of this declared validity envelope.
        let (lo,hi)=if axis=="T_b" {(knots[0].0.max(temperature_bounds.0),knots.last().unwrap().0.min(temperature_bounds.1))}
            else {(knots[0].0,knots.last().unwrap().0)};
        if lo>hi {return Err("contact and bulk temperature supports do not overlap".into());}
        validity=validity.with(axis,lo,hi);
        let mut claims=ClaimSet::new();
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(KEY,DIMS),
            value:PropertyValue::Curve {abscissa:axis.into(),abscissa_dims:axis_dims,knots,dims:DIMS},validity,
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,observations:vec![],
            provenance:Provenance {source:"synthetic interface input; not measured contact data".into(),license:"MIT fixture".into(),artifact:None}})?;
        let card=InterfaceSystemCard::assemble(surfaces[0].clone(),surfaces[1].clone(),SystemContext {
            medium:"declared dry contact".into(),third_body:None,environment:"fixed-volume equation fixture".into(),history:"prescribed constitutive chart".into()},claims,vec![])?;
        Ok(Self {heat_fraction,axis,card,surfaces,material_identity})
    }
    pub(super) fn response(&self,face:&ElectricalContact,sides:[&EnthalpyTransportSample<EquilibriumPhaseState>;2])
        ->Result<ElectricalInterfaceResponse<MaterialContact>,JouleEnthalpyError> {
        let fail=|e:String|ChtError::InvalidInput {field:"example.contact_source",reason:e};
        for side in sides {
            if side.response.material().material_card_identity()!=self.material_identity {
                return Err(JouleEnthalpyError::Invalid("contact card and receiving bulk material differ"));
            }
        }
        let mut point=QueryPoint::new().with("T_a",sides[0].temperature_k).map_err(|e|fail(e.to_string()))?
            .with("T_b",sides[1].temperature_k).map_err(|e|fail(e.to_string()))?;
        if self.axis=="fraction_b" {point=point.with("fraction_b",sides[1].response.material().liquid_mass_fraction()).map_err(|e|fail(e.to_string()))?;}
        let source=resolve_electrical_interface_state(&self.card,&point,MaterialPropertySelection::SingleClaimOnly)
            .map_err(|e|fail(e.to_string()))?;
        // MaterialContact owns the material/process/finish/order checks. Retain
        // that original object rather than copying R into anonymous metadata.
        let contact=MaterialContact::try_new(source,face.cells,self.surfaces.each_ref(),self.heat_fraction)?;
        Ok(contact.into_enthalpy_response(face.cells)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn table()->Table {Table::parse(include_str!("transport_phase_dependent.txt")).unwrap()}
    fn sample(table:&Table,c:u64,h:f64)->EnthalpyTransportSample<EquilibriumPhaseState> {
        let response=table.transport(h).unwrap();let temperature_k=response.material().temperature_k();
        EnthalpyTransportSample {cell:DepositedCell {id:CellId(c),segment_id:c,body_id:1,material_id:MaterialId(1),
            mass_kg:table.density*1e-6,specific_h_j_kg:h,geometry:CellGeometry {
                start_m:[c as f64*0.01,0.005,0.005],end_m:[(c+1) as f64*0.01,0.005,0.005],volume_m3:1e-6,surface_area_m2:0.0006},
            boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.1,reservoir:None,source_w:0.0}},temperature_k,response}
    }
    fn face()->ElectricalContact {ElectricalContact {cells:[0,1],resistance_ohm_m2:0.0,heat_fraction_to_first:0.25}}
    #[test]
    fn original_source_receipts_follow_distinct_phase_states_and_surface_temperatures() {
        let t=table();let source=ContactTable::parse(include_str!("contact_low.txt"),&t).unwrap();
        let first=sample(&t,0,0.0);let a=sample(&t,1,10000.0);let b=sample(&t,1,20000.0);
        let ra=source.response(&face(),[&first,&a]).unwrap();let rb=source.response(&face(),[&first,&b]).unwrap();
        assert_eq!(a.temperature_k,b.temperature_k);assert!(rb.resistance_ohm_m2()>ra.resistance_ohm_m2());
        let receipt=rb.material().source().resolved();
        assert_eq!(receipt.card_identity(),source.card.content_hash());
        assert_eq!(receipt.state_coordinate("T_a"),Some(300.0));assert_eq!(receipt.state_coordinate("T_b"),Some(310.0));
        assert_eq!(receipt.state_coordinate("fraction_b"),Some(0.5));
        assert_eq!(rb.resistance_ohm_m2(),rb.material().face().resistance_ohm_m2);
        assert!(rb.material().clone().into_enthalpy_response([1,0]).is_err());
        assert!(rb.material().clone().into_enthalpy_response([1,2]).is_err());
    }
    #[test]
    fn explicit_temperature_axis_and_unsupported_or_stale_sources_are_not_silently_replaced() {
        let t=table();let source=ContactTable::parse("temperature_b_K 0.25\n300 0.001\n410 0.023",&t).unwrap();
        let first=sample(&t,0,0.0);let a=sample(&t,1,30000.0);let b=sample(&t,1,50000.0);
        let ra=source.response(&face(),[&first,&a]).unwrap();let rb=source.response(&face(),[&first,&b]).unwrap();
        assert_eq!(a.response.material().liquid_mass_fraction(),b.response.material().liquid_mass_fraction());
        assert!(rb.resistance_ohm_m2()>ra.resistance_ohm_m2());
        let short=ContactTable::parse("temperature_b_K 0.25\n300 0.001\n309 0.003",&t).unwrap();
        assert!(short.response(&face(),[&first,&a]).is_err());
        let other=Table::parse(&include_str!("transport_phase_dependent.txt").replace("1000 0.01","2000 0.01")).unwrap();
        assert!(source.response(&face(),[&first,&sample(&other,1,30000.0)]).is_err());
        for bad in ["liquid_fraction_b 1.1\n0 0.001\n1 0.003", "liquid_fraction_b 0.25\n0 -0.1\n1 0.003",
            "liquid_fraction_b 0.25\n1 0.001\n0 0.003", "temperature_b_K 0.25\nNaN 0.001\n410 0.003",
            "material_name 0.25\n0 0.001\n1 0.003"] {assert!(ContactTable::parse(bad,&t).is_err());}
    }
    #[test]
    fn changing_only_contact_data_changes_melting_and_current_without_losing_pulse_history() {
        let t=table();let low=ContactTable::parse(include_str!("contact_low.txt"),&t).unwrap();
        let high=ContactTable::parse(include_str!("contact_high.txt"),&t).unwrap();
        let mut a=Vec::new();let mut b=Vec::new();
        let x=simulate_with_contact(&t,Some(&low),10.0,60.0,0.05,|r|a.push(r)).unwrap();
        let y=simulate_with_contact(&t,Some(&high),10.0,60.0,0.05,|r|b.push(r)).unwrap();
        assert_eq!(x.temperature,310.0);assert_eq!(y.temperature,310.0);
        assert!(x.liquid_fraction>y.liquid_fraction+0.2);assert!(a[29].current>b[29].current);
        // Independent three-cell current/enthalpy equations at this exact dt.
        assert!((x.liquid_fraction-0.61046337).abs()<3e-5);assert!((y.liquid_fraction-0.32950283).abs()<3e-5);
        for rows in [&a,&b] {
            assert!(rows[20].contact_resistance.unwrap()>rows[10].contact_resistance.unwrap());
            assert!(rows[30..].iter().all(|r|r.current==0.0 && r.electrode==rows[29].electrode));
            assert!(rows.iter().all(|r|r.residual.abs()<1e-6));
            assert!(rows[59].enthalpy>0.0);
        }
        let finer=simulate_with_contact(&t,Some(&low),10.0,60.0,0.025,|_|{}).unwrap();
        assert!((finer.liquid_fraction-0.61029042).abs()<3e-5);
        assert!((finer.liquid_fraction-x.liquid_fraction).abs()<0.001);
    }
    #[test]
    fn sourced_zero_contact_reproduces_the_unmodified_material_transport_route() {
        let t=table();let zero=ContactTable::parse("liquid_fraction_b 0.25\n0 0\n1 0",&t).unwrap();
        let a=simulate(&t,10.0,60.0,0.05,|_|{}).unwrap();
        let b=simulate_with_contact(&t,Some(&zero),10.0,60.0,0.05,|_|{}).unwrap();
        assert!((a.enthalpy-b.enthalpy).abs()<1e-7);assert!((a.electrode-b.electrode).abs()<1e-9);
        assert_eq!(b.contact_resistance,Some(0.0));
    }
}
