//! Uncertainty in the actual contact resistance, not an effective h.
//! One draw applies to the named matching or planar nonmatching contact
//! throughout a complete solve/trajectory. The child owns geometry admission.
use super::*;

#[derive(Debug, Clone)]
pub(super) struct Target {
    contact: String,
}

impl Target {
    pub(super) fn parse(value: &J) -> Result<Option<Self>> {
        if value.str_field("kind") != Some("contact-resistance") { return Ok(None); }
        object(value, &["kind", "contact"], "contact uncertainty target")?;
        Ok(Some(Self { contact: string(field(value, "contact")?, "contact")? }))
    }

    pub(super) fn name(&self) -> String {
        format!("contact[{}].resistance_m2_k_w", self.contact)
    }

    pub(super) fn render(&self) -> String {
        format!("{{\"kind\":\"contact-resistance\",\"contact\":{}}}", quote(&self.contact))
    }

    fn index(&self, base: &J) -> Result<usize> {
        let contacts = array(field(field(base, "solid")?, "contacts")?, "solid.contacts", 4096)?;
        let mut matches = contacts.iter().enumerate()
            .filter(|(_, row)| row.str_field("name") == Some(self.contact.as_str()));
        let (index, row) = matches.next().ok_or_else(|| bad(format!("unknown contact {}", self.contact)))?;
        if matches.next().is_some() { return Err(bad("uncertain contact name is ambiguous")); }
        let resistance = number(field(row, "resistance_m2_k_w")?, "contact resistance")?;
        admit(resistance)?;
        match (row.get("face_pairs"),row.get("nonmatching")) {
            (Some(pairs),None) => {
                if array(pairs,"contact.face_pairs",200_000)?.is_empty() {
                    return Err(bad("uncertain contact must own at least one declared face pair"));
                }
            }
            (None,Some(sides)) => {
                for key in ["side_a_faces","side_b_faces"] {
                    if array(field(sides,key)?,key,200_000)?.is_empty() {
                        return Err(bad("uncertain nonmatching contact requires both complete face sets"));
                    }
                }
            }
            _ => return Err(bad("uncertain contact requires exactly one matching or nonmatching trace declaration")),
        }
        Ok(index)
    }

    pub(super) fn validate(&self, base: &J) -> Result<()> {
        self.index(base).map(|_| ())
    }

    pub(super) fn apply(&self, base: &mut J, value: f64) -> Result<()> {
        admit(value).map_err(|error| model_failure(error.message))?;
        let index = self.index(base)?;
        let contacts = array_mut_path(base, &["solid", "contacts"])?;
        set_member_number(&mut contacts[index], "resistance_m2_k_w", value)
    }
}

fn admit(value: f64) -> Result<()> {
    if value.is_finite() && value > 0.0 && (1.0 / value).is_finite() {
        Ok(())
    } else {
        Err(bad("contact resistance must be finite and positive with representable reciprocal conductance"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/cooling-network/nonlinear-contact-pulse.json"));

    fn target() -> Target { Target { contact: "bondline".into() } }

    #[test]
    fn sampling_changes_only_the_named_physical_contact_law() {
        let mut base = J::parse(BASE).unwrap();
        let original = base.clone();
        target().validate(&base).unwrap();
        target().apply(&mut base, 0.02).unwrap();
        assert_eq!(array_path(&base, &["solid", "contacts"]).unwrap()[0]
            .f64_field("resistance_m2_k_w"), Some(0.02));
        target().apply(&mut base, 0.01).unwrap();
        assert_eq!(base, original, "no geometry, workload or material mutation");
    }

    #[test]
    fn missing_ambiguous_or_invalid_contacts_refuse() {
        let base = J::parse(BASE).unwrap();
        assert!(Target { contact: "missing".into() }.validate(&base).is_err());
        let mut duplicated = base.clone();
        let contacts = array_mut_path(&mut duplicated, &["solid", "contacts"]).unwrap();
        contacts.push(contacts[0].clone());
        assert!(target().validate(&duplicated).is_err());
        for value in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
            let mut request = base.clone();
            assert!(target().apply(&mut request, value).is_err());
            assert_eq!(request, base);
        }
        let ambiguous = J::parse(r#"{"kind":"contact-resistance","contact":"bondline","surface":"wall"}"#).unwrap();
        assert!(Target::parse(&ambiguous).is_err());
    }

    #[test]
    fn full_uq_plan_binds_contact_name_and_changes_real_request() {
        let base = J::parse(BASE).unwrap();
        let text = r#"{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":4,"wall_seconds":60,"qoi":{"kind":"transient-sampled-peak"},"correlation":{"kind":"independent"},"parameters":[{"target":{"kind":"contact-resistance","contact":"bondline"},"distribution":{"kind":"uniform","lo":0.005,"hi":0.02}}]}"#;
        let config = Config::parse(text, &base).unwrap();
        let sample = J::parse(&config.sample_request(&base, &[0.015]).unwrap()).unwrap();
        assert_eq!(array_path(&sample, &["solid", "contacts"]).unwrap()[0]
            .f64_field("resistance_m2_k_w"), Some(0.015));
        assert_eq!(sample.get("transient"), base.get("transient"));
        let parameters = config.render_parameters().unwrap();
        assert!(parameters.contains("contact[bondline].resistance_m2_k_w"));
        assert!(parameters.contains("m2 K/W"));
        assert!(Config::parse(&text.replace("\"lo\":0.005", "\"lo\":0"), &base).is_err());
    }

    #[test]
    fn nonmatching_samples_keep_both_meshes_and_all_geometry_admission_controls() {
        let original=J::parse(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
            "/../../examples/cooling-network/nonmatching-contact-hotspot.json"))).unwrap();
        let mut changed=original.clone();
        target().validate(&changed).unwrap();target().apply(&mut changed,0.025).unwrap();
        let old=&array_path(&original,&["solid","contacts"]).unwrap()[0];
        let new=&array_path(&changed,&["solid","contacts"]).unwrap()[0];
        assert_eq!(old.get("nonmatching"),new.get("nonmatching"));
        assert_eq!(new.f64_field("resistance_m2_k_w"),Some(0.025));
        target().apply(&mut changed,0.01).unwrap();assert_eq!(changed,original);
    }
}
