"""Source adapter regressions, not an alternate simulation/reference engine."""
import base64
import copy
import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("ftx_import", ROOT / "scripts/import_franken_textiles.py")
bridge = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bridge)


def record():
    return json.loads((ROOT / "data/matdb/textile-examples/synthetic-sheet.json").read_text())["materials"][0]


def render(material=None, **options):
    obj = record() if material is None else material
    raw = json.dumps(obj, ensure_ascii=False).encode()
    opts = dict(citation="synthetic source test",license_text="test-only",redistribution_terms="synthetic unit-test data",allow_estimates=True)
    opts.update(options)
    return raw, bridge.render(raw,obj["id"],**opts)


class ImportTests(unittest.TestCase):
    def test_original_bytes_qualifiers_and_measured_vs_estimated_backing(self):
        m = record();m["properties"]["areal_density"]["evidence"]="measured"
        raw,(manifest,source,report)=render(m)
        chunks=[l.split(" ",2)[2] for l in source.decode().splitlines() if l.startswith("# original_json_base64 ")]
        self.assertEqual(base64.b64decode("".join(chunks)),raw)
        self.assertIn(report["source_sha256"].encode(),manifest)
        rows=[l.split("\t") for l in source.decode().splitlines() if l.startswith("scalar\t")]
        self.assertNotEqual(next(r for r in rows if r[3]=="areal_density")[2],"-")
        self.assertEqual(next(r for r in rows if r[3]=="warp_stiffness")[2],"-")
        self.assertTrue(all(r[-1]=="constant" for r in rows))

    def test_sheet_units_do_not_become_bulk_properties(self):
        _,(_,source,_)=render();rows=[l.split("\t") for l in source.decode().splitlines() if l.startswith("scalar\t")]
        props={r[3]:(float(r[4]),r[5]) for r in rows}
        self.assertEqual(props["areal_density"],(.2,"kg/m2"))
        self.assertEqual(props["thickness"],(.001,"m"))
        self.assertEqual(props["warp_stiffness"],(8000.,"N/m"))
        self.assertEqual(props["bending_warp"],(3e-5,"N*m"))
        self.assertNotIn("density",props);self.assertNotIn("young_modulus",props)
        m=record();m["properties"]["warp_stiffness"]["unit"]="Pa"
        with self.assertRaises(ValueError):render(m)

    def test_temperature_tables_keep_domain_and_do_not_turn_sd_into_confidence(self):
        _,(_,source,_)=render(interpolate_temperature=True)
        rows=[l.split("\t") for l in source.decode().splitlines() if l.startswith("curve\t")]
        cp=next(r for r in rows if r[3]=="specific_heat_capacity")
        self.assertEqual(cp[4:7],["T","K","J/kg/K"])
        self.assertEqual(cp[7],"293.15:1200.0,333.15:1600.0")
        self.assertEqual(cp[8],"linear")
        self.assertIn(b"standard_deviation",source)
        self.assertTrue(all("\tunstated\t" in l for l in source.decode().splitlines() if l.startswith("uncertainty\t")))
        _,(_,exact,_)=render();self.assertIn(b"293.15:1200.0,333.15:1600.0\ttabulated",exact)

    def test_missing_heat_capacity_stays_missing_and_competing_laws_refuse(self):
        m=record();del m["thermal_response"]["specific_heat"]
        _,(_,source,report)=render(m,selected={"thermal_conductivity"})
        self.assertNotIn("specific_heat",report["imported"])
        self.assertFalse(any(l.startswith("curve\t") and "specific_heat_capacity" in l for l in source.decode().splitlines()))
        m["properties"]["thermal_conductivity"]=copy.deepcopy(m["properties"]["thickness"])
        with self.assertRaises(ValueError):render(m)

    def test_history_curvature_and_wrong_specimen_are_not_silently_flattened(self):
        m=record();m["membrane_relaxation"]={"model":"some nonempty model"}
        with self.assertRaises(ValueError):render(m)
        render(m,selected={"thermal_conductivity"})
        m=record();m["properties"]["preferred_curvature_warp"]={"value":1.0,"unit":"1/m"}
        with self.assertRaises(ValueError):render(m)
        m=record();m["properties"]["thickness"]["conditions"]["conditioning"]="other moisture"
        with self.assertRaises(ValueError):render(m)
        render(m,selected={"specific_heat"})

    def test_estimates_polynomials_and_wrong_source_hash_refuse(self):
        with self.assertRaises(ValueError):render(allow_estimates=False)
        with self.assertRaises(ValueError):render(expected_sha256="0"*64)
        m=record();m["thermal_response"]["conductivity"]["polynomial"]={"coefficients":[1,2]}
        with self.assertRaises(ValueError):render(m)

    def test_duplicate_json_unknown_schema_or_id_and_injection(self):
        for raw in (b'{"schema":"bad","schema":"ftx-textile-catalog/1"}',b'{"schema":"future"}'):
            with self.assertRaises(ValueError):bridge.render(raw,"x",citation="c",license_text="l",redistribution_terms="r")
        m=record();m["id"]="tab\tnew\nrecord";m["properties"]["thickness"]["method"]="another\tline\nmethod"
        _,(_,source,_)=render(m)
        self.assertTrue(all(len(l.split("\t"))==5 for l in source.decode().splitlines() if l.startswith("observation\t")))
        with self.assertRaises(ValueError):render(citation="inject\nsource\tbad")

    def test_thermal_only_catalog_record_uses_canonical_default_empty_properties(self):
        # Actual FSRI catalog records omit the serde-default properties map.
        m=record();m.pop("properties");m["form"]="fill"
        m["thermal_response"].pop("specific_heat")
        _,(_,source,report)=render(m,selected={"thermal_conductivity"})
        self.assertEqual(report["imported"],["thermal_conductivity"])
        self.assertNotIn(b"scalar\t",source)

    def test_selected_properties_are_required_not_silently_omitted(self):
        with self.assertRaises(ValueError):render(selected={"air_permeability"})
        m=record();m["form"]="fiber"
        with self.assertRaises(ValueError):render(m)
        render(m,selected={"thermal_conductivity"})


if __name__ == "__main__":unittest.main()
