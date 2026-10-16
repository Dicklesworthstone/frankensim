"""Focused product-path tests: real importer, not a substitute physics solver."""
import base64
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("ftx_import", ROOT / "scripts/import_franken_textiles.py")
IMPORTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(IMPORTER)
FIXTURE = ROOT / "data/matdb/textile-examples/synthetic-hardening-sheet.json"
TERMS = dict(citation="Synthetic test data",license_text="Synthetic numerical fixture",
             redistribution_terms="Authored synthetic data",allow_estimates=True)

class NonlinearImport(unittest.TestCase):
    def setUp(self): self.record=json.loads(FIXTURE.read_bytes())
    def render(self,record=None,**kwargs):
        raw=json.dumps(record or self.record).encode()
        return IMPORTER.render(raw,self.record["id"],**(TERMS|dict(nonlinear_membrane=True)|kwargs))
    def test_curves_keep_knots_units_semantics_and_original_bytes(self):
        raw=FIXTURE.read_bytes()
        manifest,source,report=IMPORTER.render(raw,self.record['id'],**TERMS,nonlinear_membrane=True)
        rows=[l.split('\t') for l in source.decode().splitlines() if l and not l.startswith('#')]
        curves=[r for r in rows if r[0]=='curve']
        self.assertEqual(len(curves),3)
        self.assertEqual({r[3] for r in curves},{'sheet_second_piola_warp','sheet_second_piola_weft','sheet_second_piola_shear'})
        for r in curves:
            self.assertEqual(r[4:7],['sheet_green_strain','1','N/m']);self.assertEqual(r[-1],'linear')
        self.assertIn('0.025:200.0',curves[0][7])
        encoded=''.join(l.split(' ',2)[2] for l in source.decode().splitlines() if l.startswith('# original_json_base64 '))
        self.assertEqual(base64.b64decode(encoded),raw)
        self.assertEqual(report['source_sha256'],hashlib.sha256(raw).hexdigest())
        self.assertEqual(report['query']['sheet_green_strain'],0)
        self.assertTrue(all(r[2]=='-' for r in curves)) # estimates stay estimates
        self.assertIn(b'material-tsv-v1',manifest)
    def test_cubic_default_and_history_cannot_be_silently_reinterpreted(self):
        for interpolation in [None,'monotone_cubic','unknown']:
            r=copy.deepcopy(self.record)
            if interpolation is None:r['membrane_response'].pop('interpolation')
            else:r['membrane_response']['interpolation']=interpolation
            with self.assertRaises(ValueError):self.render(r)
        r=copy.deepcopy(self.record);r['membrane_relaxation']={}
        with self.assertRaises(ValueError):self.render(r)
    def test_linear_constants_are_replaced_not_added_to_the_curve(self):
        r=copy.deepcopy(self.record);r['properties']['warp_stiffness']=dict(r['properties']['areal_density'],value=1e9,unit='N/m')
        _,source,report=self.render(r)
        self.assertNotIn('\twarp_stiffness\t',source.decode())
        self.assertIn('warp_stiffness',report['source_only_properties'])
        with self.assertRaises(ValueError):self.render(r,selected={'warp_stiffness'})
        with self.assertRaises(ValueError):self.render(r,nonlinear_membrane=False)
    def test_bad_source_wrong_stress_measure_and_partial_curve_refuse(self):
        mutations=[lambda r:r['membrane_response'].update(model='nominal_stress/engineering_strain'),
            lambda r:r['membrane_response']['warp'].update(points=[[0,0],[.1,1],[.2,2]]),
            lambda r:r['membrane_response']['shear'].update(source='missing'),
            lambda r:r['membrane_response']['coupling'].update(unit='Pa'),
            lambda r:r.update(form='fiber')]
        for mutate in mutations:
            r=copy.deepcopy(self.record);mutate(r)
            with self.assertRaises(ValueError):self.render(r)
        with self.assertRaises(ValueError):self.render(allow_estimates=False)
    def test_scalar_only_default_path_remains_unchanged(self):
        r=copy.deepcopy(self.record);r.pop('membrane_response')
        raw=json.dumps(r).encode()
        a=IMPORTER.render(raw,r['id'],**TERMS)
        b=IMPORTER.render(raw,r['id'],**TERMS,nonlinear_membrane=False)
        self.assertEqual(a,b)
        with self.assertRaises(ValueError):self.render(r)
    def test_cli_creates_complete_pack_sources_without_overwrite(self):
        # Leave created artifacts for inspection; never delete repository/user files.
        parent=Path(tempfile.mkdtemp(prefix='ftx-nonlinear-',dir='/mnt/data'))
        out=parent/'imported'
        command=[sys.executable,str(ROOT/'scripts/import_franken_textiles.py'),str(FIXTURE),self.record['id'],str(out),
            '--citation','Synthetic test','--license','Synthetic fixture','--redistribution-terms','Authored synthetic data',
            '--allow-estimates','--nonlinear-membrane']
        first=subprocess.run(command,capture_output=True,text=True)
        self.assertEqual(first.returncode,0,first.stderr)
        self.assertEqual((out/'original.json').read_bytes(),FIXTURE.read_bytes())
        before=(out/'properties.tsv').read_bytes()
        self.assertTrue((out/'manifest.tsv').is_file())
        again=subprocess.run(command,capture_output=True,text=True)
        self.assertEqual(again.returncode,2)
        self.assertEqual((out/'properties.tsv').read_bytes(),before)

if __name__=='__main__':unittest.main()
