"""Guard against false stage completion, drift and path escape."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("checks", Path(__file__).parents[1] / "scripts/check_workspace.py")
checks = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checks)

class ReceiptTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / "handoff").mkdir()
        self.policy = {"A": {"requires": [], "outputs": ["evidence.txt"], "checks": ["actual_check"]}}

    def receipt(self):
        evidence = self.root / "evidence.txt"
        evidence.write_text("Example verification evidence for validator fixture, not a scientific result.")
        data = {"stage": "A", "status": "complete", "actor": "fixture", "completed_at": "2026-11-10",
                "outputs": [{"path": "evidence.txt", "sha256": checks.digest(evidence)}],
                "checks": [{"id": "actual_check", "passed": True, "evidence": ["evidence.txt"]}]}
        self.save(data)
        return data

    def save(self, data):
        (self.root / "handoff/A.json").write_text(json.dumps(data))

    def test_missing_receipt_rejected(self):
        with self.assertRaises(ValueError):
            checks.check_stage(self.root, "A", self.policy)

    def test_valid_evidence_accepted(self):
        self.receipt()
        checks.check_stage(self.root, "A", self.policy)

    def test_changed_evidence_rejected(self):
        self.receipt()
        (self.root / "evidence.txt").write_text("changed")
        with self.assertRaises(ValueError):
            checks.check_stage(self.root, "A", self.policy)

    def test_unsupported_pass_rejected(self):
        data = self.receipt()
        data["checks"][0]["evidence"] = []
        self.save(data)
        with self.assertRaises(ValueError):
            checks.check_stage(self.root, "A", self.policy)

    def test_missing_dependency_rejected(self):
        self.receipt()
        self.policy["A"]["requires"] = ["B"]
        self.policy["B"] = {"requires": [], "outputs": [], "checks": []}
        with self.assertRaises(ValueError):
            checks.check_stage(self.root, "A", self.policy)

    def test_path_escape_rejected(self):
        with self.assertRaises(ValueError):
            checks.confined(self.root, "../outside")

    def test_e1_import_only_receipt_is_rejected(self):
        (self.root / 'environment').mkdir()
        (self.root / 'design').mkdir()
        (self.root / 'environment/qualification.json').write_text(json.dumps({'passed':True,'biological_forwards':0,'models':{}}))
        (self.root / 'design/execution-config.json').write_text(json.dumps({'status':'frozen'}))
        (self.root / 'handoff/E1.json').write_text(json.dumps({'stage':'E1','status':'complete','actor':'fixture',
            'completed_at':'2026-11-10','outputs':[],'checks':[]}))
        with self.assertRaisesRegex(ValueError,'missing core model'):
            checks.check_stage(self.root,'E1',{'E1':{'requires':[],'outputs':[],'checks':[]}})

    def test_e1_invented_tolerance_is_rejected(self):
        (self.root / 'environment').mkdir()
        (self.root / 'design').mkdir()
        model={'eta':{'1':.1},'errors':{'1':{'repeat_max':0.,'cpu_mps_max':0.}}}
        q={'passed':True,'biological_forwards':0,'models':{'esm2_650m':model,'profam1':model},'resources':{'fits_before_freeze':True}}
        (self.root / 'environment/qualification.json').write_text(json.dumps(q))
        (self.root / 'design/execution-config.json').write_text(json.dumps({'status':'frozen','models':q['models']}))
        (self.root / 'handoff/E1.json').write_text(json.dumps({'stage':'E1','status':'complete','actor':'fixture',
            'completed_at':'2026-11-10','outputs':[],'checks':[]}))
        with self.assertRaisesRegex(ValueError,'tolerance differs'):
            checks.check_stage(self.root,'E1',{'E1':{'requires':[],'outputs':[],'checks':[]}})

if __name__ == "__main__":
    unittest.main()
