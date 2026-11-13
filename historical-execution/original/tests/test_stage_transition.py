"""Synthetic intake fixtures: prevent generic approval, missing evidence and drift."""
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import check_stage_transition as gate


class TransitionTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        for p in ['design','handoff/intake']:(self.root/p).mkdir(parents=True)
        self.stages={'E1':{'requires':[]},'E2':{'requires':['E1']}}
        self.spec={'requires':['E1'],'evidence':['evidence.txt'],'topics':['visibility']}
        self.save('design/stages.json',{'stages':self.stages})
        self.save('design/stage-intake-policy.json',{'schema_version':1,'stages':{'E2':self.spec}})
        self.save('handoff/E1.json',{'fixture':'Not a real qualification receipt'})
        (self.root/'evidence.txt').write_text('Synthetic test evidence only; no real model result.')
        self.record=gate.template(self.root,'E2')
        self.record.update(status='reviewed',actor='unittest synthetic fixture',evaluated_at=datetime.now(timezone.utc).isoformat())
        self.record['assessments']=[{'topic':'visibility','observed_result':'This synthetic fixture represents an earlier-state visibility test.',
            'limitation':'The fixture establishes no biological result and is for validation only.',
            'effect_on_current_stage':'The fixture requires causal evidence to be cited before dependent execution.',
            'evidence_paths':['evidence.txt']}]
        self.record['decision']={'outcome':'proceed','scope_unchanged':True,'open_material_decisions':[],
            'rationale':'This decision belongs to a temporary synthetic validator test only.'}
    def save(self,path,data):
        (self.root/path).write_text(json.dumps(data))
    def validate(self):
        gate.validate_assessment(self.root,'E2',self.record,self.spec,self.stages)
    def test_specific_evidence_linked_assessment_accepted(self):self.validate()
    def test_unapproved_template_rejected(self):
        self.record['status']='draft'
        with self.assertRaisesRegex(ValueError,'not reviewed'):self.validate()
    def test_generic_pass_rejected(self):
        self.record['assessments'][0]['effect_on_current_stage']='All checks passed and everything is ready'
        with self.assertRaisesRegex(ValueError,'generic approval'):self.validate()
    def test_changed_output_rejected(self):
        (self.root/'evidence.txt').write_text('Changed after the intake')
        with self.assertRaisesRegex(ValueError,'Changed/missing'):self.validate()
    def test_changed_receipt_rejected(self):
        self.save('handoff/E1.json',{'changed':'upstream'})
        with self.assertRaisesRegex(ValueError,'Changed ancestor'):self.validate()
    def test_listed_but_unassessed_output_rejected(self):
        self.record['assessments'][0]['evidence_paths']=[]
        with self.assertRaisesRegex(ValueError,'Unidentified evidence'):self.validate()
    def test_material_decision_blocks_dependent_work(self):
        self.record['decision']['open_material_decisions']=['Unapproved change to the biological panel']
        with self.assertRaisesRegex(ValueError,'Unresolved material'):self.validate()
    def test_missing_intake_stops_before_dispatch(self):
        with patch.object(gate,'check_stage') as check:
            with self.assertRaisesRegex(ValueError,'Missing intake'):gate.check_transition(self.root,'E2')
            check.assert_not_called()
    def test_dispatch_rechecks_actual_dependency(self):
        self.save('handoff/intake/E2.json',self.record)
        with patch.object(gate,'check_stage',side_effect=ValueError('upstream failed')):
            with self.assertRaisesRegex(ValueError,'upstream failed'):gate.check_transition(self.root,'E2')
    def test_ancestor_intake_also_required(self):
        self.stages['W1']={'requires':['E2']}
        self.save('design/stages.json',{'stages':self.stages})
        p=gate.read(self.root/'design/stage-intake-policy.json');p['stages']['W1']={**self.spec,'requires':['E2']}
        self.save('design/stage-intake-policy.json',p)
        self.save('handoff/intake/W1.json',{'not':'a substitute for missing E2 intake'})
        with self.assertRaisesRegex(ValueError,'Missing intake.*E2'):gate.check_transition(self.root,'W1')


if __name__=='__main__':unittest.main()
