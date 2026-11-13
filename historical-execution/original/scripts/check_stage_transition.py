#!/usr/bin/env python3
"""Require an evidence-linked predecessor assessment before downstream work.

Mechanical validation cannot certify the author's semantic judgment. The
coordinator must read and assess the cited outputs, not merely fill this schema.
Does not activate a stage, load a model or change WORK_STATE/receipts.
"""
import argparse
from datetime import datetime
import json
from pathlib import Path
import sys
from check_workspace import ROOT, check_stage, check_final, confined, digest, read


def require(condition, message):
    if not condition:
        raise ValueError(message)


def explanation(value, label):
    require(isinstance(value, str) and len(value.strip()) >= 30,
            label + ': requires a specific observation and consequence, not a short approval')
    require(value.strip().lower() not in {'all checks passed and everything is ready', 'all previous outputs have been verified'},
            label + ': generic approval is not an assessment')


def ancestors(stage, stages):
    seen=set()
    def visit(current, trail):
        require(current not in trail, 'Stage dependency cycle')
        for parent in stages[current]['requires']:
            if parent not in seen:
                visit(parent, trail | {current})
                seen.add(parent)
    visit(stage,set())
    return sorted(seen)


def template(root, stage):
    stages=read(root/'design/stages.json')['stages']
    spec=read(root/'design/stage-intake-policy.json')['stages'][stage]
    def sha_if_exists(path):
        p=confined(root,path)
        return digest(p) if p.is_file() else None
    return {'schema_version':1,'stage':stage,'status':'draft','actor':None,'evaluated_at':None,
        'intake_policy_sha256':digest(root/'design/stage-intake-policy.json'),
        'upstream_receipts':[{'stage':s,'path':f'handoff/{s}.json','sha256':sha_if_exists(f'handoff/{s}.json')} for s in ancestors(stage,stages)],
        'evidence':[{'path':p,'sha256':sha_if_exists(p)} for p in spec['evidence']],
        'assessments':[{'topic':t,'observed_result':None,'limitation':None,'effect_on_current_stage':None,'evidence_paths':[]} for t in spec['topics']],
        'decision':{'outcome':'not_evaluated','rationale':None,'scope_unchanged':None,'open_material_decisions':[]}}


def validate_assessment(root, stage, assessment, spec, stages):
    require(assessment.get('schema_version')==1 and assessment.get('stage')==stage, 'Wrong assessment schema/stage')
    require(assessment.get('status')=='reviewed', 'Intake is not reviewed')
    require(bool(assessment.get('actor')), 'Missing evaluator identity')
    instant=datetime.fromisoformat(assessment.get('evaluated_at') or '')
    require(instant.tzinfo is not None, 'Intake time requires timezone')
    require(assessment.get('intake_policy_sha256')==digest(root/'design/stage-intake-policy.json'), 'Intake policy changed')
    entries=assessment.get('upstream_receipts',[])
    required=ancestors(stage,stages)
    require(len(entries)==len(required) and {r['stage'] for r in entries}==set(required), 'Missing/duplicate ancestor receipt')
    for entry in entries:
        require(entry['path']==f"handoff/{entry['stage']}.json", 'Unexpected ancestor receipt path')
        require(digest(confined(root,entry['path']))==entry['sha256'], 'Changed ancestor receipt')
    evidence=assessment.get('evidence',[])
    evidence_paths=[e['path'] for e in evidence]
    require(len(evidence_paths)==len(set(evidence_paths)), 'Duplicate evidence')
    require(set(spec['evidence']).issubset(evidence_paths), 'Required predecessor output was not assessed')
    for entry in evidence:
        path=confined(root,entry['path'])
        require(path.is_file() and digest(path)==entry['sha256'], 'Changed/missing assessed output: '+entry['path'])
    rows=assessment.get('assessments',[])
    require(len(rows)==len(spec['topics']) and {r['topic'] for r in rows}==set(spec['topics']), 'Missing/duplicate assessment topic')
    used=set()
    for row in rows:
        for field in ['observed_result','limitation','effect_on_current_stage']:
            explanation(row.get(field), row['topic']+'/'+field)
        refs=row.get('evidence_paths',[])
        require(bool(refs) and set(refs).issubset(evidence_paths), 'Unidentified evidence in assessment')
        used.update(refs)
    require(set(spec['evidence']).issubset(used), 'Required outputs listed but not used in any assessment')
    decision=assessment.get('decision',{})
    require(decision.get('outcome')=='proceed', 'Dependent stage cannot proceed')
    require(decision.get('scope_unchanged') is True, 'Material scope change needs an explicit reviewed amendment')
    require(decision.get('open_material_decisions')==[], 'Unresolved material decision blocks dependent work')
    explanation(decision.get('rationale'),'decision/rationale')


def check_transition(root, stage, complete=False, seen=None):
    seen=set() if seen is None else set(seen)
    require(stage not in seen,'Intake dependency cycle');seen.add(stage)
    policy=read(root/'design/stage-intake-policy.json')['stages']
    stages=read(root/'design/stages.json')['stages']
    spec=policy[stage]
    require(stages[stage]['requires']==spec['requires'],'Stage/intake dependency mismatch')
    # Missing intake fails before expensive ancestor rechecks. No template counts as evidence.
    path=root/'handoff/intake'/f'{stage}.json'
    require(path.is_file(),'Missing intake assessment: '+str(path.relative_to(root)))
    for dep in spec['requires']:
        if dep in policy:
            check_transition(root,dep,complete=True,seen=seen)
        else:
            check_stage(root,dep,stages)
    validate_assessment(root,stage,read(path),spec,stages)
    if complete:
        receipt=check_stage(root,stage,stages)
        output_map={r['path']:r['sha256'] for r in receipt['outputs']}
        name=f'handoff/intake/{stage}.json'
        require(output_map.get(name)==digest(path),'Intake not preserved in completion receipt')
        check=next((r for r in receipt['checks'] if r['id']=='predecessor_output_evaluated'),{})
        log=f'logs/{stage}-intake.txt'
        require({name,log}.issubset(check.get('evidence',[])),'Completion must cite intake assessment and gate log')
        if stage=='R0': check_final(root,stages)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',required=True,choices=['E2','W0','W1','R0'])
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--template',action='store_true',help='Print an unapproved draft schema; does not write or activate anything')
    mode.add_argument('--complete',action='store_true',help='Also verify the completed stage and all ancestor intake records')
    args=parser.parse_args()
    try:
        if args.template:
            print(json.dumps(template(ROOT,args.stage),indent=2))
            return 0
        check_transition(ROOT,args.stage,args.complete)
    except (ValueError,KeyError,FileNotFoundError,TypeError) as exc:
        print('NOT READY: '+str(exc),file=sys.stderr)
        return 2
    print('PASS: predecessor identities and recorded assessment checked; semantic support requires coordinator review.')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
