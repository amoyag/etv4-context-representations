#!/usr/bin/env python3
"""Frozen-model adapters. E1 uses artificial inputs; E2 requires explicit activation."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
AA = set('ACDEFGHIKLMNPQRSTVWY')
LAYERS = {'esm2_650m': [17, 22, 33], 'profam1': [8, 11, 16]}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read(name):
    return json.loads((ROOT / name).read_text())


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def validate_sequences(target, prefixes):
    for seq in [*prefixes, target]:
        if not seq or set(seq) - AA:
            raise ValueError('Only nonempty canonical amino-acid sequences are allowed')


def native_indices(native_positions, native_start, target_length, token_start=0):
    positions = list(native_positions)
    if len(set(positions)) != len(positions):
        raise ValueError('Duplicate native coordinate')
    offsets = [i - native_start for i in positions]
    if any(i < 0 or i >= target_length for i in offsets):
        raise ValueError('Native position outside target')
    return [token_start + i for i in offsets]


class Adapter:
    """Capture post-block residuals, before final model-wide normalisation."""
    def __init__(self, model_id, device='cpu', threads=4):
        import torch
        self.torch = torch
        torch.set_num_threads(threads)
        torch.set_num_interop_threads(1)
        torch.manual_seed(19001)
        self.model_id, self.device = model_id, device
        self.layers = LAYERS[model_id]
        if device == 'mps' and not torch.backends.mps.is_available():
            raise RuntimeError('MPS unavailable in this process')
        assets = read('environment/model-assets-macmini.json')['files']
        self.asset_checks = []
        for item in assets:
            relevant = item['model'].startswith('ESM') if model_id == 'esm2_650m' else item['model'].startswith('ProFam')
            if relevant:
                actual = digest(ROOT / item['path'])
                if actual != item['sha256']:
                    raise ValueError('Checkpoint/config hash mismatch: ' + item['path'])
                self.asset_checks.append({'path': item['path'], 'sha256': actual})
        start = time.perf_counter()
        if model_id == 'esm2_650m':
            import esm
            checkpoint = ROOT / 'models/esm2/esm2_t33_650M_UR50D.pt'
            # Official fair-esm checkpoint stores its architecture in argparse.Namespace.
            # Keep weights_only=True; allow only this inert configuration container.
            with torch.serialization.safe_globals([argparse.Namespace]):
                self.model, self.tokenizer = esm.pretrained.load_model_and_alphabet_local(str(checkpoint))
            self.model.eval().requires_grad_(False).to(device=device, dtype=torch.float32)
            self.blocks = self.model.layers
            self.final_norm = self.model.emb_layer_norm_after
            self.width = self.model.embed_dim
        else:
            from profam.checkpoint import load_model
            original_load = torch.nn.Module.load_state_dict
            load_checks = []
            def audited_load(module, *args, **kwargs):
                result = original_load(module, *args, **kwargs)
                state_dict = args[0] if args else kwargs['state_dict']
                alias = 'model.token_embedder.weight'
                duplicate = alias in state_dict and torch.equal(state_dict[alias], state_dict['model.model.embed_tokens.weight'])
                load_checks.append({'missing_keys': list(result.missing_keys), 'unexpected_keys': list(result.unexpected_keys),
                                    'token_embedder_is_exact_input_embedding_duplicate': duplicate})
                return result
            torch.nn.Module.load_state_dict = audited_load
            try:
                self.model = load_model(checkpoint=str(ROOT / 'models/profam-1/checkpoints/last.ckpt'),
                                        device=device, dtype='float32', attn_implementation='sdpa', auto_download=False)
            finally:
                torch.nn.Module.load_state_dict = original_load
            assert load_checks and all(not x['missing_keys'] and (not x['unexpected_keys'] or
                (x['unexpected_keys'] == ['model.token_embedder.weight'] and x['token_embedder_is_exact_input_embedding_duplicate']))
                for x in load_checks), load_checks
            self.state_dict_audit = load_checks
            self.model.requires_grad_(False)
            self.tokenizer = self.model.tokenizer
            self.blocks = self.model.model.model.layers
            self.final_norm = self.model.model.model.norm
            self.width = self.model.model.config.hidden_size
            assert self.model.model.model.rotary_emb.rope_type == 'llama3'
        self.load_seconds = time.perf_counter() - start
        assert len(self.blocks) == self.layers[-1]
        assert not self.model.training and not any(p.requires_grad for p in self.model.parameters())
        self.parameter_count = sum(p.numel() for p in self.model.parameters())

    def encode(self, target, prefixes=()):
        import numpy as np
        torch = self.torch
        prefixes = list(prefixes)
        validate_sequences(target, prefixes)
        if self.model_id == 'esm2_650m':
            if prefixes:
                raise ValueError('ESM has no family prompt')
            _, _, tokens = self.tokenizer.get_batch_converter()([('target', target)])
            expected = [self.tokenizer.cls_idx] + [self.tokenizer.get_idx(a) for a in target] + [self.tokenizer.eos_idx]
            token_start = 1
        else:
            from profam.api import _build_protein_document
            seqs = prefixes + [target]
            encoded = self.tokenizer.encode(_build_protein_document(seqs), padding=False, add_final_sep=True)
            tokens = torch.tensor(np.asarray(encoded['input_ids']), dtype=torch.long).unsqueeze(0)
            expected = [self.tokenizer.convert_tokens_to_ids(t) for t in ['[start-of-document]', '[RAW]']]
            for seq in seqs:
                expected.extend(self.tokenizer.convert_tokens_to_ids(list(seq)))
                expected.append(self.tokenizer.sep_token_id)
            token_start = 2 + sum(len(s) + 1 for s in prefixes)
        assert tokens.tolist()[0] == expected, 'Token identity/BOS/SEP/EOS mismatch'
        assert len(expected) == token_start + len(target) + 1, 'Length mismatch'
        tokens = tokens.to(self.device)
        positions = torch.arange(tokens.shape[1], device=self.device).unsqueeze(0)
        if self.model_id == 'profam1':
            actual = self.model.get_position_ids_for_model_forward(tokens, None)
            assert torch.equal(actual, positions), 'Document positions reset or differ'
        mapping = {'token_start': token_start, 'target_length': len(target), 'tokens': tokens.shape[1],
                   'token_ids': expected, 'position_ids': positions.cpu().tolist()[0],
                   'target_token_indices': list(range(token_start, token_start + len(target)))}
        return tokens, mapping

    def forward(self, target, prefixes=()):
        import numpy as np
        torch = self.torch
        tokens, mapping = self.encode(target, prefixes)
        captured, hooks, norm_input = {}, [], {}
        def capture(layer):
            def hook(module, args, output):
                x = output[0] if isinstance(output, tuple) else output
                if self.model_id == 'esm2_650m':
                    x = x.transpose(0, 1)
                captured[layer] = x.detach().clone()
            return hook
        def norm_hook(module, args):
            x = args[0]
            if self.model_id == 'esm2_650m':
                x = x.transpose(0, 1)
            norm_input['x'] = x.detach().clone()
        for layer in self.layers:
            hooks.append(self.blocks[layer - 1].register_forward_hook(capture(layer)))
        hooks.append(self.final_norm.register_forward_pre_hook(norm_hook))
        if self.device == 'mps':
            torch.mps.synchronize()
        t0 = time.perf_counter()
        try:
            with torch.inference_mode():
                if self.model_id == 'esm2_650m':
                    out = self.model(tokens, repr_layers=self.layers[:-1], return_contacts=False)
                    for layer in self.layers[:-1]:
                        assert torch.equal(captured[layer], out['representations'][layer])
                else:
                    out = self.model(input_ids=tokens, use_cache=False, output_hidden_states=False)
                assert torch.equal(captured[self.layers[-1]], norm_input['x']), 'Last-block output is not pre-normalisation'
                assert set(captured) == set(self.layers)
                if self.device == 'mps':
                    torch.mps.synchronize()
                elapsed = time.perf_counter() - t0
                start = mapping['token_start']
                arrays = {str(k): v[0, start:start + len(target)].cpu().numpy().copy() for k, v in captured.items()}
                assert all(a.shape == (len(target), self.width) and np.isfinite(a).all() for a in arrays.values())
        finally:
            for hook in hooks:
                hook.remove()
        return arrays, mapping, elapsed


def require_e2():
    state = read('WORK_STATE.json')
    if not state.get('active_execution') or 'E2' not in state.get('activation_scope', []) or state['stages']['E2']['status'] != 'in_progress':
        raise RuntimeError('E2 is not explicitly activated; biological inference prohibited')
    from check_workspace import check_stage
    check_stage(ROOT, 'E1', read('design/stages.json')['stages'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', choices=LAYERS, required=True)
    parser.add_argument('--run-id', required=True, help='New immutable E2 execution directory name')
    args = parser.parse_args()
    require_e2()
    if Path(args.run_id).name != args.run_id or args.run_id in ('.', '..'):
        raise ValueError('Invalid execution ID')
    import numpy as np
    import fcntl
    lock = (ROOT / '.model-execution.lock.tmp').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    config = read('design/execution-config.json')
    if config['status'] != 'frozen':
        raise RuntimeError('Execution configuration is not frozen')
    dest = ROOT / 'results/activations' / args.run_id
    dest.mkdir(parents=True, exist_ok=False)
    adapter = Adapter(args.model, config['models'][args.model]['device'], config['threads'])
    inputs = read('data/derived/model-inputs.json')
    prompts = read('data/derived/prompt-panel.json')
    manifest = []
    for row in read('design/run-plan.json')['runs']:
        if row['model'] != args.model:
            continue
        x = inputs[row['input_id']]
        arrays, mapping, seconds = adapter.forward(x['sequence'], prompts[row['prompt_id']]['sequences'])
        assert mapping['tokens'] == row['expected_tokens']
        path = dest / (row['run_id'] + '.npz')
        np.savez(path, **arrays)
        manifest.append({**row, 'activation_file': str(path.relative_to(ROOT)), 'sha256': digest(path),
                         'seconds': seconds, 'mapping': mapping, 'native_start': x['native_start'],
                         'execution_config_sha256': digest(ROOT / 'design/execution-config.json')})
        write_json(dest / 'manifest.json', manifest)
        print(row['run_id'], seconds, flush=True)


if __name__ == '__main__':
    main()
