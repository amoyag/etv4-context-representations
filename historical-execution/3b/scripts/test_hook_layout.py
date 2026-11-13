"""Exercise Adapter.forward on a tiny untrained CPU model; no checkpoint loading."""
import torch
import esm
import run_experiments as r
# Match the adapter's fixed width, with identity blocks to keep this test tiny.
class Block(torch.nn.Module):
    def forward(self, x, **kwargs):
        return x + 0.01, None
model = esm.model.esm2.ESM2(num_layers=0, embed_dim=2560, attention_heads=20)
model.layers = torch.nn.ModuleList([Block() for _ in range(36)])
model.eval().requires_grad_(False)
a = r.Adapter.__new__(r.Adapter)
a.torch, a.model, a.alphabet, a.device = torch, model, model.alphabet, 'cpu'
for sequence in ['ACD', 'ACDEFGH']:
    arrays, mapping, seconds = a.forward(sequence)
    assert set(arrays) == {'19','24','36'}
    assert all(x.shape == (len(sequence),2560) for x in arrays.values())
    assert mapping['target_length'] == len(sequence)
print('Adapter hook layout: two sequence lengths passed')
